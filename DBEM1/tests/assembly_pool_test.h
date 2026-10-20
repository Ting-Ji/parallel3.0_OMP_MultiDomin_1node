// Synthetic row work tests scheduling without geometry or large solves.
#include <atomic>
static pthread_mutex_t queueTestMutex;
static pthread_cond_t queueTestReady;
static bool otherDomainStarted;
static bool requireOverlap;
static bool exportFirstCleared;
static int injectedTaskFailure;
static std::atomic<int> activeSlots[2], visits[20], finishedTasks;
static const MultiDomainModel* queueTestModel;

static void ObserveAssemblyTask(const AssemblyTask& task, long slot)
{
    CHECK(slot >= 0 && slot < 2);
    CHECK(activeSlots[slot].fetch_add(1) == 0);
    if (task.kind == AssemblyTaskKind::Output) {
        if (injectedTaskFailure == 3 && task.phase == 2 && task.sourceBegin > 0) {
            RequireThreadSuccess(pthread_mutex_lock(&queueTestMutex), "test export lock");
            while (!exportFirstCleared)
                RequireThreadSuccess(pthread_cond_wait(&queueTestReady, &queueTestMutex), "test export wait");
            RequireThreadSuccess(pthread_mutex_unlock(&queueTestMutex), "test export unlock");
            throw std::bad_alloc();
        }
        return;
    }
    CHECK(&DSquareElement::ActiveDynaMat() == &queueTestModel->domains[task.domainId].mat);
    if (requireOverlap) {
        RequireThreadSuccess(pthread_mutex_lock(&queueTestMutex), "test lock");
        if (task.domainId == 1) {
            otherDomainStarted = true;
            RequireThreadSuccess(pthread_cond_broadcast(&queueTestReady), "test signal");
        }
        if (task.domainId == 0 && task.sourceBegin == 0)
            while (!otherDomainStarted)
                RequireThreadSuccess(pthread_cond_wait(&queueTestReady, &queueTestMutex), "test wait");
        RequireThreadSuccess(pthread_mutex_unlock(&queueTestMutex), "test unlock");
    }
    if (task.domainId == 0 && task.sourceBegin == 0) {
        if (injectedTaskFailure == 1) throw std::bad_alloc();
        if (injectedTaskFailure == 2) throw std::runtime_error("injected task structure failure");
    }
}
static void FinishedAssemblyTask(const AssemblyTask& task, long slot, bool restored)
{
    CHECK(activeSlots[slot].fetch_sub(1) == 1);
    CHECK(restored);
    if (task.kind == AssemblyTaskKind::Output && task.phase == 2 && task.sourceBegin == 0) {
        RequireThreadSuccess(pthread_mutex_lock(&queueTestMutex), "test export done lock");
        exportFirstCleared = true;
        RequireThreadSuccess(pthread_cond_broadcast(&queueTestReady), "test export done");
        RequireThreadSuccess(pthread_mutex_unlock(&queueTestMutex), "test export done unlock");
    }
    ++finishedTasks;
}
static void* SyntheticStep0(void* data)
{
    auto& c = *static_cast<Step0WorkerContext*>(data);
    c.ok = 1;
    CHECK(c.coefSlot == AssistCoef[c.threadId]);
    CHECK(c.uSlot == AssistUij[c.threadId] && c.tSlot == AssistTij[c.threadId]);
    CHECK(&DSquareElement::ActiveDynaMat() == &c.domain->mat);
    for (long i = c.sourceBegin; i < c.sourceEnd; ++i) {
        const long node = c.elements[c.domain->elementIds[i]].m_nodeID[0];
        CHECK(visits[node].fetch_add(1) == 0);
        const long row = c.dofMap->GetEquationRow(c.domain->id, node);
        const double block[9] = {double(node+1),0,0,0,2,0,0,0,3};
        c.unknownBuilder->AddBlock(row, c.dofMap->GetU(c.domain->id,node).globalBlock, 1, block);
    }
    return nullptr;
}
static void* SyntheticHistory(void* data)
{
    auto& c = *static_cast<MultiDomainHistoryThreadArg*>(data);
    c.ok = 1;
    CHECK(&DSquareElement::ActiveDynaMat() == &c.domain->mat);
    for (long i = c.sourceBegin; i < c.sourceEnd; ++i) {
        const long node = c.elements[c.domain->elementIds[i]].m_nodeID[0];
        CHECK(visits[node].fetch_add(1) == 0);
        const long row = c.dofMap->GetEquationRow(c.domain->id,node);
        const double block[9] = {double(node+c.step),0,0,0,2,0,0,0,3};
        c.mtsBuilder->AddBlock(row,c.dofMap->GetHistoryU(c.domain->id,node).globalBlock,1,block);
        c.mgsBuilder->AddBlock(row,c.dofMap->GetHistoryT(c.domain->id,node).globalBlock,1,block);
    }
    return nullptr;
}
static void ResetQueueVisits()
{
    for (auto& v : visits) v = 0;
    finishedTasks = 0;
    otherDomainStarted = false;
    for (auto& active : activeSlots) CHECK(active == 0);
}
static void TestAssemblyQueue()
{
    TestStep0Scratch scratch;
    DSquareElement elements[20];
    MultiDomainModel model;
    model.elementCount = model.nodeCount = 20;
    for (int d = 0; d < 2; ++d) {
        Domain domain = {};
        domain.id = d; domain.eleBegin = d ? 17 : 0; domain.eleCount = d ? 3 : 17;
        domain.mat.C1 = 10+d;
        for (long e = domain.eleBegin; e < domain.eleBegin+domain.eleCount; ++e) {
            domain.elementIds.push_back(e);
            elements[e].BCID = 456; elements[e].InfEleFlag = -1; elements[e].m_nodeID[0] = e;
        }
        model.domains.push_back(domain);
    }
    GlobalDofMap dofs; CHECK(dofs.Build(model,elements,20));
    MultiDomainInterfaceCache cache(model,elements);
    queueTestModel = &model;
    CHECK(pthread_mutex_init(&queueTestMutex,nullptr) == 0);
    CHECK(pthread_cond_init(&queueTestReady,nullptr) == 0);
    g_assemblyObserveTask = ObserveAssemblyTask;
    g_assemblyTaskFinished = FinishedAssemblyTask;
    g_step0ExecuteRows = SyntheticStep0;
    g_historyExecuteRows = SyntheticHistory;
    g_assemblyCreateThread = InjectCreate;
    failAt = -1; createCalls = 0; injectedTaskFailure = 0; requireOverlap = true;
    {
        MultiDomainAssemblyPool pool(elements,model,dofs,cache,nullptr,nullptr,2);
        CHECK(pool.Initialize()); CHECK(createCalls == 2);
        CHECK(pool.SourceTasks().size() == 9); // ceil(17/8)=3, ceil(3/8)=1.
        CHECK(pool.SourceTasks()[0].domainId == 0 && pool.SourceTasks()[1].domainId == 1);
        for (int repeat = 0; repeat < 2; ++repeat) {
            ResetQueueVisits();
            CHECK(pool.RunBatch(pool.SourceTasks()));
            for (auto& v : visits) CHECK(v == 1);
            CHECK(otherDomainStarted && finishedTasks == 9);
            MultiDomainCCSRBuilder unknown(20,20), known(20,20);
            pool.MergeInto(unknown,known);
            CHECK(unknown.NonZeroBlocks() == 20);
            for (long row=0; row<20; ++row) CHECK((*unknown.FindBlock(row,row))[0] == row+1);
        }
        for (long step : {1L,2L}) {
            ResetQueueVisits();
            MultiDomainCCSRBuilder t(20,20,false,MultiDomainCCSRBuilder::IndependentRows);
            MultiDomainCCSRBuilder g(20,20,false,MultiDomainCCSRBuilder::IndependentRows);
            CHECK(pool.RunHistory(step,t,g));
            for (auto& v : visits) CHECK(v == 1);
            CCSRMat outT,outG; SymCCSRMat outSym; bool symmetric = false;
            CHECK(BuildHistoryMatrices(t,g,outT,outSym,outG,symmetric,2,step,&pool));
            CHECK(symmetric && outT.row_ptr[20] == 20);
            CHECK(t.NonZeroBlocks() == 0 && g.NonZeroBlocks() == 0);
            CHECK(createCalls == 2);
        }
        {
            ResetQueueVisits();
            MultiDomainCCSRBuilder t(20,20,false,MultiDomainCCSRBuilder::IndependentRows);
            MultiDomainCCSRBuilder g(20,20,false,MultiDomainCCSRBuilder::IndependentRows);
            CHECK(pool.RunHistory(1,t,g));
            CCSRMat outT,outG; SymCCSRMat outSym; bool symmetric = false;
            exportFirstCleared = false; injectedTaskFailure = 3;
            CHECK(!BuildHistoryMatrices(t,g,outT,outSym,outG,symmetric,2,1,&pool));
            CHECK(pool.ResourceFailed() && exportFirstCleared);
            CHECK(t.NonZeroBlocks() > 0 && t.NonZeroBlocks() < 20);
            // A retry must reconstruct the lag, never export the partially cleared builders.
            t.Reset(20,20,false,MultiDomainCCSRBuilder::IndependentRows);
            g.Reset(20,20,false,MultiDomainCCSRBuilder::IndependentRows);
            injectedTaskFailure = 0; ResetQueueVisits();
            CHECK(pool.RunHistory(1,t,g));
            CHECK(BuildHistoryMatrices(t,g,outT,outSym,outG,symmetric,2,1,&pool));
            CHECK(outT.row_ptr[20] == 20 && createCalls == 2);
        }
        for (int failure : {1,2}) {
            injectedTaskFailure = failure; ResetQueueVisits();
            CHECK(!pool.RunBatch(pool.SourceTasks()));
            CHECK(pool.ResourceFailed() == (failure == 1));
            for (auto& active : activeSlots) CHECK(active == 0);
            MultiDomainCCSRBuilder unknown(20,20),known(20,20);
            bool rejected = false;
            try { pool.MergeInto(unknown,known); } catch (const std::runtime_error&) { rejected = true; }
            CHECK(rejected && unknown.NonZeroBlocks() == 0);
        }
        injectedTaskFailure = 0; ResetQueueVisits();
        CHECK(pool.RunBatch(pool.SourceTasks()));
        for (auto& v : visits) CHECK(v == 1);
        CHECK(createCalls == 2);
    }
    requireOverlap = false;
    createCalls = 0;
    {
        MultiDomainAssemblyPool pool(elements,model,dofs,cache,nullptr,nullptr,1);
        CHECK(pool.Initialize()); ResetQueueVisits();
        CHECK(pool.RunBatch(pool.SourceTasks()));
        for (auto& v : visits) CHECK(v == 1);
        CHECK(createCalls == 0);
        injectedTaskFailure = 2; ResetQueueVisits();
        CHECK(!pool.RunBatch(pool.SourceTasks()));
        CHECK(!pool.ResourceFailed() && finishedTasks == 1);
        for (auto& v : visits) CHECK(v == 0); // Nothing after the failed first task was claimed.
        injectedTaskFailure = 0;
    }
    g_assemblyObserveTask = nullptr; g_assemblyTaskFinished = nullptr;
    g_historyExecuteRows = MultiDomainHistoryThreadMain;
    CHECK(pthread_cond_destroy(&queueTestReady) == 0);
    CHECK(pthread_mutex_destroy(&queueTestMutex) == 0);
    printf("PASS: cross-domain queue, fixed slots, material restoration, history/export reuse, failed-batch isolation\n");
}
