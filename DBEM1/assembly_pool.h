#pragma once
// Included after the integral row kernels. Only a worker ID selects scratch.
enum class AssemblyTaskKind { Step0, History, Output };
struct AssemblyTask {
    int domainId = -1;
    long sourceBegin = 0, sourceEnd = 0;
    AssemblyTaskKind kind = AssemblyTaskKind::Step0;
    void* output = nullptr;
    void* (*exportRows)(void*) = nullptr;
    int phase = -1;
};
struct AssemblyBatch {
    long step = 0;
    MultiDomainCCSRBuilder* mts = nullptr;
    MultiDomainCCSRBuilder* mgs = nullptr;
};
static auto g_assemblyCreateThread = &pthread_create;
static auto g_step0ExecuteRows = &MultiDomainStep0ThreadMain;
static auto g_historyExecuteRows = &MultiDomainHistoryThreadMain;
#ifdef DBEM_STEP0_TEST_HOOKS
static void (*g_assemblyObserveTask)(const AssemblyTask&, long) = nullptr;
static void (*g_assemblyTaskFinished)(const AssemblyTask&, long, bool) = nullptr;
#endif

class MultiDomainAssemblyPool {
    struct Worker { MultiDomainAssemblyPool* pool; size_t index; };
    const MultiDomainModel* model = nullptr;
    const GlobalDofMap* dofMap = nullptr;
    std::vector<Step0WorkerContext> contexts;
    std::vector<Worker> workers;
    std::vector<pthread_t> threads;
    std::vector<AssemblyTask> sourceTasks;
    std::vector<AssemblyTask> historyTasks;
    struct FailureRecord { AssemblyTask task; long step = 0; int status = 1; };
    std::vector<FailureRecord> errors;
    const std::vector<AssemblyTask>* tasks = nullptr;
    AssemblyBatch batch;
    pthread_mutex_t mutex;
    pthread_cond_t ready, done;
    bool mutexReady = false, readyInitialized = false, doneInitialized = false;
    bool stopping = false, initialized = false, batchSucceeded = false;
    size_t created = 0, generation = 0, completed = 0, next = 0;
    int failure = 0; // 0 success, 1 resource failure, 2 invalid data/structure.

    int Execute(const AssemblyTask& task, size_t slot) noexcept
    {
        int status = -1;
#ifdef DBEM_STEP0_TEST_HOOKS
        const DynaMat* previousMaterial = &DSquareElement::ActiveDynaMat();
#endif
        try {
            if (task.kind != AssemblyTaskKind::Output &&
                (!model || task.domainId < 0 || (size_t)task.domainId >= model->domains.size()))
                throw std::runtime_error("invalid assembly task domain");
            if (task.kind == AssemblyTaskKind::Output) {
                if (!task.exportRows || !task.output) throw std::runtime_error("invalid export task");
#ifdef DBEM_STEP0_TEST_HOOKS
                if (g_assemblyObserveTask) g_assemblyObserveTask(task, (long)slot);
#endif
                task.exportRows(task.output);
                status = 1;
            } else {
                const auto& domain = model->domains[(size_t)task.domainId];
                if (task.sourceBegin < 0 || task.sourceEnd <= task.sourceBegin ||
                    (size_t)task.sourceEnd > domain.elementIds.size())
                    throw std::runtime_error("invalid assembly source range");
                ScopedDynaMat material(domain.mat);
#ifdef DBEM_STEP0_TEST_HOOKS
                if (g_assemblyObserveTask) g_assemblyObserveTask(task, (long)slot);
#endif
                auto& c = contexts[slot];
                if (task.kind == AssemblyTaskKind::Step0) {
                    c.domain = &domain;
                    c.sourceBegin = task.sourceBegin; c.sourceEnd = task.sourceEnd;
                    c.ok = 1;
                    g_step0ExecuteRows(&c);
                    status = c.ok;
                } else {
                    if (batch.step <= 0 || !batch.mts || !batch.mgs)
                        throw std::runtime_error("invalid history batch");
                    MultiDomainHistoryThreadArg h = {};
                    h.elements = c.elements; h.domain = &domain; h.dofMap = dofMap;
                    h.step = batch.step; h.mtsBuilder = batch.mts; h.mgsBuilder = batch.mgs;
                    h.sourceBegin = task.sourceBegin; h.sourceEnd = task.sourceEnd;
                    h.threadId = (long)slot;
                    g_historyExecuteRows(&h);
                    status = h.ok;
                }
            }
        } catch (const std::bad_alloc&) { status = 0; }
        catch (const std::exception& e) {
            fprintf(stderr, "Assembly task exception: %s\n", e.what());
        } catch (...) { fprintf(stderr, "Assembly task unknown exception\n"); }
        if (status != 1) {
            errors[slot] = {task, batch.step, status};
            fprintf(stderr, "Assembly failed kind=%d phase=%d step=%ld domain=%d range=[%ld,%ld) worker=%zu category=%s\n",
                (int)task.kind, task.phase, batch.step, task.domainId, task.sourceBegin, task.sourceEnd,
                slot, status == 0 ? "resource" : "data/structure");
        }
#ifdef DBEM_STEP0_TEST_HOOKS
        if (g_assemblyTaskFinished) {
            try { g_assemblyTaskFinished(task, (long)slot, previousMaterial == &DSquareElement::ActiveDynaMat()); }
            catch (...) { status = -1; }
        }
#endif
        return status;
    }
    static void* Work(void* data)
    {
        auto& worker = *static_cast<Worker*>(data);
        auto& pool = *worker.pool;
        size_t seen = 0;
        RequireThreadSuccess(pthread_mutex_lock(&pool.mutex), "assembly lock");
        for (;;) {
            while (!pool.stopping && seen == pool.generation)
                RequireThreadSuccess(pthread_cond_wait(&pool.ready, &pool.mutex), "assembly wait");
            if (pool.stopping) break;
            seen = pool.generation;
            while (!pool.failure && pool.next < pool.tasks->size()) {
                const size_t index = pool.next++;
                RequireThreadSuccess(pthread_mutex_unlock(&pool.mutex), "assembly unlock");
                const int status = pool.Execute((*pool.tasks)[index], worker.index);
                RequireThreadSuccess(pthread_mutex_lock(&pool.mutex), "assembly lock");
                if (status != 1) pool.failure = (std::max)(pool.failure, status == 0 ? 1 : 2);
            }
            if (++pool.completed == pool.created)
                RequireThreadSuccess(pthread_cond_signal(&pool.done), "assembly done signal");
        }
        RequireThreadSuccess(pthread_mutex_unlock(&pool.mutex), "assembly unlock");
        return nullptr;
    }
public:
    // Output-only compatibility callers need no integral workspace.
    explicit MultiDomainAssemblyPool(long count)
    {
        if (count < 1) throw std::runtime_error("invalid assembly worker count");
        contexts.resize((size_t)count); workers.resize((size_t)count); threads.resize((size_t)count);
        errors.resize((size_t)count);
        for (size_t i = 0; i < workers.size(); ++i) workers[i] = {this, i};
    }
    MultiDomainAssemblyPool(DSquareElement* elements, const MultiDomainModel& m,
        const GlobalDofMap& dofs, const MultiDomainInterfaceCache& cache,
        long** infElePid, long** elePid, long count, bool rowsValidated = false) : MultiDomainAssemblyPool(count)
    {
        if (!elements || count > thread_num || !AssistCoef || !AssistUij || !AssistTij)
            throw std::runtime_error("invalid assembly integral workspace");
        if (DSquareElement::Flagproblem != 0)
            throw std::runtime_error("parallel assembly requires dynamic finite elements");
        for (long e = 0; e < m.elementCount; ++e)
            if (elements[e].InfEleFlag != -1)
                throw std::runtime_error("parallel assembly does not support infinite elements");
        if (!rowsValidated && !HistoryRowsHaveUniqueOwners(elements, m, dofs))
            throw std::runtime_error("assembly rows must have unique owners");
        model = &m; dofMap = &dofs;
        for (long i = 0; i < count; ++i) {
            auto& c = contexts[(size_t)i];
            c.elements = elements; c.dofMap = &dofs; c.interfaceCache = &cache;
            c.infElePid = infElePid; c.elePid = elePid; c.threadId = i;
            c.coefSlot = AssistCoef[i]; c.uSlot = AssistUij[i]; c.tSlot = AssistTij[i];
            if (!c.coefSlot || !c.uSlot || !c.tSlot)
                throw std::runtime_error("missing assembly integral workspace slot");
            c.unknownBuilder = &c.unknownStorage; c.knownBuilder = &c.knownStorage;
        }
        std::vector<long> begin(m.domains.size(), 0), chunk(m.domains.size());
        for (size_t d = 0; d < m.domains.size(); ++d) {
            const size_t n = m.domains[d].elementIds.size(), parts = 4 * (size_t)count;
            chunk[d] = (long)(std::max)((size_t)1, n / parts + (n % parts != 0));
        }
        bool more = true;
        while (more) {
            more = false;
            for (size_t d = 0; d < m.domains.size(); ++d) {
                const long n = (long)m.domains[d].elementIds.size();
                if (begin[d] == n) continue;
                const long end = begin[d] + (std::min)(chunk[d], n - begin[d]);
                sourceTasks.push_back({(int)d, begin[d], end});
                begin[d] = end; more = true;
            }
        }
        historyTasks = sourceTasks;
        for (auto& task : historyTasks) task.kind = AssemblyTaskKind::History;
    }
    MultiDomainAssemblyPool(const MultiDomainAssemblyPool&) = delete;
    MultiDomainAssemblyPool& operator=(const MultiDomainAssemblyPool&) = delete;
    ~MultiDomainAssemblyPool() { Close(); }
    size_t WorkerCount() const { return workers.size(); }
    bool ResourceFailed() const { return failure == 1; }
    // Coordinator-only: no tasks may be active when recording allocation failure.
    void RecordResourceFailure() { failure = 1; batchSucceeded = false; }
    const std::vector<AssemblyTask>& SourceTasks() const { return sourceTasks; }
    bool Initialize()
    {
        if (initialized) return true;
        stopping = false; generation = completed = 0;
        if (workers.size() == 1) { initialized = true; return true; }
        int error = pthread_mutex_init(&mutex, nullptr);
        if (!error) { mutexReady = true; error = pthread_cond_init(&ready, nullptr); }
        if (!error) { readyInitialized = true; error = pthread_cond_init(&done, nullptr); }
        if (!error) {
            doneInitialized = true;
            for (; created < threads.size(); ++created) {
                error = g_assemblyCreateThread(&threads[created], nullptr, Work, &workers[created]);
                if (error) break;
            }
        }
        if (error) {
            fprintf(stderr, "Assembly pool initialization failed: code=%d created=%zu\n", error, created);
            failure = 1; Close(); return false;
        }
        initialized = true; return true;
    }
    bool RunBatch(const std::vector<AssemblyTask>& work, const AssemblyBatch& input = {})
    {
        batchSucceeded = false;
        if (!initialized) throw std::runtime_error("assembly pool is not initialized");
        failure = 0; next = 0; tasks = &work; batch = input;
        for (auto& error : errors) error.status = 1;
        try {
            if (!work.empty() && work.front().kind == AssemblyTaskKind::Step0) {
                if (!dofMap) throw std::runtime_error("Step0 batch requires integral context");
                for (auto& c : contexts) {
                    c.unknownStorage.Reset(dofMap->equationBlockCount, dofMap->unknownBlockCount);
                    c.knownStorage.Reset(dofMap->equationBlockCount, dofMap->knownBlockCount);
                }
            }
        } catch (const std::bad_alloc&) {
            failure = 1; tasks = nullptr; DiscardStep0(); return false;
        }
        if (workers.size() == 1) {
            while (!failure && next < work.size()) {
                const int status = Execute(work[next++], 0);
                if (status != 1) failure = status == 0 ? 1 : 2;
            }
        } else {
            RequireThreadSuccess(pthread_mutex_lock(&mutex), "assembly batch lock");
            completed = 0; ++generation;
            RequireThreadSuccess(pthread_cond_broadcast(&ready), "assembly batch broadcast");
            while (completed != created)
                RequireThreadSuccess(pthread_cond_wait(&done, &mutex), "assembly batch wait");
            RequireThreadSuccess(pthread_mutex_unlock(&mutex), "assembly batch unlock");
        }
        tasks = nullptr;
        batchSucceeded = failure == 0;
        if (!batchSucceeded && !work.empty() && work.front().kind == AssemblyTaskKind::Step0)
            DiscardStep0();
        return batchSucceeded;
    }
    bool RunHistory(long step, MultiDomainCCSRBuilder& mts, MultiDomainCCSRBuilder& mgs)
    {
        return RunBatch(historyTasks, {step, &mts, &mgs});
    }
    void DiscardStep0()
    {
        for (auto& c : contexts) {
            c.unknownStorage.Reset(0, 0); c.knownStorage.Reset(0, 0);
        }
        batchSucceeded = false;
    }
    void MergeInto(MultiDomainCCSRBuilder& unknown, MultiDomainCCSRBuilder& known)
    {
        if (!batchSucceeded) throw std::runtime_error("cannot publish failed assembly batch");
        try {
            for (const auto& c : contexts) {
                unknown.MergeFrom(c.unknownStorage); known.MergeFrom(c.knownStorage);
            }
        } catch (...) { DiscardStep0(); throw; }
        // The pool lives through all history lags; do not retain duplicate Step0 matrices.
        DiscardStep0();
    }
    void Close()
    {
        if (created) {
            RequireThreadSuccess(pthread_mutex_lock(&mutex), "assembly close lock");
            stopping = true;
            RequireThreadSuccess(pthread_cond_broadcast(&ready), "assembly close broadcast");
            RequireThreadSuccess(pthread_mutex_unlock(&mutex), "assembly close unlock");
            for (size_t i = 0; i < created; ++i)
                RequireThreadSuccess(pthread_join(threads[i], nullptr), "assembly join");
        }
        created = 0;
        if (doneInitialized) RequireThreadSuccess(pthread_cond_destroy(&done), "assembly destroy done");
        if (readyInitialized) RequireThreadSuccess(pthread_cond_destroy(&ready), "assembly destroy ready");
        if (mutexReady) RequireThreadSuccess(pthread_mutex_destroy(&mutex), "assembly destroy mutex");
        initialized = mutexReady = readyInitialized = doneInitialized = false;
    }
};
