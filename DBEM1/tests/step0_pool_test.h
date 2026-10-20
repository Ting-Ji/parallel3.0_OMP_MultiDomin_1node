struct TestStep0Scratch {
    double coef[2][8][9] = {}, u[2][8][9] = {}, t[2][8][9] = {};
    double* cp[2][8]; double* up[2][8]; double* tp[2][8];
    double** cpp[2]; double** upp[2]; double** tpp[2];
    double*** savedC = AssistCoef; double*** savedU = AssistUij; double*** savedT = AssistTij;
    long savedThreads = thread_num;
    int savedMethod = DSquareElement::FlagDynaMethod, savedProblem = DSquareElement::Flagproblem;
    double** savedBases = DSquareElement::m_transmat;
    double basis[9] = {1,0,0,0,1,0,0,0,1};
    double* bases[2] = {basis,basis};
    TestStep0Scratch() {
        for (int i=0;i<2;++i) {
            for (int j=0;j<8;++j) { cp[i][j]=coef[i][j]; up[i][j]=u[i][j]; tp[i][j]=t[i][j]; }
            cpp[i]=cp[i]; upp[i]=up[i]; tpp[i]=tp[i];
        }
        AssistCoef=cpp; AssistUij=upp; AssistTij=tpp; thread_num=2;
        DSquareElement::Flagproblem=0; DSquareElement::FlagDynaMethod=1;
        DSquareElement::m_transmat=bases;
    }
    ~TestStep0Scratch() {
        AssistCoef=savedC; AssistUij=savedU; AssistTij=savedT; thread_num=savedThreads;
        DSquareElement::FlagDynaMethod=savedMethod; DSquareElement::Flagproblem=savedProblem;
        DSquareElement::m_transmat=savedBases;
        g_step0ObserveContext=nullptr; g_step0ExecuteRows=&MultiDomainStep0ThreadMain;
        g_assemblyCreateThread=&pthread_create;
    }
};
static int expectedStep0Domain[2], observedStep0Batches[2];
static void ObserveStep0Context(const Step0WorkerContext& c)
{
    CHECK(&DSquareElement::ActiveDynaMat() == &c.domain->mat);
    CHECK(c.coefSlot == AssistCoef[c.threadId] && c.uSlot == AssistUij[c.threadId] && c.tSlot == AssistTij[c.threadId]);
    CHECK(c.coefSlot != AssistCoef[1-c.threadId]);
    ++observedStep0Batches[c.threadId];
}
static void* ObserveStep0Restore(void* data)
{
    const DynaMat* previous = &DSquareElement::ActiveDynaMat();
    MultiDomainStep0ThreadMain(data);
    auto& c = *static_cast<Step0WorkerContext*>(data);
    if (&DSquareElement::ActiveDynaMat() != previous) c.ok = -1;
    return nullptr;
}
static int injectedStep0Status;
static void* FailStep0AfterRows(void* data)
{
    MultiDomainStep0ThreadMain(data);
    static_cast<Step0WorkerContext*>(data)->ok = injectedStep0Status;
    return nullptr;
}
static void TestStep0Pool()
{
    TestStep0Scratch scratch;
    DSquareElement elements[2];
    MultiDomainModel model;
    model.elementCount=model.nodeCount=2;
    for (int i=0;i<2;++i) {
        elements[i].BCID=456; elements[i].InfEleFlag=-1; elements[i].m_nodeID[0]=i;
        for (int k=0;k<9;++k) {
            elements[i].m_OnElementU1ij[0][0][k]=k+1;
            elements[i].m_OnElementT1ij[0][0][k]=2*(k+1);
        }
        Domain domain = {};
        domain.id=i; domain.eleBegin=i; domain.eleCount=1; domain.elementIds={i};
        domain.mat.C1=10+i; model.domains.push_back(domain);
    }
    // Explicit override selects the right kernel without changing the element.
    for (int code : {123,456}) {
        elements[0].BCID=456;
        CHECK(UnknownSubMatrix(elements,0,0,nullptr,nullptr,0L,code));
        CHECK(elements[0].BCID==456);
        double actual[9]; CopyBlock9(AssistCoef[0][0],actual);
        for (int k=0;k<9;++k) CHECK(actual[k] == (code==123 ? -(k+1.0) : 2*(k+1.0)));
        elements[0].BCID=code; // Test-only old API reference.
        CHECK(UnknownSubMatrix(elements,0,0,nullptr,nullptr,0L));
        for (int k=0;k<9;++k) CHECK(actual[k]==AssistCoef[0][0][k]);
        elements[0].BCID=456;
        CHECK(knownSubMatrix(elements,0,0,nullptr,nullptr,0L,code));
        CHECK(elements[0].BCID==456);
        CopyBlock9(AssistCoef[0][0],actual);
        elements[0].BCID=code;
        CHECK(knownSubMatrix(elements,0,0,nullptr,nullptr,0L));
        for (int k=0;k<9;++k) CHECK(actual[k]==AssistCoef[0][0][k]);
    }
    elements[0].BCID=456;
    GlobalDofMap dofs; CHECK(dofs.Build(model,elements,2));
    const MultiDomainInterfaceCache cache(model,elements);
    g_assemblyCreateThread=InjectCreate; createCalls=0; failAt=-1;
    g_step0ObserveContext=ObserveStep0Context;
    g_step0ExecuteRows=ObserveStep0Restore;
    observedStep0Batches[0]=observedStep0Batches[1]=0;
    {
        MultiDomainAssemblyPool pool(elements,model,dofs,cache,nullptr,nullptr,2);
        CHECK(pool.Initialize()); CHECK(createCalls==2);
        expectedStep0Domain[0]=0; expectedStep0Domain[1]=1;
        CHECK(pool.RunBatch({{0,0,1},{1,0,1}}));
        expectedStep0Domain[0]=1; expectedStep0Domain[1]=0;
        CHECK(pool.RunBatch({{1,0,1},{0,0,1}}));
        CHECK(createCalls==2);
        CHECK(observedStep0Batches[0]+observedStep0Batches[1]==4);
        CHECK(elements[0].BCID==456 && elements[1].BCID==456);
    }
    g_step0ObserveContext=nullptr; g_step0ExecuteRows=&MultiDomainStep0ThreadMain;
    createCalls=0; failAt=1;
    {
        MultiDomainAssemblyPool pool(elements,model,dofs,cache,nullptr,nullptr,2);
        CHECK(!pool.Initialize()); CHECK(createCalls==2);
    }
    createCalls=0; failAt=1;
    MultiDomainCCSRBuilder unknown,known;
    CHECK(AssembleMultiDomainStep0Pthread(cache,elements,model,dofs,nullptr,nullptr,unknown,known,2)
        == (DBEMPolicy::allowSerialFallback ? 1 : 0));
    failAt=-1;
    g_step0ExecuteRows=FailStep0AfterRows;
    for (int status : {0,-1}) {
        injectedStep0Status=status;
        const bool ok=AssembleMultiDomainStep0Pthread(cache,elements,model,dofs,nullptr,nullptr,unknown,known,2)!=0;
        CHECK(ok == (status==0 && DBEMPolicy::allowSerialFallback));
        if (ok) {
            MultiDomainCCSRBuilder refUnknown,refKnown;
            CHECK(AssembleMultiDomainStep0(cache,elements,model,dofs,nullptr,nullptr,refUnknown,refKnown));
            CHECK(SparseBuilderDifference(unknown,refUnknown)==0);
            CHECK(SparseBuilderDifference(known,refKnown)==0);
        } else CHECK(unknown.NonZeroBlocks()==0 && known.NonZeroBlocks()==0);
    }
}
