// Test-only access to the real RHS implementation, in the same translation unit.
static void RhsTestRequire(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
static void RhsTestEqual(const Wvector& a, const Wvector& b, const char* message)
{
    RhsTestRequire(a.n == b.n && !memcmp(a.b, b.b, sizeof(double) * a.n), message);
}
template<class Matrix>
static void RhsTestMatrix(Matrix& matrix, long rows, long columns, int width, int seed)
{
    matrix.initial(rows, rows * columns);
    int at = 0; matrix.row_ptr[0] = 0;
    for (long row = 0; row < rows; ++row)
    {
        const long entries = row == rows / 2 ? 0 : row == 0 ? columns : (std::min)(columns, 2L);
        for (long j = 0; j < entries; ++j)
        {
            matrix.col_index[at] = (j + row) % columns;
            for (int c = 0; c < width; ++c)
                matrix.value[(size_t)at * width + c] = sin((at * width + c + seed) * .71) * .13;
            ++at;
        }
        matrix.row_ptr[row + 1] = at;
    }
}
static int g_rhsTestCreateCalls;
static int RhsTestFailCreate(pthread_t* thread, const pthread_attr_t* attr, void* (*entry)(void*), void* arg)
{
    if (g_rhsTestCreateCalls++ == 2) return EAGAIN;
    return pthread_create(thread, attr, entry, arg);
}

void RhsCpuRunTests()
{
    thread_num = 40;
    const long rows = 53, knownBlocks = 7, limit = 3, steps = 7;
    CCSRMat knownMatrix, mts[limit + 1], full[limit + 1];
    SymCCSRMat sym[limit + 1];
    int useSym[limit + 1] = { 0, 1, 0, 1 };
    RhsTestMatrix(knownMatrix, rows, knownBlocks, 9, 13);
    for (long lag = 1; lag <= limit; ++lag)
    {
        RhsTestMatrix(mts[lag], rows, rows, 9, (int)lag);
        RhsTestMatrix(full[lag], rows, rows, 9, (int)lag + 4);
        RhsTestMatrix(sym[lag], rows, rows, 6, (int)lag + 9);
    }
    std::vector<MultiDomainState> states(steps + 1);
    for (long step = 0; step <= steps; ++step)
    {
        states[step].globalT.initial(3 * rows); states[step].globalU.initial(3 * rows);
        for (long i = 0; i < 3 * rows; ++i)
        {
            states[step].globalT.b[i] = cos(i * .13 + step);
            states[step].globalU.b[i] = sin(i * .17 - step);
        }
    }
    Wvector known(3 * knownBlocks), boundary(21), rhs(3 * rows), temp(3 * rows),
        convU(3 * rows), convT(3 * rows), x0(3 * rows), x(3 * rows),
        dg(3 * rows), dt(3 * rows), dk(3 * rows), expected(3 * rows),
        referenceKnown(3 * knownBlocks), referenceG(3 * rows), referenceT(3 * rows), referenceK(3 * rows);
    for (long i = 0; i < boundary.n; ++i) boundary.b[i] = sin(i * .31);
    std::vector<std::vector<MultiDomainHistorySpMVTask<CCSRMat> > > tt(limit + 1), gf(limit + 1);
    std::vector<std::vector<MultiDomainHistorySpMVTask<SymCCSRMat> > > gs(limit + 1);
    MultiDomainHistoryThreadPool pool;
    RhsTestRequire(pool.Initialize(40), "RHS pool initialize");
    _putenv_s("DBEM_RHS_CPU_PROFILE", "1");
    RhsCpuWorkspace w(pool.WorkerCount());
    w.knownMap.assign(knownBlocks, RhsKnownSource());
    referenceKnown = 0.0;
    for (long b = 0; b < knownBlocks - 1; ++b)
    {
        RhsCpuWorkspace::StoreKnown(w.knownMap, b, 3 * b, b % 2 ? -1.0 : 1.0, boundary.n);
        for (int c = 0; c < 3; ++c) referenceKnown.b[3 * b + c] = (b % 2 ? -1.0 : 1.0) * boundary.b[3 * b + c];
    }
    RhsCpuWorkspace::StoreKnown(w.knownMap, 1, 12, -.5, boundary.n);
    for (int c = 0; c < 3; ++c) referenceKnown.b[3 + c] = -.5 * boundary.b[12 + c];
    bool caught = false;
    try { RhsCpuWorkspace::StoreKnown(w.knownMap, knownBlocks, 0, 1, boundary.n); }
    catch (const std::runtime_error&) { caught = true; }
    RhsTestRequire(caught, "RHS map bounds");
    w.knownMatrix = &knownMatrix; RhsCpuWorkspace::Partition(knownMatrix, w.count, w.knownRows);
    w.known = &known; w.boundary = &boundary; w.rhs = &rhs; w.temp = &temp;
    w.convU = &convU; w.convT = &convT; w.x0 = &x0; w.x = &x;
    w.states = &states; w.mts = mts; w.mgsFull = full; w.mgs = sym; w.symmetric = useSym;
    w.tTasks = &tt; w.gTasks = &gs; w.gFullTasks = &gf;
    w.diagKnown = &dk; w.diagG = &dg; w.diagT = &dt;
    w.maxLag = limit;
    for (long lag = 1; lag <= limit; ++lag)
    {
        PrepareRhsSpMVTasks(mts[lag], convU, temp, tt[lag], w.count);
        PrepareRhsSpMVTasks(full[lag], convT, temp, gf[lag], w.count);
        PrepareRhsSpMVTasks(sym[lag], convT, temp, gs[lag], w.count);
    }
    for (int diagnostics = 0; diagnostics < 2; ++diagnostics)
    for (long step = 1; step <= steps; ++step)
    {
        SparseMul(knownMatrix, referenceKnown, expected); referenceK = expected;
        referenceG = 0.0; referenceT = 0.0;
        for (long lag = 1; step > 1 && lag <= step && lag <= limit; ++lag)
        {
            BuildConvolutionVector(states, step, lag, 1, convT);
            if (useSym[lag]) SparseMul(sym[lag], convT, temp); else SparseMul(full[lag], convT, temp);
            AddVectorInPlace(expected, temp, 1.0); AddVectorInPlace(referenceG, temp, 1.0);
            BuildConvolutionVector(states, step, lag, 0, convU); SparseMul(mts[lag], convU, temp);
            AddVectorInPlace(expected, temp, -1.0); AddVectorInPlace(referenceT, temp, -1.0);
        }
        w.diagnostics = diagnostics != 0; w.step = step; x0 = 5.0; x = -3.0;
        pool.Submit(&w, RhsCpuWorkspace::Execute); w.Finish();
        RhsTestEqual(known, referenceKnown, "RHS known gather/sign/last assignment/unmapped");
        RhsTestEqual(rhs, expected, "RHS exact full task");
        for (long i = 0; i < x.n; ++i) RhsTestRequire(x.b[i] == 0 && x0.b[i] == 0, "RHS initial guess zero");
        if (diagnostics)
        {
            RhsTestEqual(dk, referenceK, "RHS known diagnostic");
            RhsTestEqual(dg, referenceG, "RHS G diagnostic");
            RhsTestEqual(dt, referenceT, "RHS T diagnostic");
        }
        for (long lag = 1; lag <= limit; ++lag) for (int traction = 0; traction < 2; ++traction)
        {
            Wvector reference(3 * rows);
            BuildConvolutionVector(states, step, lag, traction, reference);
            for (size_t id = 0; id < w.count; ++id) w.BuildConvolution(convU, traction != 0, lag, id);
            RhsTestEqual(convU, reference, "RHS exact convolution");
        }
    }
    boundary.b[12] = std::numeric_limits<double>::quiet_NaN(); w.step = 1;
    pool.Submit(&w, RhsCpuWorkspace::Execute);
    caught = false; try { w.Finish(); } catch (const std::runtime_error&) { caught = true; }
    RhsTestRequire(caught, "RHS nonfinite rejected after synchronized completion");
    boundary.b[12] = 1.0; pool.Submit(&w, RhsCpuWorkspace::Execute); w.Finish();
    pool.Close(); pool.Close();
    auto original = g_historyCreateThread; g_rhsTestCreateCalls = 0; g_historyCreateThread = RhsTestFailCreate;
    const bool initialized = pool.Initialize(40);
    g_historyCreateThread = original;
    RhsTestRequire(initialized == DBEMPolicy::allowSerialFallback && pool.WorkerCount() == 1, "RHS partial thread failure cleanup");
    RhsTestRequire(pool.Initialize(40), "RHS initialize after failure"); pool.Close();
    _putenv_s("DBEM_RHS_CPU_PROFILE", "0");
    _putenv_s("DBEM_RHS_CPU_MODE", "parallel"); RhsTestRequire(RhsParallelSelected(), "RHS mode parallel");
    _putenv_s("DBEM_RHS_CPU_MODE", "legacy"); RhsTestRequire(!RhsParallelSelected(), "RHS mode legacy");
    _putenv_s("DBEM_RHS_CPU_MODE", "invalid"); caught = false;
    try { RhsParallelSelected(); } catch (const std::runtime_error&) { caught = true; }
    RhsTestRequire(caught, "RHS invalid mode"); _putenv_s("DBEM_RHS_CPU_MODE", "");
    RhsTestRequire(!RhsParallelSelected(), "RHS default legacy");
}
