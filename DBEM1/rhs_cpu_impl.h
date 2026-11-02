// Included only by multidomain.cpp, after its legacy RHS operators.
// The solve session owns this workspace; no matrix-address cache is retained.
static bool RhsParallelSelected()
{
    const char* mode = std::getenv("DBEM_RHS_CPU_MODE");
    if (!mode || !*mode || !strcmp(mode, "legacy")) return false;
    if (!strcmp(mode, "parallel")) return true;
    throw std::runtime_error("DBEM_RHS_CPU_MODE must be legacy or parallel");
}

class RhsTeamBarrier
{
    pthread_mutex_t mutex;
    pthread_cond_t ready;
    size_t arrived, generation, count;
    RhsTeamBarrier(const RhsTeamBarrier&);
    RhsTeamBarrier& operator=(const RhsTeamBarrier&);
public:
    explicit RhsTeamBarrier(size_t n) : arrived(0), generation(0), count(n)
    {
        RequireThreadSuccess(pthread_mutex_init(&mutex, NULL), "RHS barrier mutex init");
        RequireThreadSuccess(pthread_cond_init(&ready, NULL), "RHS barrier cond init");
    }
    ~RhsTeamBarrier()
    {
        RequireThreadSuccess(pthread_cond_destroy(&ready), "RHS barrier cond destroy");
        RequireThreadSuccess(pthread_mutex_destroy(&mutex), "RHS barrier mutex destroy");
    }
    void Wait()
    {
        if (count <= 1) return;
        RequireThreadSuccess(pthread_mutex_lock(&mutex), "RHS barrier lock");
        const size_t seen = generation;
        if (++arrived == count)
        {
            arrived = 0; ++generation;
            RequireThreadSuccess(pthread_cond_broadcast(&ready), "RHS barrier broadcast");
        }
        else while (seen == generation)
            RequireThreadSuccess(pthread_cond_wait(&ready, &mutex), "RHS barrier wait");
        RequireThreadSuccess(pthread_mutex_unlock(&mutex), "RHS barrier unlock");
    }
};

struct RhsKnownSource
{
    long scalar;
    double sign;
    RhsKnownSource(long s = -1, double v = 1.0) : scalar(s), sign(v) {}
};

static void RhsProductRow(CCSRMat& M, Wvector& x, long i, double* values)
{
	int j;
	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

		values[0] = 0.0;
		values[1] = 0.0;
		values[2] = 0.0;
		rownum = M.row_ptr[i + 1] - M.row_ptr[i];
		rowstart = M.row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			btemp[0] = x[M.col_index[rowstart + j] * 3 + 1];
			btemp[1] = x[M.col_index[rowstart + j] * 3 + 2];
			btemp[2] = x[M.col_index[rowstart + j] * 3 + 3];
			values[0] += M.value[rowstart * 9 + j * 9] * btemp[0] + M.value[rowstart * 9 + j * 9 + 1] * btemp[1] + M.value[rowstart * 9 + j * 9 + 2] * btemp[2];
			values[1] += M.value[rowstart * 9 + j * 9 + 3] * btemp[0] + M.value[rowstart * 9 + j * 9 + 4] * btemp[1] + M.value[rowstart * 9 + j * 9 + 5] * btemp[2];
			values[2] += M.value[rowstart * 9 + j * 9 + 6] * btemp[0] + M.value[rowstart * 9 + j * 9 + 7] * btemp[1] + M.value[rowstart * 9 + j * 9 + 8] * btemp[2];
		}


}

static void RhsProductRow(SymCCSRMat& M, Wvector& x, long i, double* values)
{
	int j;
	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

		values[0] = 0.0;
		values[1] = 0.0;
		values[2] = 0.0;
		rownum = M.row_ptr[i + 1] - M.row_ptr[i];
		rowstart = M.row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			btemp[0] = x[M.col_index[rowstart + j] * 3 + 1];
			btemp[1] = x[M.col_index[rowstart + j] * 3 + 2];
			btemp[2] = x[M.col_index[rowstart + j] * 3 + 3];
			values[0] += M.value[rowstart * 6 + j * 6] * btemp[0] + M.value[rowstart * 6 + j * 6 + 1] * btemp[1] + M.value[rowstart * 6 + j * 6 + 2] * btemp[2];
			values[1] += M.value[rowstart * 6 + j * 6 + 1] * btemp[0] + M.value[rowstart * 6 + j * 6 + 3] * btemp[1] + M.value[rowstart * 6 + j * 6 + 4] * btemp[2];
			values[2] += M.value[rowstart * 6 + j * 6 + 2] * btemp[0] + M.value[rowstart * 6 + j * 6 + 4] * btemp[1] + M.value[rowstart * 6 + j * 6 + 5] * btemp[2];
		}


}

struct RhsCpuWorkspace
{
    enum Stage { KnownGather, KnownMultiply, Convolution, GMultiplyAdd,
        TMultiplyAdd, CheckZero, Wait, StageCount };
    struct WorkerStats
    {
        double seconds[StageCount];
        double lastTaskWall;
        long firstInvalid;
        WorkerStats() : lastTaskWall(0.0), firstInvalid(-1) { std::fill(seconds, seconds + StageCount, 0.0); }
    };
    const size_t count;
    const bool profile;
    RhsTeamBarrier barrier;
    std::vector<RhsKnownSource> knownMap;
    std::vector<PthArg> knownRows;
    std::vector<WorkerStats> stats;
    CCSRMat *knownMatrix, *mts, *mgsFull;
    SymCCSRMat* mgs;
    const int* symmetric;
    const std::vector<MultiDomainState>* states;
    const std::vector<std::vector<MultiDomainHistorySpMVTask<CCSRMat> > > *tTasks, *gFullTasks;
    const std::vector<std::vector<MultiDomainHistorySpMVTask<SymCCSRMat> > >* gTasks;
    Wvector *known, *boundary, *rhs, *convU, *convT, *temp, *x0, *x;
    Wvector *diagKnown, *diagG, *diagT;
    long step, maxLag;
    bool diagnostics;
    double setupSeconds;

    explicit RhsCpuWorkspace(size_t threads)
        : count(threads), profile(ProfileEnabled()), barrier(threads), stats(threads),
          knownMatrix(NULL), mts(NULL), mgsFull(NULL), mgs(NULL), symmetric(NULL), states(NULL),
          tTasks(NULL), gFullTasks(NULL), gTasks(NULL), known(NULL), boundary(NULL), rhs(NULL),
          convU(NULL), convT(NULL), temp(NULL), x0(NULL), x(NULL), diagKnown(NULL), diagG(NULL),
          diagT(NULL), step(0), maxLag(0), diagnostics(false), setupSeconds(0.0) {}
    RhsCpuWorkspace(const RhsCpuWorkspace&) = delete;
    RhsCpuWorkspace& operator=(const RhsCpuWorkspace&) = delete;

    static bool ProfileEnabled()
    {
        const char* value = std::getenv("DBEM_RHS_CPU_PROFILE");
        return value && !strcmp(value, "1");
    }
    static void StoreKnown(std::vector<RhsKnownSource>& map, long dst, long src, double sign, long sourceSize)
    {
        if (dst < 0 || dst >= (long)map.size() || src < 0 || src + 2 >= sourceSize)
            throw std::runtime_error("RHS known boundary mapping out of range");
        // Same destination may be visited repeatedly: legacy's last assignment wins.
        map[(size_t)dst] = RhsKnownSource(src, sign);
    }
    void Prepare(const GlobalDofMap& map, const MultiDomainModel& model,
        DSquareElement* elements, long boundarySize, CCSRMat& matrix)
    {
        const double start = DBEMWallTime();
        knownMap.assign((size_t)map.knownBlockCount, RhsKnownSource());
        for (size_t d = 0; d < model.domains.size(); ++d)
        {
            const Domain& domain = model.domains[d];
            for (size_t i = 0; i < domain.elementIds.size(); ++i)
            {
                const long ele = domain.elementIds[i];
                if (elements[ele].SurfaceType == SurfaceInterface) continue;
                const long node = elements[ele].m_nodeID[0];
                const VariableRef ref = elements[ele].BCID == 123
                    ? map.GetU(domain.id, node) : map.GetT(domain.id, node);
                if (ref.known) StoreKnown(knownMap, ref.globalBlock, 3 * node, ref.sign, boundarySize);
            }
        }
        Partition(matrix, count, knownRows);
        knownMatrix = &matrix;
        setupSeconds = DBEMWallTime() - start;
        printf("RHS CPU mode=parallel threads=%zu rows=%ld known_blocks=%zu known_nnz=%ld profile=%d\n",
            count, (long)matrix.n, knownMap.size(), (long)matrix.NonZeroBlocks(), profile ? 1 : 0);
    }
    template<class Matrix>
    static void Partition(Matrix& matrix, size_t threads, std::vector<PthArg>& ranges)
    {
        ranges.resize(threads);
        const unsigned long long total = (unsigned long long)matrix.row_ptr[matrix.n] + matrix.n;
        long row = 0;
        for (size_t worker = 0; worker < threads; ++worker)
        {
            ranges[worker].start = row;
            const unsigned long long target = total * (worker + 1) / threads;
            while (row < matrix.n && (unsigned long long)matrix.row_ptr[row] + row < target) ++row;
            ranges[worker].end = row;
        }
        ranges.back().end = matrix.n;
    }
    double Clock() const { return profile ? DBEMWallTime() : 0.0; }
    void EndStage(size_t id, Stage stage, double start)
    {
        if (profile) stats[id].seconds[stage] += DBEMWallTime() - start;
    }
    void Sync(size_t id)
    {
        const double start = Clock();
        barrier.Wait();
        EndStage(id, Wait, start);
    }
    void Gather(size_t id)
    {
        const long begin = (long)(knownMap.size() * id / count);
        const long end = (long)(knownMap.size() * (id + 1) / count);
        for (long b = begin; b < end; ++b)
        {
            const RhsKnownSource& source = knownMap[(size_t)b];
            for (int c = 0; c < 3; ++c)
                known->b[3 * b + c] = source.scalar < 0 ? 0.0 : source.sign * boundary->b[source.scalar + c];
        }
    }
    void BuildConvolution(Wvector& out, bool traction, long lag, size_t id)
    {
        const Wvector* source[3] = { NULL, NULL, NULL };
        for (int k = 0; k < 3; ++k)
        {
            const long previous = step - lag + 1 - k;
            if (previous > 0 && previous < step)
                source[k] = traction ? &(*states)[(size_t)previous].globalT : &(*states)[(size_t)previous].globalU;
        }
        for (long i = out.n * (long)id / (long)count; i < out.n * (long)(id + 1) / (long)count; ++i)
        {
            double value = 0.0;
            for (int k = 0; k < 3; ++k)
                if (source[k] && i < source[k]->n) value += (k == 1 ? 2.0 : 1.0) * source[k]->b[i];
            out.b[i] = value;
        }
    }
    template<class Matrix>
    void ProductAdd(Matrix& matrix, Wvector& input, const PthArg& range, double scale, Wvector* diagnostic)
    {
        // Complete each original row dot product locally before the RHS addition.
        for (long row = range.start; row < range.end; ++row)
        {
            double values[3];
            RhsProductRow(matrix, input, row, values);
            for (int c = 0; c < 3; ++c)
            {
                const long i = row * 3 + c;
                rhs->b[i] += scale * values[c];
                if (diagnostic) diagnostic->b[i] += scale * values[c];
            }
        }
    }
    template<class Matrix>
    static PthArg Range(const MultiDomainHistorySpMVTask<Matrix>& task)
    {
        PthArg range; range.start = task.start; range.end = task.end; return range;
    }
    static void Execute(void* data, size_t id) { static_cast<RhsCpuWorkspace*>(data)->RunWorker(id); }
    void RunWorker(size_t id)
    {
        const double taskStart = Clock();
        double start = Clock(); Gather(id); EndStage(id, KnownGather, start); Sync(id);
        start = Clock();
        MultiplyHistoryRows(*knownMatrix, *known, *rhs, knownRows[id].start, knownRows[id].end);
        const long begin = rhs->n * (long)id / (long)count, end = rhs->n * (long)(id + 1) / (long)count;
        // The known product partitions and vector partitions differ.
        EndStage(id, KnownMultiply, start); Sync(id);
        if (diagnostics) {
        for (long i = begin; i < end; ++i)
        {
            diagKnown->b[i] = rhs->b[i]; diagG->b[i] = 0.0; diagT->b[i] = 0.0;
        }
        // Every following phase has exactly one owner per output element.
        Sync(id);
        }
        for (long lag = 1; step > 1 && lag <= step && lag <= maxLag; ++lag)
        {
            start = Clock(); BuildConvolution(*convT, true, lag, id); BuildConvolution(*convU, false, lag, id);
            EndStage(id, Convolution, start); Sync(id);
            start = Clock();
            if (symmetric[lag]) ProductAdd(mgs[lag], *convT, Range((*gTasks)[(size_t)lag][id]), 1.0, diagnostics ? diagG : NULL);
            else ProductAdd(mgsFull[lag], *convT, Range((*gFullTasks)[(size_t)lag][id]), 1.0, diagnostics ? diagG : NULL);
            EndStage(id, GMultiplyAdd, start); Sync(id);
            start = Clock();
            ProductAdd(mts[lag], *convU, Range((*tTasks)[(size_t)lag][id]), -1.0, diagnostics ? diagT : NULL);
            EndStage(id, TMultiplyAdd, start); Sync(id);
        }
        start = Clock(); stats[id].firstInvalid = -1;
        for (long i = begin; i < end; ++i)
            if (!IsFiniteScalar(rhs->b[i]) && stats[id].firstInvalid < 0) stats[id].firstInvalid = i;
        for (long i = x0->n * (long)id / (long)count; i < x0->n * (long)(id + 1) / (long)count; ++i) x0->b[i] = 0.0;
        for (long i = x->n * (long)id / (long)count; i < x->n * (long)(id + 1) / (long)count; ++i) x->b[i] = 0.0;
        EndStage(id, CheckZero, start);
        if (profile) stats[id].lastTaskWall = DBEMWallTime() - taskStart;
    }
    void Finish()
    {
        long invalid = -1;
        for (size_t i = 0; i < count; ++i)
            if (stats[i].firstInvalid >= 0 && (invalid < 0 || stats[i].firstInvalid < invalid)) invalid = stats[i].firstInvalid;
        if (invalid >= 0) InvalidAssemblyValue("non-finite RHS before GMRES", invalid + 1, step);
    }
    void WriteProfile()
    {
        if (!profile) return;
        const char* names[StageCount] = { "rhs.cpu.known_gather", "rhs.cpu.known_spmv", "rhs.cpu.convolution",
            "rhs.cpu.g_spmv_accumulate", "rhs.cpu.t_spmv_accumulate", "rhs.cpu.check_zero", "rhs.cpu.worker_wait" };
        for (int stage = 0; stage < StageCount; ++stage)
        {
            double maximum = 0.0;
            for (size_t id = 0; id < count; ++id) maximum = (std::max)(maximum, stats[id].seconds[stage]);
            DBEMWriteTiming(names[stage], maximum);
        }
    }
};

template<class Matrix>
static void PrepareRhsSpMVTasks(Matrix& matrix, Wvector& input, Wvector& output,
    std::vector<MultiDomainHistorySpMVTask<Matrix> >& tasks, size_t count)
{
    std::vector<PthArg> ranges;
    RhsCpuWorkspace::Partition(matrix, count, ranges);
    tasks.resize(count);
    for (size_t id = 0; id < count; ++id)
    {
        tasks[id].matrix = &matrix; tasks[id].input = &input; tasks[id].output = &output;
        tasks[id].start = ranges[id].start; tasks[id].end = ranges[id].end;
    }
}

#ifdef DBEM_GMRES_CPU_TESTS
#include "rhs_cpu_tests_impl.h"
#endif
