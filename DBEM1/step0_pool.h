#pragma once
// Internal implementation. Included after Step0WorkerContext and its row kernel.
// Task IDs never select integral scratch slots: each sleeping worker owns one.
static auto g_step0CreateThread = &pthread_create;
static auto g_step0ExecuteRows = &MultiDomainStep0ThreadMain;

class Step0ThreadPool {
    struct Worker { Step0ThreadPool* pool; size_t index; };
    const MultiDomainModel& model;
    const GlobalDofMap& dofMap;
    std::vector<Step0WorkerContext> contexts;
    std::vector<Worker> workers;
    std::vector<pthread_t> threads;
    pthread_mutex_t mutex;
    pthread_cond_t ready, done;
    bool mutexReady = false, readyInitialized = false, doneInitialized = false;
    bool stopping = false, initialized = false, batchSucceeded = false;
    size_t created = 0, generation = 0, completed = 0, taskCount = 0;

    static void* Work(void* data)
    {
        auto& worker = *static_cast<Worker*>(data);
        auto& pool = *worker.pool;
        size_t seen = 0;
        RequireThreadSuccess(pthread_mutex_lock(&pool.mutex), "Step0 lock");
        for (;;) {
            while (!pool.stopping && seen == pool.generation)
                RequireThreadSuccess(pthread_cond_wait(&pool.ready, &pool.mutex), "Step0 wait");
            if (pool.stopping) break;
            seen = pool.generation;
            const bool active = worker.index < pool.taskCount;
            RequireThreadSuccess(pthread_mutex_unlock(&pool.mutex), "Step0 unlock");
            if (active) g_step0ExecuteRows(&pool.contexts[worker.index]);
            RequireThreadSuccess(pthread_mutex_lock(&pool.mutex), "Step0 lock");
            if (++pool.completed == pool.created)
                RequireThreadSuccess(pthread_cond_signal(&pool.done), "Step0 done signal");
        }
        RequireThreadSuccess(pthread_mutex_unlock(&pool.mutex), "Step0 unlock");
        return nullptr;
    }

public:
    Step0ThreadPool(DSquareElement* elements, const MultiDomainModel& m,
        const GlobalDofMap& dofs, const MultiDomainInterfaceCache& cache,
        long** infElePid, long** elePid, long count) : model(m), dofMap(dofs)
    {
        if (!elements || count < 1 || count > thread_num || !AssistCoef || !AssistUij || !AssistTij)
            throw std::runtime_error("Step0 invalid worker count/integral workspace");
        contexts.resize((size_t)count);
        workers.resize((size_t)count);
        threads.resize((size_t)count);
        for (long i = 0; i < count; ++i) {
            auto& c = contexts[(size_t)i];
            c.elements = elements;
            c.dofMap = &dofs;
            c.interfaceCache = &cache;
            c.infElePid = infElePid; c.elePid = elePid;
            c.threadId = i;
            c.coefSlot = AssistCoef[i]; c.uSlot = AssistUij[i]; c.tSlot = AssistTij[i];
            if (!c.coefSlot || !c.uSlot || !c.tSlot)
                throw std::runtime_error("Step0 missing integral workspace slot");
            c.unknownBuilder = &c.unknownStorage;
            c.knownBuilder = &c.knownStorage;
            workers[(size_t)i] = {this, (size_t)i};
        }
    }
    Step0ThreadPool(const Step0ThreadPool&) = delete;
    Step0ThreadPool& operator=(const Step0ThreadPool&) = delete;
    ~Step0ThreadPool() { Close(); }

    bool Initialize()
    {
        if (initialized) return true;
        stopping = false;
        generation = completed = 0;
        if (contexts.size() == 1) { initialized = true; return true; }
        int error = pthread_mutex_init(&mutex, nullptr);
        if (!error) { mutexReady = true; error = pthread_cond_init(&ready, nullptr); }
        if (!error) { readyInitialized = true; error = pthread_cond_init(&done, nullptr); }
        if (!error) {
            doneInitialized = true;
            for (; created < threads.size(); ++created) {
                error = g_step0CreateThread(&threads[created], nullptr, Work, &workers[created]);
                if (error) break;
            }
        }
        if (error) {
            fprintf(stderr, "Step0 pool initialization failed: code=%d created=%zu\n", error, created);
            Close();
            return false;
        }
        initialized = true;
        return true;
    }

    bool RunBatch(const std::vector<Step0Task>& tasks)
    {
        batchSucceeded = false;
        if (!initialized || tasks.size() > contexts.size())
            throw std::runtime_error("Step0 invalid task batch/pool state");
        taskCount = tasks.size();
        for (size_t i = 0; i < taskCount; ++i) {
            const auto& task = tasks[i];
            if (task.domainId < 0 || (size_t)task.domainId >= model.domains.size())
                throw std::runtime_error("Step0 invalid task domain");
            auto& c = contexts[i];
            c.domain = &model.domains[(size_t)task.domainId];
            if (task.sourceBegin < 0 || task.sourceEnd < task.sourceBegin ||
                task.sourceEnd > (long)c.domain->elementIds.size())
                throw std::runtime_error("Step0 invalid source range");
            c.sourceBegin = task.sourceBegin; c.sourceEnd = task.sourceEnd;
            c.unknownStorage.Reset(dofMap.equationBlockCount, dofMap.unknownBlockCount);
            c.knownStorage.Reset(dofMap.equationBlockCount, dofMap.knownBlockCount);
            c.ok = 1; c.blockCount = 0;
        }
        if (contexts.size() == 1) {
            if (taskCount) g_step0ExecuteRows(&contexts[0]);
        } else {
            RequireThreadSuccess(pthread_mutex_lock(&mutex), "Step0 batch lock");
            completed = 0;
            ++generation;
            RequireThreadSuccess(pthread_cond_broadcast(&ready), "Step0 batch broadcast");
            while (completed != created)
                RequireThreadSuccess(pthread_cond_wait(&done, &mutex), "Step0 batch wait");
            RequireThreadSuccess(pthread_mutex_unlock(&mutex), "Step0 batch unlock");
        }
        for (size_t i = 0; i < taskCount; ++i)
            if (contexts[i].ok != 1) return false;
        batchSucceeded = true;
        return true;
    }

    void MergeInto(MultiDomainCCSRBuilder& unknown, MultiDomainCCSRBuilder& known) const
    {
        if (!batchSucceeded) throw std::runtime_error("Step0 cannot publish a failed batch");
        for (size_t i = 0; i < taskCount; ++i) {
            unknown.MergeFrom(contexts[i].unknownStorage);
            known.MergeFrom(contexts[i].knownStorage);
        }
    }

    void Close()
    {
        if (created) {
            RequireThreadSuccess(pthread_mutex_lock(&mutex), "Step0 close lock");
            stopping = true;
            RequireThreadSuccess(pthread_cond_broadcast(&ready), "Step0 close broadcast");
            RequireThreadSuccess(pthread_mutex_unlock(&mutex), "Step0 close unlock");
            for (size_t i = 0; i < created; ++i)
                RequireThreadSuccess(pthread_join(threads[i], nullptr), "Step0 join");
        }
        created = 0;
        if (doneInitialized) RequireThreadSuccess(pthread_cond_destroy(&done), "Step0 destroy done");
        if (readyInitialized) RequireThreadSuccess(pthread_cond_destroy(&ready), "Step0 destroy ready");
        if (mutexReady) RequireThreadSuccess(pthread_mutex_destroy(&mutex), "Step0 destroy mutex");
        initialized = mutexReady = readyInitialized = doneInitialized = false;
    }
};
