// Build as a separate test translation unit, replacing multidomain.cpp and
// renaming DBEM.cpp's main to dbem_application_main (see validation runner).
#define DBEM_STEP0_TEST_HOOKS
#include "../multidomain.cpp"
#include <stdexcept>

#define CHECK(expr) do { if (!(expr)) throw std::runtime_error(#expr); } while (0)

static int createCalls = 0;
static int failAt = -1;
static int PTW32_CDECL InjectCreate(pthread_t* thread, const pthread_attr_t* attr,
	void* (*entry)(void*), void* arg)
{
	if (createCalls++ == failAt) return EAGAIN;
	return pthread_create(thread, attr, entry, arg);
}

static void SameMatrix(const CCSRMat& a, const CCSRMat& b)
{
	CHECK(a.n == b.n);
	for (long i = 0; i <= a.n; ++i) CHECK(a.row_ptr[i] == b.row_ptr[i]);
	for (long i = 0; i < a.row_ptr[a.n]; ++i)
	{
		CHECK(a.col_index[i] == b.col_index[i]);
		for (int k = 0; k < 9; ++k) CHECK(a.value[(size_t)i * 9 + k] == b.value[(size_t)i * 9 + k]);
	}
}

static void TestOutput(bool nonSymmetric, bool empty, int failedCreate)
{
	using Builder = MultiDomainCCSRBuilder;
	Builder mts(7, 7, false, Builder::IndependentRows), mgs(7, 7, false, Builder::IndependentRows);
	Builder refT(7, 7, false), refG(7, 7, false);
	const double diagonal[9] = { 1,0,0,0,2,0,0,0,3 };
	const double asymmetric[9] = { 1,2,3,4,5,6,7,8,9 };
	if (!empty)
	{
		for (long row : { 0L, 2L, 6L })
			for (long col : { 5L, 1L, 3L })
			{
				const double* g = nonSymmetric ? asymmetric : diagonal;
				mts.AddBlock(row, col, 1, diagonal); refT.AddBlock(row, col, 1, diagonal);
				mgs.AddBlock(row, col, -1, g); refG.AddBlock(row, col, -1, g);
			}
		// Repeated columns accumulate and exact cancellation removes the block.
		mts.AddBlock(2, 3, -1, diagonal); refT.AddBlock(2, 3, -1, diagonal);
		mts.AddBlock(6, 1, 2, diagonal); refT.AddBlock(6, 1, 2, diagonal);
	}
	CHECK(SameSparseBuilderStructure(mts, refT));
	CHECK(SparseBuilderDifference(mts, refT) == 0);
	CHECK(SameSparseBuilderStructure(mgs, refG));
	CHECK(SparseBuilderDifference(mgs, refG) == 0);
	CCSRMat expectedT, expectedG, actualT, actualG;
	SymCCSRMat expectedSym, actualSym;
	actualG.initial(0, 0);
	refT.Build(expectedT); refG.Build(expectedG);
	bool expectedSymmetric = true;
	refG.BuildSym(expectedSym, 1e-10, &expectedSymmetric);
	createCalls = 0; failAt = failedCreate;
	g_assemblyCreateThread = InjectCreate;
	bool symmetric = false;
	const bool ok = BuildHistoryMatrices(mts, mgs, actualT, actualSym, actualG, symmetric, 4, 1);
    if (failedCreate >= 0 && !DBEMPolicy::allowSerialFallback) {
        CHECK(!ok);
        g_assemblyCreateThread = &pthread_create;
        return;
    }
    CHECK(ok);
    CHECK(createCalls == (failedCreate >= 0 ? failedCreate + 1 : 4));
	g_assemblyCreateThread = &pthread_create;
	CHECK(symmetric == expectedSymmetric);
	SameMatrix(expectedT, actualT);
	if (!symmetric) SameMatrix(expectedG, actualG);
	else
	{
		CHECK(actualSym.n == expectedSym.n);
		for (long row = 0; row <= actualSym.n; ++row) CHECK(actualSym.row_ptr[row] == expectedSym.row_ptr[row]);
		for (long i = 0; i < actualSym.row_ptr[actualSym.n]; ++i)
		{
			CHECK(actualSym.col_index[i] == expectedSym.col_index[i]);
			for (int k = 0; k < 6; ++k) CHECK(actualSym.value[(size_t)i * 6 + k] == expectedSym.value[(size_t)i * 6 + k]);
		}
		CHECK(actualG.n == 0);
	}
	CHECK(mts.NonZeroBlocks() == 0 && mgs.NonZeroBlocks() == 0);
}

static void TestOwnershipAndPoolFallback()
{
	MultiDomainModel model;
	model.nodeCount = model.elementCount = 2;
	Domain domain = {};
	domain.id = 0; domain.eleBegin = 0; domain.eleCount = 2;
	domain.elementIds = { 0, 1 };
	model.domains.push_back(domain);
	DSquareElement elements[2];
	for (long i = 0; i < 2; ++i) { elements[i].BCID = 123; elements[i].m_nodeID[0] = i; }
	GlobalDofMap dof;
	CHECK(dof.Build(model, elements, 2));
	CHECK(HistoryRowsHaveUniqueOwners(elements, model, dof));
	elements[1].m_nodeID[0] = 0;
	CHECK(!HistoryRowsHaveUniqueOwners(elements, model, dof));

    MultiDomainCCSRBuilder mts, mgs;
    CHECK(!AssembleMultiDomainHistoryStepPthread(elements, model, dof, 1, mts, mgs, 2, false));
    createCalls = 0; failAt = 1;
    g_historyCreateThread = InjectCreate;
    MultiDomainHistoryThreadPool pool;
    CHECK(pool.Initialize(2) == DBEMPolicy::allowSerialFallback);
    pool.Close();
    g_historyCreateThread = &pthread_create;
}

static void TestBuilderCompatibility()
{
	using Builder = MultiDomainCCSRBuilder;
	Builder rows(3, 3, false, Builder::IndependentRows), legacy(3, 3, false);
	const double diagonal[9] = { 1,0,0,0,2,0,0,0,3 };
	rows.AddBlock(1, 2, 1, diagonal);
	legacy.MergeFrom(rows);
	CHECK(SparseBuilderDifference(rows, legacy) == 0);
	Builder rows2(3, 3, false, Builder::IndependentRows);
	rows2.MergeFrom(legacy);
	CHECK(SparseBuilderDifference(rows, rows2) == 0);
	std::vector<double> a, b;
	rows.BuildDense(a); legacy.BuildDense(b); CHECK(a == b);
	const double rhs[3] = { 1, 2, 3 };
	rows.AddRhsBlock(1, -1, rhs);
	Wvector v;
	rows.BuildRhs(v); CHECK(v[4] == -1 && v[5] == -2 && v[6] == -3);
    double bad[9] = { (std::numeric_limits<double>::quiet_NaN)(),0,0,0,0,0,0,0,0 };
    bool rejected = false;
    try { rows.AddBlock(0, 0, 1, bad); } catch (const std::runtime_error&) { rejected = true; }
    CHECK(rejected && !rows.FindBlock(0, 0));
	rows.Reset(3, 3, false); CHECK(!rows.HasIndependentRows() && rows.NonZeroBlocks() == 0);
}

#include "interface_cache_test.h"
#include "step0_pool_test.h"
#include "assembly_pool_test.h"

int main()
{
	setvbuf(stdout, 0, _IONBF, 0);
	try
	{
		printf("TEST interface cache\n");
		TestInterfaceCache();
		printf("TEST Step0 pool and explicit BCID\n");
		TestStep0Pool();
        TestAssemblyQueue();
		printf("TEST builder compatibility\n");
		TestBuilderCompatibility();
		printf("TEST symmetric output\n");
		TestOutput(false, false, -1);
		printf("TEST full output\n");
		TestOutput(true, false, -1);
		printf("TEST empty output\n");
		TestOutput(false, true, -1);
		printf("TEST output creation failures\n");
		// Partial pool creation failure; output phases reuse the same workers.
		TestOutput(false, false, 1);
		TestOutput(true, false, 1);
		printf("TEST ownership and assembly fallback\n");
		TestOwnershipAndPoolFallback();
		printf("PASS: history row storage, sparse export, ownership and failure recovery tests\n");
		return 0;
	}
	catch (const std::exception& error)
	{
		fprintf(stderr, "FAIL: %s\n", error.what());
		return 1;
	}
}
