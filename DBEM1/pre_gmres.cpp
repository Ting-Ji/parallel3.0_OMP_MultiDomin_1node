#include "pre_gmres.h"
#include "stdio.h"
#include "math.h"
#include "Counter.h"
#include "Gauss_solver.h"
#include <string>
#include <string.h>
#include <stdlib.h>
#include <direct.h>
#include "DBEM.h"
#include "PthreadCCSR.h"
//
//-------------------------------------------------------------------------------------------------

static const char* GMRES_CCSR_CACHE_PATH = "output\\gmres_ccsr_cache.bin";
static const int GMRES_CCSR_CACHE_VERSION = 1;
static const int GMRES_CACHE_MTS0 = 1;
static const int GMRES_CACHE_MTS = 2;
static const int GMRES_CACHE_MGS = 3;
static bool GMRES_BATCH_MEMORY_SESSION = false;
static bool GMRES_BATCH_MEMORY_SESSION_READY = false;
static long GMRES_BATCH_NODE_NUM = 0;
static long GMRES_BATCH_ELE_NUM = 0;
static long GMRES_BATCH_NSTEP = 0;
static long GMRES_BATCH_MAXN = 0;
static long GMRES_BATCH_THREAD_NUM = 0;
static char GMRES_OUTPUT_SUBDIR[260] = "";
static int RhsDiagnosticEnabled()
{
	char* value = 0;
	size_t valueLen = 0;
	_dupenv_s(&value, &valueLen, "DBEM_RHS_DIAGNOSTIC");
	int enabled = value && (strcmp(value, "1") == 0 || strcmp(value, "true") == 0 || strcmp(value, "TRUE") == 0);
	if (value)
		free(value);
	return enabled;
}

static FILE* OpenDiagnosticCsv(const char* path, const char* header)
{
	FILE* test = 0;
	fopen_s(&test, path, "r");
	int needsHeader = test ? 0 : 1;
	if (test)
		fclose(test);

	FILE* fp = 0;
	fopen_s(&fp, path, "a");
	if (fp && needsHeader && header)
		fprintf(fp, "%s\n", header);
	return fp;
}

static void WriteOldRhsBreakdownCsv(const char* solver,
	long step,
	DSquareElement* elements,
	long eleCount,
	Wvector& known,
	Wvector& historyG,
	Wvector& historyT,
	Wvector& total,
	Wvector& solution,
	CCSRMat& matrix)
{
	if (!RhsDiagnosticEnabled() || !elements)
		return;

	Wvector ax(total.n, 0);
	SparseMul(matrix, solution, ax);

	FILE* fp = OpenDiagnosticCsv("output\\rhs_breakdown.csv",
		"solver,step,rowElement,domain,surfaceType,bcid,cx,cy,cz,component,knownRhs,historyG,historyT,totalRhs,Ax,residual,solution");
	if (!fp)
		return;

	const char comp[3] = { 'x', 'y', 'z' };
	for (long ele = 0; ele < eleCount; ++ele)
	{
		for (int k = 0; k < 3; ++k)
		{
			long p = ele * 3 + k;
			double residual = ax.b[p] - total.b[p];
			fprintf(fp, "%s,%ld,%ld,%d,%d,%d,%.17g,%.17g,%.17g,%c,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
				solver,
				step,
				ele,
				elements[ele].DomainID,
				elements[ele].SurfaceType,
				elements[ele].BCID,
				elements[ele].m_center.pt[0],
				elements[ele].m_center.pt[1],
				elements[ele].m_center.pt[2],
				comp[k],
				known.b[p],
				historyG.b[p],
				historyT.b[p],
				total.b[p],
				ax.b[p],
				residual,
				solution.b[p]);
		}
	}
	fclose(fp);

	double maxAbsResidual = 0.0;
	double rhsNorm2 = 0.0;
	double residualNorm2 = 0.0;
	for (long i = 0; i < total.n; ++i)
	{
		double residual = ax.b[i] - total.b[i];
		if (fabs(residual) > maxAbsResidual)
			maxAbsResidual = fabs(residual);
		residualNorm2 += residual * residual;
		rhsNorm2 += total.b[i] * total.b[i];
	}
	FILE* rp = OpenDiagnosticCsv("output\\residual_probe.csv",
		"solver,step,maxAbsResidual,l2Residual,l2Rhs,relativeResidual");
	if (rp)
	{
		double l2Residual = sqrt(residualNorm2);
		double l2Rhs = sqrt(rhsNorm2);
		double relative = l2Rhs > 0.0 ? l2Residual / l2Rhs : l2Residual;
		fprintf(rp, "%s,%ld,%.17g,%.17g,%.17g,%.17g\n",
			solver, step, maxAbsResidual, l2Residual, l2Rhs, relative);
		fclose(rp);
	}
}
struct GmresCcsrCacheHeader
{
	char magic[16];
	int version;
	int matrixCount;
	long NodeNum;
	long EleNum;
	long NStep;
	long MaxN;
	double dt;
	double E;
	double v;
	double Rou;
	double gap;
};

struct GmresCcsrCacheEntry
{
	int type;
	long step;
	__int64 offset;
	__int64 byteSize;
};

void GmresCcsrBeginMemorySession()
{
	GMRES_BATCH_MEMORY_SESSION = true;
	GMRES_BATCH_MEMORY_SESSION_READY = false;
}

static void GmresCcsrReleaseRetainedMatrices()
{
	if (SpMVThreadT)
	{
		for (long i = 0; i <= GMRES_BATCH_MAXN; ++i)
		{
			delete[] SpMVThreadT[i];
			SpMVThreadT[i] = 0;
		}
		delete[] SpMVThreadT;
		SpMVThreadT = 0;
	}
	if (SpMVThreadG)
	{
		for (long i = 0; i <= GMRES_BATCH_MAXN; ++i)
		{
			delete[] SpMVThreadG[i];
			SpMVThreadG[i] = 0;
		}
		delete[] SpMVThreadG;
		SpMVThreadG = 0;
	}
	delete[] ConvolThreadT;
	ConvolThreadT = 0;
	delete[] ConvolThreadG;
	ConvolThreadG = 0;
	delete[] Thread;
	Thread = 0;
	delete[] Pthnum;
	Pthnum = 0;
	delete[] MGS;
	MGS = 0;
	delete[] MTS;
	MTS = 0;
	MTS0.Matdelete();
	GMRES_BATCH_MEMORY_SESSION_READY = false;
	GMRES_BATCH_NODE_NUM = 0;
	GMRES_BATCH_ELE_NUM = 0;
	GMRES_BATCH_NSTEP = 0;
	GMRES_BATCH_MAXN = 0;
	GMRES_BATCH_THREAD_NUM = 0;
}

void GmresCcsrEndMemorySession()
{
	GmresCcsrReleaseRetainedMatrices();
	GMRES_BATCH_MEMORY_SESSION = false;
	GMRES_OUTPUT_SUBDIR[0] = '\0';
}

void GmresCcsrSetOutputSubdir(const char* subdir)
{
	if (!subdir || subdir[0] == '\0')
	{
		GMRES_OUTPUT_SUBDIR[0] = '\0';
		return;
	}
	strcpy_s(GMRES_OUTPUT_SUBDIR, sizeof(GMRES_OUTPUT_SUBDIR), subdir);
}

static bool GmresCcsrCanReuseMemorySession(long NodeNum, long EleNum, long NStep, long MaxN)
{
	return GMRES_BATCH_MEMORY_SESSION && GMRES_BATCH_MEMORY_SESSION_READY &&
		GMRES_BATCH_NODE_NUM == NodeNum && GMRES_BATCH_ELE_NUM == EleNum &&
		GMRES_BATCH_NSTEP == NStep && GMRES_BATCH_MAXN == MaxN &&
		GMRES_BATCH_THREAD_NUM == thread_num;
}

static void GmresCcsrMarkMemorySessionReady(long NodeNum, long EleNum, long NStep, long MaxN)
{
	if (!GMRES_BATCH_MEMORY_SESSION)
		return;
	GMRES_BATCH_MEMORY_SESSION_READY = true;
	GMRES_BATCH_NODE_NUM = NodeNum;
	GMRES_BATCH_ELE_NUM = EleNum;
	GMRES_BATCH_NSTEP = NStep;
	GMRES_BATCH_MAXN = MaxN;
	GMRES_BATCH_THREAD_NUM = thread_num;
}
static double GmresCacheElasticE()
{
	return 2.0 * DSquareElement::m_DMat.G * (1.0 + DSquareElement::m_DMat.v);
}

static int GmresCacheMatrixCount(long MaxN)
{
	return (int)(2 * MaxN + 2);
}

static GmresCcsrCacheHeader GmresMakeCacheHeader(long NodeNum, long EleNum, long NStep, long MaxN)
{
	GmresCcsrCacheHeader header;
	memset(&header, 0, sizeof(header));
	strcpy_s(header.magic, sizeof(header.magic), "GMRESCCSRCACHE");
	header.version = GMRES_CCSR_CACHE_VERSION;
	header.matrixCount = GmresCacheMatrixCount(MaxN);
	header.NodeNum = NodeNum;
	header.EleNum = EleNum;
	header.NStep = NStep;
	header.MaxN = MaxN;
	header.dt = DSquareElement::m_DMat.Dt;
	header.E = GmresCacheElasticE();
	header.v = DSquareElement::m_DMat.v;
	header.Rou = DSquareElement::m_DMat.Rou;
	header.gap = DSquareElement::m_gap;
	return header;
}

static bool GmresCacheHeaderMatches(const GmresCcsrCacheHeader& actual, const GmresCcsrCacheHeader& expected, string& reason)
{
	if (strcmp(actual.magic, expected.magic) != 0) { reason = "magic mismatch"; return false; }
	if (actual.version != expected.version) { reason = "version mismatch"; return false; }
	if (actual.matrixCount != expected.matrixCount) { reason = "matrix count mismatch"; return false; }
	if (actual.NodeNum != expected.NodeNum) { reason = "NodeNum mismatch"; return false; }
	if (actual.EleNum != expected.EleNum) { reason = "EleNum mismatch"; return false; }
	if (actual.NStep != expected.NStep) { reason = "NStep mismatch"; return false; }
	if (actual.MaxN != expected.MaxN) { reason = "MaxN mismatch"; return false; }
	if (actual.dt != expected.dt) { reason = "dt mismatch"; return false; }
	if (actual.E != expected.E) { reason = "E mismatch"; return false; }
	if (actual.v != expected.v) { reason = "v mismatch"; return false; }
	if (actual.Rou != expected.Rou) { reason = "Rou mismatch"; return false; }
	if (actual.gap != expected.gap) { reason = "gap mismatch"; return false; }
	return true;
}

static GmresCcsrCacheEntry* GmresFindCacheEntry(GmresCcsrCacheEntry* entries, int matrixCount, int type, long step)
{
	for (int i = 0; i < matrixCount; ++i)
	{
		if (entries[i].type == type && entries[i].step == step)
			return &entries[i];
	}
	return 0;
}

static bool GmresWriteCacheEntry(FILE* fp, GmresCcsrCacheEntry* entry, int type, long step, CCSRMat& mat)
{
	entry->type = type;
	entry->step = step;
	entry->offset = _ftelli64(fp);
	if (entry->offset < 0)
		return false;
	if (!mat.WriteBinary(fp))
		return false;
	__int64 end = _ftelli64(fp);
	if (end < entry->offset)
		return false;
	entry->byteSize = end - entry->offset;
	return true;
}

static bool GmresWriteCacheEntry(FILE* fp, GmresCcsrCacheEntry* entry, int type, long step, SymCCSRMat& mat)
{
	entry->type = type;
	entry->step = step;
	entry->offset = _ftelli64(fp);
	if (entry->offset < 0)
		return false;
	if (!mat.WriteBinary(fp))
		return false;
	__int64 end = _ftelli64(fp);
	if (end < entry->offset)
		return false;
	entry->byteSize = end - entry->offset;
	return true;
}

static bool GmresReadCacheEntry(FILE* fp, GmresCcsrCacheEntry* entry, CCSRMat& mat)
{
	if (!entry || _fseeki64(fp, entry->offset, SEEK_SET) != 0)
		return false;
	if (!mat.ReadBinary(fp))
		return false;
	__int64 end = _ftelli64(fp);
	return end - entry->offset == entry->byteSize;
}

static bool GmresReadCacheEntry(FILE* fp, GmresCcsrCacheEntry* entry, SymCCSRMat& mat)
{
	if (!entry || _fseeki64(fp, entry->offset, SEEK_SET) != 0)
		return false;
	if (!mat.ReadBinary(fp))
		return false;
	__int64 end = _ftelli64(fp);
	return end - entry->offset == entry->byteSize;
}

static void GmresBuildRowCount(int* rowCount, int* rowPtr, long n)
{
	for (long i = 0; i < n; ++i)
		rowCount[i] = rowPtr[i + 1] - rowPtr[i];
}

static bool GmresRebuildSpmvTasks(long MaxN, long threadCount)
{
	if (!Thread || !SpMVThreadT || !SpMVThreadG || threadCount <= 0)
		return false;
	int* rowCount = new int[MTS0.n];
	GmresBuildRowCount(rowCount, MTS0.row_ptr, MTS0.n);
	PthDistr(rowCount, Thread, MTS0.NonZeroBlocks(), MTS0.n);
	for (long i = 0; i < threadCount; ++i)
	{
		SpMVThreadT[0][i] = Thread[i];
		SpMVThreadG[0][i] = Thread[i];
	}
	delete[] rowCount;

	for (long step = 1; step <= MaxN; ++step)
	{
		rowCount = new int[MTS[step].n];
		GmresBuildRowCount(rowCount, MTS[step].row_ptr, MTS[step].n);
		PthDistr(rowCount, Thread, MTS[step].NonZeroBlocks(), MTS[step].n);
		for (long i = 0; i < threadCount; ++i)
			SpMVThreadT[step][i] = Thread[i];
		delete[] rowCount;

		rowCount = new int[MGS[step].n];
		GmresBuildRowCount(rowCount, MGS[step].row_ptr, MGS[step].n);
		PthDistr(rowCount, Thread, MGS[step].NonZeroBlocks(), MGS[step].n);
		for (long i = 0; i < threadCount; ++i)
			SpMVThreadG[step][i] = Thread[i];
		delete[] rowCount;
	}
	return true;
}

static bool GmresWriteCcsrCache(const char* cachePath, long NodeNum, long EleNum, long NStep, long MaxN, __int64& fileSize, string& reason)
{
	_mkdir("output");
	FILE* fp = 0;
	if (fopen_s(&fp, cachePath, "wb") != 0 || !fp)
	{
		reason = "cannot open cache for write";
		return false;
	}

	GmresCcsrCacheHeader header = GmresMakeCacheHeader(NodeNum, EleNum, NStep, MaxN);
	int matrixCount = header.matrixCount;
	GmresCcsrCacheEntry* entries = new GmresCcsrCacheEntry[matrixCount];
	memset(entries, 0, sizeof(GmresCcsrCacheEntry) * matrixCount);

	bool ok = fwrite(&header, sizeof(header), 1, fp) == 1;
	ok = ok && fwrite(entries, sizeof(GmresCcsrCacheEntry), matrixCount, fp) == (size_t)matrixCount;
	int index = 0;
	ok = ok && GmresWriteCacheEntry(fp, &entries[index++], GMRES_CACHE_MTS0, 0, MTS0);
	for (long step = 0; ok && step <= MaxN; ++step)
		ok = GmresWriteCacheEntry(fp, &entries[index++], GMRES_CACHE_MTS, step, MTS[step]);
	for (long step = 1; ok && step <= MaxN; ++step)
		ok = GmresWriteCacheEntry(fp, &entries[index++], GMRES_CACHE_MGS, step, MGS[step]);

	fileSize = _ftelli64(fp);
	if (ok && _fseeki64(fp, sizeof(header), SEEK_SET) == 0)
		ok = fwrite(entries, sizeof(GmresCcsrCacheEntry), matrixCount, fp) == (size_t)matrixCount;

	delete[] entries;
	fclose(fp);
	if (!ok)
		reason = "cache write failed";
	return ok;
}

static bool GmresReadCcsrCache(const char* cachePath, long NodeNum, long EleNum, long NStep, long MaxN, long threadCount, string& reason)
{
	FILE* fp = 0;
	if (fopen_s(&fp, cachePath, "rb") != 0 || !fp)
	{
		reason = "cache file not found";
		return false;
	}

	GmresCcsrCacheHeader expected = GmresMakeCacheHeader(NodeNum, EleNum, NStep, MaxN);
	GmresCcsrCacheHeader actual;
	bool ok = fread(&actual, sizeof(actual), 1, fp) == 1;
	if (!ok)
		reason = "cannot read cache header";
	if (ok)
		ok = GmresCacheHeaderMatches(actual, expected, reason);

	GmresCcsrCacheEntry* entries = 0;
	if (ok)
	{
		entries = new GmresCcsrCacheEntry[actual.matrixCount];
		ok = fread(entries, sizeof(GmresCcsrCacheEntry), actual.matrixCount, fp) == (size_t)actual.matrixCount;
		if (!ok)
			reason = "cannot read cache index";
	}
	if (ok)
	{
		ok = GmresReadCacheEntry(fp, GmresFindCacheEntry(entries, actual.matrixCount, GMRES_CACHE_MTS0, 0), MTS0);
		if (!ok)
			reason = "cannot read MTS0";
	}
	for (long step = 0; ok && step <= MaxN; ++step)
	{
		ok = GmresReadCacheEntry(fp, GmresFindCacheEntry(entries, actual.matrixCount, GMRES_CACHE_MTS, step), MTS[step]);
		if (!ok)
			reason = "cannot read MTS";
	}
	for (long step = 1; ok && step <= MaxN; ++step)
	{
		ok = GmresReadCacheEntry(fp, GmresFindCacheEntry(entries, actual.matrixCount, GMRES_CACHE_MGS, step), MGS[step]);
		if (!ok)
			reason = "cannot read MGS";
	}
	if (ok)
	{
		ok = GmresRebuildSpmvTasks(MaxN, threadCount);
		if (!ok)
			reason = "cannot rebuild SpMV tasks";
	}

	if (entries)
		delete[] entries;
	fclose(fp);
	return ok;
}

int GMRESPreConditioner(DSquareElement* m_DSE, PreConditioner& m_Pre, long Count_PT, long MaxLeafPointCount, Wmatrix& A)
{
	Cube GlobalBox;
	Tree* m_treePre;
	PointerOfLeaf p_leafPre;

	long i, j, k, beginID;

	Wmatrix subPre(3, 3);

	long PointCount = Count_PT;

	// б
	AssignCubeSize(m_DSE[0].m_nodelist, PointCount, GlobalBox);

	// ṹ
	InitRoot(m_treePre, GlobalBox, PointCount);
	CreateTree(m_treePre, m_DSE[0].m_nodelist, MaxLeafPointCount);
	CreateLeafPointer(p_leafPre, m_treePre, PointCount);

	// Ԥ
	m_Pre.m_LeafCount = p_leafPre.LeafCount;
	m_Pre.m_LeafPointCount = new long[p_leafPre.LeafCount];
	m_Pre.m_LeafBeginID = new long[p_leafPre.LeafCount];
	m_Pre.m_RePID = new long[PointCount];
	m_Pre.IID = new long[PointCount];
	m_Pre.m_PreM = new Wmatrix[p_leafPre.LeafCount];

	for (i = 0; i < p_leafPre.LeafCount; ++i)
	{
		m_Pre.m_LeafPointCount[i] = p_leafPre.LeafPointer[i]->m_PointCount;
		m_Pre.m_LeafBeginID[i] = p_leafPre.LeafPointer[i]->m_BeginID;
	}

	// 
	RenumberPointID(p_leafPre, m_Pre.m_RePID);

	// IIDֵ
	for (i = 0; i < PointCount; ++i)
		m_Pre.IID[i] = 3 * i;

	// ʼԤ
	for (i = 0; i < p_leafPre.LeafCount; ++i)
	{
		j = 3 * m_Pre.m_LeafPointCount[i];
		m_Pre.m_PreM[i].initial(j, j);
	}

	// Ԥ
	long a1, a2;

	for (i = 0; i < m_Pre.m_LeafCount; ++i)
	{
		beginID = m_Pre.m_LeafBeginID[i];

		for (j = 0; j < m_Pre.m_LeafPointCount[i]; ++j)
		{
			a1 = 3 * m_Pre.m_RePID[beginID + j];

			for (k = 0; k < m_Pre.m_LeafPointCount[i]; ++k)
			{
				a2 = 3 * m_Pre.m_RePID[beginID + k];
				subPre(1, 1) = A(a1 + 1, a2 + 1);
				subPre(1, 2) = A(a1 + 1, a2 + 2);
				subPre(1, 3) = A(a1 + 1, a2 + 3);
				subPre(2, 1) = A(a1 + 2, a2 + 1);
				subPre(2, 2) = A(a1 + 2, a2 + 2);
				subPre(2, 3) = A(a1 + 2, a2 + 3);
				subPre(3, 1) = A(a1 + 3, a2 + 1);
				subPre(3, 2) = A(a1 + 3, a2 + 2);
				subPre(3, 3) = A(a1 + 3, a2 + 3);
				m_Pre.m_PreM[i].assign(subPre, 3 * j + 1, 3 * k + 1, 3, 3);
			}
		}
	}

	// ԤԤ
	for (i = 0; i < m_Pre.m_LeafCount; ++i)
		inv_mat(m_Pre.m_PreM[i], m_Pre.m_PreM[i].m);

	//release PreTree and p_leafPre
	DeleteTree(m_treePre);
	DeleteLeafPointer(p_leafPre);

	return 1;
}



int GMRESPreConditioner(DSquareElement* m_DSE, PreConditioner& m_Pre, long Count_PT, long MaxLeafPointCount, CCSRMat& M)
{
	Cube GlobalBox;
	Tree* m_treePre;
	PointerOfLeaf p_leafPre;

	long i, j;

	

	long PointCount = Count_PT;

	// б
	AssignCubeSize(m_DSE[0].m_nodelist, PointCount, GlobalBox);

	// ṹ
	InitRoot(m_treePre, GlobalBox, PointCount);
	CreateTree(m_treePre, m_DSE[0].m_nodelist, MaxLeafPointCount);
	CreateLeafPointer(p_leafPre, m_treePre, PointCount);

	// Ԥ
	m_Pre.m_LeafCount = p_leafPre.LeafCount;
	m_Pre.m_LeafPointCount = new long[p_leafPre.LeafCount];
	m_Pre.m_LeafBeginID = new long[p_leafPre.LeafCount];
	m_Pre.m_RePID = new long[PointCount];
	m_Pre.IID = new long[PointCount];
	m_Pre.m_PreM = new Wmatrix[p_leafPre.LeafCount];

	for (i = 0; i < p_leafPre.LeafCount; ++i)
	{
		m_Pre.m_LeafPointCount[i] = p_leafPre.LeafPointer[i]->m_PointCount;
		m_Pre.m_LeafBeginID[i] = p_leafPre.LeafPointer[i]->m_BeginID;
	}

	// 
	RenumberPointID(p_leafPre, m_Pre.m_RePID);

	// IIDֵ
	for (i = 0; i < PointCount; ++i)
		m_Pre.IID[i] = 3 * i;

	// ʼԤ
	for (i = 0; i < p_leafPre.LeafCount; ++i)
	{
		
		j = 3 * m_Pre.m_LeafPointCount[i];
		m_Pre.m_PreM[i].initial(j, j);
	}

	// Ԥ


	/*long a1, a2;
	long rownum, rowstart, beginID, tempa2,flag;*/

	//#pragma omp parallel for private (a1, a2,rownum, rowstart, beginID,tempa2,flag) num_threads(16)


	#pragma omp parallel for num_threads(thread_num)//thread_num
	for (long p = 0; p < m_Pre.m_LeafCount; ++p)
	{
		long a1, a2;
		long rownum, rowstart, beginID, tempa2, flag;
		Wmatrix subPre(3, 3);
		//printf("threadnum = %d\n", omp_get_num_threads());

		beginID = m_Pre.m_LeafBeginID[p];

		for (long l = 0; l < m_Pre.m_LeafPointCount[p]; ++l)
		{
			//a1 = 3 * m_Pre.m_RePID[beginID + j];   

			//for (k = 0; k < m_Pre.m_LeafPointCount[i]; ++k)
			//{
			//	a2 = 3 * m_Pre.m_RePID[beginID + k];
			//	subPre(1, 1) = A(a1 + 1, a2 + 1);
			//	subPre(1, 2) = A(a1 + 1, a2 + 2);
			//	subPre(1, 3) = A(a1 + 1, a2 + 3);
			//	subPre(2, 1) = A(a1 + 2, a2 + 1);
			//	subPre(2, 2) = A(a1 + 2, a2 + 2);
			//	subPre(2, 3) = A(a1 + 2, a2 + 3);
			//	subPre(3, 1) = A(a1 + 3, a2 + 1);
			//	subPre(3, 2) = A(a1 + 3, a2 + 2);
			//	subPre(3, 3) = A(a1 + 3, a2 + 3);
			//	m_Pre.m_PreM[i].assign(subPre, 3 * j + 1, 3 * k + 1, 3, 3);
			//}
			a1 = m_Pre.m_RePID[beginID + l]; //

			for (long k = 0; k < m_Pre.m_LeafPointCount[p]; ++k)
			{
				a2 = m_Pre.m_RePID[beginID + k];// 

				rownum = M.row_ptr[a1 + 1] - M.row_ptr[a1];
				rowstart = M.row_ptr[a1];
				flag = 0;
				for (long m = 0; m < rownum; ++m) {
					if (a2 == M.col_index[rowstart + m]){
						tempa2 = rowstart + m;
						flag = 1;
						break;
					}
				}
				if (flag == 1) {
					subPre(1, 1) = M.value[tempa2 * 9];
					subPre(1, 2) = M.value[tempa2 * 9 + 1];
					subPre(1, 3) = M.value[tempa2 * 9+2];
					subPre(2, 1) = M.value[tempa2 * 9+3];
					subPre(2, 2) = M.value[tempa2 * 9+4];
					subPre(2, 3) = M.value[tempa2 * 9+5];
					subPre(3, 1) = M.value[tempa2 * 9+6];
					subPre(3, 2) = M.value[tempa2 * 9+7];
					subPre(3, 3) = M.value[tempa2 * 9+8];
				}
				else {
						subPre(1, 1) = 0;
						subPre(1, 2) = 0;
						subPre(1, 3) = 0;
						subPre(2, 1) = 0;
						subPre(2, 2) = 0;
						subPre(2, 3) = 0;
						subPre(3, 1) = 0;
						subPre(3, 2) = 0;
						subPre(3, 3) = 0;
				}
				m_Pre.m_PreM[p].assign(subPre, 3 * l + 1, 3 * k + 1, 3, 3);
			}
		}
		
	}

	// ԤԤ
	for (i = 0; i < m_Pre.m_LeafCount; ++i)
		inv_mat(m_Pre.m_PreM[i], m_Pre.m_PreM[i].m);

	//release PreTree and p_leafPre
	DeleteTree(m_treePre);
	DeleteLeafPointer(p_leafPre);

	return 1;
}

int GMRES_M_X(PreConditioner& m_Pre, Wvector& X, Wvector& B)
{
	int m, n;
	long i, j, k, tj, tk;
	long a1, a2;

	B = 0.0;

	for (i = 0; i < m_Pre.m_LeafCount; ++i)
	{
		for (j = 0; j < m_Pre.m_LeafPointCount[i]; ++j)
		{
			tj = m_Pre.IID[j];
			long outPID = m_Pre.m_OutputPID ? m_Pre.m_OutputPID[m_Pre.m_LeafBeginID[i] + j] : m_Pre.m_RePID[m_Pre.m_LeafBeginID[i] + j];
			a1 = m_Pre.IID[outPID];
			for (k = 0; k < m_Pre.m_LeafPointCount[i]; ++k)
			{
				tk = m_Pre.IID[k];
				long inPID = m_Pre.m_InputPID ? m_Pre.m_InputPID[m_Pre.m_LeafBeginID[i] + k] : m_Pre.m_RePID[m_Pre.m_LeafBeginID[i] + k];
				a2 = m_Pre.IID[inPID];
				for (m = 1; m <= 3; ++m)
				{
					for (n = 1; n <= 3; ++n)
						B[a1 + m] += m_Pre.m_PreM[i](tj + m, tk + n)*X[a2 + n];
				}
			}
		}
	}

	return 1;
}

//precondictioned GMRES
int gmres(DSquareElement* m_DSE, Wmatrix& A, Wvector& b, PreConditioner& m_Pre, long int restart, double tol, long int maxit, Wvector& x0, Wvector& x, int& flag, int* iter)
{
	//A:             ϵ
	//b:             
	//restart:       ѭ
	//tol:           
	//maxit:         ѭ
	//M:             Ԥ
	//x0:            ʼ
	//x:             
	//flag:          Ƿȷ013
	//iter:              iter[0]:inner      iter[1]:outer
	long n = A.m;
	long int i, j, k, l;                   //used for loop
	double tolb;                        //tol*norm(b)
	double n2b;                         //b2
	Wvector r(n);                       //r=b-A*x
	Wvector vh(n);                      //൱r0
	Wvector u(n);
	Wvector u2(n);
	Wvector q(restart);
	double normr;                       //r2
	int stag;                           //Ƿж
	double phibar;
	double rt, c, s, temp;

	Wmatrix V(n, restart + 1);             //Arnoldi 
	Wvector h(restart + 1);               //Hessenbergһ st A*V = V*H
	Wmatrix QT(restart + 1, restart + 1);    //ת st QT*H = R
	Wmatrix R(restart, restart);         //תǾ st H = Q*R
	Wvector f(restart);                 //y = R\f => x = x0 + V*y
	Wmatrix W(n, restart);               //W = V*inv(R)

	n2b = vector_norm2(b, n);
	if (n2b == 0.0)
	{
		x = 0.0;
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	x = x0;
	flag = 1;
	tolb = tol * n2b;

	mat_vec_product(A, x, r);

	//test output for AX
//	FILE* whtB=fopen("GMRES_AX.txt","w");
//	for(i=1;i<=r.n;++i)
//		fprintf(whtB,"%4ld  %15.10lf\n",i,r[i]);
//	fclose(whtB);

	for (i = 1; i <= n; i++)
	{
		r[i] = b[i] - r[i];
	}//end for ii

	normr = vector_norm2(r, n);
	if (normr <= tolb)
	{
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	stag = 0;
	for (i = 1; i <= maxit; i++)  //ѭ
	{

		QT = 0.0;
		R = 0.0;

		// Ԥ
		GMRES_M_X(m_Pre, r, vh);

		j = 0;
		h[1] = vector_norm2(vh, n);
		for (j = 1; j <= n; j++)  //V(:,1) = vh / h(1)
			V(j, 1) = vh[j] / h[1];

		QT(1, 1) = 1.0;

		phibar = h[1];

		for (j = 1; j <= restart; j++)        //ѭ
		{
			printf("GMRES outer=%ld inner=%ld\n", i, j);

			for (k = 1; k <= n; k++)
				u[k] = V(k, j);
			mat_vec_product(A, u, u2);

			// Ԥ
			GMRES_M_X(m_Pre, u2, u);

			for (k = 1; k <= j; k++)
			{
				h[k] = 0.0;
				for (l = 1; l <= n; l++)//h(k) = V(:,k)' * u;
					h[k] += V(l, k)*u[l];
				for (l = 1; l <= n; l++)//u = u - h(k) * V(:,k)
					u[l] -= h[k] * V(l, k);
			}

			h[j + 1] = vector_norm2(u, n);

			if (h[j + 1] <= 1.0e-12)
			{
				for (k = 1; k <= n; k++)
					V(k, j + 1) = 0.0;
			}
			else
			{
				for (k = 1; k <= n; k++)//V(:,j+1) = u / h(j+1);
					V(k, j + 1) = u[k] / h[j + 1];
			}

			for (k = 1; k <= j; k++)
			{
				R(k, j) = 0.0;
				for (l = 1; l <= j; l++)
					R(k, j) += QT(k, l)*h[l];
			}//R(1:j,j) = QT(1:j,1:j) * h(1:j)

			rt = R(j, j);

			if (h[j + 1] == 0.0)
			{
				c = 1.0;
				s = 0.0;
			}
			else if (h[j + 1] * h[j + 1] > rt*rt)
			{
				temp = rt / h[j + 1];
				s = 1.0 / sqrt(1.0 + fabs(temp)*fabs(temp));
				c = -temp * s;
			}
			else
			{
				temp = h[j + 1] / rt;
				c = 1.0 / sqrt(1.0 + fabs(temp)*fabs(temp));
				s = -temp * c;
			}

			R(j, j) = c * rt - s * h[j + 1];

			for (k = 1; k <= j; k++)
				q[k] = QT(j, k);

			for (k = 1; k <= j; k++)
				QT(j, k) = c * q[k];

			for (k = 1; k <= j; k++)
				QT(j + 1, k) = s * q[k];

			QT(j, j + 1) = -s;
			QT(j + 1, j + 1) = c;
			f[j] = c * phibar;
			phibar = s * phibar;

			if (j < restart)
			{
				if (f[j] == 0.0)
					stag = 1;
				if (j == 1)
				{
					for (k = 1; k <= n; k++)
						W(k, j) = V(k, j) / R(j, j);
				}
				else
				{
					for (k = 1; k <= n; k++)
					{
						W(k, j) = V(k, j) / R(j, j);
						for (l = 1; l <= j - 1; l++)
							W(k, j) -= W(k, l)*R(l, j) / R(j, j);
					}
				}
				for (k = 1; k <= n; k++)
					x[k] += f[j] * W(k, j);
			}
			else
			{
				inv_mat(R, j);
				mat_vec_product(R, f, q);
				mat_vec_product(V, q, x);

				x += x0;
			}

			mat_vec_product(A, x, vh);
			for (k = 1; k <= n; k++)
				vh[k] = b[k] - vh[k];

			normr = vector_norm2(vh, n);

			if (normr <= tolb)
			{
				flag = 0;
				iter[0] = i;
				iter[1] = j;
				return 1;
			}

			if (stag == 1)
			{
				flag = 3;
				iter[0] = i;
				iter[1] = j;
				return -1;
			}

		}//end for j

		if (flag == 1)
		{
			x0 = x;
			mat_vec_product(A, x0, r);
			for (k = 1; k <= n; k++)
				r[k] = b[k] - r[k];
			iter[0] = i;
			iter[1] = j;
		}

	}//end for i


	return 1;
}

int gmres(DSquareElement* m_DSE, CSRMat& A, Wvector& b, PreConditioner& m_Pre, long int restart, double tol, long int maxit, Wvector& x0, Wvector& x, int& flag, int* iter)
{
	//A:             ϵ
//b:             
//restart:       ѭ
//tol:           
//maxit:         ѭ
//M:             Ԥ
//x0:            ʼ
//x:             
//flag:          Ƿȷ013
//iter:              iter[0]:inner      iter[1]:outer
	long n = A.n;
	long int i, j, k, l;                   //used for loop
	double tolb;                        //tol*norm(b)
	double n2b;                         //b2
	Wvector r(n);                       //r=b-A*x
	Wvector vh(n);                      //൱r0
	Wvector u(n);
	Wvector u2(n);
	Wvector q(restart);
	double normr;                       //r2
	int stag;                           //Ƿж
	double phibar;
	double rt, c, s, temp;

	Wmatrix V(n, restart + 1);             //Arnoldi 
	Wvector h(restart + 1);               //Hessenbergһ st A*V = V*H
	Wmatrix QT(restart + 1, restart + 1);    //ת st QT*H = R
	Wmatrix R(restart, restart);         //תǾ st H = Q*R
	Wvector f(restart);                 //y = R\f => x = x0 + V*y
	Wmatrix W(n, restart);               //W = V*inv(R)

	n2b = vector_norm2(b, n);
	if (n2b == 0.0)
	{
		x = 0.0;
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	x = x0;
	flag = 1;
	tolb = tol * n2b;

	//mat_vec_product(A, x, r);
	SparseMul(A, x, r);
	//test output for AX
//	FILE* whtB=fopen("GMRES_AX.txt","w");
//	for(i=1;i<=r.n;++i)
//		fprintf(whtB,"%4ld  %15.10lf\n",i,r[i]);
//	fclose(whtB);

	for (i = 1; i <= n; i++)
	{
		r[i] = b[i] - r[i];
	}//end for ii

	normr = vector_norm2(r, n);
	if (normr <= tolb)
	{
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	stag = 0;
	for (i = 1; i <= maxit; i++)  //ѭ
	{

		QT = 0.0;
		R = 0.0;

		// Ԥ
		GMRES_M_X(m_Pre, r, vh);

		j = 0;
		h[1] = vector_norm2(vh, n);
		for (j = 1; j <= n; j++)  //V(:,1) = vh / h(1)
			V(j, 1) = vh[j] / h[1];

		QT(1, 1) = 1.0;

		phibar = h[1];

		for (j = 1; j <= restart; j++)        //ѭ
		{
			printf("GMRES outer=%ld inner=%ld\n", i, j);

			for (k = 1; k <= n; k++)
				u[k] = V(k, j);
			//mat_vec_product(A, u, u2);
			SparseMul(A, u, u2);

			// Ԥ
			GMRES_M_X(m_Pre, u2, u);

			for (k = 1; k <= j; k++)
			{
				h[k] = 0.0;
				for (l = 1; l <= n; l++)//h(k) = V(:,k)' * u;
					h[k] += V(l, k) * u[l];
				for (l = 1; l <= n; l++)//u = u - h(k) * V(:,k)
					u[l] -= h[k] * V(l, k);
			}

			h[j + 1] = vector_norm2(u, n);

			if (h[j + 1] <= 1.0e-12)
			{
				for (k = 1; k <= n; k++)
					V(k, j + 1) = 0.0;
			}
			else
			{
				for (k = 1; k <= n; k++)//V(:,j+1) = u / h(j+1);
					V(k, j + 1) = u[k] / h[j + 1];
			}

			for (k = 1; k <= j; k++)
			{
				R(k, j) = 0.0;
				for (l = 1; l <= j; l++)
					R(k, j) += QT(k, l) * h[l];
			}//R(1:j,j) = QT(1:j,1:j) * h(1:j)

			rt = R(j, j);

			if (h[j + 1] == 0.0)
			{
				c = 1.0;
				s = 0.0;
			}
			else if (h[j + 1] * h[j + 1] > rt * rt)
			{
				temp = rt / h[j + 1];
				s = 1.0 / sqrt(1.0 + fabs(temp) * fabs(temp));
				c = -temp * s;
			}
			else
			{
				temp = h[j + 1] / rt;
				c = 1.0 / sqrt(1.0 + fabs(temp) * fabs(temp));
				s = -temp * c;
			}

			R(j, j) = c * rt - s * h[j + 1];

			for (k = 1; k <= j; k++)
				q[k] = QT(j, k);

			for (k = 1; k <= j; k++)
				QT(j, k) = c * q[k];

			for (k = 1; k <= j; k++)
				QT(j + 1, k) = s * q[k];

			QT(j, j + 1) = -s;
			QT(j + 1, j + 1) = c;
			f[j] = c * phibar;
			phibar = s * phibar;

			if (j < restart)
			{
				if (f[j] == 0.0)
					stag = 1;
				if (j == 1)
				{
					for (k = 1; k <= n; k++)
						W(k, j) = V(k, j) / R(j, j);
				}
				else
				{
					for (k = 1; k <= n; k++)
					{
						W(k, j) = V(k, j) / R(j, j);
						for (l = 1; l <= j - 1; l++)
							W(k, j) -= W(k, l) * R(l, j) / R(j, j);
					}
				}
				for (k = 1; k <= n; k++)
					x[k] += f[j] * W(k, j);
			}
			else
			{
				inv_mat(R, j);
				mat_vec_product(R, f, q);
				mat_vec_product(V, q, x);

				x += x0;
			}

			//mat_vec_product(A, x, vh);
			SparseMul(A, x, vh);
			for (k = 1; k <= n; k++)
				vh[k] = b[k] - vh[k];

			normr = vector_norm2(vh, n);

			if (normr <= tolb)
			{
				flag = 0;
				iter[0] = i;
				iter[1] = j;
				return 1;
			}

			if (stag == 1)
			{
				flag = 3;
				iter[0] = i;
				iter[1] = j;
				return -1;
			}

		}//end for j

		if (flag == 1)
		{
			x0 = x;
			//mat_vec_product(A, x0, r);
			SparseMul(A, x0, r);
			for (k = 1; k <= n; k++)
				r[k] = b[k] - r[k];
			iter[0] = i;
			iter[1] = j;
		}

	}//end for i


	return 1;
}

int gmres(DSquareElement* m_DSE, CCSRMat& A, Wvector& b, PreConditioner& m_Pre, long int restart, double tol, long int maxit, Wvector& x0, Wvector& x, int& flag, int* iter)
{
//A:             ϵ
//b:             
//restart:       ѭ
//tol:           
//maxit:         ѭ
//M:             Ԥ
//x0:            ʼ
//x:             
//flag:          Ƿȷ013
//iter:              iter[0]:inner      iter[1]:outer
	long n = A.n * 3;
	long int i, j, k, l;                   //used for loop
	double tolb;                        //tol*norm(b)
	double n2b;                         //b2
	Wvector r(n);                       //r=b-A*x
	Wvector vh(n);                      //൱r0
	Wvector u(n);
	Wvector u2(n);
	Wvector q(restart);
	double normr;                       //r2
	int stag;                           //Ƿж
	double phibar;
	double rt, c, s, temp;

	Wmatrix V(n, restart + 1);             //Arnoldi 
	Wvector h(restart + 1);               //Hessenbergһ st A*V = V*H
	Wmatrix QT(restart + 1, restart + 1);    //ת st QT*H = R
	Wmatrix R(restart, restart);         //תǾ st H = Q*R
	Wvector f(restart);                 //y = R\f => x = x0 + V*y
	Wmatrix W(n, restart);               //W = V*inv(R)


	n2b = vector_norm2(b, n);
	if (n2b == 0.0)
	{
		x = 0.0;
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	x = x0;
	flag = 1;
	tolb = tol * n2b;

	//mat_vec_product(A, x, r);
	SparseMul(A, x, r);
	//test output for AX
//	FILE* whtB=fopen("GMRES_AX.txt","w");
//	for(i=1;i<=r.n;++i)
//		fprintf(whtB,"%4ld  %15.10lf\n",i,r[i]);
//	fclose(whtB);

	for (i = 1; i <= n; i++)
	{
		r[i] = b[i] - r[i];
	}//end for ii

	normr = vector_norm2(r, n);
	if (normr <= tolb)
	{
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	stag = 0;
	for (i = 1; i <= maxit; i++)  //ѭ
	{

		QT = 0.0;
		R = 0.0;

		// Ԥ
		GMRES_M_X(m_Pre, r, vh);

		j = 0;
		h[1] = vector_norm2(vh, n);
		for (j = 1; j <= n; j++)  //V(:,1) = vh / h(1)
			V(j, 1) = vh[j] / h[1];

		QT(1, 1) = 1.0;

		phibar = h[1];

		for (j = 1; j <= restart; j++)        //ѭ
		{
			printf("GMRES outer=%ld inner=%ld\n", i, j);

			for (k = 1; k <= n; k++)
				u[k] = V(k, j);
			//mat_vec_product(A, u, u2);
			SparseMul(A, u, u2);

			// Ԥ
			GMRES_M_X(m_Pre, u2, u);

			for (k = 1; k <= j; k++)
			{
				h[k] = 0.0;
				for (l = 1; l <= n; l++)//h(k) = V(:,k)' * u;
					h[k] += V(l, k) * u[l];
				for (l = 1; l <= n; l++)//u = u - h(k) * V(:,k)
					u[l] -= h[k] * V(l, k);
			}

			h[j + 1] = vector_norm2(u, n);

			if (h[j + 1] <= 1.0e-12)
			{
				for (k = 1; k <= n; k++)
					V(k, j + 1) = 0.0;
			}
			else
			{
				for (k = 1; k <= n; k++)//V(:,j+1) = u / h(j+1);
					V(k, j + 1) = u[k] / h[j + 1];
			}

			for (k = 1; k <= j; k++)
			{
				R(k, j) = 0.0;
				for (l = 1; l <= j; l++)
					R(k, j) += QT(k, l) * h[l];
			}//R(1:j,j) = QT(1:j,1:j) * h(1:j)

			rt = R(j, j);

			if (h[j + 1] == 0.0)
			{
				c = 1.0;
				s = 0.0;
			}
			else if (h[j + 1] * h[j + 1] > rt * rt)
			{
				temp = rt / h[j + 1];
				s = 1.0 / sqrt(1.0 + fabs(temp) * fabs(temp));
				c = -temp * s;
			}
			else
			{
				temp = h[j + 1] / rt;
				c = 1.0 / sqrt(1.0 + fabs(temp) * fabs(temp));
				s = -temp * c;
			}

			R(j, j) = c * rt - s * h[j + 1];

			for (k = 1; k <= j; k++)
				q[k] = QT(j, k);

			for (k = 1; k <= j; k++)
				QT(j, k) = c * q[k];

			for (k = 1; k <= j; k++)
				QT(j + 1, k) = s * q[k];

			QT(j, j + 1) = -s;
			QT(j + 1, j + 1) = c;
			f[j] = c * phibar;
			phibar = s * phibar;

			if (j < restart)
			{
				if (f[j] == 0.0)
					stag = 1;
				if (j == 1)
				{
					for (k = 1; k <= n; k++)
						W(k, j) = V(k, j) / R(j, j);
				}
				else
				{
					for (k = 1; k <= n; k++)
					{
						W(k, j) = V(k, j) / R(j, j);
						for (l = 1; l <= j - 1; l++)
							W(k, j) -= W(k, l) * R(l, j) / R(j, j);
					}
				}
				for (k = 1; k <= n; k++)
					x[k] += f[j] * W(k, j);
			}
			else
			{
				inv_mat(R, j);
				mat_vec_product(R, f, q);
				mat_vec_product(V, q, x);

				x += x0;
			}

			//mat_vec_product(A, x, vh);
			SparseMul(A, x, vh);
			for (k = 1; k <= n; k++)
				vh[k] = b[k] - vh[k];

			normr = vector_norm2(vh, n);

			if (normr <= tolb)
			{
				flag = 0;
				iter[0] = i;
				iter[1] = j;
				return 1;
			}

			if (stag == 1)
			{
				flag = 3;
				iter[0] = i;
				iter[1] = j;
				return -1;
			}

		}//end for j

		if (flag == 1)
		{
			x0 = x;
			//mat_vec_product(A, x0, r);
			SparseMul(A, x0, r);
			for (k = 1; k <= n; k++)
				r[k] = b[k] - r[k];
			iter[0] = i;
			iter[1] = j;
		}

	}//end for i


	return 1;
}


int DynaGMRESSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// Ԥ
	PreConditioner m_Pre;

	int flag;                                      //ź
	int iter[2];

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// ȷķӦ
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf(": %ld\n", MaxN);



	pthread_t* threads = 0;

	// 1. ʼϵ
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, ΪûMT2(0)MT2(NStep)ҲڳʼλΪ0û
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)MT1[0],Ҳδ֪ MG(0)MG1[0]Ҳ֪
	// MT(1~NStep-1)MT1 + MT2MG(1~NStep-1)MG1 + MG2
	// MG(NStep)MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}

	// 2. 0ʱ̵֪δ֪˳ΪֲµλƺλΪ0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. n = 1ʱ̵ϵõ n = 1 ʱ̵Ľ

	// TEST
	printf("GMRES Solver: Step 1.\n");

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum);

	//Ԥ
	GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MT[0]);

	x0 = 0.0;

	// TEST
	printf("GMRES Solver: Step 1. Calculate T2ij1 and G2ij1...\n");

	GetMatrixDynaT2ij(m_DSE, MT[1], NodeNum, EleNum, 1.0, dt);

	// B = MG_1(0) * bd[1].luMG_1(0)ʵ֪bd[1].lu1ʱ̵֪
	mat_vec_product(MG[0], bd[1].lu, B);

	printf("GMRES Solver: Step 1. Calculate Result...\n");

	//  n = 1ʱ̵δ֪bd[1].lt
	gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[1].lt, flag, iter);

	// TEST
	printf("GMRES Solver: Step 1. Convert...\n");

	//  n = 1ʱ̵δ֪/֪תΪλƺ
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4.  n = 2 ~ NStepĽ

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("GMRES Solver: Step %ld.\n", i);

		if (i <= MaxN)
		{
			// ʹMG[i]ݴMT_1(i-1)MG_1(i-1)ֱۼӵԭMT_2(i-1)MG_2(i-1)MT(i-1)MG(i-1)
			GetMatrixDynaUij(m_DSE, MG[i - 1], NodeNum, EleNum, (i - 1)*1.0, dt);

			GetMatrixDynaT1ij(m_DSE, MG[i], NodeNum, EleNum, (i - 1)*1.0, dt);
			MT[i - 1] += MG[i];

			// MT_2(i)(i < NStep)MT(i)
			if (i < NStep)
				GetMatrixDynaT2ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt);
		}

		// B = MG_1(0) * bd[i].luMG_1(0)ʵ֪bd[i].luiʱ̵֪
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i-1
		for (j = 1; j <= i - 1; ++j)
		{
			if (j <= MaxN)
			{
				// CB1 = [MG_1(j) + MG_2(j)] * bd[i-j].ltbd[i-j].ltȷָi-jʱ̵
				mat_vec_product(MG[j], bd[i - j].lt, CB1);

				B += CB1;

				// CB1 = [MT_1(j) + MT_2(j)] * bd[i-j].lubd[i-j].luȷָi-jʱ̵λ
				mat_vec_product(MT[j], bd[i - j].lu, CB1);
				B -= CB1;
			}
		}

		//  n = iʱ̵δ֪bd[i].lt
		gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[i].lt, flag, iter);

		// TEST
//		bd[i].lt.OutputData("R",i);

		//  n = iʱ̵δ֪/֪תΪλƺ
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}

	// ͷſռ

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	char* CurrPath = _getcwd(NULL, 0);//صǰ·
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//صǰļ
	Path = StrCurrPath.substr(0, found + 1);//ϼ·

	string OutFile = "output\\";
	if (GMRES_OUTPUT_SUBDIR[0] != '\0') {
		OutFile += string(GMRES_OUTPUT_SUBDIR) + "\\";
		_mkdir("output");
		string outDir = string("output\\") + GMRES_OUTPUT_SUBDIR;
		_mkdir(outDir.c_str());
	}
	string TempPath = Path + OutFile;



	string Temp = "GMRESTime.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "ʱ:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

int DynaGMRESSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// Ԥ
	PreConditioner m_Pre;

	int flag;                                      //ź
	int iter[2];

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// ȷķӦ
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf(": %ld\n", MaxN);



	pthread_t* threads = 0;

	// 1. ʼϵ
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, ΪûMT2(0)MT2(NStep)ҲڳʼλΪ0û
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)MT1[0],Ҳδ֪ MG(0)MG1[0]Ҳ֪
	// MT(1~NStep-1)MT1 + MT2MG(1~NStep-1)MG1 + MG2
	// MG(NStep)MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}

	// 2. 

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum);

	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		GetMatrixDynaT1ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt);
		GetMatrixDynaT2ij(m_DSE, MG[MaxN], NodeNum, EleNum, i*1.0, dt);
		MT[i] += MG[MaxN];
		GetMatrixDynaUij(m_DSE, MG[i], NodeNum, EleNum, i*1.0, dt);
	}

	// 3. 0ʱ̵֪δ֪˳ΪֲµλƺλΪ0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. n = 1ʱ̵ϵõ n = 1 ʱ̵Ľ

	// TEST
	printf("GMRES Solver: Step 1.\n");

	//Ԥ
	GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MT[0]);

	x0 = 0.0;

	// B = MG_1(0) * bd[1].luMG_1(0)ʵ֪bd[1].lu1ʱ̵֪
	mat_vec_product(MG[0], bd[1].lu, B);

	//  n = 1ʱ̵δ֪bd[1].lt
	gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[1].lt, flag, iter);

	//  n = 1ʱ̵δ֪/֪תΪλƺ
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4.  n = 2 ~ NStepĽ

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("GMRES Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].luMG_1(0)ʵ֪bd[i].luiʱ̵֪
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;

				mat_vec_product(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;

				mat_vec_product(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		//  n = iʱ̵δ֪bd[i].lt
		gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[i].lt, flag, iter);

		//  n = iʱ̵δ֪/֪תΪλƺ
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}

	// ͷſռ

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	char* CurrPath = _getcwd(NULL, 0);//صǰ·
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//صǰļ
	Path = StrCurrPath.substr(0, found + 1);//ϼ·

	string OutFile = "output\\";
	if (GMRES_OUTPUT_SUBDIR[0] != '\0') {
		OutFile += string(GMRES_OUTPUT_SUBDIR) + "\\";
		_mkdir("output");
		string outDir = string("output\\") + GMRES_OUTPUT_SUBDIR;
		_mkdir(outDir.c_str());
	}
	string TempPath = Path + OutFile;



	string Temp = "GMRESTime.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "ʱ:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}
//޵Ԫ
int DynaGMRESSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// Ԥ
	PreConditioner m_Pre;

	int flag;                                      //ź
	int iter[2];

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// ȷķӦ
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf(": %ld\n", MaxN);



	pthread_t* threads = 0;

	// 1. ʼϵ
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, ΪûMT2(0)MT2(NStep)ҲڳʼλΪ0û
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)MT1[0],Ҳδ֪ MG(0)MG1[0]Ҳ֪
	// MT(1~NStep-1)MT1 + MT2MG(1~NStep-1)MG1 + MG2
	// MG(NStep)MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}

	// 2. 0ʱ̵֪δ֪˳ΪֲµλƺλΪ0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. n = 1ʱ̵ϵõ n = 1 ʱ̵Ľ

	// TEST
	printf("GMRES Solver: Step 1.\n");

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	//GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	//GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum);

	//Ԥ
	GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MT[0]);

	x0 = 0.0;

	// TEST
	printf("GMRES Solver: Step 1. Calculate T2ij1 and G2ij1...\n");

	//GetMatrixDynaT2ij(m_DSE, MT[1], NodeNum, EleNum, 1.0, dt);
	GetMatrixDynaT2ij(m_DSE, MT[1], NodeNum, EleNum, 1.0, dt, m_InfElePID, m_ElePID);

	// B = MG_1(0) * bd[1].luMG_1(0)ʵ֪bd[1].lu1ʱ̵֪
	mat_vec_product(MG[0], bd[1].lu, B);

	printf("GMRES Solver: Step 1. Calculate Result...\n");

	//  n = 1ʱ̵δ֪bd[1].lt
	gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[1].lt, flag, iter);

	// TEST
	printf("GMRES Solver: Step 1. Convert...\n");

	//  n = 1ʱ̵δ֪/֪תΪλƺ
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4.  n = 2 ~ NStepĽ

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("GMRES Solver: Step %ld.\n", i);

		if (i <= MaxN)
		{
			// ʹMG[i]ݴMT_1(i-1)MG_1(i-1)ֱۼӵԭMT_2(i-1)MG_2(i-1)MT(i-1)MG(i-1)
			GetMatrixDynaUij(m_DSE, MG[i - 1], NodeNum, EleNum, (i - 1)*1.0, dt);

			GetMatrixDynaT1ij(m_DSE, MG[i], NodeNum, EleNum, (i - 1)*1.0, dt, m_InfElePID, m_ElePID);
			MT[i - 1] += MG[i];

			// MT_2(i)(i < NStep)MT(i)
			if (i < NStep)
				GetMatrixDynaT2ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt, m_InfElePID, m_ElePID);
		}

		// B = MG_1(0) * bd[i].luMG_1(0)ʵ֪bd[i].luiʱ̵֪
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i-1
		for (j = 1; j <= i - 1; ++j)
		{
			if (j <= MaxN)
			{
				// CB1 = [MG_1(j) + MG_2(j)] * bd[i-j].ltbd[i-j].ltȷָi-jʱ̵
				mat_vec_product(MG[j], bd[i - j].lt, CB1);

				B += CB1;

				// CB1 = [MT_1(j) + MT_2(j)] * bd[i-j].lubd[i-j].luȷָi-jʱ̵λ
				mat_vec_product(MT[j], bd[i - j].lu, CB1);
				B -= CB1;
			}
		}

		//  n = iʱ̵δ֪bd[i].lt
		gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[i].lt, flag, iter);

		// TEST
//		bd[i].lt.OutputData("R",i);

		//  n = iʱ̵δ֪/֪תΪλƺ
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}

	// ͷſռ

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	char* CurrPath = _getcwd(NULL, 0);//صǰ·
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//صǰļ
	Path = StrCurrPath.substr(0, found + 1);//ϼ·

	string OutFile = "output\\";
	if (GMRES_OUTPUT_SUBDIR[0] != '\0') {
		OutFile += string(GMRES_OUTPUT_SUBDIR) + "\\";
		_mkdir("output");
		string outDir = string("output\\") + GMRES_OUTPUT_SUBDIR;
		_mkdir(outDir.c_str());
	}
	string TempPath = Path + OutFile;



	string Temp = "GMRESTime.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "ʱ:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}
//޵Ԫ
int DynaGMRESSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// Ԥ
	PreConditioner m_Pre;

	int flag;                                      //ź
	int iter[2];

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// ȷķӦ
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf(": %ld\n", MaxN);



	pthread_t* threads = 0;

	// 1. ʼϵ
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, ΪûMT2(0)MT2(NStep)ҲڳʼλΪ0û
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)MT1[0],Ҳδ֪ MG(0)MG1[0]Ҳ֪
	// MT(1~NStep-1)MT1 + MT2MG(1~NStep-1)MG1 + MG2
	// MG(NStep)MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}

	// 2. 

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum, m_InfElePID, m_ElePID);

	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		GetMatrixDynaT1ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt, m_InfElePID, m_ElePID);
		GetMatrixDynaT2ij(m_DSE, MG[MaxN], NodeNum, EleNum, i*1.0, dt, m_InfElePID, m_ElePID);
		MT[i] += MG[MaxN];
		GetMatrixDynaUij(m_DSE, MG[i], NodeNum, EleNum, i*1.0, dt);
	}

	// 3. 0ʱ̵֪δ֪˳ΪֲµλƺλΪ0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. n = 1ʱ̵ϵõ n = 1 ʱ̵Ľ

	m_Counter.EndCount();
	double DifTime1 = m_Counter.TimeDiff();

	m_Counter.StartCount();




	// TEST
	printf("GMRES Solver: Step 1.\n");

	//Ԥ
	GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MT[0]);

	x0 = 0.0;

	// B = MG_1(0) * bd[1].luMG_1(0)ʵ֪bd[1].lu1ʱ̵֪
	mat_vec_product(MG[0], bd[1].lu, B);

	//  n = 1ʱ̵δ֪bd[1].lt
	gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[1].lt, flag, iter);

	//  n = 1ʱ̵δ֪/֪תΪλƺ
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4.  n = 2 ~ NStepĽ

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("GMRES Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].luMG_1(0)ʵ֪bd[i].luiʱ̵֪
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;

				mat_vec_product(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;

				mat_vec_product(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		//  n = iʱ̵δ֪bd[i].lt
		gmres(m_DSE, MT[0], B, m_Pre, iterations, error, 1, x0, bd[i].lt, flag, iter);

		//  n = iʱ̵δ֪/֪תΪλƺ
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}



	showMemoryInfo();
	// ͷſռ

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	char* CurrPath = _getcwd(NULL, 0);//صǰ·
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//صǰļ
	Path = StrCurrPath.substr(0, found + 1);//ϼ·

	string OutFile = "output\\";
	if (GMRES_OUTPUT_SUBDIR[0] != '\0') {
		OutFile += string(GMRES_OUTPUT_SUBDIR) + "\\";
		_mkdir("output");
		string outDir = string("output\\") + GMRES_OUTPUT_SUBDIR;
		_mkdir(outDir.c_str());
	}
	string TempPath = Path + OutFile;



	string Temp = "GMRESTime.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "ʱ:\n");
	fprintf(pf, "װʱ:%lf\n", DifTime1);
	fprintf(pf, "ʱ:%lf\n", m_Counter.TimeDiff());
	fclose(pf);

	//FILE* SparseM;
	//double nnum = 9 * NodeNum * NodeNum;
	//Temp = "SparseM.txt";
	//Path = TempPath + Temp;
	//OutPath = Path.c_str();
	//fopen_s(&SparseM, OutPath, "w");
	//fprintf_s(SparseM, "MT\n");
	//for (int i = 0; i <= MaxN; i++) {
	//	fprintf_s(SparseM, "%.8lf", MTS[i].row_ptr[(MTS[i].n)] / nnum);
	//	fprintf_s(SparseM, "\n");
	//}
	//fprintf_s(SparseM, "MG\n");
	//for (int i = 0; i <= MaxN; i++) {
	//	fprintf_s(SparseM, "%.8lf", MGS[i].row_ptr[(MGS[i].n)] / nnum);
	//	fprintf_s(SparseM, "\n");
	//}
	//fclose(SparseM);




	return 1;
}

int DynaGMRESSolverNewCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;
	CSRMat* MGS;
	CSRMat* MTS;
	WmatrixCSR Mtemp;//ʱ
	WvectorCSR Btemp;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// Ԥ
	PreConditioner m_Pre;

	int flag;                                      //ź
	int iter[2];

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// ȷķӦ
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf(": %ld\n", MaxN);



	pthread_t* threads = 0;

	// 1. ʼϵ
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, ΪûMT2(0)MT2(NStep)ҲڳʼλΪ0û
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)MT1[0],Ҳδ֪ MG(0)MG1[0]Ҳ֪
	// MT(1~NStep-1)MT1 + MT2MG(1~NStep-1)MG1 + MG2
	// MG(NStep)MG2[NStep]
	MT = new Wmatrix[1];
	//MT = new Wmatrix[NStep + 1];
	//MG = new Wmatrix[NStep + 1];
	MTS = new CSRMat[NStep+1];
	MGS = new CSRMat[NStep+1];


	MT[0].initial(3 * NodeNum, 3 * NodeNum);
	Mtemp.initial(3 * NodeNum, 3 * NodeNum);//ʱʼ
	Mtemp.reinitial();//ֵ0
	Btemp.initial(3 * NodeNum);//ʱʼ



	// 2. 

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);

	GetMatrixDynaUnknown0ij(m_DSE, Mtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);//¸
	SparseConvert(Mtemp, MTS[0]);
	//GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	GetMatrixDynaKnown0ij(m_DSE, Mtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);
	SparseConvert(Mtemp, MGS[0]);




	//Ԥ
	GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MT[0]);
	delete[] MT;
	MT = 0;

	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		//GetMatrixDynaT1ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt);
		GetMatrixDynaT1ij(m_DSE, Mtemp, NodeNum, EleNum, i * 1.0, dt, m_InfElePID, m_ElePID);
		SparseConvert(Mtemp, MTS[i]);

		//GetMatrixDynaT2ij(m_DSE, MG[MaxN], NodeNum, EleNum, i*1.0, dt);
		GetMatrixDynaT2ij(m_DSE, Mtemp, NodeNum, EleNum, i * 1.0, dt, m_InfElePID, m_ElePID);
		//MT[i] += MG[MaxN];
		SparseAddition(Mtemp, MTS[i], Btemp);

		//GetMatrixDynaUij(m_DSE, MG[i], NodeNum, EleNum, i*1.0, dt);
		GetMatrixDynaUij(m_DSE, Mtemp, NodeNum, EleNum, i * 1.0, dt);
		SparseConvert(Mtemp, MGS[i]);
	}

	// 3. 0ʱ̵֪δ֪˳ΪֲµλƺλΪ0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. n = 1ʱ̵ϵõ n = 1 ʱ̵Ľ

	m_Counter.EndCount();
	double DifTime1 = m_Counter.TimeDiff();

	m_Counter.StartCount();

	// TEST
	printf("GMRES Solver: Step 1.\n");


	x0 = 0.0;

	// B = MG_1(0) * bd[1].luMG_1(0)ʵ֪bd[1].lu1ʱ̵֪
	//mat_vec_product(MG[0], bd[1].lu, B);
	SparseMul(MGS[0], bd[1].lu, B);

	//  n = 1ʱ̵δ֪bd[1].lt
	gmres(m_DSE, MTS[0], B, m_Pre, iterations, error, 1, x0, bd[1].lt, flag, iter);

	//  n = 1ʱ̵δ֪/֪תΪλƺ
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4.  n = 2 ~ NStepĽ

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("GMRES Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].luMG_1(0)ʵ֪bd[i].luiʱ̵֪
		//mat_vec_product(MG[0], bd[i].lu, B);
		SparseMul(MGS[0], bd[i].lu, B);
		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;

				SparseMul(MGS[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;

				SparseMul(MTS[j], CB1, CB2);
				B -= CB2;
			}
		}

		//  n = iʱ̵δ֪bd[i].lt
		gmres(m_DSE, MTS[0], B, m_Pre, iterations, error, 1, x0, bd[i].lt, flag, iter);

		//  n = iʱ̵δ֪/֪תΪλƺ
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}
	char* CurrPath = _getcwd(NULL, 0);//صǰ·
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//صǰļ
	Path = StrCurrPath.substr(0, found + 1);//ϼ·

	string OutFile = "output\\";
	if (GMRES_OUTPUT_SUBDIR[0] != '\0') {
		OutFile += string(GMRES_OUTPUT_SUBDIR) + "\\";
		_mkdir("output");
		string outDir = string("output\\") + GMRES_OUTPUT_SUBDIR;
		_mkdir(outDir.c_str());
	}
	string TempPath = Path + OutFile;
	//ϡȡ

	FILE* SparseM;
	double nnum = 9 * NodeNum * NodeNum;
	string Temp = "SparseM.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();
	fopen_s(&SparseM, OutPath, "w");
	fprintf_s(SparseM, "MT\n");
	for (int i = 0; i <= MaxN; i++) {
		fprintf_s(SparseM, "%.8lf", MTS[i].row_ptr[(MTS[i].n)] / nnum);
		fprintf_s(SparseM, "\n");
	}
	fprintf_s(SparseM, "MG\n");
	for (int i = 0; i <= MaxN; i++) {
		fprintf_s(SparseM, "%.8lf", MGS[i].row_ptr[(MGS[i].n)] / nnum);
		fprintf_s(SparseM, "\n");
	}
	fclose(SparseM);


	showMemoryInfo();

	//ͷſռ䡪


	delete[] MGS;
	MGS = 0;
	delete[] MTS;
	MTS = 0;
	Mtemp.finish();

	// End Counting...
	m_Counter.EndCount();
	




	Temp = "GMRESTime.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "ʱ:\n");
	fprintf(pf, "װʱ:%lf\n", DifTime1);
	fprintf(pf, "ʱ:%lf\n", m_Counter.TimeDiff());
	fclose(pf);





	return 1;
}

int DynaGMRESSolverNewCCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID)
{
	long i, j, p;

	//Wmatrix* MT;
	Wmatrix* MG;
	//SymCCSRMat* MGS;ȫ
	//CCSRMat* MTS;ȫ
	WmatrixCSR Mtemp;//ʱ
	//WmatrixCCSR MCCSRtemp;//ʱ
	//WmatrixSymCCSR MSymCtemp;//ʱ
	WvectorCSR Btemp;
	//WvectorCCSR BCCSRtemp; ȫ
	WvectorSymCCSR BSymCtemp;



	double dt = DSquareElement::m_DMat.Dt;
	dt_p = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// Ԥ
	PreConditioner m_Pre;

	int flag;                                      //ź
	int iter[2];

	Counter m_Counter;



	// ȷķӦ
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf(": %ld\n", MaxN);



	pthread_t* threads = 0;

	// 1. ʼϵ
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, ΪûMT2(0)MT2(NStep)ҲڳʼλΪ0û
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)MT1[0],Ҳδ֪ MG(0)MG1[0]Ҳ֪
	// MT(1~NStep-1)MT1 + MT2MG(1~NStep-1)MG1 + MG2
	// MG(NStep)MG2[NStep]
	MT = new Wmatrix[1];
	//MT = new Wmatrix[NStep + 1];
	//MG = new Wmatrix[NStep + 1];
	MTS = new CCSRMat[NStep + 1];
	MGS = new SymCCSRMat[NStep + 1];



	//FlagUij = new long** [MaxN + 1];
	//for (i = 0; i < MaxN + 1; i++) {
	//	FlagUij[i] = new long*[NodeNum];
	//	for (j = 0; j < NodeNum; j++) {
	//		FlagUij[i][j] = new long[NodeNum];
	//	}
	//}
	//for (i = 0; i < MaxN + 1; i++) {
	//	for (j = 0; j < NodeNum; j++) {
	//		for (p = 0; p < 3; p++) {
	//			FlagUij[i][j][p] = 0;
	//		}
	//	}
	//}



	MT[0].initial(3 * NodeNum, 3 * NodeNum);
	Mtemp.initial(3 * NodeNum, 3 * NodeNum);//ʱʼ
	Mtemp.reinitial();//ֵ0
	MCCSRtemp.initial(3 * NodeNum, 3 * NodeNum);//ʱʼ
	MCCSRtemp.reinitial();//ֵ0
	//MSymCtemp.initial(3 * NodeNum, 3 * NodeNum);//ʱʼ
	//MSymCtemp.reinitial();//ֵ0

	Btemp.initial(3 * NodeNum);//ʱʼ
	BCCSRtemp.initial(3 * NodeNum);
	BSymCtemp.initial(3 * NodeNum);

	//AssistCoef = new double[thread_num][8][9];

	Pthnum = new long[NodeNum];//ϸ
	threads = new pthread_t[thread_num];



	Thread = new PthArg[thread_num];





	for (i = 0; i < thread_num; i++) {
		Thread[i].targs = i;

	}


	// Start Counting...
	m_Counter.StartCount();



	// 2. 

	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthGetMatrixDynaUnknown0ijT0, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}
	//GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);//ݴм

	//MCCSRtemp = 0.0;
	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthGetMatrixDynaUnknown0ij, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}
	MCCSRtemp.sum();
	MTS0.initial(MCCSRtemp.m, MCCSRtemp.num);
	MTS0.row_ptr[0] = 0;
	//ʼ䡪

	PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);



	//
	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthSparseConvertCCSRMTS0, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}


	//m_Counter.EndCount();
	//double DifTime1 = m_Counter.TimeDiff();

	//printf("ʱ:%lf\n", m_Counter.TimeDiff());

	//MCCSRtemp = 0.0;
	//GetMatrixDynaUnknown0ij(m_DSE, MCCSRtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);//м
	//SparseConvert(MCCSRtemp, MTS0);

	//MCCSRtemp = 0.0;

	//GetMatrixDynaKnown0ij(m_DSE, MCCSRtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);

	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthGetMatrixDynaKnown0ij, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}
	MCCSRtemp.sum();
	MTS[0].initial(MCCSRtemp.n, MCCSRtemp.num);
	MTS[0].row_ptr[0] = 0;
	//PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);
	StepNum = 0;
	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthSparseConvertCCSR, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}

	//Ԥ
	GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MT[0]);//Сδ
	delete[] MT;
	MT = 0;


	//ʼ䡪
	//PthDistr(Thread);


	for (StepNum = 1; StepNum <= MaxN; ++StepNum)
	{
		printf("Calculate %ld-th BEM Matrices.\n", StepNum);


		//з
		//for (i = 0; i < thread_num; i++) {
		//	//targs[i] = i;
		//	pthread_create(&(threads[i]), NULL, PthAllTraverse, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		//


		//T1ij
		//GetMatrixDynaT1ij(m_DSE, MCCSRtemp, NodeNum, EleNum, StepNum * 1.0, dt_p, m_InfElePID, m_ElePID);
		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthGetMatrixDynaT2ij, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}

		MCCSRtemp.sum();
		MTS[StepNum].initial(MCCSRtemp.n, MCCSRtemp.num);
		MTS[StepNum].row_ptr[0] = 0;
		//if (StepNum == 1)
		PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);


		//SparseConvert(MCCSRtemp, MTS[StepNum]);//㲢
		for (i = 0; i < thread_num; i++) {
			pthread_create(&(threads[i]), NULL, PthSparseConvertCCSR, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		//



		//GetMatrixDynaT2ij(m_DSE, MCCSRtemp, NodeNum, EleNum, StepNum * 1.0, dt_p, m_InfElePID, m_ElePID);
		//MCCSRtemp = 0.0;
		//T2ij
		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthGetMatrixDynaT1ij, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}

		//SparseAddition(MCCSRtemp, MTS[StepNum], BCCSRtemp);

		//Ӳȫڵ
		
		//PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);
		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthSparseAdditionCCSR, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		MCCSRtemp.sum();
		MTS[StepNum].Matdelete();
		MTS[StepNum].initial(MCCSRtemp.n, MCCSRtemp.num);
		MTS[StepNum].row_ptr[0] = 0;
		PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);



		for (i = 0; i < thread_num; i++) {
			pthread_create(&(threads[i]), NULL, PthSparseConvertCCSR, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}

		//




		//Uij

		//GetMatrixDynaUij(m_DSE, MSymCtemp, NodeNum, EleNum, StepNum * 1.0, dt_p);
		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthGetMatrixDynaSymUij, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		MCCSRtemp.sum();
		MGS[StepNum].initial(MCCSRtemp.n, MCCSRtemp.num);
		MGS[StepNum].row_ptr[0] = 0;
		PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);



		for (i = 0; i < thread_num; i++) {
			pthread_create(&(threads[i]), NULL, PthSparseConvertSymCCSR, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		//








		//for (i = 0; i < thread_num; i++) {
		//	//targs[i] = i;
		//	pthread_create(&(threads[i]), NULL, PthGetMatrixDynaSymUij, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		//MSymCtemp.sum();
		//MGS[StepNum].initial(MSymCtemp.n, MSymCtemp.num);
		//MGS[StepNum].row_ptr[0] = 0;
		//PthDistr(MSymCtemp.row, Thread, MSymCtemp.num, MSymCtemp.n);
		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthSparseConvertSymCCSR, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}

		//SparseConvert(MSymCtemp, MGS[StepNum]);

	}

	// 3. 0ʱ̵֪δ֪˳ΪֲµλƺλΪ0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. n = 1ʱ̵ϵõ n = 1 ʱ̵Ľ


	delete[] threads;


	m_Counter.EndCount();
	double DifTime1 = m_Counter.TimeDiff();
	printf("װʱ:%lf\n", m_Counter.TimeDiff());

	m_Counter.StartCount();

	showMemoryInfo();

	//·á

	char* CurrPath = _getcwd(NULL, 0);//صǰ·
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//صǰļ
	Path = StrCurrPath.substr(0, found + 1);//ϼ·

	string OutFile = "output\\";
	if (GMRES_OUTPUT_SUBDIR[0] != '\0') {
		OutFile += string(GMRES_OUTPUT_SUBDIR) + "\\";
		_mkdir("output");
		string outDir = string("output\\") + GMRES_OUTPUT_SUBDIR;
		_mkdir(outDir.c_str());
	}
	string TempPath = Path + OutFile;


	//Ԫ͡

	/*string Temp1 = "EleType";
	string EtypeNum;
	string ForMatTxt = ".txt";
	FILE* Etype;

	for (i = 0; i <= MaxN; i++) {
		EtypeNum= to_string(i);
		Path = TempPath + Temp1;
		Path = Path + EtypeNum + ForMatTxt;
		const char* OutPath1 = Path.c_str();
		fopen_s(&Etype, OutPath1, "w");
		for (j = 0; j < NodeNum; j++) {
			for (p = 0; p < 3; p++) {
				fprintf_s(Etype, "%8d  ", FlagUij[i][j][p]);
			}
			fprintf_s(Etype, "\n");
		}
		fclose(Etype);
	}*/



	//зԪظ


	//FILE* MTSpare;
	//string SparseTemp = "MTSparse.txt";
	//Path = TempPath + SparseTemp;
	//const char* SparseOutPath = Path.c_str();
	//fopen_s(&MTSpare, SparseOutPath, "w");
	////fprintf_s(MTSpare, "MT\n");
	//for (i = 0; i <= 0; i++) {//  1;MaxN
	//	for (int j = 0; j < NodeNum; j++) {
	//		fprintf_s(MTSpare, "%5d  ", MTS[i].row_ptr[j + 1] - MTS[i].row_ptr[j]);
	//	}
	//	fprintf_s(MTSpare, "\n");
	//}

	//fprintf_s(MTSpare, "MG\n");
	//for (int i = 1; i <= MaxN; i++) {
	//	fprintf_s(MTSpare, "%.8lf", MGS[i].row_ptr[(MGS[i].n)]);
	//	fprintf_s(MTSpare, "\n");
	//}
	//fclose(MTSpare);



	//


	printf("GMRES Solver: Step 1.\n");




	x0 = 0.0;

	// B = MG_1(0) * bd[1].luMG_1(0)ʵ֪bd[1].lu1ʱ̵֪
	//mat_vec_product(MG[0], bd[1].lu, B);
	SparseMul(MTS[0], bd[1].lu, B);

	//  n = 1ʱ̵δ֪bd[1].lt
	gmres(m_DSE, MTS0, B, m_Pre, iterations, error, 1, x0, bd[1].lt, flag, iter);

	//  n = 1ʱ̵δ֪/֪תΪλƺ
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);



	// 4.  n = 2 ~ NStepĽ

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("GMRES Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].luMG_1(0)ʵ֪bd[i].luiʱ̵֪
		//mat_vec_product(MG[0], bd[i].lu, B);
		SparseMul(MTS[0], bd[i].lu, B);
		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;

				SparseMul(MGS[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;

				SparseMul(MTS[j], CB1, CB2);
				B -= CB2;
			}
		}

		//  n = iʱ̵δ֪bd[i].lt
		gmres(m_DSE, MTS0, B, m_Pre, iterations, error, 1, x0, bd[i].lt, flag, iter);

		//  n = iʱ̵δ֪/֪תΪλƺ
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}

	//ϡȡ

	FILE* SparseM;
	double nnum = 9 * NodeNum * NodeNum;
	string Temp = "SparseM.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();
	fopen_s(&SparseM, OutPath, "w");
	fprintf_s(SparseM, "MT\n");
	for (int i = 1; i <= MaxN; i++) {
		fprintf_s(SparseM, "%.8lf", MTS[i].row_ptr[(MTS[i].n)] / nnum);
		fprintf_s(SparseM, "\n");
	}
	fprintf_s(SparseM, "MG\n");
	for (int i = 1; i <= MaxN; i++) {
		fprintf_s(SparseM, "%.8lf", MGS[i].row_ptr[(MGS[i].n)] / nnum);
		fprintf_s(SparseM, "\n");
	}
	fclose(SparseM);



	showMemoryInfo();


	//ͷſռ䡪


	delete[] MGS;
	MGS = 0;
	delete[] MTS;
	MTS = 0;
	Mtemp.finish();

	// End Counting...
	m_Counter.EndCount();





	Temp = "GMRESTime.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "ʱ:\n");
	fprintf(pf, "װʱ:%lf\n", DifTime1);
	fprintf(pf, "ʱ:%lf\n", m_Counter.TimeDiff());
	fclose(pf);





	return 1;
}

int DynaGMRESSolverNewCCSRNewPth(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID)
{
	long i, j;

	//Wmatrix* MT;
	Wmatrix* MG;
	//SymCCSRMat* MGS;ȫ
	//CCSRMat* MTS;ȫ
	WmatrixCSR Mtemp;//ʱ
	//WmatrixCCSR MCCSRtemp;//ʱ
	//WmatrixSymCCSR MSymCtemp;//ʱ
	WvectorCSR Btemp;
	//WvectorCCSR BCCSRtemp; ȫ
	WvectorSymCCSR BSymCtemp;
	long MaxN = NStep;
	

	double dt = DSquareElement::m_DMat.Dt;
	dt_p = DSquareElement::m_DMat.Dt;

	//Wvector B(3 * NodeNum),CB1(3 * NodeNum), CB2(3 * NodeNum);
	//Wvector x0(3 * NodeNum);
	x0.initial(3 * NodeNum);
	CB1.initial(3 * NodeNum);
	CB2.initial(3 * NodeNum);
	Bsum.initial(3 * NodeNum);
	Wvector diagKnown;
	Wvector diagHistoryG;
	Wvector diagHistoryT;
	diagKnown.initial(3 * NodeNum);
	diagHistoryG.initial(3 * NodeNum);
	diagHistoryT.initial(3 * NodeNum);

	//Bvector = new Wvector[thread_num];
	//CBT = new Wvector[MaxN + 1];
	//CBG = new Wvector[MaxN + 1];

	//for (i = 0; i < thread_num; ++i) {
	//	Bvector[i].initial(3 * NodeNum);
	//	Bvector[i] = 0;
	//}
	//for (i = 0; i < MaxN + 1; ++i) {
	//	CBT[i].initial(3 * NodeNum);
	//	CBG[i].initial(3 * NodeNum);
	//}

	// Ԥ
	PreConditioner m_Pre;

	int flag;                                      //ź
	int iter[2];

	Counter m_Counter;
	Counter T_Solution;
	Counter T_GMRES;
	double GMRESTime = 0;

	// ȷķӦ


	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf(": %ld\n", MaxN);

	bool gmresMemoryReuse = GmresCcsrCanReuseMemorySession(NodeNum, EleNum, NStep, MaxN);
	if (gmresMemoryReuse)
		printf("GMRES CCSR memory session hit. Reuse in-memory matrices.\n");


	pthread_t* threads = new pthread_t[thread_num];

	// 1. ʼϵ
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, ΪûMT2(0)MT2(NStep)ҲڳʼλΪ0û
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)MT1[0],Ҳδ֪ MG(0)MG1[0]Ҳ֪
	// MT(1~NStep-1)MT1 + MT2MG(1~NStep-1)MG1 + MG2
	// MG(NStep)MG2[NStep]
	//MT = new Wmatrix[1];
	//MT = new Wmatrix[NStep + 1];
	//MG = new Wmatrix[NStep + 1];
	if (!gmresMemoryReuse) {
	MTS = new CCSRMat[NStep + 1];
	MGS = new SymCCSRMat[NStep + 1];


	//MT[0].initial(3 * NodeNum, 3 * NodeNum);
	//Mtemp.initial(3 * NodeNum, 3 * NodeNum);//ʱʼ
	//Mtemp.reinitial();//ֵ0
	// 
	// 
	// 
	//MCCSRtemp.initial(3 * NodeNum, 3 * NodeNum);//ʱʼ

	long ccsrTempColumnDofs = 3 * NodeNum;
	if (ccsrTempColumnDofs < 900)
		ccsrTempColumnDofs = 900;
	MCCSRtemp.initial(3 * NodeNum, ccsrTempColumnDofs);// Largescale capacity, block columns = max(NodeNum, 300)

	MCCSRtemp.reinitial();//ֵ0
	//MSymCtemp.initial(3 * NodeNum, 3 * NodeNum);//ʱʼ
	//MSymCtemp.reinitial();//ֵ0

	//Btemp.initial(3 * NodeNum);//ʱʼ
	BCCSRtemp.initial(3 * NodeNum);
	BSymCtemp.initial(3 * NodeNum);

	//AssistCoef = new double[thread_num][8][9];

	Pthnum= new long[NodeNum];//ϸ

	Thread = new PthArg[thread_num];
	ConvolThreadT = new PthConvol[thread_num];
	ConvolThreadG = new PthConvol[thread_num];

	SpMVThreadT = new PthArg * [MaxN + 1]; //SpMV ʼ  T
	SpMVThreadG = new PthArg * [MaxN + 1]; ////SpMV ʼ  U


	for (i = 0; i < MaxN + 1; i++) {
		SpMVThreadT[i] = new PthArg[thread_num];
		SpMVThreadG[i] = new PthArg[thread_num];
	}

	for (i = 0; i < thread_num; i++) {
		Thread[i].targs = i;

	}
	}

	double starttime, endtime;
	bool gmresCacheHit = false;
	string gmresCacheReason;
	// Start Counting...
	m_Counter.StartCount();

	if (gmresMemoryReuse)
	{
		starttime = omp_get_wtime();
		GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MTS0);
		endtime = omp_get_wtime();
		printf("Ԥʱ %lf\n", endtime - starttime);
		goto GMRES_CCSR_SOLVE_READY;
	}

	printf("GMRES CCSR cache path: %s\n", GMRES_CCSR_CACHE_PATH);
	starttime = omp_get_wtime();
	gmresCacheHit = GmresReadCcsrCache(GMRES_CCSR_CACHE_PATH, NodeNum, EleNum, NStep, MaxN, thread_num, gmresCacheReason);
	endtime = omp_get_wtime();
	if (gmresCacheHit)
	{
		printf("GMRES CCSR cache hit, read time %lf\n", endtime - starttime);
		starttime = omp_get_wtime();
		GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MTS0);
		endtime = omp_get_wtime();
		printf("Ԥʱ %lf\n", endtime - starttime);
		goto GMRES_CCSR_MATRICES_READY;
	}
	printf("GMRES CCSR cache miss: %s\n", gmresCacheReason.c_str());


	// 2. 

	//for (i = 0; i < thread_num; i++) {
	//	pthread_create(&(threads[i]), NULL, PthGetMatrixDynaUnknown0ijT0, &Thread[i]);
	//}
	//for (i = 0; i < thread_num; i++) {
	//	pthread_join(threads[i], NULL);
	//}
	//GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);//ݴм

	//MCCSRtemp = 0.0;
	printf("Aʼ\n");
	GmresThreadDiagPrintSystemInfo(thread_num);
	GmresThreadDiagBegin(thread_num);
	int* gmresThreadCreateCodes = new int[thread_num];
	long gmresThreadCreateSuccess = 0;
	starttime = omp_get_wtime();
	for (i = 0; i < thread_num; i++) {
		gmresThreadCreateCodes[i] = pthread_create(&(threads[i]), NULL, PthGetMatrixDynaUnknown0ij, &Thread[i]);
		if (gmresThreadCreateCodes[i] == 0)
			++gmresThreadCreateSuccess;
	}
	GmresThreadDiagPrintCreateResult(thread_num, gmresThreadCreateSuccess, gmresThreadCreateCodes);
	for (i = 0; i < thread_num; i++) {
		if (gmresThreadCreateCodes[i] == 0)
			pthread_join(threads[i], NULL);
	}
	GmresThreadDiagEnd();
	delete[] gmresThreadCreateCodes;
	endtime = omp_get_wtime(); // ȡʱ

	printf("Aʱ %lf\n", endtime - starttime);

	MCCSRtemp.sum();
	MTS0.initial(MCCSRtemp.m, MCCSRtemp.num);
	MTS0.row_ptr[0] = 0;
	//ʼ䡪



	starttime = omp_get_wtime();
	PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.m);
	endtime = omp_get_wtime(); // ȡʱ

	printf("Aʱ %lf\n", endtime - starttime);

	for (i = 0; i < thread_num; i++) {
		SpMVThreadT[0][i] = Thread[i];  //δ֪
		SpMVThreadG[0][i] = Thread[i];  //֪
	}

	//
	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthSparseConvertCCSRMTS0, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}


	//m_Counter.EndCount();
	//double DifTime1 = m_Counter.TimeDiff();

	//printf("ʱ:%lf\n", m_Counter.TimeDiff());

	//MCCSRtemp = 0.0;
	//GetMatrixDynaUnknown0ij(m_DSE, MCCSRtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);//м
	//SparseConvert(MCCSRtemp, MTS0);

	//MCCSRtemp = 0.0;

	//GetMatrixDynaKnown0ij(m_DSE, MCCSRtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);

	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthGetMatrixDynaKnown0ij, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}
	MCCSRtemp.sum();
	MTS[0].initial(MCCSRtemp.m, MCCSRtemp.num);//֪
	MTS[0].row_ptr[0] = 0;
	PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.m);
	StepNum = 0;
	for (i = 0; i < thread_num; i++) {
		pthread_create(&(threads[i]), NULL, PthSparseConvertCCSR, &Thread[i]);
	}
	for (i = 0; i < thread_num; i++) {
		pthread_join(threads[i], NULL);
	}

	//Ԥ
	//GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MT[0]);//Сδ
	//delete[] MT;
	//MT = 0;


	//Ԥ󡪡ϡʽ

	starttime = omp_get_wtime();
	GMRESPreConditioner(m_DSE, m_Pre, NodeNum, maxleafpointnum, MTS0);
	endtime = omp_get_wtime(); // ȡʱ

	printf("Ԥʱ %lf\n", endtime - starttime);
	

	//ʼ䡪
	//PthDistr(Thread);


	for (StepNum = 1; StepNum <= MaxN; ++StepNum)
	{
		printf("Calculate %ld-th BEM Matrices.\n", StepNum);

		
		//з
		//for (i = 0; i < thread_num; i++) {
		//	//targs[i] = i;
		//	pthread_create(&(threads[i]), NULL, PthAllTraverse, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		//


		//T1ij
		//GetMatrixDynaT1ij(m_DSE, MCCSRtemp, NodeNum, EleNum, StepNum * 1.0, dt_p, m_InfElePID, m_ElePID);



		//m_Counter.StartCount();
		//for (i = 0; i < thread_num; i++) {
		//	//targs[i] = i;
		//	pthread_create(&(threads[i]), NULL, PthGetMatrixDynaT1ij, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}

		//m_Counter.EndCount();
		//printf("ʱ:%lf\n", m_Counter.TimeDiff());
		//m_Counter.StartCount();

		//for (i = 0; i < thread_num; i++) {
		//	//targs[i] = i;
		//	pthread_create(&(threads[i]), NULL, PthGetMatrixDynaT1ijNew, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		//m_Counter.EndCount();
		//printf("ʱ:%lf\n", m_Counter.TimeDiff());





		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthGetMatrixDynaT2ijNew, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		 

		MCCSRtemp.sum();
		MTS[StepNum].initial(MCCSRtemp.m, MCCSRtemp.num);
		MTS[StepNum].row_ptr[0] = 0;
		//if (StepNum == 1) {// жϿе⣬ϸ
			//PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);
		//}
		PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.m); //ע

		//SparseConvert(MCCSRtemp, MTS[StepNum]);//㲢
		for (i = 0; i < thread_num; i++) {
			pthread_create(&(threads[i]), NULL, PthSparseConvertCCSR, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		//







		//GetMatrixDynaT2ij(m_DSE, MCCSRtemp, NodeNum, EleNum, StepNum * 1.0, dt_p, m_InfElePID, m_ElePID);
		//MCCSRtemp = 0.0;
		//T1ij
		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthGetMatrixDynaT1ijNew, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}

		//SparseAddition(MCCSRtemp, MTS[StepNum], BCCSRtemp);

		//Ӳȫڵ
		//MCCSRtemp.sum();
		//MTS[StepNum].initial(MCCSRtemp.n, MCCSRtemp.num);
		//MTS[StepNum].row_ptr[0] = 0;
		//PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);//ע




		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthSparseAdditionCCSR, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}

		MCCSRtemp.sum();
		MTS[StepNum].Matdelete();
		MTS[StepNum].initial(MCCSRtemp.m, MCCSRtemp.num);
		MTS[StepNum].row_ptr[0] = 0;
		PthDistrAdd(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.m);
		//PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.n);//ע

		for (i = 0; i < thread_num; i++) {
			SpMVThreadT[StepNum][i] = Thread[i];
		}

		for (i = 0; i < thread_num; i++) {
			//pthread_create(&(threads[i]), NULL, PthSparseConvertCCSR, &Thread[i]);
			pthread_create(&(threads[i]), NULL, PthSparseAddConvertCCSR, &Thread[i]);//ע
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}

		//

		

		//UijŻ

		//GetMatrixDynaUij(m_DSE, MSymCtemp, NodeNum, EleNum, StepNum * 1.0, dt_p);
		for (i = 0; i < thread_num; i++) {
			//targs[i] = i;
			pthread_create(&(threads[i]), NULL, PthGetMatrixDynaSymUijNew, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		MCCSRtemp.sum();
		MGS[StepNum].initial(MCCSRtemp.m, MCCSRtemp.num);
		MGS[StepNum].row_ptr[0] = 0;
		PthDistr(MCCSRtemp.row, Thread, MCCSRtemp.num, MCCSRtemp.m);//ע
		
		for (i = 0; i < thread_num; i++) {
			SpMVThreadG[StepNum][i] = Thread[i];
		}

		for (i = 0; i < thread_num; i++) {
			pthread_create(&(threads[i]), NULL, PthSparseConvertSymCCSR, &Thread[i]);
		}
		for (i = 0; i < thread_num; i++) {
			pthread_join(threads[i], NULL);
		}
		//



		//for (i = 0; i < thread_num; i++) {
		//	//targs[i] = i;
		//	pthread_create(&(threads[i]), NULL, PthGetMatrixDynaSymUij, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		//MSymCtemp.sum();
		//MGS[StepNum].initial(MSymCtemp.n, MSymCtemp.num);
		//MGS[StepNum].row_ptr[0] = 0;
		//PthDistr(MSymCtemp.row, Thread, MSymCtemp.num, MSymCtemp.n);
		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthSparseConvertSymCCSR, &Thread[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}

		//SparseConvert(MSymCtemp, MGS[StepNum]);

	}

	{
		__int64 gmresCacheFileSize = 0;
		string gmresCacheWriteReason;
		starttime = omp_get_wtime();
		bool gmresCacheWriteOk = GmresWriteCcsrCache(GMRES_CCSR_CACHE_PATH, NodeNum, EleNum, NStep, MaxN, gmresCacheFileSize, gmresCacheWriteReason);
		endtime = omp_get_wtime();
		if (gmresCacheWriteOk)
			printf("GMRES CCSR cache written, write time %lf, file size %lld bytes\n", endtime - starttime, gmresCacheFileSize);
		else
			printf("GMRES CCSR cache write failed: %s\n", gmresCacheWriteReason.c_str());
	}

GMRES_CCSR_MATRICES_READY:
	//תв

	for (i = 0; i < thread_num; i++) {
		ConvolThreadT[i].initial(MaxN + 1);
		ConvolThreadG[i].initial(MaxN + 1);
	}
	SpmvConvertConvol(SpMVThreadT, ConvolThreadT, MaxN);
	SpmvConvertConvol(SpMVThreadG, ConvolThreadG, MaxN);
	GmresCcsrMarkMemorySessionReady(NodeNum, EleNum, NStep, MaxN);

GMRES_CCSR_SOLVE_READY:
	



	// 3. 0ʱ̵֪δ֪˳ΪֲµλƺλΪ0


	StepNum = 0;   //ѭ 

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	//άȲO(n)ôС

	//for (i = 0; i < thread_num; i++) {
	//	pthread_create(&(threads[i]), NULL, PthbdConvert, &SpMVThreadG[0][i]);
	//}
	//for (i = 0; i < thread_num; i++) {
	//	pthread_join(threads[i], NULL);
	//}
	//for (i = 0; i < thread_num; i++) {
	//	pthread_create(&(threads[i]), NULL, PthbdAllLocToAbs, &SpMVThreadG[0][i]);
	//}
	//for (i = 0; i < thread_num; i++) {
	//	pthread_join(threads[i], NULL);
	//}
	//

	// 3. n = 1ʱ̵ϵõ n = 1 ʱ̵Ľ




	m_Counter.EndCount();
	double DifTime1 = m_Counter.TimeDiff();
	printf("ʱ:%lf\n", m_Counter.TimeDiff());
	
	T_Solution.StartCount();

	

	//·á

	string Path;
	_mkdir("output");
	string TempPath = "output\\";
	if (GMRES_OUTPUT_SUBDIR[0] != '\0') {
		string outDir = string("output\\") + GMRES_OUTPUT_SUBDIR;
		_mkdir(outDir.c_str());
		TempPath = outDir + "\\";
	}
	//зԪظ


	//FILE* MTSpare;
	//string SparseTemp = "MTSparse.txt";
	//Path = TempPath + SparseTemp;
	//const char* SparseOutPath = Path.c_str();
	//fopen_s(&MTSpare, SparseOutPath, "w");
	////fprintf_s(MTSpare, "MT\n");
	//for (i = 0; i <= MaxN; i++) {//  1;MaxN
	//	for (int j = 0; j < NodeNum; j++) {
	//		fprintf_s(MTSpare, "%5d  ", MTS[i].row_ptr[j + 1] - MTS[i].row_ptr[j]);
	//	}
	//	fprintf_s(MTSpare, "\n");
	//}

	////fprintf_s(MTSpare, "MG\n");
	////for (int i = 1; i <= MaxN; i++) {
	////	fprintf_s(MTSpare, "%.8lf", MGS[i].row_ptr[(MGS[i].n)]);
	////	fprintf_s(MTSpare, "\n");
	////}
	//fclose(MTSpare);



	//


	printf("GMRES Solver: Step 1.\n");
	showMemoryInfo();



	StepNum++;   

	x0 = 0.0;

	// B = MG_1(0) * bd[1].luMG_1(0)ʵ֪bd[1].lu1ʱ̵֪
	//mat_vec_product(MG[0], bd[1].lu, B);
	SparseMul(MTS[0], bd[1].lu, Bsum);
	diagKnown = Bsum;
	diagHistoryG = 0.0;
	diagHistoryT = 0.0;


	//for (i = 0; i < thread_num; i++) {
	//	pthread_create(&(threads[i]), NULL, PthSpMVKnown, &SpMVThreadG[0][i]);
	//}
	//for (i = 0; i < thread_num; i++) {
	//	pthread_join(threads[i], NULL);
	//}

	//SpMVThreadT[0][i] //δ֪
	//SpMVThreadG[0][i] //֪

	T_GMRES.StartCount();
	//  n = 1ʱ̵δ֪bd[1].lt
	gmres(m_DSE, MTS0, Bsum, m_Pre, iterations, error, 1, x0, bd[1].lt, flag, iter);
	WriteOldRhsBreakdownCsv("single", StepNum, m_DSE, EleNum,
		diagKnown, diagHistoryG, diagHistoryT, Bsum, bd[StepNum].lt, MTS0);
	T_GMRES.EndCount();
	GMRESTime += T_GMRES.TimeDiff();

	//  n = 1ʱ̵δ֪/֪תΪλƺ
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	//for (i = 0; i < thread_num; i++) {
	//	pthread_create(&(threads[i]), NULL, PthbdConvert, &SpMVThreadG[0][i]);
	//}
	//for (i = 0; i < thread_num; i++) {
	//	pthread_join(threads[i], NULL);
	//}

	//for (i = 0; i < thread_num; i++) {
	//	pthread_create(&(threads[i]), NULL, PthbdAllLocToAbs, &SpMVThreadG[0][i]);
	//}
	//for (i = 0; i < thread_num; i++) {
	//	pthread_join(threads[i], NULL);
	//}


	// 4.  n = 2 ~ NStepĽ

	for (StepNum = 2; StepNum <= NStep; ++StepNum)
	{
		// TEST
		printf("GMRES Solver: Step %ld.\n", StepNum);

		//㡪 
		//for (i = 0; i < thread_num; i++) {
		//	Bvector[i] = 0;
		//}

		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthBvectorInitial, &ConvolThreadT[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}


		Bsum = 0;

		// B = MG_1(0) * bd[i].luMG_1(0)ʵ֪bd[i].luiʱ̵֪
		//mat_vec_product(MG[0], bd[i].lu, B);
		
		

		SparseMul(MTS[0], bd[StepNum].lu, Bsum);
		diagKnown = Bsum;
		diagHistoryG = 0.0;
		diagHistoryT = 0.0;

		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthSpMVKnown, &SpMVThreadG[0][i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL); 
		//}
		
		
		// j = 1..i
		for (ConvoluNum = 1; ConvoluNum <= StepNum; ++ConvoluNum)
		{
			if (ConvoluNum <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((StepNum - ConvoluNum + 1) > 0 && (StepNum - ConvoluNum + 1) < StepNum)
					CB1 += bd[StepNum - ConvoluNum + 1].lt;
				if ((StepNum - ConvoluNum) > 0 && (StepNum - ConvoluNum) < StepNum)
					CB1.plus(bd[StepNum - ConvoluNum].lt, 2.0);
				if ((StepNum - ConvoluNum - 1) > 0 && (StepNum - ConvoluNum - 1) < StepNum)
					CB1 += bd[StepNum - ConvoluNum - 1].lt;


				//SparseMul(MGS[ConvoluNum], CB1, CB2);

				for (i = 0; i < thread_num; i++) {
					pthread_create(&(threads[i]), NULL, PthSpMVG, &SpMVThreadG[ConvoluNum][i]);
				}
				for (i = 0; i < thread_num; i++) {
					pthread_join(threads[i], NULL);
				}

				Bsum += CB2;
				diagHistoryG += CB2;

				// [T]
				CB1 = 0.0;
				if ((StepNum - ConvoluNum + 1) > 0 && (StepNum - ConvoluNum + 1) < StepNum)
					CB1 += bd[StepNum - ConvoluNum + 1].lu;
				if ((StepNum - ConvoluNum) > 0 && (StepNum - ConvoluNum) < StepNum)
					CB1.plus(bd[StepNum - ConvoluNum].lu, 2.0);
				if ((StepNum - ConvoluNum - 1) > 0 && (StepNum - ConvoluNum - 1) < StepNum)
					CB1 += bd[StepNum - ConvoluNum - 1].lu;

				//SparseMul(MTS[ConvoluNum], CB1, CB2);
				for (i = 0; i < thread_num; i++) {
					pthread_create(&(threads[i]), NULL, PthSpMVT, &SpMVThreadT[ConvoluNum][i]);
				}
				for (i = 0; i < thread_num; i++) {
					pthread_join(threads[i], NULL);
				}
				Bsum -= CB2;
				diagHistoryT -= CB2;
			}
		}


		//ConvolBvector(CBT, CBG, bd, StepNum, MaxN);//  Ҫ



		//   ȶԴ ת
		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthConvolBvector, &ConvolThreadG[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		////
		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthConvolutinG, &ConvolThreadG[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		//for (i = 0; i < thread_num; i++) {   
		//	Bsum += Bvector[i];
		//	Bvector[i] = 0;
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthConvolutinT, &ConvolThreadT[i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	Bsum -= Bvector[i];
		//}


		//  n = iʱ̵δ֪bd[i].lt
		T_GMRES.StartCount();
		
		gmres(m_DSE, MTS0, Bsum, m_Pre, iterations, error, 1, x0, bd[StepNum].lt, flag, iter);
		WriteOldRhsBreakdownCsv("single", StepNum, m_DSE, EleNum,
			diagKnown, diagHistoryG, diagHistoryT, Bsum, bd[StepNum].lt, MTS0);

		T_GMRES.EndCount();
		GMRESTime += T_GMRES.TimeDiff();


		//  n = iʱ̵δ֪/֪תΪλƺ
		bd[StepNum].Convert(m_DSE, EleNum);
		bd[StepNum].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthbdConvert, &SpMVThreadG[0][i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}

		//for (i = 0; i < thread_num; i++) {
		//	pthread_create(&(threads[i]), NULL, PthbdAllLocToAbs, &SpMVThreadG[0][i]);
		//}
		//for (i = 0; i < thread_num; i++) {
		//	pthread_join(threads[i], NULL);
		//}
	}

	//ϡȡ

	//FILE* SparseM;
	//double nnum = 9 * NodeNum * NodeNum;
	string Temp = "SparseM.txt";
	//Path = TempPath + Temp;
	const char* OutPath = Path.c_str();
	//fopen_s(&SparseM, OutPath, "w");
	//fprintf_s(SparseM, "MT\n");
	//for (int i = 1; i <= MaxN; i++) {
	//	fprintf_s(SparseM, "%.8lf", MTS[i].row_ptr[(MTS[i].n)] / nnum);
	//	fprintf_s(SparseM, "\n");
	//}
	//fprintf_s(SparseM, "MG\n");
	//for (int i = 1; i <= MaxN; i++) {
	//	fprintf_s(SparseM, "%.8lf", MGS[i].row_ptr[(MGS[i].n)] / nnum);
	//	fprintf_s(SparseM, "\n");
	//}
	//fclose(SparseM);



	showMemoryInfo();



	//ͷſռ䡪
	delete[] threads;

	if (!GMRES_BATCH_MEMORY_SESSION) {
		delete[] MGS;
		MGS = 0;
		delete[] MTS;
		MTS = 0;
	}
	//Mtemp.finish();

	// End Counting...
	T_Solution.EndCount();

	



	Temp = "GMRESTime.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "ʱ:\n");
	fprintf(pf, "װʱ:%lf\n", DifTime1);
	fprintf(pf, "ʱ:%lf\n", T_Solution.TimeDiff());
	fprintf(pf, "GMRES:%lf\n", GMRESTime);
	fclose(pf);





	return 1;
}














