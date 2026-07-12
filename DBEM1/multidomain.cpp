#include "multidomain.h"
#include "Gauss_solver.h"
#include "pre_gmres.h"
#include "output_path.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <pthread.h>

extern double*** AssistCoef;
extern double*** AssistUij;
extern double*** AssistTij;

static int IsFiniteScalar(double value)
{
	return value == value && fabs(value) <= (std::numeric_limits<double>::max)();
}

static int IsSupportedBlockBCID(int bcid)
{
	return bcid == 123 || bcid == 456;
}

static int g_validationMetricsValid = 0;
static int g_validationFinalFlag = -999;
static long g_validationSolvedSteps = 0;
static long g_validationTotalIter0 = 0;
static long g_validationTotalIter1 = 0;
static long g_validationGlobalUnknownBlocks = 0;
static long g_validationMTS0Blocks = 0;
static double g_validationAvgMTSBlocks = 0.0;
static double g_validationAvgMGSBlocks = 0.0;
static double g_validationMatrixAssemblyTime = 0.0;
static double g_validationGMRESTime = 0.0;
static int g_validationDomainCount = 0;
static int g_validationMatrixStructureValid = 0;
static long g_validationMatrixZeroRows = 0;
static long g_validationMatrixZeroCols = 0;
static long g_validationMatrixNonFiniteBlocks = 0;
static long g_validationMatrixDiagonalBlocks = 0;
static long g_validationMatrixInvertibleDiagonalBlocks = 0;
static int g_validationInterfaceDofInvariantValid = 0;
static long g_validationInterfaceDofInvariantFailures = 0;
static int g_validationMultiRateTimeMode = 0;
static double g_validationMacroDt = 0.0;
static double g_validationBetaLimit = 0.0;
static double g_validationAveElementSize = 0.0;
static long g_validationMaxSubstepsPerMacro = 0;
static std::vector<long> g_validationDomainSubsteps;
static std::vector<double> g_validationDomainStableDt;
static std::vector<double> g_validationDomainLocalDt;
static std::vector<double> g_validationDomainBetaActual;
static std::vector<double> g_validationDomainC1;
static std::vector<double> g_validationDomainC2;

static long ElementLocalNodeId(long ele, int localNode);
static int InterfaceLocalAForB(const InterfacePair& itf, int localB);
static long InterfaceNodeA(const InterfacePair& itf, int localA);
static long InterfaceNodeBForA(const InterfacePair& itf, int localA);

static int MultiDomainValidationOutputEnabled()
{
	char* value = 0;
	size_t valueLen = 0;
	if (_dupenv_s(&value, &valueLen, "DBEM_VALIDATION_OUTPUT") != 0)
		return 1;
	if (!value)
		return 1;
	int enabled = !(strcmp(value, "0") == 0 ||
		strcmp(value, "false") == 0 ||
		strcmp(value, "FALSE") == 0 ||
		strcmp(value, "off") == 0 ||
		strcmp(value, "OFF") == 0 ||
		strcmp(value, "no") == 0 ||
		strcmp(value, "NO") == 0);
	free(value);
	return enabled;
}

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

static int EnvFlagEnabled(const char* name)
{
	char* value = 0;
	size_t valueLen = 0;
	_dupenv_s(&value, &valueLen, name);
	int enabled = value && (strcmp(value, "1") == 0 ||
		strcmp(value, "true") == 0 ||
		strcmp(value, "TRUE") == 0 ||
		strcmp(value, "on") == 0 ||
		strcmp(value, "ON") == 0 ||
		strcmp(value, "yes") == 0 ||
		strcmp(value, "YES") == 0);
	if (value)
		free(value);
	return enabled;
}

static int EnvFlagDisabled(const char* name)
{
	char* value = 0;
	size_t valueLen = 0;
	_dupenv_s(&value, &valueLen, name);
	int disabled = value && (strcmp(value, "0") == 0 ||
		strcmp(value, "false") == 0 ||
		strcmp(value, "FALSE") == 0 ||
		strcmp(value, "off") == 0 ||
		strcmp(value, "OFF") == 0 ||
		strcmp(value, "no") == 0 ||
		strcmp(value, "NO") == 0 ||
		strcmp(value, "serial") == 0 ||
		strcmp(value, "SERIAL") == 0);
	if (value)
		free(value);
	return disabled;
}

static int MultiDomainPthreadAssemblyDisabled()
{
	return EnvFlagDisabled("DBEM_MULTIDOMAIN_GMRES_PTHREAD");
}

static int MultiDomainPthreadCheckEnabled()
{
	return EnvFlagEnabled("DBEM_MULTIDOMAIN_GMRES_PTHREAD_CHECK");
}

static FILE* OpenDiagnosticCsv(const char* path, const char* header)
{
	std::string actualPath = DBEMOutputPath(path);
	FILE* test = 0;
	fopen_s(&test, actualPath.c_str(), "r");
	int needsHeader = test ? 0 : 1;
	if (test)
		fclose(test);

	FILE* fp = 0;
	fopen_s(&fp, actualPath.c_str(), "a");
	if (fp && needsHeader && header)
		fprintf(fp, "%s\n", header);
	return fp;
}

struct VectorStats {
	double l2;
	double maxAbs;
	long nonFinite;
};

static VectorStats ComputeVectorStats(const Wvector& value)
{
	VectorStats stats;
	stats.l2 = 0.0;
	stats.maxAbs = 0.0;
	stats.nonFinite = 0;
	if (!value.b)
		return stats;
	for (long i = 0; i < value.n; ++i)
	{
		double v = value.b[i];
		if (!IsFiniteScalar(v))
		{
			++stats.nonFinite;
			continue;
		}
		stats.l2 += v * v;
		double av = fabs(v);
		if (av > stats.maxAbs)
			stats.maxAbs = av;
	}
	stats.l2 = sqrt(stats.l2);
	return stats;
}

static void PrintRhsSummary(long step,
	const Wvector& known,
	const Wvector& historyG,
	const Wvector& historyT,
	const Wvector& total,
	double tolerance)
{
	VectorStats knownStats = ComputeVectorStats(known);
	VectorStats historyGStats = ComputeVectorStats(historyG);
	VectorStats historyTStats = ComputeVectorStats(historyT);
	VectorStats totalStats = ComputeVectorStats(total);
	double threshold = totalStats.l2 > 0.0 ? tolerance * totalStats.l2 : tolerance;
	printf("MultiDomain RHS step %ld: knownL2=%e historyGL2=%e historyTL2=%e totalL2=%e threshold=%e maxKnown=%e maxHistoryG=%e maxHistoryT=%e maxTotal=%e nonFinite=%ld/%ld/%ld/%ld\n",
		step,
		knownStats.l2,
		historyGStats.l2,
		historyTStats.l2,
		totalStats.l2,
		threshold,
		knownStats.maxAbs,
		historyGStats.maxAbs,
		historyTStats.maxAbs,
		totalStats.maxAbs,
		knownStats.nonFinite,
		historyGStats.nonFinite,
		historyTStats.nonFinite,
		totalStats.nonFinite);
}

static void PrintLinearResidualSummary(const char* label,
	long step,
	CCSRMat& matrix,
	const Wvector& rhs,
	const Wvector& x,
	double tolerance)
{
	Wvector ax(rhs.n, 0);
	SparseMul(matrix, (Wvector&)x, ax);
	double residual2 = 0.0;
	double rhs2 = 0.0;
	double maxAbs = 0.0;
	long rhsNonFinite = 0;
	long xNonFinite = 0;
	long axNonFinite = 0;
	for (long i = 0; i < rhs.n; ++i)
	{
		double rv = rhs.b[i];
		double xv = i < x.n ? x.b[i] : 0.0;
		double av = ax.b[i];
		if (!IsFiniteScalar(rv))
			++rhsNonFinite;
		if (!IsFiniteScalar(xv))
			++xNonFinite;
		if (!IsFiniteScalar(av))
			++axNonFinite;
		if (!IsFiniteScalar(rv) || !IsFiniteScalar(av))
			continue;
		double diff = av - rv;
		residual2 += diff * diff;
		rhs2 += rv * rv;
		if (fabs(diff) > maxAbs)
			maxAbs = fabs(diff);
	}
	double residual = sqrt(residual2);
	double rhsNorm = sqrt(rhs2);
	double threshold = rhsNorm > 0.0 ? tolerance * rhsNorm : tolerance;
	printf("MultiDomain %s step %ld residual: residual=%e threshold=%e maxAbs=%e rhsNonFinite=%ld candidateNonFinite=%ld axNonFinite=%ld\n",
		label, step, residual, threshold, maxAbs, rhsNonFinite, xNonFinite, axNonFinite);
}

static void WriteRhsBreakdownCsv(const char* solver,
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
		"solver,step,rowNode,rowElement,localNode,domain,surfaceType,bcid,cx,cy,cz,component,knownRhs,historyG,historyT,totalRhs,Ax,residual,solution");
	if (!fp)
		return;

	const char comp[3] = { 'x', 'y', 'z' };
	long rowBlocks = total.n / 3;
	for (long node = 0; node < rowBlocks; ++node)
	{
		long ele = node;
		int localNode = 0;
		if (ele < 0 || ele >= eleCount)
			continue;
		for (int k = 0; k < 3; ++k)
		{
			long p = node * 3 + k;
			double residual = ax.b[p] - total.b[p];
			fprintf(fp, "%s,%ld,%ld,%ld,%d,%d,%d,%d,%.17g,%.17g,%.17g,%c,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
				solver,
				step,
				node,
				ele,
				localNode,
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

static void WriteInterfaceTransferAuditCsv(const MultiDomainModel& model,
	DSquareElement* elements,
	BoundaryValue& stepBd,
	long step)
{
	if ((!MultiDomainValidationOutputEnabled() && !RhsDiagnosticEnabled()) || !elements)
		return;
	FILE* fp = OpenDiagnosticCsv("output\\interface_transfer_audit.csv",
		"solver,step,interfaceId,domainA,domainB,normalSign,localNode,nodeA,nodeB,eleA,eleB,component,uA,uB,uJump,tA,tB,tResidual");
	if (!fp)
		return;
	const char comp[3] = { 'x', 'y', 'z' };
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		const InterfacePair& itf = model.interfaces[i];
		int tSignB = itf.normalSign == 0 ? -1 : itf.normalSign;
		for (int localNode = 0; localNode < 1; ++localNode)
		{
			long nodeA = InterfaceNodeA(itf, localNode);
			long nodeB = InterfaceNodeBForA(itf, localNode);
			long a = nodeA * 3;
			long b = nodeB * 3;
			for (int k = 0; k < 3; ++k)
			{
				double uA = stepBd.lu.b[a + k];
				double uB = stepBd.lu.b[b + k];
				double tA = stepBd.lt.b[a + k];
				double tB = stepBd.lt.b[b + k];
				fprintf(fp, "multidomain,%ld,%ld,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%c,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
					step,
					(long)i,
					itf.domainA,
					itf.domainB,
					tSignB,
					localNode,
					nodeA,
					nodeB,
					itf.eleA,
					itf.eleB,
					comp[k],
					uA,
					uB,
					uA - uB,
					tA,
					tB,
					tB - (double)tSignB * tA);
			}
		}
	}
	fclose(fp);
}

static void WriteMultiRateStateCsv(const MultiDomainModel& model,
	DSquareElement* elements,
	BoundaryValue& stepBd,
	long macroStep)
{
	if (!MultiDomainValidationOutputEnabled() || !elements || !model.IsActive())
		return;
	FILE* fp = OpenDiagnosticCsv("output\\validation_multirate_state.csv",
		"macroStep,domain,localStep,time,nodeId,element,surfaceType,bcid,x,y,z,ux,uy,uz,tx,ty,tz");
	if (!fp)
		return;
	double time = model.macroDt * (double)macroStep;
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		long localStep = macroStep * (domain.substeps > 0 ? domain.substeps : 1);
		for (size_t i = 0; i < domain.elementIds.size(); ++i)
		{
			long ele = domain.elementIds[i];
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long nodeId = elements[ele].m_nodeID[localNode];
				long base = nodeId * 3;
				Point& pt = elements[ele].m_nodelist[nodeId];
				fprintf(fp,
					"%ld,%d,%ld,%.17g,%ld,%ld,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
					macroStep,
					domain.id,
					localStep,
					time,
					nodeId,
					ele,
					elements[ele].SurfaceType,
					elements[ele].BCID,
					pt.pt[0],
					pt.pt[1],
					pt.pt[2],
					stepBd.lu.b[base + 0],
					stepBd.lu.b[base + 1],
					stepBd.lu.b[base + 2],
					stepBd.lt.b[base + 0],
					stepBd.lt.b[base + 1],
					stepBd.lt.b[base + 2]);
			}
		}
	}
	fclose(fp);
}

static int SameBlockAndSign(const VariableRef& a, const VariableRef& b, double signScale)
{
	return a.globalBlock == b.globalBlock &&
		fabs(b.sign - signScale * a.sign) <= 1.0e-12 &&
		a.known == b.known;
}

void WriteMultiDomainValidationMetrics(FILE* fp)
{
	if (!fp)
		return;
	fprintf(fp, "MultiDomainMetricsValid=%d\n", g_validationMetricsValid);
	fprintf(fp, "MultiDomainFinalFlag=%d\n", g_validationFinalFlag);
	fprintf(fp, "MultiDomainSolvedSteps=%ld\n", g_validationSolvedSteps);
	fprintf(fp, "MultiDomainTotalOuterIterations=%ld\n", g_validationTotalIter0);
	fprintf(fp, "MultiDomainTotalInnerIterations=%ld\n", g_validationTotalIter1);
	fprintf(fp, "MultiDomainAverageOuterIterations=%lf\n",
		g_validationSolvedSteps > 0 ? (double)g_validationTotalIter0 / (double)g_validationSolvedSteps : 0.0);
	fprintf(fp, "MultiDomainAverageInnerIterations=%lf\n",
		g_validationSolvedSteps > 0 ? (double)g_validationTotalIter1 / (double)g_validationSolvedSteps : 0.0);
	fprintf(fp, "MultiDomainDomainCount=%d\n", g_validationDomainCount);
	fprintf(fp, "MultiDomainGlobalUnknownBlocks=%ld\n", g_validationGlobalUnknownBlocks);
	fprintf(fp, "MultiDomainMTS0NonZeroBlocks=%ld\n", g_validationMTS0Blocks);
	fprintf(fp, "MultiDomainMTSAverageBlocks=%lf\n", g_validationAvgMTSBlocks);
	fprintf(fp, "MultiDomainMGSAverageBlocks=%lf\n", g_validationAvgMGSBlocks);
	fprintf(fp, "MultiDomainMatrixAssemblyTime=%lf\n", g_validationMatrixAssemblyTime);
	fprintf(fp, "MultiDomainGMRESTime=%lf\n", g_validationGMRESTime);
	fprintf(fp, "MultiDomainMatrixStructureValid=%d\n", g_validationMatrixStructureValid);
	fprintf(fp, "MultiDomainMatrixZeroRows=%ld\n", g_validationMatrixZeroRows);
	fprintf(fp, "MultiDomainMatrixZeroCols=%ld\n", g_validationMatrixZeroCols);
	fprintf(fp, "MultiDomainMatrixNonFiniteBlocks=%ld\n", g_validationMatrixNonFiniteBlocks);
	fprintf(fp, "MultiDomainMatrixDiagonalBlocks=%ld\n", g_validationMatrixDiagonalBlocks);
	fprintf(fp, "MultiDomainMatrixInvertibleDiagonalBlocks=%ld\n", g_validationMatrixInvertibleDiagonalBlocks);
	fprintf(fp, "MultiDomainInterfaceDofInvariantValid=%d\n", g_validationInterfaceDofInvariantValid);
	fprintf(fp, "MultiDomainInterfaceDofInvariantFailures=%ld\n", g_validationInterfaceDofInvariantFailures);
	fprintf(fp, "MultiRateTimeMode=%d\n", g_validationMultiRateTimeMode);
	fprintf(fp, "MultiRateMacroDt=%.17g\n", g_validationMacroDt);
	fprintf(fp, "MultiRateBetaLimit=%.17g\n", g_validationBetaLimit);
	fprintf(fp, "MultiRateAveElementSize=%.17g\n", g_validationAveElementSize);
	fprintf(fp, "MultiRateMaxSubstepsPerMacro=%ld\n", g_validationMaxSubstepsPerMacro);
	for (size_t i = 0; i < g_validationDomainSubsteps.size(); ++i)
	{
		fprintf(fp,
			"MultiRateDomain[%ld]=substeps %ld stableDt %.17g localDt %.17g betaActual %.17g C1 %.17g C2 %.17g\n",
			(long)i,
			g_validationDomainSubsteps[i],
			g_validationDomainStableDt[i],
			g_validationDomainLocalDt[i],
			g_validationDomainBetaActual[i],
			g_validationDomainC1[i],
			g_validationDomainC2[i]);
	}
}

static int IsZeroBlock9(const double block9[9])
{
	for (int i = 0; i < 9; ++i)
	{
		if (block9[i] != 0.0)
			return 0;
	}
	return 1;
}

static double MultiDomainAbsTolerance()
{
	return 1.0e-14;
}

static double MultiDomainRelTolerance()
{
	return 1.0e-10;
}

static double BlockMaxAbs9(const double block9[9])
{
	double maxAbs = 0.0;
	if (!block9)
		return maxAbs;
	for (int i = 0; i < 9; ++i)
	{
		double value = fabs(block9[i]);
		if (value > maxAbs)
			maxAbs = value;
	}
	return maxAbs;
}

static double ArrayBlockMaxAbs9(const std::array<double, 9>& block)
{
	double maxAbs = 0.0;
	for (int i = 0; i < 9; ++i)
	{
		double value = fabs(block[(size_t)i]);
		if (value > maxAbs)
			maxAbs = value;
	}
	return maxAbs;
}

static int IsNearZeroBlock9(const double block9[9])
{
	return BlockMaxAbs9(block9) <= MultiDomainAbsTolerance();
}

static int IsNearZeroArrayBlock9(const std::array<double, 9>& block)
{
	return ArrayBlockMaxAbs9(block) <= MultiDomainAbsTolerance();
}

static int NearlyEqualByScale(double a, double b, double scale)
{
	double threshold = MultiDomainAbsTolerance() + MultiDomainRelTolerance() * scale;
	return fabs(a - b) <= threshold;
}

static int IsFiniteBlock9(const double block9[9])
{
	for (int i = 0; i < 9; ++i)
	{
		if (!IsFiniteScalar(block9[i]))
			return 0;
	}
	return 1;
}

static int IsZeroBlock3(const double block3[3])
{
	for (int i = 0; i < 3; ++i)
	{
		if (block3[i] != 0.0)
			return 0;
	}
	return 1;
}

static int NormalizeDomainId(int id, int domainCount, int zeroBased)
{
	if (zeroBased)
	{
		if (id >= 0 && id < domainCount)
			return id;
	}
	if (id >= 1 && id <= domainCount)
		return id - 1;
	return -1;
}

static long NormalizeElementId(long id, long eleNum, int zeroBased)
{
	if (zeroBased)
	{
		if (id >= 0 && id < eleNum)
			return id;
		return -1;
	}
	if (id >= 1 && id <= eleNum)
		return id - 1;
	return -1;
}

static long ElementLocalNodeId(long ele, int localNode)
{
	(void)localNode;
	return ele;
}

static int ClampInterfaceLocalNode(int localNode)
{
	(void)localNode;
	return 0;
}

static int InterfaceLocalBForA(const InterfacePair& itf, int localA)
{
	(void)itf;
	(void)localA;
	return 0;
}

static int InterfaceLocalAForB(const InterfacePair& itf, int localB)
{
	(void)itf;
	(void)localB;
	return 0;
}

static long InterfaceNodeA(const InterfacePair& itf, int localA)
{
	return ElementLocalNodeId(itf.eleA, ClampInterfaceLocalNode(localA));
}

static long InterfaceNodeBForA(const InterfacePair& itf, int localA)
{
	return ElementLocalNodeId(itf.eleB, InterfaceLocalBForA(itf, localA));
}

static double PointDistanceSquared(const Point& a, const Point& b)
{
	double dx = a.pt[0] - b.pt[0];
	double dy = a.pt[1] - b.pt[1];
	double dz = a.pt[2] - b.pt[2];
	return dx * dx + dy * dy + dz * dz;
}

const DynaMat& DomainMaterialContext::Get(int domainId) const
{
	if (domainId >= 0 && domainId < (int)mats.size())
		return mats[domainId];
	printf("MultiDomain material lookup failed: invalid domain id %d.\n", domainId);
	if (!mats.empty())
		return mats[0];
	printf("MultiDomain material lookup failed: material table is empty.\n");
	abort();
}

MultiDomainModel::MultiDomainModel()
	: enabled(false),
	elementCount(0),
	nodeCount(0),
	multiRateTimeMode(0),
	betaLimit(0.0),
	macroDt(0.0),
	aveElementSize(0.0),
	maxSubstepsPerMacro(64)
{
}

bool MultiDomainModel::IsActive() const
{
	return enabled;
}

int MultiDomainModel::DomainCount() const
{
	return (int)domains.size();
}

void MultiDomainModel::PrintSummary(FILE* fp) const
{
	FILE* out = fp ? fp : stdout;
	fprintf(out, "MultiDomain enabled = %d\n", enabled ? 1 : 0);
	fprintf(out, "DomainCount = %d\n", DomainCount());
	fprintf(out, "InterfaceCount = %ld\n", (long)interfaces.size());
	fprintf(out, "MultiRateTimeMode = %d macroDt=%lf betaLimit=%lf aveElementSize=%lf maxSubstepsPerMacro=%ld\n",
		multiRateTimeMode,
		macroDt,
		betaLimit,
		aveElementSize,
		maxSubstepsPerMacro);
	for (size_t i = 0; i < domains.size(); ++i)
	{
		const Domain& d = domains[i];
		fprintf(out,
			"Domain %d: eleBegin=%ld eleCount=%ld nodeBegin=%ld nodeCount=%ld boundaryNodes=%ld C1=%lf C2=%lf C1Dt=%lf C2Dt=%lf stableDt=%lf localDt=%lf substeps=%ld betaActual=%lf\n",
			d.id,
			d.eleBegin,
			d.eleCount,
			d.nodeBegin,
			d.nodeCount,
			(long)d.boundaryNodeIds.size(),
			d.mat.C1,
			d.mat.C2,
			d.mat.C1Dt,
			d.mat.C2Dt,
			d.stableDt,
			d.localDt,
			d.substeps,
			d.betaActual);
	}
	for (size_t i = 0; i < interfaces.size(); ++i)
	{
		const InterfacePair& itf = interfaces[i];
		fprintf(out,
			"Interface %ld: domainA=%d eleA=%ld domainB=%d eleB=%ld normalSign=%d\n",
			(long)i,
			itf.domainA,
			itf.eleA,
			itf.domainB,
			itf.eleB,
			itf.normalSign);
	}
}

bool MultiDomainModel::Validate() const
{
	if (domains.empty())
	{
		printf("MultiDomain validation failed: no domains.\n");
		return false;
	}
	std::vector<int> owner((size_t)elementCount, -1);
	for (size_t i = 0; i < domains.size(); ++i)
	{
		const Domain& d = domains[i];
		if (d.id != (int)i)
		{
			printf("MultiDomain validation failed: non-contiguous domain id %d at index %ld.\n", d.id, (long)i);
			return false;
		}
		if (d.eleBegin < 0 || d.eleCount <= 0 || d.eleBegin + d.eleCount > elementCount)
		{
			printf("MultiDomain validation failed: invalid element range for domain %d.\n", d.id);
			return false;
		}
		for (size_t j = 0; j < d.elementIds.size(); ++j)
		{
			long ele = d.elementIds[j];
			if (ele < 0 || ele >= elementCount)
			{
				printf("MultiDomain validation failed: element %ld out of range in domain %d.\n", ele, d.id);
				return false;
			}
			if (owner[(size_t)ele] != -1)
			{
				printf("MultiDomain validation failed: element %ld assigned twice.\n", ele);
				return false;
			}
			owner[(size_t)ele] = d.id;
		}
	}
	for (long ele = 0; ele < elementCount; ++ele)
	{
		if (owner[(size_t)ele] == -1)
		{
			printf("MultiDomain validation failed: element %ld has no domain.\n", ele);
			return false;
		}
	}
	std::vector<int> interfaceUse((size_t)elementCount, 0);
	for (size_t i = 0; i < interfaces.size(); ++i)
	{
		const InterfacePair& itf = interfaces[i];
		if (itf.domainA < 0 || itf.domainA >= DomainCount() ||
			itf.domainB < 0 || itf.domainB >= DomainCount() ||
			itf.domainA == itf.domainB)
		{
			printf("MultiDomain validation failed: invalid domains in interface %ld.\n", (long)i);
			return false;
		}
		if (itf.eleA < 0 || itf.eleA >= elementCount || itf.eleB < 0 || itf.eleB >= elementCount)
		{
			printf("MultiDomain validation failed: invalid elements in interface %ld.\n", (long)i);
			return false;
		}
		if (owner[(size_t)itf.eleA] != itf.domainA || owner[(size_t)itf.eleB] != itf.domainB)
		{
			printf("MultiDomain validation failed: interface %ld element/domain mismatch.\n", (long)i);
			return false;
		}
		interfaceUse[(size_t)itf.eleA]++;
		interfaceUse[(size_t)itf.eleB]++;
		if (interfaceUse[(size_t)itf.eleA] > 1 || interfaceUse[(size_t)itf.eleB] > 1)
		{
			printf("MultiDomain validation failed: repeated interface element in interface %ld.\n", (long)i);
			return false;
		}
	}
	return true;
}

GlobalDofMap::GlobalDofMap()
	: unknownBlockCount(0),
	knownBlockCount(0),
	equationBlockCount(0),
	outerDisplacementUnknownBlocks(0),
	outerTractionUnknownBlocks(0),
	interfaceDisplacementUnknownBlocks(0),
	interfaceTractionUnknownBlocks(0),
	historyUBlockCount(0),
	historyTBlockCount(0)
{
}

static int MarkAssignedUnknownBlock(std::vector<int>& assigned,
	long block,
	long* assignedCount,
	const char* label)
{
	if (block < 0 || block >= (long)assigned.size())
	{
		printf("MultiDomain dof map failed: %s unknown block %ld is outside node block range 0..%ld.\n",
			label, block, (long)assigned.size() - 1);
		return 0;
	}
	if (assigned[(size_t)block])
	{
		printf("MultiDomain dof map failed: duplicate %s unknown block %ld.\n", label, block);
		return 0;
	}
	assigned[(size_t)block] = 1;
	if (assignedCount)
		++(*assignedCount);
	return 1;
}

static std::string MakeOuterUnknownDebugLabel(const char* quantity,
	int domainId,
	long ele,
	int localNode,
	long nodeId,
	int bcid)
{
	char buffer[256];
	snprintf(buffer, sizeof(buffer),
		"outer %s domain=%d ele=%ld localNode=%d nodeId=%ld BCID=%d",
		quantity, domainId, ele, localNode, nodeId, bcid);
	return std::string(buffer);
}

static std::string MakeInterfaceUnknownDebugLabel(const char* quantity,
	const InterfacePair& itf,
	int localNode,
	long nodeA,
	long nodeB)
{
	char buffer[320];
	snprintf(buffer, sizeof(buffer),
		"interface %s domainA=%d eleA=%ld nodeA=%ld domainB=%d eleB=%ld nodeB=%ld localNode=%d normalSign=%d",
		quantity, itf.domainA, itf.eleA, nodeA, itf.domainB, itf.eleB, nodeB, localNode, itf.normalSign);
	return std::string(buffer);
}

static long AllocateUnknownBlock(std::vector<std::string>& labels, const std::string& label)
{
	long block = (long)labels.size();
	labels.push_back(label);
	return block;
}

bool GlobalDofMap::Build(const MultiDomainModel& model, const DSquareElement* elements, long eleNum)
{
	unknownBlockCount = 0;
	knownBlockCount = 0;
	equationBlockCount = model.nodeCount;
	outerDisplacementUnknownBlocks = 0;
	outerTractionUnknownBlocks = 0;
	interfaceDisplacementUnknownBlocks = 0;
	interfaceTractionUnknownBlocks = 0;
	historyUBlockCount = 0;
	historyTBlockCount = 0;

	m_uRefs.assign((size_t)model.nodeCount, VariableRef());
	m_tRefs.assign((size_t)model.nodeCount, VariableRef());
	m_historyURefs.assign((size_t)model.nodeCount, VariableRef());
	m_historyTRefs.assign((size_t)model.nodeCount, VariableRef());
	m_equationRows.assign((size_t)model.nodeCount, -1);
	m_preferredUnknownBlockForRow.assign((size_t)model.nodeCount, -1);
	m_unknownDebugLabels.clear();
	m_domainByElement.assign((size_t)eleNum, -1);

	for (long node = 0; node < model.nodeCount; ++node)
		m_equationRows[(size_t)node] = node;

	for (size_t i = 0; i < model.domains.size(); ++i)
	{
		const Domain& d = model.domains[i];
		for (size_t j = 0; j < d.elementIds.size(); ++j)
			m_domainByElement[(size_t)d.elementIds[j]] = d.id;
	}

	std::vector<int> interfaceSide((size_t)model.nodeCount, 0);
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		const InterfacePair& itf = model.interfaces[i];
		if (!itf.nodeMapResolved)
		{
			printf("MultiDomain dof map failed: interface %ld node map is not resolved.\n", (long)i);
			return false;
		}
		int tSignB = itf.normalSign == 0 ? -1 : itf.normalSign;
		for (int localNode = 0; localNode < 1; ++localNode)
		{
			long nodeA = InterfaceNodeA(itf, localNode);
			long nodeB = InterfaceNodeBForA(itf, localNode);
			long uBlock = AllocateUnknownBlock(m_unknownDebugLabels,
				MakeInterfaceUnknownDebugLabel("U", itf, localNode, nodeA, nodeB));
			long tBlock = AllocateUnknownBlock(m_unknownDebugLabels,
				MakeInterfaceUnknownDebugLabel("T", itf, localNode, nodeA, nodeB));

			m_uRefs[(size_t)nodeA] = VariableRef(uBlock, 1.0, false);
			m_uRefs[(size_t)nodeB] = VariableRef(uBlock, 1.0, false);
			m_tRefs[(size_t)nodeA] = VariableRef(tBlock, 1.0, false);
			m_tRefs[(size_t)nodeB] = VariableRef(tBlock, (double)tSignB, false);
			m_preferredUnknownBlockForRow[(size_t)nodeA] = uBlock;
			m_preferredUnknownBlockForRow[(size_t)nodeB] = tBlock;

			long uHistoryBlock = historyUBlockCount++;
			long tHistoryBlock = historyTBlockCount++;
			m_historyURefs[(size_t)nodeA] = VariableRef(uHistoryBlock, 1.0, false);
			m_historyURefs[(size_t)nodeB] = VariableRef(uHistoryBlock, 1.0, false);
			m_historyTRefs[(size_t)nodeA] = VariableRef(tHistoryBlock, 1.0, false);
			m_historyTRefs[(size_t)nodeB] = VariableRef(tHistoryBlock, (double)tSignB, false);

			interfaceSide[(size_t)nodeA] = 1;
			interfaceSide[(size_t)nodeB] = 1;
			++interfaceDisplacementUnknownBlocks;
			++interfaceTractionUnknownBlocks;
		}
	}

	for (long ele = 0; ele < eleNum; ++ele)
	{
		if (!IsSupportedBlockBCID(elements[ele].BCID))
		{
			printf("MultiDomain dof map failed: element %ld has mixed BCID %d. First multi-domain path supports only 123 and 456.\n",
				ele, elements[ele].BCID);
			return false;
		}

		int domainId = -1;
		if (ele >= 0 && ele < (long)m_domainByElement.size())
			domainId = m_domainByElement[(size_t)ele];

		for (int localNode = 0; localNode < 1; ++localNode)
		{
			long nodeId = ElementLocalNodeId(ele, localNode);
			if (interfaceSide[(size_t)nodeId])
				continue;

			if (elements[ele].BCID == 123)
			{
				long tBlock = AllocateUnknownBlock(m_unknownDebugLabels,
					MakeOuterUnknownDebugLabel("T", domainId, ele, localNode, nodeId, elements[ele].BCID));
				m_uRefs[(size_t)nodeId] = VariableRef(knownBlockCount++, 1.0, true);
				m_tRefs[(size_t)nodeId] = VariableRef(tBlock, 1.0, false);
				m_preferredUnknownBlockForRow[(size_t)nodeId] = tBlock;
				++outerTractionUnknownBlocks;
			}
			else
			{
				long uBlock = AllocateUnknownBlock(m_unknownDebugLabels,
					MakeOuterUnknownDebugLabel("U", domainId, ele, localNode, nodeId, elements[ele].BCID));
				m_uRefs[(size_t)nodeId] = VariableRef(uBlock, 1.0, false);
				m_tRefs[(size_t)nodeId] = VariableRef(knownBlockCount++, 1.0, true);
				m_preferredUnknownBlockForRow[(size_t)nodeId] = uBlock;
				++outerDisplacementUnknownBlocks;
			}

			m_historyURefs[(size_t)nodeId] = VariableRef(historyUBlockCount++, 1.0, false);
			m_historyTRefs[(size_t)nodeId] = VariableRef(historyTBlockCount++, 1.0, false);
		}
	}

	unknownBlockCount = (long)m_unknownDebugLabels.size();
	if (unknownBlockCount != equationBlockCount)
	{
		printf("MultiDomain dof map failed: unknown blocks %ld != equation blocks %ld.\n",
			unknownBlockCount, equationBlockCount);
		return false;
	}

	return Validate();
}

VariableRef GlobalDofMap::GetU(int domainId, long localNodeOrEle) const
{
	(void)domainId;
	if (localNodeOrEle >= 0 && localNodeOrEle < (long)m_uRefs.size())
		return m_uRefs[(size_t)localNodeOrEle];
	return VariableRef();
}

VariableRef GlobalDofMap::GetT(int domainId, long localNodeOrEle) const
{
	(void)domainId;
	if (localNodeOrEle >= 0 && localNodeOrEle < (long)m_tRefs.size())
		return m_tRefs[(size_t)localNodeOrEle];
	return VariableRef();
}

VariableRef GlobalDofMap::GetHistoryU(int domainId, long localNodeOrEle) const
{
	(void)domainId;
	if (localNodeOrEle >= 0 && localNodeOrEle < (long)m_historyURefs.size())
		return m_historyURefs[(size_t)localNodeOrEle];
	return VariableRef();
}

VariableRef GlobalDofMap::GetHistoryT(int domainId, long localNodeOrEle) const
{
	(void)domainId;
	if (localNodeOrEle >= 0 && localNodeOrEle < (long)m_historyTRefs.size())
		return m_historyTRefs[(size_t)localNodeOrEle];
	return VariableRef();
}

long GlobalDofMap::GetEquationRow(int domainId, long localNodeOrEle) const
{
	(void)domainId;
	if (localNodeOrEle >= 0 && localNodeOrEle < (long)m_equationRows.size())
		return m_equationRows[(size_t)localNodeOrEle];
	return -1;
}

long GlobalDofMap::GetPreferredUnknownBlockForRow(long row) const
{
	if (row >= 0 && row < (long)m_preferredUnknownBlockForRow.size())
		return m_preferredUnknownBlockForRow[(size_t)row];
	return -1;
}

const char* GlobalDofMap::GetUnknownDebugLabel(long unknownBlock) const
{
	if (unknownBlock >= 0 && unknownBlock < (long)m_unknownDebugLabels.size())
		return m_unknownDebugLabels[(size_t)unknownBlock].c_str();
	return "<invalid unknown block>";
}

void GlobalDofMap::PrintSummary(FILE* fp) const
{
	FILE* out = fp ? fp : stdout;
	fprintf(out, "GlobalDofMap equationBlockCount = %ld\n", equationBlockCount);
	fprintf(out, "GlobalDofMap unknownBlockCount = %ld\n", unknownBlockCount);
	fprintf(out, "GlobalDofMap knownBlockCount = %ld\n", knownBlockCount);
	fprintf(out, "Outer displacement unknown blocks = %ld\n", outerDisplacementUnknownBlocks);
	fprintf(out, "Outer traction unknown blocks = %ld\n", outerTractionUnknownBlocks);
	fprintf(out, "Interface displacement unknown blocks = %ld\n", interfaceDisplacementUnknownBlocks);
	fprintf(out, "Interface traction unknown blocks = %ld\n", interfaceTractionUnknownBlocks);
	fprintf(out, "History U blocks = %ld\n", historyUBlockCount);
	fprintf(out, "History T blocks = %ld\n", historyTBlockCount);
}

bool GlobalDofMap::Validate() const
{
	if (equationBlockCount <= 0 || unknownBlockCount <= 0)
		return false;
	if (unknownBlockCount != equationBlockCount)
		return false;
	if ((long)m_unknownDebugLabels.size() != unknownBlockCount)
		return false;
	if ((long)m_preferredUnknownBlockForRow.size() != equationBlockCount)
		return false;
	std::vector<int> unknownUsed((size_t)unknownBlockCount, 0);
	std::vector<int> preferredUsed((size_t)unknownBlockCount, 0);
	for (size_t i = 0; i < m_uRefs.size(); ++i)
	{
		if (m_equationRows[i] < 0)
			return false;
		if (m_uRefs[i].globalBlock < 0 || m_tRefs[i].globalBlock < 0)
			return false;
		if (!m_uRefs[i].known && m_uRefs[i].globalBlock >= unknownBlockCount)
			return false;
		if (!m_tRefs[i].known && m_tRefs[i].globalBlock >= unknownBlockCount)
			return false;
		if (m_uRefs[i].known && m_uRefs[i].globalBlock >= knownBlockCount)
			return false;
		if (m_tRefs[i].known && m_tRefs[i].globalBlock >= knownBlockCount)
			return false;
		if (!m_uRefs[i].known)
			unknownUsed[(size_t)m_uRefs[i].globalBlock] = 1;
		if (!m_tRefs[i].known)
			unknownUsed[(size_t)m_tRefs[i].globalBlock] = 1;
		if (m_historyURefs[i].globalBlock < 0 || m_historyURefs[i].globalBlock >= historyUBlockCount)
			return false;
		if (m_historyTRefs[i].globalBlock < 0 || m_historyTRefs[i].globalBlock >= historyTBlockCount)
			return false;
		if (m_preferredUnknownBlockForRow[i] < 0 || m_preferredUnknownBlockForRow[i] >= unknownBlockCount)
			return false;
		++preferredUsed[(size_t)m_preferredUnknownBlockForRow[i]];
	}
	for (long block = 0; block < unknownBlockCount; ++block)
	{
		if (preferredUsed[(size_t)block] != 1)
		{
			printf("MultiDomain dof map validation failed: preferred unknown block %ld is used %d times (%s).\n",
				block, preferredUsed[(size_t)block], GetUnknownDebugLabel(block));
			return false;
		}
	}
	for (long block = 0; block < unknownBlockCount; ++block)
	{
		if (!unknownUsed[(size_t)block])
		{
			printf("MultiDomain dof map validation failed: unknown block %ld is unused.\n", block);
			return false;
		}
	}
	return true;
}
static int DiagnoseInterfaceDofInvariants(const GlobalDofMap& dofMap, const MultiDomainModel& model)
{
	g_validationInterfaceDofInvariantValid = 1;
	g_validationInterfaceDofInvariantFailures = 0;
	long printed = 0;
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		const InterfacePair& itf = model.interfaces[i];
		int tSignB = itf.normalSign == 0 ? -1 : itf.normalSign;
		for (int localNode = 0; localNode < 1; ++localNode)
		{
			long nodeA = InterfaceNodeA(itf, localNode);
			long nodeB = InterfaceNodeBForA(itf, localNode);
			VariableRef uA = dofMap.GetU(itf.domainA, nodeA);
			VariableRef uB = dofMap.GetU(itf.domainB, nodeB);
			VariableRef tA = dofMap.GetT(itf.domainA, nodeA);
			VariableRef tB = dofMap.GetT(itf.domainB, nodeB);
			VariableRef huA = dofMap.GetHistoryU(itf.domainA, nodeA);
			VariableRef huB = dofMap.GetHistoryU(itf.domainB, nodeB);
			VariableRef htA = dofMap.GetHistoryT(itf.domainA, nodeA);
			VariableRef htB = dofMap.GetHistoryT(itf.domainB, nodeB);
			int ok = 1;
			ok = ok && !uA.known && !uB.known && SameBlockAndSign(uA, uB, 1.0);
			ok = ok && !tA.known && !tB.known && SameBlockAndSign(tA, tB, (double)tSignB);
			ok = ok && SameBlockAndSign(huA, huB, 1.0);
			ok = ok && SameBlockAndSign(htA, htB, (double)tSignB);
			ok = ok && dofMap.GetPreferredUnknownBlockForRow(nodeA) == uA.globalBlock;
			ok = ok && dofMap.GetPreferredUnknownBlockForRow(nodeB) == tA.globalBlock;
			if (!ok)
			{
				g_validationInterfaceDofInvariantValid = 0;
				++g_validationInterfaceDofInvariantFailures;
				if (printed < 20)
				{
					printf("MultiDomain interface dof invariant failed: interface=%ld localNode=%d nodeA=%ld nodeB=%ld uA=%ld/%g uB=%ld/%g tA=%ld/%g tB=%ld/%g huA=%ld/%g huB=%ld/%g htA=%ld/%g htB=%ld/%g preferredA=%ld preferredB=%ld normalSign=%d\n",
						(long)i,
						localNode,
						nodeA,
						nodeB,
						uA.globalBlock, uA.sign,
						uB.globalBlock, uB.sign,
						tA.globalBlock, tA.sign,
						tB.globalBlock, tB.sign,
						huA.globalBlock, huA.sign,
						huB.globalBlock, huB.sign,
						htA.globalBlock, htA.sign,
						htB.globalBlock, htB.sign,
						dofMap.GetPreferredUnknownBlockForRow(nodeA),
						dofMap.GetPreferredUnknownBlockForRow(nodeB),
						tSignB);
					++printed;
				}
			}
		}
	}
	printf("MultiDomain interface dof invariants: valid=%d failures=%ld.\n",
		g_validationInterfaceDofInvariantValid,
		g_validationInterfaceDofInvariantFailures);
	if (MultiDomainValidationOutputEnabled())
	{
		FILE* metrics = 0;
		std::string validationMetricsPath = DBEMOutputPath("validation_metrics.txt");
		fopen_s(&metrics, validationMetricsPath.c_str(), "a");
		if (metrics)
		{
			fprintf(metrics, "MultiDomainInterfaceDofInvariantValid=%d\n", g_validationInterfaceDofInvariantValid);
			fprintf(metrics, "MultiDomainInterfaceDofInvariantFailures=%ld\n", g_validationInterfaceDofInvariantFailures);
			fclose(metrics);
		}
	}
	return g_validationInterfaceDofInvariantValid;
}


MultiDomainCCSRBuilder::MultiDomainCCSRBuilder(long blockRows, long blockCols)
	: m_blockRows(blockRows), m_blockCols(blockCols)
{
	m_rhs.resize((size_t)blockRows);
	for (long i = 0; i < blockRows; ++i)
		m_rhs[(size_t)i].fill(0.0);
}

void MultiDomainCCSRBuilder::Reset(long blockRows, long blockCols)
{
	m_blockRows = blockRows;
	m_blockCols = blockCols;
	m_blocks.clear();
	m_rhs.clear();
	m_rhs.resize((size_t)blockRows);
	for (long i = 0; i < blockRows; ++i)
		m_rhs[(size_t)i].fill(0.0);
}

void MultiDomainCCSRBuilder::AddBlock(long rowBlock, long colBlock, double sign, const double block9[9])
{
	if (rowBlock < 0 || rowBlock >= m_blockRows || colBlock < 0 || colBlock >= m_blockCols)
		return;
	if (!block9 || sign == 0.0 || IsZeroBlock9(block9) || IsNearZeroBlock9(block9))
		return;
	if (!IsFiniteBlock9(block9))
	{
		static long skippedNonFiniteBlocks = 0;
		if (skippedNonFiniteBlocks < 20)
			printf("MultiDomain skipped non-finite block row=%ld col=%ld\n", rowBlock, colBlock);
		++skippedNonFiniteBlocks;
		return;
	}

	std::pair<long, long> key(rowBlock, colBlock);
	std::array<double, 9>& target = m_blocks[key];
	for (int r = 0; r < 3; ++r)
	{
		for (int c = 0; c < 3; ++c)
			target[(size_t)(r * 3 + c)] += sign * block9[c * 3 + r];
	}
	if (IsNearZeroArrayBlock9(target))
		m_blocks.erase(key);
}

void MultiDomainCCSRBuilder::AddRhsBlock(long rowBlock, double sign, const double value3[3])
{
	if (rowBlock < 0 || rowBlock >= m_blockRows)
		return;
	if (!value3 || sign == 0.0 || IsZeroBlock3(value3))
		return;
	for (int i = 0; i < 3; ++i)
		m_rhs[(size_t)rowBlock][(size_t)i] += sign * value3[i];
}

void MultiDomainCCSRBuilder::MergeFrom(const MultiDomainCCSRBuilder& other)
{
	if (other.m_blockRows != m_blockRows || other.m_blockCols != m_blockCols)
		return;
	for (std::map<std::pair<long, long>, std::array<double, 9> >::const_iterator it = other.m_blocks.begin();
		it != other.m_blocks.end();
		++it)
	{
			std::array<double, 9>& target = m_blocks[it->first];
			for (int k = 0; k < 9; ++k)
				target[(size_t)k] += it->second[(size_t)k];
			if (IsNearZeroArrayBlock9(target))
				m_blocks.erase(it->first);
		}
	for (long row = 0; row < m_blockRows; ++row)
	{
		for (int c = 0; c < 3; ++c)
			m_rhs[(size_t)row][(size_t)c] += other.m_rhs[(size_t)row][(size_t)c];
	}
}

long MultiDomainCCSRBuilder::NonZeroBlocks() const
{
	return (long)m_blocks.size();
}

void MultiDomainCCSRBuilder::Build(CCSRMat& out) const
{
	std::vector<long> rowCounts((size_t)m_blockRows, 0);
	long nnz = 0;
	for (std::map<std::pair<long, long>, std::array<double, 9> >::const_iterator it = m_blocks.begin();
		it != m_blocks.end();
		++it)
	{
		++rowCounts[(size_t)it->first.first];
		++nnz;
	}

	out.initial(m_blockRows, nnz);
	out.row_ptr[0] = 0;
	for (long row = 0; row < m_blockRows; ++row)
		out.row_ptr[row + 1] = out.row_ptr[row] + (int)rowCounts[(size_t)row];

	std::vector<long> rowCursor((size_t)m_blockRows, 0);
	for (long row = 0; row < m_blockRows; ++row)
		rowCursor[(size_t)row] = out.row_ptr[row];

	for (std::map<std::pair<long, long>, std::array<double, 9> >::const_iterator it = m_blocks.begin();
		it != m_blocks.end();
		++it)
	{
		long row = it->first.first;
		long col = it->first.second;
		long pos = rowCursor[(size_t)row]++;
		out.col_index[pos] = (int)col;
		for (int k = 0; k < 9; ++k)
			out.value[pos * 9 + k] = it->second[(size_t)k];
	}
}

void MultiDomainCCSRBuilder::BuildSym(SymCCSRMat& out, double tolerance, bool* symmetric) const
{
	bool localSymmetric = true;
	std::vector<long> rowCounts((size_t)m_blockRows, 0);
	long nnz = 0;
	for (std::map<std::pair<long, long>, std::array<double, 9> >::const_iterator it = m_blocks.begin();
		it != m_blocks.end();
		++it)
	{
		const std::array<double, 9>& b = it->second;
		double scale = ArrayBlockMaxAbs9(b);
		double threshold = tolerance;
		double mixedThreshold = MultiDomainAbsTolerance() + MultiDomainRelTolerance() * scale;
		if (mixedThreshold > threshold)
			threshold = mixedThreshold;
		if (fabs(b[1] - b[3]) > threshold ||
			fabs(b[2] - b[6]) > threshold ||
			fabs(b[5] - b[7]) > threshold)
		{
			localSymmetric = false;
		}
		++rowCounts[(size_t)it->first.first];
		++nnz;
	}

	out.initial(m_blockRows, nnz);
	out.row_ptr[0] = 0;
	for (long row = 0; row < m_blockRows; ++row)
		out.row_ptr[row + 1] = out.row_ptr[row] + (int)rowCounts[(size_t)row];

	std::vector<long> rowCursor((size_t)m_blockRows, 0);
	for (long row = 0; row < m_blockRows; ++row)
		rowCursor[(size_t)row] = out.row_ptr[row];

	for (std::map<std::pair<long, long>, std::array<double, 9> >::const_iterator it = m_blocks.begin();
		it != m_blocks.end();
		++it)
	{
		long row = it->first.first;
		long col = it->first.second;
		long pos = rowCursor[(size_t)row]++;
		const std::array<double, 9>& b = it->second;
		out.col_index[pos] = (int)col;
		out.value[pos * 6 + 0] = b[0];
		out.value[pos * 6 + 1] = b[1];
		out.value[pos * 6 + 2] = b[2];
		out.value[pos * 6 + 3] = b[4];
		out.value[pos * 6 + 4] = b[5];
		out.value[pos * 6 + 5] = b[8];
	}

	if (symmetric)
		*symmetric = localSymmetric;
}

void MultiDomainCCSRBuilder::BuildRhs(Wvector& out) const
{
	out.initial(3 * m_blockRows);
	out = 0.0;
	for (long row = 0; row < m_blockRows; ++row)
	{
		out[row * 3 + 1] = m_rhs[(size_t)row][0];
		out[row * 3 + 2] = m_rhs[(size_t)row][1];
		out[row * 3 + 3] = m_rhs[(size_t)row][2];
	}
}

void MultiDomainCCSRBuilder::BuildDense(std::vector<double>& out) const
{
	out.assign((size_t)(m_blockRows * 3 * m_blockCols * 3), 0.0);
	long denseCols = m_blockCols * 3;
	for (std::map<std::pair<long, long>, std::array<double, 9> >::const_iterator it = m_blocks.begin();
		it != m_blocks.end();
		++it)
	{
		long rowBlock = it->first.first;
		long colBlock = it->first.second;
		for (int r = 0; r < 3; ++r)
		{
			for (int c = 0; c < 3; ++c)
			{
				long denseRow = rowBlock * 3 + r;
				long denseCol = colBlock * 3 + c;
				out[(size_t)(denseRow * denseCols + denseCol)] += it->second[(size_t)(r * 3 + c)];
			}
		}
	}
}

int ReadOptionalMultiDomainInput(FILE* input, MultiDomainInputConfig& config)
{
	config = MultiDomainInputConfig();
	if (!input)
		return 0;

	int domainCount = 1;
	if (fscanf_s(input, "%d", &domainCount) != 1)
		return 0;

	if (domainCount <= 0)
	{
		printf("Invalid DomainCount %d.\n", domainCount);
		return -1;
	}

	config.enabled = true;
	config.domainCount = domainCount;
	config.materials.resize((size_t)domainCount);
	for (int i = 0; i < domainCount; ++i)
	{
		DomainMaterialInput mat;
		if (fscanf_s(input, "%d %lf %lf %lf", &mat.id, &mat.E, &mat.v, &mat.Rou) != 4)
		{
			printf("Invalid material line for domain %d.\n", i);
			return -1;
		}
		config.materials[(size_t)i] = mat;
	}

	int interfaceCount = 0;
	if (fscanf_s(input, "%d", &interfaceCount) == 1)
	{
		if (interfaceCount < 0)
		{
			printf("Invalid InterfaceCount %d.\n", interfaceCount);
			return -1;
		}
		config.hasExplicitInterfaces = true;
		config.interfaces.resize((size_t)interfaceCount);
		for (int i = 0; i < interfaceCount; ++i)
		{
			InterfacePair itf;
			if (fscanf_s(input, "%d %ld %d %ld %d",
				&itf.domainA,
				&itf.eleA,
				&itf.domainB,
				&itf.eleB,
				&itf.normalSign) != 5)
			{
				printf("Invalid interface line %d.\n", i);
				return -1;
			}
			config.interfaces[(size_t)i] = itf;
		}
	}
	return 1;
}

int BuildMultiDomainModel(const MultiDomainInputConfig& config,
	long eleNum,
	long nodeNum,
	double defaultE,
	double defaultV,
	double defaultRou,
	double& dt,
	double beta,
	double aveElementSize,
	MultiDomainModel& model)
{
	model = MultiDomainModel();
	model.enabled = config.enabled;
	model.elementCount = eleNum;
	model.nodeCount = nodeNum;
	model.betaLimit = beta;
	model.aveElementSize = aveElementSize;
	model.maxSubstepsPerMacro = 64;

	int domainCount = config.enabled ? config.domainCount : 1;
	model.domains.resize((size_t)domainCount);
	model.materialContext.mats.resize((size_t)domainCount);

	bool zeroBasedDomainIds = false;
	for (size_t i = 0; i < config.materials.size(); ++i)
	{
		if (config.materials[i].id == 0)
			zeroBasedDomainIds = true;
	}

	std::vector<DomainMaterialInput> mats((size_t)domainCount);
	std::vector<int> hasDomainMat((size_t)domainCount, 0);
	for (int i = 0; i < domainCount; ++i)
	{
		mats[(size_t)i].id = i;
		mats[(size_t)i].E = defaultE;
		mats[(size_t)i].v = defaultV;
		mats[(size_t)i].Rou = defaultRou;
	}

	if (config.enabled)
	{
		for (size_t i = 0; i < config.materials.size(); ++i)
		{
			int id = NormalizeDomainId(config.materials[i].id, domainCount, zeroBasedDomainIds ? 1 : 0);
			if (id < 0)
			{
				printf("Invalid domain material id %d.\n", config.materials[i].id);
				return 0;
			}
			mats[(size_t)id] = config.materials[i];
			mats[(size_t)id].id = id;
			hasDomainMat[(size_t)id] = 1;
		}
		for (int d = 0; d < domainCount; ++d)
		{
			if (!hasDomainMat[(size_t)d])
			{
				printf("MultiDomain time grid failed: missing material row for domain %d.\n", d);
				return 0;
			}
		}
	}

	std::vector<double> c1Values((size_t)domainCount, 0.0);
	std::vector<double> c2Values((size_t)domainCount, 0.0);
	std::vector<double> stableDtValues((size_t)domainCount, dt);
	double macroDt = dt;
	if (config.enabled)
	{
		if (beta <= 0.0 || aveElementSize <= 0.0)
		{
			printf("MultiDomain time grid failed: beta=%lf aveElementSize=%lf.\n", beta, aveElementSize);
			return 0;
		}
		macroDt = 0.0;
		for (int d = 0; d < domainCount; ++d)
		{
			double G = mats[(size_t)d].E / 2.0 / (1.0 + mats[(size_t)d].v);
			DynaMat mat = BuildDynaMat(mats[(size_t)d].v, G, mats[(size_t)d].Rou, 1.0);
			c1Values[(size_t)d] = mat.C1;
			c2Values[(size_t)d] = mat.C2;
			if (mat.C1 <= 0.0)
			{
				printf("MultiDomain time grid failed: invalid C1 for domain %d.\n", d);
				return 0;
			}
			stableDtValues[(size_t)d] = beta * aveElementSize / mat.C1;
			if (stableDtValues[(size_t)d] > macroDt)
				macroDt = stableDtValues[(size_t)d];
		}
		if (macroDt <= 0.0)
		{
			printf("MultiDomain time grid failed: invalid macroDt.\n");
			return 0;
		}
		dt = macroDt;
		model.multiRateTimeMode = 1;
		model.macroDt = macroDt;
	}
	else
	{
		model.multiRateTimeMode = 0;
		model.macroDt = dt;
	}

	long eleBegin = 0;
	for (int d = 0; d < domainCount; ++d)
	{
		long base = eleNum / domainCount;
		long extra = eleNum % domainCount;
		long eleCount = base + (d < extra ? 1 : 0);

		Domain& domain = model.domains[(size_t)d];
		domain.id = d;
		domain.eleBegin = eleBegin;
		domain.eleCount = eleCount;
		domain.nodeBegin = eleBegin;
		domain.nodeCount = eleCount;
		double G = mats[(size_t)d].E / 2.0 / (1.0 + mats[(size_t)d].v);
		double localDt = dt;
		long substeps = 1;
		if (config.enabled)
		{
			substeps = (long)ceil(model.macroDt / stableDtValues[(size_t)d] - 1.0e-12);
			if (substeps < 1)
				substeps = 1;
			if (substeps > model.maxSubstepsPerMacro)
			{
				printf("MultiDomain time grid failed: domain %d requires substeps=%ld > maxSubstepsPerMacro=%ld.\n",
					d, substeps, model.maxSubstepsPerMacro);
				return 0;
			}
			localDt = model.macroDt / (double)substeps;
		}
		domain.stableDt = config.enabled ? stableDtValues[(size_t)d] : dt;
		domain.localDt = localDt;
		domain.substeps = substeps;
		domain.c1 = config.enabled ? c1Values[(size_t)d] : 0.0;
		domain.c2 = config.enabled ? c2Values[(size_t)d] : 0.0;
		domain.betaActual = (aveElementSize > 0.0 && domain.c1 > 0.0)
			? domain.c1 * localDt / aveElementSize
			: beta;
		if (config.enabled && domain.betaActual > beta * (1.0 + 1.0e-10))
		{
			printf("MultiDomain time grid failed: domain %d betaActual=%lf exceeds betaLimit=%lf.\n",
				d, domain.betaActual, beta);
			return 0;
		}
		domain.mat = BuildDynaMat(mats[(size_t)d].v, G, mats[(size_t)d].Rou, localDt);
		if (!config.enabled)
		{
			domain.c1 = domain.mat.C1;
			domain.c2 = domain.mat.C2;
		}
		model.materialContext.mats[(size_t)d] = domain.mat;
		for (long e = 0; e < eleCount; ++e)
		{
			long id = eleBegin + e;
			domain.elementIds.push_back(id);
			domain.boundaryNodeIds.push_back(ElementLocalNodeId(id, 0));
		}
		eleBegin += eleCount;
	}

	if (config.enabled && config.hasExplicitInterfaces)
	{
		bool zeroBasedEleIds = false;
		bool zeroBasedItfDomainIds = false;
		for (size_t i = 0; i < config.interfaces.size(); ++i)
		{
			if (config.interfaces[i].eleA == 0 || config.interfaces[i].eleB == 0)
				zeroBasedEleIds = true;
			if (config.interfaces[i].domainA == 0 || config.interfaces[i].domainB == 0)
				zeroBasedItfDomainIds = true;
		}

		for (size_t i = 0; i < config.interfaces.size(); ++i)
		{
			InterfacePair itf = config.interfaces[i];
			itf.domainA = NormalizeDomainId(itf.domainA, domainCount, zeroBasedItfDomainIds ? 1 : 0);
			itf.domainB = NormalizeDomainId(itf.domainB, domainCount, zeroBasedItfDomainIds ? 1 : 0);
			itf.eleA = NormalizeElementId(itf.eleA, eleNum, zeroBasedEleIds ? 1 : 0);
			itf.eleB = NormalizeElementId(itf.eleB, eleNum, zeroBasedEleIds ? 1 : 0);
			if (itf.normalSign == 0)
				itf.normalSign = -1;
			model.interfaces.push_back(itf);
		}
	}
	else if (domainCount > 1)
	{
		for (int d = 0; d < domainCount - 1; ++d)
		{
			InterfacePair itf;
			itf.domainA = d;
			itf.domainB = d + 1;
			itf.eleA = model.domains[(size_t)d].eleBegin + model.domains[(size_t)d].eleCount - 1;
			itf.eleB = model.domains[(size_t)(d + 1)].eleBegin;
			itf.normalSign = -1;
			model.interfaces.push_back(itf);
		}
	}

	if (!model.Validate())
		return 0;

	return 1;
}

int ResolveMultiDomainInterfaceNodeMaps(MultiDomainModel& model, const Point* nodeList, long nodeNum)
{
	(void)nodeList;
	(void)nodeNum;
	if (!model.IsActive() || model.interfaces.empty())
		return 1;
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		InterfacePair& itf = model.interfaces[i];
		itf.localBForA[0] = 0;
		itf.localAForB[0] = 0;
		itf.nodeMapResolved = 1;
	}
	printf("MultiDomain constant-element interface maps resolved: interfaces=%ld.\n",
		(long)model.interfaces.size());
	return 1;
}

void ApplyMultiDomainModelToElements(const MultiDomainModel& model, DSquareElement* elements, long eleNum)
{
	if (!elements)
		return;
	for (long i = 0; i < eleNum; ++i)
	{
		elements[i].DomainID = 0;
		elements[i].LocalEleID = (int)i;
		elements[i].SurfaceType = SurfaceOuter;
	}
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		for (size_t j = 0; j < domain.elementIds.size(); ++j)
		{
			long ele = domain.elementIds[j];
			if (ele >= 0 && ele < eleNum)
			{
				elements[ele].DomainID = domain.id;
				elements[ele].LocalEleID = (int)j;
			}
		}
	}
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		const InterfacePair& itf = model.interfaces[i];
		if (itf.eleA >= 0 && itf.eleA < eleNum)
			elements[itf.eleA].SurfaceType = SurfaceInterface;
		if (itf.eleB >= 0 && itf.eleB < eleNum)
			elements[itf.eleB].SurfaceType = SurfaceInterface;
	}
}

class ScopedDynaMat {
public:
	explicit ScopedDynaMat(const DynaMat& mat)
	{
		m_saved = DSquareElement::SetThreadDynaMat(&mat);
	}
	~ScopedDynaMat()
	{
		DSquareElement::SetThreadDynaMat(m_saved);
	}

private:
	const DynaMat* m_saved;
};

class ScopedBCID {
public:
	ScopedBCID(DSquareElement& element, int bcid)
		: m_element(element), m_saved(element.BCID)
	{
		m_element.BCID = bcid;
	}
	~ScopedBCID()
	{
		m_element.BCID = m_saved;
	}

private:
	DSquareElement& m_element;
	int m_saved;
};

static const DynaMat* GetDomainMatOrDefault(const MultiDomainModel& model, int domainId)
{
	if (domainId >= 0 && domainId < (int)model.domains.size())
		return &model.domains[(size_t)domainId].mat;
	printf("MultiDomain assembly failed: invalid domain id %d.\n", domainId);
	return 0;
}

static void ZeroBlock9(double block9[9])
{
	for (int i = 0; i < 9; ++i)
		block9[i] = 0.0;
}

static void CopyBlock9(const double src[9], double dst[9])
{
	for (int i = 0; i < 9; ++i)
		dst[i] = src[i];
}

static void IdentityTransform3(double transform[9])
{
	for (int i = 0; i < 9; ++i)
		transform[i] = 0.0;
	transform[0] = 1.0;
	transform[4] = 1.0;
	transform[8] = 1.0;
}

static void BuildLocalTransformFromGlobal(const double* sideTrans, double transform[9])
{
	for (int localComp = 0; localComp < 3; ++localComp)
	{
		for (int globalComp = 0; globalComp < 3; ++globalComp)
			transform[localComp * 3 + globalComp] = sideTrans[localComp * 3 + globalComp];
	}
}

static void BuildLocalTransformFromReference(const double* sideTrans,
	const double* referenceTrans,
	double transform[9])
{
	for (int sideComp = 0; sideComp < 3; ++sideComp)
	{
		for (int refComp = 0; refComp < 3; ++refComp)
		{
			double dot = 0.0;
			for (int axis = 0; axis < 3; ++axis)
				dot += sideTrans[sideComp * 3 + axis] * referenceTrans[refComp * 3 + axis];
			transform[sideComp * 3 + refComp] = dot;
		}
	}
}

static int GetElementLocalTransformFromInterfaceReference(const MultiDomainModel& model,
	int domainId,
	long ele,
	int localNode,
	double transform[9])
{
	IdentityTransform3(transform);
	if (!DSquareElement::m_transmat || localNode != 0)
		return 1;
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		const InterfacePair& itf = model.interfaces[i];
		if (itf.domainA == domainId && itf.eleA == ele)
			return 1;
		if (itf.domainB == domainId && itf.eleB == ele)
		{
			int localA = InterfaceLocalAForB(itf, localNode);
			long nodeB = ElementLocalNodeId(itf.eleB, localNode);
			long nodeA = ElementLocalNodeId(itf.eleA, localA);
			BuildLocalTransformFromReference(DSquareElement::m_transmat[nodeB],
				DSquareElement::m_transmat[nodeA],
				transform);
			return 1;
		}
	}
	return 1;
}

static int IsInterfaceReferenceSide(const MultiDomainModel& model, int domainId, long ele)
{
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		const InterfacePair& itf = model.interfaces[i];
		if (itf.domainA == domainId && itf.eleA == ele)
			return 1;
		if (itf.domainB == domainId && itf.eleB == ele)
			return 0;
	}
	return 1;
}

static void ApplyTransform3(const double transform[9], const double src[3], double dst[3])
{
	for (int r = 0; r < 3; ++r)
		dst[r] = transform[r * 3 + 0] * src[0] +
			transform[r * 3 + 1] * src[1] +
			transform[r * 3 + 2] * src[2];
}

static void TransformColumnBlock9(const double block9[9], const double transform[9], double out9[9])
{
	for (int r = 0; r < 3; ++r)
	{
		for (int q = 0; q < 3; ++q)
		{
			double value = 0.0;
			for (int c = 0; c < 3; ++c)
				value += block9[c * 3 + r] * transform[c * 3 + q];
			out9[q * 3 + r] = value;
		}
	}
}

static void AddBlockWithTransform(MultiDomainCCSRBuilder& builder,
	long row,
	const VariableRef& ref,
	const double transform[9],
	const double block9[9])
{
	double transformed[9];
	TransformColumnBlock9(block9, transform, transformed);
	builder.AddBlock(row, ref.globalBlock, ref.sign, transformed);
}

static int ComputeUnknown0Block(DSquareElement* elements,
	const MultiDomainModel& model,
	int domainId,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long** infElePid,
	long** elePid,
	double block9[9])
{
	ZeroBlock9(block9);
	const DynaMat* mat = GetDomainMatOrDefault(model, domainId);
	if (!mat)
		return 0;
	ScopedDynaMat matScope(*mat);
	long threadId = 0;
	int flag = UnknownSubMatrix(elements, sourceNode, fieldEle, infElePid, elePid, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistCoef[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeKnown0Block(DSquareElement* elements,
	const MultiDomainModel& model,
	int domainId,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long** infElePid,
	long** elePid,
	double block9[9])
{
	ZeroBlock9(block9);
	const DynaMat* mat = GetDomainMatOrDefault(model, domainId);
	if (!mat)
		return 0;
	ScopedDynaMat matScope(*mat);
	long threadId = 0;
	int flag = knownSubMatrix(elements, sourceNode, fieldEle, infElePid, elePid, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistCoef[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeT1Block(DSquareElement* elements,
	const MultiDomainModel& model,
	int domainId,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long step,
	double block9[9])
{
	ZeroBlock9(block9);
	const DynaMat* mat = GetDomainMatOrDefault(model, domainId);
	if (!mat)
		return 0;
	ScopedDynaMat matScope(*mat);
	Point& source = elements[fieldEle].m_nodelist[sourceNode];
	long threadId = 0;
	int flag = elements[fieldEle].IntDynaTijJudge(1, source, (double)step, DSquareElement::ActiveDynaMat().Dt, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistTij[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeT2Block(DSquareElement* elements,
	const MultiDomainModel& model,
	int domainId,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long step,
	double block9[9])
{
	ZeroBlock9(block9);
	const DynaMat* mat = GetDomainMatOrDefault(model, domainId);
	if (!mat)
		return 0;
	ScopedDynaMat matScope(*mat);
	int pos = -1;
	if (step <= 1 && elements[fieldEle].IsIn(sourceNode, pos))
	{
		CopyBlock9(elements[fieldEle].m_OnElementT2ij[pos][fieldLocalNode], block9);
		return 1;
	}

	Point& source = elements[fieldEle].m_nodelist[sourceNode];
	long threadId = 0;
	int flag = elements[fieldEle].IntDynaTijJudge(2, source, (double)step, DSquareElement::ActiveDynaMat().Dt, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistTij[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeUBlock(DSquareElement* elements,
	const MultiDomainModel& model,
	int domainId,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long step,
	double block9[9])
{
	ZeroBlock9(block9);
	const DynaMat* mat = GetDomainMatOrDefault(model, domainId);
	if (!mat)
		return 0;
	ScopedDynaMat matScope(*mat);
	Point& source = elements[fieldEle].m_nodelist[sourceNode];
	long threadId = 0;
	int flag = elements[fieldEle].IntDynaUijJudge(source, (double)step, DSquareElement::ActiveDynaMat().Dt, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistUij[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeUnknown0BlockCurrentMat(DSquareElement* elements,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long** infElePid,
	long** elePid,
	long threadId,
	double block9[9])
{
	ZeroBlock9(block9);
	int flag = UnknownSubMatrix(elements, sourceNode, fieldEle, infElePid, elePid, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistCoef[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeKnown0BlockCurrentMat(DSquareElement* elements,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long** infElePid,
	long** elePid,
	long threadId,
	double block9[9])
{
	ZeroBlock9(block9);
	int flag = knownSubMatrix(elements, sourceNode, fieldEle, infElePid, elePid, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistCoef[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeT1BlockCurrentMat(DSquareElement* elements,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long step,
	long threadId,
	double block9[9])
{
	ZeroBlock9(block9);
	Point& source = elements[fieldEle].m_nodelist[sourceNode];
	int flag = elements[fieldEle].IntDynaTijJudge(1, source, (double)step, DSquareElement::ActiveDynaMat().Dt, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistTij[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeT2BlockCurrentMat(DSquareElement* elements,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long step,
	long threadId,
	double block9[9])
{
	ZeroBlock9(block9);
	int pos = -1;
	if (step <= 1 && elements[fieldEle].IsIn(sourceNode, pos))
	{
		CopyBlock9(elements[fieldEle].m_OnElementT2ij[pos][fieldLocalNode], block9);
		return 1;
	}

	Point& source = elements[fieldEle].m_nodelist[sourceNode];
	int flag = elements[fieldEle].IntDynaTijJudge(2, source, (double)step, DSquareElement::ActiveDynaMat().Dt, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistTij[threadId][fieldLocalNode], block9);
	return 1;
}

static int ComputeUBlockCurrentMat(DSquareElement* elements,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long step,
	long threadId,
	double block9[9])
{
	ZeroBlock9(block9);
	Point& source = elements[fieldEle].m_nodelist[sourceNode];
	int flag = elements[fieldEle].IntDynaUijJudge(source, (double)step, DSquareElement::ActiveDynaMat().Dt, threadId);
	if (!flag)
		return 0;
	CopyBlock9(AssistUij[threadId][fieldLocalNode], block9);
	return 1;
}

static int AddKnownOrUnknownBlockCurrentMat(DSquareElement* elements,
	const GlobalDofMap& dofMap,
	long row,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long fieldNode,
	int domainId,
	int unknownPass,
	MultiDomainCCSRBuilder& builder,
	long** infElePid,
	long** elePid,
	long threadId)
{
	double block9[9];
	int flag = unknownPass
		? ComputeUnknown0BlockCurrentMat(elements, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, threadId, block9)
		: ComputeKnown0BlockCurrentMat(elements, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, threadId, block9);
	if (!flag)
		return 0;

	VariableRef ref;
	if (elements[fieldEle].BCID == 123)
		ref = unknownPass ? dofMap.GetT(domainId, fieldNode) : dofMap.GetU(domainId, fieldNode);
	else
		ref = unknownPass ? dofMap.GetU(domainId, fieldNode) : dofMap.GetT(domainId, fieldNode);

	if (ref.known == (unknownPass ? false : true))
		builder.AddBlock(row, ref.globalBlock, ref.sign, block9);
	return 1;
}

struct MultiDomainStep0ThreadArg {
	DSquareElement* elements;
	const Domain* domain;
	const GlobalDofMap* dofMap;
	long** infElePid;
	long** elePid;
	MultiDomainCCSRBuilder* unknownBuilder;
	MultiDomainCCSRBuilder* knownBuilder;
	long sourceBegin;
	long sourceEnd;
	long threadId;
	long blockCount;
	int ok;
};

struct MultiDomainHistoryThreadArg {
	DSquareElement* elements;
	const Domain* domain;
	const GlobalDofMap* dofMap;
	long step;
	MultiDomainCCSRBuilder* mtsBuilder;
	MultiDomainCCSRBuilder* mgsBuilder;
	long sourceBegin;
	long sourceEnd;
	long threadId;
	long blockCount;
	int ok;
};

static void* MultiDomainStep0ThreadMain(void* rawArg)
{
	MultiDomainStep0ThreadArg* arg = (MultiDomainStep0ThreadArg*)rawArg;
	if (!arg)
		return 0;
	arg->ok = 1;
	arg->blockCount = 0;
	if (!arg || !arg->elements || !arg->domain || !arg->dofMap ||
		!arg->unknownBuilder || !arg->knownBuilder)
	{
		arg->ok = 0;
		return 0;
	}

	DSquareElement* elements = arg->elements;
	const Domain& domain = *arg->domain;
	ScopedDynaMat matScope(domain.mat);
	const GlobalDofMap& dofMap = *arg->dofMap;
	for (long si = arg->sourceBegin; si < arg->sourceEnd; ++si)
	{
		long sourceEle = domain.elementIds[(size_t)si];
		for (int sourceLocalNode = 0; sourceLocalNode < 1; ++sourceLocalNode)
		{
			long sourceNode = elements[sourceEle].m_nodeID[sourceLocalNode];
			long row = dofMap.GetEquationRow(domain.id, sourceNode);
			if (row < 0)
				continue;

			for (size_t fj = 0; fj < domain.elementIds.size(); ++fj)
			{
				long fieldEle = domain.elementIds[fj];
				if (elements[fieldEle].SurfaceType == SurfaceInterface)
					continue;
				for (int fieldLocalNode = 0; fieldLocalNode < 1; ++fieldLocalNode)
				{
					long fieldNode = elements[fieldEle].m_nodeID[fieldLocalNode];
					arg->blockCount += AddKnownOrUnknownBlockCurrentMat(elements, dofMap, row, sourceNode,
						fieldEle, fieldLocalNode, fieldNode, domain.id, 1,
						*arg->unknownBuilder, arg->infElePid, arg->elePid, arg->threadId);
					arg->blockCount += AddKnownOrUnknownBlockCurrentMat(elements, dofMap, row, sourceNode,
						fieldEle, fieldLocalNode, fieldNode, domain.id, 0,
						*arg->knownBuilder, arg->infElePid, arg->elePid, arg->threadId);
				}
			}
		}
	}
	return 0;
}

static void* MultiDomainHistoryThreadMain(void* rawArg)
{
	MultiDomainHistoryThreadArg* arg = (MultiDomainHistoryThreadArg*)rawArg;
	if (!arg)
		return 0;
	arg->ok = 1;
	arg->blockCount = 0;
	if (!arg || !arg->elements || !arg->domain || !arg->dofMap ||
		!arg->mtsBuilder || !arg->mgsBuilder)
	{
		arg->ok = 0;
		return 0;
	}

	DSquareElement* elements = arg->elements;
	const Domain& domain = *arg->domain;
	ScopedDynaMat matScope(domain.mat);
	const GlobalDofMap& dofMap = *arg->dofMap;
	for (long si = arg->sourceBegin; si < arg->sourceEnd; ++si)
	{
		long sourceEle = domain.elementIds[(size_t)si];
		for (int sourceLocalNode = 0; sourceLocalNode < 1; ++sourceLocalNode)
		{
			long sourceNode = elements[sourceEle].m_nodeID[sourceLocalNode];
			long row = dofMap.GetEquationRow(domain.id, sourceNode);
			if (row < 0)
				continue;

			for (size_t fj = 0; fj < domain.elementIds.size(); ++fj)
			{
				long fieldEle = domain.elementIds[fj];
				double block9[9];
				for (int fieldLocalNode = 0; fieldLocalNode < 1; ++fieldLocalNode)
				{
					long fieldNode = elements[fieldEle].m_nodeID[fieldLocalNode];
					VariableRef uRef = dofMap.GetHistoryU(domain.id, fieldNode);
					if (ComputeT2BlockCurrentMat(elements, sourceNode, fieldEle, fieldLocalNode,
						arg->step, arg->threadId, block9))
					{
						arg->mtsBuilder->AddBlock(row, uRef.globalBlock, uRef.sign, block9);
						++arg->blockCount;
					}
					if (ComputeT1BlockCurrentMat(elements, sourceNode, fieldEle, fieldLocalNode,
						arg->step, arg->threadId, block9))
					{
						arg->mtsBuilder->AddBlock(row, uRef.globalBlock, uRef.sign, block9);
						++arg->blockCount;
					}

					VariableRef tRef = dofMap.GetHistoryT(domain.id, fieldNode);
					if (ComputeUBlockCurrentMat(elements, sourceNode, fieldEle, fieldLocalNode,
						arg->step, arg->threadId, block9))
					{
						arg->mgsBuilder->AddBlock(row, tRef.globalBlock, tRef.sign, block9);
						++arg->blockCount;
					}
				}
			}
		}
	}
	return 0;
}

static long EffectiveMultiDomainThreadCount(const Domain& domain, long requestedThreads)
{
	long count = requestedThreads;
	if (count <= 0)
		count = thread_num;
	if (count <= 0)
		count = 1;
	if (count > thread_num && thread_num > 0)
		count = thread_num;
	long sourceCount = (long)domain.elementIds.size();
	if (sourceCount > 0 && count > sourceCount)
		count = sourceCount;
	if (count < 1)
		count = 1;
	return count;
}

static void SetThreadSourceRange(long sourceCount, long threadId, long threadCount, long& begin, long& end)
{
	long size = sourceCount / threadCount;
	long rest = sourceCount % threadCount;
	begin = threadId * size + (threadId < rest ? threadId : rest);
	end = begin + size + (threadId < rest ? 1 : 0);
}

static int RunPthreadWorkers(pthread_t* threads, int* created, long threadCount)
{
	for (long i = 0; i < threadCount; ++i)
		created[i] = 0;
	return 1;
}

static double MaxDenseDifference(const MultiDomainCCSRBuilder& a, const MultiDomainCCSRBuilder& b)
{
	if (a.BlockRows() != b.BlockRows() || a.BlockCols() != b.BlockCols())
		return (std::numeric_limits<double>::max)();
	std::vector<double> denseA;
	std::vector<double> denseB;
	a.BuildDense(denseA);
	b.BuildDense(denseB);
	if (denseA.size() != denseB.size())
		return (std::numeric_limits<double>::max)();
	double maxDiff = 0.0;
	for (size_t i = 0; i < denseA.size(); ++i)
	{
		double diff = fabs(denseA[i] - denseB[i]);
		if (diff > maxDiff)
			maxDiff = diff;
	}
	return maxDiff;
}

static double MaxDenseMagnitude(const MultiDomainCCSRBuilder& a, const MultiDomainCCSRBuilder& b)
{
	std::vector<double> denseA;
	std::vector<double> denseB;
	a.BuildDense(denseA);
	b.BuildDense(denseB);
	double maxAbs = 0.0;
	for (size_t i = 0; i < denseA.size(); ++i)
	{
		double value = fabs(denseA[i]);
		if (value > maxAbs)
			maxAbs = value;
	}
	for (size_t i = 0; i < denseB.size(); ++i)
	{
		double value = fabs(denseB[i]);
		if (value > maxAbs)
			maxAbs = value;
	}
	return maxAbs;
}

static double BuilderCompareTolerance(const MultiDomainCCSRBuilder& a, const MultiDomainCCSRBuilder& b)
{
	return MultiDomainAbsTolerance() + MultiDomainRelTolerance() * MaxDenseMagnitude(a, b);
}

static int CompareOrFallbackStep0(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long** infElePid,
	long** elePid,
	MultiDomainCCSRBuilder& unknownBuilder,
	MultiDomainCCSRBuilder& knownBuilder)
{
	if (!MultiDomainPthreadCheckEnabled())
		return 1;
	MultiDomainCCSRBuilder serialUnknown(dofMap.equationBlockCount, dofMap.unknownBlockCount);
	MultiDomainCCSRBuilder serialKnown(dofMap.equationBlockCount, dofMap.knownBlockCount);
	if (!AssembleMultiDomainStep0(elements, model, dofMap, infElePid, elePid, serialUnknown, serialKnown))
		return 0;
	double unknownDiff = MaxDenseDifference(unknownBuilder, serialUnknown);
	double knownDiff = MaxDenseDifference(knownBuilder, serialKnown);
	double unknownTol = BuilderCompareTolerance(unknownBuilder, serialUnknown);
	double knownTol = BuilderCompareTolerance(knownBuilder, serialKnown);
	printf("MultiDomain pthread check Step0: unknownBlocks pthread=%ld serial=%ld maxDiff=%e; knownBlocks pthread=%ld serial=%ld maxDiff=%e\n",
		unknownBuilder.NonZeroBlocks(), serialUnknown.NonZeroBlocks(), unknownDiff,
		knownBuilder.NonZeroBlocks(), serialKnown.NonZeroBlocks(), knownDiff);
	if (unknownBuilder.NonZeroBlocks() != serialUnknown.NonZeroBlocks() ||
		knownBuilder.NonZeroBlocks() != serialKnown.NonZeroBlocks() ||
		unknownDiff > unknownTol || knownDiff > knownTol)
	{
		printf("MultiDomain pthread Step0 differs from serial assembly; falling back to serial matrices.\n");
		unknownBuilder = serialUnknown;
		knownBuilder = serialKnown;
	}
	return 1;
}

static int CompareOrFallbackHistory(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long step,
	MultiDomainCCSRBuilder& mtsBuilder,
	MultiDomainCCSRBuilder& mgsBuilder)
{
	if (!MultiDomainPthreadCheckEnabled())
		return 1;
	MultiDomainCCSRBuilder serialMts(dofMap.equationBlockCount, dofMap.historyUBlockCount);
	MultiDomainCCSRBuilder serialMgs(dofMap.equationBlockCount, dofMap.historyTBlockCount);
	if (!AssembleMultiDomainHistoryStep(elements, model, dofMap, step, serialMts, serialMgs))
		return 0;
	double mtsDiff = MaxDenseDifference(mtsBuilder, serialMts);
	double mgsDiff = MaxDenseDifference(mgsBuilder, serialMgs);
	double mtsTol = BuilderCompareTolerance(mtsBuilder, serialMts);
	double mgsTol = BuilderCompareTolerance(mgsBuilder, serialMgs);
	printf("MultiDomain pthread check History step=%ld: MTS pthread=%ld serial=%ld maxDiff=%e; MGS pthread=%ld serial=%ld maxDiff=%e\n",
		step,
		mtsBuilder.NonZeroBlocks(), serialMts.NonZeroBlocks(), mtsDiff,
		mgsBuilder.NonZeroBlocks(), serialMgs.NonZeroBlocks(), mgsDiff);
	if (mtsBuilder.NonZeroBlocks() != serialMts.NonZeroBlocks() ||
		mgsBuilder.NonZeroBlocks() != serialMgs.NonZeroBlocks() ||
		mtsDiff > mtsTol || mgsDiff > mgsTol)
	{
		printf("MultiDomain pthread History step=%ld differs from serial assembly; falling back to serial matrices.\n", step);
		mtsBuilder = serialMts;
		mgsBuilder = serialMgs;
	}
	return 1;
}


static int AddKnownOrUnknownBlock(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long row,
	long sourceNode,
	long fieldEle,
	int fieldLocalNode,
	long fieldNode,
	int domainId,
	int unknownPass,
	MultiDomainCCSRBuilder& builder,
	long** infElePid,
	long** elePid)
{
	double block9[9];
	int flag = unknownPass
		? ComputeUnknown0Block(elements, model, domainId, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, block9)
		: ComputeKnown0Block(elements, model, domainId, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, block9);
	if (!flag)
		return 0;

	VariableRef ref;
	if (elements[fieldEle].BCID == 123)
		ref = unknownPass ? dofMap.GetT(domainId, fieldNode) : dofMap.GetU(domainId, fieldNode);
	else
		ref = unknownPass ? dofMap.GetU(domainId, fieldNode) : dofMap.GetT(domainId, fieldNode);

	if (ref.known == (unknownPass ? false : true))
		builder.AddBlock(row, ref.globalBlock, ref.sign, block9);
	return 1;
}

int AssembleMultiDomainStep0(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long** infElePid,
	long** elePid,
	MultiDomainCCSRBuilder& unknownBuilder,
	MultiDomainCCSRBuilder& knownBuilder)
{
	if (!elements || !model.Validate() || !dofMap.Validate())
		return 0;

	unknownBuilder.Reset(dofMap.equationBlockCount, dofMap.unknownBlockCount);
	knownBuilder.Reset(dofMap.equationBlockCount, dofMap.knownBlockCount);

	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		for (size_t si = 0; si < domain.elementIds.size(); ++si)
		{
			long sourceEle = domain.elementIds[si];
			for (int sourceLocalNode = 0; sourceLocalNode < 1; ++sourceLocalNode)
			{
				long sourceNode = elements[sourceEle].m_nodeID[sourceLocalNode];
				long row = dofMap.GetEquationRow(domain.id, sourceNode);
				if (row < 0)
					continue;

				for (size_t fj = 0; fj < domain.elementIds.size(); ++fj)
				{
					long fieldEle = domain.elementIds[fj];
					int savedBCID = elements[fieldEle].BCID;
					double block9[9];

					for (int fieldLocalNode = 0; fieldLocalNode < 1; ++fieldLocalNode)
					{
						long fieldNode = elements[fieldEle].m_nodeID[fieldLocalNode];
						if (elements[fieldEle].SurfaceType == SurfaceInterface)
						{
							double transform[9];
							GetElementLocalTransformFromInterfaceReference(model, domain.id, fieldEle, fieldLocalNode, transform);
							{
								ScopedBCID bcidScope(elements[fieldEle], 456);
								if (ComputeUnknown0Block(elements, model, domain.id, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, block9))
								{
									VariableRef uRef = dofMap.GetU(domain.id, fieldNode);
									AddBlockWithTransform(unknownBuilder, row, uRef, transform, block9);
								}
							}

							{
								ScopedBCID bcidScope(elements[fieldEle], 123);
								if (ComputeUnknown0Block(elements, model, domain.id, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, block9))
								{
									VariableRef tRef = dofMap.GetT(domain.id, fieldNode);
									AddBlockWithTransform(unknownBuilder, row, tRef, transform, block9);
								}
							}
						}
						else
						{
							AddKnownOrUnknownBlock(elements, model, dofMap, row, sourceNode, fieldEle, fieldLocalNode, fieldNode, domain.id, 1, unknownBuilder, infElePid, elePid);
							AddKnownOrUnknownBlock(elements, model, dofMap, row, sourceNode, fieldEle, fieldLocalNode, fieldNode, domain.id, 0, knownBuilder, infElePid, elePid);
						}
					}

					elements[fieldEle].BCID = savedBCID;
				}
			}
		}
	}
	return 1;
}

static long MultiRateLocalHistoryStep(const Domain& domain, long macroStep)
{
	long substeps = domain.substeps > 0 ? domain.substeps : 1;
	return macroStep * substeps;
}

int AssembleMultiDomainHistoryStep(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long step,
	MultiDomainCCSRBuilder& mtsBuilder,
	MultiDomainCCSRBuilder& mgsBuilder)
{
	if (!elements || step <= 0 || !model.Validate() || !dofMap.Validate())
		return 0;

	mtsBuilder.Reset(dofMap.equationBlockCount, dofMap.historyUBlockCount);
	mgsBuilder.Reset(dofMap.equationBlockCount, dofMap.historyTBlockCount);

	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		long localStep = MultiRateLocalHistoryStep(domain, step);
		for (size_t si = 0; si < domain.elementIds.size(); ++si)
		{
			long sourceEle = domain.elementIds[si];
			for (int sourceLocalNode = 0; sourceLocalNode < 1; ++sourceLocalNode)
			{
				long sourceNode = elements[sourceEle].m_nodeID[sourceLocalNode];
				long row = dofMap.GetEquationRow(domain.id, sourceNode);
				if (row < 0)
					continue;

				for (size_t fj = 0; fj < domain.elementIds.size(); ++fj)
				{
					long fieldEle = domain.elementIds[fj];
					double block9[9];
					for (int fieldLocalNode = 0; fieldLocalNode < 1; ++fieldLocalNode)
					{
						long fieldNode = elements[fieldEle].m_nodeID[fieldLocalNode];
						VariableRef uRef = dofMap.GetHistoryU(domain.id, fieldNode);
						if (ComputeT2Block(elements, model, domain.id, sourceNode, fieldEle, fieldLocalNode, localStep, block9))
							mtsBuilder.AddBlock(row, uRef.globalBlock, uRef.sign, block9);
						if (ComputeT1Block(elements, model, domain.id, sourceNode, fieldEle, fieldLocalNode, localStep, block9))
							mtsBuilder.AddBlock(row, uRef.globalBlock, uRef.sign, block9);

						VariableRef tRef = dofMap.GetHistoryT(domain.id, fieldNode);
						if (ComputeUBlock(elements, model, domain.id, sourceNode, fieldEle, fieldLocalNode, localStep, block9))
							mgsBuilder.AddBlock(row, tRef.globalBlock, tRef.sign, block9);
					}
				}
			}
		}
	}

	return 1;
}

static long AssembleMultiDomainStep0InterfaceSerialForDomain(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	const Domain& domain,
	long** infElePid,
	long** elePid,
	MultiDomainCCSRBuilder& unknownBuilder)
{
	long interfaceBlocks = 0;
	for (size_t si = 0; si < domain.elementIds.size(); ++si)
	{
		long sourceEle = domain.elementIds[si];
		for (int sourceLocalNode = 0; sourceLocalNode < 1; ++sourceLocalNode)
		{
			long sourceNode = elements[sourceEle].m_nodeID[sourceLocalNode];
			long row = dofMap.GetEquationRow(domain.id, sourceNode);
			if (row < 0)
				continue;

			for (size_t fj = 0; fj < domain.elementIds.size(); ++fj)
			{
				long fieldEle = domain.elementIds[fj];
				if (elements[fieldEle].SurfaceType != SurfaceInterface)
					continue;
				int savedBCID = elements[fieldEle].BCID;
				double block9[9];

				for (int fieldLocalNode = 0; fieldLocalNode < 1; ++fieldLocalNode)
				{
					long fieldNode = elements[fieldEle].m_nodeID[fieldLocalNode];
					double transform[9];
					GetElementLocalTransformFromInterfaceReference(model, domain.id, fieldEle, fieldLocalNode, transform);
					{
						ScopedBCID bcidScope(elements[fieldEle], 456);
						if (ComputeUnknown0Block(elements, model, domain.id, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, block9))
						{
							VariableRef uRef = dofMap.GetU(domain.id, fieldNode);
							AddBlockWithTransform(unknownBuilder, row, uRef, transform, block9);
							++interfaceBlocks;
						}
					}

					{
						ScopedBCID bcidScope(elements[fieldEle], 123);
						if (ComputeUnknown0Block(elements, model, domain.id, sourceNode, fieldEle, fieldLocalNode, infElePid, elePid, block9))
						{
							VariableRef tRef = dofMap.GetT(domain.id, fieldNode);
							AddBlockWithTransform(unknownBuilder, row, tRef, transform, block9);
							++interfaceBlocks;
						}
					}
				}

				elements[fieldEle].BCID = savedBCID;
			}
		}
	}
	return interfaceBlocks;
}

static int AssembleMultiDomainStep0DomainPthread(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	const Domain& domain,
	long requestedThreads,
	long** infElePid,
	long** elePid,
	MultiDomainCCSRBuilder& unknownBuilder,
	MultiDomainCCSRBuilder& knownBuilder,
	long& interfaceBlocks)
{
	interfaceBlocks = 0;
	long sourceCount = (long)domain.elementIds.size();
	if (sourceCount <= 0)
		return 1;

	long localThreadCount = EffectiveMultiDomainThreadCount(domain, requestedThreads);
	clock_t domainClock = clock();

	if (localThreadCount <= 1)
	{
		MultiDomainStep0ThreadArg arg;
		arg.elements = elements;
		arg.domain = &domain;
		arg.dofMap = &dofMap;
		arg.infElePid = infElePid;
		arg.elePid = elePid;
		arg.unknownBuilder = &unknownBuilder;
		arg.knownBuilder = &knownBuilder;
		arg.sourceBegin = 0;
		arg.sourceEnd = sourceCount;
		arg.threadId = 0;
		arg.blockCount = 0;
		arg.ok = 1;
			MultiDomainStep0ThreadMain(&arg);
			interfaceBlocks = AssembleMultiDomainStep0InterfaceSerialForDomain(elements, model, dofMap,
				domain, infElePid, elePid, unknownBuilder);
			printf("MultiDomain Step0 pthread assembly domain=%d threads=1 nonInterfaceBlocks=%ld interfaceSerialBlocks=%ld time=%lf\n",
			domain.id, arg.blockCount, interfaceBlocks,
			(double)(clock() - domainClock) / (double)CLOCKS_PER_SEC);
		return arg.ok;
	}

	pthread_t* threads = new pthread_t[(size_t)localThreadCount];
	int* created = new int[(size_t)localThreadCount];
	MultiDomainStep0ThreadArg* args = new MultiDomainStep0ThreadArg[(size_t)localThreadCount];
	MultiDomainCCSRBuilder* localUnknown = new MultiDomainCCSRBuilder[(size_t)localThreadCount];
	MultiDomainCCSRBuilder* localKnown = new MultiDomainCCSRBuilder[(size_t)localThreadCount];
	RunPthreadWorkers(threads, created, localThreadCount);

	int createOk = 1;
	long createSuccess = 0;
	for (long t = 0; t < localThreadCount; ++t)
	{
		long begin = 0;
		long end = 0;
		SetThreadSourceRange(sourceCount, t, localThreadCount, begin, end);
		localUnknown[t].Reset(dofMap.equationBlockCount, dofMap.unknownBlockCount);
		localKnown[t].Reset(dofMap.equationBlockCount, dofMap.knownBlockCount);
		args[t].elements = elements;
		args[t].domain = &domain;
		args[t].dofMap = &dofMap;
		args[t].infElePid = infElePid;
		args[t].elePid = elePid;
		args[t].unknownBuilder = &localUnknown[t];
		args[t].knownBuilder = &localKnown[t];
		args[t].sourceBegin = begin;
		args[t].sourceEnd = end;
		args[t].threadId = t;
		args[t].blockCount = 0;
		args[t].ok = 1;
		int code = pthread_create(&threads[t], 0, MultiDomainStep0ThreadMain, &args[t]);
		if (code == 0)
		{
			created[t] = 1;
			++createSuccess;
		}
		else
		{
			printf("MultiDomain Step0 pthread_create failed domain=%d thread=%ld code=%d\n", domain.id, t, code);
			createOk = 0;
			break;
		}
	}

	for (long t = 0; t < localThreadCount; ++t)
	{
		if (created[t])
			pthread_join(threads[t], 0);
	}

	long nonInterfaceBlocks = 0;
	int workerOk = createOk && (createSuccess == localThreadCount);
	if (workerOk)
	{
		for (long t = 0; t < localThreadCount; ++t)
		{
			if (!args[t].ok)
				workerOk = 0;
			nonInterfaceBlocks += args[t].blockCount;
			unknownBuilder.MergeFrom(localUnknown[t]);
			knownBuilder.MergeFrom(localKnown[t]);
		}
	}

	delete[] localKnown;
	delete[] localUnknown;
	delete[] args;
	delete[] created;
	delete[] threads;

		if (workerOk)
			interfaceBlocks = AssembleMultiDomainStep0InterfaceSerialForDomain(elements, model, dofMap,
				domain, infElePid, elePid, unknownBuilder);

		printf("MultiDomain Step0 pthread assembly domain=%d threads=%ld created=%ld nonInterfaceBlocks=%ld interfaceSerialBlocks=%ld time=%lf\n",
		domain.id, localThreadCount, createSuccess, nonInterfaceBlocks, interfaceBlocks,
		(double)(clock() - domainClock) / (double)CLOCKS_PER_SEC);
	return workerOk;
}

static int AssembleMultiDomainHistoryDomainPthread(DSquareElement* elements,
	const GlobalDofMap& dofMap,
	const Domain& domain,
	long step,
	long requestedThreads,
	MultiDomainCCSRBuilder& mtsBuilder,
	MultiDomainCCSRBuilder& mgsBuilder)
{
	long sourceCount = (long)domain.elementIds.size();
	if (sourceCount <= 0)
		return 1;

	long localThreadCount = EffectiveMultiDomainThreadCount(domain, requestedThreads);
	clock_t domainClock = clock();

	if (localThreadCount <= 1)
	{
		MultiDomainHistoryThreadArg arg;
		arg.elements = elements;
		arg.domain = &domain;
		arg.dofMap = &dofMap;
		arg.step = step;
		arg.mtsBuilder = &mtsBuilder;
		arg.mgsBuilder = &mgsBuilder;
		arg.sourceBegin = 0;
		arg.sourceEnd = sourceCount;
		arg.threadId = 0;
			arg.blockCount = 0;
			arg.ok = 1;
			MultiDomainHistoryThreadMain(&arg);
			return arg.ok;
	}

	pthread_t* threads = new pthread_t[(size_t)localThreadCount];
	int* created = new int[(size_t)localThreadCount];
	MultiDomainHistoryThreadArg* args = new MultiDomainHistoryThreadArg[(size_t)localThreadCount];
	MultiDomainCCSRBuilder* localMts = new MultiDomainCCSRBuilder[(size_t)localThreadCount];
	MultiDomainCCSRBuilder* localMgs = new MultiDomainCCSRBuilder[(size_t)localThreadCount];
	RunPthreadWorkers(threads, created, localThreadCount);

	int createOk = 1;
	long createSuccess = 0;
	for (long t = 0; t < localThreadCount; ++t)
	{
		long begin = 0;
		long end = 0;
		SetThreadSourceRange(sourceCount, t, localThreadCount, begin, end);
		localMts[t].Reset(dofMap.equationBlockCount, dofMap.historyUBlockCount);
		localMgs[t].Reset(dofMap.equationBlockCount, dofMap.historyTBlockCount);
		args[t].elements = elements;
		args[t].domain = &domain;
		args[t].dofMap = &dofMap;
		args[t].step = step;
		args[t].mtsBuilder = &localMts[t];
		args[t].mgsBuilder = &localMgs[t];
		args[t].sourceBegin = begin;
		args[t].sourceEnd = end;
		args[t].threadId = t;
		args[t].blockCount = 0;
		args[t].ok = 1;
		int code = pthread_create(&threads[t], 0, MultiDomainHistoryThreadMain, &args[t]);
		if (code == 0)
		{
			created[t] = 1;
			++createSuccess;
		}
		else
		{
			printf("MultiDomain History pthread_create failed step=%ld domain=%d thread=%ld code=%d\n",
				step, domain.id, t, code);
			createOk = 0;
			break;
		}
	}

	for (long t = 0; t < localThreadCount; ++t)
	{
		if (created[t])
			pthread_join(threads[t], 0);
	}

	long blockCount = 0;
	int workerOk = createOk && (createSuccess == localThreadCount);
	if (workerOk)
	{
		for (long t = 0; t < localThreadCount; ++t)
		{
			if (!args[t].ok)
				workerOk = 0;
			blockCount += args[t].blockCount;
			mtsBuilder.MergeFrom(localMts[t]);
			mgsBuilder.MergeFrom(localMgs[t]);
		}
	}

	delete[] localMgs;
	delete[] localMts;
	delete[] args;
	delete[] created;
	delete[] threads;

		return workerOk;
	}

static int AssembleMultiDomainStep0Pthread(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long** infElePid,
	long** elePid,
	MultiDomainCCSRBuilder& unknownBuilder,
	MultiDomainCCSRBuilder& knownBuilder,
	long requestedThreads)
{
	if (!elements || !model.Validate() || !dofMap.Validate())
		return 0;
	if (requestedThreads <= 1 || MultiDomainPthreadAssemblyDisabled())
		return AssembleMultiDomainStep0(elements, model, dofMap, infElePid, elePid, unknownBuilder, knownBuilder);

	unknownBuilder.Reset(dofMap.equationBlockCount, dofMap.unknownBlockCount);
	knownBuilder.Reset(dofMap.equationBlockCount, dofMap.knownBlockCount);
	printf("MultiDomain GMRES pthread assembly enabled: requestedThreads=%ld domains=%d\n",
		requestedThreads, model.DomainCount());

	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		printf("MultiDomain GMRES material domain=%d C1=%lf C2=%lf C1Dt=%lf C2Dt=%lf elements=%ld\n",
			domain.id, domain.mat.C1, domain.mat.C2, domain.mat.C1Dt, domain.mat.C2Dt,
			(long)domain.elementIds.size());
		long interfaceBlocks = 0;
		if (!AssembleMultiDomainStep0DomainPthread(elements, model, dofMap, domain, requestedThreads,
			infElePid, elePid, unknownBuilder, knownBuilder, interfaceBlocks))
		{
			printf("MultiDomain Step0 pthread assembly failed; falling back to serial multi-domain assembly.\n");
			return AssembleMultiDomainStep0(elements, model, dofMap, infElePid, elePid, unknownBuilder, knownBuilder);
		}
	}

	return CompareOrFallbackStep0(elements, model, dofMap, infElePid, elePid, unknownBuilder, knownBuilder);
}

static int AssembleMultiDomainHistoryStepPthread(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long step,
	MultiDomainCCSRBuilder& mtsBuilder,
	MultiDomainCCSRBuilder& mgsBuilder,
	long requestedThreads)
{
	if (!elements || step <= 0 || !model.Validate() || !dofMap.Validate())
		return 0;
	if (requestedThreads <= 1 || MultiDomainPthreadAssemblyDisabled())
		return AssembleMultiDomainHistoryStep(elements, model, dofMap, step, mtsBuilder, mgsBuilder);

	mtsBuilder.Reset(dofMap.equationBlockCount, dofMap.historyUBlockCount);
	mgsBuilder.Reset(dofMap.equationBlockCount, dofMap.historyTBlockCount);
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		long localStep = MultiRateLocalHistoryStep(domain, step);
		if (!AssembleMultiDomainHistoryDomainPthread(elements, dofMap, domain, localStep, requestedThreads,
			mtsBuilder, mgsBuilder))
		{
			printf("MultiDomain History pthread assembly failed at step=%ld; falling back to serial multi-domain assembly.\n", step);
			return AssembleMultiDomainHistoryStep(elements, model, dofMap, step, mtsBuilder, mgsBuilder);
		}
	}

	return CompareOrFallbackHistory(elements, model, dofMap, step, mtsBuilder, mgsBuilder);
}

static void FillDeterministicVector(Wvector& x)
{
	for (long i = 0; i < x.n; ++i)
	{
		long v = (i * 37 + 11) % 101;
		x.b[i] = ((double)v - 50.0) / 37.0;
	}
}

static void AddVectorInPlace(Wvector& dst, const Wvector& src, double scale)
{
	for (long i = 0; i < dst.n && i < src.n; ++i)
		dst.b[i] += scale * src.b[i];
}

static int CompareMatVecResultWithTolerance(const char* name,
	const Wvector& reference,
	const Wvector& candidate,
	double relativeTolerance,
	double absoluteTolerance)
{
	double maxAbsDiff = 0.0;
	double maxRef = 0.0;
	for (long i = 0; i < reference.n && i < candidate.n; ++i)
	{
		double diff = fabs(reference.b[i] - candidate.b[i]);
		if (diff > maxAbsDiff)
			maxAbsDiff = diff;
		double refAbs = fabs(reference.b[i]);
		if (refAbs > maxRef)
			maxRef = refAbs;
	}
	double relativeDiff = maxAbsDiff / (maxRef > 1.0 ? maxRef : 1.0);
	printf("MultiDomain matrix self-check %s: maxAbsDiff=%e relativeDiff=%e\n",
		name, maxAbsDiff, relativeDiff);
	return (relativeDiff <= relativeTolerance || maxAbsDiff <= absoluteTolerance) ? 1 : 0;
}

static int CompareMatVecResult(const char* name, const Wvector& reference, const Wvector& candidate)
{
	return CompareMatVecResultWithTolerance(name, reference, candidate, 1.0e-8, 1.0e-10);
}

static void PrepareOldCCSRTemp(WmatrixCCSR& temp, long blockCount)
{
	temp.initial(3 * blockCount, 3 * blockCount);
	temp.reinitial();
}

static void PrepareOldSymCCSRTemp(WmatrixSymCCSR& temp, long blockCount)
{
	temp.initial(3 * blockCount, 3 * blockCount);
	temp.reinitial();
}

static void ConvertOldCCSRTemp(WmatrixCCSR& temp, CCSRMat& out)
{
	temp.sum();
	out.initial(temp.m, temp.num);
	out.row_ptr[0] = 0;
	SparseConvert(temp, out);
}

static void ConvertOldSymCCSRTemp(WmatrixSymCCSR& temp, SymCCSRMat& out)
{
	temp.sum();
	out.initial(temp.m, temp.num);
	out.row_ptr[0] = 0;
	SparseConvert(temp, out);
}

static void BuildOldPthEquivalentUnknown0(DSquareElement* elements,
	WmatrixCCSR& out,
	long blockCount,
	long eleCount,
	long** infElePid,
	long** elePid)
{
	long threadId = 0;
	for (long sourceNode = 0; sourceNode < blockCount; ++sourceNode)
	{
		long pos = 3 * sourceNode;
		for (long fieldEle = 0; fieldEle < eleCount; ++fieldEle)
		{
			int flag = UnknownSubMatrix(elements, sourceNode, fieldEle, infElePid, elePid, threadId);
			if (!flag)
				continue;
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long colNode = elements[fieldEle].m_nodeID[localNode];
				out.assign(AssistCoef[threadId][localNode], pos, 3 * colNode, 3, 3);
			}
		}
	}
}

static void BuildOldPthEquivalentKnown0(DSquareElement* elements,
	WmatrixCCSR& out,
	long blockCount,
	long eleCount,
	long** infElePid,
	long** elePid)
{
	long threadId = 0;
	for (long sourceNode = 0; sourceNode < blockCount; ++sourceNode)
	{
		long pos = 3 * sourceNode;
		for (long fieldEle = 0; fieldEle < eleCount; ++fieldEle)
		{
			int flag = knownSubMatrix(elements, sourceNode, fieldEle, infElePid, elePid, threadId);
			if (!flag)
				continue;
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long colNode = elements[fieldEle].m_nodeID[localNode];
				out.assign(AssistCoef[threadId][localNode], pos, 3 * colNode, 3, 3);
			}
		}
	}
}

static void BuildOldPthEquivalentT1(DSquareElement* elements,
	WmatrixCCSR& out,
	long blockCount,
	long eleCount,
	double step,
	double dt)
{
	long threadId = 0;
	for (long sourceNode = 0; sourceNode < blockCount; ++sourceNode)
	{
		long pos = 3 * sourceNode;
		for (long fieldEle = 0; fieldEle < eleCount; ++fieldEle)
		{
			int flag = elements[fieldEle].IntDynaTijJudge(1, elements[fieldEle].m_nodelist[sourceNode], step, dt, threadId);
			if (!flag)
				continue;
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long colNode = elements[fieldEle].m_nodeID[localNode];
				out.assign(AssistTij[threadId][localNode], pos, 3 * colNode, 3, 3);
			}
		}
	}
}

static void BuildOldPthEquivalentT2(DSquareElement* elements,
	WmatrixCCSR& out,
	long blockCount,
	long eleCount,
	double step,
	double dt)
{
	long threadId = 0;
	for (long sourceNode = 0; sourceNode < blockCount; ++sourceNode)
	{
		long pos = 3 * sourceNode;
		for (long fieldEle = 0; fieldEle < eleCount; ++fieldEle)
		{
			int flag = 0;
			int ipos = -1;
			if (step <= 1.01 && elements[fieldEle].IsIn(sourceNode, ipos))
			{
				for (int localNode = 0; localNode < 1; ++localNode)
					for (int k = 0; k < 9; ++k)
						AssistTij[threadId][localNode][k] = elements[fieldEle].m_OnElementT2ij[ipos][localNode][k];
				flag = 1;
			}
			else
			{
				flag = elements[fieldEle].IntDynaTijJudge(2, elements[fieldEle].m_nodelist[sourceNode], step, dt, threadId);
			}

			if (!flag)
				continue;
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long colNode = elements[fieldEle].m_nodeID[localNode];
				out.assign(AssistTij[threadId][localNode], pos, 3 * colNode, 3, 3);
			}
		}
	}
}

static void BuildOldPthEquivalentU(DSquareElement* elements,
	WmatrixSymCCSR& out,
	long blockCount,
	long eleCount,
	double step,
	double dt)
{
	long threadId = 0;
	for (long sourceNode = 0; sourceNode < blockCount; ++sourceNode)
	{
		long pos = 3 * sourceNode;
		for (long fieldEle = 0; fieldEle < eleCount; ++fieldEle)
		{
			int flag = elements[fieldEle].IntDynaUijJudge(elements[fieldEle].m_nodelist[sourceNode], step, dt, threadId);
			if (!flag)
				continue;
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long colNode = elements[fieldEle].m_nodeID[localNode];
				out.assign(AssistUij[threadId][localNode], pos, 3 * colNode, 3, 3);
			}
		}
	}
}
int RunMultiDomainMatrixSelfCheck(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long** infElePid,
	long** elePid)
{
	if (!elements || model.DomainCount() != 1)
		return 1;

	long blockCount = model.nodeCount;
	if (model.domains.empty())
	{
		printf("MultiDomain matrix self-check failed: no domain materials.\n");
		return 0;
	}
	double dt = model.domains[0].mat.Dt;
	int ok = 1;

	MultiDomainCCSRBuilder unknownBuilder(dofMap.equationBlockCount, dofMap.unknownBlockCount);
	MultiDomainCCSRBuilder knownBuilder(dofMap.equationBlockCount, dofMap.knownBlockCount);
	if (!AssembleMultiDomainStep0(elements, model, dofMap, infElePid, elePid, unknownBuilder, knownBuilder))
		return 0;
	CCSRMat newMTS0;
	CCSRMat newKnown;
	unknownBuilder.Build(newMTS0);
	knownBuilder.Build(newKnown);

	WmatrixCCSR oldTemp;
	oldTemp.a = 0;
	oldTemp.col = 0;
	oldTemp.row = 0;
	CCSRMat oldMTS0;
	PrepareOldCCSRTemp(oldTemp, blockCount);
	BuildOldPthEquivalentUnknown0(elements, oldTemp, blockCount, model.elementCount, infElePid, elePid);
	ConvertOldCCSRTemp(oldTemp, oldMTS0);

	Wvector xUnknown(3 * dofMap.unknownBlockCount, 0);
	Wvector oldY(3 * dofMap.equationBlockCount, 0);
	Wvector newY(3 * dofMap.equationBlockCount, 0);
	FillDeterministicVector(xUnknown);
	SparseMul(oldMTS0, xUnknown, oldY);
	SparseMul(newMTS0, xUnknown, newY);
	ok = CompareMatVecResult("MTS0", oldY, newY) && ok;

	CCSRMat oldKnown;
	PrepareOldCCSRTemp(oldTemp, blockCount);
	BuildOldPthEquivalentKnown0(elements, oldTemp, blockCount, model.elementCount, infElePid, elePid);
	ConvertOldCCSRTemp(oldTemp, oldKnown);

	Wvector xKnown(3 * dofMap.knownBlockCount, 0);
	FillDeterministicVector(xKnown);
	SparseMul(oldKnown, xKnown, oldY);
	SparseMul(newKnown, xKnown, newY);
	ok = CompareMatVecResult("MTS[0]", oldY, newY) && ok;

	MultiDomainCCSRBuilder mtsBuilder(dofMap.equationBlockCount, dofMap.historyUBlockCount);
	MultiDomainCCSRBuilder mgsBuilder(dofMap.equationBlockCount, dofMap.historyTBlockCount);
	if (!AssembleMultiDomainHistoryStep(elements, model, dofMap, 1, mtsBuilder, mgsBuilder))
		return 0;
	CCSRMat newMTS1;
	CCSRMat newMGS1Full;
	SymCCSRMat newMGS1;
	bool newMgsSymmetric = true;
	mtsBuilder.Build(newMTS1);
	mgsBuilder.Build(newMGS1Full);
	mgsBuilder.BuildSym(newMGS1, 1.0e-10, &newMgsSymmetric);
	if (!newMgsSymmetric)
		printf("MultiDomain matrix self-check warning: new MGS[1] has non-symmetric 3x3 blocks; full CCSR representation remains available.\n");

	CCSRMat oldT2;
	PrepareOldCCSRTemp(oldTemp, blockCount);
	BuildOldPthEquivalentT2(elements, oldTemp, blockCount, model.elementCount, 1.0, dt);
	ConvertOldCCSRTemp(oldTemp, oldT2);

	CCSRMat oldT1;
	PrepareOldCCSRTemp(oldTemp, blockCount);
	BuildOldPthEquivalentT1(elements, oldTemp, blockCount, model.elementCount, 1.0, dt);
	ConvertOldCCSRTemp(oldTemp, oldT1);

	Wvector xU(3 * dofMap.historyUBlockCount, 0);
	Wvector oldT2Y(3 * dofMap.equationBlockCount, 0);
	Wvector oldT1Y(3 * dofMap.equationBlockCount, 0);
	FillDeterministicVector(xU);
	SparseMul(oldT2, xU, oldT2Y);
	SparseMul(oldT1, xU, oldT1Y);
	oldY = oldT2Y;
	AddVectorInPlace(oldY, oldT1Y, 1.0);
	SparseMul(newMTS1, xU, newY);
	ok = CompareMatVecResult("MTS[1]", oldY, newY) && ok;

	WmatrixSymCCSR oldSymTemp;
	oldSymTemp.a = 0;
	oldSymTemp.col = 0;
	oldSymTemp.row = 0;
	SymCCSRMat oldMGS1;
	PrepareOldSymCCSRTemp(oldSymTemp, blockCount);
	BuildOldPthEquivalentU(elements, oldSymTemp, blockCount, model.elementCount, 1.0, dt);
	ConvertOldSymCCSRTemp(oldSymTemp, oldMGS1);

	Wvector xT(3 * dofMap.historyTBlockCount, 0);
	FillDeterministicVector(xT);
	SparseMul(oldMGS1, xT, oldY);
	if (newMgsSymmetric)
		SparseMul(newMGS1, xT, newY);
	else
		SparseMul(newMGS1Full, xT, newY);
	ok = CompareMatVecResult("MGS[1]", oldY, newY) && ok;

	if (!ok)
		printf("MultiDomain matrix self-check failed. Stop before multi-domain GMRES solve.\n");
	else
		printf("MultiDomain matrix self-check passed for DomainCount=1.\n");
	return ok;
}

static void BuildKnownInputVectorFromBoundary(const GlobalDofMap& dofMap,
	const MultiDomainModel& model,
	DSquareElement* elements,
	BoundaryValue& stepBd,
	Wvector& known)
{
	long expected = 3 * dofMap.knownBlockCount;
	if (!known.b)
		known.initial(expected);
	else if (known.n != expected)
	{
		known.finish();
		known.initial(expected);
	}
	known = 0.0;
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		for (size_t i = 0; i < domain.elementIds.size(); ++i)
		{
			long ele = domain.elementIds[i];
			if (elements[ele].SurfaceType == SurfaceInterface)
				continue;
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long nodeId = elements[ele].m_nodeID[localNode];
				VariableRef ref = elements[ele].BCID == 123 ? dofMap.GetU(domain.id, nodeId) : dofMap.GetT(domain.id, nodeId);
				if (!ref.known)
					continue;
				long src = nodeId * 3;
				long dst = ref.globalBlock * 3;
				known[dst + 1] = ref.sign * stepBd.lu.b[src + 0];
				known[dst + 2] = ref.sign * stepBd.lu.b[src + 1];
				known[dst + 3] = ref.sign * stepBd.lu.b[src + 2];
			}
		}
	}
}

static void ScatterUnknownToBoundary(const GlobalDofMap& dofMap,
	const MultiDomainModel& model,
	DSquareElement* elements,
	const Wvector& unknown,
	BoundaryValue& stepBd)
{
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		for (size_t i = 0; i < domain.elementIds.size(); ++i)
		{
			long ele = domain.elementIds[i];
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				double transform[9];
				GetElementLocalTransformFromInterfaceReference(model, domain.id, ele, localNode, transform);
				long nodeId = elements[ele].m_nodeID[localNode];
				VariableRef uRef = dofMap.GetU(domain.id, nodeId);
				VariableRef tRef = dofMap.GetT(domain.id, nodeId);
				long base = nodeId * 3;
				double inputU[3] = { stepBd.lu.b[base + 0], stepBd.lu.b[base + 1], stepBd.lu.b[base + 2] };
				double inputT[3] = { stepBd.lu.b[base + 0], stepBd.lu.b[base + 1], stepBd.lu.b[base + 2] };

				if (!uRef.known)
				{
					long src = uRef.globalBlock * 3;
					double refValue[3] = { unknown.b[src + 0], unknown.b[src + 1], unknown.b[src + 2] };
					double localValue[3];
					ApplyTransform3(transform, refValue, localValue);
					stepBd.lu.b[base + 0] = uRef.sign * localValue[0];
					stepBd.lu.b[base + 1] = uRef.sign * localValue[1];
					stepBd.lu.b[base + 2] = uRef.sign * localValue[2];
				}
				if (!tRef.known)
				{
					long src = tRef.globalBlock * 3;
					double refValue[3] = { unknown.b[src + 0], unknown.b[src + 1], unknown.b[src + 2] };
					double localValue[3];
					ApplyTransform3(transform, refValue, localValue);
					stepBd.lt.b[base + 0] = tRef.sign * localValue[0];
					stepBd.lt.b[base + 1] = tRef.sign * localValue[1];
					stepBd.lt.b[base + 2] = tRef.sign * localValue[2];
				}
				if (uRef.known)
				{
					stepBd.lu.b[base + 0] = uRef.sign * inputU[0];
					stepBd.lu.b[base + 1] = uRef.sign * inputU[1];
					stepBd.lu.b[base + 2] = uRef.sign * inputU[2];
				}
				if (tRef.known)
				{
					stepBd.lt.b[base + 0] = tRef.sign * inputT[0];
					stepBd.lt.b[base + 1] = tRef.sign * inputT[1];
					stepBd.lt.b[base + 2] = tRef.sign * inputT[2];
				}
			}
		}
	}
}
static void InitializeMultiDomainState(MultiDomainState& state,
	const GlobalDofMap& dofMap)
{
	if (!state.globalUnknown.b)
		state.globalUnknown.initial(3 * dofMap.unknownBlockCount);
	if (!state.globalKnown.b)
		state.globalKnown.initial(3 * dofMap.knownBlockCount);
	if (!state.globalU.b)
		state.globalU.initial(3 * dofMap.historyUBlockCount);
	if (!state.globalT.b)
		state.globalT.initial(3 * dofMap.historyTBlockCount);
	state.globalUnknown = 0.0;
	state.globalKnown = 0.0;
	state.globalU = 0.0;
	state.globalT = 0.0;
}

static void BuildPhysicalStateFromBoundary(const GlobalDofMap& dofMap,
	const MultiDomainModel& model,
	DSquareElement* elements,
	BoundaryValue& stepBd,
	MultiDomainState& state)
{
	state.globalU = 0.0;
	state.globalT = 0.0;
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		for (size_t i = 0; i < domain.elementIds.size(); ++i)
		{
			long ele = domain.elementIds[i];
			if (elements[ele].SurfaceType == SurfaceInterface &&
				!IsInterfaceReferenceSide(model, domain.id, ele))
			{
				continue;
			}
			for (int localNode = 0; localNode < 1; ++localNode)
			{
				long nodeId = elements[ele].m_nodeID[localNode];
				long src = nodeId * 3;

				VariableRef uRef = dofMap.GetHistoryU(domain.id, nodeId);
				long uDst = uRef.globalBlock * 3;
				state.globalU.b[uDst + 0] = uRef.sign * stepBd.lu.b[src + 0];
				state.globalU.b[uDst + 1] = uRef.sign * stepBd.lu.b[src + 1];
				state.globalU.b[uDst + 2] = uRef.sign * stepBd.lu.b[src + 2];

				VariableRef tRef = dofMap.GetHistoryT(domain.id, nodeId);
				long tDst = tRef.globalBlock * 3;
				state.globalT.b[tDst + 0] = tRef.sign * stepBd.lt.b[src + 0];
				state.globalT.b[tDst + 1] = tRef.sign * stepBd.lt.b[src + 1];
				state.globalT.b[tDst + 2] = tRef.sign * stepBd.lt.b[src + 2];
			}
		}
	}
}
static void BuildConvolutionVector(const std::vector<MultiDomainState>& states,
	long step,
	long convol,
	int useTraction,
	Wvector& out)
{
	out = 0.0;
	long ids[3];
	double coef[3];
	ids[0] = step - convol + 1;
	ids[1] = step - convol;
	ids[2] = step - convol - 1;
	coef[0] = 1.0;
	coef[1] = 2.0;
	coef[2] = 1.0;

	for (int k = 0; k < 3; ++k)
	{
		if (ids[k] > 0 && ids[k] < step)
		{
			const Wvector& src = useTraction ? states[(size_t)ids[k]].globalT : states[(size_t)ids[k]].globalU;
			AddVectorInPlace(out, src, coef[k]);
		}
	}
}

static void PrintInterfaceError(const GlobalDofMap& dofMap,
	const MultiDomainModel& model,
	DSquareElement* elements,
	BoundaryValue& stepBd,
	long step)
{
	double maxUDiff = 0.0;
	double maxTBalance = 0.0;
	for (size_t i = 0; i < model.interfaces.size(); ++i)
	{
		const InterfacePair& itf = model.interfaces[i];
		int tSignB = itf.normalSign == 0 ? -1 : itf.normalSign;
		for (int localNode = 0; localNode < 1; ++localNode)
		{
			long a = InterfaceNodeA(itf, localNode) * 3;
			long b = InterfaceNodeBForA(itf, localNode) * 3;
			for (int k = 0; k < 3; ++k)
			{
				double ud = fabs(stepBd.lu.b[a + k] - stepBd.lu.b[b + k]);
				if (ud > maxUDiff)
					maxUDiff = ud;
				double tb = fabs(stepBd.lt.b[b + k] - (double)tSignB * stepBd.lt.b[a + k]);
				if (tb > maxTBalance)
					maxTBalance = tb;
			}
		}
	}
	printf("MultiDomain interface error step %ld: max|uA-uB|=%e max|tB-normalSign*tA|=%e\n",
		step, maxUDiff, maxTBalance);
	(void)dofMap;
	(void)elements;
}
static void SetIdentityPreBlock(Wmatrix& block)
{
	block = 0.0;
	block(1, 1) = 1.0;
	block(2, 2) = 1.0;
	block(3, 3) = 1.0;
}


static void CopyCCSRBlockToMatrix(const CCSRMat& matrix, int blockIndex, Wmatrix& block)
{
	for (int r = 0; r < 3; ++r)
	{
		for (int c = 0; c < 3; ++c)
			block(r + 1, c + 1) = matrix.value[blockIndex * 9 + r * 3 + c];
	}
}

static int IsFiniteMatrix3(const Wmatrix& block)
{
	for (int r = 1; r <= 3; ++r)
	{
		for (int c = 1; c <= 3; ++c)
		{
			if (!IsFiniteScalar(block.a[r - 1][c - 1]))
				return 0;
		}
	}
	return 1;
}

static int TryInvertCCSRBlock(const CCSRMat& matrix, int blockIndex, Wmatrix& inverse)
{
	Wmatrix candidate(3, 3);
	CopyCCSRBlockToMatrix(matrix, blockIndex, candidate);
	if (!IsFiniteMatrix3(candidate))
		return 0;
	if (inv_mat(candidate, 3) < 0)
		return 0;
	if (!IsFiniteMatrix3(candidate))
		return 0;
	inverse = candidate;
	return 1;
}

static int FindCCSRBlockInRow(const CCSRMat& matrix, long row, long col)
{
	if (row < 0 || row >= matrix.n || col < 0 || col >= matrix.n || !matrix.row_ptr || !matrix.col_index)
		return -1;
	for (int p = matrix.row_ptr[row]; p < matrix.row_ptr[row + 1]; ++p)
	{
		if (matrix.col_index[p] == col)
			return p;
	}
	return -1;
}

static int DiagnoseMultiDomainMatrixStructure(const CCSRMat& matrix, const GlobalDofMap& dofMap)
{
	g_validationMatrixStructureValid = 0;
	g_validationMatrixZeroRows = 0;
	g_validationMatrixZeroCols = 0;
	g_validationMatrixNonFiniteBlocks = 0;
	g_validationMatrixDiagonalBlocks = 0;
	g_validationMatrixInvertibleDiagonalBlocks = 0;

	if (matrix.n <= 0 || !matrix.row_ptr || !matrix.col_index || !matrix.value)
	{
		printf("MultiDomain matrix structure invalid: empty or uninitialized matrix.\n");
		return 0;
	}

	std::vector<long> colCounts((size_t)matrix.n, 0);
	std::vector<long> zeroColumnExamples;
	long printedZeroRows = 0;
	long printedBadBlocks = 0;
	long printedBadCols = 0;
	for (long row = 0; row < matrix.n; ++row)
	{
		int begin = matrix.row_ptr[row];
		int end = matrix.row_ptr[row + 1];
		if (begin == end)
		{
			++g_validationMatrixZeroRows;
			if (printedZeroRows < 20)
			{
				printf("MultiDomain matrix zero row block %ld.\n", row);
				++printedZeroRows;
			}
		}
		for (int pos = begin; pos < end; ++pos)
		{
			int col = matrix.col_index[pos];
			if (col >= 0 && col < matrix.n)
				++colCounts[(size_t)col];
			else
			{
				++g_validationMatrixNonFiniteBlocks;
				if (printedBadBlocks < 20)
				{
					printf("MultiDomain matrix invalid column row=%ld blockIndex=%d col=%d.\n", row, pos, col);
					++printedBadBlocks;
				}
				continue;
			}

			int finiteBlock = 1;
			for (int k = 0; k < 9; ++k)
			{
				if (!IsFiniteScalar(matrix.value[pos * 9 + k]))
				{
					finiteBlock = 0;
					break;
				}
			}
			if (!finiteBlock)
			{
				++g_validationMatrixNonFiniteBlocks;
				if (printedBadBlocks < 20)
				{
					printf("MultiDomain matrix non-finite block row=%ld col=%d blockIndex=%d.\n", row, col, pos);
					++printedBadBlocks;
				}
			}
			if (col == row)
			{
				++g_validationMatrixDiagonalBlocks;
				Wmatrix inverse(3, 3);
				if (TryInvertCCSRBlock(matrix, pos, inverse))
					++g_validationMatrixInvertibleDiagonalBlocks;
			}
		}
	}

	for (long col = 0; col < matrix.n; ++col)
	{
		if (colCounts[(size_t)col] == 0)
		{
			++g_validationMatrixZeroCols;
			if ((long)zeroColumnExamples.size() < 50)
				zeroColumnExamples.push_back(col);
			if (printedBadCols < 20)
			{
				printf("MultiDomain matrix zero column block %ld: %s\n", col, dofMap.GetUnknownDebugLabel(col));
				++printedBadCols;
			}
		}
	}

	printf("MultiDomain matrix structure: rows=%ld zeroRows=%ld zeroCols=%ld nonFiniteBlocks=%ld diagonalBlocks=%ld invertibleDiagonalBlocks=%ld.\n",
		matrix.n,
		g_validationMatrixZeroRows,
		g_validationMatrixZeroCols,
		g_validationMatrixNonFiniteBlocks,
		g_validationMatrixDiagonalBlocks,
		g_validationMatrixInvertibleDiagonalBlocks);

	g_validationMatrixStructureValid =
		(g_validationMatrixZeroRows == 0 &&
			g_validationMatrixZeroCols == 0 &&
			g_validationMatrixNonFiniteBlocks == 0) ? 1 : 0;

	if (MultiDomainValidationOutputEnabled())
	{
		FILE* metrics = 0;
		std::string validationMetricsPath = DBEMOutputPath("validation_metrics.txt");
		fopen_s(&metrics, validationMetricsPath.c_str(), "a");
		if (metrics)
		{
			fprintf(metrics, "MultiDomainMatrixStructureValid=%d\n", g_validationMatrixStructureValid);
			fprintf(metrics, "MultiDomainMatrixZeroRows=%ld\n", g_validationMatrixZeroRows);
			fprintf(metrics, "MultiDomainMatrixZeroCols=%ld\n", g_validationMatrixZeroCols);
			fprintf(metrics, "MultiDomainMatrixNonFiniteBlocks=%ld\n", g_validationMatrixNonFiniteBlocks);
			fprintf(metrics, "MultiDomainMatrixDiagonalBlocks=%ld\n", g_validationMatrixDiagonalBlocks);
			fprintf(metrics, "MultiDomainMatrixInvertibleDiagonalBlocks=%ld\n", g_validationMatrixInvertibleDiagonalBlocks);
			for (size_t i = 0; i < zeroColumnExamples.size(); ++i)
			{
				long col = zeroColumnExamples[i];
				fprintf(metrics, "MultiDomainMatrixZeroColumnLabel[%ld]=block %ld %s\n",
					(long)i, col, dofMap.GetUnknownDebugLabel(col));
			}
			fclose(metrics);
		}
	}
	return g_validationMatrixStructureValid;
}

static int TryUsePreconditionerColumn(const CCSRMat& matrix,
	long row,
	long col,
	Wmatrix& inverse,
	long* existingBlockCount,
	long* invertedBlockCount)
{
	int blockIndex = FindCCSRBlockInRow(matrix, row, col);
	if (blockIndex < 0)
		return 0;
	if (existingBlockCount)
		++(*existingBlockCount);
	if (!TryInvertCCSRBlock(matrix, blockIndex, inverse))
		return 0;
	if (invertedBlockCount)
		++(*invertedBlockCount);
	return 1;
}

static void CaptureMultiRateValidationMetrics(const MultiDomainModel& model)
{
	g_validationMultiRateTimeMode = model.multiRateTimeMode;
	g_validationMacroDt = model.macroDt;
	g_validationBetaLimit = model.betaLimit;
	g_validationAveElementSize = model.aveElementSize;
	g_validationMaxSubstepsPerMacro = model.maxSubstepsPerMacro;
	g_validationDomainSubsteps.clear();
	g_validationDomainStableDt.clear();
	g_validationDomainLocalDt.clear();
	g_validationDomainBetaActual.clear();
	g_validationDomainC1.clear();
	g_validationDomainC2.clear();
	for (size_t i = 0; i < model.domains.size(); ++i)
	{
		const Domain& d = model.domains[i];
		g_validationDomainSubsteps.push_back(d.substeps);
		g_validationDomainStableDt.push_back(d.stableDt);
		g_validationDomainLocalDt.push_back(d.localDt);
		g_validationDomainBetaActual.push_back(d.betaActual);
		g_validationDomainC1.push_back(d.c1);
		g_validationDomainC2.push_back(d.c2);
	}
}


static void AssignCCSRBlockToLocalPreconditioner(const CCSRMat& matrix,
	long row,
	long col,
	Wmatrix& local,
	long localRow,
	long localCol)
{
	Wmatrix subPre(3, 3);
	int blockIndex = FindCCSRBlockInRow(matrix, row, col);
	if (blockIndex >= 0)
	{
		subPre(1, 1) = matrix.value[blockIndex * 9 + 0];
		subPre(1, 2) = matrix.value[blockIndex * 9 + 1];
		subPre(1, 3) = matrix.value[blockIndex * 9 + 2];
		subPre(2, 1) = matrix.value[blockIndex * 9 + 3];
		subPre(2, 2) = matrix.value[blockIndex * 9 + 4];
		subPre(2, 3) = matrix.value[blockIndex * 9 + 5];
		subPre(3, 1) = matrix.value[blockIndex * 9 + 6];
		subPre(3, 2) = matrix.value[blockIndex * 9 + 7];
		subPre(3, 3) = matrix.value[blockIndex * 9 + 8];
	}
	else
	{
		subPre = 0.0;
	}
	local.assign(subPre, 3 * localRow + 1, 3 * localCol + 1, 3, 3);
}

static int BuildMappedLeafPreConditioner(DSquareElement* elements,
	PreConditioner& pre,
	long pointCount,
	long maxLeafPointCount,
	const CCSRMat& matrix,
	const GlobalDofMap& dofMap)
{
	pre.Finish();
	if (!elements || pointCount <= 0 || maxLeafPointCount <= 0 || matrix.n != pointCount)
		return 0;

	Cube globalBox;
	Tree* treePre;
	PointerOfLeaf leafPre;
	AssignCubeSize(elements[0].m_nodelist, pointCount, globalBox);
	InitRoot(treePre, globalBox, pointCount);
	CreateTree(treePre, elements[0].m_nodelist, maxLeafPointCount);
	CreateLeafPointer(leafPre, treePre, pointCount);

	pre.m_LeafCount = leafPre.LeafCount;
	pre.m_LeafPointCount = new long[leafPre.LeafCount];
	pre.m_LeafBeginID = new long[leafPre.LeafCount];
	pre.m_RePID = new long[pointCount];
	pre.IID = new long[pointCount];
	pre.m_InputPID = new long[pointCount];
	pre.m_OutputPID = new long[pointCount];
	pre.m_PreM = new Wmatrix[leafPre.LeafCount];

	for (long i = 0; i < leafPre.LeafCount; ++i)
	{
		pre.m_LeafPointCount[i] = leafPre.LeafPointer[i]->m_PointCount;
		pre.m_LeafBeginID[i] = leafPre.LeafPointer[i]->m_BeginID;
	}
	RenumberPointID(leafPre, pre.m_RePID);
	for (long i = 0; i < pointCount; ++i)
	{
		pre.IID[i] = 3 * i;
		pre.m_InputPID[i] = pre.m_RePID[i];
		pre.m_OutputPID[i] = dofMap.GetPreferredUnknownBlockForRow(pre.m_RePID[i]);
	}

	long missingBlocks = 0;
	long localBlocks = 0;
	for (long leaf = 0; leaf < pre.m_LeafCount; ++leaf)
	{
		long localCount = pre.m_LeafPointCount[leaf];
		pre.m_PreM[leaf].initial(3 * localCount, 3 * localCount);
		long begin = pre.m_LeafBeginID[leaf];
		for (long r = 0; r < localCount; ++r)
		{
			long row = pre.m_InputPID[begin + r];
			for (long c = 0; c < localCount; ++c)
			{
				long col = pre.m_OutputPID[begin + c];
				if (FindCCSRBlockInRow(matrix, row, col) < 0)
					++missingBlocks;
				else
					++localBlocks;
				AssignCCSRBlockToLocalPreconditioner(matrix, row, col, pre.m_PreM[leaf], r, c);
			}
		}
	}

	for (long leaf = 0; leaf < pre.m_LeafCount; ++leaf)
		inv_mat(pre.m_PreM[leaf], pre.m_PreM[leaf].m);

	printf("MultiDomain mapped leaf preconditioner stats: leaves=%ld pointCount=%ld maxLeafPointCount=%ld localBlocks=%ld missingLocalBlocks=%ld.\n",
		pre.m_LeafCount, pointCount, maxLeafPointCount, localBlocks, missingBlocks);

	DeleteTree(treePre);
	DeleteLeafPointer(leafPre);
	return 1;
}

static int BuildBlockDiagonalPreConditioner(PreConditioner& pre,
	const CCSRMat& matrix,
	const GlobalDofMap& dofMap)
{
	pre.Finish();
	long blockCount = matrix.n;
	if (blockCount <= 0)
		return 0;
	pre.m_LeafCount = blockCount;
	pre.m_LeafPointCount = new long[blockCount];
	pre.m_LeafBeginID = new long[blockCount];
	pre.m_RePID = new long[blockCount];
	pre.IID = new long[blockCount];
	pre.m_InputPID = new long[blockCount];
	pre.m_OutputPID = new long[blockCount];
	pre.m_PreM = new Wmatrix[blockCount];
	for (long i = 0; i < blockCount; ++i)
	{
		pre.m_LeafPointCount[i] = 1;
		pre.m_LeafBeginID[i] = i;
		pre.m_RePID[i] = i;
		pre.IID[i] = 3 * i;
		pre.m_InputPID[i] = i;
		pre.m_OutputPID[i] = i;
		pre.m_PreM[i].initial(3, 3);
		SetIdentityPreBlock(pre.m_PreM[i]);
	}

	long preferredBlockCount = 0;
	long preferredInverseCount = 0;
	long alternateBlockCount = 0;
	long alternateInverseCount = 0;
	long identityBlockCount = 0;
	for (long row = 0; row < blockCount; ++row)
	{
		long preferred = dofMap.GetPreferredUnknownBlockForRow(row);
		if (TryUsePreconditionerColumn(matrix, row, preferred, pre.m_PreM[row],
			&preferredBlockCount, &preferredInverseCount))
		{
			pre.m_OutputPID[row] = preferred;
			continue;
		}

		VariableRef uRef = dofMap.GetU(0, row);
		if (!uRef.known && uRef.globalBlock != preferred &&
			TryUsePreconditionerColumn(matrix, row, uRef.globalBlock, pre.m_PreM[row],
				&alternateBlockCount, &alternateInverseCount))
		{
			pre.m_OutputPID[row] = uRef.globalBlock;
			continue;
		}

		VariableRef tRef = dofMap.GetT(0, row);
		if (!tRef.known && tRef.globalBlock != preferred &&
			TryUsePreconditionerColumn(matrix, row, tRef.globalBlock, pre.m_PreM[row],
				&alternateBlockCount, &alternateInverseCount))
		{
			pre.m_OutputPID[row] = tRef.globalBlock;
			continue;
		}

		SetIdentityPreBlock(pre.m_PreM[row]);
		++identityBlockCount;
	}
	printf("MultiDomain block preconditioner stats: rows=%ld preferredBlocks=%ld preferredInverted=%ld alternateBlocks=%ld alternateInverted=%ld identityBlocks=%ld.\n",
		blockCount,
		preferredBlockCount,
		preferredInverseCount,
		alternateBlockCount,
		alternateInverseCount,
		identityBlockCount);
	return 1;
}

static void PrintPreconditionerSelfCheck(const CCSRMat& matrix, PreConditioner& pre)
{
	Wvector x(3 * matrix.n, 0);
	Wvector ax(3 * matrix.n, 0);
	Wvector pax(3 * matrix.n, 0);
	FillDeterministicVector(x);
	SparseMul((CCSRMat&)matrix, x, ax);
	GMRES_M_X(pre, ax, pax);
	double maxAbs = 0.0;
	double ref = 0.0;
	for (long i = 0; i < x.n; ++i)
	{
		double diff = fabs(pax.b[i] - x.b[i]);
		if (diff > maxAbs)
			maxAbs = diff;
		if (fabs(x.b[i]) > ref)
			ref = fabs(x.b[i]);
	}
	double rel = ref > 0.0 ? maxAbs / ref : maxAbs;
	printf("MultiDomain preconditioner self-check M^-1*A*x: maxAbsDiff=%e relativeDiff=%e\n", maxAbs, rel);
}


int DynaGMRESSolverMultiDomainCCSR(DSquareElement* elements,
	BoundaryValue* bd,
	const MultiDomainModel& model,
	long iterations,
	double error,
	long maxleafpointnum,
	long NStep,
	double MaxLength,
	long** infElePid,
	long** elePid)
{
	(void)maxleafpointnum;

	if (!elements || !bd)
		return -1;
	if (!model.IsActive())
	{
		printf("MultiDomain solver called with inactive model.\n");
		return -1;
	}
	g_validationMetricsValid = 0;
	g_validationFinalFlag = -999;
	g_validationSolvedSteps = 0;
	g_validationTotalIter0 = 0;
	g_validationTotalIter1 = 0;
	g_validationGlobalUnknownBlocks = 0;
	g_validationMTS0Blocks = 0;
	g_validationAvgMTSBlocks = 0.0;
	g_validationAvgMGSBlocks = 0.0;
	g_validationMatrixAssemblyTime = 0.0;
	g_validationGMRESTime = 0.0;
	g_validationDomainCount = model.DomainCount();
	g_validationMatrixStructureValid = 0;
	g_validationMatrixZeroRows = 0;
	g_validationMatrixZeroCols = 0;
	g_validationMatrixNonFiniteBlocks = 0;
	g_validationMatrixDiagonalBlocks = 0;
	g_validationMatrixInvertibleDiagonalBlocks = 0;
	g_validationInterfaceDofInvariantValid = 0;
	g_validationInterfaceDofInvariantFailures = 0;
	g_validationMultiRateTimeMode = 0;
	g_validationMacroDt = 0.0;
	g_validationBetaLimit = 0.0;
	g_validationAveElementSize = 0.0;
	g_validationMaxSubstepsPerMacro = 0;
	g_validationDomainSubsteps.clear();
	g_validationDomainStableDt.clear();
	g_validationDomainLocalDt.clear();
	g_validationDomainBetaActual.clear();
	g_validationDomainC1.clear();
	g_validationDomainC2.clear();
	CaptureMultiRateValidationMetrics(model);

	GlobalDofMap dofMap;
	if (!dofMap.Build(model, elements, model.elementCount))
		return -1;
	dofMap.PrintSummary(stdout);
	if (!DiagnoseInterfaceDofInvariants(dofMap, model))
		return -1;

	if (model.DomainCount() == 1 &&
		!RunMultiDomainMatrixSelfCheck(elements, model, dofMap, infElePid, elePid))
	{
		return -1;
	}

	long MaxN = NStep;
	double minC2MacroDt = (std::numeric_limits<double>::max)();
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		long substeps = model.domains[d].substeps > 0 ? model.domains[d].substeps : 1;
		double c2MacroDt = model.domains[d].mat.C2Dt * (double)substeps;
		if (c2MacroDt > 0.0 && c2MacroDt < minC2MacroDt)
			minC2MacroDt = c2MacroDt;
	}
	if (minC2MacroDt == (std::numeric_limits<double>::max)())
	{
		printf("MultiDomain MaxN failed: no valid domain macro C2Dt.\n");
		return -1;
	}
	for (long i = 0; i <= NStep; ++i)
	{
		if (minC2MacroDt * i > MaxLength)
		{
			MaxN = i;
			break;
		}
	}
	if (MaxN > NStep)
		MaxN = NStep;
	printf("MultiDomain MaxN = %ld\n", MaxN);

	MultiDomainCCSRBuilder unknownBuilder(dofMap.equationBlockCount, dofMap.unknownBlockCount);
	MultiDomainCCSRBuilder knownBuilder(dofMap.equationBlockCount, dofMap.knownBlockCount);
	clock_t matrixClock = clock();
	if (!AssembleMultiDomainStep0Pthread(elements, model, dofMap, infElePid, elePid, unknownBuilder, knownBuilder, thread_num))
		return -1;

	printf("MultiDomain MTS0 blocks = %ld\n", unknownBuilder.NonZeroBlocks());
	printf("MultiDomain known blocks = %ld\n", knownBuilder.NonZeroBlocks());

	CCSRMat A;
	CCSRMat KnownM;
	unknownBuilder.Build(A);
	knownBuilder.Build(KnownM);
	if (!DiagnoseMultiDomainMatrixStructure(A, dofMap))
	{
		printf("MultiDomain matrix structure check failed. Stop before GMRES solve.\n");
		return -1;
	}

	CCSRMat* MTS = new CCSRMat[MaxN + 1];
	SymCCSRMat* MGS = new SymCCSRMat[MaxN + 1];
	CCSRMat* MGSFull = new CCSRMat[MaxN + 1];
	int* MGSUseSym = new int[MaxN + 1];
	for (long i = 0; i <= MaxN; ++i)
		MGSUseSym[i] = 1;

	long totalHistoryTBlocks = 0;
	long totalHistoryGBlocks = 0;
	for (long step = 1; step <= MaxN; ++step)
	{
		MultiDomainCCSRBuilder mtsBuilder(dofMap.equationBlockCount, dofMap.historyUBlockCount);
		MultiDomainCCSRBuilder mgsBuilder(dofMap.equationBlockCount, dofMap.historyTBlockCount);
		if (!AssembleMultiDomainHistoryStepPthread(elements, model, dofMap, step, mtsBuilder, mgsBuilder, thread_num))
		{
			delete[] MTS;
			delete[] MGS;
			delete[] MGSFull;
			delete[] MGSUseSym;
			return -1;
		}
		mtsBuilder.Build(MTS[step]);
		bool symmetric = true;
		mgsBuilder.BuildSym(MGS[step], 1.0e-10, &symmetric);
		MGSUseSym[step] = symmetric ? 1 : 0;
		if (!symmetric)
		{
			printf("MultiDomain MGS[%ld] has non-symmetric 3x3 blocks; using full CCSR representation for this step.\n", step);
			mgsBuilder.Build(MGSFull[step]);
		}
		totalHistoryTBlocks += mtsBuilder.NonZeroBlocks();
		totalHistoryGBlocks += mgsBuilder.NonZeroBlocks();
		printf("MultiDomain history step %ld: MTS blocks=%ld MGS blocks=%ld\n",
			step, mtsBuilder.NonZeroBlocks(), mgsBuilder.NonZeroBlocks());
	}
	double matrixAssemblyTime = (double)(clock() - matrixClock) / (double)CLOCKS_PER_SEC;

	PreConditioner pre;
	int preOk = 0;
	if (model.DomainCount() == 1)
		preOk = GMRESPreConditioner(elements, pre, model.nodeCount, maxleafpointnum, A);
	else
		preOk = BuildMappedLeafPreConditioner(elements, pre, model.nodeCount, maxleafpointnum, A, dofMap);
	if (!preOk)
	{
		delete[] MTS;
		delete[] MGS;
		delete[] MGSFull;
		delete[] MGSUseSym;
		return -1;
	}
	if (model.DomainCount() == 1)
		PrintPreconditionerSelfCheck(A, pre);

	std::vector<MultiDomainState> states((size_t)NStep + 1);
	for (long step = 0; step <= NStep; ++step)
		InitializeMultiDomainState(states[(size_t)step], dofMap);

	bd[0].lt = 0.0;
	bd[0].Convert(elements, model.elementCount);
	bd[0].AllLocToAbs(model.nodeCount, DSquareElement::m_transmat);
	BuildPhysicalStateFromBoundary(dofMap, model, elements, bd[0], states[0]);
	WriteMultiRateStateCsv(model, elements, bd[0], 0);

	Wvector rhs(3 * dofMap.equationBlockCount, 0);
	Wvector x0(3 * dofMap.unknownBlockCount, 0);
	Wvector x(3 * dofMap.unknownBlockCount, 0);
	Wvector temp(3 * dofMap.equationBlockCount, 0);
	Wvector diagKnown(3 * dofMap.equationBlockCount, 0);
	Wvector diagHistoryG(3 * dofMap.equationBlockCount, 0);
	Wvector diagHistoryT(3 * dofMap.equationBlockCount, 0);
	Wvector convU(3 * dofMap.historyUBlockCount, 0);
	Wvector convT(3 * dofMap.historyTBlockCount, 0);

	double gmresTime = 0.0;
	long totalIter0 = 0;
	long totalIter1 = 0;
	long solvedSteps = 0;
	int finalFlag = 0;

	for (long step = 1; step <= NStep; ++step)
	{
		printf("MultiDomain GMRES Solver: Step %ld.\n", step);
		BuildKnownInputVectorFromBoundary(dofMap, model, elements, bd[step], states[(size_t)step].globalKnown);
		SparseMul(KnownM, states[(size_t)step].globalKnown, rhs);
		diagKnown = rhs;
		diagHistoryG = 0.0;
		diagHistoryT = 0.0;

		for (long convol = 1; step > 1 && convol <= step; ++convol)
		{
			if (convol > MaxN)
				break;

			BuildConvolutionVector(states, step, convol, 1, convT);
			if (MGSUseSym[convol])
				SparseMul(MGS[convol], convT, temp);
			else
				SparseMul(MGSFull[convol], convT, temp);
			AddVectorInPlace(rhs, temp, 1.0);
			AddVectorInPlace(diagHistoryG, temp, 1.0);

			BuildConvolutionVector(states, step, convol, 0, convU);
			SparseMul(MTS[convol], convU, temp);
			AddVectorInPlace(rhs, temp, -1.0);
			AddVectorInPlace(diagHistoryT, temp, -1.0);
		}

		PrintRhsSummary(step, diagKnown, diagHistoryG, diagHistoryT, rhs, error);
		x0 = 0.0;
		x = 0.0;
		int flag = 1;
		int iter[2] = { 0, 0 };
		clock_t gmresClock = clock();
		gmres(elements, A, rhs, pre, iterations, error, 1, x0, x, flag, iter);
		gmresTime += (double)(clock() - gmresClock) / (double)CLOCKS_PER_SEC;
		PrintLinearResidualSummary("GMRES candidate", step, A, rhs, x, error);
		totalIter0 += iter[0];
		totalIter1 += iter[1];
		finalFlag = flag;
		++solvedSteps;
		printf("MultiDomain GMRES step %ld: flag=%d outer=%d inner=%d\n", step, flag, iter[0], iter[1]);
		WriteRhsBreakdownCsv("multidomain", step, elements, model.elementCount,
			diagKnown, diagHistoryG, diagHistoryT, rhs, x, A);
		if (flag != 0)
			break;

		states[(size_t)step].globalUnknown = x;
		ScatterUnknownToBoundary(dofMap, model, elements, x, bd[step]);
		bd[step].AllLocToAbs(model.nodeCount, DSquareElement::m_transmat);
		BuildPhysicalStateFromBoundary(dofMap, model, elements, bd[step], states[(size_t)step]);
		WriteInterfaceTransferAuditCsv(model, elements, bd[step], step);
		WriteMultiRateStateCsv(model, elements, bd[step], step);
		PrintInterfaceError(dofMap, model, elements, bd[step], step);
	}

	double avgMTSBlocks = MaxN > 0 ? (double)totalHistoryTBlocks / (double)MaxN : 0.0;
	double avgMGSBlocks = MaxN > 0 ? (double)totalHistoryGBlocks / (double)MaxN : 0.0;
	printf("MultiDomain performance: DomainCount=%d GlobalUnknownBlocks=%ld MTS0NonZeroBlocks=%ld MTSAverageBlocks=%lf MGSAverageBlocks=%lf GMRESAverageOuter=%lf GMRESAverageInner=%lf MatrixAssemblyTime=%lf GMRESTime=%lf\n",
		model.DomainCount(),
		dofMap.unknownBlockCount,
		A.NonZeroBlocks(),
		avgMTSBlocks,
		avgMGSBlocks,
		NStep > 0 ? (double)totalIter0 / (double)NStep : 0.0,
		NStep > 0 ? (double)totalIter1 / (double)NStep : 0.0,
		matrixAssemblyTime,
		gmresTime);

	g_validationMetricsValid = 1;
	g_validationFinalFlag = finalFlag;
	g_validationSolvedSteps = solvedSteps;
	g_validationTotalIter0 = totalIter0;
	g_validationTotalIter1 = totalIter1;
	g_validationGlobalUnknownBlocks = dofMap.unknownBlockCount;
	g_validationMTS0Blocks = A.NonZeroBlocks();
	g_validationAvgMTSBlocks = avgMTSBlocks;
	g_validationAvgMGSBlocks = avgMGSBlocks;
	g_validationMatrixAssemblyTime = matrixAssemblyTime;
	g_validationGMRESTime = gmresTime;
	g_validationDomainCount = model.DomainCount();

	delete[] MTS;
	delete[] MGS;
	delete[] MGSFull;
	delete[] MGSUseSym;
	return finalFlag == 0 ? 1 : -1;
}








