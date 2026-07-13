#include "stdio.h"
#include "DSquareElement.h"
#include "Gauss_solver.h"
#include "pre_gmres.h"
#include "aca.h"
#include "write_tec.h"
#include "quadinfo.h"
#include <string>
#include<time.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "DBEM.h"
#include <direct.h>
#include "dsquareelement.h"

#include "PthreadCCSR.h"
#include "multidomain.h"
#include "output_path.h"
//#include<iostream>


#pragma comment(lib, "pthreadVC2.lib")

//全局变量
WmatrixCCSR MCCSRtemp;//临时矩阵
WmatrixSymCCSR MSymCtemp;//临时矩阵
Wmatrix* MT;


DSquareElement* m_DSE;    // 单元类列表
long EleNum;              // 单元数
long NodeNum;             // 物理节点数	
long** m_ElePID;          // 单元几何节点编号
long InfEleNum;
long** m_InfElePID;
long StepNum,ConvoluNum;
double dt_p;
//long* targs;
PthArg* Thread;
PthArg** SpMVThreadT,** SpMVThreadG;
PthConvol* ConvolThreadT,* ConvolThreadG;

long thread_num;
double*** AssistCoef;
double*** AssistUij;
double*** AssistTij;
long* Pthnum;//各行非零子块个数数组    并行负载均衡

WvectorCCSR BCCSRtemp;
SymCCSRMat* MGS;
CCSRMat* MTS;
CCSRMat MTS0;
BoundaryValue* bd;
Wvector* Bvector, CB1, CB2, x0, * CBG, * CBT, Bsum;
long*** FlagUij;

static int DBEMFinite(double value)
{
	return _finite(value) != 0;
}

static void DBEMZeroStress(double stress[3][3])
{
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			stress[i][j] = 0.0;
}

static int DBEMStressComponentsFinite(double stress[3][3])
{
	return DBEMFinite(stress[0][0]) &&
		DBEMFinite(stress[1][1]) &&
		DBEMFinite(stress[2][2]) &&
		DBEMFinite(stress[0][1]) &&
		DBEMFinite(stress[0][2]) &&
		DBEMFinite(stress[1][2]);
}

static void DBEMCleanTecValues(double* values, int count)
{
	for (int i = 0; i < count; ++i)
		if (!DBEMFinite(values[i]))
			values[i] = 0.0;
}

static const char DBEM_BINARY_MAGIC[8] = { 'D','B','E','M','B','I','N','\0' };
static const uint32_t DBEM_BINARY_VERSION = 1;
static const uint32_t DBEM_BINARY_KIND_MODEL = 1;
static const uint32_t DBEM_BINARY_KIND_BD = 2;
static const uint32_t DBEM_BINARY_ENDIAN = 0x01020304;

static void BuildBinaryPath(const char* textPath, char* binPath, size_t binPathSize)
{
	sprintf_s(binPath, binPathSize, "%s.bin", textPath);
}

static int BinaryFileExists(const char* filename)
{
	FILE* fp = 0;
	if (fopen_s(&fp, filename, "rb") != 0 || !fp)
		return 0;
	fclose(fp);
	return 1;
}

static int ReadBinaryHeader(FILE* fp, uint32_t expectedKind, const char* filename)
{
	char magic[8];
	uint32_t version = 0, kind = 0, endian = 0, reserved = 0;
	if (fread(magic, 1, 8, fp) != 8 ||
		fread(&version, sizeof(uint32_t), 1, fp) != 1 ||
		fread(&kind, sizeof(uint32_t), 1, fp) != 1 ||
		fread(&endian, sizeof(uint32_t), 1, fp) != 1 ||
		fread(&reserved, sizeof(uint32_t), 1, fp) != 1)
	{
		printf("Invalid binary header in %s\n", filename);
		return 0;
	}
	if (memcmp(magic, DBEM_BINARY_MAGIC, 8) != 0 || version != DBEM_BINARY_VERSION || kind != expectedKind || endian != DBEM_BINARY_ENDIAN)
	{
		printf("Binary header mismatch in %s\n", filename);
		return 0;
	}
	return 1;
}

static int ReadModelBinary(const char* filename, Point*& pointList, long& pointNum, long**& elePid, char*& eleflag, long& eleNum, long**& infElePid, long& infEleNum, double m_min[3], double m_max[3])
{
	FILE* input = 0;
	if (fopen_s(&input, filename, "rb") != 0 || !input)
	{
		printf("Cannot open binary model file: %s\n", filename);
		return 0;
	}
	if (!ReadBinaryHeader(input, DBEM_BINARY_KIND_MODEL, filename))
	{
		fclose(input);
		return 0;
	}

	int64_t count64 = 0;
	if (fread(&count64, sizeof(int64_t), 1, input) != 1 || count64 <= 0)
	{
		fclose(input);
		printf("Invalid PointNum in binary model file: %s\n", filename);
		return 0;
	}
	pointNum = (long)count64;
	pointList = new Point[pointNum];
	m_max[0] = -10000000000, m_max[1] = -10000000000, m_max[2] = -10000000000;
	m_min[0] = 10000000000, m_min[1] = 10000000000, m_min[2] = 10000000000;
	for (long i = 0; i < pointNum; ++i)
	{
		for (long j = 0; j < 3; ++j)
		{
			if (fread(&pointList[i].pt[j], sizeof(double), 1, input) != 1)
			{
				fclose(input);
				printf("Invalid point data in binary model file: %s\n", filename);
				return 0;
			}
			if (m_max[j] < pointList[i].pt[j])
				m_max[j] = pointList[i].pt[j];
			if (m_min[j] > pointList[i].pt[j])
				m_min[j] = pointList[i].pt[j];
		}
	}

	if (fread(&count64, sizeof(int64_t), 1, input) != 1 || count64 <= 0)
	{
		fclose(input);
		printf("Invalid EleNum in binary model file: %s\n", filename);
		return 0;
	}
	eleNum = (long)count64;
	elePid = new long* [eleNum];
	eleflag = new char[eleNum];
	for (long i = 0; i < eleNum; ++i)
	{
		elePid[i] = new long[8];
		for (long j = 0; j < 8; ++j)
		{
			int64_t nodeId = 0;
			if (fread(&nodeId, sizeof(int64_t), 1, input) != 1)
			{
				fclose(input);
				printf("Invalid element data in binary model file: %s\n", filename);
				return 0;
			}
			elePid[i][j] = (long)nodeId;
		}
		eleflag[i] = char('T');
	}

	if (fread(&count64, sizeof(int64_t), 1, input) != 1 || count64 < 0)
	{
		fclose(input);
		printf("Invalid InfEleNum in binary model file: %s\n", filename);
		return 0;
	}
	infEleNum = (long)count64;
	infElePid = new long* [infEleNum];
	for (long i = 0; i < infEleNum; ++i)
	{
		infElePid[i] = new long[4];
		for (long j = 0; j < 4; ++j)
		{
			int64_t eleId = 0;
			if (fread(&eleId, sizeof(int64_t), 1, input) != 1)
			{
				fclose(input);
				printf("Invalid infinite element data in binary model file: %s\n", filename);
				return 0;
			}
			infElePid[i][j] = (long)eleId;
		}
	}

	fclose(input);
	printf("Read binary model file: %s\n", filename);
	return 1;
}

static int ReadBoundaryBinary(const char* filename, BoundaryValue* bd, int* bdid, long NStep, long EleNum, int checkBoundaryType)
{
	FILE* input = 0;
	if (fopen_s(&input, filename, "rb") != 0 || !input)
	{
		printf("Cannot open binary boundary file: %s\n", filename);
		return 0;
	}
	if (!ReadBinaryHeader(input, DBEM_BINARY_KIND_BD, filename))
	{
		fclose(input);
		return 0;
	}

	int64_t fileEleNum = 0, fileNStep = 0;
	if (fread(&fileEleNum, sizeof(int64_t), 1, input) != 1 ||
		fread(&fileNStep, sizeof(int64_t), 1, input) != 1 ||
		fileEleNum != EleNum || fileNStep != NStep)
	{
		fclose(input);
		printf("Boundary binary size mismatch in %s. file EleNum=%lld NStep=%lld, expected EleNum=%ld NStep=%ld\n", filename, (long long)fileEleNum, (long long)fileNStep, EleNum, NStep);
		return 0;
	}

	int32_t currentBdid = 0;
	for (long i = 0; i < EleNum; ++i)
	{
		if (fread(&currentBdid, sizeof(int32_t), 1, input) != 1)
		{
			fclose(input);
			printf("Invalid boundary type section in %s\n", filename);
			return 0;
		}
		if (checkBoundaryType && currentBdid != bdid[i])
		{
			fclose(input);
			printf("Boundary type changed in %s at element %ld. Matrix cache cannot be reused.\n", filename, i);
			return 0;
		}
		if (!checkBoundaryType)
			bdid[i] = currentBdid;
	}

	for (long i = 0; i <= NStep; ++i)
	{
		bd[i].lu = 0.0;
		bd[i].lt = 0.0;
		for (long j = 0; j < EleNum; ++j)
		{
			double tbd[3];
			if (fread(tbd, sizeof(double), 3, input) != 3)
			{
				fclose(input);
				printf("Invalid boundary value section in %s\n", filename);
				return 0;
			}
				for (long k = 0; k < 8; ++k)
				{
					long j3 = 3 * (j * 8 + k);
					bd[i].lu.b[j3 + 0] = tbd[0];
					bd[i].lu.b[j3 + 1] = tbd[1];
					bd[i].lu.b[j3 + 2] = tbd[2];
				}
		}
	}

	fclose(input);
	printf("Read binary boundary file: %s\n", filename);
	return 1;
}

static int ReadBoundaryValuesFile(const char* filename, BoundaryValue* bd, int* expectedBdid, long NStep, long EleNum)
{
	FILE* input = 0;
	if (fopen_s(&input, filename, "r") != 0 || !input)
	{
		printf("Cannot open boundary file: %s\n", filename);
		return 0;
	}

	int currentBdid = 0;
	for (long i = 0; i < EleNum; ++i)
	{
		if (fscanf_s(input, "%d", &currentBdid) != 1)
		{
			fclose(input);
			printf("Invalid boundary type section in %s\n", filename);
			return 0;
		}
		if (expectedBdid && currentBdid != expectedBdid[i])
		{
			fclose(input);
			printf("Boundary type changed in %s at element %ld. Matrix cache cannot be reused.\n", filename, i);
			return 0;
		}
	}

	double tbd[3];
	for (long i = 0; i <= NStep; ++i)
	{
		bd[i].lu = 0.0;
		bd[i].lt = 0.0;
		for (long j = 0; j < EleNum; ++j)
		{
			if (fscanf_s(input, "%lf", &tbd[0]) != 1 ||
				fscanf_s(input, "%lf", &tbd[1]) != 1 ||
				fscanf_s(input, "%lf", &tbd[2]) != 1)
			{
				fclose(input);
				printf("Invalid boundary value section in %s\n", filename);
				return 0;
			}
				for (long k = 0; k < 8; ++k)
				{
					long j3 = 3 * (j * 8 + k);
					bd[i].lu.b[j3 + 0] = tbd[0];
					bd[i].lu.b[j3 + 1] = tbd[1];
					bd[i].lu.b[j3 + 2] = tbd[2];
				}
		}
	}

	fclose(input);
	return 1;
}

static int ReadBoundaryValuesAny(const char* filename, BoundaryValue* bd, int* expectedBdid, long NStep, long EleNum, int inputBinaryMode)
{
	char binaryPath[300];
	BuildBinaryPath(filename, binaryPath, sizeof(binaryPath));
	if (inputBinaryMode == 1 || inputBinaryMode == 2)
	{
		if (BinaryFileExists(binaryPath))
			return ReadBoundaryBinary(binaryPath, bd, expectedBdid, NStep, EleNum, 1);
		if (inputBinaryMode == 1)
		{
			printf("Missing required binary boundary file: %s\n", binaryPath);
			return 0;
		}
	}
	return ReadBoundaryValuesFile(filename, bd, expectedBdid, NStep, EleNum);
}

static char** AllocateBatchBoundaryFiles(long count)
{
	if (count <= 0)
		return 0;

	char** files = new char* [count];
	for (long i = 0; i < count; ++i)
	{
		files[i] = new char[260];
		files[i][0] = '\0';
	}

	return files;
}

static void FreeBatchBoundaryFiles(char** files, long count)
{
	if (!files)
		return;

	for (long i = 0; i < count; ++i)
	{
		delete[] files[i];
		files[i] = 0;
	}

	delete[] files;
}




using namespace std;

static int ValidationOutputEnabled()
{
	char* value = 0;
	size_t valueSize = 0;
	if (_dupenv_s(&value, &valueSize, "DBEM_VALIDATION_OUTPUT") != 0)
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

static void GetValidationElementCenter(Point* pointList, long** elePid, long ele, double center[3])
{
	center[0] = center[1] = center[2] = 0.0;
	if (!pointList || !elePid || !elePid[ele])
		return;
	for (int i = 0; i < 8; ++i)
	{
		long pid = elePid[ele][i];
		center[0] += pointList[pid].pt[0];
		center[1] += pointList[pid].pt[1];
		center[2] += pointList[pid].pt[2];
	}
	center[0] /= 8.0;
	center[1] /= 8.0;
	center[2] /= 8.0;
}

static void WriteValidationOutputs(DSquareElement* elements,
	Point* pointList,
	long** elePid,
	int* bdid,
	BoundaryValue* bdValue,
	long eleCount,
	long NStep,
	double dt,
	int solver,
	const MultiDomainModel& model)
{
	if (!ValidationOutputEnabled() || !elements || !pointList || !elePid || !bdValue)
		return;

	FILE* state = 0;
	std::string validationStatePath = DBEMOutputPath("validation_state.csv");
	fopen_s(&state, validationStatePath.c_str(), "w");
	if (state)
	{
		fprintf(state, "step,element,localNode,nodeId,domain,surfaceType,bcid,x,y,z,ux,uy,uz,tx,ty,tz\n");
		for (long step = 0; step <= NStep; ++step)
		{
			for (long ele = 0; ele < eleCount; ++ele)
			{
				int bcid = bdid ? bdid[ele] : elements[ele].BCID;
					for (int localNode = 0; localNode < 8; ++localNode)
					{
						long nodeId = elements[ele].m_nodeID[localNode];
						long base = nodeId * 3;
						Point& pt = elements[ele].m_nodelist[nodeId];
					fprintf(state,
						"%ld,%ld,%d,%ld,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
						step,
						ele,
						localNode,
						nodeId,
						elements[ele].DomainID,
						elements[ele].SurfaceType,
						bcid,
						pt.pt[0],
						pt.pt[1],
						pt.pt[2],
						bdValue[step].lu.b[base + 0],
						bdValue[step].lu.b[base + 1],
						bdValue[step].lu.b[base + 2],
						bdValue[step].lt.b[base + 0],
						bdValue[step].lt.b[base + 1],
						bdValue[step].lt.b[base + 2]);
				}
			}
		}
		fclose(state);
	}
	else
	{
		printf("Cannot open %s for writing.\n", validationStatePath.c_str());
	}

	FILE* metrics = 0;
	std::string validationMetricsPath = DBEMOutputPath("validation_metrics.txt");
	fopen_s(&metrics, validationMetricsPath.c_str(), "w");
	if (!metrics)
	{
		printf("Cannot open %s for writing.\n", validationMetricsPath.c_str());
		return;
	}

	fprintf(metrics, "ValidationOutput=1\n");
	fprintf(metrics, "Solver=%d\n", solver);
	fprintf(metrics, "NStep=%ld\n", NStep);
	fprintf(metrics, "Dt=%.17g\n", dt);
	fprintf(metrics, "ElementCount=%ld\n", eleCount);
		fprintf(metrics, "NodeCount=%ld\n", eleCount * 8);
	fprintf(metrics, "MultiDomainActive=%d\n", model.IsActive() ? 1 : 0);
	fprintf(metrics, "DomainCount=%d\n", model.DomainCount());
	fprintf(metrics, "InterfaceCount=%ld\n", (long)model.interfaces.size());
	for (size_t d = 0; d < model.domains.size(); ++d)
	{
		const Domain& domain = model.domains[d];
		fprintf(metrics,
			"DomainMaterial[%d]=C1 %.17g C2 %.17g C1Dt %.17g C2Dt %.17g\n",
			domain.id,
			domain.mat.C1,
			domain.mat.C2,
			domain.mat.C1Dt,
			domain.mat.C2Dt);
	}

	double overallMaxU = 0.0;
	double overallMaxT = 0.0;
	for (long step = 0; step <= NStep; ++step)
	{
		double maxU = 0.0;
		double maxT = 0.0;
		for (size_t i = 0; i < model.interfaces.size(); ++i)
		{
			const InterfacePair& itf = model.interfaces[i];
			int tSignB = itf.normalSign == 0 ? -1 : itf.normalSign;
				for (int localNode = 0; localNode < 8; ++localNode)
				{
					int localB = itf.localBForA[localNode];
					if (localB < 0 || localB >= 8)
						localB = localNode;
					long a = (itf.eleA * 8 + localNode) * 3;
					long b = (itf.eleB * 8 + localB) * 3;
					for (int k = 0; k < 3; ++k)
					{
					double du = fabs(bdValue[step].lu.b[a + k] - bdValue[step].lu.b[b + k]);
					double tb = fabs(bdValue[step].lt.b[b + k] - (double)tSignB * bdValue[step].lt.b[a + k]);
					if (du > maxU)
						maxU = du;
					if (tb > maxT)
						maxT = tb;
				}
			}
		}
		if (maxU > overallMaxU)
			overallMaxU = maxU;
		if (maxT > overallMaxT)
			overallMaxT = maxT;
		fprintf(metrics, "InterfaceStep[%ld].MaxAbsU=%.17g\n", step, maxU);
		fprintf(metrics, "InterfaceStep[%ld].MaxAbsTBalance=%.17g\n", step, maxT);
	}
	fprintf(metrics, "InterfaceOverallMaxAbsU=%.17g\n", overallMaxU);
	fprintf(metrics, "InterfaceOverallMaxAbsTBalance=%.17g\n", overallMaxT);
	if (model.IsActive())
		WriteMultiDomainValidationMetrics(metrics);
	else
		fprintf(metrics, "MultiDomainMetricsValid=0\n");
	fclose(metrics);
}
//int DBEM_main() { return 0; };
int main()
{
	//--------------变量定义----------------
	int m_flag = 0;
	long PointNum;            // 几何节点数
	//long EleNum;              // 单元数   全局
	//long NodeNum;             // 物理节点数	  全局
	Point* m_PointList = 0;       // 几何节点列表
	//long** m_ElePID;          // 单元几何节点编号   全局
	Point* m_NodeList = 0;        // 物理节点列表
	long** m_EleNID = 0;          // 单元物理节点编号

	double dt;                // 时间步长
	double Beta;              // Beta = c1 * dt / AveElementSize
	long NStep;               // 时间步数 n = 0 .. NStep
	int* bdid;                // 单元边界条件代码
	//BoundaryValue* bd;        // 单元边界条件值    并行改全局

	//DSquareElement* m_DSE;    // 单元类列表 全局

	double gap;               // 单元不连续参数 0 < gap <= 1
	double v;                 // 泊松比
	double G;                 // 剪切模量
	double E;                 // 弹性模量
	double Rou;               // 密度
	double amplitude;         // 输出变形图时的变形放大系数
	int solvers;              // 求解器类型
	double error;             // 迭代收敛残差
	long iterations;          // 最大迭代步
	long maxleafelenum;       // ACA叶子节点允许最大element数量
	long premaxleafpointnum;  // 预处理ACA叶子节点允许最大node数量

	int Flag_ReadWrite;       // 0: ACA不进行读写  1: ACA进行读写

	int FlagDyna;             // 1: Standard  2: Interpolation Method

	long i, j, k;               // 循环变量

	double errorvalue = 0.00000001;
	double errorstop = 0.0001;
	char* eleflag;//面类型
	char temptitlename[50], titlename[260], modelname[260], bdname[260], YN, SD;//文件名
	int InputBinaryMode = 0;
	int BatchBoundaryMode = 0;
	long BatchBoundaryCaseCount = 1;
	char** BatchBoundaryFiles = 0;
	MultiDomainInputConfig MultiDomainConfig;
	MultiDomainModel MultiDomain;
	FILE* logfile;
	//—————————返回路径—————————
	char* CurrPath = _getcwd(NULL, 0);//返回当前路径
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//返回当前文件夹名称
	Path = StrCurrPath.substr(0, found + 1);//返回上级路径

	string OutFile = "output\\";
	string TempPath = Path + OutFile;

	

	string Temp = "log.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();
	//——————————————————————————
	clock_t starttime, endtime, temptime;
	//--------------参数卡片----------------

	fopen_s(&logfile, OutPath, "w");
	printf_s("读入参数卡片 BEM_DATACARD.DAT...\n\n");
	FILE* input;
	starttime = clock();

	fopen_s(&input, "BEM_DATACARD.DAT", "r");
	fscanf_s(input, "%s", temptitlename, 19);
	strcpy_s(titlename, ".//input//");
	strcat_s(titlename, temptitlename);
	strcpy_s(modelname, titlename);
	strcat_s(modelname, ".model");
	strcpy_s(bdname, titlename);
	strcat_s(bdname, ".bd");
	printf_s("Titlename is %s .\n", titlename);
	fprintf_s(logfile,"Titlename is %s .\n", titlename);
	fscanf_s(input, "%d", &thread_num);
	fscanf_s(input, "%d", &FlagDyna);
	fscanf_s(input, "%lf", &Beta);
	fscanf_s(input, "%ld", &NStep);
	DBEMInitializeOutputDirectory(temptitlename, NStep);
	fscanf_s(input, "%lf", &gap);
	fscanf_s(input, "%lf", &E);
	fscanf_s(input, "%lf", &v);
	fscanf_s(input, "%lf", &Rou);
	fscanf_s(input, "%lf", &amplitude);
	fscanf_s(input, "%d", &solvers);
	if (solvers == 1)
	{
		long legacyPos = ftell(input);
		long legacyIterations = 0;
		double legacyError = 0.0;
		long legacyPremax = 0;
		int legacyInputMode = 0;
		int legacyBatchMode = 0;
		int legacyRead = fscanf_s(input, "%ld", &legacyIterations);
		legacyRead += fscanf_s(input, "%lf", &legacyError);
		legacyRead += fscanf_s(input, "%ld", &legacyPremax);
		legacyRead += fscanf_s(input, "%d", &legacyInputMode);
		legacyRead += fscanf_s(input, "%d", &legacyBatchMode);
		if (legacyRead == 5 && legacyIterations > 0 && legacyError > 0.0 && legacyPremax >= 0 &&
			legacyInputMode >= 0 && legacyInputMode <= 2 && legacyBatchMode >= 0 && legacyBatchMode <= 2)
		{
			printf("Gauss solver ignored legacy GMRES parameter block.\n");
			fprintf_s(logfile, "Gauss solver ignored legacy GMRES parameter block.\n");
		}
		else
		{
			fseek(input, legacyPos, SEEK_SET);
		}
	}
	if (solvers == 2)
	{
		fscanf_s(input, "%ld", &iterations);
		fscanf_s(input, "%lf", &error);
		fscanf_s(input, "%ld", &premaxleafpointnum);
		long PreReadBatchCaseCount = -1;
		if (fscanf_s(input, "%d", &InputBinaryMode) != 1)
			InputBinaryMode = 0;
		if (InputBinaryMode < 0 || InputBinaryMode > 2)
		{
			BatchBoundaryMode = InputBinaryMode;
			InputBinaryMode = 0;
		}
		else if (fscanf_s(input, "%d", &BatchBoundaryMode) != 1)
		{
			BatchBoundaryMode = 0;
		}
		else if (BatchBoundaryMode < 0 || BatchBoundaryMode > 2)
		{
			PreReadBatchCaseCount = BatchBoundaryMode;
			BatchBoundaryMode = InputBinaryMode;
			InputBinaryMode = 0;
		}

		if (BatchBoundaryMode == 1)
		{
			if (PreReadBatchCaseCount > 0)
				BatchBoundaryCaseCount = PreReadBatchCaseCount;
			else
				fscanf_s(input, "%ld", &BatchBoundaryCaseCount);
			BatchBoundaryFiles = AllocateBatchBoundaryFiles(BatchBoundaryCaseCount);
			for (i = 0; i < BatchBoundaryCaseCount; ++i)
				fscanf_s(input, "%s", BatchBoundaryFiles[i], 260);
		}
		else if (BatchBoundaryMode == 2)
		{
			char BatchBoundaryDir[260] = { 0 };
			char BatchBoundaryPrefix[260] = { 0 };
			char BatchBoundaryExt[32] = { 0 };
			long BatchBoundaryStartIndex = 0;
			long BatchBoundaryEndIndex = 0;

			if (PreReadBatchCaseCount > 0)
				BatchBoundaryCaseCount = PreReadBatchCaseCount;
			else
				fscanf_s(input, "%ld", &BatchBoundaryCaseCount);
			fscanf_s(input, "%s", BatchBoundaryDir, (unsigned)_countof(BatchBoundaryDir));
			fscanf_s(input, "%s", BatchBoundaryPrefix, (unsigned)_countof(BatchBoundaryPrefix));
			fscanf_s(input, "%ld", &BatchBoundaryStartIndex);
			fscanf_s(input, "%ld", &BatchBoundaryEndIndex);
			fscanf_s(input, "%s", BatchBoundaryExt, (unsigned)_countof(BatchBoundaryExt));

			long RangeCount = BatchBoundaryEndIndex - BatchBoundaryStartIndex + 1;
			if (BatchBoundaryCaseCount <= 0 || RangeCount <= 0)
			{
				printf("Invalid GMRES batch boundary directory config.\n");
				BatchBoundaryMode = 0;
				BatchBoundaryCaseCount = 0;
			}
			else
			{
				if (RangeCount != BatchBoundaryCaseCount)
					printf("GMRES batch boundary count mismatch: count=%ld, range=%ld. Use count.\n", BatchBoundaryCaseCount, RangeCount);
				BatchBoundaryFiles = AllocateBatchBoundaryFiles(BatchBoundaryCaseCount);
				for (i = 0; i < BatchBoundaryCaseCount; ++i)
					sprintf_s(BatchBoundaryFiles[i], 260, "%s%s%ld%s", BatchBoundaryDir, BatchBoundaryPrefix, BatchBoundaryStartIndex + i, BatchBoundaryExt);
			}
		}
		else
		{
			BatchBoundaryMode = 0;
			BatchBoundaryCaseCount = 0;
		}
	}
	if (solvers == 3)
	{
		fscanf_s(input, "%ld", &iterations);
		fscanf_s(input, "%lf", &error);
		fscanf_s(input, "%ld", &maxleafelenum);
		fscanf_s(input, "%ld", &premaxleafpointnum);
		fscanf_s(input, "%lf", &errorvalue);
		fscanf_s(input, "%lf", &errorstop);
		fscanf_s(input, "%d", &Flag_ReadWrite);
	}

	if ((BatchBoundaryMode == 1 || BatchBoundaryMode == 2) && BatchBoundaryCaseCount > 0 && BatchBoundaryFiles)
	{
		printf("Input binary mode %d\n", InputBinaryMode);
		printf("GMRES batch boundary mode %d, cases %ld\n", BatchBoundaryMode, BatchBoundaryCaseCount);
		printf("GMRES batch first case: %s\n", BatchBoundaryFiles[0]);
		printf("GMRES batch last case: %s\n", BatchBoundaryFiles[BatchBoundaryCaseCount - 1]);
		fprintf_s(logfile, "Input binary mode %d\n", InputBinaryMode);
		fprintf_s(logfile, "GMRES batch boundary mode %d, cases %ld\n", BatchBoundaryMode, BatchBoundaryCaseCount);
		fprintf_s(logfile, "GMRES batch first case: %s\n", BatchBoundaryFiles[0]);
		fprintf_s(logfile, "GMRES batch last case: %s\n", BatchBoundaryFiles[BatchBoundaryCaseCount - 1]);
	}

	int MultiDomainRead = ReadOptionalMultiDomainInput(input, MultiDomainConfig);
	if (MultiDomainRead < 0)
	{
		FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
		return -1;
	}
	if (MultiDomainConfig.enabled)
	{
		printf("MultiDomain config enabled, DomainCount = %d\n", MultiDomainConfig.domainCount);
		fprintf_s(logfile, "MultiDomain config enabled, DomainCount = %d\n", MultiDomainConfig.domainCount);
	}

	fclose(input);

	G = E / 2.0 / (1 + v);




	AssistCoef = new double** [thread_num];
	AssistUij = new double** [thread_num];
	AssistTij = new double** [thread_num];
	for (i = 0; i < thread_num; ++i) {
		AssistCoef[i] = new double* [8];
		AssistUij[i] = new double* [8];
		AssistTij[i] = new double* [8];
		for (j = 0; j < 8; ++j) {
			AssistCoef[i][j] = new double[9];
			AssistUij[i][j] = new double[9];
			AssistTij[i][j] = new double[9];
		}

	}

	for (i = 0; i < thread_num; ++i) {
		for (j = 0; j < 8; ++j) {
			for (k = 0; k < 9; ++k) {
				AssistCoef[i][j][k] = 0.0;
				AssistUij[i][j][k] = 0.0;
				AssistTij[i][j][k] = 0.0;
			}
		}
	}



	//--------------读入BEM_INP.DAT----------------
	double m_max[3];
	double m_min[3];
	char binaryModelName[300];
	BuildBinaryPath(modelname, binaryModelName, sizeof(binaryModelName));
	if (InputBinaryMode == 1 || (InputBinaryMode == 2 && BinaryFileExists(binaryModelName)))
	{
		if (!ReadModelBinary(binaryModelName, m_PointList, PointNum, m_ElePID, eleflag, EleNum, m_InfElePID, InfEleNum, m_min, m_max))
		{
			if (InputBinaryMode == 1)
			{
				FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
				return -1;
			}
		}
	}
	if (!m_PointList)
	{
		if (InputBinaryMode == 1)
		{
			printf("Missing required binary model file: %s\n", binaryModelName);
			FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
			return -1;
		}
		fopen_s(&input, modelname, "r");
		// 节点
		fscanf_s(input, "%ld", &PointNum);
		m_PointList = new Point[PointNum];
		//取得区域大小
		m_max[0] = -10000000000, m_max[1] = -10000000000, m_max[2] = -10000000000;
		m_min[0] = 10000000000, m_min[1] = 10000000000, m_min[2] = 10000000000;
		for (i = 0; i < PointNum; ++i)
			for (j = 0; j < 3; ++j)
			{
				fscanf_s(input, "%lf", &m_PointList[i].pt[j]);
				//km  m_PointList[i].pt[j] *= 1000;
				if (m_max[j] < m_PointList[i].pt[j])
					m_max[j] = m_PointList[i].pt[j];
				if (m_min[j] > m_PointList[i].pt[j])
					m_min[j] = m_PointList[i].pt[j];
			}


		// 单元
		fscanf_s(input, "%ld", &EleNum);
		m_ElePID = new long*[EleNum];
		eleflag = new char[EleNum];
		//	printf_s("%d",eleflag.size());
		for (i = 0; i < EleNum; ++i)
		{
			m_ElePID[i] = new long[8];
			for (j = 0; j < 8; ++j)
			{
				fscanf_s(input, "%ld", &m_ElePID[i][j]);
				m_ElePID[i][j] --;//从1开始编节点,如果输入文件从0开始则注释掉
			}
			//fscanf_s(input, " %c", &eleflag[i]);//根据需求设置,暂删除
			eleflag[i] = char('T');//根据需求设置,暂添加
		}

		//printf_s("%d", EleNum);
		//读入无限单元
		//long InfEleNum; 全局
		//long ** m_InfElePID;全局
		fscanf_s(input, "%ld", &InfEleNum);//无限单元个数
		//InfEleNum = 0;//暂设置
		m_InfElePID = new long*[InfEleNum];
		for (i = 0; i < InfEleNum; i++) {
			m_InfElePID[i] = new long[4];
			for (j = 0; j < 4; j++) {
				fscanf_s(input, "%ld", &m_InfElePID[i][j]);
				m_InfElePID[i][j]--; //王老师程序为从1开始编节点,如果输入文件从0开始则注释掉
			}
		}

		fclose(input);
	}
	double m_z[3];
	//自动设置区域中心为z点
	for (i = 0; i < 3; i++) {
		m_z[i] = (m_max[i] + m_min[i]) / 2.0;
	}
	//手动读入z点坐标，暂删去
	//for (i = 0; i < 3; i++) {
	//	fscanf_s(input, "%ld", &m_z[i]);
	//}
	// 获取单元平均尺寸
	double AveEleSize = GetEleSize(m_PointList, m_ElePID, PointNum, EleNum);
	double c1 = sqrt((2.0 - 2.0 * v) / (1.0 - 2.0 * v) * G / Rou);
	dt = Beta * AveEleSize / c1;//0.078125000000000*2;// Beta * AveEleSize / c1;//调整一步大小*****
	//dt = Beta * AveEleSize / c1;
	double tempbeta = c1 * dt / AveEleSize;
	printf("\nBeta change to %lf\n\n", tempbeta);
	//printf("请检查beta是否合适(Y,N)？\n");
	//scanf_s(" %c", &YN, 1);
	//if (YN == 'n' || YN == 'N')
	//{
	//	printf("请调整...\n");
	//	return -1;
	//}
	//else if (YN == 'y' || YN == 'Y')
	//{
	//	printf("正常运行...\n");
	//}
	//else {
	//	printf("不可接收的参数,程序退出...\n");
	//	return -1;
	//}

	printf("Element Size: %lf\n\n", AveEleSize);
	fprintf_s(logfile, "Element Size: %lf\n\n", AveEleSize);
	printf("dt: %2.15lf\n\n", dt);
	fprintf_s(logfile, "dt: %2.15lf\n\n", dt);
	FILE* dtpf;
	Temp = "dt.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();
	fopen_s(&dtpf, OutPath, "w");
	fprintf(dtpf, "%2.15lf", dt);
	fclose(dtpf);

	if (!BuildMultiDomainModel(MultiDomainConfig, EleNum, EleNum * 8, E, v, Rou, dt, MultiDomain))
	{
		FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
		return -1;
	}
	if (MultiDomain.IsActive())
	{
		MultiDomain.PrintSummary(stdout);
		MultiDomain.PrintSummary(logfile);
	}

	printf("请输入求解问题类型(S静力学,D动力学):\n");
	//暂删去scanf_s(" %c", &SD, 1);
	SD = 'd';//暂使用
	if (SD == 's' || SD == 'S') {
		FlagDyna = 1;
		DSquareElement::FlagDynaMethod = -1;
		DSquareElement::Flagproblem = 1;
		printf("静力学问题\n");
		fprintf_s(logfile, "静力学问题\n");
	}
	else if (SD == 'd' || SD == 'D')
	{
		DSquareElement::Flagproblem = 0;
		printf("动力学问题\n");
		fprintf_s(logfile, "动力学问题\n");
	}
	else {
		printf("不可接收的参数,程序退出...\n");
		fprintf_s(logfile, "不可接收的参数,程序退出...\n");
		FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
		return -1;
	}

	printf("\n输入信号文件是否完成(Y完成,N未完成)？\n");
	//scanf_s(" %c", &YN, 1);
	YN = 'y';//暂使用
	if (YN == 'n' || YN == 'N')
	{
		printf("请调用bdmaker...\n");
		FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
		return -1;
	}
	else if (YN == 'y' || YN == 'Y')
	{
		printf("读入输入信号文件...\n");
	}
	else {
		printf("不可接收的参数,程序退出...\n");
		FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
		return -1;
	}

	// 边界条件	
//	fscanf_s(input,"%lf",&dt);
//	fscanf_s(input,"%ld",&NStep);
	if ((BatchBoundaryMode == 1 || BatchBoundaryMode == 2) && BatchBoundaryCaseCount > 0)
		strcpy_s(bdname, BatchBoundaryFiles[0]);
	//王老师要保证是偶数步,不知道是为啥,暂删去
	//	if ((NStep / 2) * 2 != NStep)
	//	{
	//		NStep++;
	//	}
	printf("NStep = %ld\n",NStep);
	bdid = new int[EleNum];



	bd = new BoundaryValue[NStep + 1];
	for (i = 0; i <= NStep; ++i)
		bd[i].init(EleNum * 8); // 物理节点边界条件变量初始化

	double tbd[3];
	long j24, k3;

	// 注意：给定边界条件是单元上的，而计算结果是节点的

	char binaryBdName[300];
	BuildBinaryPath(bdname, binaryBdName, sizeof(binaryBdName));
	if (InputBinaryMode == 1 || (InputBinaryMode == 2 && BinaryFileExists(binaryBdName)))
	{
		if (!ReadBoundaryBinary(binaryBdName, bd, bdid, NStep, EleNum, 0))
		{
			if (InputBinaryMode == 1)
			{
				FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
				return -1;
			}
		}
	}
	else
	{
		if (InputBinaryMode == 1)
		{
			printf("Missing required binary boundary file: %s\n", binaryBdName);
			FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
			return -1;
		}
		fopen_s(&input, bdname, "r");
		for (i = 0; i < EleNum; ++i)
		{
			fscanf_s(input, "%d", &bdid[i]);

		}
		//自动对应时间步长给边界条件,如不用则注释
		//double m_A0 = 1.0, m_f0=1.0,m_omega = 2.0*pi * 1.0;
		
		FILE* sourcedata;
		//fopen_s(&sourcedata, "source.dat", "r");
		//fscanf_s(sourcedata, "%lf", &m_A0);
		//fscanf_s(sourcedata, "%lf", &m_f0);
		//m_omega = 2.0*pi * m_f0;
		//fclose(sourcedata);
		  

		FILE* TIMEDATE;

		Temp = "timedate.txt";
		Path = TempPath + Temp;
		OutPath = Path.c_str();
		fopen_s(&TIMEDATE, OutPath, "w");
		int flag_m = 0;

		for (i = 0; i <= NStep; ++i)
		{
			flag_m = 0;
			for (j = 0; j < EleNum; ++j)
			{
				fscanf_s(input, "%lf", &tbd[0]);
				fscanf_s(input, "%lf", &tbd[1]);
				fscanf_s(input, "%lf", &tbd[2]);
				j24 = 24 * j;
				for (long k = 0; k < 8; ++k)
				{
					k3 = j24 + 3 * k;
					bd[i].lu.b[k3 + 0] = tbd[0];
					bd[i].lu.b[k3 + 1] = tbd[1];
					bd[i].lu.b[k3 + 2] = tbd[2];
				}
			}
		}
		fclose(TIMEDATE);
		fclose(input);
	}

	
	//-------------单元初始化----------------

	// 生成单元类列表
	m_DSE = new DSquareElement[EleNum];

	// 单元静态初始化
	DSquareElement::StaticInit(gap, v, G, Rou, dt, FlagDyna);
	ApplyMultiDomainModelToElements(MultiDomain, m_DSE, EleNum);

	// TEST
	const DynaMat& logMat = MultiDomain.IsActive() ? MultiDomain.materialContext.Get(0) : DSquareElement::m_DMat;
	double logE = MultiDomain.IsActive() ? 2.0 * logMat.G * (1.0 + logMat.v) : E;
	double logV = logMat.v;
	double dudt = logMat.C1 * sqrt((1.0 + logV)*(1.0 - 2.0*logV) / (1.0 - logV)) / logE * dt / dt;
	double C_1D = logMat.C1 * sqrt((1.0 + logV)*(1.0 - 2.0*logV) / (1.0 - logV));
	printf("c1: %lf   c2: %lf  c1*dt: %lf\n", logMat.C1, logMat.C2, logMat.C1Dt);
	printf("c: %lf u = %lf\n", C_1D, dudt * dt);
	fprintf_s(logfile, "c1: %lf   c2: %lf  c1*dt: %lf\n", logMat.C1, logMat.C2, logMat.C1Dt);
	fprintf_s(logfile, "c: %lf u = %lf\n", C_1D, dudt * dt);
		fflush(stdout); fflush(logfile);
	// 单元几何初始化

	#pragma omp parallel for num_threads(thread_num)
	for (long p = 0; p < EleNum; ++p)
	{
		m_DSE[p].GeoInit(m_PointList, m_ElePID[p]);   //      计算各单元雅可比、法向量、积分点坐标。   占空间    需改为积分时再计算
		fprintf_s(logfile, "Trace: GeoInit loop reached\n"); fflush(logfile);
	}



	//标记无限单元的相邻单元
	for (i = 0; i < InfEleNum; i++) {
		m_DSE[m_InfElePID[i][0]].InfEleFlag = i;  //无限单元暂不修改
	}
	// 生成物理节点和物理单元节点编号列表
	m_EleNID = new long*[EleNum];
	for (i = 0; i < EleNum; ++i)
		m_EleNID[i] = new long[8];

	// 根据单元几何信息确定不连续信息
	GetNodeInfo(m_PointList, m_NodeList, m_ElePID, m_EleNID, PointNum, EleNum, NodeNum, m_DSE);
		fprintf_s(logfile, "Trace: GetNodeInfo finished NodeNum=%ld EleNum=%ld\n", NodeNum, EleNum); fflush(logfile);
	if (MultiDomain.IsActive() && !ResolveMultiDomainInterfaceNodeMaps(MultiDomain, m_NodeList, NodeNum))
	{
		FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
		return -1;
	}

	// 计算节点所属单元信息
	DSquareElement::StaticGetBeInfo(NodeNum, EleNum, m_EleNID);
		fprintf_s(logfile, "Trace: StaticGetBeInfo finished\n"); fflush(logfile);
	
	// 单元物理初始化
	if (MultiDomain.IsActive())
	{
		DynaMat savedStaticMat = DSquareElement::m_DMat;
		for (long p = 0; p < EleNum; ++p) {
			DSquareElement::m_DMat = MultiDomain.materialContext.Get(m_DSE[p].DomainID);
			m_DSE[p].PthPhyInit(m_NodeList, m_EleNID[p], bdid[p]);     //计算各单元奇异积分并存储。  减少存储量
		}
		DSquareElement::m_DMat = savedStaticMat;
	}
	else
	{
		#pragma omp parallel for num_threads(thread_num)
		for (long p = 0; p < EleNum; ++p) {
			m_DSE[p].PthPhyInit(m_NodeList, m_EleNID[p], bdid[p]);     //计算各单元奇异积分并存储。  减少存储量
		}
	}


		fprintf_s(logfile, "Trace: PthPhyInit finished\n"); fflush(logfile);
	//无限单元的参考点
	DSquareElement::InfCenter[0] =0.0;// m_z[0]
	DSquareElement::InfCenter[1] = 0.0;//m_z[1]
	DSquareElement::InfCenter[2] = 1.5;//m_z[2]
	printf_s("Center z is (%lf,%lf,%lf)\n", m_z[0], m_z[1], m_z[2]);
	fprintf_s(logfile, "Center z is (%lf,%lf,%lf)\n", m_z[0], m_z[1], m_z[2]);
	// 把不连续单元画出来
//	Plot(m_NodeList,m_EleNID,"dis-element.dat",NodeNum,EleNum,2,amplitude);
	
	endtime = clock();
	printf_s("读入文件用时 %ld\n", endtime - starttime);
	fprintf_s(logfile, "读入文件用时 %ld\n", endtime - starttime);
	temptime = endtime;




	//OutPoint* OutP;
	//OutEle* OutE;
	//long numOut = PointNum;
	//OutP = new OutPoint[numOut];
	//OutE = new OutEle[EleNum];
	//PreOutput(OutP, OutE, m_DSE, EleNum, PointNum, m_PointList,numOut);
	
	//-------------选择求解问题类型------------------

	if (SD == 'D' || SD == 'd')
	{
		printf("调用动力学求解...\n");
		fprintf_s(logfile, "Trace: before dynamic solve\n"); fflush(logfile);
		//-------------选择求解器进行求解----------------

		Point* m_EleCenterList;

		double MaxLe = 1000;//MaxLength_adapt(m_PointList, PointNum);//大规模算例默认给出

		printf("Select Solvers:\n1: Gauss Elimination\n2: GMRES\n3: ACA\n\nPlease Enter your choice:");

		printf("solver = %d\n", solvers);

		if (solvers == 1)
		{
			DSquareElement::FlagDynaMethod == 2;//强制稳定性处理
			if (DSquareElement::FlagDynaMethod == 1)
			{
				DynaGaussSolver(m_DSE, bd, NodeNum, EleNum, NStep, MaxLe, m_InfElePID, m_ElePID);
				//DynaGaussSolver(m_DSE, bd, NodeNum, EleNum, NStep, MaxLe);
			}

			else//修改中
			{
				DynaGaussSolverNew(m_DSE, bd, NodeNum, EleNum, NStep, MaxLe, m_InfElePID, m_ElePID, m_flag);//输出各矩阵
				//DynaGaussSolverNewCSR(m_DSE, bd, NodeNum, EleNum, NStep, MaxLe, m_InfElePID, m_ElePID, m_flag);
				//DynaGaussSolverNewCCSR(m_DSE, bd, NodeNum, EleNum, NStep, MaxLe, m_InfElePID, m_ElePID, m_flag);
				
				//DynaGaussSolverNew(m_DSE, bd, NodeNum, EleNum, NStep, MaxLe);
			}
				
			printf("end gauss\n");
		}
		else if (solvers == 2)
		{
			DSquareElement::FlagDynaMethod = 2;//强制稳定性处理
			if (MultiDomain.IsActive())
			{
				if ((BatchBoundaryMode == 1 || BatchBoundaryMode == 2) && BatchBoundaryCaseCount > 0)
				{
					printf("MultiDomain GMRES does not support batch boundary mode yet.\n");
					FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
					return -1;
				}
				if (DynaGMRESSolverMultiDomainCCSR(m_DSE, bd, MultiDomain, iterations, error, premaxleafpointnum, NStep, MaxLe, m_InfElePID, m_ElePID) < 0)
				{
					FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
					return -1;
				}
			}
			else if ((BatchBoundaryMode == 1 || BatchBoundaryMode == 2) && BatchBoundaryCaseCount > 0 && DSquareElement::FlagDynaMethod != 1)
			{
				GmresCcsrBeginMemorySession();
				for (long caseID = 0; caseID < BatchBoundaryCaseCount; ++caseID)
				{
					if (caseID > 0 && !ReadBoundaryValuesAny(BatchBoundaryFiles[caseID], bd, bdid, NStep, EleNum, InputBinaryMode))
					{
						GmresCcsrEndMemorySession();
						FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
						return -1;
					}
					char caseDir[64];
					sprintf_s(caseDir, "case_%03ld", caseID + 1);
					GmresCcsrSetOutputSubdir(caseDir);
					printf("GMRES batch boundary case %ld/%ld: %s\n", caseID + 1, BatchBoundaryCaseCount, BatchBoundaryFiles[caseID]);
					DynaGMRESSolverNewCCSRNewPth(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe, m_InfElePID, m_ElePID);
				}
				GmresCcsrSetOutputSubdir("");
				GmresCcsrEndMemorySession();
			}
			else if (DSquareElement::FlagDynaMethod == 1)
				//DynaGMRESSolver(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe);
				DynaGMRESSolverNewCCSRNewPth(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe, m_InfElePID, m_ElePID);

			else
				//DynaGMRESSolverNew(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe);
				//DynaGMRESSolverNew(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe, m_InfElePID, m_ElePID);
				////DynaGMRESSolverNewCSR(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe, m_InfElePID, m_ElePID);
				//DynaGMRESSolverNewCCSR(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe, m_InfElePID, m_ElePID);//普通并行
				DynaGMRESSolverNewCCSRNewPth(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe, m_InfElePID, m_ElePID);//改进并行
			printf("end GMRES\n");
		}
		else if (solvers == 3)
		{
			//ACA
			DSquareElement::FlagDynaMethod = 2;//强制稳定性处理
			m_EleCenterList = new Point[EleNum];
			for (i = 0; i < EleNum; ++i)
			{
				m_DSE[i].GetFromLocal(0.0, 0.0, m_EleCenterList[i]);
			}

			if (DSquareElement::FlagDynaMethod == 1)
				DynaACASolver(m_DSE, m_NodeList, m_EleCenterList, bd, NodeNum, EleNum, iterations, maxleafelenum, premaxleafpointnum, error, errorvalue, errorstop, NStep, MaxLe, Flag_ReadWrite);
			else if (DSquareElement::FlagDynaMethod == 2)
				DynaACASolverNew(m_DSE, m_NodeList, m_EleCenterList, bd, NodeNum, EleNum, iterations, maxleafelenum, premaxleafpointnum, error, errorvalue, errorstop, NStep, MaxLe, Flag_ReadWrite);
			else
				DynaACASolverMem(m_DSE, m_NodeList, m_EleCenterList, bd, NodeNum, EleNum, iterations, maxleafelenum, premaxleafpointnum, error, errorvalue, errorstop, NStep, MaxLe, Flag_ReadWrite);


			delete[] m_EleCenterList;
			m_EleCenterList = 0;

			printf("end ACA\n");
		}

		endtime = clock();	
		printf_s("求解用时 %ld\n", endtime - temptime);
		fprintf_s(logfile, "求解用时 %ld\n", endtime - temptime);
		temptime = endtime;

		//tecplot
		for(i=0;i<=NStep;++i)
			ResultPlotNodeAverage_14VARS(m_DSE,m_PointList,m_ElePID,bdid,bd[i],i,PointNum,NodeNum,EleNum,DBEMOutputPath("TecValueFile").c_str(),amplitude);


		// ------输出选定节点的位移/面力随时间的变化信息--------
		//GetInfoFromPoint_All(m_NodeList, bd, m_EleNID, EleNum, eleflag, NStep);
		//GetInfoFromPoint_onenode(m_NodeList, bd, NodeNum, NStep, dt, m_EleNID);
		//GetInfoFromPoint_OutPut(OutP, OutE, bd, m_EleNID, numOut, EleNum, NStep);
		//GetInfoFromPoint_Uppersurface(m_NodeList, bd, m_EleNID, EleNum, eleflag, NStep);
		//GetInfoFromPoint_SingleColumn(m_NodeList, bd, NodeNum, NStep, dt, dudt, C_1D, 4.0, m_EleNID);
		
		endtime = clock();
		printf_s("后处理用时 %ld\n", endtime - temptime);
		fprintf_s(logfile, "后处理用时 %ld\n", endtime - temptime);
		temptime = endtime;
	}
	else if (SD == 'S' || SD == 's')
	{
		printf("调用静力学求解...\n");

		//-------------选择求解器进行求解----------------

	//		Point* m_EleCenterList;

		double MaxLe = MaxLength(m_PointList, PointNum);

		printf("Select Solvers:\n1: Gauss Elimination\n2: GMRES\n3: ACA\n\nPlease Enter your choice:");

		printf("solver = %d\n", solvers);

		if (solvers == 1)
		{
			StaticGaussSolver(m_DSE, bd, NodeNum, EleNum, InfEleNum, m_InfElePID, m_ElePID, MaxLe);//增加适用无限单元
			//StaticGaussSolver(m_DSE, bd, NodeNum, EleNum, MaxLe);
			printf("end gauss\n");
		}
		else if (solvers == 2)
		{
			//GMRESSolver(m_DSE, bd, NodeNum, EleNum, iterations, error, premaxleafpointnum, NStep, MaxLe);
			printf("end GMRES\n");
		}
		else if (solvers == 3)
		{
			//ACA


			printf("end ACA\n");
		}

		//tecplot
		for(i=0;i<=NStep;++i)
			ResultPlotNodeAverage_14VARS(m_DSE,m_PointList,m_ElePID,bdid,bd[i],i,PointNum,NodeNum,EleNum,DBEMOutputPath("TecValueFile").c_str(),amplitude);


		// ----输出选定节点的位移/面力随时间的变化信息-----
		//GetInfoFromPoint_UppersurfaceforStatic(m_NodeList, bd, m_EleNID, EleNum, eleflag, gap);
		//GetInfoFromPoint_AllforStatic(m_NodeList, bd, m_EleNID, EleNum);
		//GetInfoFromPoint_SingleColumn(m_NodeList, bd, NodeNum, 1, dt, dudt, C_1D, 1.0);
	}
	else
	{
		printf("不可接收的参数,程序退出...\n");

	}

	WriteValidationOutputs(m_DSE, m_PointList, m_ElePID, bdid, bd, EleNum, NStep, dt, solvers, MultiDomain);

	//---------------释放分配空间---------------------

	DSquareElement::StaticRelease(NodeNum);

	delete[] m_PointList;
	m_PointList = 0;

	for (i = 0; i < EleNum; ++i)
	{
		delete[] m_ElePID[i];
		m_ElePID[i] = 0;
		delete[] m_EleNID[i];
		m_EleNID[i] = 0;
	}
	delete[] bd;
	bd = 0;
	delete[] m_ElePID;
	m_ElePID = 0;
	delete[] m_EleNID;
	m_EleNID = 0;

	delete[] bdid;
	bdid = 0;

	delete[] m_DSE;
	m_DSE = 0;

	delete m_NodeList;
	m_NodeList = 0;
	
	endtime = clock();
	printf_s("总用时 %ld\n", endtime - starttime);
	fprintf_s(logfile, "求解用时 %ld\n", endtime - starttime);

	FreeBatchBoundaryFiles(BatchBoundaryFiles, BatchBoundaryCaseCount);
	fclose(logfile);

	return 1;

}

// 根据几何节点信息获取物理点信息
int GetNodeInfo(Point* m_pointlist, Point* &m_NodeList, long** m_ElePID, long** m_EleNID, long PointNum, long EleNum, long& NodeNum, DSquareElement* m_DSE)
{
	// 前提：单元已经几何初始化


	double gap = DSquareElement::m_gap;

	NodeNum = EleNum * 8;
	m_NodeList = new Point[NodeNum];

	#pragma omp parallel for num_threads(thread_num)
	for (long i = 0; i < EleNum; ++i)
	{
		m_EleNID[i][0] = i * 8 + 0;
		m_EleNID[i][1] = i * 8 + 1;
		m_EleNID[i][2] = i * 8 + 2;
		m_EleNID[i][3] = i * 8 + 3;
		m_EleNID[i][4] = i * 8 + 4;
		m_EleNID[i][5] = i * 8 + 5;
		m_EleNID[i][6] = i * 8 + 6;
		m_EleNID[i][7] = i * 8 + 7;

		m_DSE[i].GetFromLocal(-gap, -gap, m_NodeList[i * 8 + 0]);
		m_DSE[i].GetFromLocal(gap, -gap, m_NodeList[i * 8 + 1]);
		m_DSE[i].GetFromLocal(gap, gap, m_NodeList[i * 8 + 2]);
		m_DSE[i].GetFromLocal(-gap, gap, m_NodeList[i * 8 + 3]);
		m_DSE[i].GetFromLocal(0.0, -gap, m_NodeList[i * 8 + 4]);
		m_DSE[i].GetFromLocal(gap, 0.0, m_NodeList[i * 8 + 5]);
		m_DSE[i].GetFromLocal(0.0, gap, m_NodeList[i * 8 + 6]);
		m_DSE[i].GetFromLocal(-gap, 0.0, m_NodeList[i * 8 + 7]);
	}

	return 1;
}
// 画网格函数
int Plot(Point* &m_point, long** &m_ElePID, char* name, long PN, long EN, int flag, double amplitude)
{
	// flag
	// 1: 三角形
	// 2: 四边形线性
	// 3: 四边形二次

	int i, j;
	int TPointNum, TEleNum;
	double ** Tvalue;
	int** TEID;

	int eletype;

	switch (flag)
	{
	case 1:
	{
		eletype = 3;
		break;
	}
	case 2:
	{
		eletype = 4;
		break;
	}
	case 3:
	{
		eletype = 8;
		break;
	}
	}

	TPointNum = PN;
	TEleNum = EN;

	TEID = new int*[EN];
	for (i = 0; i < EN; ++i)
	{
		TEID[i] = new int[eletype];
		for (j = 0; j < eletype; ++j)
			TEID[i][j] = m_ElePID[i][j];
	}

	Tvalue = new double*[TPointNum];
	for (i = 0; i < TPointNum; ++i)
	{
		Tvalue[i] = new double[1];
		Tvalue[i][0] = 0.0;
	}
	for (i = 0; i < TEleNum; ++i)
	{
		for (j = 0; j < eletype; ++j)
			TEID[i][j] ++;
	}

	write_tecplot(TPointNum, TEleNum, eletype, 1, m_point, TEID, Tvalue, 1, 0, name, amplitude);

	for (i = 0; i < TEleNum; ++i)
		delete[] TEID[i];
	delete[] TEID;
	for (i = 0; i < TPointNum; ++i)
		delete[] Tvalue[i];
	delete[] Tvalue;

	return 1;
}



static void DBEMWriteTecplot14Row(FILE* out, const double coord[3], const double* values)
{
	for (int j = 0; j < 3; ++j)
		fprintf(out, "%.10g\t", coord[j]);
	for (int j = 0; j < 14; ++j)
		fprintf(out, "%.10g\t", values[j]);
	fprintf(out, "\n");
}

static void DBEMWriteTecplot14Zone(FILE* out,
	const char* title,
	long nodeCount,
	long eleCount,
	Point* nodes,
	int** elePid,
	double** values,
	Point* centerNodes,
	double** centerValues,
	double** displayDisp,
	double** centerDisplayDisp,
	int deformed)
{
	fprintf(out, "ZONE T=\"%s\", N=%ld, E=%ld, F=FEPOINT, ET=QUADRILATERAL\n", title, nodeCount + eleCount, 4 * eleCount);
	for (long i = 0; i < nodeCount; ++i)
	{
		double coord[3];
		for (int j = 0; j < 3; ++j)
			coord[j] = nodes[i].pt[j] + (deformed && displayDisp ? displayDisp[i][j] : 0.0);
		DBEMWriteTecplot14Row(out, coord, values[i]);
	}
	for (long i = 0; i < eleCount; ++i)
	{
		double coord[3];
		for (int j = 0; j < 3; ++j)
			coord[j] = centerNodes[i].pt[j] + (deformed && centerDisplayDisp ? centerDisplayDisp[i][j] : 0.0);
		DBEMWriteTecplot14Row(out, coord, centerValues[i]);
	}
	for (long i = 0; i < eleCount; ++i)
	{
		long centerId = nodeCount + i + 1;
		fprintf(out, "%d\t%d\t%ld\t%d\n", elePid[i][0], elePid[i][4], centerId, elePid[i][7]);
		fprintf(out, "%d\t%d\t%d\t%ld\n", elePid[i][4], elePid[i][1], elePid[i][5], centerId);
		fprintf(out, "%ld\t%d\t%d\t%d\n", centerId, elePid[i][5], elePid[i][2], elePid[i][6]);
		fprintf(out, "%d\t%ld\t%d\t%d\n", elePid[i][7], centerId, elePid[i][6], elePid[i][3]);
	}
}

static int DBEMWriteTecplot14WithDisplayDisp(long nodeCount,
	long eleCount,
	Point* nodes,
	int** elePid,
	double** values,
	double** displayDisp,
	char* filename)
{
	FILE* out = 0;
	if (fopen_s(&out, filename, "w") != 0 || !out)
		return 0;

	Point* centerNodes = new Point[eleCount];
	double** centerValues = new double* [eleCount];
	double** centerDisplayDisp = new double* [eleCount];
	for (long ele = 0; ele < eleCount; ++ele)
	{
		centerValues[ele] = new double[14];
		centerDisplayDisp[ele] = new double[3];
		for (int j = 0; j < 3; ++j)
		{
			centerNodes[ele].pt[j] = 0.0;
			centerDisplayDisp[ele][j] = 0.0;
		}
		for (int j = 0; j < 14; ++j)
			centerValues[ele][j] = 0.0;

		for (int localNode = 0; localNode < 8; ++localNode)
		{
			long node = elePid[ele][localNode] - 1;
			double weight = localNode < 4 ? -0.25 : 0.5;
			for (int j = 0; j < 3; ++j)
			{
				centerNodes[ele].pt[j] += weight * nodes[node].pt[j];
				centerDisplayDisp[ele][j] += weight * displayDisp[node][j];
			}
			for (int j = 0; j < 14; ++j)
				centerValues[ele][j] += weight * values[node][j];
		}
	}

	fprintf(out, "TITLE = \"BEM: 3D SURFACE ZONES\"\n");
	fprintf(out, "VARIABLES = \"X\", \"Y\", \"Z\"");
	for (int j = 0; j < 14; ++j)
		fprintf(out, ", \"var%d\"", j);
	fprintf(out, "\n");

	DBEMWriteTecplot14Zone(out, "original", nodeCount, eleCount, nodes, elePid, values, centerNodes, centerValues, displayDisp, centerDisplayDisp, 0);
	DBEMWriteTecplot14Zone(out, "deformed", nodeCount, eleCount, nodes, elePid, values, centerNodes, centerValues, displayDisp, centerDisplayDisp, 1);

	for (long ele = 0; ele < eleCount; ++ele)
	{
		delete[] centerValues[ele];
		delete[] centerDisplayDisp[ele];
	}
	delete[] centerValues;
	delete[] centerDisplayDisp;
	delete[] centerNodes;
	fclose(out);
	return 1;
}

static int WriteElementGeometryTecplot14VARS(DSquareElement* m_DSE,
	Point* m_PointList,
	long** m_ElePID,
	BoundaryValue& bd,
	long n,
	long PointNum,
	long NodeNum,
	long EleNum,
	const char* fp,
	double amplitude)
{
	if (!m_DSE || !m_PointList || !m_ElePID || !bd.lu.b || !bd.lt.b || NodeNum <= 0 || EleNum <= 0)
		return 0;

	char fpn[512];
	sprintf_s(fpn, "%s_%ld.dat", fp, n);

	long displayNodeCount = EleNum * 8;
	Point* tecNodes = new Point[displayNodeCount];
	int** tecElePid = new int* [EleNum];
	double** tecvalue = new double* [displayNodeCount];
	double** displayDisp = new double* [displayNodeCount];
	long* tecGeometryNode = new long[displayNodeCount];
	double* geometryDispSum = new double[PointNum * 3];
	long* geometryDispCount = new long[PointNum];

	for (long point = 0; point < PointNum; ++point)
	{
		geometryDispCount[point] = 0;
		geometryDispSum[3 * point + 0] = 0.0;
		geometryDispSum[3 * point + 1] = 0.0;
		geometryDispSum[3 * point + 2] = 0.0;
	}

	for (long node = 0; node < displayNodeCount; ++node)
	{
		tecNodes[node].pt[0] = 0.0;
		tecNodes[node].pt[1] = 0.0;
		tecNodes[node].pt[2] = 0.0;
		tecGeometryNode[node] = -1;
		tecvalue[node] = new double[14];
		displayDisp[node] = new double[3];
		for (int k = 0; k < 14; ++k)
			tecvalue[node][k] = 0.0;
		for (int k = 0; k < 3; ++k)
			displayDisp[node][k] = 0.0;
	}

	for (long ele = 0; ele < EleNum; ++ele)
	{
		tecElePid[ele] = new int[8];
		for (int localNode = 0; localNode < 8; ++localNode)
		{
			long displayNode = ele * 8 + localNode;
			long geometryNode = m_ElePID[ele][localNode];
			long physicalNode = m_DSE[ele].m_nodeID[localNode];
			tecElePid[ele][localNode] = (int)displayNode + 1;
			tecGeometryNode[displayNode] = geometryNode;
			if (geometryNode >= 0 && geometryNode < PointNum)
				tecNodes[displayNode] = m_PointList[geometryNode];
			if (physicalNode >= 0 && physicalNode < NodeNum)
			{
				tecvalue[displayNode][0] = bd.lu.b[3 * physicalNode + 0];
				tecvalue[displayNode][1] = bd.lu.b[3 * physicalNode + 1];
				tecvalue[displayNode][2] = bd.lu.b[3 * physicalNode + 2];
				if (geometryNode >= 0 && geometryNode < PointNum)
				{
					geometryDispSum[3 * geometryNode + 0] += tecvalue[displayNode][0];
					geometryDispSum[3 * geometryNode + 1] += tecvalue[displayNode][1];
					geometryDispSum[3 * geometryNode + 2] += tecvalue[displayNode][2];
					geometryDispCount[geometryNode]++;
				}
			}
		}
	}

	for (long node = 0; node < displayNodeCount; ++node)
	{
		long geometryNode = tecGeometryNode[node];
		for (int k = 0; k < 3; ++k)
		{
			displayDisp[node][k] = tecvalue[node][k];
			if (geometryNode >= 0 && geometryNode < PointNum && geometryDispCount[geometryNode] > 0)
				displayDisp[node][k] = geometryDispSum[3 * geometryNode + k] / (double)geometryDispCount[geometryNode];
		}
	}

	double stresses[3][3];
	double absdisp[8][3];
	double abstrac[8][3];
	for (long ele = 0; ele < EleNum; ++ele)
	{
		for (int localNode = 0; localNode < 8; ++localNode)
		{
			long physicalNode = m_DSE[ele].m_nodeID[localNode];
			if (physicalNode >= 0 && physicalNode < NodeNum)
			{
				absdisp[localNode][0] = bd.lu.b[3 * physicalNode + 0];
				absdisp[localNode][1] = bd.lu.b[3 * physicalNode + 1];
				absdisp[localNode][2] = bd.lu.b[3 * physicalNode + 2];
				abstrac[localNode][0] = bd.lt.b[3 * physicalNode + 0];
				abstrac[localNode][1] = bd.lt.b[3 * physicalNode + 1];
				abstrac[localNode][2] = bd.lt.b[3 * physicalNode + 2];
			}
			else
			{
				absdisp[localNode][0] = absdisp[localNode][1] = absdisp[localNode][2] = 0.0;
				abstrac[localNode][0] = abstrac[localNode][1] = abstrac[localNode][2] = 0.0;
			}
		}

		for (int localNode = 0; localNode < 8; ++localNode)
		{
			long displayNode = ele * 8 + localNode;
			long physicalNode = m_DSE[ele].m_nodeID[localNode];
			if (physicalNode < 0 || physicalNode >= NodeNum)
				continue;
			if (!m_DSE[ele].ObtainStress(DSquareElement::m_quadinfo.m_LocalCoord[localNode][0],
				DSquareElement::m_quadinfo.m_LocalCoord[localNode][1], stresses, absdisp, abstrac))
				DBEMZeroStress(stresses);

			tecvalue[displayNode][3] = stresses[0][0];
			tecvalue[displayNode][4] = stresses[1][1];
			tecvalue[displayNode][5] = stresses[2][2];
			tecvalue[displayNode][6] = stresses[0][1];
			tecvalue[displayNode][7] = stresses[0][2];
			tecvalue[displayNode][8] = stresses[1][2];

			if (DBEMStressComponentsFinite(stresses))
				m_DSE[ele].MainStress(stresses, tecvalue[displayNode][9], tecvalue[displayNode][10],
					tecvalue[displayNode][11], tecvalue[displayNode][12], tecvalue[displayNode][13]);
			else
			{
				tecvalue[displayNode][9] = 0.0;
				tecvalue[displayNode][10] = 0.0;
				tecvalue[displayNode][11] = 0.0;
				tecvalue[displayNode][12] = 0.0;
				tecvalue[displayNode][13] = 0.0;
			}
			DBEMCleanTecValues(tecvalue[displayNode], 14);
		}
	}

	(void)amplitude;
	if (!DBEMWriteTecplot14WithDisplayDisp(displayNodeCount, EleNum, tecNodes, tecElePid, tecvalue, displayDisp, fpn))
	{
		for (long ele = 0; ele < EleNum; ++ele)
			delete[] tecElePid[ele];
		delete[] tecElePid;
		for (long node = 0; node < displayNodeCount; ++node)
		{
			delete[] tecvalue[node];
			delete[] displayDisp[node];
		}
		delete[] tecvalue;
		delete[] displayDisp;
		delete[] tecGeometryNode;
		delete[] geometryDispSum;
		delete[] geometryDispCount;
		delete[] tecNodes;
		return 0;
	}

	for (long ele = 0; ele < EleNum; ++ele)
		delete[] tecElePid[ele];
	delete[] tecElePid;
	for (long node = 0; node < displayNodeCount; ++node)
	{
		delete[] tecvalue[node];
		delete[] displayDisp[node];
	}
	delete[] tecvalue;
	delete[] displayDisp;
	delete[] tecGeometryNode;
	delete[] geometryDispSum;
	delete[] geometryDispCount;
	delete[] tecNodes;

	return 1;
}

int ResultPlotNodeAverage_14VARS(DSquareElement* m_DSE, Point* m_PointList, long** m_ElePID, int* bdid, BoundaryValue& bd, long n, long PointNum, long NodeNum, long EleNum, const char* fp, double amplitude)
{
	(void)bdid;
	return WriteElementGeometryTecplot14VARS(m_DSE, m_PointList, m_ElePID, bd, n, PointNum, NodeNum, EleNum, fp, amplitude);
}

int TResultPlotNodeAverage_14VARS(DSquareElement* m_DSE, Point* m_PointList, long** m_ElePID, int* bdid, BoundaryValue& bd, long n, long PointNum, long NodeNum, long EleNum, const char* fp, double amplitude, long PointID, double& DispR)
{
	long i, j, NodeID;

	char fpn[512];
	sprintf_s(fpn, "%s_%ld.dat", fp, n);

	//Tecplot arrays
	Wvector DispX(NodeNum, 0);
	Wvector DispY(NodeNum, 0);
	Wvector DispZ(NodeNum, 0);

	//stress
	Wvector Stress11(NodeNum, 0);
	Wvector Stress22(NodeNum, 0);
	Wvector Stress33(NodeNum, 0);
	Wvector Stress12(NodeNum, 0);
	Wvector Stress13(NodeNum, 0);
	Wvector Stress23(NodeNum, 0);

	Wvector Stress1(NodeNum, 0);
	Wvector Stress2(NodeNum, 0);
	Wvector Stress3(NodeNum, 0);
	Wvector Mises(NodeNum, 0);
	Wvector Tresca(NodeNum, 0);

	//flag
	int* NodeFlag = new int[NodeNum];

	int * TecPointAdd = new int[PointNum];
	int **tecElePid = new int*[EleNum];
	for (i = 0; i < EleNum; ++i)
		tecElePid[i] = new int[8];
	double **tecvalue = new double*[PointNum];
	for (i = 0; i < PointNum; ++i)
		tecvalue[i] = new double[14];

	//Tecplot arrays
	for (i = 0; i < EleNum; ++i)
	{
		for (j = 0; j < 8; ++j)
			tecElePid[i][j] = m_ElePID[i][j] + 1;
	}
	for (i = 0; i < PointNum; ++i)
		TecPointAdd[i] = 0;

	// 根据求解向量和边界条件确定局部坐标值
	for (i = 0; i < NodeNum; ++i)
		NodeFlag[i] = 1;

	for (i = 0; i < EleNum; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			NodeID = m_DSE[i].m_nodeID[j];
			DispX[NodeID + 1] = bd.lu.b[3 * NodeID + 0];
			DispY[NodeID + 1] = bd.lu.b[3 * NodeID + 1];
			DispZ[NodeID + 1] = bd.lu.b[3 * NodeID + 2];
		}
	}

	//应力计算
	double stresses[3][3];
	double absdisp[8][3];
	double abstrac[8][3];
	//	double ts1, ts2, ts3, es1, es2;

	for (i = 0; i < EleNum; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			NodeID = m_DSE[i].m_nodeID[j];
			absdisp[j][0] = bd.lu.b[3 * NodeID + 0];
			absdisp[j][1] = bd.lu.b[3 * NodeID + 1];
			absdisp[j][2] = bd.lu.b[3 * NodeID + 2];
			abstrac[j][0] = bd.lt.b[3 * NodeID + 0];
			abstrac[j][1] = bd.lt.b[3 * NodeID + 1];
			abstrac[j][2] = bd.lt.b[3 * NodeID + 2];
		}

		for (j = 0; j < 8; ++j)
		{
			NodeID = m_DSE[i].m_nodeID[j];
			if (!m_DSE[i].ObtainStress(DSquareElement::m_quadinfo.m_LocalCoord[j][0], DSquareElement::m_quadinfo.m_LocalCoord[j][1], stresses, absdisp, abstrac))
				DBEMZeroStress(stresses);

			Stress11[NodeID + 1] = stresses[0][0];
			Stress22[NodeID + 1] = stresses[1][1];
			Stress33[NodeID + 1] = stresses[2][2];
			Stress12[NodeID + 1] = stresses[0][1];
			Stress13[NodeID + 1] = stresses[0][2];
			Stress23[NodeID + 1] = stresses[1][2];
		}
	}

	for (i = 0; i < PointNum; ++i)
	{
		tecvalue[i][0] = 0.0;
		tecvalue[i][1] = 0.0;
		tecvalue[i][2] = 0.0;
		tecvalue[i][3] = 0.0;
		tecvalue[i][4] = 0.0;
		tecvalue[i][5] = 0.0;
		tecvalue[i][6] = 0.0;
		tecvalue[i][7] = 0.0;
		tecvalue[i][8] = 0.0;
		tecvalue[i][9] = 0.0;
		tecvalue[i][10] = 0.0;
		tecvalue[i][11] = 0.0;
		tecvalue[i][12] = 0.0;
		tecvalue[i][13] = 0.0;
	}

	Wvector tecinput(8), tecoutput(8);
	Wvector* nodalValues[9] = { &DispX, &DispY, &DispZ, &Stress11, &Stress22, &Stress33, &Stress12, &Stress13, &Stress23 };

	for (i = 0; i < EleNum; ++i)
	{
		for (int valueId = 0; valueId < 9; ++valueId)
		{
			for (j = 0; j < 8; ++j)
				tecinput[j + 1] = (*nodalValues[valueId])[m_DSE[i].m_nodeID[j] + 1];
			m_DSE[i].Approximate(tecinput.b, tecoutput.b);
			for (j = 0; j < 8; ++j)
				tecvalue[tecElePid[i][j] - 1][valueId] += tecoutput.b[j];
		}


		// 节点所属的单元数量
		TecPointAdd[tecElePid[i][0] - 1]++;
		TecPointAdd[tecElePid[i][1] - 1]++;
		TecPointAdd[tecElePid[i][2] - 1]++;
		TecPointAdd[tecElePid[i][3] - 1]++;
		TecPointAdd[tecElePid[i][4] - 1]++;
		TecPointAdd[tecElePid[i][5] - 1]++;
		TecPointAdd[tecElePid[i][6] - 1]++;
		TecPointAdd[tecElePid[i][7] - 1]++;

	}

	for (i = 0; i < PointNum; ++i)
	{
		if (TecPointAdd[i] > 0)
		{
			tecvalue[i][0] /= TecPointAdd[i];
			tecvalue[i][1] /= TecPointAdd[i];
			tecvalue[i][2] /= TecPointAdd[i];
			tecvalue[i][3] /= TecPointAdd[i];
			tecvalue[i][4] /= TecPointAdd[i];
			tecvalue[i][5] /= TecPointAdd[i];
			tecvalue[i][6] /= TecPointAdd[i];
			tecvalue[i][7] /= TecPointAdd[i];
			tecvalue[i][8] /= TecPointAdd[i];
		}
		DBEMCleanTecValues(tecvalue[i], 14);
	}

	for (i = 0; i < PointNum; ++i)
	{
		stresses[0][0] = tecvalue[i][3];
		stresses[1][1] = tecvalue[i][4];
		stresses[2][2] = tecvalue[i][5];
		stresses[0][1] = tecvalue[i][6];
		stresses[0][2] = tecvalue[i][7];
		stresses[1][2] = tecvalue[i][8];
		stresses[1][0] = stresses[0][1];
		stresses[2][0] = stresses[0][2];
		stresses[2][1] = stresses[1][2];

		if (DBEMStressComponentsFinite(stresses))
			m_DSE[1].MainStress(stresses, tecvalue[i][9], tecvalue[i][10], tecvalue[i][11], tecvalue[i][12], tecvalue[i][13]);
		else
		{
			tecvalue[i][9] = 0.0;
			tecvalue[i][10] = 0.0;
			tecvalue[i][11] = 0.0;
			tecvalue[i][12] = 0.0;
			tecvalue[i][13] = 0.0;
		}
		DBEMCleanTecValues(tecvalue[i], 14);
	}



	// 判断中间节点所属的单元信息

	int* InterPointCount = new int[PointNum];           // 记录中间节点所属的单元数量
	long(*InterElementID)[4] = new long[PointNum][4];   // 记录中间节点所属单元的ID，最多属于4个单元
	int(*InterElementPos)[4] = new int[PointNum][4];    // 记录中间节点在所属单元的位置

	for (i = 0; i < PointNum; ++i)
		InterPointCount[i] = 0;

	for (i = 0; i < EleNum; ++i)
	{
		InterElementID[m_ElePID[i][4]][InterPointCount[m_ElePID[i][4]]] = i;
		InterElementID[m_ElePID[i][5]][InterPointCount[m_ElePID[i][5]]] = i;
		InterElementID[m_ElePID[i][6]][InterPointCount[m_ElePID[i][6]]] = i;
		InterElementID[m_ElePID[i][7]][InterPointCount[m_ElePID[i][7]]] = i;

		InterElementPos[m_ElePID[i][4]][InterPointCount[m_ElePID[i][4]]] = 4;
		InterElementPos[m_ElePID[i][5]][InterPointCount[m_ElePID[i][5]]] = 5;
		InterElementPos[m_ElePID[i][6]][InterPointCount[m_ElePID[i][6]]] = 6;
		InterElementPos[m_ElePID[i][7]][InterPointCount[m_ElePID[i][7]]] = 7;

		InterPointCount[m_ElePID[i][4]] ++;
		InterPointCount[m_ElePID[i][5]] ++;
		InterPointCount[m_ElePID[i][6]] ++;
		InterPointCount[m_ElePID[i][7]] ++;
	}

	write_tecplot(PointNum, EleNum, 8, 14, m_PointList, tecElePid, tecvalue, 1, 1, fpn, amplitude);
	//Tecplot arrays
	//tecplot

	//---------------release space---------------------------



	DispR = tecvalue[PointID][0];




	//tecplot arrays
	delete[] TecPointAdd;
	TecPointAdd = 0;
	for (i = 0; i < EleNum; ++i)
	{
		delete[] tecElePid[i];
		tecElePid[i] = 0;
	}
	delete[] tecElePid;
	tecElePid = 0;
	for (i = 0; i < PointNum; ++i)
	{
		delete[] tecvalue[i];
		tecvalue[i] = 0;
	}
	delete[] tecvalue;
	tecvalue = 0;
	//tecplot arrays

	delete[] InterPointCount;
	InterPointCount = 0;
	delete[] InterElementID;
	InterElementID = 0;
	delete[] InterElementPos;
	InterElementPos = 0;

	delete[] NodeFlag;

	return 1;
}

int GetInfoFromPoint_SingleColumn(Point* m_NodeList, BoundaryValue* bd, long NodeNum, long NStep, double dt, double dudt, double C_1D, double ColumnLength, long ** m_EleNID)
{
	long i;
	Point DispPoint, TracPoint;
	double mindisp, mintrac;
	long MinDispNodeID, MinTracNodeID;

	// 解析解定义
	int FlagTurn, FlagOne;
	double TimeSum;
	double TimeLength = ColumnLength / C_1D;

	double MaxDisp = dudt * 2.0 * TimeLength;

	double* Ana_Disp = new double[10 * NStep + 10];
	double* Ana_Trac = new double[10 * NStep + 10];

	TimeSum = 0.0;
	FlagTurn = 1;

	// 位移解析解

	for (i = 0; i < 10 * NStep + 10; ++i)
	{
		if (TimeSum > 2.0 * TimeLength)
		{
			TimeSum -= 2.0 * TimeLength;
			if (FlagTurn == 1)
				FlagTurn = -1;
			else
				FlagTurn = 1;
		}

		if (FlagTurn == 1)
			Ana_Disp[i] = -(dudt * TimeSum);
		else
			Ana_Disp[i] = -(MaxDisp - dudt * TimeSum);

		TimeSum += dt / 10.0;
	}

	// 面力解析解

	TimeSum = 0.0;
	FlagTurn = 1;
	FlagOne = 0;

	for (i = 0; i < 10 * NStep + 10; ++i)
	{
		if (FlagOne == 0)
		{
			if (TimeSum > TimeLength)
			{
				TimeSum -= TimeLength;
				if (FlagTurn == 1)
					FlagTurn = -1;
				else
					FlagTurn = 1;

				FlagOne = 1;
			}
		}
		else
		{
			if (TimeSum > 2 * TimeLength)
			{
				TimeSum -= 2 * TimeLength;
				if (FlagTurn == 1)
					FlagTurn = -1;
				else
					FlagTurn = 1;

				FlagOne = 1;
			}
		}

		if (FlagTurn == 1)
			Ana_Trac[i] = 0.0;
		else
			Ana_Trac[i] = 2.0;

		TimeSum += dt / 10.0;
	}

	// 获取距离提取点最近的点的ID

	mindisp = 100000.0;
	mintrac = 100000.0;

	DispPoint.pt[0] = 0.0;
	DispPoint.pt[1] = 0.0;
	DispPoint.pt[2] = ColumnLength;//ColumnLength;
	TracPoint.pt[0] = 0.0;
	TracPoint.pt[1] = 0.0;
	TracPoint.pt[2] = 0.0;

	for (i = 0; i < NodeNum; ++i)
	{
		if (mindisp > Dist(DispPoint, m_NodeList[i]))
		{
			mindisp = Dist(DispPoint, m_NodeList[i]);
			MinDispNodeID = i;
		}

		if (mintrac > Dist(TracPoint, m_NodeList[i]))
		{
			mintrac = Dist(TracPoint, m_NodeList[i]);
			MinTracNodeID = i;
		}
	}

	//FILE* fp;
	//char nameofstep[10];
	//int itemp,jtemp;
	//for (i = 0; i <= NStep; i=i+1)
	//{
	//	_itoa_s(i,nameofstep,10);
	//	fopen_s(&fp, nameofstep, "w");
	//	for (itemp = 0; itemp < 648; itemp++)
	//	{
	//		fprintf(fp, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", m_NodeList[itemp].pt[0], m_NodeList[itemp].pt[1], m_NodeList[itemp].pt[2], bd[i].lu.b[3*itemp + 0], bd[i].lu.b[3 * itemp + 1], bd[i].lu.b[3 * itemp + 2], bd[i].lt.b[3 * itemp + 0], bd[i].lt.b[3 * itemp + 1], bd[i].lt.b[3 * itemp + 2]);
	//	}
	//	fclose(fp);
	////	fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e\n", i * dt, bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2]/*[3 * MinDispNodeID + 2]*/, bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);
	//}

	
	long PointID, tempNodeID;
	double outputcheck=0.0;
	double res[8];
	FILE* fp2;
	std::string dispTracPath2 = DBEMOutputPath("Disp-Trac-Versus-Time-2.dat");
	fopen_s(&fp2, dispTracPath2.c_str(), "w");
	FILE* fp;
	std::string dispTracPath = DBEMOutputPath("Disp-Trac-Versus-Time.dat");
	fopen_s(&fp, dispTracPath.c_str(), "w");
	for (i = 0; i <= NStep; ++i)
	{
		//输出中心点位移
		//计算形函数
		DSquareElement::m_quadinfo.N(res, 1.0, 1.0);//局部坐标
		//计算位移
		for (PointID = 0; PointID < 8; ++PointID)
		{
			tempNodeID = m_EleNID[0][PointID];
			outputcheck += bd[i].lu.b[3 * tempNodeID + 2] * res[PointID];
		}
		fprintf_s(fp2, "%20.10e   %20.10e\n", i * dt, outputcheck);

		fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e\n",
			i * dt, bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2],
			bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);
	}
	fclose(fp);
	fclose(fp2);
	

	//FILE *fp;
	//fopen_s(&fp,"Ana-Disp-Trac-Versus-Time.dat", "w");
	//for (i = 0; i < 10 * NStep + 10; ++i)
	//{
	//	fprintf(fp, "%20.10e   %20.10e   %20.10e\n", i * dt / 10.0, Ana_Disp[i], Ana_Trac[i]);
	//}

	delete[] Ana_Disp;
	Ana_Disp = 0;
	delete[] Ana_Trac;
	Ana_Trac = 0;

	return 1;
}

double GetEleSize(Point* m_PointList, long** m_ElePID, long PointNum, long EleNum)
{
	long i;
	double avesize = 0.0;

	for (i = 0; i < EleNum; ++i)
	{
		avesize += Dist(m_PointList[m_ElePID[i][0]], m_PointList[m_ElePID[i][1]]);
	}

	avesize /= EleNum;

	return avesize;
}

double GetRadiusDispFromSphere(double rou, double a, double v, double t, double c1)
{
	double gamma, s, f, dfdt;

	gamma = (1.0 - 2.0 * v) / (1.0 - v) * c1 / a;
	s = sqrt(1.0 / (1.0 - 2.0 * v));
	f = -a * a / (2.0 * rou * gamma * c1) * (1.0 - exp(-gamma * t) * (cos(gamma*s*t) + 1.0 / s * sin(gamma*s*t)));
	dfdt = -a * a / (2.0 * rou * c1) * exp(-gamma * t) * (1.0 / s + s) * sin(gamma*s*t);

	return (-1.0 / c1 / a * dfdt - f / a / a);
}

long GetPoint100ID(Point* m_pointlist, long PointNum)
{
	long i, ID;
	double dist = 1000000.0;
	Point ID100;
	ID100.pt[0] = 1.0;
	ID100.pt[1] = 0.0;
	ID100.pt[2] = 0.0;

	for (i = 0; i < PointNum; ++i)
	{
		if (dist > Dist(ID100, m_pointlist[i]))
		{
			dist = Dist(ID100, m_pointlist[i]);
			ID = i;
		}
	}

	return ID;
}

int GetInfoFromPoint_Uppersurface(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag, long NStep)
{
	long i;
	FILE* fp;
	char nameofstep[10];
	int itemp, jtemp;
	int tempNodeID;
	for (i = 0; i <= NStep; i = i + 1)
	{

		_itoa_s(i, nameofstep, 10);
		fopen_s(&fp, nameofstep, "w");
		if (i < NStep)//最后一步有问题，暂复制前一步结果
		{
			for (itemp = 0; itemp < EleNum; itemp++)
			{

				if (eleflag[itemp] == 'T')
				{

					for (jtemp = 0; jtemp < 8; jtemp++)
					{
						tempNodeID = m_EleNID[itemp][jtemp];
						fprintf(fp, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", m_NodeList[tempNodeID].pt[0], m_NodeList[tempNodeID].pt[1], m_NodeList[tempNodeID].pt[2], bd[i].lu.b[3 * tempNodeID + 0], bd[i].lu.b[3 * tempNodeID + 1], bd[i].lu.b[3 * tempNodeID + 2], bd[i].lt.b[3 * tempNodeID + 0], bd[i].lt.b[3 * tempNodeID + 1], bd[i].lt.b[3 * tempNodeID + 2]);

					}

				}
			}
		}
		else {
			for (itemp = 0; itemp < EleNum; itemp++)
			{

				if (eleflag[itemp] == 'T')
				{

					for (jtemp = 0; jtemp < 8; jtemp++)
					{
						tempNodeID = m_EleNID[itemp][jtemp];
						fprintf(fp, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", m_NodeList[tempNodeID].pt[0], m_NodeList[tempNodeID].pt[1], m_NodeList[tempNodeID].pt[2], bd[i-1].lu.b[3 * tempNodeID + 0], bd[i-1].lu.b[3 * tempNodeID + 1], bd[i-1].lu.b[3 * tempNodeID + 2], bd[i-1].lt.b[3 * tempNodeID + 0], bd[i-1].lt.b[3 * tempNodeID + 1], bd[i-1].lt.b[3 * tempNodeID + 2]);

					}

				}
			}
		}
		
		fclose(fp);



	}
	return 1;
}

int GetInfoFromPoint_All(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag, long NStep)
{
	long i;
	FILE* fp;
	char nameofstep[10],openfile[30],tempchar[20]=".//result//";
	int itemp, jtemp;
	int tempNodeID;
	char* CurrPath = _getcwd(NULL, 0);//返回当前路径
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//返回当前文件夹名称
	Path = StrCurrPath.substr(0, found + 1);//返回上级路径
	string OutFile = "output\\";
	string TempPath= Path + OutFile;
	string TempNum;
	string ForMatTxt = ".txt";//文件后缀    Git test
	const char* OutPath;



	for (i = 0; i <= NStep; i = i + 1)
	{
		FILE* fp;
		TempNum = to_string(i);

		Path = TempPath + TempNum + ForMatTxt;
		OutPath = Path.c_str();
		fopen_s(&fp, OutPath, "w");
		//strcpy_s(openfile, tempchar);
		//_itoa_s(i, nameofstep, 10);
		//strcat_s(openfile, nameofstep);
		//fopen_s(&fp, openfile, "w");

		if (i <= NStep)//最后一步有问题，暂复制前一步结果
		{
			for (itemp = 0; itemp < EleNum; itemp++)
			{

				for (jtemp = 0; jtemp < 8; jtemp++)
				{
					tempNodeID = m_EleNID[itemp][jtemp];
					fprintf(fp, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", m_NodeList[tempNodeID].pt[0], m_NodeList[tempNodeID].pt[1], m_NodeList[tempNodeID].pt[2], bd[i].lu.b[3 * tempNodeID + 0], bd[i].lu.b[3 * tempNodeID + 1], bd[i].lu.b[3 * tempNodeID + 2], bd[i].lt.b[3 * tempNodeID + 0], bd[i].lt.b[3 * tempNodeID + 1], bd[i].lt.b[3 * tempNodeID + 2]);

				}
			}
		}
		
		fclose(fp);

	}
	return 1;
}

int GetInfoFromPoint_UppersurfaceforStatic(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag)
{
	FILE* fp;
	int itemp, jtemp;
	int tempNodeID;


	fopen_s(&fp, "static_result.txt", "w");
	for (itemp = 0; itemp < EleNum; itemp++)
	{

		if (eleflag[itemp] == 'T')
		{

			for (jtemp = 0; jtemp < 8; jtemp++)
			{
				tempNodeID = m_EleNID[itemp][jtemp];
				fprintf(fp, "%d %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", itemp, m_NodeList[tempNodeID].pt[0], m_NodeList[tempNodeID].pt[1], m_NodeList[tempNodeID].pt[2], bd[1].lu.b[3 * tempNodeID + 0], bd[1].lu.b[3 * tempNodeID + 1], bd[1].lu.b[3 * tempNodeID + 2], bd[1].lt.b[3 * tempNodeID + 0], bd[1].lt.b[3 * tempNodeID + 1], bd[1].lt.b[3 * tempNodeID + 2]);

			}

		}
	}
	fclose(fp);
	return 1;
}

int GetInfoFromPoint_UppersurfaceforStatic(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag, double gap)
{
	FILE* fp;
	int itemp, jtemp, PointID;
	int tempNodeID;
	double res[8];
	double outputcheck = 0.0;

	fopen_s(&fp, "static_result.txt", "w");
	fprintf_s(fp, "gap = %lf\n\n", gap);
	fprintf_s(fp, "z = %lf\n\n", DSquareElement::InfCenter[2]);
	//输出中心点位移
	//计算形函数
	DSquareElement::m_quadinfo.N(res, 1.0, 1.0);//局部坐标
	//计算位移
	for (PointID = 0; PointID < 8; ++PointID)
	{
		tempNodeID = m_EleNID[5][PointID];
		outputcheck += bd[1].lu.b[3 * tempNodeID + 2] * res[PointID];

		fprintf_s(fp, "N = %20.10e\n\n", res[PointID]);
	}
	fprintf_s(fp, "center = %20.10e\n\n", outputcheck);
	//输出角点位移
	DSquareElement::m_quadinfo.N(res, -1.0, -1.0);//局部坐标
	//计算位移
	outputcheck = 0.0;
	for (PointID = 0; PointID < 8; ++PointID)
	{
		tempNodeID = m_EleNID[0][PointID];
		outputcheck += bd[1].lu.b[3 * tempNodeID + 2] * res[PointID];

		fprintf_s(fp, "N = %20.10e\n\n", res[PointID]);
	}
	fprintf_s(fp, "corner = %20.10e\n\n", outputcheck);

	for (itemp = 0; itemp < EleNum; itemp++)
	{

		if (eleflag[itemp] == 'T')
		{

			for (jtemp = 0; jtemp < 8; jtemp++)
			{
				tempNodeID = m_EleNID[itemp][jtemp];
				fprintf(fp, "%d %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", itemp, m_NodeList[tempNodeID].pt[0], m_NodeList[tempNodeID].pt[1], m_NodeList[tempNodeID].pt[2], bd[1].lu.b[3 * tempNodeID + 0], bd[1].lu.b[3 * tempNodeID + 1], bd[1].lu.b[3 * tempNodeID + 2], bd[1].lt.b[3 * tempNodeID + 0], bd[1].lt.b[3 * tempNodeID + 1], bd[1].lt.b[3 * tempNodeID + 2]);

			}

		}
	}
	fclose(fp);
	return 1;
}

int GetInfoFromPoint_AllforStatic(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum)
{
	FILE* fp;
	int itemp, jtemp;
	int tempNodeID;


	fopen_s(&fp, "static_result_AllPoint.txt", "w");
	for (itemp = 0; itemp < EleNum; itemp++)
	{
		for (jtemp = 0; jtemp < 8; jtemp++)
		{
			tempNodeID = m_EleNID[itemp][jtemp];
			fprintf(fp, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
				m_NodeList[tempNodeID].pt[0], m_NodeList[tempNodeID].pt[1], m_NodeList[tempNodeID].pt[2],
				bd[1].lu.b[3 * tempNodeID + 0], bd[1].lu.b[3 * tempNodeID + 1], bd[1].lu.b[3 * tempNodeID + 2],
				bd[1].lt.b[3 * tempNodeID + 0], bd[1].lt.b[3 * tempNodeID + 1], bd[1].lt.b[3 * tempNodeID + 2]);

		}
	}
	fclose(fp);
	return 1;
}

//int GetInfoFromPoint_onenode(Point* m_NodeList, BoundaryValue* bd, long NodeNum, long NStep, double dt, long ** m_EleNID)
//{
//	long i;
//	Point DispPoint, TracPoint;
//	double mindisp, mintrac;
//	long MinDispNodeID, MinTracNodeID;
//
//	// 获取距离提取点最近的点的ID
//
//	mindisp = 100000.0;
//	mintrac = 100000.0;
//
//	DispPoint.pt[0] = 0;
//	DispPoint.pt[1] = 0;
//	DispPoint.pt[2] = 0;
//	TracPoint.pt[0] = 0;
//	TracPoint.pt[1] =0;
//	TracPoint.pt[2] = 0;
//
//	for (i = 0; i < NodeNum; ++i)
//	{
//		if (mindisp > Dist(DispPoint, m_NodeList[i]))
//		{
//			mindisp = Dist(DispPoint, m_NodeList[i]);
//			MinDispNodeID = i;
//		}
//
//		if (mintrac > Dist(TracPoint, m_NodeList[i]))
//		{
//			mintrac = Dist(TracPoint, m_NodeList[i]);
//			MinTracNodeID = i;
//		}
//	}
//
//	//FILE* fp;
//	//char nameofstep[10];
//	//int itemp,jtemp;
//	//for (i = 0; i <= NStep; i=i+1)
//	//{
//	//	_itoa_s(i,nameofstep,10);
//	//	fopen_s(&fp, nameofstep, "w");
//	//	for (itemp = 0; itemp < 648; itemp++)
//	//	{
//	//		fprintf(fp, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", m_NodeList[itemp].pt[0], m_NodeList[itemp].pt[1], m_NodeList[itemp].pt[2], bd[i].lu.b[3*itemp + 0], bd[i].lu.b[3 * itemp + 1], bd[i].lu.b[3 * itemp + 2], bd[i].lt.b[3 * itemp + 0], bd[i].lt.b[3 * itemp + 1], bd[i].lt.b[3 * itemp + 2]);
//	//	}
//	//	fclose(fp);
//	////	fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e\n", i * dt, bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2]/*[3 * MinDispNodeID + 2]*/, bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);
//	//}
//
//	double x, y, z,ux,uy,uz,tx,ty,tz,u,t;
//	
//	long PointID, tempNodeID;
//	double outputcheck = 0.0;
//	double res[8];
//	FILE* fp2;
//	fopen_s(&fp2, "Disp-Trac-Versus-Time-onenode.dat", "w");
//	FILE* fp;
//	fopen_s(&fp, "Disp-Trac-Versus-Time.dat", "w");
//	for (i = 0; i <= NStep; ++i)
//	{
//		x = 0.0;
//		y = 0.0;
//		z = 0.0;
//		ux = uy = uz = 0;
//		tx = ty = tz = 0;
//		//输出中心点位移
//		//计算形函数
//		DSquareElement::m_quadinfo.N(res, 0.0, 0.0);//局部坐标
//		double ttm = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			ttm += res[j];
//			printf("%lf\n", res[j]);
//		}
//		printf("%lf!!!\n", ttm);
//		ttm = 0;
//		//计算位移
//		for (PointID = 0; PointID < 8; ++PointID)
//		{
//			tempNodeID = m_EleNID[0][PointID];//取单元此处修改
//			
//			x += m_NodeList[tempNodeID].pt[0] *res[PointID];
//			y += m_NodeList[tempNodeID].pt[1] *res[PointID];
//			z += m_NodeList[tempNodeID].pt[2] *res[PointID];
//			ux += bd[i].lu.b[3 * tempNodeID + 0] *res[PointID];
//			uy += bd[i].lu.b[3 * tempNodeID + 1] *res[PointID];
//			uz += bd[i].lu.b[3 * tempNodeID + 2] *res[PointID];
//			tx += bd[i].lt.b[3 * tempNodeID + 0]*res[PointID];
//			ty += bd[i].lt.b[3 * tempNodeID + 1] *res[PointID];
//			tz += bd[i].lt.b[3 * tempNodeID + 2] *res[PointID];
//			
////			printf("%.16lf %.16lf %.16lf\n", sqrt(bd[i].lu.b[3 * tempNodeID + 0]* bd[i].lu.b[3 * tempNodeID + 0]+bd[i].lu.b[3 * tempNodeID + 1]* bd[i].lu.b[3 * tempNodeID + 1]+ bd[i].lu.b[3 * tempNodeID + 2]* bd[i].lu.b[3 * tempNodeID + 2]),res[PointID], sqrt(bd[i].lu.b[3 * tempNodeID + 0] * bd[i].lu.b[3 * tempNodeID + 0] + bd[i].lu.b[3 * tempNodeID + 1] * bd[i].lu.b[3 * tempNodeID + 1] + bd[i].lu.b[3 * tempNodeID + 2] * bd[i].lu.b[3 * tempNodeID + 2])*res[PointID]);
//			
//		}
//		u = sqrt(ux*ux + uy * uy + uz * uz);
//		t = sqrt(tx*tx + ty * ty + tz * tz);
////		printf("%20.10e  %20.10e  %20.10e %20.10e  %20.10e %20.10e\n", i * dt, x, y, z, u, t);
//		//fprintf_s(fp2, "%20.10e  %20.10e  %20.10e %20.10e  %20.10e  %20.10e %20.10e  %20.10e  %20.10e %20.10e\n", i * dt, x,y,z,ux,uy,uz,tx,ty,tz);
//		fprintf_s(fp2, "%20.10e  %20.10e  %20.10e %20.10e  %20.10e %20.10e\n", i * dt, x, y, z, u,t);
//
//		//fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e\n",
//		//	i * dt, m_NodeList[MinDispNodeID].pt[0], m_NodeList[MinDispNodeID].pt[1], m_NodeList[MinDispNodeID].pt[2], sqrt(bd[i].lu.b[3 * MinDispNodeID + 0]* bd[i].lu.b[3 * MinDispNodeID + 0]+ bd[i].lu.b[3 * MinDispNodeID + 1] * bd[i].lu.b[3 * MinDispNodeID + 1] + bd[i].lu.b[3 * MinDispNodeID + 2] * bd[i].lu.b[3 * MinDispNodeID + 2]),
//		//	sqrt(bd[i].lt.b[3 * MinDispNodeID + 0] * bd[i].lt.b[3 * MinDispNodeID + 0] + bd[i].lt.b[3 * MinDispNodeID + 1] * bd[i].lt.b[3 * MinDispNodeID + 1] + bd[i].lt.b[3 * MinDispNodeID + 2] * bd[i].lt.b[3 * MinDispNodeID + 2]));
//		fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
//			i * dt, m_NodeList[MinDispNodeID].pt[0], m_NodeList[MinDispNodeID].pt[1], m_NodeList[MinDispNodeID].pt[2], bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2],
//			bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);
//		
//	}
//	fclose(fp);
//	fclose(fp2);
//
//	return 1;
//}

int GetInfoFromPoint_onenode(Point* m_NodeList, BoundaryValue* bd, long NodeNum, long NStep, double dt, long ** m_EleNID)
{
	long i;
	Point DispPoint, TracPoint;
	double mindisp, mintrac;
	long MinDispNodeID, MinTracNodeID;
	double m_X, m_Y, m_Z;

	m_X = 0;
	m_Y = 0;
	m_Z = 0;
	// 获取距离提取点最近的点的ID

	mindisp = 100000.0;
	mintrac = 100000.0;

	DispPoint.pt[0] = 0.001+ m_X;
	DispPoint.pt[1] = 0.001 + m_Y;
	DispPoint.pt[2] = 0+m_Z;
	TracPoint.pt[0] = 0.001 + m_X;
	TracPoint.pt[1] = 0.001 + m_Y;
	TracPoint.pt[2] = 0+m_Z;

	for (i = 0; i < NodeNum; ++i)
	{
		if (mindisp > Dist(DispPoint, m_NodeList[i]))
		{
			mindisp = Dist(DispPoint, m_NodeList[i]);
			MinDispNodeID = i;
		}

		if (mintrac > Dist(TracPoint, m_NodeList[i]))
		{
			mintrac = Dist(TracPoint, m_NodeList[i]);
			MinTracNodeID = i;
		}
	}
	char* CurrPath = _getcwd(NULL, 0);//返回当前路径
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//返回当前文件夹名称
	Path = StrCurrPath.substr(0, found + 1);//返回上级路径
	string OutFile = "output\\";
	string TempPath = Path + OutFile;
	string TempNum = to_string(i);
	string ForMatTxt = ".txt";
	//Path = TempPath + TempNum + ForMatTxt;
	//const char* OutPath = Path.c_str();

	MinDispNodeID = 38;//86;//指定节点

	FILE* fp;
	string Temp = "Disp-Trac-Versus-Time_1.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();
	fopen_s(&fp, OutPath, "w");
	for (i = 0; i <= NStep; ++i)
	{
			fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
			i * dt, m_NodeList[MinDispNodeID].pt[0], m_NodeList[MinDispNodeID].pt[1], m_NodeList[MinDispNodeID].pt[2], bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2],
			bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);

	}
	fclose(fp);

	// 获取距离提取点最近的点的ID

	mindisp = 100000.0;
	mintrac = 100000.0;

	DispPoint.pt[0] = 0.001 + m_X;
	DispPoint.pt[1] = -0.001 + m_Y;
	DispPoint.pt[2] = 0+m_Z;
	TracPoint.pt[0] = 0.001 + m_X;
	TracPoint.pt[1] = -0.001 + m_Y;
	TracPoint.pt[2] = 0+ m_Z;

	for (i = 0; i < NodeNum; ++i)
	{
		if (mindisp > Dist(DispPoint, m_NodeList[i]))
		{
			mindisp = Dist(DispPoint, m_NodeList[i]);
			MinDispNodeID = i;
		}

		if (mintrac > Dist(TracPoint, m_NodeList[i]))
		{
			mintrac = Dist(TracPoint, m_NodeList[i]);
			MinTracNodeID = i;
		}
	}

	Temp = "Disp-Trac-Versus-Time_2.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();
	fopen_s(&fp, OutPath, "w");
	for (i = 0; i <= NStep; ++i)
	{
		fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
			i * dt, m_NodeList[MinDispNodeID].pt[0], m_NodeList[MinDispNodeID].pt[1], m_NodeList[MinDispNodeID].pt[2], bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2],
			bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);

	}
	fclose(fp);
	// 获取距离提取点最近的点的ID

	mindisp = 100000.0;
	mintrac = 100000.0;

	DispPoint.pt[0] = -0.001 + m_X;
	DispPoint.pt[1] = -0.001 + m_Y;
	DispPoint.pt[2] = 0+m_Z;
	TracPoint.pt[0] = -0.001 + m_X;
	TracPoint.pt[1] = -0.001 + m_Y;
	TracPoint.pt[2] = 0+ m_Z;

	for (i = 0; i < NodeNum; ++i)
	{
		if (mindisp > Dist(DispPoint, m_NodeList[i]))
		{
			mindisp = Dist(DispPoint, m_NodeList[i]);
			MinDispNodeID = i;
		}

		if (mintrac > Dist(TracPoint, m_NodeList[i]))
		{
			mintrac = Dist(TracPoint, m_NodeList[i]);
			MinTracNodeID = i;
		}
	}
	Temp = "Disp-Trac-Versus-Time_3.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();

	fopen_s(&fp, OutPath, "w");
	for (i = 0; i <= NStep; ++i)
	{
		fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
			i * dt, m_NodeList[MinDispNodeID].pt[0], m_NodeList[MinDispNodeID].pt[1], m_NodeList[MinDispNodeID].pt[2], bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2],
			bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);

	}
	fclose(fp);

	// 获取距离提取点最近的点的ID

	mindisp = 100000.0;
	mintrac = 100000.0;

	DispPoint.pt[0] = -0.001 + m_X;
	DispPoint.pt[1] = 0.001 + m_Y;
	DispPoint.pt[2] = 0+ m_Z;
	TracPoint.pt[0] = -0.001 + m_X;
	TracPoint.pt[1] = 0.001 + m_Y;
	TracPoint.pt[2] = 0+ m_Z;

	for (i = 0; i < NodeNum; ++i)
	{
		if (mindisp > Dist(DispPoint, m_NodeList[i]))
		{
			mindisp = Dist(DispPoint, m_NodeList[i]);
			MinDispNodeID = i;
		}

		if (mintrac > Dist(TracPoint, m_NodeList[i]))
		{
			mintrac = Dist(TracPoint, m_NodeList[i]);
			MinTracNodeID = i;
		}
	}
	Temp = "Disp-Trac-Versus-Time_4.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();
	fopen_s(&fp, OutPath, "w");

	for (i = 0; i <= NStep; ++i)
	{
		fprintf(fp, "%20.10e   %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
			i * dt, m_NodeList[MinDispNodeID].pt[0], m_NodeList[MinDispNodeID].pt[1], m_NodeList[MinDispNodeID].pt[2], bd[i].lu.b[3 * MinDispNodeID + 0], bd[i].lu.b[3 * MinDispNodeID + 1], bd[i].lu.b[3 * MinDispNodeID + 2],
			bd[i].lt.b[3 * MinTracNodeID + 0], bd[i].lt.b[3 * MinTracNodeID + 1], bd[i].lt.b[3 * MinTracNodeID + 2]);

	}
	fclose(fp);

	return 1;
}

void showMemoryInfo() {//输出内存使用量
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	//std::string strs;
	printf("memory used %d M/ %d M  +  usage %d M/ %d M\n", pmc.WorkingSetSize / 1024 / 1024, pmc.PeakWorkingSetSize / 1024 / 1024, pmc.PagefileUsage / 1024 / 1024, pmc.PeakPagefileUsage / 1024 / 1024);
}

int GetInfoFromPoint_OutPut(OutPoint* OutP, OutEle* OutE, BoundaryValue* bd, long** m_EleNID, long OutPoNum, long EleNum, long NStep)
{
	long i, j;
	FILE* fp;
	char nameofstep[10], openfile[30], tempchar[20] = ".//result//";
	int itemp, jtemp;
	int tempNodeID;
	char* CurrPath = _getcwd(NULL, 0);//返回当前路径
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//返回当前文件夹名称
	Path = StrCurrPath.substr(0, found + 1);//返回上级路径
	string TempPath = DBEMOutputDirectory() + "\\";
	string TempNum;
	string ForMatTxt = ".dat";//文件后缀    Git test
	const char* OutPath;
	long tempf, tempEle;
	double tempValue;


	for (i = 0; i <= NStep; i = i + 1)
	{
		FILE* fp;
		TempNum = to_string(i);

		Path = TempPath + TempNum + ForMatTxt;
		OutPath = Path.c_str();
		fopen_s(&fp, OutPath, "w");
		//strcpy_s(openfile, tempchar);
		//_itoa_s(i, nameofstep, 10);
		//strcat_s(openfile, nameofstep);
		//fopen_s(&fp, openfile, "w");

		if (i <= NStep)//最后一步有问题，暂复制前一步结果
		{
			fprintf(fp, "Title=\"Exapmle 3D-U  %5d\"\n", i);
			fprintf(fp, "\"x\", \"y\", \"z\", \"U\"\n");
			fprintf(fp, "zone N =%6d, E =%6d, f = fepoint, ET = QUADRILATERAL\n", OutPoNum, EleNum);

			for (itemp = 0; itemp < OutPoNum; itemp++)
			{

				//for (jtemp = 0; jtemp < 8; jtemp++)
				//{
				//	tempNodeID = m_EleNID[itemp][jtemp];
				//	fprintf(fp, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", m_NodeList[tempNodeID].pt[0], m_NodeList[tempNodeID].pt[1], m_NodeList[tempNodeID].pt[2], bd[i].lu.b[3 * tempNodeID + 0], bd[i].lu.b[3 * tempNodeID + 1], bd[i].lu.b[3 * tempNodeID + 2], bd[i].lt.b[3 * tempNodeID + 0], bd[i].lt.b[3 * tempNodeID + 1], bd[i].lt.b[3 * tempNodeID + 2]);
				//	
				//}


				//输出物理量
				tempValue = 0;
				tempf = OutP[itemp].flag;
				for (j = 0; j < tempf; j++) {
					tempEle = OutP[itemp].EleList[j];
					tempNodeID = m_EleNID[tempEle][0];
					tempValue = tempValue + sqrt(bd[i].lu.b[3 * tempNodeID + 0] * bd[i].lu.b[3 * tempNodeID + 0] + bd[i].lu.b[3 * tempNodeID + 1] * bd[i].lu.b[3 * tempNodeID + 1] + bd[i].lu.b[3 * tempNodeID + 2] * bd[i].lu.b[3 * tempNodeID + 2]);
				}
				tempValue = tempValue / tempf;




				fprintf(fp, "%10.5f %10.5f %10.5f %10.5f\n", OutP[itemp].pt[0], OutP[itemp].pt[1], OutP[itemp].pt[2], tempValue);
			}

			for (itemp = 0; itemp < EleNum; itemp++) {
				fprintf(fp, "%10d %10d %10d %10d\n", OutE[itemp].PointList[0], OutE[itemp].PointList[1], OutE[itemp].PointList[2], OutE[itemp].PointList[3]);
			}
		}

		fclose(fp);

	}
	return 1;
}
