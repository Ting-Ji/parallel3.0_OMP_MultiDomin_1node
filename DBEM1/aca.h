#ifndef _ACA
#define _ACA 

//--------------------------------------------------------------
//
//      ACA模板类的定义，可用于各种数据类型
//
//--------------------------------------------------------------

#include "geo.h"
#include "tree.h"
#include "mat_vec.h"
#include "lin_comp.h"
#include "stdio.h"
#include "dsquareelement.h"
#include "DBEM.h"

//--------------------------------------------------------------
//
//       Declaration of ACA classes and functions
//
//--------------------------------------------------------------


class ACAEntry
{
public:
	// 变量定义

	int flag_zero;       // 是否完全是零矩阵  0: 是  1: 非零

	int flag;            // 是否ACA运算
	int flag_overflow;   // 是否溢出
	long mbeginID;
	long nbeginID;
	long m;
	long n;
	long num;
	Wmatrix Entry;
	Wvector* u;
	Wvector* v;

	static long* IID;     // IID = T_DIMENSION * i
	static long* VID;     // VID = i / T_DIMENSION
	static int* VID_RC;
	static long MAX_M;
	static long MAX_N;

	static double errorvalue;
	static double errorstop;

	// 函数定义

	ACAEntry();
	~ACAEntry();

	double ACAMemory();
	int FullACA();
	int WriteToFile(FILE* fp,long i);
	int ReadFromFile(FILE* fp);

	ACAEntry& operator+=(ACAEntry& o_ACAEntry);

	static int StaticInit(long PointCount, double o_errorvalue, double o_errorstop);
	static int StaticRelease();
};

class ACAMatrix
{
public:
	long EntryCount;
	ACAEntry* m_ACAEntry;

	char MatrixName[100];

	static long* RePID;
	static int Flag_ReadWrite;    // 0 不进行读写  1 进行读写

	// 函数定义

	ACAMatrix()
	{
		m_ACAEntry = 0;
	}

	~ACAMatrix()
	{
		if (m_ACAEntry)
		{
			delete[] m_ACAEntry;
			m_ACAEntry = 0;
		}

		if (Flag_ReadWrite)
		{
			FILE* fp;
			fopen_s(&fp, MatrixName, "w");
			fprintf(fp, "\n");
			fclose(fp);
		}
	}

	ACAMatrix& operator+=(ACAMatrix& o_ACAMatrix);

	int Memory(double& ACAMem, double& ClassicMem);

	int WriteToFile(char* filename, long ID);
	int WriteToFile();
	int ReadFromFile();

	static int StaticRelease()
	{
		if (RePID)
		{
			delete[] RePID;
			RePID = 0;
		}

		return 1;
	}
};

class ACAPreConditioner
{
public:
	// variables
	long m_LeafCount;
	long* m_LeafPointCount;
	long* m_LeafBeginID;
	long* m_RePID;
	long* IID;
	Wmatrix* m_PreM;

	// functions
	ACAPreConditioner()
	{
		m_LeafPointCount = 0;
		m_LeafBeginID = 0;
		m_RePID = 0;
		IID = 0;
		m_PreM = 0;
	}

	~ACAPreConditioner()
	{
		if (m_LeafPointCount)
		{
			delete[] m_LeafPointCount;
			m_LeafPointCount = 0;
		}
		if (m_LeafBeginID)
		{
			delete[] m_LeafBeginID;
			m_LeafBeginID = 0;
		}
		if (m_RePID)
		{
			delete[] m_RePID;
			m_RePID = 0;
		}
		if (IID)
		{
			delete[] IID;
			IID = 0;
		}
		if (m_PreM)
		{
			delete[] m_PreM;
			m_PreM = 0;
		}
	}

	int Finish()
	{
		if (m_LeafPointCount)
		{
			delete[] m_LeafPointCount;
			m_LeafPointCount = 0;
		}
		if (m_LeafBeginID)
		{
			delete[] m_LeafBeginID;
			m_LeafBeginID = 0;
		}
		if (m_RePID)
		{
			delete[] m_RePID;
			m_RePID = 0;
		}
		if (IID)
		{
			delete[] IID;
			IID = 0;
		}
		if (m_PreM)
		{
			delete[] m_PreM;
			m_PreM = 0;
		}

		return 1;
	}

};

// 计算预处理矩阵
int ACAGMRESPreConditioner(DSquareElement* m_DSE, ACAPreConditioner& m_Pre, Point* m_PointList, long MaxLeafPointCount);

// 计算MX
int ACA_MX(ACAPreConditioner& m_Pre, Wvector& X, Wvector& B);

// 计算树中的ACA矩阵数量
void ACAMatrixCount(Tree* &m_tree, long& sum);

// 确定ACA矩阵信息
void ACAInfo(Tree* &m_tree, long& sum, ACAMatrix& m_ACAMatrix);

// 初始化ACA矩阵
void InitACAEntry(ACAMatrix& m_ACAMatrix);

// 设置ACA矩阵和RePID
int Setup_ACAEntry_RePID(ACAMatrix& m_ACAMatrix, Point* m_PointList, long Count_PT, long MaxLeafPointCount);

// 利用一个ACA矩阵初始化另外一个ACA矩阵
int GetInfoFromACAEntry(ACAMatrix& m_ACAMatrix1, ACAMatrix& m_ACAMatrix2);

// 计算基于树的ACA矩阵
void CalACAEntry(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, int flag);

// 计算基于树的不同n的ACA矩阵
void CalACAEntryUij(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, double n, double dt);
void CalACAEntryT1ij(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, double n, double dt);
void CalACAEntryT2ij(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, double n, double dt);

// 释放基于树的ACA的U和V向量
void ReleaseUV(ACAMatrix& m_ACAMatrix);

// ---------------------以下函数有待进一步检查---------------------------

// 计算AX
int ACA_BuildAX(ACAMatrix& m_ACAMatrix, Wvector& X, Wvector& AX);

// GMRES求解器
int ACA_GMRES(ACAMatrix& m_ACAMatrix, Wvector& b, Wvector& x0, Wvector& x, ACAPreConditioner& m_Pre, long restart, double tol, long maxit, int& flag, int* iter);

// Dyna ACA求解器
int DynaACASolver(DSquareElement* m_DSE, Point* m_NodeList, Point* m_EleCenterList, BoundaryValue* bd, long NodeNum, long EleNum, long MaxIterations, long MaxLeafPointCount, long MaxLeafPointCountPre, double Error, double ACAErrorValue, double ACAErrorStop, long NStep, double MaxLength, int Flag_ReadWrite);
int DynaACASolverNew(DSquareElement* m_DSE, Point* m_NodeList, Point* m_EleCenterList, BoundaryValue* bd, long NodeNum, long EleNum, long MaxIterations, long MaxLeafPointCount, long MaxLeafPointCountPre, double Error, double ACAErrorValue, double ACAErrorStop, long NStep, double MaxLength, int Flag_ReadWrite);
int DynaACASolverMem(DSquareElement* m_DSE, Point* m_NodeList, Point* m_EleCenterList, BoundaryValue* bd, long NodeNum, long EleNum, long MaxIterations, long MaxLeafPointCount, long MaxLeafPointCountPre, double Error, double ACAErrorValue, double ACAErrorStop, long NStep, double MaxLength, int Flag_ReadWrite);


//--------------------------------------------------------------

#endif

