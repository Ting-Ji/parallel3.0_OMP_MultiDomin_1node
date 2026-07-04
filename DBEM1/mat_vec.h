#ifndef MAT_VEC_
#define MAT_VEC_
#include"vector"
#include <omp.h>
#include <stdio.h>
using namespace std;
//---------------------------------------------------------------------------------------
const int IDENTITY = 1;
const int ZERO = 0;

class Wmatrix
{
public:
	//variables
	long m, n;
	double** a;
	//functions
	Wmatrix();
	Wmatrix(long mm, long nn);
	Wmatrix(long mm, long nn, int flag);
	Wmatrix(const Wmatrix& aa);
	~Wmatrix();
	void initial(long mm, long nn);
	void initial(long mm, long nn, int flag);
	void finish();
	Wmatrix& operator=(const Wmatrix& aa);
	Wmatrix& operator=(double* aa);
	Wmatrix& operator=(const double aa);
	double& operator()(long i, long j) { return a[i - 1][j - 1]; }
	void assign(Wmatrix& sub, long br, long bc, long m, long n);
	void assign(double* aa, long br, long bc, long m, long n);
	void Trans();
	void Minus();
	void Add(Wmatrix& aa);
	void Add(Wmatrix& sub, long br, long bc, long mm, long nn);

	Wmatrix& operator+=(double v);
	Wmatrix& operator-=(double v);
	Wmatrix& operator+=(Wmatrix& v);
	Wmatrix& operator-=(Wmatrix& v);
	Wmatrix& operator*=(Wmatrix& v);
	Wmatrix& operator*=(double v);
	Wmatrix& operator/=(Wmatrix& v);
	Wmatrix& operator/=(double v);
};

class Wvector
{
public:
	//variables
	long n;
	double* b;
	//functions
	Wvector();
	Wvector(long nn);
	Wvector(long nn, int flag);
	Wvector(const Wvector& bb);
	~Wvector();
	void initial(long nn);
	void finish();
	Wvector& operator=(const Wvector& bb);
	Wvector& operator=(const double bb);
	double& operator[](long i) { return b[i - 1]; }
	double& operator[](int i) { return b[i - 1]; }
	void assign(Wvector& sub, long br, long n);

	Wvector& operator+=(double v);
	Wvector& operator-=(double v);
	Wvector& operator+=(Wvector& v);
	Wvector& operator-=(Wvector& v);
	Wvector& operator*=(Wvector& v);
	Wvector& operator*=(double v);
	Wvector& operator/=(Wvector& v);
	Wvector& operator/=(double v);

	int plus(Wvector& v, double c);
	int minus(Wvector& v, double c);

	int OutputData(char* Name, long ID);
	int OutputData(char* Name, long ID, double* Value);
};

class WmatrixCSR
{
public:
	//variables
	long m, n;
	long num;//非零元素个数
	int* row;//每行非零元素
	double** a;
	int** col;//非零元素对应列号
	//functions
	WmatrixCSR();
	WmatrixCSR(long mm, long nn);
	WmatrixCSR(long mm, long nn, int flag);
	WmatrixCSR(const WmatrixCSR& aa);
	~WmatrixCSR();
	void initial(long mm, long nn);//新加—>对num赋值0（初始化）
	void sum();
	void reinitial();//重新初始化：赋初值0
	void initial(long mm, long nn, int flag);
	void finish();
	WmatrixCSR& operator=(const WmatrixCSR& aa);
	WmatrixCSR& operator=(double* aa);
	WmatrixCSR& operator=(const double aa);
	//WmatrixCSR& operator=(const int aa);
	double& operator()(long i, long j) { return a[i - 1][j - 1]; }
	void assign(WmatrixCSR& sub, long br, long bc, long m, long n);
	void assign(double* aa, long br, long bc, long m, long n);
	void Trans();
	void Minus();
	void Add(WmatrixCSR& aa);
	void Add(WmatrixCSR& sub, long br, long bc, long mm, long nn);

	WmatrixCSR& operator+=(double v);
	WmatrixCSR& operator-=(double v);
	WmatrixCSR& operator+=(WmatrixCSR& v);
	WmatrixCSR& operator-=(WmatrixCSR& v);
	WmatrixCSR& operator*=(WmatrixCSR& v);
	WmatrixCSR& operator*=(double v);
	WmatrixCSR& operator/=(WmatrixCSR& v);
	WmatrixCSR& operator/=(double v);
};

class WvectorCSR
{
public:
	//variables
	long n;
	double* value;
	int* col;
	//functions
	void initial(long nn);
	void reinitial();
	void finish();

};

class CSRMat
{
public:
	//variables
	long n;//数组行数
	CSRMat();
	//long num;//非零元素个数
	int* row_ptr;//矩阵每行第一列在value中对应的位置，最后一个数为value元素总数
	int* col_index;//矩阵列序号  从0开始
	double* value;//矩阵值     
	~CSRMat();
	//functions 
	void initial(long nn, long nnum);//初始化
	void Matdelete();
	//void assign(double a, int ncol, int nnum);//赋值 
	//void assign(double a, int ncol, int nnum, int nrow);//赋值 —>第nrow行第一个非零元素

};

class WmatrixCCSR
{
public:
	//variables
	long m, n;  //细胞行列		总行号为 原始行/3；总列号为 原始列号/3
	long num;//非零元素个数
	int* row;//每细胞行非零细胞个数		总行号为 原始行/3
	double** a;
	int** col;//非零细胞对应细胞列号    总列号为 原始列号/3
	//functions
	WmatrixCCSR();
	WmatrixCCSR(long mm, long nn);
	WmatrixCCSR(long mm, long nn, int flag);
	WmatrixCCSR(const WmatrixCCSR& aa);
	~WmatrixCCSR();
	void initial(long mm, long nn);//新加—>对num赋值0（初始化）
	void sum();//计算非零的细胞矩阵个数
	void reinitial();//重新初始化：赋初值0
	void initial(long mm, long nn, int flag);
	void finish();
	WmatrixCCSR& operator=(const WmatrixCCSR& aa);
	WmatrixCCSR& operator=(double* aa);
	WmatrixCCSR& operator=(const double aa);
	//WmatrixCCSR& operator=(const int aa);
	double& operator()(long i, long j) { return a[i - 1][j - 1]; }
	void assign(WmatrixCCSR& sub, long br, long bc, long m, long n);
	void assign(double* aa, long br, long bc, long m, long n);
	void assignSym(double* aa, long br, long bc, long m, long n);
	void Trans();
	void Minus();
	void Add(WmatrixCCSR& aa);
	void Add(WmatrixCCSR& sub, long br, long bc, long mm, long nn);

	WmatrixCCSR& operator+=(double v);
	WmatrixCCSR& operator-=(double v);
	WmatrixCCSR& operator+=(WmatrixCCSR& v);
	WmatrixCCSR& operator-=(WmatrixCCSR& v);
	WmatrixCCSR& operator*=(WmatrixCCSR& v);
	WmatrixCCSR& operator*=(double v);
	WmatrixCCSR& operator/=(WmatrixCCSR& v);
	WmatrixCCSR& operator/=(double v);
};

class WvectorCCSR
{
public:
	//variables
	long n;
	double* value;
	int* col;
	//functions
	WvectorCCSR();
	void initial(long nn);
	void reinitial();
	void finish();
	~WvectorCCSR();

};




class CCSRMat
{
public:
	//variables
	long n;//“细胞”行数
	//long num;//非零元素个数
	int* row_ptr;//矩阵每超行第一个超列对应的位置，最后一个数为非零子矩阵总数
	int* col_index;//矩阵列序号  从0开始
	double* value;//矩阵值      
	//functions 
	CCSRMat();//构造！！数组需初始化！！
	~CCSRMat();
	void initial(long nn, long nnum);//初始化  nnum非零细胞矩阵个数
	void Matdelete();
	long NonZeroBlocks() const;
	bool WriteBinary(FILE* fp) const;
	bool ReadBinary(FILE* fp);
	//void assign(double a, int ncol, int nnum);//赋值 
	//void assign(double a, int ncol, int nnum, int nrow);//赋值 —>第nrow行第一个非零元素

};

class WmatrixSymCCSR
{
public:
	//variables
	long m, n;  //细胞行列		总行号为 原始行/3；总列号为 原始列号/3
	long num;//非零元素个数
	int* row;//每细胞行非零细胞个数		总行号为 原始行/3
	double** a;
	int** col;//非零细胞对应细胞列号    总列号为 原始列号/3
	//functions
	WmatrixSymCCSR();
	~WmatrixSymCCSR();
	void initial(long mm, long nn);//新加—>对num赋值0（初始化）
	void sum();//计算非零的细胞矩阵个数
	void reinitial();//重新初始化：赋初值0
	void finish();
	//WmatrixSymCCSR& operator=(const int aa);
	void assign(WmatrixCCSR& sub, long br, long bc, long m, long n);
	void assign(double* aa, long br, long bc, long m, long n);
};

class WvectorSymCCSR
{
public:
	//variables
	long n;
	double* value;
	int* col;
	//functions
	WvectorSymCCSR();
	void initial(long nn);
	void reinitial();
	void finish();
	~WvectorSymCCSR();
};

class SymCCSRMat
{
public:
	//variables
	long n;//“细胞”行数
	//long num;//非零元素个数
	SymCCSRMat();//析构！！数组需初始化！！
	int* row_ptr;//矩阵每细胞行第一列在value中对应的位置，最后一个数为value元素总数
	int* col_index;//矩阵列序号  从0开始
	double* value;//矩阵值      
	//functions 
	~SymCCSRMat();
	void initial(long nn, long nnum);//初始化  nnum非零细胞矩阵个数
	void Matdelete();
	long NonZeroBlocks() const;
	bool WriteBinary(FILE* fp) const;
	bool ReadBinary(FILE* fp);
	//void assign(double a, int ncol, int nnum);//赋值 
	//void assign(double a, int ncol, int nnum, int nrow);//赋值 —>第nrow行第一个非零元素

};



class PthArg
{
public:
	long targs;
	long start;
	long end;
	long sum;
	long AddStrat;
	long AddEnd;
	long AddSum;
	PthArg();
	PthArg& operator=(PthArg& M);
	~PthArg();
	//long Distr;
};

class PthConvol //卷积并行任务划分
{
public:
	long targs;
	long PthStep;
	long MaxN;
	long* start;
	long* end;
	PthConvol();
	void initial(long nn);
	~PthConvol();
	//long Distr;
};



//稀疏存储转化
void SparseConvert(WmatrixCCSR& A, CCSRMat& M);// CCSR
void SparseConvert(WmatrixSymCCSR& A, SymCCSRMat& M);//对称CCSR

//稀疏矩阵相加
void SparseAddition(WmatrixCCSR& A, CCSRMat& M, WvectorCCSR& b);//A为稀疏矩阵1，M为稀疏矩阵2，M2=M1+M2，返回稀疏矩阵2
void SparseAddition(WmatrixSymCCSR& A, SymCCSRMat& M, WvectorSymCCSR& b);

//稀疏矩阵乘向量
int SparseMul(CCSRMat& M, Wvector& x, Wvector& b);
int SparseMul(SymCCSRMat& M, Wvector& x, Wvector& b);


//稀疏存储转化
void SparseConvert(WmatrixCSR& A, CSRMat& M);


//稀疏矩阵相加
void SparseAddition(WmatrixCSR& A, CSRMat& M, WvectorCSR& b);//A为稀疏矩阵1，M为稀疏矩阵2，M2=M1+M2，返回稀疏矩阵2


//稀疏矩阵乘向量
int SparseMul(CSRMat& M, Wvector& x, Wvector& b);







//---------------------------------------------------------------------------------------
#endif
