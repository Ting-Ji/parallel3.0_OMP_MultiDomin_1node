#ifndef _lIN_COMP
#define _LIN_COMP

#include "mat_vec.h"
//-----------------------------------------------------
//高斯消去法，带右手向量
int gauss_mat(Wmatrix& A,Wvector& b,long n);

//矩阵求逆
int inv_mat(Wmatrix& A,long n);

//求向量内积
double d_dot(Wvector& a,Wvector& b,long n);

//求向量2范数
double vector_norm2(Wvector& a,long n);

//求矩阵向量相乘
int mat_vec_product(Wmatrix& A,Wvector& x,Wvector& b);

//求矩阵-矩阵相乘
int MMMul(Wmatrix& A,Wmatrix& B,Wmatrix& C);

//取矩阵的最大值和位置
int ArgMax(Wmatrix& A,long& i,long& j,double& value);

//取矩阵的最大值和位置
int ArgMax(Wvector& b,long& i,double& value);

//取矩阵的F范数
int normF(Wmatrix& A,double& value);

//取向量的F范数
int normF(Wvector& b,double& value);

//-----------------------------------------------------
#endif