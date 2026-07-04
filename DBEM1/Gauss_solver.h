#ifndef _GAUSSSOLVER_
#define _GAUSSSOLVER_

#include "mat_vec.h"
#include "DSquareElement.h"
#include "DBEM.h"
//---------------------------------------------------------------------------

int StaticGaussSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, double MaxLength);
int StaticGaussSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum,long InfEleNum,long **m_InfElePID,long **m_ElePID, double MaxLength);
int DynaGaussSolver(DSquareElement* m_DSE,BoundaryValue* bd,long NodeNum,long EleNum,long NStep,double MaxLength);
int DynaGaussSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID);
int DynaGaussSolverNew(DSquareElement* m_DSE,BoundaryValue* bd,long NodeNum,long EleNum,long NStep,double MaxLength);
int DynaGaussSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID);
int DynaGaussSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID,int m_flag);
int DynaGaussSolverNewCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID, int m_flag);
int DynaGaussSolverNewCCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID, int m_flag);
// 获取动态基本解对应矩阵，矩阵系数转换到了局部坐标系

// (已逐行检查)
int GetMatrixDynaUij(DSquareElement* m_DSE,Wmatrix& A,long NodeNum,long EleNum,double n,double dt);
int GetMatrixDynaUij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, double n, double dt);
int GetMatrixDynaUij(DSquareElement* m_DSE, WmatrixSymCCSR& A, long NodeNum, long EleNum, double n, double dt);

int GetMatrixDynaT1ij(DSquareElement* m_DSE,Wmatrix& A,long NodeNum,long EleNum,double n,double dt);
int GetMatrixDynaT1ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID);
int GetMatrixDynaT1ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID);
int GetMatrixDynaT1ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID);

int GetMatrixDynaT2ij(DSquareElement* m_DSE,Wmatrix& A,long NodeNum,long EleNum,double n,double dt);
int GetMatrixDynaT2ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID);
int GetMatrixDynaT2ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID);
int GetMatrixDynaT2ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID);

int GetMatrixDynaTij(DSquareElement* m_DSE,Wmatrix& A,long NodeNum,long EleNum,double n,double dt);

// (已逐行检查)
int GetMatrixStaticUnknown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum);
int GetMatrixStaticKnown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum);
int GetMatrixStaticUnknown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum,long** InfElePID, long** m_ElePID);
int GetMatrixStaticKnown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum,long** InfElePID, long** m_ElePID);

int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE,Wmatrix& A,long NodeNum,long EleNum);
int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID);

int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID);
int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID);

int GetMatrixDynaKnown0ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID);
int GetMatrixDynaKnown0ij(DSquareElement* m_DSE,Wmatrix& A,long NodeNum,long EleNum);
int GetMatrixDynaKnown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID);
int GetMatrixDynaKnown0ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID);

//---------------------------------------------------------------------------
int GetInfoFromMatrixNstep(Wmatrix& MG, Wmatrix& MT, long Step, long NodeNum);
#endif