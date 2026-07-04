#ifndef _PTHREADCCSR_
#define _PTHREADCCSR_
#include "Gauss_solver.h"
#include "dsquareelement.h"
#include "DBEM.h"
using namespace std;



void *PthGetMatrixDynaUnknown0ij(void *arg);

void* PthGetMatrixDynaUnknown0ijT0(void* arg);

void *PthGetMatrixDynaKnown0ij(void* arg);

void *PthGetMatrixDynaT1ij(void* arg);
void* PthGetMatrixDynaT1ijNew(void* arg);

void *PthGetMatrixDynaT2ij(void* arg);
void* PthGetMatrixDynaT2ijNew(void* arg);

void *PthGetMatrixDynaUij(void* arg);

void* PthGetMatrixDynaSymUij(void* arg);
void* PthGetMatrixDynaSymUijNew(void* arg);

void* PthSparseConvertCCSR(void* arg);
void* PthSparseAddConvertCCSR(void* arg);

void* PthSparseConvertCCSRMTS0(void* arg);

void* PthSparseConvertSymCCSR(void* arg);

void* PthSparseAdditionCCSR(void* arg);

void* PthDistr(int* row,PthArg* Thread, long num,long n);
void* PthDistrAdd(int* row, PthArg* ThreadAdd, long num, long n);
void* PthDistr(PthArg* Thread);

void* PthAllTraverse(void* arg);


void* PthSpMVKnown(void* arg);
void* PthSpMVUnkown(void* arg);
void* PthSpMVT(void* arg);
void* PthSpMVG(void* arg);
void* PthNorm(void* arg);

void* PthbdConvert(void* arg);
void* PthbdAllLocToAbs(void* arg);


void* SpmvConvertConvol(PthArg** Spmvthread, PthConvol* Convolthread, long MaxN);
void* PthConvolutinT(void* arg);
void* PthConvolutinG(void* arg);
void* ConvolBvector(Wvector* CBT, Wvector* CBG, BoundaryValue* bd, long StepNum, long MaxN);
void* PthConvolBvector(void* arg);
void* PthCollecBA(void* arg);
void* PthCollecBS(void* arg);
void* PthBvectorInitial(void* arg);


void* PthDSEPhyInit(void* arg);

void GmresThreadDiagPrintSystemInfo(long requestedThreadNum);
void GmresThreadDiagBegin(long requestedThreadNum);
void GmresThreadDiagRecord(long threadID);
void GmresThreadDiagPrintCreateResult(long requestedThreadNum, long successCount, const int* returnCodes);
void GmresThreadDiagEnd();

#endif
