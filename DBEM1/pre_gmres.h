#ifndef _PREGMRES
#define _PREGMRES

#include "DSquareElement.h"
#include "mat_vec.h"
#include "lin_comp.h"
#include "tree.h"
#include "DBEM.h"
//-------------------------------------------------------------------------------------------------
//precondictioned GMRES迭代法

class PreConditioner
{
public:
	// variables
	long m_LeafCount;
	long* m_LeafPointCount;
	long* m_LeafBeginID;
	long* m_RePID;
	long* IID;
	long* m_InputPID;
	long* m_OutputPID;
	Wmatrix* m_PreM;
	
	// functions
	PreConditioner()
	{
		m_LeafPointCount = 0; 
		m_LeafBeginID = 0;
		m_RePID = 0;
		IID = 0;
		m_InputPID = 0;
		m_OutputPID = 0;
		m_PreM = 0;
	}
	
	~PreConditioner()
	{
		if(m_LeafPointCount)
		{
			delete[] m_LeafPointCount;
			m_LeafPointCount = 0;
		}
		if(m_LeafBeginID)
		{
			delete[] m_LeafBeginID;
			m_LeafBeginID = 0;
		}
		if(m_RePID)
		{
			delete[] m_RePID;
			m_RePID = 0;
		}
		if(IID)
		{
			delete[] IID;
			IID = 0;
		}
		if(m_InputPID)
		{
			delete[] m_InputPID;
			m_InputPID = 0;
		}
		if(m_OutputPID)
		{
			delete[] m_OutputPID;
			m_OutputPID = 0;
		}
		if(m_PreM)
		{
			delete[] m_PreM;
			m_PreM = 0;
		}
	}

	int Finish()
	{
		if(m_LeafPointCount)
		{
			delete[] m_LeafPointCount;
			m_LeafPointCount = 0;
		}
		if(m_LeafBeginID)
		{
			delete[] m_LeafBeginID;
			m_LeafBeginID = 0;
		}
		if(m_RePID)
		{
			delete[] m_RePID;
			m_RePID = 0;
		}
		if(IID)
		{
			delete[] IID;
			IID = 0;
		}
		if(m_InputPID)
		{
			delete[] m_InputPID;
			m_InputPID = 0;
		}
		if(m_OutputPID)
		{
			delete[] m_OutputPID;
			m_OutputPID = 0;
		}
		if(m_PreM)
		{
			delete[] m_PreM;
			m_PreM = 0;
		}

		return 1;
	}

};

int GMRESPreConditioner(DSquareElement* m_DSE,PreConditioner& m_Pre,long Count_PT,long MaxLeafPointCount,Wmatrix& A);
int GMRESPreConditioner(DSquareElement* m_DSE, PreConditioner& m_Pre, long Count_PT, long MaxLeafPointCount, CCSRMat& MT);

int GMRES_M_X(PreConditioner& m_Pre,Wvector& X,Wvector& B);

int gmres(DSquareElement* m_DSE,Wmatrix& A,Wvector& b,PreConditioner& m_Pre,long int restart,double tol,long int maxit,Wvector& x0,Wvector& x,int& flag,int* iter);
int gmres(DSquareElement* m_DSE, CSRMat& A, Wvector& b, PreConditioner& m_Pre, long int restart, double tol, long int maxit, Wvector& x0, Wvector& x, int& flag, int* iter);
int gmres(DSquareElement* m_DSE, CCSRMat& A, Wvector& b, PreConditioner& m_Pre, long int restart, double tol, long int maxit, Wvector& x0, Wvector& x, int& flag, int* iter);

int DynaGMRESSolver(DSquareElement* m_DSE,BoundaryValue* bd,long NodeNum,long EleNum,long iterations,double error,long maxleafpointnum,long NStep,double MaxLength);
int DynaGMRESSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID);


int DynaGMRESSolverNew(DSquareElement* m_DSE,BoundaryValue* bd,long NodeNum,long EleNum,long iterations,double error,long maxleafpointnum,long NStep,double MaxLength);
int DynaGMRESSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID);
int DynaGMRESSolverNewCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID);
int DynaGMRESSolverNewCCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID);
int DynaGMRESSolverNewCCSRNewPth(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long iterations, double error, long maxleafpointnum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID);
void GmresCcsrBeginMemorySession();
void GmresCcsrEndMemorySession();
void GmresCcsrSetOutputSubdir(const char* subdir);
//-------------------------------------------------------------------------------------------------
#endif
