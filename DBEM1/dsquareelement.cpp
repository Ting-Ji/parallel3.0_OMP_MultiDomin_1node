#include "DSquareElement.h"
#include "lin_comp.h"
#include "stdio.h"
#include <float.h>
#include "dsquareelement.h"
#include "DBEM.h"
//realization of class: DSquareElement

//1.Initialization

quadinfo DSquareElement::m_quadinfo;
FValue DSquareElement::m_FValue;
SecInfo DSquareElement::m_SecInfo;

double DSquareElement::m_gap = 0.5;

GaussNumeMatrix DSquareElement::m_gnm;

double DSquareElement::m_OnEleN[8][4][SINGAUSSPOINT2][3];
double DSquareElement::m_OnEleJ[8][4][SINGAUSSPOINT2];
double DSquareElement::m_OnEleR[8][4][SINGAUSSPOINT2][4];

double DSquareElement::m_AssistUij[8][9];
double DSquareElement::m_AssistTij[8][9];
double DSquareElement::m_AssistCoef[8][9];

double DSquareElement::m_StaticTij[8][8][9];

Counter DSquareElement::m_counter;

long* DSquareElement::BeEleID;
int* DSquareElement::BeElePos;
double** DSquareElement::m_transmat;  //移至全局变量

DynaMat DSquareElement::m_DMat;

int DSquareElement::FlagDynaMethod;
int DSquareElement::Flagproblem;
double DSquareElement::InfCenter[3];

static int DSquareFinite(double value)
{
	return _finite(value) != 0;
}

static void DSquareZeroStress(double stress[3][3])
{
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			stress[i][j] = 0.0;
}

static int DSquareStressFinite(double stress[3][3])
{
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (!DSquareFinite(stress[i][j]))
				return 0;
	return 1;
}

DynaMat BuildDynaMat(double o_v, double o_G, double o_Rou, double o_dt)
{
	DynaMat mat;
	mat.v = o_v;
	mat.G = o_G;
	mat.Rou = o_Rou;
	mat.C1 = sqrt((2.0 - 2.0 * o_v) / (1.0 - 2.0 * o_v) * o_G / o_Rou);
	mat.C2 = sqrt(o_G / o_Rou);

	mat.Dt = o_dt;
	mat.Dt2 = 2.0 * o_dt;
	mat.Dt_2 = o_dt * o_dt;
	mat.Dt_3 = o_dt * o_dt * o_dt;
	mat.OneOverFourPaiRou = ONE_OVER_4PI / o_Rou;

	mat.C1_2 = mat.C1 * mat.C1;
	mat.C2_2 = mat.C2 * mat.C2;
	mat.OneOverC1_2 = 1.0 / mat.C1_2;
	mat.OneOverC2_2 = 1.0 / mat.C2_2;
	mat.OneOverC1_3 = 1.0 / mat.C1_2 / mat.C1;
	mat.OneOverC2_3 = 1.0 / mat.C2_2 / mat.C2;
	mat.C2_2OverC1_2 = mat.C2_2 / mat.C1_2;
	mat.C2_2OverC1_3 = mat.C2_2OverC1_2 / mat.C1;

	mat.C1Dt = mat.C1 * mat.Dt;
	mat.C2Dt = mat.C2 * mat.Dt;

	mat.STCoef1 = -1.0 / (8.0 * pi * (1.0 - mat.v));
	mat.STCoef2 = 1.0 - 2.0 * mat.v;

	mat.OneOverThree_1OC23_1OC13 = 1.0 / 3.0 * (mat.OneOverC2_3 - mat.OneOverC1_3);
	mat.Half_1OC22_1OC12 = 0.5 * (mat.OneOverC2_2 - mat.OneOverC1_2);

	return mat;
}

int BoundaryValue::Convert(int bdid, long n)
{
	double tlu[3], tlt[3];
	long n3 = 3 * n;

	switch (bdid)
	{
	case 123:
	{
		break;
	}
	case 126:
	{
		tlu[2] = lt.b[n3 + 2];
		tlt[2] = lu.b[n3 + 2];
		lu.b[n3 + 2] = tlu[2];
		lt.b[n3 + 2] = tlt[2];
		break;
	}
	case 153:
	{
		tlu[1] = lt.b[n3 + 1];
		tlt[1] = lu.b[n3 + 1];
		lu.b[n3 + 1] = tlu[1];
		lt.b[n3 + 1] = tlt[1];
		break;
	}
	case 156:
	{
		tlu[1] = lt.b[n3 + 1];
		tlt[1] = lu.b[n3 + 1];
		lu.b[n3 + 1] = tlu[1];
		lt.b[n3 + 1] = tlt[1];

		tlu[2] = lt.b[n3 + 2];
		tlt[2] = lu.b[n3 + 2];
		lu.b[n3 + 2] = tlu[2];
		lt.b[n3 + 2] = tlt[2];

		break;
	}
	case 423:
	{
		tlu[0] = lt.b[n3 + 0];
		tlt[0] = lu.b[n3 + 0];
		lu.b[n3 + 0] = tlu[0];
		lt.b[n3 + 0] = tlt[0];
		break;
	}
	case 426:
	{
		tlu[0] = lt.b[n3 + 0];
		tlt[0] = lu.b[n3 + 0];
		lu.b[n3 + 0] = tlu[0];
		lt.b[n3 + 0] = tlt[0];

		tlu[2] = lt.b[n3 + 2];
		tlt[2] = lu.b[n3 + 2];
		lu.b[n3 + 2] = tlu[2];
		lt.b[n3 + 2] = tlt[2];

		break;
	}
	case 453:
	{
		tlu[0] = lt.b[n3 + 0];
		tlt[0] = lu.b[n3 + 0];
		lu.b[n3 + 0] = tlu[0];
		lt.b[n3 + 0] = tlt[0];

		tlu[1] = lt.b[n3 + 1];
		tlt[1] = lu.b[n3 + 1];
		lu.b[n3 + 1] = tlu[1];
		lt.b[n3 + 1] = tlt[1];

		break;
	}
	case 456:
	{
		tlu[0] = lt.b[n3 + 0];
		tlt[0] = lu.b[n3 + 0];
		lu.b[n3 + 0] = tlu[0];
		lt.b[n3 + 0] = tlt[0];

		tlu[1] = lt.b[n3 + 1];
		tlt[1] = lu.b[n3 + 1];
		lu.b[n3 + 1] = tlu[1];
		lt.b[n3 + 1] = tlt[1];

		tlu[2] = lt.b[n3 + 2];
		tlt[2] = lu.b[n3 + 2];
		lu.b[n3 + 2] = tlu[2];
		lt.b[n3 + 2] = tlt[2];

		break;
	}
	}

	return 1;
}

int BoundaryValue::Convert(DSquareElement* m_DSE, long EleNum)
{
	long i;

	for (i = 0; i < EleNum; ++i)
		Convert(m_DSE[i].BCID, i);

	return 1;
}

int BoundaryValue::AllLocToAbs(long NodeNum, double** Trans)
{
	long i;

	for (i = 0; i < NodeNum; ++i)
		LocToAbs(i, Trans[i]);

	return 1;
}

int BoundaryValue::ToZero(double error)
{
	long i;
	double MaxLu, MaxLt;

	MaxLu = 0.0;
	MaxLt = 0.0;

	for (i = 0; i < lu.n; ++i)
	{
		if (MaxLu < fabs(lu.b[i]))
			MaxLu = fabs(lu.b[i]);

		if (MaxLt < fabs(lt.b[i]))
			MaxLt = fabs(lt.b[i]);
	}

	MaxLu *= error;
	MaxLt *= error;

	for (i = 0; i < lu.n; ++i)
	{
		if (fabs(lu.b[i]) < MaxLu)
			lu.b[i] = 0.0;

		if (fabs(lt.b[i]) < MaxLt)
			lt.b[i] = 0.0;
	}

	return 1;
}

void DSquareElement::StaticInit(double gap, double o_v, double o_G, double o_Rou, double o_dt, int o_FlagDynaMethod)
{
	m_gap = gap;
	FlagDynaMethod = o_FlagDynaMethod;

	// 材料参数赋值
	m_DMat = BuildDynaMat(o_v, o_G, o_Rou, o_dt);

	m_quadinfo.init(gap);
	m_quadinfo.AssignSecInfo(m_SecInfo);
	m_quadinfo.AssignFValue(m_FValue);
}

void DSquareElement::StaticGetBeInfo(long NodeNum, long EleNum, long** m_EleNID)
{
	long i, j;
	int pos;

	BeEleID = new long[NodeNum];
	BeElePos = new int[NodeNum];

	for (i = 0; i < EleNum; ++i)
	{
		for (j = 0; j < 1; ++j)
		{
			pos = (int)m_EleNID[i][j];
			BeEleID[pos] = i;
			BeElePos[pos] = (int)j;
		}
	}

	m_transmat = new double*[NodeNum];
	for (i = 0; i < NodeNum; ++i)
		m_transmat[i] = new double[9];
}

void DSquareElement::StaticRelease(long NodeNum)
{
	long i;

	if (BeEleID)
	{
		delete[] BeEleID;
		BeEleID = 0;
	}

	if (BeElePos)
	{
		delete[] BeElePos;
		BeElePos = 0;
	}

	if (m_transmat)
	{
		for (i = 0; i < NodeNum; ++i)
		{
			if (m_transmat[i])
			{
				delete[] m_transmat[i];
				m_transmat[i] = 0;
			}
		}
		delete[] m_transmat;
		m_transmat = 0;
	}

}

void DSquareElement::GeoInit(Point* pointlist, long* pointID)
{
	int i, j;

	m_pointlist = pointlist;
	for (i = 0; i < 8; ++i)
		m_pointID[i] = pointID[i];
	InfEleFlag = -1;
	//New Definations
	GetCoordCoef();

	//Jacobi values and normal vectors on area gauss position
	//for (i = 0; i < GAUSSPOINT2; ++i)
	//{
	//	//  ———————— 改为积分时计算———————
	//	m_Jacobi[i] = Jacobi(m_gnm.rp[i][0], m_gnm.rp[i][1]);
	//	Normal(m_gnm.rp[i][0], m_gnm.rp[i][1], m_normal[i][0], m_normal[i][1], m_normal[i][2]);
	//	GetFromLocal(m_gnm.rp[i][0], m_gnm.rp[i][1], m_intpt[i]);
	//	//  ———————— 改为积分时计算———————
	//}


	// 分片单元积分信息
	for (i = 0; i < NUMPW2; ++i)
	{
		//for (j = 0; j < GAUSSPOINTPW2; ++j)
		//{//  ———————— 改为积分时计算———————
		//	m_JacobiPW[i][j] = Jacobi(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1]);
		//	Normal(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1], m_normalPW[i][j][0], m_normalPW[i][j][1], m_normalPW[i][j][2]);
		//	GetFromLocal(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1], m_intptPW[i][j]);
		//	//  ———————— 改为积分时计算———————
		//}

		GetFromLocal(m_quadinfo.m_CENPW[i][0], m_quadinfo.m_CENPW[i][1], m_rcenPW[i]);// 预先生成
	}

	// 计算距离判断信息
	double m_dist;
	GetFromLocal(0.0, 0.0, m_center);
	m_radius = -1.0;

	for (i = 0; i < 8; ++i)
	{
		m_dist = Dist(m_center, m_pointlist[m_pointID[i]]);
		if (m_radius < m_dist)
			m_radius = m_dist;
	}

	m_radiusPW = m_radius / NUMPW;

}

void DSquareElement::PhyInit(Point* nodelist, long* nodeID, int boundaryID)
{
	int i;

	BCID = boundaryID;
	m_nodelist = nodelist;
	for (i = 0; i < 1; ++i)
		m_nodeID[i] = nodeID[i];

	for (i = 0; i < 1; ++i)
		TranslateMatrix(m_quadinfo.m_LocalCoord[i][0], m_quadinfo.m_LocalCoord[i][1], m_transmat[m_nodeID[i]]);

	if (Flagproblem == 1)
		InElementIntStatic();
	else
		InElementIntDynamic();
}

void DSquareElement::PthPhyInit(Point* nodelist, long* nodeID, int boundaryID)
{
	int i;

	BCID = boundaryID;
	m_nodelist = nodelist;
	for (i = 0; i < 1; ++i)
		m_nodeID[i] = nodeID[i];

	for (i = 0; i < 1; ++i)
		TranslateMatrix(m_quadinfo.m_LocalCoord[i][0], m_quadinfo.m_LocalCoord[i][1], m_transmat[m_nodeID[i]]);

	if (Flagproblem == 1)
		InElementIntStatic();
	else
		PthInElementIntDynamic();
}


double DSquareElement::N(double s1, double s2, int ID)
{
	switch (ID)
	{
	case 0:
		return 0.25*(1.0 - s2)*(s1 - 1.0)*(s1 + s2 + 1.0);
	case 1:
		return 0.25*(1.0 - s2)*(1.0 + s1)*(s1 - 1.0 - s2);
	case 2:
		return 0.25*(1.0 + s2)*(1.0 + s1)*(s1 - 1.0 + s2);
	case 3:
		return 0.25*(1.0 + s2)*(s1 - 1.0)*(s1 + 1.0 - s2);
	case 4:
		return 0.5*(s1 - 1.0)*(1.0 + s1)*(s2 - 1.0);
	case 5:
		return 0.5*(1.0 - s2)*(1.0 + s2)*(1.0 + s1);
	case 6:
		return 0.5*(1.0 - s1)*(1.0 + s1)*(1.0 + s2);
	case 7:
		return 0.5*(s2 - 1.0)*(1.0 + s2)*(s1 - 1.0);
	default:
		return 0.0;
	}
}

double DSquareElement::NDiff1(double s1, double s2, int ID)
{
	switch (ID)
	{
	case 0:
		return 0.25*(1.0 - s2)*(s2 + s1 + s1);
	case 1:
		return 0.25*(1.0 - s2)*(s1 + s1 - s2);
	case 2:
		return 0.25*(1.0 + s2)*(s2 + s1 + s1);
	case 3:
		return 0.25*(1.0 + s2)*(s1 + s1 - s2);
	case 4:
		return s1 * (s2 - 1.0);
	case 5:
		return 0.5*(1.0 - s2 * s2);
	case 6:
		return -s1 * (1.0 + s2);
	case 7:
		return 0.5*(s2 - 1.0)*(1 + s2);
	default:
		return 0.0;
	}
}

double DSquareElement::NDiff2(double s1, double s2, int ID)
{
	switch (ID)
	{
	case 0:
		return 0.25*(1.0 - s1)*(s1 + s2 + s2);
	case 1:
		return 0.25*(1.0 + s1)*(s2 + s2 - s1);
	case 2:
		return 0.25*(1.0 + s1)*(s1 + s2 + s2);
	case 3:
		return 0.25*(s1 - 1.0)*(s1 - s2 - s2);
	case 4:
		return 0.5*(s1*s1 - 1.0);
	case 5:
		return -(1.0 + s1)*s2;
	case 6:
		return 0.5*(1.0 - s1 * s1);
	case 7:
		return (s1 - 1.0)*s2;
	default:
		return 0.0;
	}
}

void DSquareElement::GetTij(double* T, double* R, double* n)
{
	//coef1=-1/(8*pi*(1-v))
	//coef2=1-2*v

	double t1 = R[0] * R[0];
	double t2 = m_DMat.STCoef1 / t1;
	double t3 = m_DMat.STCoef2 / R[0];
	double drdn = (R[1] * n[0] + R[2] * n[1] + R[3] * n[2]) / R[0];
	double t4, t5, t6, t7, t8, t9, t10;

	t5 = t3 * (R[2] * n[0] - R[1] * n[1]);
	t6 = t3 * (R[3] * n[0] - R[1] * n[2]);
	t7 = t3 * (R[3] * n[1] - R[2] * n[2]);
	t8 = t2 * drdn;
	t9 = t8 * m_DMat.STCoef2;
	t10 = t8 * 3.0 / t1;
	t4 = t10 * R[1];
	T[0] = t9 + t4 * R[1];
	T[1] = t4 * R[2] - t2 * t5;
	T[2] = t4 * R[3] - t2 * t6;
	T[3] = t4 * R[2] + t2 * t5;
	T[4] = t9 + t10 * R[2] * R[2];
	t4 = t10 * R[3];
	T[5] = t4 * R[2] - t2 * t7;
	T[6] = t4 * R[1] + t2 * t6;
	T[7] = t4 * R[2] + t2 * t7;
	T[8] = t9 + t10 * R[3] * R[3];
	//	double t1=R[0]*R[0];
	//	double t2=coef1/t1;
	//	double t3=coef2/R[0];
	//	double drdn=(R[1]*n[0]+R[2]*n[1]+R[3]*n[2])/R[0];
	//	double t4,t5,t6,t7;
	//	t4=drdn*3.0*R[1]/t1;
	//	t5=t3*(R[2]*n[0]-R[1]*n[1]);
	//	t6=t3*(R[3]*n[0]-R[1]*n[2]);
	//	t7=t3*(R[3]*n[1]-R[2]*n[2]);
	//	T[0]=t2*drdn*(coef2+3.0*R[1]*R[1]/t1);
	//	T[1]=t2*(t4*R[2]-t5);
	//	T[2]=t2*(t4*R[3]-t6);
	//	T[3]=t2*(t4*R[2]+t5);
	//	T[4]=t2*drdn*(coef2+3.0*R[2]*R[2]/t1);
	//	t4=drdn*3.0*R[3]/t1;
	//	T[5]=t2*(t4*R[2]-t7);
	//	T[6]=t2*(t4*R[1]+t6);
	//	T[7]=t2*(t4*R[2]+t7);
	//	T[8]=t2*drdn*(coef2+3.0*R[3]*R[3]/t1);
}
//没使用
int DSquareElement::GetUij(double* U, double* R)
{
	U[0] = 0.0;
	U[1] = 0.0;
	U[2] = 0.0;
	U[3] = 0.0;
	U[4] = 0.0;
	U[5] = 0.0;
	U[6] = 0.0;
	U[7] = 0.0;
	U[8] = 0.0;
	double t1 = 1.0 / (16.0 *pi*(1.0 - m_DMat.v)*m_DMat.G)/R[0];
	double t2 = t1 * (3.0 - 4.0*m_DMat.v);
	double t3 = t1 / (R[0] * R[0]);
	/*
	[0 1 2
	 3 4 5
	 6 7 8]
	*/
	U[0] += t3 * R[1] * R[1] + t2;
	U[1] += t3 * R[1] * R[2];
	U[2] += t3 * R[1] * R[3];
	U[3] += U[1];
	U[4] += t3 * R[2] * R[2] + t2;
	U[5] += t3 * R[2] * R[3];
	U[6] += U[2];
	U[7] += U[5];
	U[8] += t3 * R[3] * R[3] + t2;

	return 1;
}
//静力学奇异积分UIJ,动力学无需使用
int DSquareElement::GetStaticUij(double* U, double* R)
{
	U[0] = 0.0;
	U[1] = 0.0;
	U[2] = 0.0;
	U[3] = 0.0;
	U[4] = 0.0;
	U[5] = 0.0;
	U[6] = 0.0;
	U[7] = 0.0;
	U[8] = 0.0;
	double t1 = 1.0 / (16.0 *pi*(1.0-m_DMat.v)*m_DMat.G)/R[0];
	double t2 = t1*(3.0 - 4.0*m_DMat.v);
	double t3 = t1/(R[0] * R[0]);
	/*
	[0 1 2
	 3 4 5
	 6 7 8]
	*/
	U[0] = t3 * R[1] * R[1] + t2;
	U[1] = t3 * R[1] * R[2];
	U[2] = t3 * R[1] * R[3];
	U[3] = U[1];
	U[4] = t3 * R[2] * R[2] + t2;
	U[5] = t3 * R[2] * R[3];
	U[6] = U[2];
	U[7] = U[5];
	U[8] = t3 * R[3] * R[3] + t2;

	return 1;
}

void DSquareElement::GetFijMinusOne(double* FijMinusOne, double* A, double* J, double coef1, double coef2)
{
	//coef1=-1/(8*pi*(1-v))
	//coef2=1-2*v
	double AkJk0 = A[1] * J[0] + A[2] * J[1] + A[3] * J[2];
	double t1, t2, t3, t4;
	t1 = 3.0*A[1] * AkJk0 / A[5];
	t2 = coef2 * (A[2] * J[0] - A[1] * J[1]) / A[4];
	t3 = coef2 * (A[3] * J[0] - A[1] * J[2]) / A[4];
	t4 = coef2 * (A[3] * J[1] - A[2] * J[2]) / A[4];
	double temp = coef2 * AkJk0 / A[4];
	//F(-1)11
	FijMinusOne[0] = coef1 * (temp + t1 * A[1]);
	//F(-1)21
	FijMinusOne[1] = coef1 * (t1*A[2] - t2);
	//F(-1)31
	FijMinusOne[2] = coef1 * (t1*A[3] - t3);
	//F(-1)12
	FijMinusOne[3] = coef1 * (t1*A[2] + t2);
	//F(-1)22
	FijMinusOne[4] = coef1 * (temp + 3.0*A[2] * A[2] * AkJk0 / A[5]);
	t1 = 3.0*A[3] * AkJk0 / A[5];
	//F(-1)32
	FijMinusOne[5] = coef1 * (t1*A[2] - t4);
	//F(-1)13
	FijMinusOne[6] = coef1 * (t1*A[1] + t3);
	//F(-1)23
	FijMinusOne[7] = coef1 * (t1*A[2] + t4);
	//F(-1)33
	FijMinusOne[8] = coef1 * (temp + t1 * A[3]);
}

int DSquareElement::Approximate(double* input, double* output)
{
	//input[8],output[8]
	double res[8];
	m_quadinfo.N(res, -1.0, -1.0);
	output[0] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 1.0, -1.0);
	output[1] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 1.0, 1.0);
	output[2] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, -1.0, 1.0);
	output[3] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 0.0, -1.0);
	output[4] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 1.0, 0.0);
	output[5] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 0.0, 1.0);
	output[6] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, -1.0, 0.0);
	output[7] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];

	return 1;
}

int DSquareElement::Approximate9P(double* input, double* output)
{
	//input[8],output[8]
	double res[8];
	m_quadinfo.N(res, -1.0, -1.0);
	output[0] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 1.0, -1.0);
	output[1] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 1.0, 1.0);
	output[2] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, -1.0, 1.0);
	output[3] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 0.0, -1.0);
	output[4] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 1.0, 0.0);
	output[5] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 0.0, 1.0);
	output[6] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, -1.0, 0.0);
	output[7] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];
	m_quadinfo.N(res, 0.0, 0.0);
	output[8] = res[0] * input[0] + res[1] * input[1] + res[2] * input[2] + res[3] * input[3] + res[4] * input[4] + res[5] * input[5] + res[6] * input[6] + res[7] * input[7];

	return 1;
}
//Tecplot Arrays


//calculate stress distributions
int DSquareElement::ObtainStress(double s1, double s2, double stress[3][3], double disp[8][3], double trac[8][3])
{
	// disp和trac是绝对坐标系

	double absdisp[8][3];
	double abstrac[8][3];
	double abst[3];
	double duds1[3], duds2[3];
	double dxds1[3], dxds2[3];
	double dN[8], dNds1[8], dNds2[8];
	double n[3];
	Wmatrix Y(9, 9, 0);
	Wvector P(9);
	int i;
	int r, c;

	DSquareZeroStress(stress);

	for (i = 0; i < 8; ++i)
	{
		absdisp[i][0] = disp[i][0];
		absdisp[i][1] = disp[i][1];
		absdisp[i][2] = disp[i][2];

		abstrac[i][0] = trac[i][0];
		abstrac[i][1] = trac[i][1];
		abstrac[i][2] = trac[i][2];
	}
	for (i = 0; i < 3; ++i)
	{
		abst[i] = 0.0;
		duds1[i] = 0.0;
		duds2[i] = 0.0;
	}

	m_quadinfo.N(dN, s1, s2);
	m_quadinfo.NDiffs1(dNds1, s1, s2);
	m_quadinfo.NDiffs2(dNds2, s1, s2);
	RDiffs1(s1, s2, dxds1);
	RDiffs2(s1, s2, dxds2);
	Normal(s1, s2, n[0], n[1], n[2]);
	for (i = 0; i < 8; ++i)
	{
		abst[0] += dN[i] * abstrac[i][0];
		abst[1] += dN[i] * abstrac[i][1];
		abst[2] += dN[i] * abstrac[i][2];

		duds1[0] += dNds1[i] * absdisp[i][0];
		duds1[1] += dNds1[i] * absdisp[i][1];
		duds1[2] += dNds1[i] * absdisp[i][2];

		duds2[0] += dNds2[i] * absdisp[i][0];
		duds2[1] += dNds2[i] * absdisp[i][1];
		duds2[2] += dNds2[i] * absdisp[i][2];
	}

	P[1] = abst[0];
	P[2] = abst[1];
	P[3] = abst[2];
	P[4] = duds1[0];
	P[5] = duds2[0];
	P[6] = duds1[1];
	P[7] = duds2[1];
	P[8] = duds1[2];
	P[9] = duds2[2];

	Y(1, 1) = 2.0*(1 - m_DMat.v)*m_DMat.G*n[0] / (1 - 2.0*m_DMat.v);
	Y(1, 2) = m_DMat.G*n[1];
	Y(1, 3) = m_DMat.G*n[2];
	Y(1, 4) = m_DMat.G*n[1];
	Y(1, 5) = 2.0*m_DMat.v*m_DMat.G*n[0] / (1 - 2.0*m_DMat.v);
	Y(1, 7) = m_DMat.G*n[2];
	Y(1, 9) = 2.0*m_DMat.v*m_DMat.G*n[0] / (1 - 2.0*m_DMat.v);

	Y(2, 1) = 2.0*m_DMat.v*m_DMat.G*n[1] / (1 - 2.0*m_DMat.v);
	Y(2, 2) = m_DMat.G*n[0];
	Y(2, 4) = m_DMat.G*n[0];
	Y(2, 5) = 2.0*(1 - m_DMat.v)*m_DMat.G*n[1] / (1 - 2.0*m_DMat.v);
	Y(2, 6) = m_DMat.G*n[2];
	Y(2, 8) = m_DMat.G*n[2];
	Y(2, 9) = 2.0*m_DMat.v*m_DMat.G*n[1] / (1 - 2.0*m_DMat.v);

	Y(3, 1) = 2.0*m_DMat.v*m_DMat.G*n[2] / (1 - 2.0*m_DMat.v);
	Y(3, 3) = m_DMat.G*n[0];
	Y(3, 5) = 2.0*m_DMat.v*m_DMat.G*n[2] / (1 - 2.0*m_DMat.v);
	Y(3, 6) = m_DMat.G*n[1];
	Y(3, 7) = m_DMat.G*n[0];
	Y(3, 8) = m_DMat.G*n[1];
	Y(3, 9) = 2.0*(1 - m_DMat.v)*m_DMat.G*n[2] / (1 - 2.0*m_DMat.v);

	Y(4, 1) = dxds1[0];
	Y(4, 2) = dxds1[1];
	Y(4, 3) = dxds1[2];

	Y(5, 1) = dxds2[0];
	Y(5, 2) = dxds2[1];
	Y(5, 3) = dxds2[2];

	Y(6, 4) = dxds1[0];
	Y(6, 5) = dxds1[1];
	Y(6, 6) = dxds1[2];

	Y(7, 4) = dxds2[0];
	Y(7, 5) = dxds2[1];
	Y(7, 6) = dxds2[2];

	Y(8, 7) = dxds1[0];
	Y(8, 8) = dxds1[1];
	Y(8, 9) = dxds1[2];

	Y(9, 7) = dxds2[0];
	Y(9, 8) = dxds2[1];
	Y(9, 9) = dxds2[2];

	for (r = 1; r <= 9; ++r)
	{
		if (!DSquareFinite(P[r]))
			return 0;
		for (c = 1; c <= 9; ++c)
			if (!DSquareFinite(Y(r, c)))
				return 0;
	}

	if (!gauss_mat(Y, P, 9))
		return 0;

	for (r = 1; r <= 9; ++r)
		if (!DSquareFinite(P[r]))
			return 0;

	double ukk = P[1] + P[5] + P[9];
	stress[0][0] = 2.0*m_DMat.v*m_DMat.G / (1 - 2.0*m_DMat.v)*ukk + m_DMat.G*(P[1] + P[1]);
	stress[1][0] = m_DMat.G*(P[2] + P[4]);
	stress[2][0] = m_DMat.G*(P[3] + P[7]);
	stress[0][1] = stress[1][0];
	stress[1][1] = 2.0*m_DMat.v*m_DMat.G / (1 - 2.0*m_DMat.v)*ukk + m_DMat.G*(P[5] + P[5]);
	stress[2][1] = m_DMat.G*(P[6] + P[8]);
	stress[0][2] = stress[2][0];
	stress[1][2] = stress[2][1];
	stress[2][2] = 2.0*m_DMat.v*m_DMat.G / (1 - 2.0*m_DMat.v)*ukk + m_DMat.G*(P[9] + P[9]);

	if (!DSquareStressFinite(stress))
	{
		DSquareZeroStress(stress);
		return 0;
	}

	return 1;
}

int DSquareElement::MainStress(double stress[3][3], double& stress1, double&stress2, double& stress3, double& mises, double& tresca)
{

	/*
	// 老方法计算主应力、Mises和Tresca，可能有问题
	double Pstress[3][3];
	double stress0;
	int i,j,k;
	double J2,J3;
	double theta;
	double ts;

	stress0 = (stress[0][0]+stress[1][1]+stress[2][2])/3.0;
	Pstress[0][0] = stress[0][0] - stress0;
	Pstress[1][0] = stress[1][0];
	Pstress[2][0] = stress[2][0];
	Pstress[0][1] = stress[0][1];
	Pstress[1][1] = stress[1][1] - stress0;
	Pstress[2][1] = stress[2][1];
	Pstress[0][2] = stress[0][2];
	Pstress[1][2] = stress[1][2];
	Pstress[2][2] = stress[2][2] - stress0;

	J2 = 0.0;
	J3 = 0.0;
	for(i=0;i<3;++i)
	{
		for(j=0;j<3;++j)
		{
			J2 += Pstress[i][j]*Pstress[i][j];
			for(k=0;k<3;++k)
				J3 += Pstress[i][j]*Pstress[j][k]*Pstress[k][i];
		}
	}
	J2 = J2/2.0;
	J3 = J3/3.0;
	mises = sqrt(2.0*J2/3.0);
	theta = 1.0/3.0*acos(sqrt(2.0)*J3/mises/mises/mises);
	stress1 = stress0 + sqrt(2.0)*mises*cos(theta);
	stress2 = stress0 + sqrt(2.0)*mises*cos(theta+2.0*pi/3.0);
	stress3 = stress0 + sqrt(2.0)*mises*cos(theta-2.0*pi/3.0);

	if(stress1 < stress2)
	{
		ts = stress1;
		stress1 = stress2;
		stress2 = ts;
	}
	if(stress1 < stress3)
	{
		ts = stress1;
		stress1 = stress3;
		stress3 = ts;
	}
	if(stress2 < stress3)
	{
		ts = stress2;
		stress2 = stress3;
		stress3 = ts;
	}

	tresca = stress1 - stress3;

	*/

	// 新方法计算主应力、Mises和Tresca，源自Graphitecode

	double sum = fabs(stress[0][0]) + fabs(stress[1][1]) + fabs(stress[2][2]) + 2.0*(fabs(stress[0][1]) + fabs(stress[0][2]) + fabs(stress[1][2]));

	if (sum < ERRORMAC)
	{
		stress1 = 0.0;
		stress2 = 0.0;
		stress3 = 0.0;
		mises = 0.0;
		tresca = 0.0;

		return 1;
	}

	double x0, x1, xsig0, xshear0, xi1, xi2, xi3, xj3, xsita;
	double s1, s2, s3, s4, s5, s6, ts;
	// s1: Stresses[0][0]
	// s2: Stresses[1][1]
	// s3: Stresses[2][2]
	// s4: Stresses[0][1]
	// s5: Stresses[1][2]
	// s6: Stresses[0][2]
	s1 = stress[0][0];
	s2 = stress[1][1];
	s3 = stress[2][2];
	s4 = stress[0][1];
	s5 = stress[1][2];
	s6 = stress[0][2];

	x0 = (s1 - s2)*(s1 - s2) + (s2 - s3)*(s2 - s3) + (s1 - s3)*(s1 - s3);
	x1 = 6.0*(s4*s4 + s5 * s5 + s6 * s6);
	xsig0 = (s1 + s2 + s3) / 3.0;
	xshear0 = sqrt(fabs(x0 + x1)) / 3.0;

	xi1 = xsig0 * 3.0;
	xi2 = s1 * s2 + s2 * s3 + s3 * s1 - s4 * s4 - s5 * s5 - s6 * s6;
	xi3 = s1 * s2*s3 + 2.0*s4*s5*s6 - s1 * s5*s5 - s2 * s6*s6 - s3 * s4*s4;
	xj3 = xi3 - xi1 * xi2 / 3.0 + xi1 * xi1*xi1 * 2 / 27.0;

	double over = xj3 / (xshear0*xshear0*xshear0)*sqrt(2.0);
	if (over > 1.0)
		over = 1.0;
	else if (over < -1.0)
		over = -1.0;

	xsita = acos(over) / 3.0;

	stress1 = xsig0 + xshear0 * cos(xsita)*sqrt(2.0);
	stress2 = xsig0 + xshear0 * cos(xsita + pi * 2 / 3.0)*sqrt(2.0);
	stress3 = xsig0 + xshear0 * cos(xsita - pi * 2 / 3.0)*sqrt(2.0);

	if (stress1 < stress2)
	{
		ts = stress1;
		stress1 = stress2;
		stress2 = ts;
	}
	if (stress1 < stress3)
	{
		ts = stress1;
		stress1 = stress3;
		stress3 = ts;
	}
	if (stress2 < stress3)
	{
		ts = stress2;
		stress2 = stress3;
		stress3 = ts;
	}

	mises = sqrt(3.0 / 2.0*fabs(stress1*stress1 + stress2 * stress2 + stress3 * stress3 - xi1 * xi1 / 3.0));
	tresca = stress1 - stress3;



	return 1;
}

// New Definations

int DSquareElement::GetCoordCoef()
{
	// coordinate coef arrays
	long i, j;
	double coord[3][8];
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 8; ++j)
			coord[i][j] = m_pointlist[m_pointID[j]].pt[i];
		GetCoef(m_CoordCoef[i], coord[i]);
	}

	return 1;
}

int DSquareElement::IsIn(long NodeID, int& Pos)
{
	int i, flag;
	flag = 0;
	Pos = -1;
	for (i = 0; i < 1; ++i)
	{
		if (m_nodeID[i] == NodeID)
		{
			flag = 1;
			Pos = i;
			break;
		}
	}
	return flag;
}

int DSquareElement::FuncI(double n, double dt, double* R, double& ResIn)
{
	double ndt = n * dt;

	if (ndt < R[1] / m_DMat.C1)
	{
		ResIn = 0.0;
		return 0;
	}
	else if (ndt >= R[1] / m_DMat.C1 && ndt <= R[1] / m_DMat.C2)
	{
		ResIn = 0.5 * (ndt * ndt / R[2] - m_DMat.OneOverC1_2);
		return 1;
	}
	else
	{
		// ResIn = 0.5 * (m_DMat.OneOverC2_2 - m_DMat.OneOverC1_2);
		ResIn = m_DMat.Half_1OC22_1OC12;
		return 2;
	}

	return 1;
}

int DSquareElement::FuncJ(double n, double dt, double* R, double& ResJn)
{
	double ndt = n * dt;

	if (ndt < R[1] / m_DMat.C1)
	{
		ResJn = 0.0;
		return 0;
	}
	else if (ndt >= R[1] / m_DMat.C1 && ndt <= R[1] / m_DMat.C2)
	{
		ResJn = 0.333333333333 * (ndt * ndt * ndt / R[3] - m_DMat.OneOverC1_3);
		return 1;
	}
	else
	{
		// ResJn = 1.0 / 3.0 * (m_DMat.OneOverC2_3 - m_DMat.OneOverC1_3);
		ResJn = m_DMat.OneOverThree_1OC23_1OC13;
		return 2;
	}

	return 1;
}

int DSquareElement::ISeries(double n, double dt, double* R, double& ResISeries)
{
	double ndt = n * dt;
	double np1dt = ndt + dt;

	double ROverC1 = R[1] / m_DMat.C1;
	double ROverC2 = R[1] / m_DMat.C2;

	double t1;

	if (np1dt <= ROverC1 || ndt >= ROverC2)
	{
		ResISeries = 0.0;
		return 0;
	}
	else if (ndt >= ROverC1 && np1dt <= ROverC2)
	{
		ResISeries = (n + 0.5) * m_DMat.Dt_2 / R[2];
		return 1;
	}
	else if (ndt <= ROverC1 && (np1dt >= ROverC1 && np1dt <= ROverC2))
	{
		t1 = np1dt / R[1];
		ResISeries = 0.5 * (t1 * t1 - m_DMat.OneOverC1_2);
		return 1;
	}
	else if (ndt <= ROverC1 && np1dt >= ROverC2)
	{
		ResISeries = m_DMat.Half_1OC22_1OC12;
		return 1;
	}
	else if (np1dt >= ROverC2 && (ndt >= ROverC1 && ndt <= ROverC2))
	{
		t1 = ndt / R[1];
		ResISeries = 0.5 * (m_DMat.OneOverC2_2 - t1 * t1);
		return 1;
	}

	ResISeries = 0.0;
	return 0;
}

int DSquareElement::JSeries(double n, double dt, double* R, double& ResJSeries)
{
	double ndt = n * dt;
	double np1dt = ndt + dt;

	double ROverC1 = R[1] / m_DMat.C1;
	double ROverC2 = R[1] / m_DMat.C2;

	double t1;

	if (np1dt <= ROverC1 || ndt >= ROverC2)
	{
		ResJSeries = 0.0;
		return 0;
	}
	else if (ndt >= ROverC1 && np1dt <= ROverC2)
	{
		t1 = m_DMat.Dt_3 / R[3];
		ResJSeries = (n * n + n + 0.33333333333333333) * t1;
		return 1;
	}
	else if (ndt <= ROverC1 && (np1dt >= ROverC1 && np1dt <= ROverC2))
	{
		t1 = np1dt / R[1];
		ResJSeries = 0.33333333333333 * (t1 * t1 * t1 - m_DMat.OneOverC1_3);
		return 1;
	}
	else if (ndt <= ROverC1 && np1dt >= ROverC2)
	{
		ResJSeries = m_DMat.OneOverThree_1OC23_1OC13;
		return 1;
	}
	else if (np1dt >= ROverC2 && (ndt >= ROverC1 && ndt <= ROverC2))
	{
		t1 = ndt / R[1];
		ResJSeries = 0.333333333333333 * (m_DMat.OneOverC2_3 - t1 * t1 * t1);
		return 1;
	}

	ResJSeries = 0.0;
	return 0;

}

int DSquareElement::NeedCal(int type, double n, Point& o_Point)
{
	// 该函数用来判断源点与场点单元的距离是否满足计算动态系数的条件
	// 0: 没有收到波动影响，不参与计算
	// 1: 部分收到波动影响，调用PW计算
	// 2: 完全收到波动影响，调用常规计算

	// new method

	double nn = ((type == 2) ? n - 1 : n);

	double ndt = nn * m_DMat.Dt;
	double np1dt = (nn + 1) * m_DMat.Dt;
	double r = Dist(m_center, o_Point);

	double Rmax = r + m_radius;
	double Rmin = r - m_radius;

	if (Rmin < 0.0)
		Rmin = 0.0;

	double RmaxOverC1 = Rmax / m_DMat.C1;
	double RmaxOverC2 = Rmax / m_DMat.C2;

	double RminOverC1 = Rmin / m_DMat.C1;
	double RminOverC2 = Rmin / m_DMat.C2;

	if (ndt >= RmaxOverC2 || np1dt <= RminOverC1)
		return 0;
	else if (ndt <= RminOverC1 && (np1dt >= RmaxOverC1 && np1dt <= RminOverC2))
		return 2;
	else if (ndt <= RminOverC1 && np1dt >= RmaxOverC2)
		return 2;
	else if (ndt >= RmaxOverC1 && np1dt <= RminOverC2)
		return 2;
	else if ((ndt >= RmaxOverC1 && ndt <= RminOverC2) && np1dt >= RmaxOverC2)
		return 2;
	else
		return 1;

	// old method
	/*
	double ndt = n * m_DMat.Dt;
	double np1dt = ndt + m_DMat.Dt;

	double r;
	double ROverC1,ROverC2;

	int i,flag,count;

	count = 0;

	for(i=0;i<NUMPW2;++i)
	{
		flag = 0;

		r = Dist(m_rcenPW[i],o_Point);
		ROverC1 = r / m_DMat.C1;
		ROverC2 = r / m_DMat.C2;

		if(ROverC1 >= ndt && ROverC1 < np1dt)
			flag = 1;

		if(ROverC2 >= ndt && ROverC2 < np1dt)
			flag = 1;

		if(!(np1dt <= ROverC1 || ndt >= ROverC2))
			flag = 1;

		if(flag)
			count ++;
	}

	if(count == 0)
		return 0;
	else if(count < NUMPW2)
		return 1;
	else
		return 2;
	*/
}

int DSquareElement::GetDynaUij(double* U, double n, double dt, double* R, double* RI)
{
	// n = 0..M-1

	U[0] = 0.0;
	U[1] = 0.0;
	U[2] = 0.0;
	U[3] = 0.0;
	U[4] = 0.0;
	U[5] = 0.0;
	U[6] = 0.0;
	U[7] = 0.0;
	U[8] = 0.0;

	double Tij[9], ResISeries;
	double t1;

	double RIRJ[6];
	int Flag_RIRJ = 0;

	int Flag = 0;

	if (Fai(n, m_DMat.Dt, R[1], m_DMat.C1))
	{
		Flag++;

		Flag_RIRJ++;
		RIRJ[0] = RI[1] * RI[1];
		RIRJ[1] = RI[1] * RI[2];
		RIRJ[2] = RI[1] * RI[3];
		RIRJ[3] = RI[2] * RI[2];
		RIRJ[4] = RI[2] * RI[3];
		RIRJ[5] = RI[3] * RI[3];

		GetAij(Tij, R, RI, RIRJ);
		t1 = m_DMat.OneOverFourPaiRou;
		U[0] += Tij[0] * t1;
		U[1] += Tij[1] * t1;
		U[2] += Tij[2] * t1;
		U[3] += Tij[3] * t1;
		U[4] += Tij[4] * t1;
		U[5] += Tij[5] * t1;
		U[6] += Tij[6] * t1;
		U[7] += Tij[7] * t1;
		U[8] += Tij[8] * t1;
	}

	if (Fai(n, m_DMat.Dt, R[1], m_DMat.C2))
	{
		Flag++;

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetBij(Tij, R, RI, RIRJ);
		t1 = m_DMat.OneOverFourPaiRou;
		U[0] += Tij[0] * t1;
		U[1] += Tij[1] * t1;
		U[2] += Tij[2] * t1;
		U[3] += Tij[3] * t1;
		U[4] += Tij[4] * t1;
		U[5] += Tij[5] * t1;
		U[6] += Tij[6] * t1;
		U[7] += Tij[7] * t1;
		U[8] += Tij[8] * t1;
	}

	if (ISeries(n, m_DMat.Dt, R, ResISeries))
	{
		Flag++;

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetCij(Tij, R, RI, RIRJ);
		t1 = m_DMat.OneOverFourPaiRou * ResISeries;
		U[0] += Tij[0] * t1;
		U[1] += Tij[1] * t1;
		U[2] += Tij[2] * t1;
		U[3] += Tij[3] * t1;
		U[4] += Tij[4] * t1;
		U[5] += Tij[5] * t1;
		U[6] += Tij[6] * t1;
		U[7] += Tij[7] * t1;
		U[8] += Tij[8] * t1;
	}

	return Flag;
}

int DSquareElement::GetDynaT1ij(double* T, double n, double dt, double* R, double* RI, double* Nor)
{
	// n = 0..M-1

	T[0] = 0.0;
	T[1] = 0.0;
	T[2] = 0.0;
	T[3] = 0.0;
	T[4] = 0.0;
	T[5] = 0.0;
	T[6] = 0.0;
	T[7] = 0.0;
	T[8] = 0.0;

	double T1ij[9], T2ij[9], ResISeries, ResJSeries;
	double t1, t2;

	double RIRJ[6];
	int Flag_RIRJ = 0;

	int Flag = 0;

	if (Fai(n, m_DMat.Dt, R[1], m_DMat.C1))
	{
		Flag++;

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetDij(T1ij, R, RI, Nor, RIRJ);
		GetGij(T2ij, R, RI, Nor, RIRJ);

		t1 = ONE_OVER_4PI * (n + 1.0 - R[1] / m_DMat.C1Dt);
		//		t2 = - ONE_OVER_4PI / m_DMat.Dt;   // 我推导
		t2 = ONE_OVER_4PI / m_DMat.Dt;     // 英文书结果
		T[0] += T1ij[0] * t1 + T2ij[0] * t2;
		T[1] += T1ij[1] * t1 + T2ij[1] * t2;
		T[2] += T1ij[2] * t1 + T2ij[2] * t2;
		T[3] += T1ij[3] * t1 + T2ij[3] * t2;
		T[4] += T1ij[4] * t1 + T2ij[4] * t2;
		T[5] += T1ij[5] * t1 + T2ij[5] * t2;
		T[6] += T1ij[6] * t1 + T2ij[6] * t2;
		T[7] += T1ij[7] * t1 + T2ij[7] * t2;
		T[8] += T1ij[8] * t1 + T2ij[8] * t2;
	}

	if (Fai(n, m_DMat.Dt, R[1], m_DMat.C2))
	{
		Flag++;

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetEij(T1ij, R, RI, Nor, RIRJ);
		GetHij(T2ij, R, RI, Nor, RIRJ);

		t1 = ONE_OVER_4PI * (n + 1.0 - R[1] / m_DMat.C2Dt);
		//		t2 = - ONE_OVER_4PI / m_DMat.Dt;   // 我推导
		t2 = ONE_OVER_4PI / m_DMat.Dt;   // 英文书结果
		T[0] += T1ij[0] * t1 + T2ij[0] * t2;
		T[1] += T1ij[1] * t1 + T2ij[1] * t2;
		T[2] += T1ij[2] * t1 + T2ij[2] * t2;
		T[3] += T1ij[3] * t1 + T2ij[3] * t2;
		T[4] += T1ij[4] * t1 + T2ij[4] * t2;
		T[5] += T1ij[5] * t1 + T2ij[5] * t2;
		T[6] += T1ij[6] * t1 + T2ij[6] * t2;
		T[7] += T1ij[7] * t1 + T2ij[7] * t2;
		T[8] += T1ij[8] * t1 + T2ij[8] * t2;
	}

	if (ISeries(n, m_DMat.Dt, R, ResISeries))
	{
		Flag++;

		JSeries(n, m_DMat.Dt, R, ResJSeries);

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetFij(T1ij, R, RI, Nor, RIRJ);

		t1 = ONE_OVER_4PI * ((n + 1.0) * ResISeries - R[1] * ResJSeries / m_DMat.Dt);
		T[0] += T1ij[0] * t1;
		T[1] += T1ij[1] * t1;
		T[2] += T1ij[2] * t1;
		T[3] += T1ij[3] * t1;
		T[4] += T1ij[4] * t1;
		T[5] += T1ij[5] * t1;
		T[6] += T1ij[6] * t1;
		T[7] += T1ij[7] * t1;
		T[8] += T1ij[8] * t1;
	}

	return Flag;
}

int DSquareElement::GetDynaT2ij(double* T, double n, double dt, double* R, double* RI, double* Nor)
{
	// n = 1..M

	T[0] = 0.0;
	T[1] = 0.0;
	T[2] = 0.0;
	T[3] = 0.0;
	T[4] = 0.0;
	T[5] = 0.0;
	T[6] = 0.0;
	T[7] = 0.0;
	T[8] = 0.0;

	double T1ij[9], T2ij[9], ResISeries, ResJSeries;
	double t1, t2;

	double RIRJ[6];
	int Flag_RIRJ = 0;

	int Flag = 0;

	if (Fai(n - 1, m_DMat.Dt, R[1], m_DMat.C1))
	{
		Flag++;

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetDij(T1ij, R, RI, Nor, RIRJ);
		GetGij(T2ij, R, RI, Nor, RIRJ);

		t1 = ONE_OVER_4PI * (1.0 - n + R[1] / m_DMat.C1Dt);
		//		t2 = ONE_OVER_4PI / m_DMat.Dt;    // 我推导
		t2 = -ONE_OVER_4PI / m_DMat.Dt;    // 英文书结果
		T[0] += T1ij[0] * t1 + T2ij[0] * t2;
		T[1] += T1ij[1] * t1 + T2ij[1] * t2;
		T[2] += T1ij[2] * t1 + T2ij[2] * t2;
		T[3] += T1ij[3] * t1 + T2ij[3] * t2;
		T[4] += T1ij[4] * t1 + T2ij[4] * t2;
		T[5] += T1ij[5] * t1 + T2ij[5] * t2;
		T[6] += T1ij[6] * t1 + T2ij[6] * t2;
		T[7] += T1ij[7] * t1 + T2ij[7] * t2;
		T[8] += T1ij[8] * t1 + T2ij[8] * t2;
	}

	if (Fai(n - 1, m_DMat.Dt, R[1], m_DMat.C2))
	{
		Flag++;

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetEij(T1ij, R, RI, Nor, RIRJ);
		GetHij(T2ij, R, RI, Nor, RIRJ);

		t1 = ONE_OVER_4PI * (1.0 - n + R[1] / m_DMat.C2Dt);
		//		t2 = ONE_OVER_4PI / m_DMat.Dt;    // 我推导
		t2 = -ONE_OVER_4PI / m_DMat.Dt;    // 英文书结果
		T[0] += T1ij[0] * t1 + T2ij[0] * t2;
		T[1] += T1ij[1] * t1 + T2ij[1] * t2;
		T[2] += T1ij[2] * t1 + T2ij[2] * t2;
		T[3] += T1ij[3] * t1 + T2ij[3] * t2;
		T[4] += T1ij[4] * t1 + T2ij[4] * t2;
		T[5] += T1ij[5] * t1 + T2ij[5] * t2;
		T[6] += T1ij[6] * t1 + T2ij[6] * t2;
		T[7] += T1ij[7] * t1 + T2ij[7] * t2;
		T[8] += T1ij[8] * t1 + T2ij[8] * t2;
	}

	if (ISeries(n - 1, m_DMat.Dt, R, ResISeries))
	{
		Flag++;

		JSeries(n - 1, m_DMat.Dt, R, ResJSeries);

		if (!Flag_RIRJ)
		{
			Flag_RIRJ++;
			RIRJ[0] = RI[1] * RI[1];
			RIRJ[1] = RI[1] * RI[2];
			RIRJ[2] = RI[1] * RI[3];
			RIRJ[3] = RI[2] * RI[2];
			RIRJ[4] = RI[2] * RI[3];
			RIRJ[5] = RI[3] * RI[3];
		}

		GetFij(T1ij, R, RI, Nor, RIRJ);

		t1 = ONE_OVER_4PI * ((1.0 - n) * ResISeries + R[1] * ResJSeries / m_DMat.Dt);
		T[0] += T1ij[0] * t1;
		T[1] += T1ij[1] * t1;
		T[2] += T1ij[2] * t1;
		T[3] += T1ij[3] * t1;
		T[4] += T1ij[4] * t1;
		T[5] += T1ij[5] * t1;
		T[6] += T1ij[6] * t1;
		T[7] += T1ij[7] * t1;
		T[8] += T1ij[8] * t1;
	}

	return Flag;
}

int DSquareElement::GetDynaTij(int typeT, double* T, double n, double dt, double* R, double* RI, double* Nor)
{
	if (typeT == 1)
		return GetDynaT1ij(T, n, dt, R, RI, Nor);
	else
		return GetDynaT2ij(T, n, dt, R, RI, Nor);
}

int DSquareElement::IntUij(Point& source)
{
	// n = 0..M-1

	double RI[4], UT[9], temp=0;
	int i,  PointID;

	int Flag = 0;

	for (i = 0; i < 8; ++i)
	{
		m_AssistUij[i][0] = 0.0;
		m_AssistUij[i][1] = 0.0;
		m_AssistUij[i][2] = 0.0;
		m_AssistUij[i][4] = 0.0;
		m_AssistUij[i][5] = 0.0;
		m_AssistUij[i][8] = 0.0;
	}

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], RI);
		GetStaticUij(UT, RI);

		for (PointID = 0; PointID < 8; ++PointID)
		{
			//temp = m_quadinfo.m_NRGV[PointID][i] * m_Jacobi[i];
			m_AssistUij[PointID][0] += UT[0] * temp;
			m_AssistUij[PointID][1] += UT[1] * temp;
			m_AssistUij[PointID][2] += UT[2] * temp;
			m_AssistUij[PointID][4] += UT[4] * temp;
			m_AssistUij[PointID][5] += UT[5] * temp;
			m_AssistUij[PointID][8] += UT[8] * temp;
		}
	}

	for (PointID = 0; PointID < 8; ++PointID)
	{
		m_AssistUij[PointID][3] = m_AssistUij[PointID][1];
		m_AssistUij[PointID][6] = m_AssistUij[PointID][2];
		m_AssistUij[PointID][7] = m_AssistUij[PointID][5];
	}

	return 1;
}

int DSquareElement::IntTij(Point& source)
{
	// n = 0..M-1

	double  RI[4], UT[9], temp = 0;
	int i,  PointID;

	int Flag = 0;

	for (i = 0; i < 8; ++i)
	{
		m_AssistTij[i][0] = 0.0;
		m_AssistTij[i][1] = 0.0;
		m_AssistTij[i][2] = 0.0;
		m_AssistTij[i][3] = 0.0;
		m_AssistTij[i][4] = 0.0;
		m_AssistTij[i][5] = 0.0;
		m_AssistTij[i][6] = 0.0;
		m_AssistTij[i][7] = 0.0;
		m_AssistTij[i][8] = 0.0;
	}

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], RI);
		//GetTij(UT, RI, m_normal[i]);

		for (PointID = 0; PointID < 8; ++PointID)
		{
			//temp = m_quadinfo.m_NRGV[PointID][i] * m_Jacobi[i];
			m_AssistTij[PointID][0] += UT[0] * temp;
			m_AssistTij[PointID][1] += UT[1] * temp;
			m_AssistTij[PointID][2] += UT[2] * temp;
			m_AssistTij[PointID][3] += UT[3] * temp;
			m_AssistTij[PointID][4] += UT[4] * temp;
			m_AssistTij[PointID][5] += UT[5] * temp;
			m_AssistTij[PointID][6] += UT[6] * temp;
			m_AssistTij[PointID][7] += UT[7] * temp;
			m_AssistTij[PointID][8] += UT[8] * temp;
		}
	}

	return 1;
}

int DSquareElement::IntUijforInfEle()
{
	int i;

	for (i = 0; i < 8; ++i)
	{
		m_AssistUij[i][0] = 0.0;
		m_AssistUij[i][1] = 0.0;
		m_AssistUij[i][2] = 0.0;
		m_AssistUij[i][3] = 0.0;
		m_AssistUij[i][4] = 0.0;
		m_AssistUij[i][5] = 0.0;
		m_AssistUij[i][6] = 0.0;
		m_AssistUij[i][7] = 0.0;
		m_AssistUij[i][8] = 0.0;
	}

	return 1;
}

int DSquareElement::IntTijforInfEle(Point& source,int type)
{
	double  RI[4], UT[9], temp,arpha,beta,gamma;
	double ru[3], rv[3],drda[3],drdb[3],E,F,G,Jacobi_Inf;
	double res[8],m_normalforinf[3];
	int i, j, PointID;
	Point R_bd,R_gp;
	for (i = 0; i < 8; ++i)
	{
		m_AssistTij[i][0] = 0.0;
		m_AssistTij[i][1] = 0.0;
		m_AssistTij[i][2] = 0.0;
		m_AssistTij[i][3] = 0.0;
		m_AssistTij[i][4] = 0.0;
		m_AssistTij[i][5] = 0.0;
		m_AssistTij[i][6] = 0.0;
		m_AssistTij[i][7] = 0.0;
		m_AssistTij[i][8] = 0.0;
	}

	if (type == 1)
	{
		for (i = 0; i < GAUSSPOINT2; ++i)
		{
			arpha = m_gnm.rp[i][1];
			gamma = m_gnm.rp[i][2];
			beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
			//计算R
			GetFromLocal(-arpha, -1.0, R_bd);//局部坐标
			R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
			R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
			R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
			GetR(source, R_gp, RI);
			////以交界的n为无限域n
			Normal(-arpha, -1.0, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
			GetTij(UT, RI, m_normalforinf);
			//Jacobi
			RDiffs1(-arpha, -1.0, ru);//局部坐标
			RDiffs2(-arpha, -1.0, rv);//局部坐标
			for (j = 0; j < 3; j++) {
				drda[j] = -(3.0 + beta) / 2.0*ru[j];//偏导对象ru rv (-1)
				drdb[j] = 1.0 / 2.0*(R_bd.pt[j] - DSquareElement::InfCenter[j]);
			}
			E = drda[0] * drda[0] + drda[1] * drda[1] + drda[2] * drda[2];
			F = drda[0] * drdb[0] + drda[1] * drdb[1] + drda[2] * drdb[2];
			G = drdb[0] * drdb[0] + drdb[1] * drdb[1] + drdb[2] * drdb[2];
			Jacobi_Inf = sqrt(fabs(E*G - F * F));
			//计算物理形函数
			DSquareElement::m_quadinfo.N(res, -arpha, -1.0);//局部坐标
			//计算积分
			for (PointID = 0; PointID < 8; ++PointID)
			{
				temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i];
				m_AssistTij[PointID][0] += UT[0] * temp;
				m_AssistTij[PointID][1] += UT[1] * temp;
				m_AssistTij[PointID][2] += UT[2] * temp;
				m_AssistTij[PointID][3] += UT[3] * temp;
				m_AssistTij[PointID][4] += UT[4] * temp;
				m_AssistTij[PointID][5] += UT[5] * temp;
				m_AssistTij[PointID][6] += UT[6] * temp;
				m_AssistTij[PointID][7] += UT[7] * temp;
				m_AssistTij[PointID][8] += UT[8] * temp;
			}
		}
	}
	else if (type == 2) 
	{
		for (i = 0; i < GAUSSPOINT2; ++i)
		{
			arpha = m_gnm.rp[i][1];
			gamma = m_gnm.rp[i][2];
			beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
			//计算R
			GetFromLocal(1.0, -arpha, R_bd);//局部坐标
			R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
			R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
			R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
			GetR(source, R_gp, RI);
			////以交界的n为无限域n
			Normal(1.0, -arpha, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
			GetTij(UT, RI, m_normalforinf);
			//Jacobi
			RDiffs1(1.0, -arpha, ru);//局部坐标
			RDiffs2(1.0, -arpha, rv);//局部坐标
			for (j = 0; j < 3; j++) {
				drda[j] = -(3.0 + beta) / 2.0*rv[j];//偏导对象ru rv (-1)
				drdb[j] = 1.0 / 2.0*(R_bd.pt[j] - DSquareElement::InfCenter[j]);
			}
			E = drda[0] * drda[0] + drda[1] * drda[1] + drda[2] * drda[2];
			F = drda[0] * drdb[0] + drda[1] * drdb[1] + drda[2] * drdb[2];
			G = drdb[0] * drdb[0] + drdb[1] * drdb[1] + drdb[2] * drdb[2];
			Jacobi_Inf = sqrt(fabs(E*G - F * F));
			//计算形函数
			DSquareElement::m_quadinfo.N(res, 1.0, -arpha);//局部坐标
			//计算积分
			for (PointID = 0; PointID < 8; ++PointID)
			{
				temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i];
				m_AssistTij[PointID][0] += UT[0] * temp;
				m_AssistTij[PointID][1] += UT[1] * temp;
				m_AssistTij[PointID][2] += UT[2] * temp;
				m_AssistTij[PointID][3] += UT[3] * temp;
				m_AssistTij[PointID][4] += UT[4] * temp;
				m_AssistTij[PointID][5] += UT[5] * temp;
				m_AssistTij[PointID][6] += UT[6] * temp;
				m_AssistTij[PointID][7] += UT[7] * temp;
				m_AssistTij[PointID][8] += UT[8] * temp;
			}
		}
	}
	else if (type == 3)
	{
		for (i = 0; i < GAUSSPOINT2; ++i)
		{
			arpha = m_gnm.rp[i][1];
			gamma = m_gnm.rp[i][2];
			beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
			//计算R
			GetFromLocal(arpha, 1.0, R_bd);//局部坐标
			R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
			R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
			R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
			GetR(source, R_gp, RI);
			////以交界的n为无限域n
			Normal(arpha, 1.0, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
			GetTij(UT, RI, m_normalforinf);
			//Jacobi
			RDiffs1(arpha, 1.0, ru);//局部坐标
			RDiffs2(arpha, 1.0, rv);//局部坐标
			for (j = 0; j < 3; j++) {
				drda[j] = (3.0 + beta) / 2.0*ru[j];//偏导对象ru rv
				drdb[j] = 1.0 / 2.0*(R_bd.pt[j] - DSquareElement::InfCenter[j]);
			}
			E = drda[0] * drda[0] + drda[1] * drda[1] + drda[2] * drda[2];
			F = drda[0] * drdb[0] + drda[1] * drdb[1] + drda[2] * drdb[2];
			G = drdb[0] * drdb[0] + drdb[1] * drdb[1] + drdb[2] * drdb[2];
			Jacobi_Inf = sqrt(fabs(E*G - F * F));
			//计算形函数
			DSquareElement::m_quadinfo.N(res, arpha, 1.0);//局部坐标
			//计算积分
			for (PointID = 0; PointID < 8; ++PointID)
			{
				temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i];
				m_AssistTij[PointID][0] += UT[0] * temp;
				m_AssistTij[PointID][1] += UT[1] * temp;
				m_AssistTij[PointID][2] += UT[2] * temp;
				m_AssistTij[PointID][3] += UT[3] * temp;
				m_AssistTij[PointID][4] += UT[4] * temp;
				m_AssistTij[PointID][5] += UT[5] * temp;
				m_AssistTij[PointID][6] += UT[6] * temp;
				m_AssistTij[PointID][7] += UT[7] * temp;
				m_AssistTij[PointID][8] += UT[8] * temp;
			}
		}
	}
	else if (type == 4)
	{
		for (i = 0; i < GAUSSPOINT2; ++i)
		{
			arpha = m_gnm.rp[i][1];
			gamma = m_gnm.rp[i][2];
			beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
			//计算R
			GetFromLocal(-1.0, arpha, R_bd);//局部坐标
			R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
			R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
			R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
			GetR(source, R_gp, RI);
			////以交界的n为无限域n
			Normal(-1.0, arpha, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
			GetTij(UT, RI, m_normalforinf);
			//Jacobi
			RDiffs1(-1.0, arpha, ru);//局部坐标
			RDiffs2(-1.0, arpha, rv);//局部坐标
			for (j = 0; j < 3; j++) {
				drda[j] = (3.0 + beta) / 2.0*rv[j];//偏导对象ru rv
				drdb[j] = 1.0 / 2.0*(R_bd.pt[j] - DSquareElement::InfCenter[j]);
			}
			E = drda[0] * drda[0] + drda[1] * drda[1] + drda[2] * drda[2];
			F = drda[0] * drdb[0] + drda[1] * drdb[1] + drda[2] * drdb[2];
			G = drdb[0] * drdb[0] + drdb[1] * drdb[1] + drdb[2] * drdb[2];
			Jacobi_Inf = sqrt(fabs(E*G - F * F));
			//计算形函数
			DSquareElement::m_quadinfo.N(res, -1.0, arpha);//局部坐标
			//计算积分
			for (PointID = 0; PointID < 8; ++PointID)
			{
				temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i];
				m_AssistTij[PointID][0] += UT[0] * temp;
				m_AssistTij[PointID][1] += UT[1] * temp;
				m_AssistTij[PointID][2] += UT[2] * temp;
				m_AssistTij[PointID][3] += UT[3] * temp;
				m_AssistTij[PointID][4] += UT[4] * temp;
				m_AssistTij[PointID][5] += UT[5] * temp;
				m_AssistTij[PointID][6] += UT[6] * temp;
				m_AssistTij[PointID][7] += UT[7] * temp;
				m_AssistTij[PointID][8] += UT[8] * temp;
			}
		}
	}
	return 1;
}

int DSquareElement::IntDynaUij(Point& source, double n, double dt)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp = 0;
	int i,  PointID;

	int Flag = 0;

	for (i = 0; i < 8; ++i)
	{
		m_AssistUij[i][0] = 0.0;
		m_AssistUij[i][1] = 0.0;
		m_AssistUij[i][2] = 0.0;
		m_AssistUij[i][4] = 0.0;
		m_AssistUij[i][5] = 0.0;
		m_AssistUij[i][8] = 0.0;
	}

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], R, RI);
		Flag = GetDynaUij(UT, n, dt, R, RI);

		if (Flag)
		{
			for (PointID = 0; PointID < 8; ++PointID)
			{
				//temp = m_quadinfo.m_NRGV[PointID][i] * m_Jacobi[i];
				m_AssistUij[PointID][0] += UT[0] * temp;
				m_AssistUij[PointID][1] += UT[1] * temp;
				m_AssistUij[PointID][2] += UT[2] * temp;
				m_AssistUij[PointID][4] += UT[4] * temp;
				m_AssistUij[PointID][5] += UT[5] * temp;
				m_AssistUij[PointID][8] += UT[8] * temp;
			}
		}
	}

	for (PointID = 0; PointID < 8; ++PointID)
	{
		m_AssistUij[PointID][3] = m_AssistUij[PointID][1];
		m_AssistUij[PointID][6] = m_AssistUij[PointID][2];
		m_AssistUij[PointID][7] = m_AssistUij[PointID][5];
	}

	return 1;
}

int DSquareElement::IntDynaUij(Point& source, double n, double dt, long T_ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp;
	int i, PointID;
	Point intpt;
	double JacobiT;


	int Flag = 0;

	for (i = 0; i < 8; ++i)
	{
		AssistUij[T_ID][i][0] = 0.0;
		AssistUij[T_ID][i][1] = 0.0;
		AssistUij[T_ID][i][2] = 0.0;
		AssistUij[T_ID][i][4] = 0.0;
		AssistUij[T_ID][i][5] = 0.0;
		AssistUij[T_ID][i][8] = 0.0;
	}

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		GetFromLocal(m_gnm.rp[i][0], m_gnm.rp[i][1], intpt);
		GetR(source, intpt, R, RI);
		//GetR(source, m_intpt[i], R, RI);
		Flag = GetDynaUij(UT, n, dt, R, RI);

		if (Flag)
		{
			for (PointID = 0; PointID < 8; ++PointID)
			{
				JacobiT = Jacobi(m_gnm.rp[i][0], m_gnm.rp[i][1]);
				temp = m_quadinfo.m_NRGV[PointID][i] * JacobiT;
				//temp = m_quadinfo.m_NRGV[PointID][i] * m_Jacobi[i];
				AssistUij[T_ID][PointID][0] += UT[0] * temp;
				AssistUij[T_ID][PointID][1] += UT[1] * temp;
				AssistUij[T_ID][PointID][2] += UT[2] * temp;
				AssistUij[T_ID][PointID][4] += UT[4] * temp;
				AssistUij[T_ID][PointID][5] += UT[5] * temp;
				AssistUij[T_ID][PointID][8] += UT[8] * temp;
			}
		}
	}

	for (PointID = 0; PointID < 8; ++PointID)
	{
		AssistUij[T_ID][PointID][3] = AssistUij[T_ID][PointID][1];
		AssistUij[T_ID][PointID][6] = AssistUij[T_ID][PointID][2];
		AssistUij[T_ID][PointID][7] = AssistUij[T_ID][PointID][5];
	}

	return 1;
}

int DSquareElement::IntDynaUij(Point& source, double n, double dt, int ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp = 0;
	int i;

	int Flag = 0;

	m_AssistUij[ID][0] = 0.0;
	m_AssistUij[ID][1] = 0.0;
	m_AssistUij[ID][2] = 0.0;
	m_AssistUij[ID][4] = 0.0;
	m_AssistUij[ID][5] = 0.0;
	m_AssistUij[ID][8] = 0.0;

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], R, RI);
		Flag = GetDynaUij(UT, n, dt, R, RI);

		if (Flag)
		{
			//temp = m_quadinfo.m_NRGV[ID][i] * m_Jacobi[i];
			m_AssistUij[ID][0] += UT[0] * temp;
			m_AssistUij[ID][1] += UT[1] * temp;
			m_AssistUij[ID][2] += UT[2] * temp;
			m_AssistUij[ID][4] += UT[4] * temp;
			m_AssistUij[ID][5] += UT[5] * temp;
			m_AssistUij[ID][8] += UT[8] * temp;
		}
	}

	m_AssistUij[ID][3] = m_AssistUij[ID][1];
	m_AssistUij[ID][6] = m_AssistUij[ID][2];
	m_AssistUij[ID][7] = m_AssistUij[ID][5];

	return 1;
}

int DSquareElement::IntDynaTij(int typeT, Point& source, double n, double dt)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp = 0;
	int i,  PointID;

	int Flag = 0;

	for (i = 0; i < 8; ++i)
	{
		m_AssistTij[i][0] = 0.0;
		m_AssistTij[i][1] = 0.0;
		m_AssistTij[i][2] = 0.0;
		m_AssistTij[i][3] = 0.0;
		m_AssistTij[i][4] = 0.0;
		m_AssistTij[i][5] = 0.0;
		m_AssistTij[i][6] = 0.0;
		m_AssistTij[i][7] = 0.0;
		m_AssistTij[i][8] = 0.0;
	}

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], R, RI);
		//Flag = GetDynaTij(typeT, UT, n, dt, R, RI, m_normal[i]);

		if (Flag)
		{
			for (PointID = 0; PointID < 8; ++PointID)
			{
			//	temp = m_quadinfo.m_NRGV[PointID][i] * m_Jacobi[i];
				m_AssistTij[PointID][0] += UT[0] * temp;
				m_AssistTij[PointID][1] += UT[1] * temp;
				m_AssistTij[PointID][2] += UT[2] * temp;
				m_AssistTij[PointID][3] += UT[3] * temp;
				m_AssistTij[PointID][4] += UT[4] * temp;
				m_AssistTij[PointID][5] += UT[5] * temp;
				m_AssistTij[PointID][6] += UT[6] * temp;
				m_AssistTij[PointID][7] += UT[7] * temp;
				m_AssistTij[PointID][8] += UT[8] * temp;
			}
		}
	}

	return 1;
}

int DSquareElement::IntDynaTij(int typeT, Point& source, double n, double dt, int ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp=0;
	int i;

	int Flag = 0;

	m_AssistTij[ID][0] = 0.0;
	m_AssistTij[ID][1] = 0.0;
	m_AssistTij[ID][2] = 0.0;
	m_AssistTij[ID][3] = 0.0;
	m_AssistTij[ID][4] = 0.0;
	m_AssistTij[ID][5] = 0.0;
	m_AssistTij[ID][6] = 0.0;
	m_AssistTij[ID][7] = 0.0;
	m_AssistTij[ID][8] = 0.0;

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], R, RI);
		//Flag = GetDynaTij(typeT, UT, n, dt, R, RI, m_normal[i]);

		if (Flag)
		{
			//temp = m_quadinfo.m_NRGV[ID][i] * m_Jacobi[i];
			m_AssistTij[ID][0] += UT[0] * temp;
			m_AssistTij[ID][1] += UT[1] * temp;
			m_AssistTij[ID][2] += UT[2] * temp;
			m_AssistTij[ID][3] += UT[3] * temp;
			m_AssistTij[ID][4] += UT[4] * temp;
			m_AssistTij[ID][5] += UT[5] * temp;
			m_AssistTij[ID][6] += UT[6] * temp;
			m_AssistTij[ID][7] += UT[7] * temp;
			m_AssistTij[ID][8] += UT[8] * temp;
		}
	}

	return 1;
}

int DSquareElement::IntDynaTij(Point& source, double n, double dt)
{
	// n = 0..M-1

	double R[6], RI[4], UT1[9], UT2[9], temp=0;
	int i,  PointID;

	int Flag = 0;

	for (i = 0; i < 8; ++i)
	{
		m_AssistTij[i][0] = 0.0;
		m_AssistTij[i][1] = 0.0;
		m_AssistTij[i][2] = 0.0;
		m_AssistTij[i][3] = 0.0;
		m_AssistTij[i][4] = 0.0;
		m_AssistTij[i][5] = 0.0;
		m_AssistTij[i][6] = 0.0;
		m_AssistTij[i][7] = 0.0;
		m_AssistTij[i][8] = 0.0;
	}

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], R, RI);
		//Flag = GetDynaTij(1, UT1, n, dt, R, RI, m_normal[i]);
		//Flag += GetDynaTij(2, UT2, n, dt, R, RI, m_normal[i]);

		if (Flag)
		{
			for (PointID = 0; PointID < 8; ++PointID)
			{
			//	temp = m_quadinfo.m_NRGV[PointID][i] * m_Jacobi[i];
				m_AssistTij[PointID][0] += (UT1[0] + UT2[0]) * temp;
				m_AssistTij[PointID][1] += (UT1[1] + UT2[1]) * temp;
				m_AssistTij[PointID][2] += (UT1[2] + UT2[2]) * temp;
				m_AssistTij[PointID][3] += (UT1[3] + UT2[3]) * temp;
				m_AssistTij[PointID][4] += (UT1[4] + UT2[4]) * temp;
				m_AssistTij[PointID][5] += (UT1[5] + UT2[5]) * temp;
				m_AssistTij[PointID][6] += (UT1[6] + UT2[6]) * temp;
				m_AssistTij[PointID][7] += (UT1[7] + UT2[7]) * temp;
				m_AssistTij[PointID][8] += (UT1[8] + UT2[8]) * temp;
			}
		}
	}

	return 1;
}

int DSquareElement::IntDynaTij(Point& source, double n, double dt, int ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT1[9], UT2[9], temp=0;
	int i;

	int Flag = 0;

	m_AssistTij[ID][0] = 0.0;
	m_AssistTij[ID][1] = 0.0;
	m_AssistTij[ID][2] = 0.0;
	m_AssistTij[ID][3] = 0.0;
	m_AssistTij[ID][4] = 0.0;
	m_AssistTij[ID][5] = 0.0;
	m_AssistTij[ID][6] = 0.0;
	m_AssistTij[ID][7] = 0.0;
	m_AssistTij[ID][8] = 0.0;

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		//GetR(source, m_intpt[i], R, RI);
		//Flag = GetDynaTij(1, UT1, n, dt, R, RI, m_normal[i]);
		//Flag += GetDynaTij(2, UT2, n, dt, R, RI, m_normal[i]);

		if (Flag)
		{
			//temp = m_quadinfo.m_NRGV[ID][i] * m_Jacobi[i];
			m_AssistTij[ID][0] += (UT1[0] + UT2[0]) * temp;
			m_AssistTij[ID][1] += (UT1[1] + UT2[1]) * temp;
			m_AssistTij[ID][2] += (UT1[2] + UT2[2]) * temp;
			m_AssistTij[ID][3] += (UT1[3] + UT2[3]) * temp;
			m_AssistTij[ID][4] += (UT1[4] + UT2[4]) * temp;
			m_AssistTij[ID][5] += (UT1[5] + UT2[5]) * temp;
			m_AssistTij[ID][6] += (UT1[6] + UT2[6]) * temp;
			m_AssistTij[ID][7] += (UT1[7] + UT2[7]) * temp;
			m_AssistTij[ID][8] += (UT1[8] + UT2[8]) * temp;
		}
	}

	return 1;
}

int DSquareElement::IntDynaTij(int typeT, Point& source, double n, double dt, long T_ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT1[9], UT2[9], temp;
	int i, PointID;
	Point intpt;
	double JacobiT;
	double normalT[3];

	int Flag = 0;

	for (i = 0; i < 8; ++i)
	{
		AssistTij[T_ID][i][0] = 0.0;
		AssistTij[T_ID][i][1] = 0.0;
		AssistTij[T_ID][i][2] = 0.0;
		AssistTij[T_ID][i][3] = 0.0;
		AssistTij[T_ID][i][4] = 0.0;
		AssistTij[T_ID][i][5] = 0.0;
		AssistTij[T_ID][i][6] = 0.0;
		AssistTij[T_ID][i][7] = 0.0;
		AssistTij[T_ID][i][8] = 0.0;
	}

	for (i = 0; i < GAUSSPOINT2; ++i)
	{
		GetFromLocal(m_gnm.rp[i][0], m_gnm.rp[i][1], intpt);
		GetR(source, intpt, R, RI);
		Normal(m_gnm.rp[i][0], m_gnm.rp[i][1], normalT[0], normalT[1], normalT[2]);
		Flag = GetDynaTij(1, UT1, n, dt, R, RI, normalT);
		Flag += GetDynaTij(2, UT2, n, dt, R, RI, normalT);
		//GetR(source, m_intpt[i], R, RI);
		//Flag = GetDynaTij(1, UT1, n, dt, R, RI, m_normal[i]);
		//Flag += GetDynaTij(2, UT2, n, dt, R, RI, m_normal[i]);

		if (Flag)
		{
			for (PointID = 0; PointID < 8; ++PointID)
			{
				JacobiT = Jacobi(m_gnm.rp[i][0], m_gnm.rp[i][1]);
				temp = m_quadinfo.m_NRGV[PointID][i] * JacobiT;
				//temp = m_quadinfo.m_NRGV[PointID][i] * m_Jacobi[i];
				AssistTij[T_ID][PointID][0] += (UT1[0] + UT2[0]) * temp;
				AssistTij[T_ID][PointID][1] += (UT1[1] + UT2[1]) * temp;
				AssistTij[T_ID][PointID][2] += (UT1[2] + UT2[2]) * temp;
				AssistTij[T_ID][PointID][3] += (UT1[3] + UT2[3]) * temp;
				AssistTij[T_ID][PointID][4] += (UT1[4] + UT2[4]) * temp;
				AssistTij[T_ID][PointID][5] += (UT1[5] + UT2[5]) * temp;
				AssistTij[T_ID][PointID][6] += (UT1[6] + UT2[6]) * temp;
				AssistTij[T_ID][PointID][7] += (UT1[7] + UT2[7]) * temp;
				AssistTij[T_ID][PointID][8] += (UT1[8] + UT2[8]) * temp;
			}
		}
	}

	return 1;
}

int DSquareElement::InElementStaticTij()
{
	// 计算本单元静态奇异积分Tij，存于m_StaticTij公共变量中
	double coef1, coef2, UT[9], temp, s1, s2;
	double FVS1, FVL1;
	int i, j, k, l, AreaID;

	//获取8个物理点，每个物理点对应4个三角形的共8*4*SINGULARGAUSSPOINT*SINGULARGAUSSPOINT个积分点上的J、normal、R
	for (i = 0; i < 1; ++i)
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT2; ++j)
			{
				s1 = m_SecInfo.m_SecConstPos[AreaID][0][j];
				s2 = m_SecInfo.m_SecConstPos[AreaID][1][j];
				m_OnEleJ[i][AreaID][j] = Jacobi(s1, s2);
				Normal(s1, s2, m_OnEleN[i][AreaID][j][0], m_OnEleN[i][AreaID][j][1], m_OnEleN[i][AreaID][j][2]);
				GetR(m_nodelist[m_nodeID[i]], s1, s2, m_OnEleR[i][AreaID][j]);
			}
		}

	// 计算Tij的奇异积分：对于强奇异，需要在角坐标系下加一项减一项；对于弱奇异，只需要在角坐标系下积分即可

	coef1 = -1.0 / (8.0*pi*(1.0 - m_DMat.v));
	coef2 = 1.0 - 2.0*m_DMat.v;
	for (i = 0; i < 1; ++i)
		for (j = 0; j < 8; ++j)
			for (k = 0; k < 9; ++k)
				m_StaticTij[i][j][k] = 0.0;

	// 强奇异

	//variables for assissant of FijMinusOne
	double m_FijMinusOne[9], A[6], Jh[3], DrDs1[3], DrDs2[3];
	double Beta;   // 对多个单元共享一个节点的情况要用到
	for (i = 0; i < 1; ++i)
	{
		s1 = 0.0;
		s2 = 0.0;
		GetJi0(i, Jh);
		RDiffs1(s1, s2, DrDs1);
		RDiffs2(s1, s2, DrDs2);

		//Tij
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT2; ++j)
			{					//Tij
				GetTij(UT, m_OnEleR[i][AreaID][j], m_OnEleN[i][AreaID][j]);
				temp = m_SecInfo.m_SecConstVal[AreaID][j] * m_OnEleJ[i][AreaID][j];
				for (l = 0; l < 9; ++l)
				{
					m_StaticTij[i][i][l] += UT[l] * temp;
				}
			}
		}
		//FijMinusOne
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT; ++j)
			{
				GetAi(m_SecInfo.m_ConstLinePos[AreaID][j], DrDs1, DrDs2, A);

				GetBeta(Beta, A);
				Beta = fabs(Beta);

				GetFijMinusOne(m_FijMinusOne, A, Jh, coef1, coef2);
				for (l = 0; l < SINGAUSSPOINT; ++l)
				{
					FVS1 = m_FValue.m_ConstFMinusOneSecVal[AreaID][l * SINGAUSSPOINT + j];  // [l][j] TEST
					m_StaticTij[i][i][0] -= m_FijMinusOne[0] * FVS1;
					m_StaticTij[i][i][1] -= m_FijMinusOne[1] * FVS1;
					m_StaticTij[i][i][2] -= m_FijMinusOne[2] * FVS1;
					m_StaticTij[i][i][3] -= m_FijMinusOne[3] * FVS1;
					m_StaticTij[i][i][4] -= m_FijMinusOne[4] * FVS1;
					m_StaticTij[i][i][5] -= m_FijMinusOne[5] * FVS1;
					m_StaticTij[i][i][6] -= m_FijMinusOne[6] * FVS1;
					m_StaticTij[i][i][7] -= m_FijMinusOne[7] * FVS1;
					m_StaticTij[i][i][8] -= m_FijMinusOne[8] * FVS1;
				}

				// 对单个节点由若干个单元共享的情况，需要用到Beta
				FVL1 = m_FValue.m_ConstFMinusOneLineVal[AreaID][j]; // 连续单元用： + m_FValue.m_FMinusOneLineVal_2[i][AreaID][j]*log(Beta);
				m_StaticTij[i][i][0] += m_FijMinusOne[0] * FVL1;
				m_StaticTij[i][i][1] += m_FijMinusOne[1] * FVL1;
				m_StaticTij[i][i][2] += m_FijMinusOne[2] * FVL1;
				m_StaticTij[i][i][3] += m_FijMinusOne[3] * FVL1;
				m_StaticTij[i][i][4] += m_FijMinusOne[4] * FVL1;
				m_StaticTij[i][i][5] += m_FijMinusOne[5] * FVL1;
				m_StaticTij[i][i][6] += m_FijMinusOne[6] * FVL1;
				m_StaticTij[i][i][7] += m_FijMinusOne[7] * FVL1;
				m_StaticTij[i][i][8] += m_FijMinusOne[8] * FVL1;
			}

		}
	}

	//	for(i=0;i<8;++i)
	//	{
	//		m_StaticTij[i][i][0] += 0.5;
	//		m_StaticTij[i][i][4] += 0.5;
	//		m_StaticTij[i][i][8] += 0.5;
	//
	//	}

		// 弱奇异

	for (i = 0; i < 1; ++i)// recycle for source point
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINGAUSSPOINT2; ++k)
			{
				GetTij(UT, m_OnEleR[i][AreaID][k], m_OnEleN[i][AreaID][k]);

				for (j = 0; j < 1; ++j)// recycle for field point
				{
					if (i != j)
					{
						temp = m_SecInfo.m_SecConstVal[AreaID][k] * m_OnEleJ[i][AreaID][k];
						m_StaticTij[i][j][0] += UT[0] * temp;
						m_StaticTij[i][j][1] += UT[1] * temp;
						m_StaticTij[i][j][2] += UT[2] * temp;
						m_StaticTij[i][j][3] += UT[3] * temp;
						m_StaticTij[i][j][4] += UT[4] * temp;
						m_StaticTij[i][j][5] += UT[5] * temp;
						m_StaticTij[i][j][6] += UT[6] * temp;
						m_StaticTij[i][j][7] += UT[7] * temp;
						m_StaticTij[i][j][8] += UT[8] * temp;
					}
				}
			}

		}

	return 1;
}

int DSquareElement::InElementStaticTij(double* StaticTij)
{
	// 计算本单元静态奇异积分Tij，存于m_StaticTij公共变量中
	double coef1, coef2, UT[9], temp, s1, s2;
	double FVS1, FVL1;
	int i, j, k, l, AreaID;
	double** OnEleJ;
	double*** OnEleN, ***OnEleR;
	//double OnEleJ[4][SINGAUSSPOINT2];
	//double OnEleN[4][SINGAUSSPOINT2][3];
	//double OnEleR[4][SINGAUSSPOINT2][3];

	OnEleJ = new double* [4];
	OnEleN = new double** [4];
	OnEleR = new double** [4];
	for (i = 0; i < 4; i++) {
		OnEleJ[i] = new double [SINGAUSSPOINT2];
		OnEleN[i] = new double* [SINGAUSSPOINT2];
		OnEleR[i] = new double* [SINGAUSSPOINT2];
		for (j = 0; j < SINGAUSSPOINT2; j++) {
			OnEleN[i][j] = new double [3];
			OnEleR[i][j] = new double [4];
		}
	}

	//获取8个物理点，每个物理点对应4个三角形的共8*4*SINGULARGAUSSPOINT*SINGULARGAUSSPOINT个积分点上的J、normal、R
	for (i = 0; i < 1; ++i)
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT2; ++j)
			{
				s1 = m_SecInfo.m_SecConstPos[AreaID][0][j];
				s2 = m_SecInfo.m_SecConstPos[AreaID][1][j];
				//m_OnEleJ[i][AreaID][j] = Jacobi(s1, s2);
				//Normal(s1, s2, m_OnEleN[i][AreaID][j][0], m_OnEleN[i][AreaID][j][1], m_OnEleN[i][AreaID][j][2]);
				//GetR(m_nodelist[m_nodeID[i]], s1, s2, m_OnEleR[i][AreaID][j]);
				OnEleJ[AreaID][j] = Jacobi(s1, s2);//m_OnEleJ[i][AreaID][j] = Jacobi(s1, s2);
				Normal(s1, s2, OnEleN[AreaID][j][0], OnEleN[AreaID][j][1], OnEleN[AreaID][j][2]);
				GetR(m_nodelist[m_nodeID[i]], s1, s2, OnEleR[AreaID][j]);
			}
		}

	// 计算Tij的奇异积分：对于强奇异，需要在角坐标系下加一项减一项；对于弱奇异，只需要在角坐标系下积分即可

	coef1 = -1.0 / (8.0 * pi * (1.0 - m_DMat.v));
	coef2 = 1.0 - 2.0 * m_DMat.v;
	for (i = 0; i < 1; ++i)
		for (j = 0; j < 8; ++j)
			for (k = 0; k < 9; ++k)
				StaticTij[k] = 0.0;//m_StaticTij[i][j][k] = 0.0;

	// 强奇异

	//variables for assissant of FijMinusOne
	double m_FijMinusOne[9], A[6], Jh[3], DrDs1[3], DrDs2[3];
	double Beta;   // 对多个单元共享一个节点的情况要用到
	for (i = 0; i < 1; ++i)
	{
		s1 = 0.0;
		s2 = 0.0;
		GetJi0(i, Jh);
		RDiffs1(s1, s2, DrDs1);
		RDiffs2(s1, s2, DrDs2);

		//Tij
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT2; ++j)
			{					//Tij
				//GetTij(UT, m_OnEleR[i][AreaID][j], m_OnEleN[i][AreaID][j]);
				//temp = m_SecInfo.m_SecConstVal[AreaID][j] * m_OnEleJ[i][AreaID][j];
				GetTij(UT, OnEleR[AreaID][j], OnEleN[AreaID][j]);
				temp = m_SecInfo.m_SecConstVal[AreaID][j] * OnEleJ[AreaID][j];

				for (l = 0; l < 9; ++l)
				{
					//m_StaticTij[i][i][l] += UT[l] * temp;
					StaticTij[l] += UT[l] * temp;
				}
			}
		}
		//FijMinusOne
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT; ++j)
			{
				GetAi(m_SecInfo.m_ConstLinePos[AreaID][j], DrDs1, DrDs2, A);

				GetBeta(Beta, A);
				Beta = fabs(Beta);

				GetFijMinusOne(m_FijMinusOne, A, Jh, coef1, coef2);
				for (l = 0; l < SINGAUSSPOINT; ++l)
				{
					FVS1 = m_FValue.m_ConstFMinusOneSecVal[AreaID][l * SINGAUSSPOINT + j];  // [l][j] TEST

					/*m_StaticTij[i][i][0] -= m_FijMinusOne[0] * FVS1;
					m_StaticTij[i][i][1] -= m_FijMinusOne[1] * FVS1;
					m_StaticTij[i][i][2] -= m_FijMinusOne[2] * FVS1;
					m_StaticTij[i][i][3] -= m_FijMinusOne[3] * FVS1;
					m_StaticTij[i][i][4] -= m_FijMinusOne[4] * FVS1;
					m_StaticTij[i][i][5] -= m_FijMinusOne[5] * FVS1;
					m_StaticTij[i][i][6] -= m_FijMinusOne[6] * FVS1;
					m_StaticTij[i][i][7] -= m_FijMinusOne[7] * FVS1;
					m_StaticTij[i][i][8] -= m_FijMinusOne[8] * FVS1;*/
					StaticTij[0] -= m_FijMinusOne[0] * FVS1;
					StaticTij[1] -= m_FijMinusOne[1] * FVS1;
					StaticTij[2] -= m_FijMinusOne[2] * FVS1;
					StaticTij[3] -= m_FijMinusOne[3] * FVS1;
					StaticTij[4] -= m_FijMinusOne[4] * FVS1;
					StaticTij[5] -= m_FijMinusOne[5] * FVS1;
					StaticTij[6] -= m_FijMinusOne[6] * FVS1;
					StaticTij[7] -= m_FijMinusOne[7] * FVS1;
					StaticTij[8] -= m_FijMinusOne[8] * FVS1;
				}

				// 对单个节点由若干个单元共享的情况，需要用到Beta
				FVL1 = m_FValue.m_ConstFMinusOneLineVal[AreaID][j]; // 连续单元用： + m_FValue.m_FMinusOneLineVal_2[i][AreaID][j]*log(Beta);
				//m_StaticTij[i][i][0] += m_FijMinusOne[0] * FVL1;
				//m_StaticTij[i][i][1] += m_FijMinusOne[1] * FVL1;
				//m_StaticTij[i][i][2] += m_FijMinusOne[2] * FVL1;
				//m_StaticTij[i][i][3] += m_FijMinusOne[3] * FVL1;
				//m_StaticTij[i][i][4] += m_FijMinusOne[4] * FVL1;
				//m_StaticTij[i][i][5] += m_FijMinusOne[5] * FVL1;
				//m_StaticTij[i][i][6] += m_FijMinusOne[6] * FVL1;
				//m_StaticTij[i][i][7] += m_FijMinusOne[7] * FVL1;
				//m_StaticTij[i][i][8] += m_FijMinusOne[8] * FVL1;
				StaticTij[0] += m_FijMinusOne[0] * FVL1;
				StaticTij[1] += m_FijMinusOne[1] * FVL1;
				StaticTij[2] += m_FijMinusOne[2] * FVL1;
				StaticTij[3] += m_FijMinusOne[3] * FVL1;
				StaticTij[4] += m_FijMinusOne[4] * FVL1;
				StaticTij[5] += m_FijMinusOne[5] * FVL1;
				StaticTij[6] += m_FijMinusOne[6] * FVL1;
				StaticTij[7] += m_FijMinusOne[7] * FVL1;
				StaticTij[8] += m_FijMinusOne[8] * FVL1;
			}

		}
	}

	//	for(i=0;i<8;++i)
	//	{
	//		m_StaticTij[i][i][0] += 0.5;
	//		m_StaticTij[i][i][4] += 0.5;
	//		m_StaticTij[i][i][8] += 0.5;
	//
	//	}

		// 弱奇异

	for (i = 0; i < 1; ++i)// recycle for source point
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINGAUSSPOINT2; ++k)
			{
				//GetTij(UT, m_OnEleR[i][AreaID][k], m_OnEleN[i][AreaID][k]);
				GetTij(UT, OnEleR[AreaID][k], OnEleN[AreaID][k]);

				for (j = 0; j < 1; ++j)// recycle for field point
				{
					if (i != j)
					{
						//temp = m_SecInfo.m_SecConstVal[AreaID][k] * m_OnEleJ[i][AreaID][k];
						temp = m_SecInfo.m_SecConstVal[AreaID][k] * OnEleJ[AreaID][k];
						//m_StaticTij[i][j][0] += UT[0] * temp;
						//m_StaticTij[i][j][1] += UT[1] * temp;
						//m_StaticTij[i][j][2] += UT[2] * temp;
						//m_StaticTij[i][j][3] += UT[3] * temp;
						//m_StaticTij[i][j][4] += UT[4] * temp;
						//m_StaticTij[i][j][5] += UT[5] * temp;
						//m_StaticTij[i][j][6] += UT[6] * temp;
						//m_StaticTij[i][j][7] += UT[7] * temp;
						//m_StaticTij[i][j][8] += UT[8] * temp;
						StaticTij[0] += UT[0] * temp;
						StaticTij[1] += UT[1] * temp;
						StaticTij[2] += UT[2] * temp;
						StaticTij[3] += UT[3] * temp;
						StaticTij[4] += UT[4] * temp;
						StaticTij[5] += UT[5] * temp;
						StaticTij[6] += UT[6] * temp;
						StaticTij[7] += UT[7] * temp;
						StaticTij[8] += UT[8] * temp;
					}
				}
			}

		}
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < SINGAUSSPOINT2; j++) {
			delete[] OnEleR[i][j];
			OnEleR[i][j] = NULL;
			delete[] OnEleN[i][j];
			OnEleN[i][j] = NULL;
		}
		delete[] OnEleR[i];
		delete[] OnEleN[i];
		delete[] OnEleJ[i];
		OnEleR[i] = NULL;
		OnEleN[i] = NULL;
		OnEleJ[i] = NULL;
	}
	delete[] OnEleR;
	delete[] OnEleN;
	delete[] OnEleJ;
	OnEleR = NULL;
	OnEleN = NULL;
	OnEleJ = NULL;

	return 1;
}

int DSquareElement::ConvertUnknownCoef()
{
	int i, j, k, PointID;
	for (i = 0; i < 9; ++i)
	{
		for (int PointID = 0; PointID < 1; ++PointID)
			m_AssistCoef[PointID][i] = 0.0;
		//m_AssistCoef[1][i] = 0.0;
		//m_AssistCoef[2][i] = 0.0;
		//m_AssistCoef[3][i] = 0.0;
		//m_AssistCoef[4][i] = 0.0;
		//m_AssistCoef[5][i] = 0.0;
		//m_AssistCoef[6][i] = 0.0;
		//m_AssistCoef[7][i] = 0.0;
	}
	switch (BCID)
	{
		case 123:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
						for (k = 0; k < 3; ++k)
							m_AssistCoef[PointID][3 * j + i] -= m_AssistUij[PointID][3 * k + i] * m_transmat[m_nodeID[PointID]][3 * j + k];
			break;
		}
		case 126:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
					{
						m_AssistCoef[PointID][i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
						m_AssistCoef[PointID][3 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
						m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
					}
			break;
		}
		case 153:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
					{
						m_AssistCoef[PointID][i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
						m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
						m_AssistCoef[PointID][6 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
					}
			break;
		}
		case 156:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
					{
						m_AssistCoef[PointID][i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
						m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
						m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
					}
			break;
		}
		case 423:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
					{
						m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
						m_AssistCoef[PointID][3 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
						m_AssistCoef[PointID][6 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
					}
			break;
		}
		case 426:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
					{
						m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
						m_AssistCoef[PointID][3 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
						m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
					}
			break;
		}
		case 453:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
					{
						m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
						m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
						m_AssistCoef[PointID][6 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
					}
			break;
		}
		case 456:
		{
			for (PointID = 0; PointID < 1; ++PointID)
				for (i = 0; i < 3; ++i)
					for (j = 0; j < 3; ++j)
					{
						m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
						m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
						m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
					}
			break;
		}
	}
	return 1;
}

int DSquareElement::ConvertUnknownCoef(long T_ID)
{
	int i, j, k, PointID;
	for (i = 0; i < 9; ++i)
	{
		for (int PointID = 0; PointID < 1; ++PointID)
			AssistCoef[T_ID][PointID][i] = 0.0;
		//m_AssistCoef[1][i] = 0.0;
		//m_AssistCoef[2][i] = 0.0;
		//m_AssistCoef[3][i] = 0.0;
		//m_AssistCoef[4][i] = 0.0;
		//m_AssistCoef[5][i] = 0.0;
		//m_AssistCoef[6][i] = 0.0;
		//m_AssistCoef[7][i] = 0.0;
	}
	switch (BCID)
	{
	case 123:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
					for (k = 0; k < 3; ++k)
						AssistCoef[T_ID][PointID][3 * j + i] -= AssistUij[T_ID][PointID][3 * k + i] * m_transmat[m_nodeID[PointID]][3 * j + k];
		break;
	}
	case 126:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 153:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 156:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 423:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 426:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 453:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] -= AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 456:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] += AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	}
	return 1;
}

int DSquareElement::ConvertUnknownCoef(int ID)
{
	int i, j, k;
	for (i = 0; i < 9; ++i)
		m_AssistCoef[ID][i] = 0.0;

	switch (BCID)
	{
	case 123:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
				for (k = 0; k < 3; ++k)
					m_AssistCoef[ID][3 * j + i] -= m_AssistUij[ID][3 * k + i] * m_transmat[m_nodeID[ID]][3 * j + k];
		break;
	}
	case 126:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 153:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 156:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 423:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 426:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 453:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] -= m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 456:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] += m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	}
	return 1;
}

int DSquareElement::ConvertUnknownCoefforInf()
{
	int i, j, k, PointID;

	switch (BCID)
	{
	case 123:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
					for (k = 0; k < 3; ++k)
						m_AssistCoef[PointID][3 * j + i] -= m_AssistUij[PointID][3 * k + i] * m_transmat[m_nodeID[PointID]][3 * j + k];
		break;
	}
	case 126:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 153:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 156:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 423:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 426:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 453:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 456:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	}
	return 1;
}

int DSquareElement::ConvertKnownCoef()
{
	int i, j, k, PointID;
	for (i = 0; i < 9; ++i)
	{
		for (int PointID = 0; PointID < 1; ++PointID)
			m_AssistCoef[PointID][i] = 0.0;
		//m_AssistCoef[1][i] = 0.0;
		//m_AssistCoef[2][i] = 0.0;
		//m_AssistCoef[3][i] = 0.0;
		//m_AssistCoef[4][i] = 0.0;
		//m_AssistCoef[5][i] = 0.0;
		//m_AssistCoef[6][i] = 0.0;
		//m_AssistCoef[7][i] = 0.0;
	}
	switch (BCID)
	{
	case 456:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
					for (k = 0; k < 3; ++k)
						m_AssistCoef[PointID][3 * j + i] += m_AssistUij[PointID][3 * k + i] * m_transmat[m_nodeID[PointID]][3 * j + k];
		break;
	}
	case 453:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 426:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 423:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 156:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 153:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 126:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 123:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	}
	return 1;
}

int DSquareElement::ConvertKnownCoef(int ID)
{
	int i, j, k;
	for (i = 0; i < 9; ++i)
		m_AssistCoef[ID][i] = 0.0;

	switch (BCID)
	{
	case 456:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
				for (k = 0; k < 3; ++k)
					m_AssistCoef[ID][3 * j + i] += m_AssistUij[ID][3 * k + i] * m_transmat[m_nodeID[ID]][3 * j + k];
		break;
	}
	case 453:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 426:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 423:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 156:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 153:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 126:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] += m_AssistUij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	case 123:
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				m_AssistCoef[ID][i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][j];
				m_AssistCoef[ID][3 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][3 + j];
				m_AssistCoef[ID][6 + i] -= m_AssistTij[ID][3 * j + i] * m_transmat[m_nodeID[ID]][6 + j];
			}
		break;
	}
	}
	return 1;
}

int DSquareElement::ConvertKnownCoef(long T_ID)
{
	int i, j, k, PointID;
	for (i = 0; i < 9; ++i)
	{
		for (int PointID = 0; PointID < 1; ++PointID)
			AssistCoef[T_ID][PointID][i] = 0.0;
		//m_AssistCoef[1][i] = 0.0;
		//m_AssistCoef[2][i] = 0.0;
		//m_AssistCoef[3][i] = 0.0;
		//m_AssistCoef[4][i] = 0.0;
		//m_AssistCoef[5][i] = 0.0;
		//m_AssistCoef[6][i] = 0.0;
		//m_AssistCoef[7][i] = 0.0;
	}
	switch (BCID)
	{
	case 456:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
					for (k = 0; k < 3; ++k)
						AssistCoef[T_ID][PointID][3 * j + i] += AssistUij[T_ID][PointID][3 * k + i] * m_transmat[m_nodeID[PointID]][3 * j + k];
		break;
	}
	case 453:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 426:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 423:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 156:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 153:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 126:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] += AssistUij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 123:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					AssistCoef[T_ID][PointID][i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					AssistCoef[T_ID][PointID][3 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					AssistCoef[T_ID][PointID][6 + i] -= AssistTij[T_ID][PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	}
	return 1;
}

int DSquareElement::ConvertKnownCoefforInf()
{
	int i, j, k, PointID;

	switch (BCID)
	{
	case 456:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
					for (k = 0; k < 3; ++k)
						m_AssistCoef[PointID][3 * j + i] += m_AssistUij[PointID][3 * k + i] * m_transmat[m_nodeID[PointID]][3 * j + k];
		break;
	}
	case 453:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 426:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 423:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 156:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 153:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 126:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] += m_AssistUij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	case 123:
	{
		for (PointID = 0; PointID < 1; ++PointID)
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					m_AssistCoef[PointID][i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][j];
					m_AssistCoef[PointID][3 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][3 + j];
					m_AssistCoef[PointID][6 + i] -= m_AssistTij[PointID][3 * j + i] * m_transmat[m_nodeID[PointID]][6 + j];
				}
		break;
	}
	}
	return 1;
}

//静力学使用
int DSquareElement::RegularIntUnknownCoefforStatic(Point& source)
{
	
	if (BCID == 123)//Uij
	{
		IntUij(source);
	}
	else if (BCID == 456)//Tij
	{
		IntTij(source);
	}
	else//both Uij and Tij
	{
		IntUij(source);
		IntTij(source);
	}

	ConvertUnknownCoef();

	return 1;
}
//无限单元未知系数矩阵
int DSquareElement::IntUnknownCoefforInf(DSquareElement* m_elelist, long NodeID, long FieldEleID, long ** InfElePID, long** ElePID,int m_InfEleFlag)
{
	int i, j;

	if (DSquareElement::Flagproblem == 0)
	{
		if (DSquareElement::FlagDynaMethod == 1)
		{
			m_elelist[FieldEleID].IntUijforInfEle();
			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 1);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 2);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 3);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 4);

			}

		}
		else {
			//稳定性处理

		// 常速度插值方法
			double  TT[8][9];

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] = 0.0;
				}
			}
			// G[0], T[1,0] and T[2,1]
			m_elelist[FieldEleID].IntUijforInfEle();//远处条件，赋值为0
		

			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 1);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 2);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 3);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 4);

			}

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += 4.0 * DSquareElement::m_AssistTij[i][j];
				}
			}

			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 1);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 2);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 3);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 4);

			}

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += DSquareElement::m_AssistTij[i][j];
				}
			}


			// T[1,1]
			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 1);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 2);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 3);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 4);

			}

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += DSquareElement::m_AssistTij[i][j];
				}
			}

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{

					DSquareElement::m_AssistTij[i][j] = TT[i][j];
				}
			}
		}
	}
	else if (DSquareElement::Flagproblem == 1)
	{
		//判断无限扩展方向
		m_elelist[FieldEleID].IntUijforInfEle();
		if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 1);

		}
		else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 2);

		}
		else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 3);

		}
		else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 4);

		}

	}

	m_elelist[FieldEleID].ConvertUnknownCoefforInf();//集成

	return 0;
}
//无限单元已知系数矩阵
int DSquareElement::IntKnownCoefforInf(DSquareElement* m_elelist, long NodeID, long FieldEleID, long ** InfElePID, long** ElePID, int m_InfEleFlag)
{

	int i, j;

	if (DSquareElement::Flagproblem == 0)
	{
		if (DSquareElement::FlagDynaMethod == 1)
		{
			m_elelist[FieldEleID].IntUijforInfEle();
			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {
				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 1);
			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {
				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 2);
			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {
				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 3);
			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {
				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 4);
			}
		}
		else
		{
			//稳定性处理
			// 常速度插值方法
			double TT[8][9];
			
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] = 0.0;
				}
			}

			// G[0], T[1,0] and T[2,1]
			m_elelist[FieldEleID].IntUijforInfEle();//远处条件，赋值为0

			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 1);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 2);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 3);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 0.0, DSquareElement::m_DMat.Dt, 4);

			}

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += 4.0 * m_AssistTij[i][j];
				}
			}

			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 1);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 2);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 3);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(2, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 4);

			}

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += m_AssistTij[i][j];
				}
			}

			// G[1], T[1,1]
			if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 1);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 2);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 3);

			}
			else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {

				m_elelist[FieldEleID].IntDynaTijPWforInfEle(1, m_elelist[FieldEleID].m_nodelist[NodeID], 1.0, DSquareElement::m_DMat.Dt, 4);
			}

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += m_AssistTij[i][j];
				}
			}
			
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					
					m_AssistTij[i][j] = TT[i][j];
				}
			}
		}
	}
	else if(DSquareElement::Flagproblem == 1)
	{
		if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][4]) {
			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 1);
		}
		else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][5]) {
			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 2);
		}
		else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][6]) {
			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 3);
		}
		else if (InfElePID[m_InfEleFlag][2] == ElePID[FieldEleID][7]) {
			m_elelist[FieldEleID].IntTijforInfEle(m_elelist[FieldEleID].m_nodelist[NodeID], 4);
		}
	}
	m_elelist[FieldEleID].ConvertKnownCoefforInf();//集成

	return 1;
}

int DSquareElement::RegularIntKnownCoefforStatic(Point& source)
{
	if (BCID == 456)//Uij
	{
		IntUij(source);
	}
	else if (BCID == 123)//Tij
	{
		IntTij(source);
	}
	else//both Uij and Tij
	{
		IntUij(source);
		IntTij( source);
	}
	ConvertKnownCoef();


	return 2;
}

int DSquareElement::RegularIntUnknownCoef(Point& source)
{
	int Flag, FlagSum;
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)
	{
		// 常规方法

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt);
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt);
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt);
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);
				IntDynaTij(1, source, 0.0, m_DMat.Dt);
			}
		}

		if (FlagSum)
		{
			ConvertUnknownCoef();
		}
	}
	else
	{
		// 常速度插值方法

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}

				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}

				IntDynaTij(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}

		// G[1], T[1,1]

		Flag = NeedCal(1, 1.0, source);

		FlagSum += Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[1]
				IntDynaUijPW(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}

				IntDynaTijPW(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}

				IntDynaTij(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}

		if (FlagSum)
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistUij[i][j] = TU[i][j];
					m_AssistTij[i][j] = TT[i][j];
				}
			}

			ConvertUnknownCoef();
		}
	}

	return FlagSum;
}

int DSquareElement::RegularIntUnknownCoef(Point& source, long T_ID)
{
	int Flag, FlagSum;
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)
	{
		// 常规方法   不做改动

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt, T_ID);
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, T_ID);
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, T_ID);
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, T_ID);
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);
				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);
			}
		}

		if (FlagSum)
		{
			ConvertUnknownCoef(T_ID);
		}
	}
	else
	{
		// 常速度插值方法 并行   做相关调整     高斯点、雅可比、法向量  —— 使用时计算

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt,T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt,T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt,T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}

				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}

				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}

		// G[1], T[1,1]

		Flag = NeedCal(1, 1.0, source);

		FlagSum += Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[1]
				IntDynaUijPW(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}

				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}

				IntDynaTij(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}

		if (FlagSum)
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistUij[T_ID][i][j] = TU[i][j];
					AssistTij[T_ID][i][j] = TT[i][j];
				}
			}

			ConvertUnknownCoef(T_ID);
		}
	}

	return FlagSum;
}

int DSquareElement::RegularIntUnknownCoef(Point& source, int ID)
{
	int Flag, FlagSum;
	int i;
	double TU[9], TT[9];

	if (FlagDynaMethod == 1)
	{
		// 标准方法

		Flag = NeedCal(1, 0.0, source);
		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);
				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);
			}
		}

		if (FlagSum)
		{
			ConvertUnknownCoef(ID);
		}
	}
	else
	{
		// 常速度插值方法

		for (i = 0; i < 9; ++i)
		{
			TU[i] = 0.0;
			TT[i] = 0.0;
		}


		// G[0], T[1,0] and T[2,1]

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}

				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}

				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}

		// G[1], T[1,1]

		Flag = NeedCal(1, 1.0, source);

		FlagSum += Flag;

		if (Flag == 1)
		{
			if (BCID == 123)//Uij
			{
				// G[1]
				IntDynaUijPW(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}
			}
			else if (BCID == 456) // Tij
			{
				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}

				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 123)//Uij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}
			}
			else if (BCID == 456)//Tij
			{
				IntDynaTij(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}

				IntDynaTij(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}

		if (FlagSum)
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistUij[ID][i] = TU[i];
				m_AssistTij[ID][i] = TT[i];
			}

			ConvertUnknownCoef(ID);
		}
	}

	return FlagSum;
}

int DSquareElement::RegularIntKnownCoef(Point& source)
{
	int Flag, FlagSum;
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)
	{
		// 常规方法

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt);
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt);
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt);
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);
				IntDynaTij(1, source, 0.0, m_DMat.Dt);
			}
		}

		if (FlagSum)
		{
			ConvertKnownCoef();
		}
	}
	else
	{
		// 常速度插值方法

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}

				IntDynaTijPW(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * m_AssistUij[i][j];
					}
				}

				IntDynaTij(1, source, 0.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * m_AssistTij[i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}

		// G[1], T[1,1]

		Flag = NeedCal(1, 1.0, source);

		FlagSum += Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				// G[1]
				IntDynaUijPW(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}

				IntDynaTijPW(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += m_AssistUij[i][j];
					}
				}

				IntDynaTij(1, source, 1.0, m_DMat.Dt);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += m_AssistTij[i][j];
					}
				}
			}
		}

		if (FlagSum)
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistUij[i][j] = TU[i][j];
					m_AssistTij[i][j] = TT[i][j];
				}
			}

			ConvertKnownCoef();
		}
	}

	return FlagSum;
}

int DSquareElement::RegularIntKnownCoef(Point& source, long T_ID)
{
	int Flag, FlagSum;
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)
	{
		// 常规方法

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, T_ID);
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, T_ID);
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, T_ID);
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, T_ID);
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);
				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);
			}
		}

		if (FlagSum)
		{
			ConvertKnownCoef(T_ID);
		}
	}
	else
	{
		// 常速度插值方法

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}

				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += 4.0 * AssistUij[T_ID][i][j];
					}
				}

				IntDynaTij(1, source, 0.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += 4.0 * AssistTij[T_ID][i][j];
					}
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}

		// G[1], T[1,1]

		Flag = NeedCal(1, 1.0, source);

		FlagSum += Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				// G[1]
				IntDynaUijPW(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}

				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TU[i][j] += AssistUij[T_ID][i][j];
					}
				}

				IntDynaTij(1, source, 1.0, m_DMat.Dt, T_ID);

				for (i = 0; i < 8; ++i)
				{
					for (j = 0; j < 9; ++j)
					{
						TT[i][j] += AssistTij[T_ID][i][j];
					}
				}
			}
		}

		if (FlagSum)
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistUij[T_ID][i][j] = TU[i][j];
					AssistTij[T_ID][i][j] = TT[i][j];
				}
			}

			ConvertKnownCoef(T_ID);
		}
	}

	return FlagSum;
}

int DSquareElement::RegularIntKnownCoef(Point& source, int ID)
{
	int Flag, FlagSum;
	int i;
	double TU[9], TT[9];

	if (FlagDynaMethod == 1)
	{
		// 常规方法

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);
				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);
			}
		}

		if (FlagSum)
		{
			ConvertKnownCoef(ID);
		}
	}
	else
	{
		// 常速度插值方法

		for (i = 0; i < 9; ++i)
		{
			TU[i] = 0.0;
			TT[i] = 0.0;
		}


		// G[0], T[1,0] and T[2,1]

		Flag = NeedCal(1, 0.0, source);

		FlagSum = Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				// G[0]
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}

				IntDynaTijPW(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTijPW(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += 4.0 * m_AssistUij[ID][i];
				}

				IntDynaTij(1, source, 0.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += 4.0 * m_AssistTij[ID][i];
				}

				IntDynaTij(2, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}

		// G[1], T[1,1]

		Flag = NeedCal(1, 1.0, source);

		FlagSum += Flag;

		if (Flag == 1)
		{
			if (BCID == 456)//Uij
			{
				// G[1]
				IntDynaUijPW(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}
			}
			else if (BCID == 123) // Tij
			{
				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else // both Uij and Tij
			{
				IntDynaUijPW(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}

				IntDynaTijPW(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}
		else if (Flag == 2)
		{
			if (BCID == 456)//Uij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}
			}
			else if (BCID == 123)//Tij
			{
				IntDynaTij(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
			else//both Uij and Tij
			{
				IntDynaUij(source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TU[i] += m_AssistUij[ID][i];
				}

				IntDynaTij(1, source, 1.0, m_DMat.Dt, ID);

				for (i = 0; i < 9; ++i)
				{
					TT[i] += m_AssistTij[ID][i];
				}
			}
		}

		if (FlagSum)
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistUij[ID][i] = TU[i];
				m_AssistTij[ID][i] = TT[i];
			}

			ConvertKnownCoef(ID);
		}
	}

	return FlagSum;
}

int DSquareElement::InElementIntUnknownCoef(int sourceID)
{
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)//常单元未修改
	{
		// 标准方法

		if (BCID == 123)//Uij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistUij[i][j] = m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 456) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistTij[i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistUij[i][j] = m_OnElementU1ij[sourceID][i][j];
					m_AssistTij[i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
	}
	else// 常单元已修改
	{
		// 常速度插值方法

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		if (BCID == 123)//Uij
		{
			// G[0]

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 456) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}

		// G[1], T[1,1]

		if (BCID == 123)//Uij
		{
			// G[1]
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += m_AssistUij[i][j];
				}
			}
		}
		else if (BCID == 456) // Tij
		{
			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += m_AssistTij[i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += m_AssistUij[i][j];
				}
			}

			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += m_AssistTij[i][j];
				}
			}
		}

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				m_AssistUij[i][j] = TU[i][j];
				m_AssistTij[i][j] = TT[i][j];
			}
		}
	}

	ConvertUnknownCoef();

	return 1;
}

int DSquareElement::InElementIntUnknownCoef(int sourceID, long T_ID)
{
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)//常单元未修改
	{
		// 标准方法

		if (BCID == 123)//Uij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistUij[T_ID][i][j] = m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 456) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistTij[T_ID][i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistUij[T_ID][i][j] = m_OnElementU1ij[sourceID][i][j];
					AssistTij[T_ID][i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
	}
	else// 常单元已修改
	{
		// 常速度插值方法

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		if (BCID == 123)//Uij
		{
			// G[0]

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 456) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}

		// G[1], T[1,1]

		if (BCID == 123)//Uij
		{
			// G[1]
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt,T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += AssistUij[T_ID][i][j];
				}
			}
		}
		else if (BCID == 456) // Tij
		{
			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += AssistTij[T_ID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += AssistUij[T_ID][i][j];
				}
			}

			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += AssistTij[T_ID][i][j];
				}
			}
		}

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				AssistUij[T_ID][i][j] = TU[i][j];
				AssistTij[T_ID][i][j] = TT[i][j];
			}
		}
	}

	ConvertUnknownCoef(T_ID);

	return 1;
}

int DSquareElement::InElementIntUnknownCoef(int sourceID, int ID)
{
	int i;
	double TU[9], TT[9];

	if (FlagDynaMethod == 1)
	{
		// 常规方法

		if (BCID == 123)//Uij
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistUij[ID][i] = m_OnElementU1ij[sourceID][ID][i];
			}
		}
		else if (BCID == 456) // Tij
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistTij[ID][i] = m_OnElementT1ij[sourceID][ID][i];
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistUij[ID][i] = m_OnElementU1ij[sourceID][ID][i];
				m_AssistTij[ID][i] = m_OnElementT1ij[sourceID][ID][i];
			}
		}
	}
	else
	{
		// 常速度插值法

		for (i = 0; i < 9; ++i)
		{
			TU[i] = 0.0;
			TT[i] = 0.0;
		}


		// G[0], T[1,0] and T[2,1]

		if (BCID == 123)//Uij
		{
			// G[0]

			for (i = 0; i < 9; ++i)
			{
				TU[i] += 4.0 * m_OnElementU1ij[sourceID][ID][i];
			}
		}
		else if (BCID == 456) // Tij
		{
			for (i = 0; i < 9; ++i)
			{
				TT[i] += 4.0 * m_OnElementT1ij[sourceID][ID][i];
				TT[i] += m_OnElementT2ij[sourceID][ID][i];
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 9; ++i)
			{
				TU[i] += 4.0 * m_OnElementU1ij[sourceID][ID][i];
				TT[i] += 4.0 * m_OnElementT1ij[sourceID][ID][i];
				TT[i] += m_OnElementT2ij[sourceID][ID][i];
			}
		}

		// G[1], T[1,1]

		if (BCID == 123)//Uij
		{
			// G[1]
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TU[i] += m_AssistUij[ID][i];
			}
		}
		else if (BCID == 456) // Tij
		{
			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TT[i] += m_AssistTij[ID][i];
			}
		}
		else // both Uij and Tij
		{
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TU[i] += m_AssistUij[ID][i];
			}

			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TT[i] += m_AssistTij[ID][i];
			}
		}

		for (i = 0; i < 9; ++i)
		{
			m_AssistUij[ID][i] = TU[i];
			m_AssistTij[ID][i] = TT[i];
		}
	}

	ConvertUnknownCoef(ID);

	return 1;
}

int DSquareElement::InElementIntKnownCoef(int sourceID)
{
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)
	{

		if (BCID == 456)//Uij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistUij[i][j] = m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 123) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistTij[i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					m_AssistUij[i][j] = m_OnElementU1ij[sourceID][i][j];
					m_AssistTij[i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
	}
	else
	{

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		if (BCID == 456)//Uij
		{
			// G[0]

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 123) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}

		// G[1], T[1,1]

		if (BCID == 456)//Uij
		{
			// G[1]
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += m_AssistUij[i][j];
				}
			}
		}
		else if (BCID == 123) // Tij
		{
			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += m_AssistTij[i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += m_AssistUij[i][j];
				}
			}

			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += m_AssistTij[i][j];
				}
			}
		}

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				m_AssistUij[i][j] = TU[i][j];
				m_AssistTij[i][j] = TT[i][j];
			}
		}
	}

	ConvertKnownCoef();

	return 1;
}

int DSquareElement::InElementIntKnownCoef(int sourceID, int ID)
{
	int i;
	double TU[9], TT[9];

	if (FlagDynaMethod == 1)
	{
		if (BCID == 456)//Uij
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistUij[ID][i] = m_OnElementU1ij[sourceID][ID][i];
			}
		}
		else if (BCID == 123) // Tij
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistTij[ID][i] = m_OnElementT1ij[sourceID][ID][i];
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 9; ++i)
			{
				m_AssistUij[ID][i] = m_OnElementU1ij[sourceID][ID][i];
				m_AssistTij[ID][i] = m_OnElementT1ij[sourceID][ID][i];
			}
		}
	}
	else
	{
		for (i = 0; i < 9; ++i)
		{
			TU[i] = 0.0;
			TT[i] = 0.0;
		}


		// G[0], T[1,0] and T[2,1]

		if (BCID == 456)//Uij
		{
			// G[0]

			for (i = 0; i < 9; ++i)
			{
				TU[i] += 4.0 * m_OnElementU1ij[sourceID][ID][i];
			}
		}
		else if (BCID == 123) // Tij
		{
			for (i = 0; i < 9; ++i)
			{
				TT[i] += 4.0 * m_OnElementT1ij[sourceID][ID][i];
				TT[i] += m_OnElementT2ij[sourceID][ID][i];
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 9; ++i)
			{
				TU[i] += 4.0 * m_OnElementU1ij[sourceID][ID][i];
				TT[i] += 4.0 * m_OnElementT1ij[sourceID][ID][i];
				TT[i] += m_OnElementT2ij[sourceID][ID][i];
			}
		}

		// G[1], T[1,1]

		if (BCID == 456)//Uij
		{
			// G[1]
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TU[i] += m_AssistUij[ID][i];
			}
		}
		else if (BCID == 123) // Tij
		{
			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TT[i] += m_AssistTij[ID][i];
			}
		}
		else // both Uij and Tij
		{
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TU[i] += m_AssistUij[ID][i];
			}

			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, ID);

			for (i = 0; i < 9; ++i)
			{
				TT[i] += m_AssistTij[ID][i];
			}
		}

		for (i = 0; i < 9; ++i)
		{
			m_AssistUij[ID][i] = TU[i];
			m_AssistTij[ID][i] = TT[i];
		}
	}

	ConvertKnownCoef(ID);

	return 1;
}

int DSquareElement::InElementIntKnownCoef(int sourceID, long T_ID)
{
	int i, j;
	double TU[8][9], TT[8][9];

	if (FlagDynaMethod == 1)
	{

		if (BCID == 456)//Uij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistUij[T_ID][i][j] = m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 123) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistTij[T_ID][i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					AssistUij[T_ID][i][j] = m_OnElementU1ij[sourceID][i][j];
					AssistTij[T_ID][i][j] = m_OnElementT1ij[sourceID][i][j];
				}
			}
		}
	}
	else
	{

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				TU[i][j] = 0.0;
				TT[i][j] = 0.0;
			}
		}


		// G[0], T[1,0] and T[2,1]

		if (BCID == 456)//Uij
		{
			// G[0]

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
				}
			}
		}
		else if (BCID == 123) // Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += 4.0 * m_OnElementU1ij[sourceID][i][j];
					TT[i][j] += 4.0 * m_OnElementT1ij[sourceID][i][j];
					TT[i][j] += m_OnElementT2ij[sourceID][i][j];
				}
			}
		}

		// G[1], T[1,1]

		if (BCID == 456)//Uij
		{
			// G[1]
			//IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt);

			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += AssistUij[T_ID][i][j];
				}
			}
		}
		else if (BCID == 123) // Tij
		{
			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += AssistTij[T_ID][i][j];
				}
			}
		}
		else // both Uij and Tij
		{
			IntDynaUijPW(m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TU[i][j] += AssistUij[T_ID][i][j];
				}
			}

			IntDynaTijPW(1, m_nodelist[m_nodeID[sourceID]], 1.0, m_DMat.Dt, T_ID);

			for (i = 0; i < 1; ++i)
			{
				for (j = 0; j < 9; ++j)
				{
					TT[i][j] += AssistTij[T_ID][i][j];
				}
			}
		}

		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < 9; ++j)
			{
				AssistUij[T_ID][i][j] = TU[i][j];
				AssistTij[T_ID][i][j] = TT[i][j];
			}
		}
	}

	ConvertKnownCoef(T_ID);

	return 1;
}



//-----------------END-------------------------------------------------------

//Functions of operation for DSquareElement

int UnknownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID)
{
	int pos = -1;

	if (m_elelist[FieldEleID].IsIn(NodeID, pos))
	{
		m_elelist[FieldEleID].InElementIntUnknownCoef(pos);
		return 1;
	}
	else
	{
		if (DSquareElement::Flagproblem == 0) {
			return m_elelist[FieldEleID].RegularIntUnknownCoef(m_elelist[FieldEleID].m_nodelist[NodeID]);
		}
		else if (DSquareElement::Flagproblem == 1) {
			return m_elelist[FieldEleID].RegularIntUnknownCoefforStatic(m_elelist[FieldEleID].m_nodelist[NodeID]);
		}
		
	}

	return 1;
}
int UnknownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID, long** InfElePID, long** ElePID)
{
	int pos = -1, needcalflag = 0;
	int i;
	long m_InfEleFlag;

	//8个点3*3矩阵
	for (i = 0; i < 9; ++i)
	{

		for (int PointID = 0; PointID < 8; ++PointID)
			DSquareElement::m_AssistCoef[PointID][i] = 0.0;
		//DSquareElement::m_AssistCoef[1][i] = 0.0;
		//DSquareElement::m_AssistCoef[2][i] = 0.0;
		//DSquareElement::m_AssistCoef[3][i] = 0.0;
		//DSquareElement::m_AssistCoef[4][i] = 0.0;
		//DSquareElement::m_AssistCoef[5][i] = 0.0;
		//DSquareElement::m_AssistCoef[6][i] = 0.0;
		//DSquareElement::m_AssistCoef[7][i] = 0.0;
	}


	if (m_elelist[FieldEleID].IsIn(NodeID, pos))
	{
		m_elelist[FieldEleID].InElementIntUnknownCoef(pos);
		needcalflag++;

	}
	else
	{

		if (DSquareElement::Flagproblem == 0) {
			needcalflag += m_elelist[FieldEleID].RegularIntUnknownCoef(m_elelist[FieldEleID].m_nodelist[NodeID]);
		}
		else if (DSquareElement::Flagproblem == 1)
		{
			m_elelist[FieldEleID].RegularIntUnknownCoefforStatic(m_elelist[FieldEleID].m_nodelist[NodeID]);
			needcalflag++;
		}
		//无限单元
		//m_InfEleFlag = m_elelist[FieldEleID].InfEleFlag;
		//if (m_InfEleFlag != (-1))
		//{
		//	m_elelist[FieldEleID].IntUnknownCoefforInf(m_elelist, NodeID, FieldEleID, InfElePID, ElePID, m_InfEleFlag);
		//	needcalflag++;
		//}
	}
	//无限单元
	m_InfEleFlag = m_elelist[FieldEleID].InfEleFlag;
	if (m_InfEleFlag != (-1))
	{
		m_elelist[FieldEleID].IntUnknownCoefforInf(m_elelist, NodeID, FieldEleID, InfElePID, ElePID, m_InfEleFlag);
		needcalflag++;
	}
	return  needcalflag;
}
//增加无限单元适用
int UnknownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID,long ** InfElePID,long** ElePID,long T_ID)
{
	int pos = -1,needcalflag=0;
	int i;
	long m_InfEleFlag;

	//8个点3*3矩阵
	for (i = 0; i < 9; ++i)
	{

		for (int PointID = 0; PointID < 8; ++PointID)
			AssistCoef[T_ID][PointID][i] = 0.0;
		//DSquareElement::m_AssistCoef[1][i] = 0.0;
		//DSquareElement::m_AssistCoef[2][i] = 0.0;
		//DSquareElement::m_AssistCoef[3][i] = 0.0;
		//DSquareElement::m_AssistCoef[4][i] = 0.0;
		//DSquareElement::m_AssistCoef[5][i] = 0.0;
		//DSquareElement::m_AssistCoef[6][i] = 0.0;
		//DSquareElement::m_AssistCoef[7][i] = 0.0;
	}


	if (m_elelist[FieldEleID].IsIn(NodeID, pos))
	{
		m_elelist[FieldEleID].InElementIntUnknownCoef(pos, T_ID);
		needcalflag++;
		
	}
	else
	{
		
		if (DSquareElement::Flagproblem == 0) {//仅修改动力学
			needcalflag+=m_elelist[FieldEleID].RegularIntUnknownCoef(m_elelist[FieldEleID].m_nodelist[NodeID],T_ID);		
		}
		else if (DSquareElement::Flagproblem == 1) //
		{
			m_elelist[FieldEleID].RegularIntUnknownCoefforStatic(m_elelist[FieldEleID].m_nodelist[NodeID]);
			needcalflag++;
		}
		//无限单元
		//m_InfEleFlag = m_elelist[FieldEleID].InfEleFlag;
		//if (m_InfEleFlag != (-1))
		//{
		//	m_elelist[FieldEleID].IntUnknownCoefforInf(m_elelist, NodeID, FieldEleID, InfElePID, ElePID, m_InfEleFlag);
		//	needcalflag++;
		//}
	}
	//无限单元
	m_InfEleFlag = m_elelist[FieldEleID].InfEleFlag;
	if (m_InfEleFlag != (-1))
	{
		m_elelist[FieldEleID].IntUnknownCoefforInf(m_elelist, NodeID, FieldEleID, InfElePID, ElePID, m_InfEleFlag);//无线单元暂未并行
		needcalflag++;
	}
	return  needcalflag;
}

int NodeUnknownSubMatrix(DSquareElement* m_elelist, long SourceNodeID, long FieldNodeID)
{
	// 使用DSquareElement::m_AssistCoef[0]存储系数子矩阵
	int pos, EleID, posID, Flag;
	double var[9];

	{
		EleID = DSquareElement::BeEleID[FieldNodeID];
		posID = DSquareElement::BeElePos[FieldNodeID];

		pos = -1;

		if (m_elelist[EleID].IsIn(SourceNodeID, pos))
		{
			m_elelist[EleID].InElementIntUnknownCoef(pos, posID);
		}
		else
		{
			Flag = m_elelist[EleID].RegularIntUnknownCoef(m_elelist[EleID].m_nodelist[SourceNodeID], posID);

			if (Flag == 0)
				return 0;
		}

		var[0] = DSquareElement::m_AssistCoef[posID][0];
		var[1] = DSquareElement::m_AssistCoef[posID][1];
		var[2] = DSquareElement::m_AssistCoef[posID][2];
		var[3] = DSquareElement::m_AssistCoef[posID][3];
		var[4] = DSquareElement::m_AssistCoef[posID][4];
		var[5] = DSquareElement::m_AssistCoef[posID][5];
		var[6] = DSquareElement::m_AssistCoef[posID][6];
		var[7] = DSquareElement::m_AssistCoef[posID][7];
		var[8] = DSquareElement::m_AssistCoef[posID][8];
	}

	DSquareElement::m_AssistCoef[0][0] = var[0];
	DSquareElement::m_AssistCoef[0][1] = var[1];
	DSquareElement::m_AssistCoef[0][2] = var[2];
	DSquareElement::m_AssistCoef[0][3] = var[3];
	DSquareElement::m_AssistCoef[0][4] = var[4];
	DSquareElement::m_AssistCoef[0][5] = var[5];
	DSquareElement::m_AssistCoef[0][6] = var[6];
	DSquareElement::m_AssistCoef[0][7] = var[7];
	DSquareElement::m_AssistCoef[0][8] = var[8];

	return 1;
}

int knownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID)
{
	int pos = -1;

	if (m_elelist[FieldEleID].IsIn(NodeID, pos))
	{
		m_elelist[FieldEleID].InElementIntKnownCoef(pos);
		return 1;
	}
	else
	{
		if (DSquareElement::Flagproblem == 0) {
			return m_elelist[FieldEleID].RegularIntKnownCoef(m_elelist[FieldEleID].m_nodelist[NodeID]);
		}
		else if (DSquareElement::Flagproblem == 1) {
			return m_elelist[FieldEleID].RegularIntKnownCoefforStatic(m_elelist[FieldEleID].m_nodelist[NodeID]);
		}
	}

	return 1;
}
//增加无限单元适用
int knownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID, long ** InfElePID, long** ElePID)
{
	int pos = -1, needcalflag = 0;
	int i;
	long m_InfEleFlag;

	for (i = 0; i < 9; ++i)//8个点3*3矩阵
	{
		for (int PointID = 0; PointID < 8; ++PointID)
			DSquareElement::m_AssistCoef[PointID][i] = 0.0;
	/*	DSquareElement::m_AssistCoef[1][i] = 0.0;
		DSquareElement::m_AssistCoef[2][i] = 0.0;
		DSquareElement::m_AssistCoef[3][i] = 0.0;
		DSquareElement::m_AssistCoef[4][i] = 0.0;
		DSquareElement::m_AssistCoef[5][i] = 0.0;
		DSquareElement::m_AssistCoef[6][i] = 0.0;
		DSquareElement::m_AssistCoef[7][i] = 0.0;*/
	}

	if (m_elelist[FieldEleID].IsIn(NodeID, pos))
	{
		m_elelist[FieldEleID].InElementIntKnownCoef(pos);
		needcalflag++;
	}
	else
	{
		if (DSquareElement::Flagproblem == 0) {
			needcalflag += m_elelist[FieldEleID].RegularIntKnownCoef(m_elelist[FieldEleID].m_nodelist[NodeID]);
		}
		else if (DSquareElement::Flagproblem == 1) {
			m_elelist[FieldEleID].RegularIntKnownCoefforStatic(m_elelist[FieldEleID].m_nodelist[NodeID]);
			needcalflag++;
		}

	}
	//无限单元部分
	m_InfEleFlag = m_elelist[FieldEleID].InfEleFlag;
	if (m_InfEleFlag != (-1))
	{
		m_elelist[FieldEleID].IntKnownCoefforInf(m_elelist, NodeID, FieldEleID, InfElePID, ElePID, m_InfEleFlag);
		needcalflag++;
	}
	return needcalflag;
}

int knownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID, long** InfElePID, long** ElePID, long T_ID)
{
	int pos = -1;

	if (m_elelist[FieldEleID].IsIn(NodeID, pos))
	{
		m_elelist[FieldEleID].InElementIntKnownCoef(pos, T_ID);
		return 1;
	}
	else
	{
		if (DSquareElement::Flagproblem == 0) {
			return m_elelist[FieldEleID].RegularIntKnownCoef(m_elelist[FieldEleID].m_nodelist[NodeID], T_ID);
		}
		else if (DSquareElement::Flagproblem == 1) {
			return m_elelist[FieldEleID].RegularIntKnownCoefforStatic(m_elelist[FieldEleID].m_nodelist[NodeID]);
		}
	}

	return 1;
}

int NodeknownSubMatrix(DSquareElement* m_elelist, long SourceNodeID, long FieldNodeID)
{
	// 使用DSquareElement::m_AssistCoef[0]存储系数子矩阵
	int pos, EleID, posID, Flag;
	double var[9];

	{
		EleID = DSquareElement::BeEleID[FieldNodeID];
		posID = DSquareElement::BeElePos[FieldNodeID];

		pos = -1;

		if (m_elelist[EleID].IsIn(SourceNodeID, pos))
		{
			m_elelist[EleID].InElementIntKnownCoef(pos, posID);
		}
		else
		{
			Flag = m_elelist[EleID].RegularIntKnownCoef(m_elelist[EleID].m_nodelist[SourceNodeID], posID);

			if (Flag == 0)
				return 0;
		}

		var[0] = DSquareElement::m_AssistCoef[posID][0];
		var[1] = DSquareElement::m_AssistCoef[posID][1];
		var[2] = DSquareElement::m_AssistCoef[posID][2];
		var[3] = DSquareElement::m_AssistCoef[posID][3];
		var[4] = DSquareElement::m_AssistCoef[posID][4];
		var[5] = DSquareElement::m_AssistCoef[posID][5];
		var[6] = DSquareElement::m_AssistCoef[posID][6];
		var[7] = DSquareElement::m_AssistCoef[posID][7];
		var[8] = DSquareElement::m_AssistCoef[posID][8];
	}

	DSquareElement::m_AssistCoef[0][0] = var[0];
	DSquareElement::m_AssistCoef[0][1] = var[1];
	DSquareElement::m_AssistCoef[0][2] = var[2];
	DSquareElement::m_AssistCoef[0][3] = var[3];
	DSquareElement::m_AssistCoef[0][4] = var[4];
	DSquareElement::m_AssistCoef[0][5] = var[5];
	DSquareElement::m_AssistCoef[0][6] = var[6];
	DSquareElement::m_AssistCoef[0][7] = var[7];
	DSquareElement::m_AssistCoef[0][8] = var[8];

	return 1;
}
//////
