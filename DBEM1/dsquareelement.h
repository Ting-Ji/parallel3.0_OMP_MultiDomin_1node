#ifndef _DSQUAREELEMENT_
#define _DSQUAREELEMENT_

#include "Geo.h"
#include "GaussNumeMatrix.h"
#include "mat_vec.h"
#include "quadinfo.h"
//time count
#include "Counter.h"
//#include "DBEM.h"
//end time count


//--------------------------------------------------------
//  Definition of discontinuous square element for 3D
//  elastostatics BEM
//
//  Author: Haitao Wang
//--------------------------------------------------------

// 动态材料属性类
class DynaMat
{
public:
	double G;          // 剪切模量
	double v;          // 泊松比
	double Rou;        // 密度
	double C1;         // 纵波波速
	double C2;         // 横波波速
	
	double Dt;         // 时间步长
	double Dt2;        // 2倍时间步长
	double Dt_2;
	double Dt_3;

	// 其他系数
	double C1_2;
	double C2_2;
	double OneOverC1_2;
	double OneOverC2_2;
	double OneOverC1_3;
	double OneOverC2_3;
	double C2_2OverC1_2;
	double C2_2OverC1_3;

	double C1Dt;
	double C2Dt;

	double OneOverFourPaiRou;

	double OneOverThree_1OC23_1OC13;
	double Half_1OC22_1OC12;

	double STCoef1;
	double STCoef2;
};

DynaMat BuildDynaMat(double o_v, double o_G, double o_Rou, double o_dt);


class DSquareElement;

// 边界条件类
class BoundaryValue
{
public:
	Wvector lu;
	Wvector lt;

	BoundaryValue(){;}
	~BoundaryValue(){;}

	int init(long NodeNum)
	{
		lu.initial(3 * NodeNum);
		lt.initial(3 * NodeNum);

		return 1;
	}

	// (已验证) 根据边界条件，把第n个节点的未知量/已知量信息自我转化到位移/面力信息上
	int Convert(int bdid,long n);

	// (已验证) 根据边界条件，把所有节点的未知量/已知量信息自我转化到位移/面力信息上
	int Convert(DSquareElement* m_DSE,long EleNum);

	// (已验证) 将第n个节点的局部坐标值转化为绝对坐标值
	int LocToAbs(long n,double* Trans)
	{
		double tv[3],atv[3];
		long i;

		i = n * 3;
		
		// lu -> u
		tv[0] = lu.b[i];
		tv[1] = lu.b[i+1];
		tv[2] = lu.b[i+2];

		atv[0] = Trans[0] * tv[0] + Trans[3] * tv[1] + Trans[6] * tv[2];
		atv[1] = Trans[1] * tv[0] + Trans[4] * tv[1] + Trans[7] * tv[2];
		atv[2] = Trans[2] * tv[0] + Trans[5] * tv[1] + Trans[8] * tv[2];

		lu.b[i] = atv[0];
		lu.b[i+1] = atv[1];
		lu.b[i+2] = atv[2];

		// lt -> t
		tv[0] = lt.b[i];
		tv[1] = lt.b[i+1];
		tv[2] = lt.b[i+2];

		atv[0] = Trans[0] * tv[0] + Trans[3] * tv[1] + Trans[6] * tv[2];
		atv[1] = Trans[1] * tv[0] + Trans[4] * tv[1] + Trans[7] * tv[2];
		atv[2] = Trans[2] * tv[0] + Trans[5] * tv[1] + Trans[8] * tv[2];

		lt.b[i] = atv[0];
		lt.b[i+1] = atv[1];
		lt.b[i+2] = atv[2];

		return 1;
	}

	// (已验证) 将所有节点的局部坐标值转化为绝对坐标值
	int AllLocToAbs(long NodeNum,double** Trans);

	// 将lu和lt数组中低于最大值10-6的数置0
	int ToZero(double error);

};

class DSquareElement
{
public:
	
	// --------------- 变量 ---------------

	// 几何信息
	
	Point* m_pointlist;
	long m_pointID[8];

	Point* m_nodelist;

	//————根据单元类型调整————   默认常值单元   需改
	long m_nodeID[8];

	// 边界条件代码

	int BCID;
	int DomainID;
	int LocalEleID;
	int SurfaceType;
	//BCID=123: 给定U1t,U2t,Un
	//BCID=126: 给定U1t,U2t,Tn
	//BCID=153: 给定U1t,T2t,Un
	//BCID=156: 给定U1t,T2t,Tn
	//BCID=423: 给定T1t,U2t,Un
	//BCID=426: 给定T1t,U2t,Tn
	//BCID=453: 给定T1t,T2t,Un
	//BCID=456: 给定T1t,T2t,Tn

	// 坐标系数      形函数+绝对坐标
	
	double m_CoordCoef[3][8];

	// 单元内多次调用变量

	////————根据单元类型调整————   默认常值单元    需改
	//double m_OnElementU1ij[8][8][9];                       // DU1ij奇异积分,用于MU1(0)
	//double m_OnElementT1ij[8][8][9];                      // DT1ij奇异积分,用于MT1(0)  静力学
	//double m_OnElementT2ij[8][8][9];                      // DT2ij奇异积分,用于MT2(1)
	

	double m_OnElementU1ij[8][8][9];                       // DU1ij奇异积分,用于MU1(0)
	double m_OnElementT1ij[8][8][9];                      // DT1ij奇异积分,用于MT1(0)  静力学
	double m_OnElementT2ij[8][8][9];                      // DT2ij奇异积分,用于MT2(1)

	 //————

	//————————————大规模算例占空间————————
	//     //  ———————— 改为积分时计算———————


	//double m_Jacobi[GAUSSPOINT2];                         // 雅克比
	//double m_normal[GAUSSPOINT2][3];                      // 法向量
	//Point m_intpt[GAUSSPOINT2];                           // 常规高斯积分点坐标

	// 分片积分多次调用变量 - 5000单元即对应 106 MB

	//double m_JacobiPW[NUMPW2][GAUSSPOINTPW2];              
	//double m_normalPW[NUMPW2][GAUSSPOINTPW2][3];  
	//Point m_intptPW[NUMPW2][GAUSSPOINTPW2];

	// 单元变量总字节数：8 * (43 + 64 * 9 * 4 + GAUSSPOINT2 * 7 + NUMPW2 * GAUSSPOINTPW2 * 7 + NUMPW2 * 3)
	// 如果5000单元，总字节数：1123MB (GAUSSPOINT2=100,NUMPW2=400,GAUSSPOINTPW2=9)
	//                          368MB (GAUSSPOINT2=100,NUMPW2=100,GAUSSPOINTPW2=9)
	//                          589MB (GAUSSPOINT2=100,NUMPW2=400,GAUSSPOINTPW2=4)
	// 
	// 
	//————————————————————//  ———————— 改为积分时计算———————


	// ————————略占空间，不做修改，暂存
	// 用于判断是否分片积分   
	Point m_rcenPW[NUMPW2];
	Point m_center;
	double m_radius,m_radiusPW;
	//标识是否为无限单元,及对应第几个无限单元
	long InfEleFlag;


	// 静态变量
	


	static double m_gap;
	static GaussNumeMatrix m_gnm;

	

	// 用于单元内奇异积分的辅助静态变量
	static double m_OnEleN[8][4][SINGAUSSPOINT2][3];
	static double m_OnEleJ[8][4][SINGAUSSPOINT2];
	static double m_OnEleR[8][4][SINGAUSSPOINT2][4];

	static double m_AssistUij[8][9];
	static double m_AssistTij[8][9];
	static double m_AssistCoef[8][9];

	static double m_StaticTij[8][8][9];

	// 计时静态变量
    static Counter m_counter;

	// 单元几何静态变量
	static quadinfo m_quadinfo;
	static FValue m_FValue;
	static SecInfo m_SecInfo;

	// 节点所属单元信息变量
	static long* BeEleID;              // 节点所属的单元ID
	static int* BeElePos;              // 节点所属的单元中节点在单元的编号
	static double** m_transmat;         // 节点局部坐标转换矩阵，从局部坐标到整体坐标转换

	static DynaMat m_DMat;
	
	static double InfCenter[3];//无限单元中心
	static int FlagDynaMethod;//是否进行稳定性处理 1否2是
	static int Flagproblem;//1静态 0动态
	// --------------- 函数 ---------------

	// 初始化
	
	DSquareElement() : BCID(0), DomainID(0), LocalEleID(0), SurfaceType(0) {;}
	~DSquareElement(){;}
	static void StaticInit(double gap,double o_v,double o_G,double o_Rou,double o_dt,int o_FlagDynaMethod);
	static void StaticGetBeInfo(long NodeNum,long EleNum,long** m_EleNID);
	static void StaticRelease(long NodeNum);
	void GeoInit(Point* pointlist,long* pointID);
	void PhyInit(Point* nodelist,long* nodeID,int boundaryID);
	void PthPhyInit(Point* nodelist, long* nodeID, int boundaryID);

	// Lagrange函数
	
	static double N(double s1,double s2,int ID);
	static void N(double s1,double s2,double *res)
	{
//		res[0] = 0.25*(1.0-s2)*(s1-1.0)*(s1+s2+1.0);
//		res[1] = 0.25*(1.0-s2)*(1.0+s1)*(s1-1.0-s2);
//		res[2] = 0.25*(1.0+s2)*(1.0+s1)*(s1-1.0+s2);
//		res[3] = 0.25*(1.0+s2)*(s1-1.0)*(s1+1.0-s2);
//		res[4] = 0.5*(s1-1.0)*(1.0+s1)*(s2-1.0);
//		res[5] = 0.5*(1.0-s2)*(1.0+s2)*(1.0+s1);
//		res[6] = 1.0-s1*s1-res[4];
//		res[7] = 1.0-s2*s2-res[5];
		double s1s1=0.25*s1*s1;
		double s1s2=0.25*s1*s2;
		double s2s2=0.25*s2*s2;
		double s12s2=s1s1*s2;
		double s1s22=s1*s2s2;
		double halfs1=0.5*s1;
		double halfs2=0.5*s2;
		res[0]=-0.25+s1s2+s1s1-s12s2+s2s2-s1s22;
		res[1]=-0.25-s1s2+s1s1-s12s2+s2s2+s1s22;
		res[2]=-0.25+s1s2+s2s2+s1s22+s1s1+s12s2;
		res[3]=-0.25-s1s2+s1s1+s12s2+s2s2-s1s22;
		res[4]=0.5-halfs2-s1s1+s12s2-s1s1+s12s2;
		res[5]=0.5+halfs1-s2s2-s1s22-s2s2-s1s22;
		res[6]=0.5+halfs2-s1s1-s12s2-s1s1-s12s2;
		res[7]=0.5-halfs1-s2s2+s1s22-s2s2+s1s22;
	}

	// 拉氏函数求导
	
	static double NDiff1(double s1,double s2,int ID);
	static void NDiff1(double s1,double s2,double *res)
	{
		res[0] = 0.25*(1.0-s2)*(s2+s1+s1);
		res[1] = 0.25*(1.0-s2)*(s1+s1-s2);
		res[2] = 0.5*s2+s1-res[0];
		res[3] = s1-0.5*s2-res[1];
		res[4] = s1*(s2-1.0);
		res[5] = 0.5*(1.0-s2*s2);
		res[6] = -s1-s1-res[4];
		res[7] = -res[5];
	}
	static double NDiff2(double s1,double s2,int ID);
	static void NDiff2(double s1,double s2,double *res)
	{
		res[0] = 0.25*(1.0-s1)*(s1+s2+s2);
		res[1] = 0.25*(1.0+s1)*(s2+s2-s1);
		res[2] = 0.5*s1+s2-res[0];
		res[3] = -0.5*s1+s2-res[1];
		res[4] = 0.5*(s1*s1-1.0);
		res[5] = -(1.0+s1)*s2;
		res[6] = -res[4];
		res[7] = -s2-s2-res[5];
	}
	static void NDiff11(double s1,double s2,double *res)
	{
		res[0] = 0.5*(1.0-s2);
		res[1] = 0.5*(1.0-s2);
		res[2] = 0.5*(1.0+s2);
		res[3] = 0.5*(1.0+s2);
		res[4] = s2-1.0;
		res[5] = 0.0;
		res[6] = -1.0-s2;
		res[7] = 0.0;
	}
	static void NDiff12(double s1,double s2,double *res)
	{
		res[0] = -0.5*(s1+s2)+0.25;
		res[1] = 0.5*(s2-s1)-0.25;
		res[2] = 0.5*(s1+s2)+0.25;
		res[3] = 0.5*(s1-s2)-0.25;
		res[4] = s1;
		res[5] = -s2;
		res[6] = -s1;
		res[7] = s2;
	}
	static void NDiff22(double s1,double s2,double *res)
	{
		res[0] = 0.5*(1.0-s1);
		res[1] = 0.5*(1.0+s1);
		res[2] = 0.5*(1.0+s1);
		res[3] = 0.5*(1.0-s1);
		res[4] = 0.0;
		res[5] = -1.0-s1;
		res[6] = 0.0;
		res[7] = s1-1.0;
	}

	// 由局部坐标获得绝对坐标

	double GetFromLocal(double s1,double s2,int ID)
	{
		return GetValueFromCoefs(m_CoordCoef[ID],s1,s2);
	}
	void GetFromLocal(double s1,double s2,Point& pt)
	{
		double s12,s22,s1s22,s12s2;
		s12 = s1*s1;
		s22 = s2*s2;
		s1s22 = s1*s22;
		s12s2 = s12*s2;

		pt.pt[0] = m_CoordCoef[0][0] + m_CoordCoef[0][1]*s1 + m_CoordCoef[0][2]*s2 + m_CoordCoef[0][3]*s1*s2 + m_CoordCoef[0][4]*s12 + m_CoordCoef[0][5]*s22 + m_CoordCoef[0][6]*s12s2 + m_CoordCoef[0][7]*s1s22;
		pt.pt[1] = m_CoordCoef[1][0] + m_CoordCoef[1][1]*s1 + m_CoordCoef[1][2]*s2 + m_CoordCoef[1][3]*s1*s2 + m_CoordCoef[1][4]*s12 + m_CoordCoef[1][5]*s22 + m_CoordCoef[1][6]*s12s2 + m_CoordCoef[1][7]*s1s22;
		pt.pt[2] = m_CoordCoef[2][0] + m_CoordCoef[2][1]*s1 + m_CoordCoef[2][2]*s2 + m_CoordCoef[2][3]*s1*s2 + m_CoordCoef[2][4]*s12 + m_CoordCoef[2][5]*s22 + m_CoordCoef[2][6]*s12s2 + m_CoordCoef[2][7]*s1s22;
	}
	
	// 距离函数r的求导
	
	double RDiffs1(double s1,double s2,int ID)
	{
		return m_CoordCoef[ID][1] + m_CoordCoef[ID][3]*s2 + 2.0*m_CoordCoef[ID][4]*s1 + 2.0*m_CoordCoef[ID][6]*s1*s2 + m_CoordCoef[ID][7]*s2*s2;
	}
	void RDiffs1(double s1,double s2,double* res)
	{
		double s1_2 = 2.0*s1;
		double s1s2_2 = 2.0*s1*s2;
		double s22 = s2*s2;
		res[0] = m_CoordCoef[0][1] + m_CoordCoef[0][3]*s2 + m_CoordCoef[0][4]*s1_2 + m_CoordCoef[0][6]*s1s2_2 + m_CoordCoef[0][7]*s22;
		res[1] = m_CoordCoef[1][1] + m_CoordCoef[1][3]*s2 + m_CoordCoef[1][4]*s1_2 + m_CoordCoef[1][6]*s1s2_2 + m_CoordCoef[1][7]*s22;
		res[2] = m_CoordCoef[2][1] + m_CoordCoef[2][3]*s2 + m_CoordCoef[2][4]*s1_2 + m_CoordCoef[2][6]*s1s2_2 + m_CoordCoef[2][7]*s22;
	}
	double RDiffs2(double s1,double s2,int ID)
	{
		return m_CoordCoef[ID][2] + m_CoordCoef[ID][3]*s1 + 2.0*m_CoordCoef[ID][5]*s2 + m_CoordCoef[ID][6]*s1*s1 + 2.0*m_CoordCoef[ID][7]*s1*s2;
	}
	void RDiffs2(double s1,double s2,double* res)
	{
		double s2_2 = 2.0*s2;
		double s1s2_2 = 2.0*s1*s2;
		double s12 = s1*s1;
		res[0] = m_CoordCoef[0][2] + m_CoordCoef[0][3]*s1 + m_CoordCoef[0][5]*s2_2 + m_CoordCoef[0][6]*s12 + m_CoordCoef[0][7]*s1s2_2;
		res[1] = m_CoordCoef[1][2] + m_CoordCoef[1][3]*s1 + m_CoordCoef[1][5]*s2_2 + m_CoordCoef[1][6]*s12 + m_CoordCoef[1][7]*s1s2_2;
		res[2] = m_CoordCoef[2][2] + m_CoordCoef[2][3]*s1 + m_CoordCoef[2][5]*s2_2 + m_CoordCoef[2][6]*s12 + m_CoordCoef[2][7]*s1s2_2;
	}
	double RDiffs11(double s1,double s2,int ID)
	{
		return 2.0*(m_CoordCoef[ID][4] + m_CoordCoef[ID][6]*s2);
	}
	void RDiffs11(double s1,double s2,double* res)
	{
		res[0] = 2.0*(m_CoordCoef[0][4] + m_CoordCoef[0][6]*s2);
		res[1] = 2.0*(m_CoordCoef[1][4] + m_CoordCoef[1][6]*s2);
		res[2] = 2.0*(m_CoordCoef[2][4] + m_CoordCoef[2][6]*s2);
	}
	double RDiffs12(double s1,double s2,int ID)
	{
		return m_CoordCoef[ID][3] + 2.0*(m_CoordCoef[ID][6]*s1 + m_CoordCoef[ID][7]*s2);
	}
	void RDiffs12(double s1,double s2,double* res)
	{
		res[0] = m_CoordCoef[0][3] + 2.0*(m_CoordCoef[0][6]*s1 + m_CoordCoef[0][7]*s2);
		res[1] = m_CoordCoef[1][3] + 2.0*(m_CoordCoef[1][6]*s1 + m_CoordCoef[1][7]*s2);
		res[2] = m_CoordCoef[2][3] + 2.0*(m_CoordCoef[2][6]*s1 + m_CoordCoef[2][7]*s2);
	}
	double RDiffs22(double s1,double s2,int ID)
	{
		return 2.0*(m_CoordCoef[ID][5] + m_CoordCoef[ID][7]*s1);
	}
	void RDiffs22(double s1,double s2,double* res)
	{
		res[0] = 2.0*(m_CoordCoef[0][5] + m_CoordCoef[0][7]*s1);
		res[1] = 2.0*(m_CoordCoef[1][5] + m_CoordCoef[1][7]*s1);
		res[2] = 2.0*(m_CoordCoef[2][5] + m_CoordCoef[2][7]*s1);
	}

	// 根据局部坐标输出切向量1、2和法向量
	
	void TanS1(double s1,double s2,double& t1,double& t2,double& t3)
	{
		double res[3];
		RDiffs1(s1,s2,res);
		double length;
		length = sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);
		t1 = res[0]/length;
		t2 = res[1]/length;
		t3 = res[2]/length;
	}
	void TanS2(double s1,double s2,double& t1,double& t2,double& t3)
	{
		double res[3];
		RDiffs2(s1,s2,res);
		double length;
		length = sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);
		t1 = res[0]/length;
		t2 = res[1]/length;
		t3 = res[2]/length;
	}
	void Normal(double s1,double s2,double& n1,double& n2,double& n3)
	{
		double r1[3],r2[3];
		RDiffs1(s1,s2,r1);
		RDiffs2(s1,s2,r2);
		n1 = r1[1]*r2[2]-r2[1]*r1[2];
		n2 = r1[2]*r2[0]-r2[2]*r1[0];
		n3 = r1[0]*r2[1]-r2[0]*r1[1];
		double temp = sqrt(n1*n1+n2*n2+n3*n3);
		n1 = n1/temp;
		n2 = n2/temp;
		n3 = n3/temp;
	}


	
	// 根据局部坐标求雅克比系数

	double Jacobi(double s1,double s2)
	{
		double ru[3],rv[3],E,F,G;
		RDiffs1(s1,s2,ru);
		RDiffs2(s1,s2,rv);
		E = ru[0]*ru[0]+ru[1]*ru[1]+ru[2]*ru[2];
		F = ru[0]*rv[0]+ru[1]*rv[1]+ru[2]*rv[2];
		G = rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2];
		return sqrt(fabs(E*G-F*F));
	}
	
	// 根据局部坐标求转换矩阵, 从局部坐标到整体坐标转换

	void TranslateMatrix(double s1,double s2,double* TransMat)
	{
		// 方向1：由TanS1得到
		// 方向2：由垂直和方向1正交得到
		TanS1(s1,s2,TransMat[0],TransMat[1],TransMat[2]);
		TanS2(s1,s2,TransMat[3],TransMat[4],TransMat[5]);
		TransMat[6] = TransMat[1]*TransMat[5]-TransMat[4]*TransMat[2];
		TransMat[7] = TransMat[3]*TransMat[2]-TransMat[0]*TransMat[5];
		TransMat[8] = TransMat[0]*TransMat[4]-TransMat[3]*TransMat[1];

		double temp = sqrt(TransMat[6]*TransMat[6]+TransMat[7]*TransMat[7]+TransMat[8]*TransMat[8]);

		TransMat[6] /= temp;
		TransMat[7] /= temp;
		TransMat[8] /= temp;

		TransMat[3] = TransMat[7]*TransMat[2]-TransMat[8]*TransMat[1];
		TransMat[4] = TransMat[8]*TransMat[0]-TransMat[6]*TransMat[2];
		TransMat[5] = TransMat[6]*TransMat[1]-TransMat[7]*TransMat[0];
	}
	
	// 根据源点和场点所在单元的局部坐标获取距离矢量

	void GetR(const Point& source,double s1,double s2,double* R)
	{
		//R[0]: r
		//R[1]: dr/dx
		//R[2]: dr/dy
		//R[3]: dr/dz
		Point temp;
		GetFromLocal(s1,s2,temp);
		R[1]=temp.pt[0]-source.pt[0];
		R[2]=temp.pt[1]-source.pt[1];
		R[3]=temp.pt[2]-source.pt[2];
		R[0]=sqrt(R[1]*R[1]+R[2]*R[2]+R[3]*R[3]);
	}
	void GetR(const Point& source,Point& intpt,double* R)
	{
		//R[0]: r
		//R[1]: dr/dx
		//R[2]: dr/dy
		//R[3]: dr/dz
		R[1]=intpt.pt[0]-source.pt[0];
		R[2]=intpt.pt[1]-source.pt[1];
		R[3]=intpt.pt[2]-source.pt[2];
		R[0]=sqrt(R[1]*R[1]+R[2]*R[2]+R[3]*R[3]);
	}
	
	// 获取基本解函数和积分
	
	void GetTij(double* T,double* R,double* n);
	int GetUij(double* T, double* R);
	int GetStaticUij(double* T, double* R);
	// 奇异积分用辅助变量A，B

	void GetAi(double theta,double* DrDs1,double* DrDs2,double* A)
	{
		A[1]=DrDs1[0]*cos(theta)+DrDs2[0]*sin(theta);
		A[2]=DrDs1[1]*cos(theta)+DrDs2[1]*sin(theta);
		A[3]=DrDs1[2]*cos(theta)+DrDs2[2]*sin(theta);
		A[0]=sqrt(A[1]*A[1]+A[2]*A[2]+A[3]*A[3]);
		A[4]=A[0]*A[0]*A[0];   //A^3
		A[5]=A[4]*A[0]*A[0];   //A^5
	}
	static void GetBeta(double& Beta,double* A)
	{
		Beta = 1.0/A[0];
	}
	// J*ni = (J*ni)(theta)[rou = 0] + (J*ni)'(theta)[rou = 0]*rou + O(rou^2)
	// herein Ji0 = (J*ni)(theta)[rou = 0]
	void GetJi0(int NodeID,double* J)
	{
		double s1 = m_quadinfo.m_LocalCoord[NodeID][0];
		double s2 = m_quadinfo.m_LocalCoord[NodeID][1];
		double tJ=Jacobi(s1,s2);
		double n[3];
		Normal(s1,s2,n[0],n[1],n[2]);
		J[0]=tJ*n[0];
		J[1]=tJ*n[1];
		J[2]=tJ*n[2];
	}
	void GetFijMinusOne(double* FijMinusOne,double* A,double* J,double coef1,double coef2);
	
	// 从物理点向几何点的物理量外插

	int Approximate(double* input,double* output);
	int Approximate9P(double* input,double* output);

	// 计算应力相关函数
	int ObtainStress(double s1,double s2,double stress[3][3],double disp[8][3],double trac[8][3]);
	int MainStress(double stress[3][3],double& stress1,double&stress2,double& stress3,double& mises,double& tresca);

	// 计算多项式系数以加快计算的函数

	static int GetCoef(double* m_coef,double* m_pvar)
	{
			// (1/2*V5+1/4*V4+1/4*V3-1/2*V7-1/4*V2-1/4*V1)*s1^2*s2
			// +(-1/2*V5+1/4*V4+1/4*V3-1/2*V7+1/4*V2+1/4*V1)*s1^2
			// +(1/2*V8-1/2*V6+1/4*V2-1/4*V4+1/4*V3-1/4*V1)*s1*s2^2
			// +(1/4*V1-1/4*V2+1/4*V3-1/4*V4)*s1*s2
			// +(1/2*V6-1/2*V8)*s1
			// +(-1/2*V6+1/4*V2-1/2*V8+1/4*V3+1/4*V1+1/4*V4)*s2^2
			// +(-1/2*V5+1/2*V7)*s2
			// -1/4*V1-1/4*V4+1/2*V6-1/4*V3-1/4*V2+1/2*V5+1/2*V8+1/2*V7

			// 对应常数项
			m_coef[0] = 0.25*(-m_pvar[0] - m_pvar[3] + 2*m_pvar[5] - m_pvar[2] - m_pvar[1] + 2*m_pvar[4] + 2*m_pvar[7] + 2*m_pvar[6]);
			// 对应s1
			m_coef[1] = 0.5*(m_pvar[5] - m_pvar[7]);
			// 对应s2
			m_coef[2] = 0.5*(-m_pvar[4] + m_pvar[6]);
			// 对应s1*s2
			m_coef[3] = 0.25*(m_pvar[0] - m_pvar[1] + m_pvar[2] - m_pvar[3]);
			// 对应s1^2
			m_coef[4] = 0.25*(-2*m_pvar[4] + m_pvar[3] + m_pvar[2] - 2*m_pvar[6] + m_pvar[1] + m_pvar[0]);
			// 对应s2^2
			m_coef[5] = 0.25*(-2*m_pvar[5] + m_pvar[1] - 2*m_pvar[7] + m_pvar[2] + m_pvar[0] + m_pvar[3]);
			// 对应s1^2*s2
			m_coef[6] = 0.25*(2*m_pvar[4] + m_pvar[3] + m_pvar[2] - 2*m_pvar[6] - m_pvar[1] - m_pvar[0]);
			// 对应s2^2*s1
			m_coef[7] = 0.25*(2*m_pvar[7] - 2*m_pvar[5] + m_pvar[1] - m_pvar[3] + m_pvar[2] - m_pvar[0]);

		return 1;
	}

	static double GetValueFromCoefs(double* m_coef,double s1,double s2)
	{
		double s12,s22;
		s12 = s1*s1;
		s22 = s2*s2;

		return m_coef[0] + m_coef[1]*s1 + m_coef[2]*s2 + m_coef[3]*s1*s2 + m_coef[4]*s12 + m_coef[5]*s22 + m_coef[6]*s12*s2 + m_coef[7]*s1*s22;
	}

	int IsIn(long NodeID,int& Pos);

	int GetCoordCoef();



	// 动力学函数

	// (已验证) 获取距离函数
	void GetR(const Point& source,Point& intpt,double* R,double* RI)
	{
		RI[1] = intpt.pt[0] - source.pt[0];
		RI[2] = intpt.pt[1] - source.pt[1];
		RI[3] = intpt.pt[2] - source.pt[2];
		

		R[1] = sqrt(RI[1] * RI[1] + RI[2] * RI[2] + RI[3] * RI[3]);
		R[2] = R[1] * R[1];
		R[3] = R[2] * R[1];
		R[4] = R[3] * R[1];
		R[5] = R[4] * R[1];

		RI[0] = R[1];
	}

	// (已验证: 与英文教科书一致) Aij
	int GetAij(double* aij,double* R,double* RI,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Aij = 1 / c1^2 * ri * rj / r^3

		double C1;
		C1 = m_DMat.OneOverC1_2 / R[3];

		aij[0] = C1 * RIRJ[0]; // RI[1] * RI[1];
		aij[1] = C1 * RIRJ[1]; // RI[2] * RI[1];
		aij[2] = C1 * RIRJ[2]; // RI[3] * RI[1];
		aij[3] = aij[1];
		aij[4] = C1 * RIRJ[3]; // RI[2] * RI[2];
		aij[5] = C1 * RIRJ[4]; // RI[3] * RI[2];
		aij[6] = aij[2];
		aij[7] = aij[5];
		aij[8] = C1 * RIRJ[5]; // RI[3] * RI[3];

		return 1;
	}

	// (已验证: 与英文教科书一致) Bij
	int GetBij(double* bij,double* R,double* RI,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Bij = 1 / c2^2 * (deltaij / r - ri * rj / r^3)

		bij[0] = m_DMat.OneOverC2_2 * (1.0 / R[1] - RIRJ[0] / R[3]);
		bij[1] = m_DMat.OneOverC2_2 * (- RIRJ[1] / R[3]);
		bij[2] = m_DMat.OneOverC2_2 * (- RIRJ[2] / R[3]);
		bij[3] = bij[1];
		bij[4] = m_DMat.OneOverC2_2 * (1.0 / R[1] - RIRJ[3] / R[3]);
		bij[5] = m_DMat.OneOverC2_2 * (- RIRJ[4] / R[3]);
		bij[6] = bij[2];
		bij[7] = bij[5];
		bij[8] = m_DMat.OneOverC2_2 * (1.0 / R[1] - RIRJ[5] / R[3]);

		return 1;
	}

	// (已验证: 与英文教科书一致) Cij
	int GetCij(double* cij,double* R,double* RI,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Cij = 3 * ri * rj / r^3 - deltaij / r

		double C1,C2;
		C1 = 3.0 / R[3];
		C2 = 1.0 / R[1];

		cij[0] = C1 * RIRJ[0] - C2;
		cij[1] = C1 * RIRJ[1];
		cij[2] = C1 * RIRJ[2];
		cij[3] = cij[1];
		cij[4] = C1 * RIRJ[3] - C2;
		cij[5] = C1 * RIRJ[4];
		cij[6] = cij[2];
		cij[7] = cij[5];
		cij[8] = C1 * RIRJ[5] - C2;

		return 1;
	}

	// (已验证: 与英文教科书一致, 前提ij交换) Dij
	int GetDij(double* dij,double* R,double* RI,double* Nor,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Dij = (2 * c2^2/c1^2 - 1) * ri * nj / r^3 - c2^2/c1^2 * (12 * ri * rj * rm * nm / r^5 - 2 * (ri * nj + rj * ni + deltaij * rm * nm) / r^3)
		//     = (4 * c2^2/c1^2 - 1) * ri * nj / r^3 - c2^2/c1^2 * (12 * ri * rj * rm * nm / r^5 - 2 * (rj * ni + deltaij * rm * nm) / r^3)
		
		double C1,C2,C3,RmNm;
		
		RmNm = RI[1] * Nor[0] + RI[2] * Nor[1] + RI[3] * Nor[2];
		C1 = (4.0 * m_DMat.C2_2OverC1_2 - 1.0) / R[3];
		C2 = m_DMat.C2_2OverC1_2 * 12.0 * RmNm / R[5];
		C3 = m_DMat.C2_2OverC1_2 * 2.0 / R[3];

		dij[0] = C1 * RI[1] * Nor[0] - (C2 * RIRJ[0] - C3 * (RI[1] * Nor[0] + RmNm));
		dij[1] = C1 * RI[2] * Nor[0] - (C2 * RIRJ[1] - C3 * (RI[1] * Nor[1]));
		dij[2] = C1 * RI[3] * Nor[0] - (C2 * RIRJ[2] - C3 * (RI[1] * Nor[2]));
		dij[3] = C1 * RI[1] * Nor[1] - (C2 * RIRJ[1] - C3 * (RI[2] * Nor[0]));
		dij[4] = C1 * RI[2] * Nor[1] - (C2 * RIRJ[3] - C3 * (RI[2] * Nor[1] + RmNm));
		dij[5] = C1 * RI[3] * Nor[1] - (C2 * RIRJ[4] - C3 * (RI[2] * Nor[2]));
		dij[6] = C1 * RI[1] * Nor[2] - (C2 * RIRJ[2] - C3 * (RI[3] * Nor[0]));
		dij[7] = C1 * RI[2] * Nor[2] - (C2 * RIRJ[4] - C3 * (RI[3] * Nor[1]));
		dij[8] = C1 * RI[3] * Nor[2] - (C2 * RIRJ[5] - C3 * (RI[3] * Nor[2] + RmNm));		

		return 1;
	}

	// (已验证: 与英文教科书一致, 前提ij交换) Eij
	int GetEij(double* eij,double* R,double* RI,double* Nor,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Eij = 12 * ri * rj * rm * nm / r^5 - 2 * ri * nj / r^3 - 3 * (rj * ni + deltaij * rm * nm) / r^3
		
		double C1,C2,C3,RmNm;
		
		RmNm = RI[1] * Nor[0] + RI[2] * Nor[1] + RI[3] * Nor[2];
		C1 = 12.0 * RmNm / R[5];
		C2 = 2.0 / R[3];
		C3 = 3.0 / R[3];

		eij[0] = C1 * RIRJ[0] - C2 * RI[1] * Nor[0] - C3 * (RI[1] * Nor[0] + RmNm);
		eij[1] = C1 * RIRJ[1] - C2 * RI[2] * Nor[0] - C3 * (RI[1] * Nor[1]);
		eij[2] = C1 * RIRJ[2] - C2 * RI[3] * Nor[0] - C3 * (RI[1] * Nor[2]);
		eij[3] = C1 * RIRJ[1] - C2 * RI[1] * Nor[1] - C3 * (RI[2] * Nor[0]);
		eij[4] = C1 * RIRJ[3] - C2 * RI[2] * Nor[1] - C3 * (RI[2] * Nor[1] + RmNm);
		eij[5] = C1 * RIRJ[4] - C2 * RI[3] * Nor[1] - C3 * (RI[2] * Nor[2]);
		eij[6] = C1 * RIRJ[2] - C2 * RI[1] * Nor[2] - C3 * (RI[3] * Nor[0]);
		eij[7] = C1 * RIRJ[4] - C2 * RI[2] * Nor[2] - C3 * (RI[3] * Nor[1]);
		eij[8] = C1 * RIRJ[5] - C2 * RI[3] * Nor[2] - C3 * (RI[3] * Nor[2] + RmNm);

		return 1;
	}

	// (已验证: 与英文教科书一致, 前提ij交换) Fij
	int GetFij(double* fij,double* R,double* RI,double* Nor,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Fij = -6 * c2^2 * (5 * ri * rj * rm * nm / r^5 - (ri * nj + rj * ni + deltaij * rm * nm) / r^3)
		
		double C1,C2,RmNm;
		
		RmNm = RI[1] * Nor[0] + RI[2] * Nor[1] + RI[3] * Nor[2];
		C1 = - 6.0 * m_DMat.C2_2;
		C2 = 5.0 * RmNm / R[5];

		fij[0] = C1 * (C2 * RIRJ[0] - (RI[1] * Nor[0] + RI[1] * Nor[0] + RmNm) / R[3]);
		fij[1] = C1 * (C2 * RIRJ[1] - (RI[2] * Nor[0] + RI[1] * Nor[1]) / R[3]);
		fij[2] = C1 * (C2 * RIRJ[2] - (RI[3] * Nor[0] + RI[1] * Nor[2]) / R[3]);
		fij[3] = fij[1]; // C1 * (C2 * RI[1] * RI[2] - (RI[1] * Nor[1] + RI[2] * Nor[0]) / R[3]);
		fij[4] = C1 * (C2 * RIRJ[3] - (RI[2] * Nor[1] + RI[2] * Nor[1] + RmNm) / R[3]);
		fij[5] = C1 * (C2 * RIRJ[4] - (RI[3] * Nor[1] + RI[2] * Nor[2]) / R[3]);
		fij[6] = fij[2]; // C1 * (C2 * RI[1] * RI[3] - (RI[1] * Nor[2] + RI[3] * Nor[0]) / R[3]);
		fij[7] = fij[5]; // C1 * (C2 * RI[2] * RI[3] - (RI[2] * Nor[2] + RI[3] * Nor[1]) / R[3]);
		fij[8] = C1 * (C2 * RIRJ[5] - (RI[3] * Nor[2] + RI[3] * Nor[2] + RmNm) / R[3]);

		return 1;
	}

	// (已验证: 与英文教科书一致, 前提ij交换) Gij
	int GetGij(double* gij,double* R,double* RI,double* Nor,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Gij = -2 * c2^2/c1^3 * ri * rj * rm * nm / r^4  + 1 / c1 * (2 * c2^2/c1^2 - 1) * ri * nj / r^2
		
		double C1,C2,RmNm;
		
		RmNm = RI[1] * Nor[0] + RI[2] * Nor[1] + RI[3] * Nor[2];
		C1 = - 2.0 * m_DMat.C2_2OverC1_3 * RmNm / R[4];
		C2 = 1.0 / m_DMat.C1 * (2.0 * m_DMat.C2_2OverC1_2 - 1.0) / R[2];

		gij[0] = C1 * RIRJ[0] + C2 * RI[1] * Nor[0];
		gij[1] = C1 * RIRJ[1] + C2 * RI[2] * Nor[0];
		gij[2] = C1 * RIRJ[2] + C2 * RI[3] * Nor[0];
		gij[3] = C1 * RIRJ[1] + C2 * RI[1] * Nor[1];
		gij[4] = C1 * RIRJ[3] + C2 * RI[2] * Nor[1];
		gij[5] = C1 * RIRJ[4] + C2 * RI[3] * Nor[1];
		gij[6] = C1 * RIRJ[2] + C2 * RI[1] * Nor[2];
		gij[7] = C1 * RIRJ[4] + C2 * RI[2] * Nor[2];
		gij[8] = C1 * RIRJ[5] + C2 * RI[3] * Nor[2];

		return 1;
	}

	// (已验证: 与英文教科书一致, 前提ij交换) Hij
	int GetHij(double* hij,double* R,double* RI,double* Nor,double* RIRJ)
	{
		// R[0]: blank
		// R[1]: r
		// R[2]: r^2;
		// R[3]: r^3;
		// R[4]: r^4
		// R[5]: r^5
		// RI[1]~RI[3]: r1,r2,r3

		// Hij = 2 / c2 * ri * rj * rm * nm / r^4  - 1 / c2 * (rj * ni + deltaij * rm * nm) / r^2
		
		double C1,C2,RmNm;
		
		RmNm = RI[1] * Nor[0] + RI[2] * Nor[1] + RI[3] * Nor[2];
		C1 = 2.0 / m_DMat.C2 * RmNm / R[4];
		C2 = - 1.0 / m_DMat.C2 / R[2];

		hij[0] = C1 * RIRJ[0] + C2 * (RI[1] * Nor[0] + RmNm);
		hij[1] = C1 * RIRJ[1] + C2 * RI[1] * Nor[1];
		hij[2] = C1 * RIRJ[2] + C2 * RI[1] * Nor[2];
		hij[3] = C1 * RIRJ[1] + C2 * RI[2] * Nor[0];
		hij[4] = C1 * RIRJ[3] + C2 * (RI[2] * Nor[1] + RmNm);
		hij[5] = C1 * RIRJ[4] + C2 * RI[2] * Nor[2];
		hij[6] = C1 * RIRJ[2] + C2 * RI[3] * Nor[0];
		hij[7] = C1 * RIRJ[4] + C2 * RI[3] * Nor[1];
		hij[8] = C1 * RIRJ[5] + C2 * (RI[3] * Nor[2] + RmNm);

		return 1;
	}

	// (已验证，与我推导公式一致)
	int Fai(double n,double dt,double r,double c)
	{
		// c: 膨胀波速或剪切波速
		double ndt = n * dt;
		double roverc = r / c;

		return ((roverc >= ndt) && (roverc < (ndt + dt))) ? 1 : 0;
	}

	int FuncI(double n,double dt,double* R,double& ResIn);

	int FuncJ(double n,double dt,double* R,double& ResJn);

	// (经验证与 FuncI(n+1) - FuncI(n) 相同) 计算I(n-1) - I(n)
	int ISeries(double n,double dt,double* R,double& ResISeries);

	// (经验证与 FuncJ(n+1) - FuncJ(n) 相同) 计算J(n-1) - J(n)
	int JSeries(double n,double dt,double* R,double& ResJSeries);

	// 判断是否需要计算
	int NeedCal(int type,double n,Point& o_Point);

	int PointNeedCal(int type,double n,Point& s_Point,Point& o_Point)
	{
		double nn = (type == 2 ? n - 1 : n);
		double r = Dist(s_Point,o_Point);
		return (r >= nn * m_DMat.C2Dt && r < (nn + 1) * m_DMat.C1Dt) ? 1 : 0;
	}

	int PWNeedCal(int type,double n,Point& PW_center,Point& o_Point)
	{
	// 该函数用来判断分片中心点与场点单元的距离是否满足计算动态系数的条件
	// 0: 没有收到波动影响，不参与计算
	// 1: 收到波动影响，调用PW计算
	
		double nn = (type == 2 ? n - 1 : n);

		double c2ndt = nn * m_DMat.C2Dt;
		double c1np1dt = (nn + 1) * m_DMat.C1Dt;
		double r = Dist(PW_center,o_Point);

		double rpr = r + m_radiusPW;
		double rmr = r - m_radiusPW;

		return (rmr >= c1np1dt || rpr < c2ndt) ? 0 : 1;
	}

	// // (已验证，与我推导公式一致) 计算动态基本解U
	int GetDynaUij(double* U,double n,double dt,double* R,double* RI);

	int GetDynaT1ij(double* T,double n,double dt,double* R,double* RI,double* Nor);
	int GetDynaT2ij(double* T,double n,double dt,double* R,double* RI,double* Nor);
	int GetDynaTij(int typeT,double* T,double n,double dt,double* R,double* RI,double* Nor);
	//静力学
	int IntUij(Point & source);
	int IntTij(Point & source);
	int IntTijforInfEle(Point & source,int type);
	int IntUijforInfEle();
	

	// (已验证) 计算动态基本解积分,得到的积分值存在公共变量m_AssistUij或者m_AssistTij，绝对坐标
	int IntDynaUij(Point& source,double n,double dt);
	int IntDynaUij(Point& source, double n, double dt, long T_ID);
	int IntDynaTij(int typeT,Point& source,double n,double dt);
	int IntDynaTij(Point& source,double n,double dt);

	int IntDynaUij(Point& source,double n,double dt,int ID);
	int IntDynaTij(int typeT,Point& source,double n,double dt,int ID);
	int IntDynaTij(Point& source,double n,double dt,int ID);
	int IntDynaTij(int typeT, Point& source, double n, double dt, long T_ID);

	// 分片积分方案
	int IntDynaUijPW(Point& source,double n,double dt);
	int IntDynaUijPW(Point& source, double n, double dt, long T_ID);
	int IntDynaTijPW(int typeT,Point& source,double n,double dt);
	int IntDynaTijPW(int typeT, Point& source, double n, double dt,long T_ID);
	int IntDynaTijPW(Point& source,double n,double dt);

	int IntDynaUijPW(Point& source,double n,double dt,int ID);
	int IntDynaTijPW(int typeT,Point& source,double n,double dt,int ID);
	int IntDynaTijPW(Point& source,double n,double dt,int ID);
	int IntDynaTijPWforInfEle(Point& source, double n, double dt, int type);
	int IntDynaTijPWforInfEle(int Ttype, Point& source, double n, double dt, int type);

	// 选择积分方案：根据距离判断是分片还是常规
	int IntDynaUijJudge(Point& source,double n,double dt);
	int IntDynaTijJudge(int typeT,Point& source,double n,double dt);
	int IntDynaTijJudge(Point& source,double n,double dt);
	int IntDynaTijJudge(int typeT, Point& source, double n, double dt, long T_ID);



	int IntDynaUijJudge(Point& source,double n,double dt,int ID);
	int IntDynaUijJudge(Point& source, double n, double dt, long T_ID);
	int IntDynaJudge(Point& source, double n, double dt, long T_ID);
	int IntDynaTijJudge(int typeT,Point& source,double n,double dt,int ID);
	int IntDynaTijJudge(Point& source,double n,double dt,int ID);
	int IntDynaTijJudge(Point& source, double n, double dt,long T_ID);

	// 计算n=0时的单元内奇异积分项，存在单元变量m_OnElementUij和m_OnElementTij，绝对坐标
	    
		// (已验证) 静态奇异基本解
	int InElementStaticTij();
	int InElementStaticTij(double* StaticTij);

		// (已逐行检查) 动态奇异基本解, PW解法
	int InElementIntDynamic();
	int PthInElementIntDynamic();
	//添加普通静力学部分
	int InElementIntStatic();
	int PthInElementIntStatic(long threadID);
	
	// (已验证与原有三维弹性力学程序一致) 调整已知量和未知量,仅限U1ij(n=0)和T1ij(n=0)
	
	int ConvertUnknownCoef(); 
	int ConvertUnknownCoef(long T_ID);
	int ConvertKnownCoef();
	int ConvertUnknownCoef(int ID);
	int ConvertKnownCoef(int ID);
	int ConvertKnownCoef(long T_ID);
	int ConvertUnknownCoefforInf();
	int ConvertKnownCoefforInf();

	
	// 分别进行已知量和未知量系数的计算,仅限U1ij(n=0)和T1ij(n=0)

		// (已逐行检查) 常规积分计算，包括坐标和已知未知量转换，适用于U1ij(n=0)和T1ij(n=0)
		// 可根据距离判断调用PW还是常规积分，返回0, 1, 2
	int RegularIntUnknownCoef(Point& source);
	int RegularIntUnknownCoef(Point& source,long T_ID);
	int RegularIntKnownCoef(Point& source);
	int RegularIntKnownCoef(Point& source, long T_ID);
	int RegularIntUnknownCoef(Point& source,int ID);
	int RegularIntKnownCoef(Point& source,int ID);
	int RegularIntUnknownCoefforStatic(Point & source);
	int IntUnknownCoefforInf(DSquareElement* m_elelist, long NodeID, long FieldEleID, long ** InfElePID, long** ElePID,int m_InfEleFlag);
	int IntKnownCoefforInf(DSquareElement* m_elelist, long NodeID, long FieldEleID, long ** InfElePID, long** ElePID, int m_InfEleFlag);
	int RegularIntKnownCoefforStatic(Point & source);
		// (已逐行检查) 奇异积分计算，包括坐标和已知未知量转换，U1ij(n=0)和T1ij(n=0)
	int InElementIntUnknownCoef(int sourceID);
	int InElementIntUnknownCoef(int sourceID,long T_ID);
	int InElementIntKnownCoef(int sourceID);
	int InElementIntUnknownCoef(int sourceID,int ID);
	int InElementIntKnownCoef(int sourceID,int ID);	
	int InElementIntKnownCoef(int sourceID, long T_ID);

};

// 矩阵系数运算

// 获取未知量子矩阵，仅适用于DU1ij(n=0)和DT1ij(n=0)
int UnknownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID);
int UnknownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID, long** InfElePID, long** ElePID);
int UnknownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID, long** InfElePID, long** ElePID,long T_ID);
int NodeUnknownSubMatrix(DSquareElement* m_elelist, long SourceNodeID, long FieldNodeID);

// 获取已知量子矩阵，仅适用于DU1ij(n=0)和DT1ij(n=0)
int knownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID);
int knownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID, long** InfElePID, long** ElePID);
int knownSubMatrix(DSquareElement* m_elelist, long NodeID, long FieldEleID, long** InfElePID, long** ElePID, long T_ID);
int NodeknownSubMatrix(DSquareElement* m_elelist, long SourceNodeID, long FieldNodeID);

//--------------------------------------------------------
#endif
