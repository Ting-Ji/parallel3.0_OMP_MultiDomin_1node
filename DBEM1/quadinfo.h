#ifndef _QUADINFO_
#define _QUADINFO_
//---------------------------------------------------------------------------------------
#include "precompiler.h"
#include "gaussnumematrix.h"
#include "math.h"

class SecInfo
{
public:
	// 用于将笛卡尔坐标积分转换到角坐标系下的积分点坐标，共8个物理点，每个物理点对应4个三角形
	double m_SecPos[8][4][2][SINGAUSSPOINT2];
	// 用于将笛卡尔坐标积分转换到角坐标系下的积分权值（物理量直接和该值相乘再乘以雅克比系数即可得到积分值），8个物理点，8*8
	double m_SecVal[8][8][4][SINGAUSSPOINT2];
	// 用于角坐标系下在线坐标上的积分点，共8个物理点，每个物理点对应4个三角形
	double m_LinePos[8][4][SINGAUSSPOINT];

	// 分片定义
	double m_SecPosPW[8][4][2][SINNUMPW2][SINGAUSSPOINTPW2];
	double m_SecValPW[8][8][4][SINNUMPW2][SINGAUSSPOINTPW2];

	SecInfo();
	~SecInfo();
};

class FValue
{
	// 专门用于强奇异积分时在角坐标系下与F(-1)ij相乘用的权值函数
public:
	//values of 1/rou*f(theta)/2*(theta2-theta1)/2*hij on gauss points of 4 areas(polar coordinates)
	double m_FMinusOneSecVal[8][4][SINGAUSSPOINT2];
	//gauss values of log(|rou~~(theta)|)*0.5*(theta2-theta1) on four lines
	double m_FMinusOneLineVal[8][4][SINGAUSSPOINT];
	//gauss values of 0.5*(theta2-theta1) on four lines, used for log(beta)
	double m_FMinusOneLineVal_2[8][4][SINGAUSSPOINT];

	FValue();
	~FValue();
};

class quadinfo
{
public:
	// 变量
	double m_gap;
	double m_LocalCoord[8][2];

	GaussNumeMatrix m_gnm;

	//multiplication of N and Gauss value at Gauss Point
	double m_NRGV[8][GAUSSPOINT2];               // 用于辅助笛卡尔坐标系下的常规高斯积分

	double m_NPWGV[NUMPW2][8][GAUSSPOINTPW2];    // 用于辅助笛卡尔坐标下下的分片高斯积分
	double m_CENPW[NUMPW2][2];                   // 分片中心局部坐标
	double m_LPOSPW[NUMPW2][GAUSSPOINTPW2][2];   // 分片高斯坐标局部坐标，对应在-1~1的范围的分片值

	double m_SINNPWGV[SINGAUSSPOINTPW2];                  // 用于辅助笛卡尔坐标下下的分片高斯积分，不含插值函数N，只有高斯权重*(L/2)^2
	double m_SINCENPW[SINNUMPW2][2];                      // 分片中心局部坐标
	double m_SINLPOSPW[SINNUMPW2][SINGAUSSPOINTPW2][2];   // 分片高斯坐标局部坐标，对应在-1~1的范围的分片值

	double m_CornerTheta[8][4];
	double m_CornerThetaDiff[8][4];

	// 函数
	quadinfo();
	~quadinfo();

	int init(double gap);

	static int GetLocalCoord(double LocalCoord[][2],double& gap)
	{
		LocalCoord[0][0] = -gap;
		LocalCoord[0][1] = -gap;

		LocalCoord[1][0] = gap;
		LocalCoord[1][1] = -gap;

		LocalCoord[2][0] = gap;
		LocalCoord[2][1] = gap;

		LocalCoord[3][0] = -gap;
		LocalCoord[3][1] = gap;

		LocalCoord[4][0] = 0.0;
		LocalCoord[4][1] = -gap;

		LocalCoord[5][0] = gap;
		LocalCoord[5][1] = 0.0;

		LocalCoord[6][0] = 0.0;
		LocalCoord[6][1] = gap;

		LocalCoord[7][0] = -gap;
		LocalCoord[7][1] = 0.0;
		
		return 1;
	}

	int AssignTheta();

	int AssignSecInfo(SecInfo& m_SecInfo);

	int AssignFValue(FValue& m_FValue);

	double FAlfpaTheta(double theta,int NodeID,int AreaID);

	double GetThetaFromLocalTheta(double LocalTheta,int NodeID,int AreaID)
	{
		//localtheta: -1~1
		//result: 在该三角形内，从初始角到结束角之间的角，根据localtheta确定
		//适用于0面积的三角形
		return m_CornerTheta[NodeID][AreaID] + 0.5*(LocalTheta + 1.0)*m_CornerThetaDiff[NodeID][AreaID];
	}

	double GetPolarFromLocalPolar(double localrou,double theta,int NodeID,int AreaID)
	{
		//localrou: -1~1
		//result: 在该三角形内，从0到边之间的rou，根据localrou确定
		//适用于0面积的三角形
		return 0.5*(localrou+1.0)*FAlfpaTheta(theta,NodeID,AreaID);
	}
	
	void GetPolCoordFromLocPolCoord(double &rou,double &theta,double localrou,double localtheta,int NodeID,int AreaID)
	{
		// 在某个三角形内，有局部的-1~1分布的theta和rou计算得到实际的在三角形内的绝对坐标theta和rou
		//适用于0面积的三角形
		theta = GetThetaFromLocalTheta(localtheta,NodeID,AreaID);
		rou = GetPolarFromLocalPolar(localrou,theta,NodeID,AreaID);

	}
	void GetFromPolar(double &s1,double &s2,double rou,double theta,int NodeID)
	{
		// 根据三角形内的绝对极坐标rou和theta，计算得到对应的笛卡尔坐标s1和s2

		s1 = m_LocalCoord[NodeID][0]+rou*cos(theta);
		s2 = m_LocalCoord[NodeID][1]+rou*sin(theta);
	}

	int N(double* res,double s1,double s2)
	{
		double lx = m_LocalCoord[1][0] - m_LocalCoord[0][0];
		double ly = m_LocalCoord[3][1] - m_LocalCoord[0][1];
		double lxly = lx*ly;
		double hx = 0.25*lx*lxly;
		double hy = 0.25*ly*lxly;
		
		res[0] = (s1 - m_LocalCoord[1][0]) * (s2 - m_LocalCoord[3][1])/lxly;
		res[1] = (s1 - m_LocalCoord[0][0]) * (s2 - m_LocalCoord[2][1])/(-lxly);
		res[2] = (s1 - m_LocalCoord[3][0]) * (s2 - m_LocalCoord[1][1])/(lxly);
		res[3] = (s1 - m_LocalCoord[2][0]) * (s2 - m_LocalCoord[0][1])/(-lxly);
		res[4] = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
		res[5] = (s1 - m_LocalCoord[7][0])*(s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
		res[6] = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
		res[7] = (s1 - m_LocalCoord[5][0])*(s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);

		res[0] -= 0.5*(res[4] + res[7]);
		res[1] -= 0.5*(res[4] + res[5]);
		res[2] -= 0.5*(res[5] + res[6]);
		res[3] -= 0.5*(res[6] + res[7]);

		return 1;
	}

	int NDiffs1(double* res,double s1,double s2)
	{
		double lx = m_LocalCoord[1][0] - m_LocalCoord[0][0];
		double ly = m_LocalCoord[3][1] - m_LocalCoord[0][1];
		double lxly = lx*ly;
		double hx = 0.25*lx*lxly;
		double hy = 0.25*ly*lxly;
		
		res[0] = (s2 - m_LocalCoord[3][1])/(lxly);
		res[1] = (s2 - m_LocalCoord[2][1])/(-lxly);
		res[2] = (s2 - m_LocalCoord[1][1])/(lxly);
		res[3] = (s2 - m_LocalCoord[0][1])/(-lxly);
		res[4] = (2.0*s1 - m_LocalCoord[0][0] - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
		res[5] = (s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
		res[6] = (2.0*s1 - m_LocalCoord[2][0] - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
		res[7] = (s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);

		res[0] -= 0.5*(res[4] + res[7]);
		res[1] -= 0.5*(res[4] + res[5]);
		res[2] -= 0.5*(res[5] + res[6]);
		res[3] -= 0.5*(res[6] + res[7]);

		return 1;
	}

	int NDiffs2(double* res,double s1,double s2)
	{
		double lx = m_LocalCoord[1][0] - m_LocalCoord[0][0];
		double ly = m_LocalCoord[3][1] - m_LocalCoord[0][1];
		double lxly = lx*ly;
		double hx = 0.25*lx*lxly;
		double hy = 0.25*ly*lxly;
		
		res[0] = (s1 - m_LocalCoord[1][0])/lxly;
		res[1] = (s1 - m_LocalCoord[0][0])/(-lxly);
		res[2] = (s1 - m_LocalCoord[3][0])/(lxly);
		res[3] = (s1 - m_LocalCoord[2][0])/(-lxly);
		res[4] = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])/(hx);
		res[5] = (s1 - m_LocalCoord[7][0])*(2.0*s2 - m_LocalCoord[1][1] - m_LocalCoord[2][1])/(-hy);
		res[6] = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])/(-hx);
		res[7] = (s1 - m_LocalCoord[5][0])*(2.0*s2 - m_LocalCoord[0][1] - m_LocalCoord[3][1])/(hy);

		res[0] -= 0.5*(res[4] + res[7]);
		res[1] -= 0.5*(res[4] + res[5]);
		res[2] -= 0.5*(res[5] + res[6]);
		res[3] -= 0.5*(res[6] + res[7]);

		return 1;
	}

	int N(double& res,double s1,double s2,int ID)
	{
		double lx = m_LocalCoord[1][0] - m_LocalCoord[0][0];
		double ly = m_LocalCoord[3][1] - m_LocalCoord[0][1];
		double hx = 0.25*lx*lx*ly;
		double hy = 0.25*ly*ly*lx;
		double t1,t2;

		switch(ID)
			{
			case 0:
				{
					res = (s1 - m_LocalCoord[1][0])/(-lx) * (s2 - m_LocalCoord[3][1])/(-ly);
					t1 = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
					t2 = (s1 - m_LocalCoord[5][0])*(s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);
					res -= 0.5*(t1 + t2);
					break;
				}
			case 1:
				{
					res = (s1 - m_LocalCoord[0][0])/(lx) * (s2 - m_LocalCoord[2][1])/(-ly);
					t1 = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
					t2 = (s1 - m_LocalCoord[7][0])*(s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
					res -= 0.5*(t1 + t2);
					break;
				}
			case 2:
				{
					res = (s1 - m_LocalCoord[3][0])/(lx) * (s2 - m_LocalCoord[1][1])/(ly);
					t1 = (s1 - m_LocalCoord[7][0])*(s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
					t2 = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
					res -= 0.5*(t1 + t2);
					break;
				}
			case 3:
				{
					res = (s1 - m_LocalCoord[2][0])/(-lx) * (s2 - m_LocalCoord[0][1])/(ly);
					t1 = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
					t2 = (s1 - m_LocalCoord[5][0])*(s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);
					res -= 0.5*(t1 + t2);
					break;
				}
			case 4:
				{
					res = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
					break;
				}
			case 5:
				{
					res = (s1 - m_LocalCoord[7][0])*(s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
					break;
				}
			case 6:
				{
					res = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
					break;
				}
			case 7:
				{
					res = (s1 - m_LocalCoord[5][0])*(s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);
					break;
				}
			}

		return 1;
	}

	int NDiffs1(double& res,double s1,double s2,int ID)
	{
		double lx = m_LocalCoord[1][0] - m_LocalCoord[0][0];
		double ly = m_LocalCoord[3][1] - m_LocalCoord[0][1];
		double lxly = lx*ly;
		double hx = 0.25*lx*lxly;
		double hy = 0.25*ly*lxly;
		double t1,t2;

		switch(ID)
		{
		case 0:
			{
				res = (s2 - m_LocalCoord[3][1])/(lxly);
				t1 = (2.0*s1 - m_LocalCoord[0][0] - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
				t2 = (s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);
				res -= 0.5*(t1 + t2);
				break;
			}
		case 1:
			{
				res = (s2 - m_LocalCoord[2][1])/(-lxly);
				t1 = (2.0*s1 - m_LocalCoord[0][0] - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
				t2 = (s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
				res -= 0.5*(t1 + t2);
				break;
			}
		case 2:
			{
				res = (s2 - m_LocalCoord[1][1])/(lxly);
				t1 = (s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
				t2 = (2.0*s1 - m_LocalCoord[2][0] - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
				res -= 0.5*(t1 + t2);
				break;
			}
		case 3:
			{
				res = (s2 - m_LocalCoord[0][1])/(-lxly);
				t1 = (2.0*s1 - m_LocalCoord[2][0] - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
				t2 = (s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);
				res -= 0.5*(t1 + t2);
				break;
			}
		case 4:
			{
				res = (2.0*s1 - m_LocalCoord[0][0] - m_LocalCoord[1][0])*(s2 - m_LocalCoord[6][1])/(hx);
				break;
			}
		case 5:
			{
				res = (s2 - m_LocalCoord[1][1])*(s2 - m_LocalCoord[2][1])/(-hy);
				break;
			}
		case 6:
			{
				res = (2.0*s1 - m_LocalCoord[2][0] - m_LocalCoord[3][0])*(s2 - m_LocalCoord[4][1])/(-hx);
				break;
			}
		case 7:
			{
				res = (s2 - m_LocalCoord[0][1])*(s2 - m_LocalCoord[3][1])/(hy);
				break;
			}
		}

		return 1;
	}

	int NDiffs2(double& res,double s1,double s2,int ID)
	{
		double lx = m_LocalCoord[1][0] - m_LocalCoord[0][0];
		double ly = m_LocalCoord[3][1] - m_LocalCoord[0][1];
		double lxly = lx*ly;
		double hx = 0.25*lx*lxly;
		double hy = 0.25*ly*lxly;
		double t1,t2;

		switch(ID)
		{
		case 0:
			{
				res = (s1 - m_LocalCoord[1][0])/lxly;
				t1 = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])/(hx);
				t2 = (s1 - m_LocalCoord[5][0])*(2.0*s2 - m_LocalCoord[0][1] - m_LocalCoord[3][1])/(hy);

				res -= 0.5*(t1 + t2);
				break;
			}
		case 1:
			{
				res = (s1 - m_LocalCoord[0][0])/(-lxly);
				t1 = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])/(hx);
				t2 = (s1 - m_LocalCoord[7][0])*(2.0*s2 - m_LocalCoord[1][1] - m_LocalCoord[2][1])/(-hy);
				res -= 0.5*(t1 + t2);
				break;
			}
		case 2:
			{
				res = (s1 - m_LocalCoord[3][0])/(lxly);
				t1 = (s1 - m_LocalCoord[7][0])*(2.0*s2 - m_LocalCoord[1][1] - m_LocalCoord[2][1])/(-hy);
				t2 = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])/(-hx);
				res -= 0.5*(t1 + t2);
				break;
			}
		case 3:
			{
				res = (s1 - m_LocalCoord[2][0])/(-lxly);
				t1 = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])/(-hx);
				t2 = (s1 - m_LocalCoord[5][0])*(2.0*s2 - m_LocalCoord[0][1] - m_LocalCoord[3][1])/(hy);
				res -= 0.5*(t1 + t2);
				break;
			}
		case 4:
			{
				res = (s1 - m_LocalCoord[0][0])*(s1 - m_LocalCoord[1][0])/(hx);
				break;
			}
		case 5:
			{
				res = (s1 - m_LocalCoord[7][0])*(2.0*s2 - m_LocalCoord[1][1] - m_LocalCoord[2][1])/(-hy);
				break;
			}
		case 6:
			{
				res = (s1 - m_LocalCoord[2][0])*(s1 - m_LocalCoord[3][0])/(-hx);
				break;
			}
		case 7:
			{
				res = (s1 - m_LocalCoord[5][0])*(2.0*s2 - m_LocalCoord[0][1] - m_LocalCoord[3][1])/(hy);
				break;
			}
		}

		return 1;
	}
};

//---------------------------------------------------------------------------------------
#endif
