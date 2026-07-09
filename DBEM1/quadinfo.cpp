#include "quadinfo.h"

namespace
{
	double ConstThetaStart(int area)
	{
		switch (area)
		{
		case 0: return -0.75 * pi;
		case 1: return -0.25 * pi;
		case 2: return 0.25 * pi;
		case 3: return 0.75 * pi;
		default: return 0.0;
		}
	}

	double ConstThetaDiff()
	{
		return 0.5 * pi;
	}

	double ConstThetaFromLocal(double localTheta, int area)
	{
		return ConstThetaStart(area) + 0.5 * (localTheta + 1.0) * ConstThetaDiff();
	}

	double ConstAlphaTheta(double theta, int area)
	{
		switch (area)
		{
		case 0: return -1.0 / sin(theta);
		case 1: return 1.0 / cos(theta);
		case 2: return 1.0 / sin(theta);
		case 3: return -1.0 / cos(theta);
		default: return 0.0;
		}
	}

	double ConstPolarFromLocal(double localRou, double theta, int area)
	{
		return 0.5 * (localRou + 1.0) * ConstAlphaTheta(theta, area);
	}

	void ConstFromPolar(double& s1, double& s2, double rou, double theta)
	{
		s1 = rou * cos(theta);
		s2 = rou * sin(theta);
	}
}

//------------------------------------------------------------------------------

SecInfo::SecInfo()
{
	;
}

SecInfo::~SecInfo()
{
	;
}

FValue::FValue()
{
	;
}

FValue::~FValue()
{
	;
}

quadinfo::quadinfo()
{
	;
}

quadinfo::~quadinfo()
{
	;
}

int quadinfo::init(double gap)
{
	long i, j, m,  p;
	double res[8];

	m_gap = gap;
	GetLocalCoord(m_LocalCoord, gap);

	// 计算差值函数N和高斯权值的乘积
	for (m = 0; m < GAUSSPOINT2; ++m)
	{
		N(res, m_gnm.rp[m][0], m_gnm.rp[m][1]);
		for (p = 0; p < 8; ++p)
			m_NRGV[p][m] = res[p] * m_gnm.ra[m];
		m_RGV[m] = m_gnm.ra[m];
	}

	// 计算分片参数

	double LocLength = 2.0 / NUMPW;

	// 分片中心和分片高斯点坐标
	for (i = 0; i < NUMPW; ++i)
	{
		for (j = 0; j < NUMPW; ++j)
		{
			p = i * NUMPW + j; // 编号

			m_CENPW[p][0] = -1.0 + LocLength / 2.0 + LocLength * j;
			m_CENPW[p][1] = -1.0 + LocLength / 2.0 + LocLength * i;

			for (m = 0; m < GAUSSPOINTPW2; ++m)
			{
				m_LPOSPW[p][m][0] = m_CENPW[p][0] + m_gnm.pwrp[m][0] / NUMPW;
				m_LPOSPW[p][m][1] = m_CENPW[p][1] + m_gnm.pwrp[m][1] / NUMPW;
			}
		}
	}

	// 分片权重：N * 高斯值 * (LocLength / 2) ^ 2
	for (i = 0; i < NUMPW2; ++i)
	{
		for (m = 0; m < GAUSSPOINTPW2; ++m)
		{
			N(res, m_LPOSPW[i][m][0], m_LPOSPW[i][m][1]);
			m_NPWGV[i][0][m] = res[0] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_NPWGV[i][1][m] = res[1] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_NPWGV[i][2][m] = res[2] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_NPWGV[i][3][m] = res[3] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_NPWGV[i][4][m] = res[4] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_NPWGV[i][5][m] = res[5] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_NPWGV[i][6][m] = res[6] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_NPWGV[i][7][m] = res[7] * m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
			m_PWGV[i][m] = m_gnm.pwra[m] * (LocLength / 2.0) * (LocLength / 2.0);
		}
	}

	// 奇异积分计算分片参数

	LocLength = 2.0 / SINNUMPW;

	// 奇异积分分片中心和分片高斯点坐标
	for (i = 0; i < SINNUMPW; ++i)
	{
		for (j = 0; j < SINNUMPW; ++j)
		{
			p = i * SINNUMPW + j; // 编号

			m_SINCENPW[p][0] = -1.0 + LocLength / 2.0 + LocLength * j;
			m_SINCENPW[p][1] = -1.0 + LocLength / 2.0 + LocLength * i;

			for (m = 0; m < SINGAUSSPOINTPW2; ++m)
			{
				m_SINLPOSPW[p][m][0] = m_SINCENPW[p][0] + m_gnm.pwsp[m][0] / SINNUMPW;
				m_SINLPOSPW[p][m][1] = m_SINCENPW[p][1] + m_gnm.pwsp[m][1] / SINNUMPW;
			}
		}
	}

	// 奇异积分分片权重：高斯值 * (LocLength / 2) ^ 2
	for (m = 0; m < SINGAUSSPOINTPW2; ++m)
		m_SINNPWGV[m] = m_gnm.pwsa[m] * (LocLength / 2.0) * (LocLength / 2.0);


	// 给角度赋值
	AssignTheta();

	return 1;
}

int quadinfo::AssignTheta()
{
	int ID;

	for (ID = 0; ID < 8; ++ID)
	{
		m_CornerTheta[ID][0] = atan((m_LocalCoord[ID][1] + 1.0) / (m_LocalCoord[ID][0] + 1.0)) - pi;
		m_CornerTheta[ID][1] = atan((-1.0 - m_LocalCoord[ID][1]) / (1.0 - m_LocalCoord[ID][0]));
		m_CornerTheta[ID][2] = atan((1.0 - m_LocalCoord[ID][1]) / (1.0 - m_LocalCoord[ID][0]));
		m_CornerTheta[ID][3] = pi - atan((1.0 - m_LocalCoord[ID][1]) / (m_LocalCoord[ID][0] + 1.0));
	}




	// 计算角度差
	for (ID = 0; ID < 8; ++ID)
	{
		m_CornerThetaDiff[ID][0] = m_CornerTheta[ID][1] - m_CornerTheta[ID][0];
		m_CornerThetaDiff[ID][1] = m_CornerTheta[ID][2] - m_CornerTheta[ID][1];
		m_CornerThetaDiff[ID][2] = m_CornerTheta[ID][3] - m_CornerTheta[ID][2];
		m_CornerThetaDiff[ID][3] = m_CornerTheta[ID][0] - m_CornerTheta[ID][3] + 2.0*pi;
	}

	return 1;
}

int quadinfo::AssignSecInfo(SecInfo& m_SecInfo)
{

	int ID, AreaID;
	long i, m, n;
	double NValue[8];
	double s1, s2, rou, theta, temp, tempFA;

	// 计算面上的积分点信息和权值信息
	for (ID = 0; ID < 8; ++ID)
	{
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (n = 0; n < SINGAUSSPOINT2; ++n) // recycle for localtheta
			{
				theta = GetThetaFromLocalTheta(m_gnm.sp[n][0], ID, AreaID);
				tempFA = 0.25*FAlfpaTheta(theta, ID, AreaID)*m_CornerThetaDiff[ID][AreaID];

				rou = GetPolarFromLocalPolar(m_gnm.sp[n][1], theta, ID, AreaID);
				GetFromPolar(s1, s2, rou, theta, ID);
				m_SecInfo.m_SecPos[ID][AreaID][0][n] = s1;
				m_SecInfo.m_SecPos[ID][AreaID][1][n] = s2;
				temp = rou * tempFA*m_gnm.sa[n];
				N(NValue, s1, s2);
				for (i = 0; i < 8; ++i)
					m_SecInfo.m_SecVal[ID][i][AreaID][n] = NValue[i] * temp;
			}
		}
	}

	for (AreaID = 0; AreaID < 4; ++AreaID)
	{
		for (n = 0; n < SINGAUSSPOINT2; ++n)
		{
			theta = ConstThetaFromLocal(m_gnm.sp[n][0], AreaID);
			tempFA = 0.25 * ConstAlphaTheta(theta, AreaID) * ConstThetaDiff();
			rou = ConstPolarFromLocal(m_gnm.sp[n][1], theta, AreaID);
			ConstFromPolar(s1, s2, rou, theta);
			m_SecInfo.m_SecConstPos[AreaID][0][n] = s1;
			m_SecInfo.m_SecConstPos[AreaID][1][n] = s2;
			m_SecInfo.m_SecConstVal[AreaID][n] = rou * tempFA * m_gnm.sa[n];
		}
	}

	// 计算线上的积分点信息
	for (ID = 0; ID < 8; ++ID)
	{
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (n = 0; n < SINGAUSSPOINT; ++n)//recycle for localtheta
			{
				theta = GetThetaFromLocalTheta(m_gnm.spos[n], ID, AreaID);
				m_SecInfo.m_LinePos[ID][AreaID][n] = theta;
			}
		}
	}

	for (AreaID = 0; AreaID < 4; ++AreaID)
	{
		for (n = 0; n < SINGAUSSPOINT; ++n)
		{
			m_SecInfo.m_ConstLinePos[AreaID][n] = ConstThetaFromLocal(m_gnm.spos[n], AreaID);
		}
	}

	// 分片计算(非常可能出错，确认检查)

	for (ID = 0; ID < 8; ++ID)
	{
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (n = 0; n < SINNUMPW2; ++n)//recycle for localtheta
			{
				for (m = 0; m < SINGAUSSPOINTPW2; ++m)//recycle for localrou
				{
					theta = GetThetaFromLocalTheta(m_SINLPOSPW[n][m][0], ID, AreaID);
					tempFA = 0.25*FAlfpaTheta(theta, ID, AreaID) * m_CornerThetaDiff[ID][AreaID];

					rou = GetPolarFromLocalPolar(m_SINLPOSPW[n][m][1], theta, ID, AreaID);
					GetFromPolar(s1, s2, rou, theta, ID);
					m_SecInfo.m_SecPosPW[ID][AreaID][0][n][m] = s1;
					m_SecInfo.m_SecPosPW[ID][AreaID][1][n][m] = s2;
					temp = rou * tempFA * m_SINNPWGV[m];
					N(NValue, s1, s2);
					for (i = 0; i < 8; ++i)
						m_SecInfo.m_SecValPW[ID][i][AreaID][n][m] = NValue[i] * temp;
				}
			}
		}
	}

	for (AreaID = 0; AreaID < 4; ++AreaID)
	{
		for (n = 0; n < SINNUMPW2; ++n)
		{
			for (m = 0; m < SINGAUSSPOINTPW2; ++m)
			{
				theta = ConstThetaFromLocal(m_SINLPOSPW[n][m][0], AreaID);
				tempFA = 0.25 * ConstAlphaTheta(theta, AreaID) * ConstThetaDiff();
				rou = ConstPolarFromLocal(m_SINLPOSPW[n][m][1], theta, AreaID);
				ConstFromPolar(s1, s2, rou, theta);
				m_SecInfo.m_SecConstPosPW[AreaID][0][n][m] = s1;
				m_SecInfo.m_SecConstPosPW[AreaID][1][n][m] = s2;
				m_SecInfo.m_SecConstValPW[AreaID][n][m] = rou * tempFA * m_SINNPWGV[m];
			}
		}
	}

	return 1;
}

int quadinfo::AssignFValue(FValue& m_FValue)
{

	int ID, AreaID;
	long n;
	double rou, theta, temp, tempFA;

	// 计算面上的积分点信息和权值信息
	for (ID = 0; ID < 8; ++ID)
	{
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (n = 0; n < SINGAUSSPOINT2; ++n)//recycle for localtheta
			{
				theta = GetThetaFromLocalTheta(m_gnm.sp[n][0], ID, AreaID);
				tempFA = 0.25*FAlfpaTheta(theta, ID, AreaID)*m_CornerThetaDiff[ID][AreaID];

				rou = GetPolarFromLocalPolar(m_gnm.sp[n][1], theta, ID, AreaID);
				temp = rou * tempFA*m_gnm.sa[n];
				m_FValue.m_FMinusOneSecVal[ID][AreaID][n] = temp / rou / rou;
			}
		}
	}

	for (AreaID = 0; AreaID < 4; ++AreaID)
	{
		for (n = 0; n < SINGAUSSPOINT2; ++n)
		{
			theta = ConstThetaFromLocal(m_gnm.sp[n][0], AreaID);
			tempFA = 0.25 * ConstAlphaTheta(theta, AreaID) * ConstThetaDiff();
			rou = ConstPolarFromLocal(m_gnm.sp[n][1], theta, AreaID);
			temp = rou * tempFA * m_gnm.sa[n];
			m_FValue.m_ConstFMinusOneSecVal[AreaID][n] = temp / rou / rou;
		}
	}

	// 计算线上的积分点信息
	for (ID = 0; ID < 8; ++ID)
	{
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (n = 0; n < SINGAUSSPOINT; ++n)//recycle for localtheta
			{
				theta = GetThetaFromLocalTheta(m_gnm.spos[n], ID, AreaID);
				rou = GetPolarFromLocalPolar(1.0, theta, ID, AreaID);
				m_FValue.m_FMinusOneLineVal[ID][AreaID][n] = 0.5*m_CornerThetaDiff[ID][AreaID] * log(rou)*m_gnm.sb[n];

				// 根据Guiggiani的书
				// 对单一单元，节点完全在单元内时，m_FMinusOneLineVal_2无须计算；
				// 对节点属于多个单元时，m_FMinusOneLineVal_2必须计算，与log(beta)有关；
				m_FValue.m_FMinusOneLineVal_2[ID][AreaID][n] = -0.5*m_CornerThetaDiff[ID][AreaID] * m_gnm.sb[n];

			}
		}
	}

	for (AreaID = 0; AreaID < 4; ++AreaID)
	{
		for (n = 0; n < SINGAUSSPOINT; ++n)
		{
			theta = ConstThetaFromLocal(m_gnm.spos[n], AreaID);
			rou = ConstPolarFromLocal(1.0, theta, AreaID);
			m_FValue.m_ConstFMinusOneLineVal[AreaID][n] = 0.5 * ConstThetaDiff() * log(rou) * m_gnm.sb[n];
			m_FValue.m_ConstFMinusOneLineVal_2[AreaID][n] = -0.5 * ConstThetaDiff() * m_gnm.sb[n];
		}
	}

	return 1;
}

double quadinfo::FAlfpaTheta(double theta, int NodeID, int AreaID)
{
	//返回 rou~(theta),theta在角度的积分点上
	//适用于0面积的三角形
	switch (AreaID)
	{
	case 0:
		return (m_LocalCoord[NodeID][1] + 1.0) / sin(pi + theta);
	case 1:
		return (1.0 - m_LocalCoord[NodeID][0]) / cos(theta);
	case 2:
		return (1.0 - m_LocalCoord[NodeID][1]) / sin(theta);
	case 3:
		return (m_LocalCoord[NodeID][0] + 1.0) / cos(pi - theta);
	default:
		return 0.0;
	}
}
