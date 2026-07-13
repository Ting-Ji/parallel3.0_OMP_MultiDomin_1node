#include "DSquareElement.h"
#include "lin_comp.h"
#include "stdio.h"
#include "DBEM.h"

//realization of class: DSquareElement dyna version


int DSquareElement::InElementIntDynamic()
{
	// 计算n=0时的单元奇异积分，只有n=0时才有奇异积分

	double ST[9], temp, s1, s2;
	double DU1[9], DT1[9], DT2[9];
	int i, j, k, l, AreaID;

	double R[6], RI[4];

	double Jac, Nor[3], SR[6];

	// 首先计算静态Tij奇异积分
	InElementStaticTij();

	// 计算DUij的奇异积分，在角坐标系下完成

	for (i = 0; i < 8; ++i)
		for (j = 0; j < 8; ++j)
		{
			m_OnElementU1ij[i][j][0] = 0.0;
			m_OnElementU1ij[i][j][1] = 0.0;
			m_OnElementU1ij[i][j][2] = 0.0;
			m_OnElementU1ij[i][j][4] = 0.0;
			m_OnElementU1ij[i][j][5] = 0.0;
			m_OnElementU1ij[i][j][8] = 0.0;
		}

	for (i = 0; i < 8; ++i)// recycle for source point
	{
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINNUMPW2; ++k)
			{
				for (l = 0; l < SINGAUSSPOINTPW2; ++l)
				{
					s1 = m_SecInfo.m_SecConstPosPW[AreaID][0][k][l];
					s2 = m_SecInfo.m_SecConstPosPW[AreaID][1][k][l];
					Jac = Jacobi(s1, s2);
					GetR(m_nodelist[m_nodeID[i]], s1, s2, SR);

					R[1] = SR[0];
					R[2] = R[1] * R[1];
					R[3] = R[2] * R[1];
					R[4] = R[3] * R[1];
					R[5] = R[4] * R[1];
					RI[1] = SR[1];
					RI[2] = SR[2];
					RI[3] = SR[3];

					GetDynaUij(DU1, 0.0, ActiveDynaMat().Dt, R, RI);

					for (j = 0; j < 8; ++j)// recycle for field point
					{
						temp = m_SecInfo.m_SecConstValPW[AreaID][k][l] * Jac;

						m_OnElementU1ij[i][j][0] += DU1[0] * temp;
						m_OnElementU1ij[i][j][1] += DU1[1] * temp;
						m_OnElementU1ij[i][j][2] += DU1[2] * temp;
						m_OnElementU1ij[i][j][4] += DU1[4] * temp;
						m_OnElementU1ij[i][j][5] += DU1[5] * temp;
						m_OnElementU1ij[i][j][8] += DU1[8] * temp;
					}

				}
			}
		}
	}

	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			m_OnElementU1ij[i][j][3] = m_OnElementU1ij[i][j][1];
			m_OnElementU1ij[i][j][6] = m_OnElementU1ij[i][j][2];
			m_OnElementU1ij[i][j][7] = m_OnElementU1ij[i][j][5];
		}
	}


	// 计算DTij的奇异积分，也在角坐标系下完成

	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			m_OnElementT1ij[i][j][0] = 0.0;
			m_OnElementT1ij[i][j][1] = 0.0;
			m_OnElementT1ij[i][j][2] = 0.0;
			m_OnElementT1ij[i][j][3] = 0.0;
			m_OnElementT1ij[i][j][4] = 0.0;
			m_OnElementT1ij[i][j][5] = 0.0;
			m_OnElementT1ij[i][j][6] = 0.0;
			m_OnElementT1ij[i][j][7] = 0.0;
			m_OnElementT1ij[i][j][8] = 0.0;

			m_OnElementT2ij[i][j][0] = 0.0;
			m_OnElementT2ij[i][j][1] = 0.0;
			m_OnElementT2ij[i][j][2] = 0.0;
			m_OnElementT2ij[i][j][3] = 0.0;
			m_OnElementT2ij[i][j][4] = 0.0;
			m_OnElementT2ij[i][j][5] = 0.0;
			m_OnElementT2ij[i][j][6] = 0.0;
			m_OnElementT2ij[i][j][7] = 0.0;
			m_OnElementT2ij[i][j][8] = 0.0;
		}
	}

	for (i = 0; i < 8; ++i)// recycle for source point
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINNUMPW2; ++k)
			{
				for (l = 0; l < SINGAUSSPOINTPW2; ++l)
				{
					s1 = m_SecInfo.m_SecConstPosPW[AreaID][0][k][l];
					s2 = m_SecInfo.m_SecConstPosPW[AreaID][1][k][l];
					Jac = Jacobi(s1, s2);
					Normal(s1, s2, Nor[0], Nor[1], Nor[2]);
					GetR(m_nodelist[m_nodeID[i]], s1, s2, SR);

					R[1] = SR[0];
					R[2] = R[1] * R[1];
					R[3] = R[2] * R[1];
					R[4] = R[3] * R[1];
					R[5] = R[4] * R[1];
					RI[1] = SR[1];
					RI[2] = SR[2];
					RI[3] = SR[3];

					GetDynaT1ij(DT1, 0.0, ActiveDynaMat().Dt, R, RI, Nor);
					GetDynaT2ij(DT2, 1.0, ActiveDynaMat().Dt, R, RI, Nor);

					GetTij(ST, SR, Nor);

					for (j = 0; j < 8; ++j)// recycle for field point
					{
						temp = m_SecInfo.m_SecConstValPW[AreaID][k][l] * Jac;
						if (i == j)
						{
							m_OnElementT1ij[i][j][0] += (DT1[0] - ST[0]) * temp;
							m_OnElementT1ij[i][j][1] += (DT1[1] - ST[1]) * temp;
							m_OnElementT1ij[i][j][2] += (DT1[2] - ST[2]) * temp;
							m_OnElementT1ij[i][j][3] += (DT1[3] - ST[3]) * temp;
							m_OnElementT1ij[i][j][4] += (DT1[4] - ST[4]) * temp;
							m_OnElementT1ij[i][j][5] += (DT1[5] - ST[5]) * temp;
							m_OnElementT1ij[i][j][6] += (DT1[6] - ST[6]) * temp;
							m_OnElementT1ij[i][j][7] += (DT1[7] - ST[7]) * temp;
							m_OnElementT1ij[i][j][8] += (DT1[8] - ST[8]) * temp;

							m_OnElementT2ij[i][j][0] += (DT2[0] - ST[0]) * temp;
							m_OnElementT2ij[i][j][1] += (DT2[1] - ST[1]) * temp;
							m_OnElementT2ij[i][j][2] += (DT2[2] - ST[2]) * temp;
							m_OnElementT2ij[i][j][3] += (DT2[3] - ST[3]) * temp;
							m_OnElementT2ij[i][j][4] += (DT2[4] - ST[4]) * temp;
							m_OnElementT2ij[i][j][5] += (DT2[5] - ST[5]) * temp;
							m_OnElementT2ij[i][j][6] += (DT2[6] - ST[6]) * temp;
							m_OnElementT2ij[i][j][7] += (DT2[7] - ST[7]) * temp;
							m_OnElementT2ij[i][j][8] += (DT2[8] - ST[8]) * temp;
						}
						else
						{
							m_OnElementT1ij[i][j][0] += DT1[0] * temp;
							m_OnElementT1ij[i][j][1] += DT1[1] * temp;
							m_OnElementT1ij[i][j][2] += DT1[2] * temp;
							m_OnElementT1ij[i][j][3] += DT1[3] * temp;
							m_OnElementT1ij[i][j][4] += DT1[4] * temp;
							m_OnElementT1ij[i][j][5] += DT1[5] * temp;
							m_OnElementT1ij[i][j][6] += DT1[6] * temp;
							m_OnElementT1ij[i][j][7] += DT1[7] * temp;
							m_OnElementT1ij[i][j][8] += DT1[8] * temp;

							m_OnElementT2ij[i][j][0] += DT2[0] * temp;
							m_OnElementT2ij[i][j][1] += DT2[1] * temp;
							m_OnElementT2ij[i][j][2] += DT2[2] * temp;
							m_OnElementT2ij[i][j][3] += DT2[3] * temp;
							m_OnElementT2ij[i][j][4] += DT2[4] * temp;
							m_OnElementT2ij[i][j][5] += DT2[5] * temp;
							m_OnElementT2ij[i][j][6] += DT2[6] * temp;
							m_OnElementT2ij[i][j][7] += DT2[7] * temp;
							m_OnElementT2ij[i][j][8] += DT2[8] * temp;

						}
					}
				}
			}
		}

	// 叠加静态Tij
	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			if (i == j)
			{
				m_OnElementT1ij[i][j][0] += m_StaticTij[i][j][0];
				m_OnElementT1ij[i][j][1] += m_StaticTij[i][j][1];
				m_OnElementT1ij[i][j][2] += m_StaticTij[i][j][2];
				m_OnElementT1ij[i][j][3] += m_StaticTij[i][j][3];
				m_OnElementT1ij[i][j][4] += m_StaticTij[i][j][4];
				m_OnElementT1ij[i][j][5] += m_StaticTij[i][j][5];
				m_OnElementT1ij[i][j][6] += m_StaticTij[i][j][6];
				m_OnElementT1ij[i][j][7] += m_StaticTij[i][j][7];
				m_OnElementT1ij[i][j][8] += m_StaticTij[i][j][8];

				m_OnElementT2ij[i][j][0] += m_StaticTij[i][j][0];
				m_OnElementT2ij[i][j][1] += m_StaticTij[i][j][1];
				m_OnElementT2ij[i][j][2] += m_StaticTij[i][j][2];
				m_OnElementT2ij[i][j][3] += m_StaticTij[i][j][3];
				m_OnElementT2ij[i][j][4] += m_StaticTij[i][j][4];
				m_OnElementT2ij[i][j][5] += m_StaticTij[i][j][5];
				m_OnElementT2ij[i][j][6] += m_StaticTij[i][j][6];
				m_OnElementT2ij[i][j][7] += m_StaticTij[i][j][7];
				m_OnElementT2ij[i][j][8] += m_StaticTij[i][j][8];
			}
		}
	}

	// MT1(0)叠加常数项

	for (i = 0; i < 8; ++i)
	{
		m_OnElementT1ij[i][i][0] += 0.5;
		m_OnElementT1ij[i][i][4] += 0.5;
		m_OnElementT1ij[i][i][8] += 0.5;
	}

	return 1;

}
int DSquareElement::PthInElementIntDynamic()
{
	// 计算n=0时的单元奇异积分，只有n=0时才有奇异积分


	double ST[9], temp, s1, s2;
	double DU1[9], DT1[9], DT2[9];
	int i, j, k, l, AreaID;

	double R[6], RI[4];

	double Jac, Nor[3], SR[6];

	// 首先计算静态Tij奇异积分
	InElementStaticTij();

	// 计算DUij的奇异积分，在角坐标系下完成

	for (i = 0; i < 8; ++i)
		for (j = 0; j < 8; ++j)
		{
			m_OnElementU1ij[i][j][0] = 0.0;
			m_OnElementU1ij[i][j][1] = 0.0;
			m_OnElementU1ij[i][j][2] = 0.0;
			m_OnElementU1ij[i][j][4] = 0.0;
			m_OnElementU1ij[i][j][5] = 0.0;
			m_OnElementU1ij[i][j][8] = 0.0;
		}

	for (i = 0; i < 8; ++i)// recycle for source point
	{
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINNUMPW2; ++k)
			{
				for (l = 0; l < SINGAUSSPOINTPW2; ++l)
				{
					s1 = m_SecInfo.m_SecConstPosPW[AreaID][0][k][l];
					s2 = m_SecInfo.m_SecConstPosPW[AreaID][1][k][l];
					Jac = Jacobi(s1, s2);
					GetR(m_nodelist[m_nodeID[i]], s1, s2, SR);

					R[1] = SR[0];
					R[2] = R[1] * R[1];
					R[3] = R[2] * R[1];
					R[4] = R[3] * R[1];
					R[5] = R[4] * R[1];
					RI[1] = SR[1];
					RI[2] = SR[2];
					RI[3] = SR[3];

					GetDynaUij(DU1, 0.0, ActiveDynaMat().Dt, R, RI);

					for (j = 0; j < 8; ++j)// recycle for field point
					{
						temp = m_SecInfo.m_SecConstValPW[AreaID][k][l] * Jac;

						m_OnElementU1ij[i][j][0] += DU1[0] * temp;
						m_OnElementU1ij[i][j][1] += DU1[1] * temp;
						m_OnElementU1ij[i][j][2] += DU1[2] * temp;
						m_OnElementU1ij[i][j][4] += DU1[4] * temp;
						m_OnElementU1ij[i][j][5] += DU1[5] * temp;
						m_OnElementU1ij[i][j][8] += DU1[8] * temp;
					}

				}
			}
		}
	}

	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			m_OnElementU1ij[i][j][3] = m_OnElementU1ij[i][j][1];
			m_OnElementU1ij[i][j][6] = m_OnElementU1ij[i][j][2];
			m_OnElementU1ij[i][j][7] = m_OnElementU1ij[i][j][5];
		}
	}


	// 计算DTij的奇异积分，也在角坐标系下完成

	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			m_OnElementT1ij[i][j][0] = 0.0;
			m_OnElementT1ij[i][j][1] = 0.0;
			m_OnElementT1ij[i][j][2] = 0.0;
			m_OnElementT1ij[i][j][3] = 0.0;
			m_OnElementT1ij[i][j][4] = 0.0;
			m_OnElementT1ij[i][j][5] = 0.0;
			m_OnElementT1ij[i][j][6] = 0.0;
			m_OnElementT1ij[i][j][7] = 0.0;
			m_OnElementT1ij[i][j][8] = 0.0;

			m_OnElementT2ij[i][j][0] = 0.0;
			m_OnElementT2ij[i][j][1] = 0.0;
			m_OnElementT2ij[i][j][2] = 0.0;
			m_OnElementT2ij[i][j][3] = 0.0;
			m_OnElementT2ij[i][j][4] = 0.0;
			m_OnElementT2ij[i][j][5] = 0.0;
			m_OnElementT2ij[i][j][6] = 0.0;
			m_OnElementT2ij[i][j][7] = 0.0;
			m_OnElementT2ij[i][j][8] = 0.0;
		}
	}

	for (i = 0; i < 8; ++i)// recycle for source point
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINNUMPW2; ++k)
			{
				for (l = 0; l < SINGAUSSPOINTPW2; ++l)
				{
					s1 = m_SecInfo.m_SecConstPosPW[AreaID][0][k][l];
					s2 = m_SecInfo.m_SecConstPosPW[AreaID][1][k][l];
					Jac = Jacobi(s1, s2);
					Normal(s1, s2, Nor[0], Nor[1], Nor[2]);
					GetR(m_nodelist[m_nodeID[i]], s1, s2, SR);

					R[1] = SR[0];
					R[2] = R[1] * R[1];
					R[3] = R[2] * R[1];
					R[4] = R[3] * R[1];
					R[5] = R[4] * R[1];
					RI[1] = SR[1];
					RI[2] = SR[2];
					RI[3] = SR[3];

					GetDynaT1ij(DT1, 0.0, ActiveDynaMat().Dt, R, RI, Nor);
					GetDynaT2ij(DT2, 1.0, ActiveDynaMat().Dt, R, RI, Nor);

					GetTij(ST, SR, Nor);

					for (j = 0; j < 8; ++j)// recycle for field point
					{
						temp = m_SecInfo.m_SecConstValPW[AreaID][k][l] * Jac;
						if (i == j)
						{
							m_OnElementT1ij[i][j][0] += (DT1[0] - ST[0]) * temp;
							m_OnElementT1ij[i][j][1] += (DT1[1] - ST[1]) * temp;
							m_OnElementT1ij[i][j][2] += (DT1[2] - ST[2]) * temp;
							m_OnElementT1ij[i][j][3] += (DT1[3] - ST[3]) * temp;
							m_OnElementT1ij[i][j][4] += (DT1[4] - ST[4]) * temp;
							m_OnElementT1ij[i][j][5] += (DT1[5] - ST[5]) * temp;
							m_OnElementT1ij[i][j][6] += (DT1[6] - ST[6]) * temp;
							m_OnElementT1ij[i][j][7] += (DT1[7] - ST[7]) * temp;
							m_OnElementT1ij[i][j][8] += (DT1[8] - ST[8]) * temp;

							m_OnElementT2ij[i][j][0] += (DT2[0] - ST[0]) * temp;
							m_OnElementT2ij[i][j][1] += (DT2[1] - ST[1]) * temp;
							m_OnElementT2ij[i][j][2] += (DT2[2] - ST[2]) * temp;
							m_OnElementT2ij[i][j][3] += (DT2[3] - ST[3]) * temp;
							m_OnElementT2ij[i][j][4] += (DT2[4] - ST[4]) * temp;
							m_OnElementT2ij[i][j][5] += (DT2[5] - ST[5]) * temp;
							m_OnElementT2ij[i][j][6] += (DT2[6] - ST[6]) * temp;
							m_OnElementT2ij[i][j][7] += (DT2[7] - ST[7]) * temp;
							m_OnElementT2ij[i][j][8] += (DT2[8] - ST[8]) * temp;
						}
						else
						{
							m_OnElementT1ij[i][j][0] += DT1[0] * temp;
							m_OnElementT1ij[i][j][1] += DT1[1] * temp;
							m_OnElementT1ij[i][j][2] += DT1[2] * temp;
							m_OnElementT1ij[i][j][3] += DT1[3] * temp;
							m_OnElementT1ij[i][j][4] += DT1[4] * temp;
							m_OnElementT1ij[i][j][5] += DT1[5] * temp;
							m_OnElementT1ij[i][j][6] += DT1[6] * temp;
							m_OnElementT1ij[i][j][7] += DT1[7] * temp;
							m_OnElementT1ij[i][j][8] += DT1[8] * temp;

							m_OnElementT2ij[i][j][0] += DT2[0] * temp;
							m_OnElementT2ij[i][j][1] += DT2[1] * temp;
							m_OnElementT2ij[i][j][2] += DT2[2] * temp;
							m_OnElementT2ij[i][j][3] += DT2[3] * temp;
							m_OnElementT2ij[i][j][4] += DT2[4] * temp;
							m_OnElementT2ij[i][j][5] += DT2[5] * temp;
							m_OnElementT2ij[i][j][6] += DT2[6] * temp;
							m_OnElementT2ij[i][j][7] += DT2[7] * temp;
							m_OnElementT2ij[i][j][8] += DT2[8] * temp;

						}
					}
				}
			}
		}

	// 叠加静态Tij
	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			if (i == j)
			{
			/*	m_OnElementT1ij[i][j][0] += m_StaticTij[i][j][0];
				m_OnElementT1ij[i][j][1] += m_StaticTij[i][j][1];
				m_OnElementT1ij[i][j][2] += m_StaticTij[i][j][2];
				m_OnElementT1ij[i][j][3] += m_StaticTij[i][j][3];
				m_OnElementT1ij[i][j][4] += m_StaticTij[i][j][4];
				m_OnElementT1ij[i][j][5] += m_StaticTij[i][j][5];
				m_OnElementT1ij[i][j][6] += m_StaticTij[i][j][6];
				m_OnElementT1ij[i][j][7] += m_StaticTij[i][j][7];
				m_OnElementT1ij[i][j][8] += m_StaticTij[i][j][8];

				m_OnElementT2ij[i][j][0] += m_StaticTij[i][j][0];
				m_OnElementT2ij[i][j][1] += m_StaticTij[i][j][1];
				m_OnElementT2ij[i][j][2] += m_StaticTij[i][j][2];
				m_OnElementT2ij[i][j][3] += m_StaticTij[i][j][3];
				m_OnElementT2ij[i][j][4] += m_StaticTij[i][j][4];
				m_OnElementT2ij[i][j][5] += m_StaticTij[i][j][5];
				m_OnElementT2ij[i][j][6] += m_StaticTij[i][j][6];
				m_OnElementT2ij[i][j][7] += m_StaticTij[i][j][7];
				m_OnElementT2ij[i][j][8] += m_StaticTij[i][j][8];*/

				m_OnElementT1ij[i][j][0] += m_StaticTij[i][j][0];
				m_OnElementT1ij[i][j][1] += m_StaticTij[i][j][1];
				m_OnElementT1ij[i][j][2] += m_StaticTij[i][j][2];
				m_OnElementT1ij[i][j][3] += m_StaticTij[i][j][3];
				m_OnElementT1ij[i][j][4] += m_StaticTij[i][j][4];
				m_OnElementT1ij[i][j][5] += m_StaticTij[i][j][5];
				m_OnElementT1ij[i][j][6] += m_StaticTij[i][j][6];
				m_OnElementT1ij[i][j][7] += m_StaticTij[i][j][7];
				m_OnElementT1ij[i][j][8] += m_StaticTij[i][j][8];

				m_OnElementT2ij[i][j][0] += m_StaticTij[i][j][0];
				m_OnElementT2ij[i][j][1] += m_StaticTij[i][j][1];
				m_OnElementT2ij[i][j][2] += m_StaticTij[i][j][2];
				m_OnElementT2ij[i][j][3] += m_StaticTij[i][j][3];
				m_OnElementT2ij[i][j][4] += m_StaticTij[i][j][4];
				m_OnElementT2ij[i][j][5] += m_StaticTij[i][j][5];
				m_OnElementT2ij[i][j][6] += m_StaticTij[i][j][6];
				m_OnElementT2ij[i][j][7] += m_StaticTij[i][j][7];
				m_OnElementT2ij[i][j][8] += m_StaticTij[i][j][8];
			}
		}
	}

	// MT1(0)叠加常数项

	for (i = 0; i < 8; ++i)
	{
		m_OnElementT1ij[i][i][0] += 0.5;
		m_OnElementT1ij[i][i][4] += 0.5;
		m_OnElementT1ij[i][i][8] += 0.5;
	}

	return 1;

	return 1;
}
//静力学使用
int DSquareElement::InElementIntStatic()
{
	// 计算单元奇异积分

	double  temp, s1, s2;
	double DU1[9];
	int i, j, k, AreaID;

	double  RI[4];

	double Jac;

	// 首先计算静态Tij奇异积分
	InElementStaticTij();
	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			m_OnElementT1ij[i][j][0] = 0.0;
			m_OnElementT1ij[i][j][1] = 0.0;
			m_OnElementT1ij[i][j][2] = 0.0;
			m_OnElementT1ij[i][j][3] = 0.0;
			m_OnElementT1ij[i][j][4] = 0.0;
			m_OnElementT1ij[i][j][5] = 0.0;
			m_OnElementT1ij[i][j][6] = 0.0;
			m_OnElementT1ij[i][j][7] = 0.0;
			m_OnElementT1ij[i][j][8] = 0.0;
		}
	}
	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{

			m_OnElementT1ij[i][j][0] += m_StaticTij[i][j][0];
			m_OnElementT1ij[i][j][1] += m_StaticTij[i][j][1];
			m_OnElementT1ij[i][j][2] += m_StaticTij[i][j][2];
			m_OnElementT1ij[i][j][3] += m_StaticTij[i][j][3];
			m_OnElementT1ij[i][j][4] += m_StaticTij[i][j][4];
			m_OnElementT1ij[i][j][5] += m_StaticTij[i][j][5];
			m_OnElementT1ij[i][j][6] += m_StaticTij[i][j][6];
			m_OnElementT1ij[i][j][7] += m_StaticTij[i][j][7];
			m_OnElementT1ij[i][j][8] += m_StaticTij[i][j][8];
		}
	}

	//叠加常数项

	for (i = 0; i < 8; ++i)
	{
		m_OnElementT1ij[i][i][0] += 0.5;
		m_OnElementT1ij[i][i][4] += 0.5;
		m_OnElementT1ij[i][i][8] += 0.5;
	}
	
	// 计算静态Uij的奇异积分，在角坐标系下完成
	for (i = 0; i < 8; ++i)
		for (j = 0; j < 8; ++j)
		{
			m_OnElementU1ij[i][j][0] = 0.0;
			m_OnElementU1ij[i][j][1] = 0.0;
			m_OnElementU1ij[i][j][2] = 0.0;
			m_OnElementU1ij[i][j][4] = 0.0;
			m_OnElementU1ij[i][j][5] = 0.0;
			m_OnElementU1ij[i][j][8] = 0.0;
		}
	//
	

	for (i = 0; i < 8; ++i)//recycle for source point
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINGAUSSPOINT2; ++k)
			{
				s1 = m_SecInfo.m_SecPos[i][AreaID][0][k];
				s2 = m_SecInfo.m_SecPos[i][AreaID][1][k];
				Jac = Jacobi(s1, s2);
				GetR(m_nodelist[m_nodeID[i]], s1, s2, RI);
				
				GetStaticUij(DU1,RI);// m_OnEleR[i][AreaID][k] 

				for (j = 0; j < 8; ++j)//recycle for field point
				{
					temp = m_SecInfo.m_SecVal[i][j][AreaID][k] * Jac;//m_OnEleJ[i][AreaID][k]
					m_OnElementU1ij[i][j][0] += DU1[0] * temp;
					m_OnElementU1ij[i][j][1] += DU1[1] * temp;
					m_OnElementU1ij[i][j][2] += DU1[2] * temp;
					m_OnElementU1ij[i][j][4] += DU1[4] * temp;
					m_OnElementU1ij[i][j][5] += DU1[5] * temp;
					m_OnElementU1ij[i][j][8] += DU1[8] * temp;
				}
			}
		}

	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < 8; ++j)
		{
			m_OnElementU1ij[i][j][3] = m_OnElementU1ij[i][j][1];
			m_OnElementU1ij[i][j][6] = m_OnElementU1ij[i][j][2];
			m_OnElementU1ij[i][j][7] = m_OnElementU1ij[i][j][5];
		}
	}


	return 1;

}

int DSquareElement::PthInElementIntStatic(long threadID)
{
	//未处理    暂留
	
	
	
	// 计算本单元静态奇异积分Tij，存于m_StaticTij公共变量中
	double coef1, coef2, UT[9], temp, s1, s2;
	double FVS1, FVL1;
	int i, j, k, l, AreaID;

	//获取8个物理点，每个物理点对应4个三角形的共8*4*SINGULARGAUSSPOINT*SINGULARGAUSSPOINT个积分点上的J、normal、R
	for (i = 0; i < 8; ++i)
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT2; ++j)
			{
				s1 = m_SecInfo.m_SecPos[i][AreaID][0][j];
				s2 = m_SecInfo.m_SecPos[i][AreaID][1][j];
				m_OnEleJ[i][AreaID][j] = Jacobi(s1, s2);
				Normal(s1, s2, m_OnEleN[i][AreaID][j][0], m_OnEleN[i][AreaID][j][1], m_OnEleN[i][AreaID][j][2]);
				GetR(m_nodelist[m_nodeID[i]], s1, s2, m_OnEleR[i][AreaID][j]);
			}
		}

	// 计算Tij的奇异积分：对于强奇异，需要在角坐标系下加一项减一项；对于弱奇异，只需要在角坐标系下积分即可

	coef1 = -1.0 / (8.0 * pi * (1.0 - ActiveDynaMat().v));
	coef2 = 1.0 - 2.0 * ActiveDynaMat().v;
	for (i = 0; i < 8; ++i)
		for (j = 0; j < 8; ++j)
			for (k = 0; k < 9; ++k)
				m_StaticTij[i][j][k] = 0.0;

	// 强奇异

	//variables for assissant of FijMinusOne
	double m_FijMinusOne[9], A[6], Jh[3], DrDs1[3], DrDs2[3];
	double Beta;   // 对多个单元共享一个节点的情况要用到
	for (i = 0; i < 8; ++i)
	{
		s1 = m_quadinfo.m_LocalCoord[i][0];
		s2 = m_quadinfo.m_LocalCoord[i][1];
		GetJi0(i, Jh);
		RDiffs1(s1, s2, DrDs1);
		RDiffs2(s1, s2, DrDs2);

		//Tij
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (j = 0; j < SINGAUSSPOINT2; ++j)
			{					//Tij
				GetTij(UT, m_OnEleR[i][AreaID][j], m_OnEleN[i][AreaID][j]);
				temp = m_SecInfo.m_SecVal[i][i][AreaID][j] * m_OnEleJ[i][AreaID][j];
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
				GetAi(m_SecInfo.m_LinePos[i][AreaID][j], DrDs1, DrDs2, A);

				GetBeta(Beta, A);
				Beta = fabs(Beta);

				GetFijMinusOne(m_FijMinusOne, A, Jh, coef1, coef2);
				for (l = 0; l < SINGAUSSPOINT; ++l)
				{
					FVS1 = m_FValue.m_FMinusOneSecVal[i][AreaID][l * SINGAUSSPOINT + j];  // [l][j] TEST
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
				FVL1 = m_FValue.m_FMinusOneLineVal[i][AreaID][j]; // 连续单元用： + m_FValue.m_FMinusOneLineVal_2[i][AreaID][j]*log(Beta);
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

	for (i = 0; i < 8; ++i)// recycle for source point
		for (AreaID = 0; AreaID < 4; ++AreaID)
		{
			for (k = 0; k < SINGAUSSPOINT2; ++k)
			{
				GetTij(UT, m_OnEleR[i][AreaID][k], m_OnEleN[i][AreaID][k]);

				for (j = 0; j < 8; ++j)// recycle for field point
				{
					if (i != j)
					{

						temp = m_SecInfo.m_SecVal[i][j][AreaID][k] * m_OnEleJ[i][AreaID][k];
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
	return 1;
}

int DSquareElement::IntDynaUijPW(Point& source, double n, double dt)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp = 0;
	int i, j, PointID;

	for (i = 0; i < 8; ++i)
	{
		m_AssistUij[i][0] = 0.0;
		m_AssistUij[i][1] = 0.0;
		m_AssistUij[i][2] = 0.0;
		m_AssistUij[i][4] = 0.0;
		m_AssistUij[i][5] = 0.0;
		m_AssistUij[i][8] = 0.0;
	}

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(1, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				//GetR(source, m_intptPW[i][j], R, RI);
				GetDynaUij(UT, n, dt, R, RI);

				for (PointID = 0; PointID < 8; ++PointID)
				{
					//temp = m_quadinfo.m_NPWGV[i][PointID][j] * m_JacobiPW[i][j];
					m_AssistUij[PointID][0] += UT[0] * temp;
					m_AssistUij[PointID][1] += UT[1] * temp;
					m_AssistUij[PointID][2] += UT[2] * temp;
					m_AssistUij[PointID][4] += UT[4] * temp;
					m_AssistUij[PointID][5] += UT[5] * temp;
					m_AssistUij[PointID][8] += UT[8] * temp;
				}
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

int DSquareElement::IntDynaUijPW(Point& source, double n, double dt, long T_ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp;
	int i, j, PointID;
	Point inputPW;
	double JacobiPW;

	for (i = 0; i < 8; ++i)
	{
		AssistUij[T_ID][i][0] = 0.0;
		AssistUij[T_ID][i][1] = 0.0;
		AssistUij[T_ID][i][2] = 0.0;
		AssistUij[T_ID][i][4] = 0.0;
		AssistUij[T_ID][i][5] = 0.0;
		AssistUij[T_ID][i][8] = 0.0;
	}

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(1, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				GetFromLocal(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1], inputPW);
				GetR(source, inputPW, R, RI);
				
				
				//GetR(source, m_intptPW[i][j], R, RI);	//计算高斯积分点绝对坐标
				GetDynaUij(UT, n, dt, R, RI);

				JacobiPW = Jacobi(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1]);
				for (PointID = 0; PointID < 8; ++PointID)
				{
					temp = m_quadinfo.m_NPWGV[i][PointID][j] * JacobiPW;

					AssistUij[T_ID][PointID][0] += UT[0] * temp;
					AssistUij[T_ID][PointID][1] += UT[1] * temp;
					AssistUij[T_ID][PointID][2] += UT[2] * temp;
					AssistUij[T_ID][PointID][4] += UT[4] * temp;
					AssistUij[T_ID][PointID][5] += UT[5] * temp;
					AssistUij[T_ID][PointID][8] += UT[8] * temp;
				}
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

int DSquareElement::IntDynaUijPW(Point& source, double n, double dt, int ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp = 0;
	int i, j;

	m_AssistUij[ID][0] = 0.0;
	m_AssistUij[ID][1] = 0.0;
	m_AssistUij[ID][2] = 0.0;
	m_AssistUij[ID][4] = 0.0;
	m_AssistUij[ID][5] = 0.0;
	m_AssistUij[ID][8] = 0.0;

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(1, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				//GetR(source, m_intptPW[i][j], R, RI);
				GetDynaUij(UT, n, dt, R, RI);

				//temp = m_quadinfo.m_NPWGV[i][ID][j] * m_JacobiPW[i][j];
				m_AssistUij[ID][0] += UT[0] * temp;
				m_AssistUij[ID][1] += UT[1] * temp;
				m_AssistUij[ID][2] += UT[2] * temp;
				m_AssistUij[ID][4] += UT[4] * temp;
				m_AssistUij[ID][5] += UT[5] * temp;
				m_AssistUij[ID][8] += UT[8] * temp;
			}
		}
	}

	m_AssistUij[ID][3] = m_AssistUij[ID][1];
	m_AssistUij[ID][6] = m_AssistUij[ID][2];
	m_AssistUij[ID][7] = m_AssistUij[ID][5];

	return 1;
}

int DSquareElement::IntDynaTijPW(int typeT, Point& source, double n, double dt)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp = 0;
	int i, j, PointID;

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

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(typeT, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				//GetR(source, m_intptPW[i][j], R, RI);
				//GetDynaTij(typeT, UT, n, dt, R, RI, m_normalPW[i][j]);

				for (PointID = 0; PointID < 8; ++PointID)
				{
					//temp = m_quadinfo.m_NPWGV[i][PointID][j] * m_JacobiPW[i][j];
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
	}

	return 1;
}

int DSquareElement::IntDynaTijPW(int typeT, Point& source, double n, double dt, long T_ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp;
	int i, j, PointID;
	Point inputPW;
	double JacobiPW;
	double normalPW[3];

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

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(typeT, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				GetFromLocal(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1], inputPW);
				GetR(source, inputPW, R, RI);
				Normal(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1], normalPW[0], normalPW[1], normalPW[2]);
				GetDynaTij(typeT, UT, n, dt, R, RI, normalPW);

				//GetR(source, m_intptPW[i][j], R, RI);
				//GetDynaTij(typeT, UT, n, dt, R, RI, m_normalPW[i][j]);  //计算法向量

				JacobiPW = Jacobi(m_quadinfo.m_LPOSPW[i][j][0], m_quadinfo.m_LPOSPW[i][j][1]);
				for (PointID = 0; PointID < 8; ++PointID)
				{
					temp = m_quadinfo.m_NPWGV[i][PointID][j] * JacobiPW;
					AssistTij[T_ID][PointID][0] += UT[0] * temp;
					AssistTij[T_ID][PointID][1] += UT[1] * temp;
					AssistTij[T_ID][PointID][2] += UT[2] * temp;
					AssistTij[T_ID][PointID][3] += UT[3] * temp;
					AssistTij[T_ID][PointID][4] += UT[4] * temp;
					AssistTij[T_ID][PointID][5] += UT[5] * temp;
					AssistTij[T_ID][PointID][6] += UT[6] * temp;
					AssistTij[T_ID][PointID][7] += UT[7] * temp;
					AssistTij[T_ID][PointID][8] += UT[8] * temp;
				}
			}
		}
	}

	return 1;
}

int DSquareElement::IntDynaTijPW(int typeT, Point& source, double n, double dt, int ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT[9], temp = 0;
	int i, j;

	m_AssistTij[ID][0] = 0.0;
	m_AssistTij[ID][1] = 0.0;
	m_AssistTij[ID][2] = 0.0;
	m_AssistTij[ID][3] = 0.0;
	m_AssistTij[ID][4] = 0.0;
	m_AssistTij[ID][5] = 0.0;
	m_AssistTij[ID][6] = 0.0;
	m_AssistTij[ID][7] = 0.0;
	m_AssistTij[ID][8] = 0.0;

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(typeT, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				//GetR(source, m_intptPW[i][j], R, RI);
				//GetDynaTij(typeT, UT, n, dt, R, RI, m_normalPW[i][j]);

				//temp = m_quadinfo.m_NPWGV[i][ID][j] * m_JacobiPW[i][j];
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
	}

	return 1;
}

int DSquareElement::IntDynaTijPW(Point& source, double n, double dt)
{
	// n = 0..M-1

	double R[6], RI[4], UT1[9], UT2[9], temp=0;
	int i, j, PointID;

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

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(1, n, m_rcenPW[i], source) || PWNeedCal(2, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				//GetR(source, m_intptPW[i][j], R, RI);
				//GetDynaTij(1, UT1, n, dt, R, RI, m_normalPW[i][j]);
				//GetDynaTij(2, UT2, n, dt, R, RI, m_normalPW[i][j]);

				for (PointID = 0; PointID < 8; ++PointID)
				{
					//temp = m_quadinfo.m_NPWGV[i][PointID][j] * m_JacobiPW[i][j];
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
	}

	return 1;
}

int DSquareElement::IntDynaTijPW(Point& source, double n, double dt, int ID)
{
	// n = 0..M-1

	double R[6], RI[4], UT1[9], UT2[9], temp=0;
	int i, j;

	m_AssistTij[ID][0] = 0.0;
	m_AssistTij[ID][1] = 0.0;
	m_AssistTij[ID][2] = 0.0;
	m_AssistTij[ID][3] = 0.0;
	m_AssistTij[ID][4] = 0.0;
	m_AssistTij[ID][5] = 0.0;
	m_AssistTij[ID][6] = 0.0;
	m_AssistTij[ID][7] = 0.0;
	m_AssistTij[ID][8] = 0.0;

	for (i = 0; i < NUMPW2; ++i)
	{
		if (PWNeedCal(1, n, m_rcenPW[i], source) || PWNeedCal(2, n, m_rcenPW[i], source))
		{
			for (j = 0; j < GAUSSPOINTPW2; ++j)
			{
				//GetR(source, m_intptPW[i][j], R, RI);
				//GetDynaTij(1, UT1, n, dt, R, RI, m_normalPW[i][j]);
				//GetDynaTij(2, UT2, n, dt, R, RI, m_normalPW[i][j]);

				//temp = m_quadinfo.m_NPWGV[i][ID][j] * m_JacobiPW[i][j];
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
	}

	return 1;
}
//无限单元,没用这个
int DSquareElement::IntDynaTijPWforInfEle(Point& source, double n, double dt, int type)
{
	// n = 0..M-1
	double  RI[4], UT1[9], UT2[9], temp, arpha, beta, gamma;
	double ru[3], rv[3], drda[3], drdb[3], E, F, G, Jacobi_Inf;
	double res[8], m_normalforinf[3];
	int i, ii,j, k, PointID;
	Point R_bd, R_gp;

	double R[6];
	double m_CENPW[2], m_LPOSPW[2];

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

	double m_LocLength = 2.0 / NUMPW;//考虑施加一个给定分片数专门针对无限单元，暂未实现


	for (i = 0; i < NUMPW; ++i)
	{
		for (ii= 0; ii < NUMPW; ++ii)
		{
			m_CENPW[0] = -1.0 + m_LocLength / 2.0 + m_LocLength * ii;
			m_CENPW[1] = -1.0 + m_LocLength / 2.0 + m_LocLength * i;
			for (k = 0; k < GAUSSPOINTPW2; k++)//可以优化增加块的判断，还未实现
			{
				m_LPOSPW[0] = m_CENPW[0] + m_gnm.pwrp[k][0] / NUMPW;
				m_LPOSPW[1] = m_CENPW[1] + m_gnm.pwrp[k][1] / NUMPW;
				if (type == 1)
				{
					//获得高斯点
					arpha = m_LPOSPW[0];
					gamma = m_LPOSPW[1];
					beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
					//计算R
					GetFromLocal(-arpha, -1.0, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(-arpha, -1.0, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(1, UT1, n, dt, R, RI, m_normalforinf);
					//GetDynaTij(2, UT2, n, dt, R, RI, m_normalforinf);
					//计算J
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
					//计算形函数
					DSquareElement::m_quadinfo.N(res, -arpha, -1.0);//局部坐标
					//计算积分
					for (PointID = 0; PointID < 8; ++PointID)
					{
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0]) * temp;
						m_AssistTij[PointID][1] += (UT1[1] ) * temp;
						m_AssistTij[PointID][2] += (UT1[2] ) * temp;
						m_AssistTij[PointID][3] += (UT1[3] ) * temp;
						m_AssistTij[PointID][4] += (UT1[4] ) * temp;
						m_AssistTij[PointID][5] += (UT1[5] ) * temp;
						m_AssistTij[PointID][6] += (UT1[6]) * temp;
						m_AssistTij[PointID][7] += (UT1[7]) * temp;
						m_AssistTij[PointID][8] += (UT1[8] ) * temp;
					}
				}
				else if (type == 2)
				{
					//获得高斯点
					arpha = m_LPOSPW[0];
					gamma = m_LPOSPW[1];
					beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
					//计算R
					GetFromLocal(1.0, -arpha, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(1.0, -arpha, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(1, UT1, n, dt, R, RI, m_normalforinf);
					GetDynaTij(2, UT2, n, dt, R, RI, m_normalforinf);
					//计算J
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
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0] ) * temp;
						m_AssistTij[PointID][1] += (UT1[1] ) * temp;
						m_AssistTij[PointID][2] += (UT1[2]) * temp;
						m_AssistTij[PointID][3] += (UT1[3] ) * temp;
						m_AssistTij[PointID][4] += (UT1[4] ) * temp;
						m_AssistTij[PointID][5] += (UT1[5] ) * temp;
						m_AssistTij[PointID][6] += (UT1[6] ) * temp;
						m_AssistTij[PointID][7] += (UT1[7] ) * temp;
						m_AssistTij[PointID][8] += (UT1[8] ) * temp;
					}
				}
				else if (type == 3)
				{
					//获得高斯点
					arpha = m_LPOSPW[0];
					gamma = m_LPOSPW[1];
					beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
					//计算R
					GetFromLocal(arpha, 1.0, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(arpha, 1.0, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(1, UT1, n, dt, R, RI, m_normalforinf);
					//GetDynaTij(2, UT2, n, dt, R, RI, m_normalforinf);
					//计算J
					RDiffs1(arpha, 1.0, ru);//局部坐标
					RDiffs2(arpha, 1.0, rv);//局部坐标
					for (j = 0; j < 3; j++) {
						drda[j] = (3.0 + beta) / 2.0*ru[j];//偏导对象ru rv (-1)
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
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0] ) * temp;
						m_AssistTij[PointID][1] += (UT1[1] ) * temp;
						m_AssistTij[PointID][2] += (UT1[2] ) * temp;
						m_AssistTij[PointID][3] += (UT1[3] ) * temp;
						m_AssistTij[PointID][4] += (UT1[4] ) * temp;
						m_AssistTij[PointID][5] += (UT1[5] ) * temp;
						m_AssistTij[PointID][6] += (UT1[6] ) * temp;
						m_AssistTij[PointID][7] += (UT1[7] ) * temp;
						m_AssistTij[PointID][8] += (UT1[8] ) * temp;
					}
				}
				else if (type == 4)
				{
					//获得高斯点
					arpha = m_LPOSPW[0];
					gamma = m_LPOSPW[1];
					beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
					//计算R
					GetFromLocal(-1.0, arpha, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(-1.0, arpha, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(1, UT1, n, dt, R, RI, m_normalforinf);
					//GetDynaTij(2, UT2, n, dt, R, RI, m_normalforinf);
					//计算J
					RDiffs1(-1.0, arpha, ru);//局部坐标
					RDiffs2(-1.0, arpha, rv);//局部坐标
					for (j = 0; j < 3; j++) {
						drda[j] = (3.0 + beta) / 2.0*rv[j];//偏导对象ru rv (-1)
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
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0]) * temp;
						m_AssistTij[PointID][1] += (UT1[1]) * temp;
						m_AssistTij[PointID][2] += (UT1[2] ) * temp;
						m_AssistTij[PointID][3] += (UT1[3] ) * temp;
						m_AssistTij[PointID][4] += (UT1[4]) * temp;
						m_AssistTij[PointID][5] += (UT1[5] ) * temp;
						m_AssistTij[PointID][6] += (UT1[6]) * temp;
						m_AssistTij[PointID][7] += (UT1[7] ) * temp;
						m_AssistTij[PointID][8] += (UT1[8] ) * temp;
					}
				}
			}
		}
	}


	return 1;
}
//无限单元用
int DSquareElement::IntDynaTijPWforInfEle(int Ttype, Point& source, double n, double dt, int type)
{
	// n = 0..M-1
	double  RI[4], UT1[9], temp, arpha, beta, gamma;
	double ru[3], rv[3], drda[3], drdb[3], E, F, G, Jacobi_Inf;
	double res[8], m_normalforinf[3];
	int i,ii, j, k, PointID;
	Point R_bd, R_gp;

	double R[6];
	double m_CENPW[2], m_LPOSPW[2];

	for (i = 0; i < 9; i++) {
		UT1[i] = 0.0;

	}

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

	double m_LocLength = 2.0 / NUMPW;//考虑施加一个给定分片数专门针对无限单元，暂未实现

	for (i = 0; i < NUMPW; ++i)
	{
		for (ii = 0; ii < NUMPW; ++ii)
		{
			m_CENPW[0] = -1.0 + m_LocLength / 2.0 + m_LocLength * ii;
			m_CENPW[1] = -1.0 + m_LocLength / 2.0 + m_LocLength * i;
			for (k = 0; k < GAUSSPOINTPW2; k++)//可以优化增加块的判断，还未实现
			{
				m_LPOSPW[0] = m_CENPW[0] + m_gnm.pwrp[k][0] / NUMPW;
				m_LPOSPW[1] = m_CENPW[1] + m_gnm.pwrp[k][1] / NUMPW;
				//获得高斯点
				arpha = m_LPOSPW[0];
				gamma = m_LPOSPW[1];
				beta = (1.0 + 3.0*gamma) / (1.0 - gamma);
				if (type == 1)
				{
					//计算R
					GetFromLocal(-arpha, -1.0, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(-arpha, -1.0, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(Ttype, UT1, n, dt, R, RI, m_normalforinf);
					//计算J
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
					//计算形函数
					DSquareElement::m_quadinfo.N(res, -arpha, -1.0);//局部坐标
					//计算积分
					for (PointID = 0; PointID < 8; ++PointID)
					{
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0]) * temp;
						m_AssistTij[PointID][1] += (UT1[1]) * temp;
						m_AssistTij[PointID][2] += (UT1[2]) * temp;
						m_AssistTij[PointID][3] += (UT1[3]) * temp;
						m_AssistTij[PointID][4] += (UT1[4]) * temp;
						m_AssistTij[PointID][5] += (UT1[5]) * temp;
						m_AssistTij[PointID][6] += (UT1[6]) * temp;
						m_AssistTij[PointID][7] += (UT1[7]) * temp;
						m_AssistTij[PointID][8] += (UT1[8]) * temp;
					}
				}
				else if (type == 2)
				{
					//计算R
					GetFromLocal(1.0, -arpha, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(1.0, -arpha, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(Ttype, UT1, n, dt, R, RI, m_normalforinf);
					//计算J
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
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0]) * temp;
						m_AssistTij[PointID][1] += (UT1[1]) * temp;
						m_AssistTij[PointID][2] += (UT1[2]) * temp;
						m_AssistTij[PointID][3] += (UT1[3]) * temp;
						m_AssistTij[PointID][4] += (UT1[4]) * temp;
						m_AssistTij[PointID][5] += (UT1[5]) * temp;
						m_AssistTij[PointID][6] += (UT1[6]) * temp;
						m_AssistTij[PointID][7] += (UT1[7]) * temp;
						m_AssistTij[PointID][8] += (UT1[8]) * temp;
					}
				}
				else if (type == 3)
				{
					//计算R
					GetFromLocal(arpha, 1.0, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(arpha, 1.0, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(Ttype, UT1, n, dt, R, RI, m_normalforinf);
					//计算J
					RDiffs1(arpha, 1.0, ru);//局部坐标
					RDiffs2(arpha, 1.0, rv);//局部坐标
					for (j = 0; j < 3; j++) {
						drda[j] = (3.0 + beta) / 2.0*ru[j];//偏导对象ru rv (-1)
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
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0]) * temp;
						m_AssistTij[PointID][1] += (UT1[1]) * temp;
						m_AssistTij[PointID][2] += (UT1[2]) * temp;
						m_AssistTij[PointID][3] += (UT1[3]) * temp;
						m_AssistTij[PointID][4] += (UT1[4]) * temp;
						m_AssistTij[PointID][5] += (UT1[5]) * temp;
						m_AssistTij[PointID][6] += (UT1[6]) * temp;
						m_AssistTij[PointID][7] += (UT1[7]) * temp;
						m_AssistTij[PointID][8] += (UT1[8]) * temp;
					}
				}
				else if (type == 4)
				{
					//计算R
					GetFromLocal(-1.0, arpha, R_bd);//局部坐标
					R_gp.pt[0] = (3.0 + beta) / 2.0*(R_bd.pt[0] - DSquareElement::InfCenter[0]) + DSquareElement::InfCenter[0];
					R_gp.pt[1] = (3.0 + beta) / 2.0*(R_bd.pt[1] - DSquareElement::InfCenter[1]) + DSquareElement::InfCenter[1];
					R_gp.pt[2] = (3.0 + beta) / 2.0*(R_bd.pt[2] - DSquareElement::InfCenter[2]) + DSquareElement::InfCenter[2];
					GetR(source, R_gp, R, RI);
					//计算n
					Normal(-1.0, arpha, m_normalforinf[0], m_normalforinf[1], m_normalforinf[2]);//局部坐标
					//计算T
					GetDynaTij(Ttype, UT1, n, dt, R, RI, m_normalforinf);
					//计算J
					RDiffs1(-1.0, arpha, ru);//局部坐标
					RDiffs2(-1.0, arpha, rv);//局部坐标
					for (j = 0; j < 3; j++) {
						drda[j] = (3.0 + beta) / 2.0*rv[j];//偏导对象ru rv (-1)
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
						temp = 2.0 / (1.0 - gamma)*res[PointID] * Jacobi_Inf*m_gnm.ra[i] * (m_LocLength / 2.0) * (m_LocLength / 2.0);
						m_AssistTij[PointID][0] += (UT1[0]) * temp;
						m_AssistTij[PointID][1] += (UT1[1]) * temp;
						m_AssistTij[PointID][2] += (UT1[2]) * temp;
						m_AssistTij[PointID][3] += (UT1[3]) * temp;
						m_AssistTij[PointID][4] += (UT1[4]) * temp;
						m_AssistTij[PointID][5] += (UT1[5]) * temp;
						m_AssistTij[PointID][6] += (UT1[6]) * temp;
						m_AssistTij[PointID][7] += (UT1[7]) * temp;
						m_AssistTij[PointID][8] += (UT1[8]) * temp;
					}
				}
			}
		}
	}
	return 1;
}
				
				

int DSquareElement::IntDynaUijJudge(Point& source, double n, double dt)
{
	int Flag = NeedCal(1, n, source);

	if (Flag == 0)
		return 0;

	if (Flag == 1)
		IntDynaUijPW(source, n, ActiveDynaMat().Dt);
	else if (Flag == 2)
		IntDynaUij(source, n, ActiveDynaMat().Dt);

	return Flag;
}

int DSquareElement::IntDynaUijJudge(Point& source, double n, double dt, int ID)
{
	int Flag = NeedCal(1, n, source);

	if (Flag == 0)
		return 0;

	if (Flag == 1)
		IntDynaUijPW(source, n, ActiveDynaMat().Dt, ID);
	else if (Flag == 2)
		IntDynaUij(source, n, ActiveDynaMat().Dt, ID);

	return Flag;
}

int DSquareElement::IntDynaUijJudge(Point& source, double n, double dt, long T_ID)
{
	int Flag = NeedCal(1, n, source);

	if (Flag == 0)
		return 0;

	if (Flag == 1)
		IntDynaUijPW(source, n, ActiveDynaMat().Dt, T_ID);
	else if (Flag == 2)
		IntDynaUij(source, n, ActiveDynaMat().Dt, T_ID);

	return Flag;
}

int DSquareElement::IntDynaJudge(Point& source, double n, double dt, long T_ID)
{
	int Flag = NeedCal(1, n, source);

	if (Flag == 0) {
		return 0;
	}
	else {
		return 1;
	}

}

int DSquareElement::IntDynaTijJudge(int typeT, Point& source, double n, double dt)
{
	int Flag = NeedCal(typeT, n, source);

	if (Flag == 0)
		return 0;

	if (Flag == 1)
		IntDynaTijPW(typeT, source, n, ActiveDynaMat().Dt);
	else if (Flag == 2)
		IntDynaTij(typeT, source, n, ActiveDynaMat().Dt);

	return Flag;
}

int DSquareElement::IntDynaTijJudge(int typeT, Point& source, double n, double dt, int ID)
{
	int Flag = NeedCal(typeT, n, source);

	if (Flag == 0)
		return 0;

	if (Flag == 1)
		IntDynaTijPW(typeT, source, n, ActiveDynaMat().Dt, ID);
	else if (Flag == 2)
		IntDynaTij(typeT, source, n, ActiveDynaMat().Dt, ID);

	return Flag;
}

int DSquareElement::IntDynaTijJudge(Point& source, double n, double dt)
{
	int Flag1 = NeedCal(1, n, source);
	int Flag2 = NeedCal(2, n, source);

	if (Flag1 == 0 && Flag2 == 0)
		return 0;

	if (Flag1 == 1 || Flag2 == 1)
		IntDynaTijPW(source, n, ActiveDynaMat().Dt);
	else
		IntDynaTij(source, n, ActiveDynaMat().Dt);

	return (Flag1 || Flag2);
}

int DSquareElement::IntDynaTijJudge(int typeT, Point& source, double n, double dt, long T_ID)
{
	int Flag = NeedCal(typeT, n, source);

	if (Flag == 0)
		return 0;

	if (Flag == 1)
		IntDynaTijPW(typeT, source, n, ActiveDynaMat().Dt, T_ID);
	else if (Flag == 2)
		IntDynaTij(typeT, source, n, ActiveDynaMat().Dt, T_ID);

	return Flag;
}

int DSquareElement::IntDynaTijJudge(Point& source, double n, double dt, int ID)
{
	int Flag1 = NeedCal(1, n, source);
	int Flag2 = NeedCal(2, n, source);

	if (Flag1 == 0 && Flag2 == 0)
		return 0;

	if (Flag1 == 1 || Flag2 == 1)
		IntDynaTijPW(source, n, ActiveDynaMat().Dt, ID);
	else
		IntDynaTij(source, n, ActiveDynaMat().Dt, ID);

	return (Flag1 || Flag2);
}

int DSquareElement::IntDynaTijJudge(Point& source, double n, double dt, long T_ID)
{
	return 0;
}
