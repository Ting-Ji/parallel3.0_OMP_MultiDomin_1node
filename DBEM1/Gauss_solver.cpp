#include "Gauss_solver.h"
#include "lin_comp.h"
#include "stdio.h"
#include "Counter.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include "output_path.h"
using namespace std;

int GetMatrixDynaUij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt)
{

	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaUijJudge(m_DSE[j].m_nodelist[i], n, dt);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistUij[l];
					A.assign(DSquareElement::m_AssistUij[l], pos + 1, 3 * k + 1, 3, 3);
				}
			}
		}
	}

	return 1;
}

int GetMatrixDynaUij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, double n, double dt)
{
	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	//A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaUijJudge(m_DSE[j].m_nodelist[i], n, dt);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistUij[l];
					A.assign(DSquareElement::m_AssistUij[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}

	return 1;
}

int GetMatrixDynaUij(DSquareElement* m_DSE, WmatrixSymCCSR& A, long NodeNum, long EleNum, double n, double dt)
{
	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	//A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaUijJudge(m_DSE[j].m_nodelist[i], n, dt);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistUij[l];
					A.assign(DSquareElement::m_AssistUij[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}

	return 1;
}

int GetMatrixDynaT1ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt)
{

	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaTijJudge(1, m_DSE[j].m_nodelist[i], n, dt);

			if (Flag)
			{
				for (l = 0; l < 8; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistTij[l];
					A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
				}
			}
		}
	}

	return 1;
}
//无限单元
int GetMatrixDynaT1ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID)
{

	long i, j, k, l, pos;//
	Wmatrix sub(3, 3);
	long m_InfEleFlag;
	int Flag;
	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaTijJudge(1, m_DSE[j].m_nodelist[i], n, dt);
			//
			//
			if (Flag)
			{

				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistTij[l];
					A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
				}
			}

			//无限单元部分
			m_InfEleFlag = m_DSE[j].InfEleFlag;
			if (m_InfEleFlag != (-1))
			{
				if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
				}
				else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
				}
				else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
				}
				else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
				}
				//集成



				//for (l = 0; l < 8; ++l)
				//{
				//	k = m_DSE[j].m_nodeID[l];
				//	// sub = DSquareElement::m_AssistTij[l];
				//	A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
				//}
				//for (l = 0; l < 8; ++l)
				//{
				//	k = m_DSE[j].m_nodeID[l];
				//	// sub = DSquareElement::m_AssistTij[l];
				//	A.assign(temp[l], pos + 1, 3 * k + 1, 3, 3);
				//}
			}
			//for (l = 0; l < 8; ++l)
			//{
			//	for (jj = 0; jj < 9; ++jj) {
			//		temp[l][jj] = tempa[l][jj] + tempb[l][jj];
			//	}
			//}
			for (l = 0; l < 1; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				// sub = DSquareElement::m_AssistTij[l];
				A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
			}


				
		}
	}

	return 1;
}

int GetMatrixDynaT1ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos, ii, jj;
	Wmatrix sub(3, 3);
	long m_InfEleFlag;
	int Flag;
	double temp[8][9];

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;


	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			for (ii = 0; ii < 1; ++ii) {
				for (jj = 0; jj < 9; ++jj) {
					temp[ii][jj] = 0.0;
				}
			}

			Flag = m_DSE[j].IntDynaTijJudge(1, m_DSE[j].m_nodelist[i], n, dt);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistTij[l];
					//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);

					for (jj = 0; jj < 9; ++jj) {
						temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
					}
				}
			}

			//无限单元部分
			m_InfEleFlag = m_DSE[j].InfEleFlag;
			if (m_InfEleFlag != (-1))
			{
				if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
				}
				else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
				}
				else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
				}
				else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
					m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
				}
				//集成
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistTij[l];
					for (jj = 0; jj < 9; ++jj) {
						temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						
					}
					//system("pause");
					//A.assign(temp[l], pos, 3 * k, 3, 3);
					//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					
				}
			}
			
			for (l = 0; l < 1; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				A.assign(temp[l], pos, 3 * k, 3, 3);
			}
		}
	}

	return 1;
}

int GetMatrixDynaT1ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos, ii, jj;
	Wmatrix sub(3, 3);
	long m_InfEleFlag;
	int Flag;
	double temp[8][9];

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;


	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			for (ii = 0; ii < 1; ++ii) {
				for (jj = 0; jj < 9; ++jj) {
					temp[ii][jj] = 0.0;
				}
			}

			Flag = m_DSE[j].IntDynaTijJudge(1, m_DSE[j].m_nodelist[i], n, dt);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistTij[l];
					//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);

					for (jj = 0; jj < 9; ++jj) {
						temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
					}
				}
			}

			////无限单元部分
			//m_InfEleFlag = m_DSE[j].InfEleFlag;
			//if (m_InfEleFlag != (-1))
			//{
			//	if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
			//	}
			//	//集成
			//	for (l = 0; l < 1; ++l)
			//	{
			//		k = m_DSE[j].m_nodeID[l];
			//		// sub = DSquareElement::m_AssistTij[l];
			//		for (jj = 0; jj < 9; ++jj) {
			//			temp[l][jj] += DSquareElement::m_AssistTij[l][jj];

			//		}
			//		//system("pause");
			//		//A.assign(temp[l], pos, 3 * k, 3, 3);
			//		//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);

			//	}
			//}

			for (l = 0; l < 1; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				A.assign(temp[l], pos, 3 * k, 3, 3);
			}
		}
	}

	return 1;
}

int GetMatrixDynaT2ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt)
{

	long i, j, k, l, pos;
	int ipos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	if (n > 1.01) // 常规积分
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 8; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}
			}
		}
	}
	else // n = 1
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				if (m_DSE[j].IsIn(i, ipos))
				{
					for (l = 0; l < 8; ++l)
					{
						DSquareElement::m_AssistTij[l][0] = m_DSE[j].m_OnElementT2ij[ipos][l][0];
						DSquareElement::m_AssistTij[l][1] = m_DSE[j].m_OnElementT2ij[ipos][l][1];
						DSquareElement::m_AssistTij[l][2] = m_DSE[j].m_OnElementT2ij[ipos][l][2];
						DSquareElement::m_AssistTij[l][3] = m_DSE[j].m_OnElementT2ij[ipos][l][3];
						DSquareElement::m_AssistTij[l][4] = m_DSE[j].m_OnElementT2ij[ipos][l][4];
						DSquareElement::m_AssistTij[l][5] = m_DSE[j].m_OnElementT2ij[ipos][l][5];
						DSquareElement::m_AssistTij[l][6] = m_DSE[j].m_OnElementT2ij[ipos][l][6];
						DSquareElement::m_AssistTij[l][7] = m_DSE[j].m_OnElementT2ij[ipos][l][7];
						DSquareElement::m_AssistTij[l][8] = m_DSE[j].m_OnElementT2ij[ipos][l][8];
					}

					Flag = 1;
				}
				else
					Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 8; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}
			}
		}
	}


	return 1;
}
//无限单元
int GetMatrixDynaT2ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID)
{

	long i, j, k, l, pos;
	int ipos;
	Wmatrix sub(3, 3);

	int Flag;
	long m_InfEleFlag;
	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	if (n > 1.01) // 常规积分
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}
				//无限单元部分
				m_InfEleFlag = m_DSE[j].InfEleFlag;
				if (m_InfEleFlag != (-1))
				{
					if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
					}
					//集成
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}

			}
		}
	}
	else // n = 1
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				if (m_DSE[j].IsIn(i, ipos))
				{
					for (l = 0; l < 1; ++l)
					{
						DSquareElement::m_AssistTij[l][0] = m_DSE[j].m_OnElementT2ij[ipos][l][0];
						DSquareElement::m_AssistTij[l][1] = m_DSE[j].m_OnElementT2ij[ipos][l][1];
						DSquareElement::m_AssistTij[l][2] = m_DSE[j].m_OnElementT2ij[ipos][l][2];
						DSquareElement::m_AssistTij[l][3] = m_DSE[j].m_OnElementT2ij[ipos][l][3];
						DSquareElement::m_AssistTij[l][4] = m_DSE[j].m_OnElementT2ij[ipos][l][4];
						DSquareElement::m_AssistTij[l][5] = m_DSE[j].m_OnElementT2ij[ipos][l][5];
						DSquareElement::m_AssistTij[l][6] = m_DSE[j].m_OnElementT2ij[ipos][l][6];
						DSquareElement::m_AssistTij[l][7] = m_DSE[j].m_OnElementT2ij[ipos][l][7];
						DSquareElement::m_AssistTij[l][8] = m_DSE[j].m_OnElementT2ij[ipos][l][8];
					}

					Flag = 1;
				}
				else
					Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}
				//无限单元部分
				m_InfEleFlag = m_DSE[j].InfEleFlag;
				if (m_InfEleFlag != (-1))
				{
					if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
					}
					//集成
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}



			}
		}
	}


	return 1;
}

int GetMatrixDynaT2ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos, ii, jj;
	int ipos;
	Wmatrix sub(3, 3);
	double temp[8][9];
	int Flag;
	long m_InfEleFlag;
	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	if (n > 1.01) // 常规积分
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 1; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}

				Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						// sub = DSquareElement::m_AssistTij[l];
						//A.assign(DSquareElement::m_AssistTij[l], pos, 3 * k, 3, 3);
					}
				}
				//无限单元部分
				m_InfEleFlag = m_DSE[j].InfEleFlag;
				if (m_InfEleFlag != (-1))
				{
					if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
					}
					//集成
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						//A.assign(temp[l], pos, 3 * k, 3, 3);
					}
					//system("pause");
				}
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					A.assign(temp[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}
	else // n = 1
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 1; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}
				if (m_DSE[j].IsIn(i, ipos))
				{
					for (l = 0; l < 1; ++l)
					{
						DSquareElement::m_AssistTij[l][0] = m_DSE[j].m_OnElementT2ij[ipos][l][0];
						DSquareElement::m_AssistTij[l][1] = m_DSE[j].m_OnElementT2ij[ipos][l][1];
						DSquareElement::m_AssistTij[l][2] = m_DSE[j].m_OnElementT2ij[ipos][l][2];
						DSquareElement::m_AssistTij[l][3] = m_DSE[j].m_OnElementT2ij[ipos][l][3];
						DSquareElement::m_AssistTij[l][4] = m_DSE[j].m_OnElementT2ij[ipos][l][4];
						DSquareElement::m_AssistTij[l][5] = m_DSE[j].m_OnElementT2ij[ipos][l][5];
						DSquareElement::m_AssistTij[l][6] = m_DSE[j].m_OnElementT2ij[ipos][l][6];
						DSquareElement::m_AssistTij[l][7] = m_DSE[j].m_OnElementT2ij[ipos][l][7];
						DSquareElement::m_AssistTij[l][8] = m_DSE[j].m_OnElementT2ij[ipos][l][8];
					}

					Flag = 1;
				}
				else
					Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						// sub = DSquareElement::m_AssistTij[l];
						//A.assign(DSquareElement::m_AssistTij[l], pos, 3 * k, 3, 3);
					}
				}
				//无限单元部分
				m_InfEleFlag = m_DSE[j].InfEleFlag;
				if (m_InfEleFlag != (-1))
				{
					if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
					}
					//集成
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						//A.assign(DSquareElement::m_AssistTij[l], pos , 3 * k , 3, 3);
					}
				}

				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					A.assign(temp[l], pos, 3 * k, 3, 3);
				}


			}
		}
	}


	return 1;
}

int GetMatrixDynaT2ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, double n, double dt, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos, ii, jj;
	int ipos;
	Wmatrix sub(3, 3);
	double temp[8][9];
	int Flag;
	long m_InfEleFlag;
	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	if (n > 1.01) // 常规积分
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 1; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}

				Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						// sub = DSquareElement::m_AssistTij[l];
						//A.assign(DSquareElement::m_AssistTij[l], pos, 3 * k, 3, 3);
					}
				}
				//无限单元部分
				m_InfEleFlag = m_DSE[j].InfEleFlag;
				if (m_InfEleFlag != (-1))
				{
					if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
					}
					//集成
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						//A.assign(temp[l], pos, 3 * k, 3, 3);
					}
					//system("pause");
				}
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					A.assign(temp[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}
	else // n = 1
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 1; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}
				if (m_DSE[j].IsIn(i, ipos))
				{
					for (l = 0; l < 1; ++l)
					{
						DSquareElement::m_AssistTij[l][0] = m_DSE[j].m_OnElementT2ij[ipos][l][0];
						DSquareElement::m_AssistTij[l][1] = m_DSE[j].m_OnElementT2ij[ipos][l][1];
						DSquareElement::m_AssistTij[l][2] = m_DSE[j].m_OnElementT2ij[ipos][l][2];
						DSquareElement::m_AssistTij[l][3] = m_DSE[j].m_OnElementT2ij[ipos][l][3];
						DSquareElement::m_AssistTij[l][4] = m_DSE[j].m_OnElementT2ij[ipos][l][4];
						DSquareElement::m_AssistTij[l][5] = m_DSE[j].m_OnElementT2ij[ipos][l][5];
						DSquareElement::m_AssistTij[l][6] = m_DSE[j].m_OnElementT2ij[ipos][l][6];
						DSquareElement::m_AssistTij[l][7] = m_DSE[j].m_OnElementT2ij[ipos][l][7];
						DSquareElement::m_AssistTij[l][8] = m_DSE[j].m_OnElementT2ij[ipos][l][8];
					}

					Flag = 1;
				}
				else
					Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						// sub = DSquareElement::m_AssistTij[l];
						//A.assign(DSquareElement::m_AssistTij[l], pos, 3 * k, 3, 3);
					}
				}
				//无限单元部分
				m_InfEleFlag = m_DSE[j].InfEleFlag;
				if (m_InfEleFlag != (-1))
				{
					if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
					}
					else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
						m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
					}
					//集成
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
						}
						//A.assign(DSquareElement::m_AssistTij[l], pos , 3 * k , 3, 3);
					}
				}

				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					A.assign(temp[l], pos, 3 * k, 3, 3);
				}


			}
		}
	}


	return 1;
}

int GetMatrixDynaTij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, double n, double dt)
{

	long i, j, k, l, pos;
	int ipos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	if (n > 1.01) // 常规积分
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				Flag = m_DSE[j].IntDynaTijJudge(m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 8; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}
			}
		}
	}
	else // n = 1
	{
		for (i = 0; i < NodeNum; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				if (m_DSE[j].IsIn(i, ipos))
				{
					Flag = m_DSE[j].IntDynaTijJudge(1, m_DSE[j].m_nodelist[i], n, dt);

					if (Flag)
					{
						for (l = 0; l < 8; ++l)
						{
							DSquareElement::m_AssistTij[l][0] += m_DSE[j].m_OnElementT2ij[ipos][l][0];
							DSquareElement::m_AssistTij[l][1] += m_DSE[j].m_OnElementT2ij[ipos][l][1];
							DSquareElement::m_AssistTij[l][2] += m_DSE[j].m_OnElementT2ij[ipos][l][2];
							DSquareElement::m_AssistTij[l][3] += m_DSE[j].m_OnElementT2ij[ipos][l][3];
							DSquareElement::m_AssistTij[l][4] += m_DSE[j].m_OnElementT2ij[ipos][l][4];
							DSquareElement::m_AssistTij[l][5] += m_DSE[j].m_OnElementT2ij[ipos][l][5];
							DSquareElement::m_AssistTij[l][6] += m_DSE[j].m_OnElementT2ij[ipos][l][6];
							DSquareElement::m_AssistTij[l][7] += m_DSE[j].m_OnElementT2ij[ipos][l][7];
							DSquareElement::m_AssistTij[l][8] += m_DSE[j].m_OnElementT2ij[ipos][l][8];
						}
					}
					else
					{
						for (l = 0; l < 8; ++l)
						{
							DSquareElement::m_AssistTij[l][0] = m_DSE[j].m_OnElementT2ij[ipos][l][0];
							DSquareElement::m_AssistTij[l][1] = m_DSE[j].m_OnElementT2ij[ipos][l][1];
							DSquareElement::m_AssistTij[l][2] = m_DSE[j].m_OnElementT2ij[ipos][l][2];
							DSquareElement::m_AssistTij[l][3] = m_DSE[j].m_OnElementT2ij[ipos][l][3];
							DSquareElement::m_AssistTij[l][4] = m_DSE[j].m_OnElementT2ij[ipos][l][4];
							DSquareElement::m_AssistTij[l][5] = m_DSE[j].m_OnElementT2ij[ipos][l][5];
							DSquareElement::m_AssistTij[l][6] = m_DSE[j].m_OnElementT2ij[ipos][l][6];
							DSquareElement::m_AssistTij[l][7] = m_DSE[j].m_OnElementT2ij[ipos][l][7];
							DSquareElement::m_AssistTij[l][8] = m_DSE[j].m_OnElementT2ij[ipos][l][8];
						}
					}

					Flag = 1;
				}
				else
					Flag = m_DSE[j].IntDynaTijJudge(m_DSE[j].m_nodelist[i], n, dt);

				if (Flag)
				{
					for (l = 0; l < 8; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						// sub = DSquareElement::m_AssistTij[l];
						A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);
					}
				}
			}
		}
	}


	return 1;
}

int GetMatrixStaticUnknown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum)
{

	long i, j, k, l, pos;
	
	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;


	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			UnknownSubMatrix(m_DSE, i, j);
			for (l = 0; l < 8; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				// sub = DSquareElement::m_AssistCoef[l];
				A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
			}
		}
	}

	return 1;
}
//增加无限单元适用
int GetMatrixStaticUnknown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{

	long i, j, k, l, pos;



	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;


	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			UnknownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);
			for (l = 0; l < 8; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				// sub = DSquareElement::m_AssistCoef[l];
				A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
			}
		}
	}

	return 1;
}

int GetMatrixStaticKnown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum)
{

	long i, j, k, l, pos;



	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			knownSubMatrix(m_DSE, i, j);
			for (l = 0; l < 8; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				// sub = DSquareElement::m_AssistCoef[l];
				A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
			}
		}
	}

	return 1;
}
//增加无限单元适用
int GetMatrixStaticKnown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{

	long i, j, k, l, pos;



	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			knownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);
			for (l = 0; l < 8; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				// sub = DSquareElement::m_AssistCoef[l];
				A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
			}
		}
	}

	return 1;
}

int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum)
{

	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = UnknownSubMatrix(m_DSE, i, j);

			if (Flag)
			{
				for (l = 0; l < 8; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
				}
			}
		}
	}

	return 1;
}

//增加无限单元适用
int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{

	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = UnknownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
				}
			}
		}
	}

	return 1;
}


int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = UnknownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}

	return 1;
}

int GetMatrixDynaUnknown0ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = UnknownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);

			if (Flag)
			{
				for (l = 0; l < 8; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos, 3 * k, 3, 3);

				}
			}
		}
	}

	return 1;
}


int GetMatrixDynaKnown0ij(DSquareElement* m_DSE, WmatrixCCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = knownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}


	return 0;
}

int GetMatrixDynaKnown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum)
{

	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = knownSubMatrix(m_DSE, i, j);

			if (Flag)
			{
				for (l = 0; l < 8; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
				}
			}
		}
	}

	return 1;
}
//增加无限单元
int GetMatrixDynaKnown0ij(DSquareElement* m_DSE, Wmatrix& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{

	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = knownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos + 1, 3 * k + 1, 3, 3);
				}
			}
		}
	}

	return 1;
}

int GetMatrixDynaKnown0ij(DSquareElement* m_DSE, WmatrixCSR& A, long NodeNum, long EleNum, long** InfElePID, long** m_ElePID)
{
	long i, j, k, l, pos;
	Wmatrix sub(3, 3);

	int Flag;

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	A = 0.0;

	for (i = 0; i < NodeNum; ++i)//recycle for source element
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = knownSubMatrix(m_DSE, i, j, InfElePID, m_ElePID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					A.assign(DSquareElement::m_AssistCoef[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}

	return 1;
}


int StaticGaussSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, double MaxLength)
{
	Wmatrix* MT;
	Wmatrix* MG;

	Wvector B(3 * NodeNum);

	Counter m_Counter;
	clock_t time1;
	// Start Counting...
	m_Counter.StartCount();
	time1 = clock();
	// 1. 初始化矩阵系列

	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵

	MT = new Wmatrix[1];
	MG = new Wmatrix[1];
	printf("node=%d\n", 3 * NodeNum);
	MT[0].initial(3 * NodeNum, 3 * NodeNum);
	MG[0].initial(3 * NodeNum, 3 * NodeNum);
	/* 我认为没用到*/
		// 2. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	printf_s("\ntime1=%ld\n", clock() - time1);
	// 3. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	GetMatrixStaticUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	GetMatrixStaticKnown0ij(m_DSE, MG[0], NodeNum, EleNum);
	//GetMatrixStaticUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	//GetMatrixStaticKnown0ij(m_DSE, MG[0], NodeNum, EleNum);
	//test
	FILE *fptemp;
	fopen_s(&fptemp, "MT_S.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MT[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n输出未知系数矩阵元素完成!\n");
	fopen_s(&fptemp, "MG_S.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MG[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n输出已知系数矩阵元素完成!\n");
	// 把MT[0]求逆
	printf_s("\ntime2=%ld\n", clock() - time1);
	int InvFlag = inv_mat(MT[0], MT[0].n);
	printf_s("\ntime3=%ld\n", clock() - time1);

	// TEST
	printf("Gauss Solver:   Calculate B...\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	mat_vec_product(MG[0], bd[1].lu, B);//原本是bd1,直接改成bd0

	printf("Gauss Solver:   Calculate Result...\n");

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);//原本是bd1,直接改成bd0

	// TEST
	printf("Gauss Solver:   Convert...\n");

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);//原本是bd1,直接改成bd0
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);//原本是bd1,直接改成bd0



	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();
	printf_s("\ntime4=%ld\n", clock() - time1);
	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}
//增加无限单元适用
int StaticGaussSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum,long InfEleNum,long **m_InfElePID, long **m_ElePID, double MaxLength)
{
	Wmatrix* MT;
	Wmatrix* MG;

	Wvector B(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 1. 初始化矩阵系列

	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵

	MT = new Wmatrix[1];
	MG = new Wmatrix[1];

	MT[0].initial(3 * NodeNum, 3 * NodeNum);
	MG[0].initial(3 * NodeNum, 3 * NodeNum);
/* 我认为没用到*/
	// 2. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	GetMatrixStaticUnknown0ij(m_DSE, MT[0], NodeNum, EleNum,m_InfElePID, m_ElePID);
	GetMatrixStaticKnown0ij(m_DSE, MG[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	//GetMatrixStaticUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	//GetMatrixStaticKnown0ij(m_DSE, MG[0], NodeNum, EleNum);
	//test
	FILE *fptemp;
	fopen_s(&fptemp, "MT_S.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MT[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n输出未知系数矩阵元素完成!\n");
	fopen_s(&fptemp, "MG_S.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MG[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n输出已知系数矩阵元素完成!\n");
	// 把MT[0]求逆
	int InvFlag = inv_mat(MT[0], MT[0].n);

	// TEST
	printf("Gauss Solver:   Calculate B...\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	mat_vec_product(MG[0], bd[1].lu, B);//原本是bd1,直接改成bd0

	printf("Gauss Solver:   Calculate Result...\n");

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);//原本是bd1,直接改成bd0

	// TEST
	printf("Gauss Solver:   Convert...\n");

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);//原本是bd1,直接改成bd0
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);//原本是bd1,直接改成bd0


	
	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

int DynaGaussSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 确定最大的非零矩阵对应序号
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf("最大非零矩阵迭代步: %ld\n", MaxN);


	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}
	//	MG[NStep].initial(3 * NodeNum, 3 * NodeNum);


		// 2. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum);


	//test
	FILE *fptemp;
	fopen_s(&fptemp, "MT_D.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MT[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n!!!\n");
	fopen_s(&fptemp, "MG_D.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MG[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n!!!\n");
//	FILE* fp16;
//	fopen_s(&fp16, "UnA.txt", "w");
//	for (i = 0; i < MT[0].m; i++)
//	{
//		for (j = 0; j < MT[0].n; j++)
//		{
//			fprintf(fp16, "%20.10e ", MT[0].a[i][j]);
//		}
//		fprintf(fp16, "\n");
//	}
//	fclose(fp16);

	// 把MT[0]求逆
	int InvFlag = inv_mat(MT[0], MT[0].n);

	// TEST
	printf("Gauss Solver: Step 1. Calculate T2ij1...\n");

	GetMatrixDynaT2ij(m_DSE, MT[1], NodeNum, EleNum, 1.0, dt);

	// TEST
	printf("Gauss Solver: Step 1. Calculate B...\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	mat_vec_product(MG[0], bd[1].lu, B);

	printf("Gauss Solver: Step 1. Calculate Result...\n");

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);

	// TEST
	printf("Gauss Solver: Step 1. Convert...\n");

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("Gauss Solver: Step %ld.\n", i);

		// 使用MG[i]暂存MT_1(i-1)，并累加到原来存贮MT_2(i-1)的MT(i-1)
		if (i <= MaxN)
		{
			GetMatrixDynaUij(m_DSE, MG[i - 1], NodeNum, EleNum, (i - 1)*1.0, dt);

			GetMatrixDynaT1ij(m_DSE, MG[i], NodeNum, EleNum, (i - 1)*1.0, dt);
			MT[i - 1] += MG[i];

			// 计算MT_2(i)(如果i < NStep)，存在MT(i)中
			if (i < NStep)
				GetMatrixDynaT2ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt);
		}

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i-1
		for (j = 1; j <= i - 1; ++j)
		{
			if (j <= MaxN)
			{
				// CB1 = [MG_1(j) + MG_2(j)] * bd[i-j].lt，这里的bd[i-j].lt明确指的是i-j时刻的面力
				mat_vec_product(MG[j], bd[i - j].lt, CB1);
				B += CB1;

				// CB1 = [MT_1(j) + MT_2(j)] * bd[i-j].lu，这里的bd[i-j].lu明确指的是i-j时刻的位移
				mat_vec_product(MT[j], bd[i - j].lu, CB1);
				B -= CB1;
			}
		}

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		mat_vec_product(MT[0], B, bd[i].lt);

		// 将 n = i时刻的未知量/已知量转化为位移和面力

		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	}

	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

//增加无限单元适用
int DynaGaussSolver(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 确定最大的非零矩阵对应序号
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf("最大非零矩阵迭代步: %ld\n", MaxN);


	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}
	//	MG[NStep].initial(3 * NodeNum, 3 * NodeNum);


		// 2. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	//GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	//GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum);

	//test
	FILE *fptemp;
	fopen_s(&fptemp, "MT_D.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MT[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n输出未知系数矩阵元素完成!\n");
	fopen_s(&fptemp, "MG_D.txt", "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fptemp, " %.8lf ", MG[0](itemp + 1, jtemp + 1));
		}
		fprintf_s(fptemp, "\n");
	}
	fclose(fptemp);
	printf("\n输出已知系数矩阵元素完成!\n");
	// 把MT[0]求逆
	int InvFlag = inv_mat(MT[0], MT[0].n);

	// TEST
	printf("Gauss Solver: Step 1. Calculate T2ij1...\n");

	//GetMatrixDynaT2ij(m_DSE, MT[1], NodeNum, EleNum, 1.0, dt);
	GetMatrixDynaT2ij(m_DSE, MT[1], NodeNum, EleNum, 1.0, dt, m_InfElePID, m_ElePID);
	// TEST
	printf("Gauss Solver: Step 1. Calculate B...\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	mat_vec_product(MG[0], bd[1].lu, B);

	printf("Gauss Solver: Step 1. Calculate Result...\n");

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);

	// TEST
	printf("Gauss Solver: Step 1. Convert...\n");

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("Gauss Solver: Step %ld.\n", i);

		// 使用MG[i]暂存MT_1(i-1)，并累加到原来存贮MT_2(i-1)的MT(i-1)
		if (i <= MaxN)
		{
			GetMatrixDynaUij(m_DSE, MG[i - 1], NodeNum, EleNum, (i - 1)*1.0, dt);//无限单元不需要修改U

			GetMatrixDynaT1ij(m_DSE, MG[i], NodeNum, EleNum, (i - 1)*1.0, dt, m_InfElePID, m_ElePID);
			MT[i - 1] += MG[i];

			// 计算MT_2(i)(如果i < NStep)，存在MT(i)中
			if (i < NStep)
				GetMatrixDynaT2ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt, m_InfElePID, m_ElePID);
		}

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i-1
		for (j = 1; j <= i - 1; ++j)
		{
			if (j <= MaxN)
			{
				// CB1 = [MG_1(j) + MG_2(j)] * bd[i-j].lt，这里的bd[i-j].lt明确指的是i-j时刻的面力
				mat_vec_product(MG[j], bd[i - j].lt, CB1);
				B += CB1;

				// CB1 = [MT_1(j) + MT_2(j)] * bd[i-j].lu，这里的bd[i-j].lu明确指的是i-j时刻的位移
				mat_vec_product(MT[j], bd[i - j].lu, CB1);
				B -= CB1;
			}
		}

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		mat_vec_product(MT[0], B, bd[i].lt);

		// 将 n = i时刻的未知量/已知量转化为位移和面力

		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	}

	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}


int DynaGaussSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 确定最大的非零矩阵对应序号
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf("最大非零矩阵迭代步: %ld\n", MaxN);


	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}

	// 2. 计算矩阵

	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum);

	// 把MT[0]求逆
	int InvFlag = inv_mat(MT[0], MT[0].n);

	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		GetMatrixDynaT1ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt);
		GetMatrixDynaT2ij(m_DSE, MG[MaxN], NodeNum, EleNum, i*1.0, dt);
		MT[i] += MG[MaxN];
		GetMatrixDynaUij(m_DSE, MG[i], NodeNum, EleNum, i*1.0, dt);
	}

	// 3. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	mat_vec_product(MG[0], bd[1].lu, B);

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("Gauss Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;

				mat_vec_product(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;

				mat_vec_product(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		mat_vec_product(MT[0], B, bd[i].lt);

		// 将 n = i时刻的未知量/已知量转化为位移和面力

		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	}

	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}
//增加无限单元适用
int DynaGaussSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 确定最大的非零矩阵对应序号
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf("最大非零矩阵迭代步: %ld\n", MaxN);


	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}

	// 2. 计算矩阵
	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	
	clock_t start = clock();
	// 把MT[0]求逆
	int InvFlag = inv_mat(MT[0], MT[0].n);
	
	clock_t end = clock();
	printf("矩阵求逆完成,用时 %d ms\n",end-start);
	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		GetMatrixDynaT1ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt,m_InfElePID,m_ElePID);
		GetMatrixDynaT2ij(m_DSE, MG[MaxN], NodeNum, EleNum, i*1.0, dt, m_InfElePID, m_ElePID);
		MT[i] += MG[MaxN];
		GetMatrixDynaUij(m_DSE, MG[i], NodeNum, EleNum, i*1.0, dt);
	}

	// 3. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	mat_vec_product(MG[0], bd[1].lu, B);

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("Gauss Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;

				mat_vec_product(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;

				mat_vec_product(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		mat_vec_product(MT[0], B, bd[i].lt);

		// 将 n = i时刻的未知量/已知量转化为位移和面力

		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	}


	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

//增加无限单元适用 增加矩阵输入输出
int DynaGaussSolverNew(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long **m_InfElePID, long **m_ElePID,int m_flag)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 确定最大的非零矩阵对应序号
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf("最大非零矩阵迭代步: %ld\n", MaxN);


	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为	0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	MT = new Wmatrix[NStep + 1];
	MG = new Wmatrix[NStep + 1];

	for (i = 0; i <= MaxN; ++i)
	{
		MT[i].initial(3 * NodeNum, 3 * NodeNum);
		MG[i].initial(3 * NodeNum, 3 * NodeNum);
	}
	int InvFlag;
	clock_t start = clock();
	clock_t end;
	if (m_flag != 0)
	{
		start = clock();
		ifstream MTfile;
		MTfile.open(".\\output\\MT0", ios::in);
		
		for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
			for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
				MTfile>> MT[0](itemp + 1, jtemp + 1) ;
			}
		}
		MTfile.close();
		MTfile.clear();


		MTfile.open(".\\output\\MG0", ios::in);
		
		for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
			for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
				MTfile >> MG[0](itemp + 1, jtemp + 1) ;

			}
			
		}
		MTfile.close();
		MTfile.clear();
		cout << "读入0矩阵完成" << endl;

/*
		char nameofstep[10], openfile[30], tempchar1[20] = ".\\output\\MT", tempchar2[20] = ".\\output\\MG";

		for (i = 1; i <= MaxN; ++i)
		{
			strcpy_s(openfile, tempchar1);
			_itoa_s(i, nameofstep, 10);
			strcat_s(openfile, nameofstep);
			MTfile.open(openfile, ios::in);
			
			for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
				for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
					MTfile >> MT[i](itemp + 1, jtemp + 1) ;
				}
				
			}
			MTfile.close();
			MTfile.clear();

			strcpy_s(openfile, tempchar2);
			_itoa_s(i, nameofstep, 10);
			strcat_s(openfile, nameofstep);
			MTfile.open(openfile, ios::in);
			
			for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
				for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
					MTfile >> MG[i](itemp + 1, jtemp + 1) ;
				}
				
			}
			MTfile.close();
			MTfile.clear();

			cout<<"读入"<<i<<"矩阵完成"<<endl;
		}*/
		clock_t end = clock();
		printf("矩阵读入完成,用时 %d ms\n", end - start);

	}
	else
	{
		// 2. 计算矩阵
		GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
		GetMatrixDynaKnown0ij(m_DSE, MG[0], NodeNum, EleNum, m_InfElePID, m_ElePID);

		//GetInfoFromMatrixNstep(MG[0], MT[0], 0, NodeNum);    //输出第0步的MG、MT矩阵 


		clock_t start = clock();
		// 把MT[0]求逆
		InvFlag = inv_mat(MT[0], MT[0].n);

		clock_t end = clock();
		printf("矩阵求逆完成,用时 %d ms\n", end - start);
		
	}
	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		GetMatrixDynaT1ij(m_DSE, MT[i], NodeNum, EleNum, i*1.0, dt, m_InfElePID, m_ElePID);
		GetMatrixDynaT2ij(m_DSE, MG[MaxN], NodeNum, EleNum, i*1.0, dt, m_InfElePID, m_ElePID);
		MT[i] += MG[MaxN];
		GetMatrixDynaUij(m_DSE, MG[i], NodeNum, EleNum, i*1.0, dt);
	}
	end = clock();
	printf("矩阵计算完成,用时 %d ms\n", end - start);

	// 3. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	mat_vec_product(MG[0], bd[1].lu, B);

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("Gauss Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		mat_vec_product(MG[0], bd[i].lu, B);

		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;

				mat_vec_product(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;

				mat_vec_product(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		mat_vec_product(MT[0], B, bd[i].lt);

		// 将 n = i时刻的未知量/已知量转化为位移和面力

		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	}

	//输出矩阵
//	if (m_flag == 0)
//	{
//
//		ofstream MTfile;
//		MTfile.open(".\\output\\MT0", ios::trunc);
//		MTfile.precision(15);
//		for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
//			for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
//				MTfile << MT[0](itemp + 1, jtemp + 1) << "\t";
//			}
//			MTfile << endl;
//		}
//		MTfile.close();
//		MTfile.clear();
//
//
//		MTfile.open(".\\output\\MG0", ios::trunc);
//		MTfile.precision(15);
//		for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
//			for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
//				MTfile << MG[0](itemp + 1, jtemp + 1) << "\t";
//
//			}
//			MTfile << endl;
//		}
//		MTfile.close();
//		MTfile.clear();
//
///*
//
//		char nameofstep[10], openfile[30], tempchar1[20] = ".\\output\\MT", tempchar2[20] = ".\\output\\MG";
//
//		for (i = 1; i <= MaxN; ++i)
//		{
//			strcpy_s(openfile, tempchar1);
//			_itoa_s(i, nameofstep, 10);
//			strcat_s(openfile, nameofstep);
//			MTfile.open(openfile, ios::trunc);
//			MTfile.precision(15);
//			for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
//				for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
//					MTfile << MT[i](itemp + 1, jtemp + 1) << "\t";
//				}
//				MTfile << endl;
//			}
//			MTfile.close();
//			MTfile.clear();
//
//			strcpy_s(openfile, tempchar2);
//			_itoa_s(i, nameofstep, 10);
//			strcat_s(openfile, nameofstep);
//			MTfile.open(openfile, ios::trunc);
//			MTfile.precision(15);
//			for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
//				for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
//					MTfile << MG[i](itemp + 1, jtemp + 1) << "\t";
//				}
//				MTfile << endl;
//			}
//			MTfile.close();
//			MTfile.clear();
//
//		}
//*/
//	}
	for (long i = 0; i <= MaxN; ++i)
	{
		GetInfoFromMatrixNstep(MG[i], MT[i], i, NodeNum);    //输出1~Ntsep步的MG、MT矩阵 
	}
	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

int DynaGaussSolverNewCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID, int m_flag)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;
	CSRMat* MGS;
	CSRMat* MTS;
	WmatrixCSR Mtemp;//临时矩阵
	WvectorCSR Btemp;


	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 确定最大的非零矩阵对应序号
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf("最大非零矩阵迭代步: %ld\n", MaxN);


	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	MT = new Wmatrix[1];
	//MG = new Wmatrix[NStep + 1];
	MTS = new CSRMat[NStep + 1];
	MGS = new CSRMat[NStep + 1];


	MT[0].initial(3 * NodeNum, 3 * NodeNum);
	Mtemp.initial(3 * NodeNum, 3 * NodeNum);//临时矩阵初始化
	Mtemp.reinitial();//赋值0
	Btemp.initial(3 * NodeNum);//临时向量初始化


	// 2. 计算矩阵
	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	GetMatrixDynaKnown0ij(m_DSE, Mtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);
	SparseConvert(Mtemp, MGS[0]);


	//FILE* fptemp;
	//fopen_s(&fptemp, "MT_D.txt", "w");
	//for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
	//	for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
	//		fprintf_s(fptemp, " %.8lf ", MT[0](itemp + 1, jtemp + 1));
	//	}
	//	fprintf_s(fptemp, "\n");
	//}
	//fclose(fptemp);
	//printf("\n输出未知系数矩阵元素完成!\n");
	//fopen_s(&fptemp, "MG_D.txt", "w");
	//for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
	//	for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
	//		fprintf_s(fptemp, " %.8lf ", MG[0](itemp + 1, jtemp + 1));
	//	}
	//	fprintf_s(fptemp, "\n");
	//}
	//fclose(fptemp);
	//printf("\n输出已知系数矩阵元素完成!\n");

	// 把MT[0]求逆
	int InvFlag = inv_mat(MT[0], MT[0].n);


	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		GetMatrixDynaT1ij(m_DSE, Mtemp, NodeNum, EleNum, i * 1.0, dt, m_InfElePID, m_ElePID);
		SparseConvert(Mtemp, MTS[i]);

		GetMatrixDynaT2ij(m_DSE, Mtemp, NodeNum, EleNum, i * 1.0, dt, m_InfElePID, m_ElePID);
		SparseAddition(Mtemp, MTS[i], Btemp);

		GetMatrixDynaUij(m_DSE, Mtemp, NodeNum, EleNum, i * 1.0, dt);
		SparseConvert(Mtemp, MGS[i]);
	}

	// 3. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	//mat_vec_product(MG[0], bd[1].lu, B);
	SparseMul(MGS[0], bd[1].lu, B);

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("Gauss Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		//mat_vec_product(MG[0], bd[i].lu, B);
		SparseMul(MGS[0], bd[i].lu, B);

		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;
				SparseMul(MGS[j], CB1, CB2);
				//mat_vec_product(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;
				SparseMul(MTS[j], CB1, CB2);
				//mat_vec_product(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		mat_vec_product(MT[0], B, bd[i].lt);

		// 将 n = i时刻的未知量/已知量转化为位移和面力

		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	}
	//for (long i = 1; i <= MaxN; ++i)
	//{
	//	GetInfoFromMatrixNstep(MG[i], MT[i], i, NodeNum);    //输出1~Ntsep步的MG、MT矩阵 
	//}


	// 释放空间

	delete[] MT;
	MT = 0;
	//delete[] MG;
	//MG = 0;
	delete[] MGS;
	MGS = 0;
	delete[] MTS;
	MTS = 0;
	Mtemp.finish();


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

int DynaGaussSolverNewCCSR(DSquareElement* m_DSE, BoundaryValue* bd, long NodeNum, long EleNum, long NStep, double MaxLength, long** m_InfElePID, long** m_ElePID, int m_flag)
{
	long i, j;

	Wmatrix* MT;
	Wmatrix* MG;
	SymCCSRMat* MGS;
	CCSRMat* MTS;
	WmatrixCSR Mtemp;//临时矩阵
	WmatrixCCSR MCCSRtemp;//临时矩阵
	WmatrixSymCCSR MSymCtemp;//临时矩阵
	WvectorCSR Btemp;
	WvectorCCSR BCCSRtemp;
	WvectorSymCCSR BSymCtemp;


	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	// 确定最大的非零矩阵对应序号
	long MaxN = NStep;

	for (i = 0; i <= NStep; ++i)
	{
		if (DSquareElement::m_DMat.C2Dt * i > MaxLength)
		{
			MaxN = i + 1;
			break;
		}
	}

	if (MaxN > NStep)
		MaxN = NStep;

	printf("最大非零矩阵迭代步: %ld\n", MaxN);


	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	MT = new Wmatrix[1];
	//MG = new Wmatrix[NStep + 1];
	MTS = new CCSRMat[NStep+1];
	MGS = new SymCCSRMat[NStep+1];


	MT[0].initial(3 * NodeNum, 3 * NodeNum);
	Mtemp.initial(3 * NodeNum, 3 * NodeNum);//临时矩阵初始化
	Mtemp.reinitial();//赋值0
	MCCSRtemp.initial(3 * NodeNum, 3 * NodeNum);//临时矩阵初始化
	MCCSRtemp.reinitial();//赋值0
	MSymCtemp.initial(3 * NodeNum, 3 * NodeNum);//临时矩阵初始化
	MSymCtemp.reinitial();//赋值0

	Btemp.initial(3 * NodeNum);//临时向量初始化
	BCCSRtemp.initial(3 * NodeNum);
	BSymCtemp.initial(3 * NodeNum);


	// 2. 计算矩阵
	GetMatrixDynaUnknown0ij(m_DSE, MT[0], NodeNum, EleNum, m_InfElePID, m_ElePID);
	GetMatrixDynaKnown0ij(m_DSE, MCCSRtemp, NodeNum, EleNum, m_InfElePID, m_ElePID);

	SparseConvert(MCCSRtemp, MTS[0]);

	//FILE* fptemp;
	//fopen_s(&fptemp, "MT_D.txt", "w");
	//for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
	//	for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
	//		fprintf_s(fptemp, " %.8lf ", MT[0](itemp + 1, jtemp + 1));
	//	}
	//	fprintf_s(fptemp, "\n");
	//}
	//fclose(fptemp);
	//printf("\n输出未知系数矩阵元素完成!\n");
	//fopen_s(&fptemp, "MG_D.txt", "w");
	//for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
	//	for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
	//		fprintf_s(fptemp, " %.8lf ", MG[0](itemp + 1, jtemp + 1));
	//	}
	//	fprintf_s(fptemp, "\n");
	//}
	//fclose(fptemp);
	//printf("\n输出已知系数矩阵元素完成!\n");

	// ————把MT[0]求逆————
	int InvFlag = inv_mat(MT[0], MT[0].n);
	

	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th BEM Matrices.\n", i);
		GetMatrixDynaT1ij(m_DSE, MCCSRtemp, NodeNum, EleNum, i * 1.0, dt, m_InfElePID, m_ElePID);
		SparseConvert(MCCSRtemp, MTS[i]);

		GetMatrixDynaT2ij(m_DSE, MCCSRtemp, NodeNum, EleNum, i * 1.0, dt, m_InfElePID, m_ElePID);
		SparseAddition(MCCSRtemp, MTS[i], BCCSRtemp);

		GetMatrixDynaUij(m_DSE, MSymCtemp, NodeNum, EleNum, i * 1.0, dt);
		SparseConvert(MSymCtemp, MGS[i]);
	}

	// 3. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;

	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("Gauss Solver: Step 1.\n");

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	//mat_vec_product(MG[0], bd[1].lu, B);
	SparseMul(MTS[0], bd[1].lu, B);

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中
	mat_vec_product(MT[0], B, bd[1].lt);

	// 将 n = 1时刻的未知量/已知量转化为位移和面力

	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("Gauss Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		//mat_vec_product(MG[0], bd[i].lu, B);
		SparseMul(MTS[0], bd[i].lu, B);

		// j = 1..i
		for (j = 1; j <= i; ++j)
		{
			if (j <= MaxN)
			{
				// [G]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lt;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lt, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lt;
				SparseMul(MGS[j], CB1, CB2);
				//mat_vec_product(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;
				SparseMul(MTS[j], CB1, CB2);
				//mat_vec_product(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		mat_vec_product(MT[0], B, bd[i].lt);

		// 将 n = i时刻的未知量/已知量转化为位移和面力

		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);

	}
	//for (long i = 1; i <= MaxN; ++i)
	//{
	//	GetInfoFromMatrixNstep(MG[i], MT[i], i, NodeNum);    //输出1~Ntsep步的MG、MT矩阵 
	//}


	// 释放空间

	delete[] MT;
	MT = 0;
	//delete[] MG;
	//MG = 0;
	delete[] MGS;
	MGS = 0;
	delete[] MTS;
	MTS = 0;
	Mtemp.finish();
	MCCSRtemp.finish();
	MSymCtemp.finish();


	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	std::string gaussTimePath = DBEMOutputPath("GaussTime.dat");
	fopen_s(&pf, gaussTimePath.c_str(), "w");
	if (InvFlag == 1)
		fprintf(pf, "高斯求解成功!\n");
	else
		fprintf(pf, "高斯求解失败!\n");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

int GetInfoFromMatrixNstep(Wmatrix& MG, Wmatrix& MT, long Step, long NodeNum) {
	//输出Step步的MG、MT矩阵
	FILE* fpMG, * fpMT;
	char nameofstep[10], openfileMG[30], openfileMT[30], tempMG[10] = "MG", tempMT[10] = "MT", tempchar[30] = ".//result//";
	strcpy_s(openfileMG, tempchar);
	strcpy_s(openfileMT, tempchar);
	strcpy_s(openfileMG, tempMG);
	strcpy_s(openfileMT, tempMT);
	_itoa_s(Step, nameofstep, 10);
	strcat_s(openfileMG, nameofstep);
	strcat_s(openfileMT, nameofstep);
	fopen_s(&fpMG, openfileMG, "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fpMG, " %.8lf ", MG(itemp + 1, jtemp + 1));
		}
		fprintf_s(fpMG, "\n");
	}
	fclose(fpMG);


	fopen_s(&fpMT, openfileMT, "w");
	for (int itemp = 0; itemp < 3 * NodeNum; itemp++) {
		for (int jtemp = 0; jtemp < 3 * NodeNum; jtemp++) {
			fprintf_s(fpMT, " %.8lf ", MT(itemp + 1, jtemp + 1));
		}
		fprintf_s(fpMT, "\n");
	}
	fclose(fpMT);
	return 1;
}
