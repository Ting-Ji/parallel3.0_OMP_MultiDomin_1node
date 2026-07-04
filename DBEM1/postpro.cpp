#include "postpro.h"


int PreOutput(OutPoint* OutP, OutEle* OutE, DSquareElement* m_DSE,long EleNum,long PointNum, Point* m_PointList,long& numOut)
{
	long i, j, k;
	long OutTemp = 0;
	long tflag;//
	long TList[10];
	long TLocList[10];
	long tempNum=0;
	//！！！！！！！！OutP兜兵晒！！！！！！！！
	for (i = 0; i < PointNum; i++) {
		tflag = 0;
		for (j = 0; j < 10; j++) {
			TList[j] = 0;
			TLocList[j] = 0;
		}
		for (j = 0; j < EleNum; j++) {
			//if (temp = 4) 
				//break;
			for (k = 0; k < 4; k++) {
				if (i == m_DSE[j].m_pointID[k]) {
					TList[tflag] = j;
					TLocList[tflag] = k;
					tflag++;
					
				}
			}
		}
		if (tflag != 0) {
			OutP[OutTemp].pt[0] = m_PointList[i].pt[0];
			OutP[OutTemp].pt[1] = m_PointList[i].pt[1];
			OutP[OutTemp].pt[2] = m_PointList[i].pt[2];
			for (j = 0; j < tflag; j++) {
				OutP[OutTemp].EleList[j] = TList[j];
				OutP[OutTemp].LocList[j] = TLocList[j];
			}
			OutP[OutTemp].flag = tflag;
			OutTemp++;
			tempNum++;
		}
	}

	//！！！！！！！OutE兜兵晒！！！！！！！！
	long tempE,tempL;
	for (i = 0; i < tempNum; i++) {
		for (j = 0; j < OutP[i].flag; j++) {
			tempE = OutP[i].EleList[j];
			tempL = OutP[i].LocList[j];
			OutE[tempE].PointList[tempL] = i + 1;
		}
	}
	numOut = tempNum;
	return 0;
}
