#ifndef _POSTPRO_
#define _POSTPRO_
#include "Precompiler.h"
#include "geo.h"
#include "DSquareElement.h"

class OutPoint //总输出节点数=输入节点数-单元数*2
{
public:
	double pt[3] = {0,0,0};

	int EleList[5] = { 0,0,0,0,0 };//程序从0编号！
	int LocList[5] = { 0,0,0,0,0 };//左下->右上（顺时针）
	int flag = 0;

};

class OutEle
{
public:
	int PointList[4] = { 0,0,0,0 };//tecplot 节点从1开始编号！！！

};


int PreOutput(OutPoint* OutP, OutEle* OutE, DSquareElement* m_DSE, long EleNum, long PointNum,Point* m_PointList,long& numOut);

#endif