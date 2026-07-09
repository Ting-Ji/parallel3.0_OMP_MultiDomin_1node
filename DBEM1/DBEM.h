#ifndef _DBEM
#define _DBEM

#include "geo.h"
#include "DSquareElement.h"
#include <windows.h>//输出内存使用量 头文件  Windows API
#include <psapi.h> //输出内存使用量 头文件
#include "mat_vec.h"
#include <pthread.h> // pthread support

#include"postpro.h"
//#include "PthreadCCSR.h"

// 根据几何节点信息获取物理点信息
int GetNodeInfo(Point* m_pointlist, Point* &m_NodeList, long** m_ElePID, long** m_EleNID, long PointNum, long EleNum, long& NodeNum, DSquareElement* m_DSE);
// 画网格函数
int Plot(Point* &m_point, long** &m_ElePID, char* name, long PN, long EN, int flag, double amplitude);
// 画结果函数
int ResultPlotNodeAverage_14VARS(DSquareElement* m_DSE, Point* m_PointList, long** m_ElePID, int* bdid, BoundaryValue& bd, long n, long PointNum, long NodeNum, long EleNum, const char* fp, double amplitude);
// 输出选定节点的位移/面力随时间变化信息
int GetInfoFromPoint_SingleColumn(Point* m_NodeList, BoundaryValue* bd, long NodeNum, long NStep, double dt, double dudt, double C_1D, double ColumnLength, long** m_ElePID);
// 获取单元平均尺寸：根据单元第1、2个几何节点距离平均得到
double GetEleSize(Point* m_PointList, long** m_ElePID, long PointNum, long EleNum);
// 获取无限大域内圆球表面径向位移
double GetRadiusDispFromSphere(double rou, double a, double v, double t, double c1);
// 获取圆球(1,0,0)位置的X方向位移，以代替径向位移
int TResultPlotNodeAverage_14VARS(DSquareElement* m_DSE, Point* m_PointList, long** m_ElePID, int* bdid, BoundaryValue& bd, long n, long PointNum, long NodeNum, long EleNum, const char* fp, double amplitude, long PointID, double& DispR);
// 获取圆球(1,0,0)位置节点编号
long GetPoint100ID(Point* m_pointlist, long PointNum);
// 获取上表面节点的位移/面力随时间变化信息
int GetInfoFromPoint_Uppersurface(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag, long NStep);
//获取节点的位移/面力随时间变化信息
int GetInfoFromPoint_All(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag, long NStep);
//获取上表面节点的位移/面力静力学问题解
int GetInfoFromPoint_UppersurfaceforStatic(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag);
//
int GetInfoFromPoint_onenode(Point* m_NodeList, BoundaryValue* bd, long NodeNum, long NStep, double dt, long** m_ElePID);
//获取节点的位移/面力静力学问题解
int GetInfoFromPoint_AllforStatic(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum);
int GetInfoFromPoint_UppersurfaceforStatic(Point* m_NodeList, BoundaryValue* bd, long** m_EleNID, long EleNum, char* eleflag, double gap);
void showMemoryInfo();
int GetInfoFromPoint_OutPut(OutPoint* OutP, OutEle* OutE, BoundaryValue* bd, long** m_EleNID, long OutPoNum, long EleNum, long NStep);


extern WmatrixCCSR MCCSRtemp;//临时矩阵
extern WmatrixSymCCSR MSymCtemp;//临时矩阵


extern DSquareElement* m_DSE;    // 单元类列表
extern long EleNum;              // 单元数
extern long NodeNum;             // 物理节点数	
extern long** m_ElePID;          // 单元几何节点编号
extern long InfEleNum;
extern long** m_InfElePID;
extern long StepNum, ConvoluNum;
extern double dt_p;
extern long thread_num;
//extern long* targs;
extern double*** AssistCoef;
extern double*** AssistUij;
extern double*** AssistTij;
extern Wmatrix* MT;
extern WvectorCCSR BCCSRtemp;
extern SymCCSRMat* MGS;
extern CCSRMat* MTS;
extern CCSRMat MTS0;
extern BoundaryValue* bd;
extern Wvector* Bvector, CB1, CB2, x0, * CBG, * CBT, Bsum;

extern PthArg* Thread;
extern PthArg** SpMVThreadT, ** SpMVThreadG;
extern PthConvol* ConvolThreadT, *ConvolThreadG;
extern long* Pthnum;
extern long*** FlagUij;



int DBEM_main(int m_flag);


#endif
