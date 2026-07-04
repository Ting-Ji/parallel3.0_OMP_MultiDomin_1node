#ifndef _GEO_
#define _GEO_
#include "Precompiler.h"
//-------------------------------------------------------------------------------------------

//defination of a point
struct Point
{
	double pt[DIMENSION];
	Point& operator = (const Point& a)
	{pt[0]=a.pt[0];pt[1]=a.pt[1];pt[2]=a.pt[2];return *this;}
};

//defination of cube
struct Cube
{
	Point center;           //center of the cube
	double length;          //edge length of the cube
};

struct AniCube
{
	Point center;                   //center of the cube
	double length[DIMENSION];       //edge length of the cube
	AniCube& operator = (const AniCube& a)
	{
		center=a.center;
		for(int i=0;i<DIMENSION;++i)
			length[i] = a.length[i];
		return *this;
	}
};


//distance of two points
inline double Dist(Point& pt1,Point& pt2)
{
	return sqrt((pt1.pt[0]-pt2.pt[0])*(pt1.pt[0]-pt2.pt[0]) + (pt1.pt[1]-pt2.pt[1])*(pt1.pt[1]-pt2.pt[1]) + (pt1.pt[2]-pt2.pt[2])*(pt1.pt[2]-pt2.pt[2]));
}

//if point is inside cube. return 1 means yes，0 means no
int Isin(Point& pt,Cube& cb);

// 将节点用立方体包围
int AssignBox(Cube& m_cube,Point* m_point,long PointNum);

//determain cube position and size according to point list
int AssignCubeSize(Point* PtList,long PtNum,Cube& cb);

//determain cube position and size according to point list
int AssignCubeSize(Point* PtList,long PtNum,AniCube& cb); 

// get max length of the point list
double MaxLength(Point* PtList,long PtNum); 
double MaxLength_adapt(Point* PtList, long PtNum);

//create sub-cubes of a cube according to the position
void CreateSubCube(const Cube& father,Cube& child,int position);

//-------------------------------------------------------------------------------------------
#endif