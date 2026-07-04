#include "math.h"
#include "Geo.h"
//--------------------------------------------------------------------------------------



//if point is inside cube. return 1 means yes，0 means no
int Isin(Point& pt, Cube& cb)
{
	double ds[DIMENSION], l;
	int flag = 1;
	int i;
	for (i = 0; i < DIMENSION; i++)
		ds[i] = cb.center.pt[i] - pt.pt[i];
	l = 0.5*cb.length;
	for (i = 0; i < DIMENSION; i++)
	{
		if (ds[i] * ds[i] > l*l*1.0001)
		{
			flag = 0; break;
		}
	}
	return flag;
}

//create sub-cubes of a cube according to the position
void CreateSubCube(const Cube& father, Cube& child, int position)
{
	double half;
	int i, bit;
	child.length = 0.5*father.length;
	half = 0.5*child.length;
	bit = 1;
	for (i = 0; i < DIMENSION; i++)
	{
		child.center.pt[i] = father.center.pt[i] + ((bit&position) == 0 ? -1 : 1)*half;
		bit <<= 1;
	}
}

// 将节点用立方体包围
int AssignBox(Cube& m_cube, Point* m_point, long PointNum)
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double maxlength;
	long i;
	xmin = m_point[0].pt[0];
	ymin = m_point[0].pt[1];
	zmin = m_point[0].pt[2];
	xmax = xmin;
	ymax = ymin;
	zmax = zmin;
	for (i = 0; i < PointNum; ++i)
	{
		if (xmin > m_point[i].pt[0])
			xmin = m_point[i].pt[0];
		if (ymin > m_point[i].pt[1])
			ymin = m_point[i].pt[1];
		if (zmin > m_point[i].pt[2])
			zmin = m_point[i].pt[2];
		if (xmax < m_point[i].pt[0])
			xmax = m_point[i].pt[0];
		if (ymax < m_point[i].pt[1])
			ymax = m_point[i].pt[1];
		if (zmax < m_point[i].pt[2])
			zmax = m_point[i].pt[2];
	}
	maxlength = xmax - xmin;
	if (maxlength < (ymax - ymin))
		maxlength = ymax - ymin;
	if (maxlength < (zmax - zmin))
		maxlength = zmax - zmin;
	m_cube.center.pt[0] = xmin + 0.5*maxlength;
	m_cube.center.pt[1] = ymin + 0.5*maxlength;
	m_cube.center.pt[2] = zmin + 0.5*maxlength;
	m_cube.length = 1.02*maxlength;

	return 1;
}

//determain cube position and size according to point list
int AssignCubeSize(Point* PtList, long PtNum, Cube& cb)
{
	double rmin[DIMENSION], rmax[DIMENSION];
	double maxlength;
	long i, j;
	for (i = 0; i < DIMENSION; ++i)
	{
		rmin[i] = PtList[0].pt[i];
		rmax[i] = rmin[i];
	}
	for (i = 1; i < PtNum; ++i)
	{
		for (j = 0; j < DIMENSION; ++j)
		{
			if (rmin[j] > PtList[i].pt[j])
				rmin[j] = PtList[i].pt[j];
			if (rmax[j] < PtList[i].pt[j])
				rmax[j] = PtList[i].pt[j];
		}
	}
	maxlength = rmax[0] - rmin[0];
	for (i = 0; i < DIMENSION; ++i)
	{
		if (maxlength < (rmax[i] - rmin[i]))
			maxlength = rmax[i] - rmin[i];
		cb.center.pt[i] = rmin[i] + 0.5*(rmax[i] - rmin[i]);
	}
	cb.length = 1.02*maxlength;

	return 1;
}

//determain cube position and size according to point list
int AssignCubeSize(Point* PtList, long PtNum, AniCube& cb)
{
	double rmin[DIMENSION], rmax[DIMENSION];
	long i, j;
	for (i = 0; i < DIMENSION; ++i)
	{
		rmin[i] = PtList[0].pt[i];
		rmax[i] = rmin[i];
	}
	for (i = 1; i < PtNum; ++i)
	{
		for (j = 0; j < DIMENSION; ++j)
		{
			if (rmin[j] > PtList[i].pt[j])
				rmin[j] = PtList[i].pt[j];
			if (rmax[j] < PtList[i].pt[j])
				rmax[j] = PtList[i].pt[j];
		}
	}
	for (i = 0; i < DIMENSION; ++i)
	{
		cb.center.pt[i] = rmin[i] + 0.5*(rmax[i] - rmin[i]);
		cb.length[i] = 1.0*(rmax[i] - rmin[i]);
	}

	return 1;
}

// get max length of the point list
double MaxLength(Point* PtList, long PtNum)
{
	AniCube m_cube;
	AssignCubeSize(PtList, PtNum, m_cube);

	return sqrt(m_cube.length[0] * m_cube.length[0] + m_cube.length[1] * m_cube.length[1] + m_cube.length[2] * m_cube.length[2]);
}

// get max length of the point list 改进版
double MaxLength_adapt(Point* PtList, long PtNum)
{
	int i, j,k;
	double maxl=-1,templ=0;
	for (i = 0; i < PtNum; i++)
	{
		for (j = i + 1; j < PtNum; j++)
		{
			templ = 0.0;
			for (k = 0; k < DIMENSION; k++)
			{
				templ+=(PtList[i].pt[k] - PtList[j].pt[k])*(PtList[i].pt[k] - PtList[j].pt[k]);
			}
			if (maxl < templ)
				maxl = templ;
		}
	}

	return sqrt(maxl);
}

//--------------------------------------------------------------------------------------