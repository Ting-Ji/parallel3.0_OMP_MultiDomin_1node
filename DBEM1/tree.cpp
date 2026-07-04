#include "tree.h"

//------------------------------------------------------------------
//
//      树类的定义实现
//
//      该类与FMM-Tree类的最大区别在于：叶子是自身的邻居，用于ACA运算
//
//------------------------------------------------------------------

// initial the root of tree
void InitRoot(Tree* &m_tree, const Cube& m_cube, long totalpoint)
{
	//m_tree:        root waiting for initialization
	//m_cube:        geometry information of root
	//totalpoint:    total number of points contained in root

	long i;
	PointList* tempt;

	m_tree = new Tree;

	m_tree->m_Level = 1;
	m_tree->m_Father = 0;
	m_tree->Flag = 0;
	m_tree->m_Cube = m_cube;

	m_tree->m_PointList = new PointList;
	m_tree->m_PointList->PointID = 0;
	m_tree->m_PointList->next = 0;
	tempt = m_tree->m_PointList;

	for (i = 1; i < totalpoint; ++i)
	{
		tempt->next = new PointList;
		tempt->next->PointID = i;
		tempt = tempt->next;
		tempt->next = 0;
	}

	tempt = 0;
	m_tree->m_PointCount = totalpoint;

	//用于节点重新编号
	m_tree->m_BeginID = 0;
}

// create oct-tree structure according to the list of points
void CreateTree(Tree* &m_tree, Point* m_plist, long max_point)
{
	//m_tree:    root of tree
	//m_plist:   list of all points
	//max_point: the maxium number of points contained in a leaf

	//if a blank tree ,return
	if (m_tree == 0)
		return;

	long i, j;

	//if a leaf point
	if (m_tree->m_PointCount <= max_point)
	{
		m_tree->Flag = 1;
		for (i = 0; i < CHILDNUMBER; ++i)
			m_tree->m_Children[i] = 0;
		m_tree->m_Neibor = 0;
		m_tree->m_Interaction = 0;
		return;
	}

	m_tree->Flag = 0;
	//create children
	for (i = 0; i < CHILDNUMBER; ++i)
	{
		m_tree->m_Children[i] = new Tree;
		CreateSubCube(m_tree->m_Cube, m_tree->m_Children[i]->m_Cube, i);
		m_tree->m_Children[i]->m_Level = m_tree->m_Level + 1;
		m_tree->m_Children[i]->m_PointList = 0;
		m_tree->m_Children[i]->m_Father = m_tree;
		m_tree->m_Children[i]->m_PointCount = 0;
	}

	//localize the points contained to children
	PointList* tempt;
	int Inflag;
	double PtDist[CHILDNUMBER];
	double tmpvalue;
	int NearChild;
	while (m_tree->m_PointList)
	{
		Inflag = 0;
		for (i = 0; i < CHILDNUMBER; ++i)
		{
			if (Isin(m_plist[m_tree->m_PointList->PointID], m_tree->m_Children[i]->m_Cube))
			{
				Inflag = 1;
				NearChild = i;
				break;
			}
		}
		if (!Inflag)
		{
			//the point is not allocated successfully
			for (i = 0; i < CHILDNUMBER; ++i)
				PtDist[i] = Dist(m_plist[m_tree->m_PointList->PointID], m_tree->m_Children[i]->m_Cube.center);
			tmpvalue = PtDist[0];
			NearChild = 0;
			for (i = 1; i < CHILDNUMBER; ++i)
			{
				if (tmpvalue > PtDist[i])
				{
					tmpvalue = PtDist[i];
					NearChild = i;
				}
			}
		}
		//allocate point
		i = NearChild;
		if (m_tree->m_Children[i]->m_PointList)
		{
			tempt = m_tree->m_Children[i]->m_PointList;
			m_tree->m_Children[i]->m_PointList = m_tree->m_PointList;
			m_tree->m_PointList = m_tree->m_PointList->next;
			m_tree->m_Children[i]->m_PointList->next = tempt;
			m_tree->m_Children[i]->m_PointCount++;
		}
		else
		{
			m_tree->m_Children[i]->m_PointList = m_tree->m_PointList;
			m_tree->m_PointList = m_tree->m_PointList->next;
			m_tree->m_Children[i]->m_PointList->next = 0;
			m_tree->m_Children[i]->m_PointCount++;
		}

		if (m_tree->m_PointList == 0)
			break;
	}//end while

	tempt = 0;

	//delete children which contain no points
	j = m_tree->m_BeginID;
	for (i = 0; i < CHILDNUMBER; ++i)
	{
		if (m_tree->m_Children[i]->m_PointCount == 0)
		{
			delete m_tree->m_Children[i];
			m_tree->m_Children[i] = 0;
		}
		else //calculate m_beginID
		{
			m_tree->m_Children[i]->m_BeginID = j;
			j += m_tree->m_Children[i]->m_PointCount;
		}
	}

	//reset the neighbor and interaction list
	m_tree->m_Neibor = 0;
	m_tree->m_Interaction = 0;

	for (i = 0; i < CHILDNUMBER; ++i)
		CreateTree(m_tree->m_Children[i], m_plist, max_point);
}

//reset moments of a treenode
void ResetTree(Tree* &m_tree)
{
	if (m_tree == 0)
		return;

	// add code of resetting
	long i;
	for (i = 0; i < CHILDNUMBER; ++i)
		ResetTree(m_tree->m_Children[i]);
}

//delete the total oct-tree
void DeleteTree(Tree* &m_tree)
{
	if (m_tree == 0)
		return;
	int i;
	//delete children first
	for (i = 0; i < CHILDNUMBER; ++i)
		DeleteTree(m_tree->m_Children[i]);

	//delete list of points contained
	PointList* tempt;
	while (m_tree->m_PointList)
	{
		tempt = m_tree->m_PointList;
		m_tree->m_PointList = m_tree->m_PointList->next;
		delete tempt;
		tempt = 0;
	}

	//delete list of neighbors
	TreeNeighbor* nei;
	while (m_tree->m_Neibor)
	{
		nei = m_tree->m_Neibor;
		m_tree->m_Neibor = m_tree->m_Neibor->next;
		delete nei;
		nei = 0;
	}

	//delete interaction list
	while (m_tree->m_Interaction)
	{
		nei = m_tree->m_Interaction;
		m_tree->m_Interaction = m_tree->m_Interaction->next;
		delete nei;
		nei = 0;
	}

	delete m_tree;
	m_tree = 0;
}

//calculate the total number of leaves of the tree
void CalLeafNum(Tree* &m_tree, long &sum)
{
	//precondition: sum=0
	if (m_tree == 0)
		return;
	int i, flag;
	flag = 0;
	for (i = 0; i < CHILDNUMBER; ++i)
	{
		if (m_tree->m_Children[i] != 0)
		{
			CalLeafNum(m_tree->m_Children[i], sum);
			flag = 1;
		}
	}
	if (flag == 0)
		sum++;
}

//initialize the leaf pointer
void InitLeafPointer(PointerOfLeaf& p_leaf, Tree* &m_tree, long &p)
{
	//precondition:p=0
	if (m_tree == 0)
		return;
	int i, flag;
	PointList* tempt;
	flag = 0;
	for (i = 0; i < CHILDNUMBER; ++i)
	{
		if (m_tree->m_Children[i] != 0)
		{
			InitLeafPointer(p_leaf, m_tree->m_Children[i], p);
			flag = 1;
		}
	}
	//flag=0 means this treenode is a leaf
	if (flag == 0)
	{
		p_leaf.LeafPointer[p] = m_tree;
		p++;

		tempt = m_tree->m_PointList;
		while (tempt)
		{
			p_leaf.PointPointer[tempt->PointID] = m_tree;
			tempt = tempt->next;
		}
	}
	tempt = 0;
}

//Create the leaf pointer
void CreateLeafPointer(PointerOfLeaf& p_leaf, Tree* &m_tree, long totalcount)
{
	//Attention: p_leaf.result_pointer[i] might be zero, which means that the center
	//           of the ith res_p is beyond any leaf.

	long sum;
	sum = 0;
	CalLeafNum(m_tree, sum);
	p_leaf.LeafCount = sum;
	p_leaf.PointCount = totalcount;
	p_leaf.LeafPointer = new Tree*[p_leaf.LeafCount];
	p_leaf.PointPointer = new Tree*[p_leaf.PointCount];
	sum = 0;
	InitLeafPointer(p_leaf, m_tree, sum);
}

// delete the leaf_pointer
void DeleteLeafPointer(PointerOfLeaf& p_leaf)
{
	if (p_leaf.LeafPointer)
	{
		delete[] p_leaf.LeafPointer;
		p_leaf.LeafPointer = 0;
	}
	if (p_leaf.PointPointer)
	{
		delete[] p_leaf.PointPointer;
		p_leaf.PointPointer = 0;
	}
}

//Renumber the elements
void RenumberPointID(PointerOfLeaf& p_leaf, long* RePID)
{
	long i, j, beginID;
	PointList* m_list;

	for (i = 0; i < p_leaf.LeafCount; ++i)
	{
		m_list = p_leaf.LeafPointer[i]->m_PointList;
		beginID = p_leaf.LeafPointer[i]->m_BeginID;

		for (j = 0; j < p_leaf.LeafPointer[i]->m_PointCount; ++j)
		{
			RePID[beginID + j] = m_list->PointID;
			m_list = m_list->next;
		}
	}
}

// add treenode 2 to list of treenode 1 as neighbor or interaction list
void AddList(Tree* &m_tree1, Tree* m_tree2, int position)
{
	//m_tree1: treenode which add m_tree2 to its list
	//m_tree2: treenode which is added to list of m_tree1
	//position:list number,meaning:
	//0:  add to neighbor list
	//>0: add to interaction list
	if (m_tree1 == 0 || m_tree2 == 0)
		return;

	TreeNeighbor* tempt;

	if (position == 0)
	{
		if (m_tree1->m_Neibor == 0)
		{
			m_tree1->m_Neibor = new TreeNeighbor;
			m_tree1->m_Neibor->nei = m_tree2;
			m_tree1->m_Neibor->next = 0;
		}
		else
		{
			tempt = m_tree1->m_Neibor;
			m_tree1->m_Neibor = 0;
			m_tree1->m_Neibor = new TreeNeighbor;
			m_tree1->m_Neibor->nei = m_tree2;
			m_tree1->m_Neibor->next = tempt;
		}
	}
	else
	{
		if (m_tree1->m_Interaction == 0)
		{
			m_tree1->m_Interaction = new TreeNeighbor;
			m_tree1->m_Interaction->nei = m_tree2;
			m_tree1->m_Interaction->next = 0;
		}
		else
		{
			tempt = m_tree1->m_Interaction;
			m_tree1->m_Interaction = 0;
			m_tree1->m_Interaction = new TreeNeighbor;
			m_tree1->m_Interaction->nei = m_tree2;
			m_tree1->m_Interaction->next = tempt;
		}
	}
	tempt = 0;
}

//Judge if m_tree 2 is neighbor of m_tree1. If yes, return 0, else, return 1
int JudgePosition(Tree* m_tree1, Tree* m_tree2)
{
	//precondition: m_tree1!=0, m_tree2!=0, and the levels of this
	//two treenodes are equal to each other
	//return value means:
	//-1:  error
	//0:   2 is neighbor of 1
	//1:   2 is interaction list of 1

	//m_tree1: a treenode containing various kinds of lists
	//m_tree2: a treenode which will be added to m_tree1's lists

	if (m_tree1 == 0 || m_tree2 == 0 || m_tree1->m_Level != m_tree2->m_Level)
		return -1;
	double distance = Dist(m_tree1->m_Cube.center, m_tree2->m_Cube.center);
	//m_tree2 is neighbor of m_tree1
	if (distance <= 1.1*sqrt(1.0*DIMENSION)*m_tree1->m_Cube.length)
		return 0;
	else
		return 1;
}

void CreateNeiInter(Tree* &m_tree)
{
	if (m_tree == 0)
		return;

	int i, j, flag, position;
	TreeNeighbor* tempt;

	Tree* t_tree;

	if (m_tree->m_Father)//父节点的邻居的子节点是否与本节点相邻
	{
		tempt = m_tree->m_Father->m_Neibor;
		while (tempt)
		{
			flag = 0;
			for (i = 0; i < CHILDNUMBER; ++i)
			{
				if (tempt->nei->m_Children[i] != 0)
				{
					flag = 1;
					position = JudgePosition(m_tree, tempt->nei->m_Children[i]);
					if (position >= 0)
						AddList(m_tree, tempt->nei->m_Children[i], position);
				}
			}
			if (flag == 0)
				AddList(m_tree, tempt->nei, 0);
			tempt = tempt->next;
		}
		tempt = 0;
	}

	//让子节点相互为邻居
	for (i = 0; i < CHILDNUMBER; ++i)
	{
		if (m_tree->m_Children[i] != 0)
		{
			for (j = 0; j < CHILDNUMBER; ++j)
			{
				if (i != j)
					AddList(m_tree->m_Children[i], m_tree->m_Children[j], 0);
			}
		}
	}

	// 与FMM-Tree的不同：叶子自身是自身的邻居
	if (m_tree->Flag)
	{
		t_tree = m_tree;
		AddList(m_tree, t_tree, 0);
	}

	for (i = 0; i < CHILDNUMBER; ++i)
		CreateNeiInter(m_tree->m_Children[i]);
}

//return the actual maxium number of points contained in leaf
long MaxLeafPointnum(PointerOfLeaf& p_leaf)
{
	//p_leaf: information structure of leaf pointers
	long maxnum, i;
	maxnum = 0;

	for (i = 0; i < p_leaf.LeafCount; ++i)
	{
		if (maxnum <= p_leaf.LeafPointer[i]->m_PointCount)
			maxnum = p_leaf.LeafPointer[i]->m_PointCount;
	}
	return maxnum;
}

//return the actual maxium level of tree
long MaxTreeLevel(PointerOfLeaf& p_leaf)
{
	//p_leaf: information structure of leaf pointers
	long maxlevel, i;
	maxlevel = 0;

	for (i = 0; i < p_leaf.LeafCount; ++i)
	{
		if (maxlevel <= p_leaf.LeafPointer[i]->m_Level)
			maxlevel = p_leaf.LeafPointer[i]->m_Level;
	}
	return maxlevel;
}



//-----------------------------------------------------------------------------------------