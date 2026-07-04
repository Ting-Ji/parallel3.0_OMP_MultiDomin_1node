#ifndef _TREE
#define _TREE

#include "geo.h"

//-----------------------------------------------
//
//      树类的声明
//
//-----------------------------------------------


// list of point
struct PointList
{
	long PointID;
	PointList* next;
};

struct TreeNeighbor;

// tree structure
struct Tree
{
	Tree* m_Father;                        // pointer to father
	Tree* m_Children[CHILDNUMBER];         // pointers to children
	int Flag;                              // flag=1: leaf; flag=0: non_leaf
	Cube m_Cube;                           // geometry information
	PointList* m_PointList;                // list of points contained
	long m_PointCount;                     // total number of points contained
	long m_BeginID;                        // begin point ID
	TreeNeighbor* m_Neibor;                // list of neighbors
	TreeNeighbor* m_Interaction;           // interaction list
    long m_Level;                          // level
};

// tree neighbor or interaction list structure
struct TreeNeighbor
{
	Tree* nei;
	TreeNeighbor* next;
};

// information of pointer to leaves
struct PointerOfLeaf
{
	long LeafCount;
	long PointCount;
	long MaxLeafPointCount;
	Tree** LeafPointer;
	Tree** PointPointer;
};



//------------------------------------------------------------------------------------------
//     This part is the declaration of operations of Tree code.
//------------------------------------------------------------------------------------------

// initial the root of tree
void InitRoot(Tree* &m_tree,const Cube& m_cube,long totalpoint);

// create oct-tree structure according to the list of points
void CreateTree(Tree* &m_tree,Point* m_plist,long max_point);

// reset moments of a treenode
void ResetTree(Tree* &m_tree);

// delete the total oct-tree
void DeleteTree(Tree* &m_tree);

// calculate the total number of leaves of the tree
void CalLeafCount(Tree* &m_tree,long &sum);

// initialize the leaf pointer
void InitLeafPointer(PointerOfLeaf& p_leaf,Tree* &m_tree,long &p);

// Create the leaf pointer
void CreateLeafPointer(PointerOfLeaf& p_leaf,Tree* &m_tree,long totalcount);

// delete the leaf_pointer
void DeleteLeafPointer(PointerOfLeaf& p_leaf);

// add treenode 2 to list of treenode 1 as neighbor or interaction list
void AddList(Tree* &m_tree1,Tree* m_tree2,int position);

// Judge if m_tree 2 is neighbor of m_tree1. If yes, return 0, else, return interaction list code of m_tree2 to m_tree1
int JudgePosition(Tree* m_tree1,Tree* m_tree2);

//Create neighbor and interaction lists of a tree structure
void CreateNeiInter(Tree* &m_tree);

//return the actual maxium number of points contained in leaf
long MaxLeafPointCount(PointerOfLeaf& p_leaf);

//return the actual maxium level of tree
long MaxTreeLevel(PointerOfLeaf& p_leaf);

//Renumber the points
void RenumberPointID(PointerOfLeaf& p_leaf,long* RePID);


//-----------------------------------------------------------------------------------------
#endif