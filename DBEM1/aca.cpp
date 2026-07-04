#include "aca.h"
#include "DBEM.h"
#include <string>
#include "Counter.h"
#include <direct.h>
long* ACAEntry::IID;
long* ACAEntry::VID;
int* ACAEntry::VID_RC;
long ACAEntry::MAX_M;
long ACAEntry::MAX_N;

double ACAEntry::errorvalue;
double ACAEntry::errorstop;

long* ACAMatrix::RePID;
int ACAMatrix::Flag_ReadWrite;


ACAEntry::ACAEntry()
{
	flag_zero = 1; // 默认是非零矩阵

	flag = 0;  // 默认是邻居项
	flag_overflow = 1;
	mbeginID = 0;
	nbeginID = 0;
	m = 0;
	n = 0;
	num = 0;
	u = 0;
	v = 0;
}

ACAEntry::~ACAEntry()
{
	if (u)
	{
		delete[] u;
		u = 0;
	}

	if (v)
	{
		delete[] v;
		v = 0;
	}

	Entry.finish();
}

double ACAEntry::ACAMemory()
{
	double Mem;

	if (flag_zero == 0) // 零矩阵
	{
		Mem = 0.0;
	}
	else if (flag == 1 && flag_overflow == 1)
	{
		Mem = num * (m + n) * 3 * 8.0;
	}
	else
	{
		Mem = m * n * 3 * 3 * 8.0;
	}

	Mem /= (1024 * 1024);

	return Mem;
}

ACAEntry& ACAEntry::operator+=(ACAEntry& o_ACAEntry)
{
	long i, j, k;

	Wvector* tu;
	Wvector* tv;

	if (flag_zero == 0) // 零矩阵
	{
		if (o_ACAEntry.flag_zero == 0)
		{
			;
		}
		else if (o_ACAEntry.flag == 1 && o_ACAEntry.flag_overflow == 1)
		{
			flag_zero = 1;
			flag_overflow = 1;
			num = o_ACAEntry.num;

			u = new Wvector[num];
			v = new Wvector[num];

			for (k = 0; k < num; ++k)
			{
				u[k].initial(o_ACAEntry.u[k].n);
				u[k] = o_ACAEntry.u[k];
				v[k].initial(o_ACAEntry.v[k].n);
				v[k] = o_ACAEntry.v[k];
			}
		}
		else
		{
			flag_zero = 1;
			flag_overflow = 0;

			Entry.initial(o_ACAEntry.Entry.m, o_ACAEntry.Entry.n);
			Entry = o_ACAEntry.Entry;
		}
	}
	else if (flag == 1 && flag_overflow == 1) // ACA矩阵
	{
		if (o_ACAEntry.flag_zero == 0)
		{
			;
		}
		else if (o_ACAEntry.flag == 1 && o_ACAEntry.flag_overflow == 1)
		{
			// 方案1：采取加长u和v列表的方法
			/*
			tu = u;
			tv = v;

			u = new Wvector[num + o_ACAEntry.num];
			v = new Wvector[num + o_ACAEntry.num];

			for(k = 0; k < num; ++k)
			{
				u[k].initial(tu[k].n);
				u[k] = tu[k];
				v[k].initial(tv[k].n);
				v[k] = tv[k];
			}
			for(k = 0; k < o_ACAEntry.num; ++k)
			{
				u[k + num].initial(o_ACAEntry.u[k].n);
				u[k + num] = o_ACAEntry.u[k];
				v[k + num].initial(o_ACAEntry.v[k].n);
				v[k + num] = o_ACAEntry.v[k];
			}

			num += o_ACAEntry.num;

			if(tu)
			{
				delete[] tu;
				tu = 0;
			}
			if(tv)
			{
				delete[] tv;
				tv = 0;
			}
			*/
			// 方案2：重新计算ACA
			Entry.initial(m * 3, n * 3, 0);

			for (k = 0; k < num; ++k)
			{
				for (i = 0; i < Entry.m; ++i)
				{
					for (j = 0; j < Entry.n; ++j)
						Entry.a[i][j] += u[k].b[i] * v[k].b[j];
				}
			}

			for (k = 0; k < o_ACAEntry.num; ++k)
			{
				for (i = 0; i < Entry.m; ++i)
				{
					for (j = 0; j < Entry.n; ++j)
						Entry.a[i][j] += o_ACAEntry.u[k].b[i] * o_ACAEntry.v[k].b[j];
				}
			}

			if (u)
			{
				delete[] u;
				u = 0;
			}
			if (v)
			{
				delete[] v;
				v = 0;
			}

			FullACA();
		}
		else
		{
			Entry.initial(o_ACAEntry.Entry.m, o_ACAEntry.Entry.n);
			Entry = o_ACAEntry.Entry;

			for (k = 0; k < num; ++k)
			{
				for (i = 0; i < Entry.m; ++i)
				{
					for (j = 0; j < Entry.n; ++j)
						Entry.a[i][j] += u[k].b[i] * v[k].b[j];
				}
			}

			if (u)
			{
				delete[] u;
				u = 0;
			}
			if (v)
			{
				delete[] v;
				v = 0;
			}

			flag_overflow = 0;
		}

	}
	else // 常规矩阵
	{
		if (o_ACAEntry.flag_zero == 0)
		{
			;
		}
		else if (o_ACAEntry.flag == 1 && o_ACAEntry.flag_overflow == 1)
		{
			for (k = 0; k < o_ACAEntry.num; ++k)
			{
				for (i = 0; i < Entry.m; ++i)
				{
					for (j = 0; j < Entry.n; ++j)
						Entry.a[i][j] += o_ACAEntry.u[k].b[i] * o_ACAEntry.v[k].b[j];
				}
			}
		}
		else
			Entry += o_ACAEntry.Entry;
	}

	return *this;
}

int ACAEntry::FullACA()
{
	// 前提：m_ACAEntry.Entry计算完毕

	long i, j, k, ik, jk, min, mk, nk;
	double max, FB, FR;
	double comp;

	Wmatrix R(Entry);

	mk = R.m;
	nk = R.n;

	comp = mk * nk*1.0 / (mk + nk);
	flag_overflow = 1;   // 默认值得进行ACA

	min = mk > nk ? nk : mk;

	Wvector* tu = new Wvector[min];
	Wvector* tv = new Wvector[min];

	for (i = 0; i < min; ++i)
	{
		tu[i].initial(R.m);
		tv[i].initial(R.n);
	}

	normF(R, FB);

	num = 0;

	for (k = 0; k < min; ++k)
	{
		ArgMax(R, ik, jk, max);

		if (fabs(max) < errorvalue)    //矩阵为0
			break;

		num = k + 1;

		if (num > comp)  // 不值得用ACA
		{
			flag_overflow = 0;
			break;
		}

		for (i = 0; i < mk; ++i)
		{
			tu[k].b[i] = R.a[i][jk];
		}

		for (i = 0; i < nk; ++i)
		{
			tv[k].b[i] = R.a[ik][i] / max;
		}

		FR = 0.0;

		for (i = 0; i < mk; ++i)
		{
			for (j = 0; j < nk; ++j)
			{
				R.a[i][j] -= tu[k].b[i] * tv[k].b[j];

				FR += R.a[i][j] * R.a[i][j];
			}
		}

		FR = sqrt(FR);

		if (FR <= errorstop * FB)
			break;

	}

	//	printf("Fully: m:%ld  n:%ld  num:%ld\n",R.m,R.n,m_ACAEntry.num);

	if (num > 0 && flag_overflow == 1)
	{
		u = new Wvector[num];
		v = new Wvector[num];

		for (i = 0; i < num; ++i)
		{
			u[i].initial(tu[i].n);
			u[i] = tu[i];

			v[i].initial(tv[i].n);
			v[i] = tv[i];
		}

		Entry.finish();
	}

	if (tu)
	{
		delete[] tu;
		tu = 0;
	}
	if (tv)
	{
		delete[] tv;
		tv = 0;
	}

	return 1;

}

int ACAEntry::WriteToFile(FILE* fp,long k)
{
	long i, j;

	long ZeNum;//零元素个数
	long NoZeNum;//非零元素个数

	ZeNum = NoZeNum = 0;


	//fprintf(fp, "%d %d %d %ld %ld %ld\n", flag_zero, flag, flag_overflow, m, n, k);

	if (flag_zero == 0) {
		//ZeNum = m * n * 9;
		fprintf(fp, "%d %d \n", ZeNum, NoZeNum);
		return 1;
	}

		

	if (flag == 1 && flag_overflow == 1)
	{
		//fprintf(fp, "%ld\n", num);

		for (i = 0; i < num; ++i)
		{
			for (j = 0; j < u[0].n; ++j) {
				if (u[i].b[j] == 0.0) {
					ZeNum++;
				}
				else {
					NoZeNum++;
				}
			}
				//fprintf(fp, "%20.10e ", u[i].b[j]);
			//fprintf(fp, "\n");
			for (j = 0; j < v[0].n; ++j) {
				if (v[i].b[j] == 0.0) {
					ZeNum++;
				}
				else {
					NoZeNum++;
				}
			}
				//fprintf(fp, "%20.10e ", v[i].b[j]);
			//fprintf(fp, "\n");
		}

		fprintf(fp, "%d %d \n", ZeNum, NoZeNum);
	}
	else
	{
		for (i = 0; i < Entry.m; ++i)
		{
			for (j = 0; j < Entry.n; ++j)
			{
				if (Entry.a[i][j] == 0.0) {
					ZeNum++;
				}
				else {
					NoZeNum++;
				}
				//fprintf(fp, "%20.10e ", Entry.a[i][j]);
			}
			
		}

		fprintf(fp, "%d %d \n", ZeNum, NoZeNum);
	}

	/*fprintf(fp, "%d %d %d %ld %ld %ld\n", flag_zero, flag, flag_overflow, m, n,k);

	if (flag_zero == 0)
		return 1;

	if (flag == 1 && flag_overflow == 1)
	{
		fprintf(fp, "%ld\n", num);

		for (i = 0; i < num; ++i)
		{
			for (j = 0; j < u[0].n; ++j)
				fprintf(fp, "%20.10e ", u[i].b[j]);
			fprintf(fp, "\n");
			for (j = 0; j < v[0].n; ++j)
				fprintf(fp, "%20.10e ", v[i].b[j]);
			fprintf(fp, "\n");
		}
	}
	else
	{
		for (i = 0; i < Entry.m; ++i)
		{
			for (j = 0; j < Entry.n; ++j)
			{
				fprintf(fp, "%20.10e ", Entry.a[i][j]);
			}
			fprintf(fp, "\n");
		}
	}*/

	if (u)
	{
		delete[] u;
		u = 0;
	}

	if (v)
	{
		delete[] v;
		v = 0;
	}

	Entry.finish();

	return 1;
}

int ACAEntry::ReadFromFile(FILE* fp)
{
	// 前提：Entry, u, v均没有初始化
	long i, j;

	fscanf_s(fp, "%d", &flag_zero);
	fscanf_s(fp, "%d", &flag);
	fscanf_s(fp, "%d", &flag_overflow);
	fscanf_s(fp, "%ld", &m);
	fscanf_s(fp, "%ld", &n);

	if (flag_zero == 0)
		return 1;

	if (flag == 1 && flag_overflow == 1)
	{
		fscanf_s(fp, "%ld", &num);

		u = new Wvector[num];
		v = new Wvector[num];

		for (i = 0; i < num; ++i)
		{
			u[i].initial(m * 3);
			v[i].initial(n * 3);

			for (j = 0; j < u[0].n; ++j)
				fscanf_s(fp, "%lf", &u[i].b[j]);
			for (j = 0; j < v[0].n; ++j)
				fscanf_s(fp, "%lf", &v[i].b[j]);
		}
	}
	else
	{
		Entry.initial(m * 3, n * 3);

		for (i = 0; i < Entry.m; ++i)
		{
			for (j = 0; j < Entry.n; ++j)
			{
				fscanf_s(fp, "%lf", &Entry.a[i][j]);
			}
		}
	}

	return 1;
}

int ACAEntry::StaticInit(long PointCount, double o_errorvalue, double o_errorstop)
{
	long i;

	MAX_M = 1;
	MAX_N = 1;

	IID = new long[ PointCount];
	VID = new long[ PointCount];
	VID_RC = new int[ PointCount];

	for (i = 0; i < 1 * PointCount; ++i)
	{
		IID[i] = 3 * i;
		VID[i] = i / 3;
		VID_RC[i] = i - VID[i] * 3;
	}

	errorvalue = o_errorvalue;
	errorstop = o_errorstop;

	return 1;
}

int ACAEntry::StaticRelease()
{
	if (IID)
	{
		delete[] IID;
		IID = 0;
	}
	if (VID)
	{
		delete[] VID;
		VID = 0;
	}
	if (VID_RC)
	{
		delete[] VID_RC;
		VID_RC = 0;
	}

	return 1;
}

ACAMatrix& ACAMatrix::operator+=(ACAMatrix& o_ACAMatrix)
{
	long i;

	if (Flag_ReadWrite == 1) // 首先读入矩阵
	{
		ReadFromFile();
		o_ACAMatrix.ReadFromFile();
	}

	for (i = 0; i < EntryCount; ++i)
		m_ACAEntry[i] += o_ACAMatrix.m_ACAEntry[i];

	if (Flag_ReadWrite == 1) // 首先读入矩阵
	{
		WriteToFile();
		o_ACAMatrix.WriteToFile();
	}

	return *this;
}

int ACAMatrix::Memory(double& ACAMem, double& ClassicMem)
{
	long i;

	ACAMem = 0.0;
	ClassicMem = 0.0;

	for (i = 0; i < EntryCount; ++i)
	{
		ACAMem += m_ACAEntry[i].ACAMemory();
		ClassicMem += m_ACAEntry[i].m * m_ACAEntry[i].n * 3.0 * 3.0 * 8.0 / 1024.0 / 1024.0;
	}

	return 1;
}

int ACAMatrix::WriteToFile(char* filename, long ID)
{
	long i;

	sprintf_s(MatrixName, "%s_%ld.dat", filename, ID);

	FILE* fp;
	fopen_s(&fp, MatrixName, "w");

	for (i = 0; i < EntryCount; ++i)
		m_ACAEntry[i].WriteToFile(fp,i);

	fclose(fp);

	return 1;
}

int ACAMatrix::WriteToFile()
{
	long i;

	// TEST
	printf("Write %s\n", MatrixName);

	FILE* fp;
	fopen_s(&fp, MatrixName, "w");

	for (i = 0; i < EntryCount; ++i)
		m_ACAEntry[i].WriteToFile(fp,i);

	fclose(fp);

	return 1;
}

int ACAMatrix::ReadFromFile()
{
	long i;

	// TEST
	printf("Read %s\n", MatrixName);

	FILE* fp;
	fopen_s(&fp, MatrixName, "r");

	for (i = 0; i < EntryCount; ++i)
		m_ACAEntry[i].ReadFromFile(fp);

	fclose(fp);

	return 1;
}

int ACAGMRESPreConditioner(DSquareElement* m_DSE, ACAPreConditioner& m_Pre, Point* m_PointList, long Count_PT, long MaxLeafPointCount)
{
	Cube GlobalBox;
	Tree* m_treePre;
	PointerOfLeaf p_leafPre;

	long i, j, tj, k, beginID;

	Wmatrix subPre(3, 3);

	long PointCount = Count_PT;

	// 计算包含点列表的立方体
	AssignCubeSize(m_PointList, PointCount, GlobalBox);

	// 生成树结构
	InitRoot(m_treePre, GlobalBox, PointCount);
	CreateTree(m_treePre, m_PointList, MaxLeafPointCount);
	CreateLeafPointer(p_leafPre, m_treePre, PointCount);

	// 预处理类
	m_Pre.m_LeafCount = p_leafPre.LeafCount;
	m_Pre.m_LeafPointCount = new long[p_leafPre.LeafCount];
	m_Pre.m_LeafBeginID = new long[p_leafPre.LeafCount];
	m_Pre.m_RePID = new long[PointCount];
	m_Pre.IID = new long[PointCount];
	m_Pre.m_PreM = new Wmatrix[p_leafPre.LeafCount];

	for (i = 0; i < p_leafPre.LeafCount; ++i)
	{
		m_Pre.m_LeafPointCount[i] = p_leafPre.LeafPointer[i]->m_PointCount;
		m_Pre.m_LeafBeginID[i] = p_leafPre.LeafPointer[i]->m_BeginID;
	}

	// 编号重新排序
	RenumberPointID(p_leafPre, m_Pre.m_RePID);

	// 给IID赋值
	for (i = 0; i < PointCount; ++i)
		m_Pre.IID[i] = 3 * i;

	// 初始化预处理矩阵
	for (i = 0; i < p_leafPre.LeafCount; ++i)
	{
		j = 3 * m_Pre.m_LeafPointCount[i];
		m_Pre.m_PreM[i].initial(j, j);
	}

	// 计算预处理矩阵
	for (i = 0; i < m_Pre.m_LeafCount; ++i)
	{
		beginID = m_Pre.m_LeafBeginID[i];

		for (j = 0; j < m_Pre.m_LeafPointCount[i]; ++j)
		{
			tj = m_Pre.m_RePID[beginID + j];

			for (k = 0; k < m_Pre.m_LeafPointCount[i]; ++k)
			{
				NodeUnknownSubMatrix(m_DSE, tj, m_Pre.m_RePID[beginID + k]);
				subPre = DSquareElement::m_AssistCoef[0];
				m_Pre.m_PreM[i].assign(subPre, 3 * j + 1, 3 * k + 1, 3, 3);
			}
		}
	}

	// 预处理类的预处理矩阵求逆运算
	for (i = 0; i < m_Pre.m_LeafCount; ++i)
		inv_mat(m_Pre.m_PreM[i], m_Pre.m_PreM[i].m);

	//release PreTree and p_leafPre
	DeleteTree(m_treePre);
	DeleteLeafPointer(p_leafPre);

	return 1;
}

int ACA_MX(ACAPreConditioner& m_Pre, Wvector& X, Wvector& B)
{
	int m, n;
	long i, j, k, tj, tk;
	long a1, a2;

	B = 0.0;

	for (i = 0; i < m_Pre.m_LeafCount; ++i)
	{
		for (j = 0; j < m_Pre.m_LeafPointCount[i]; ++j)
		{
			tj = m_Pre.IID[j];
			a1 = m_Pre.IID[m_Pre.m_RePID[m_Pre.m_LeafBeginID[i] + j]];
			for (k = 0; k < m_Pre.m_LeafPointCount[i]; ++k)
			{
				tk = m_Pre.IID[k];
				a2 = m_Pre.IID[m_Pre.m_RePID[m_Pre.m_LeafBeginID[i] + k]];
				for (m = 1; m <= 3; ++m)
				{
					for (n = 1; n <= 3; ++n)
						B[a1 + m] += m_Pre.m_PreM[i](tj + m, tk + n)*X[a2 + n];
				}
			}
		}
	}

	return 1;
}

void ACAInfo(Tree* &m_tree, long& sum, ACAMatrix& m_ACAMatrix)
{
	// 给ACA矩阵列表赋值：基本属性
	// 前提：ACA矩阵列表已经生成相应空间

	if (m_tree == 0)
		return;

	int i;

	TreeNeighbor* tempt;

	if (m_tree->Flag) //叶子节点
	{
		tempt = m_tree->m_Neibor;
		while (tempt)
		{
			m_ACAMatrix.m_ACAEntry[sum].flag = 0;
			m_ACAMatrix.m_ACAEntry[sum].m = m_tree->m_PointCount;
			m_ACAMatrix.m_ACAEntry[sum].mbeginID = m_tree->m_BeginID;
			m_ACAMatrix.m_ACAEntry[sum].n = tempt->nei->m_PointCount;
			m_ACAMatrix.m_ACAEntry[sum].nbeginID = tempt->nei->m_BeginID;
			sum++;
			tempt = tempt->next;
		}
	}

	tempt = m_tree->m_Interaction;
	while (tempt)
	{
		m_ACAMatrix.m_ACAEntry[sum].flag = 1;
		m_ACAMatrix.m_ACAEntry[sum].flag_overflow = 1;
		m_ACAMatrix.m_ACAEntry[sum].m = m_tree->m_PointCount;
		m_ACAMatrix.m_ACAEntry[sum].mbeginID = m_tree->m_BeginID;
		m_ACAMatrix.m_ACAEntry[sum].n = tempt->nei->m_PointCount;
		m_ACAMatrix.m_ACAEntry[sum].nbeginID = tempt->nei->m_BeginID;
		sum++;
		tempt = tempt->next;
	}

	tempt = 0;

	for (i = 0; i < CHILDNUMBER; i++)
		ACAInfo(m_tree->m_Children[i], sum, m_ACAMatrix);
}

void InitACAEntry(ACAMatrix& m_ACAMatrix)
{
	// 给ACA矩阵初始化
	// 前提：ACA矩阵列表已经生成相应空间
	long i;

	for (i = 0; i < m_ACAMatrix.EntryCount; ++i)
	{

		/*
		if(m_ACAMatrix.m_ACAEntry[i].flag == 0) // 邻居
		{
			m_ACAMatrix.m_ACAEntry[i].Entry.initial(m_ACAMatrix.m_ACAEntry[i].m*24,m_ACAMatrix.m_ACAEntry[i].n*24);
		}
		else // 远节点
		{
			m_ACAMatrix.m_ACAEntry[i].Entry.m = m_ACAMatrix.m_ACAEntry[i].m*24;
			m_ACAMatrix.m_ACAEntry[i].Entry.n = m_ACAMatrix.m_ACAEntry[i].n*24;
		}
		*/

		m_ACAMatrix.m_ACAEntry[i].Entry.m = m_ACAMatrix.m_ACAEntry[i].m * 3;
		m_ACAMatrix.m_ACAEntry[i].Entry.n = m_ACAMatrix.m_ACAEntry[i].n * 3;

		if (ACAEntry::MAX_M < m_ACAMatrix.m_ACAEntry[i].Entry.m)
			ACAEntry::MAX_M = m_ACAMatrix.m_ACAEntry[i].Entry.m;
		if (ACAEntry::MAX_N < m_ACAMatrix.m_ACAEntry[i].Entry.n)
			ACAEntry::MAX_N = m_ACAMatrix.m_ACAEntry[i].Entry.n;
	}
}

int Setup_ACAEntry_RePID(ACAMatrix& m_ACAMatrix, Point* m_PointList, long Count_PT, long MaxLeafPointCount)
{
	// 设置ACA矩阵列表和RePID单元重新编号数组，进行相应的初始化和分配空间

	long i;
	Cube GlobalBox;
	Tree* m_tree;
	PointerOfLeaf p_leaf;

	// 计算包含点列表的立方体
	AssignCubeSize(m_PointList, Count_PT, GlobalBox);

	// 生成树结构
	InitRoot(m_tree, GlobalBox, Count_PT);
	CreateTree(m_tree, m_PointList, MaxLeafPointCount);
	CreateLeafPointer(p_leaf, m_tree, Count_PT);
	CreateNeiInter(m_tree);

	//重新编号
	ACAMatrix::RePID = new long[Count_PT];
	RenumberPointID(p_leaf, ACAMatrix::RePID);

	//计算树中的ACA矩阵数量
	m_ACAMatrix.EntryCount = 0;
	ACAMatrixCount(m_tree, m_ACAMatrix.EntryCount);

	//确定ACA矩阵信息
	m_ACAMatrix.m_ACAEntry = new ACAEntry[m_ACAMatrix.EntryCount];
	i = 0;
	ACAInfo(m_tree, i, m_ACAMatrix);

	//初始化ACA矩阵
	InitACAEntry(m_ACAMatrix);

	DeleteTree(m_tree);
	DeleteLeafPointer(p_leaf);

	return 1;
}

int GetInfoFromACAEntry(ACAMatrix& m_ACAMatrix1, ACAMatrix& m_ACAMatrix2)
{
	// 利用ACAMatrix1对ACAMatrix2进行相应的初始化和分配空间

	long i;

	//计算树中的ACA矩阵数量
	m_ACAMatrix2.EntryCount = m_ACAMatrix1.EntryCount;

	//确定ACA矩阵信息
	m_ACAMatrix2.m_ACAEntry = new ACAEntry[m_ACAMatrix2.EntryCount];
	for (i = 0; i < m_ACAMatrix2.EntryCount; ++i)
	{
		m_ACAMatrix2.m_ACAEntry[i].flag = m_ACAMatrix1.m_ACAEntry[i].flag;
		m_ACAMatrix2.m_ACAEntry[i].flag_overflow = m_ACAMatrix1.m_ACAEntry[i].flag_overflow;
		m_ACAMatrix2.m_ACAEntry[i].m = m_ACAMatrix1.m_ACAEntry[i].m;
		m_ACAMatrix2.m_ACAEntry[i].mbeginID = m_ACAMatrix1.m_ACAEntry[i].mbeginID;
		m_ACAMatrix2.m_ACAEntry[i].n = m_ACAMatrix1.m_ACAEntry[i].n;
		m_ACAMatrix2.m_ACAEntry[i].nbeginID = m_ACAMatrix1.m_ACAEntry[i].nbeginID;
	}

	//初始化ACA矩阵
	InitACAEntry(m_ACAMatrix2);

	return 1;
}

void CalACAEntry(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, int flag)
{
	// 计算ACA矩阵列表，选择Fully ACA
	// 前提：1. ACA矩阵列表初始化完毕
	//       2. 如果ACA矩阵直接计算，则直接计算矩阵已经分配空间
	//          如果ACA矩阵使用ACA分解，则其中的矩阵和ACA向量均没有分配空间
	//       3. RePID数组初始化和计算完毕

	long num, i, i24, j, j3, k, k24, it, jt, mbegin, nbegin, ik, jk, mk, nk;

	Wmatrix sub(3, 3);

	int Cal_Flag, Sum_Flag;

	for (num = 0; num < m_ACAMatrix.EntryCount; ++num)
	{
		mk = m_ACAMatrix.m_ACAEntry[num].m;
		nk = m_ACAMatrix.m_ACAEntry[num].n;

		mbegin = m_ACAMatrix.m_ACAEntry[num].mbeginID;
		nbegin = m_ACAMatrix.m_ACAEntry[num].nbeginID;

		Sum_Flag = 0;  // 用于判断是否是零矩阵

		// 预先计算矩阵初始化
		m_ACAMatrix.m_ACAEntry[num].Entry.initial(mk * 3, nk * 3, 0);

		if (flag == 0)
		{
			for (i = 0; i < mk; ++i) // mk个源点单元
			{
				it = ACAMatrix::RePID[mbegin + i]; // 源点单元实际编号
				i24 = 3 * i;
				for (j = 0; j < 1; ++j)
				{
					ik = 1 * it + j; // 第ik个源点编号
					j3 = 3 * j;

					for (k = 0; k < nk; ++k) // nk的场点单元
					{
						jt = ACAMatrix::RePID[nbegin + k]; // 场点单元使劲编号
						k24 = 3 * k;
						Cal_Flag = UnknownSubMatrix(m_DSE, ik, jt);

						Sum_Flag += Cal_Flag;

						if (Cal_Flag)
						{
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[0], i24 + j3 + 1, k24 + 1, 3, 3);
						/*	m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[1], i24 + j3 + 1, k24 + 4, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[2], i24 + j3 + 1, k24 + 7, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[3], i24 + j3 + 1, k24 + 10, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[4], i24 + j3 + 1, k24 + 13, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[5], i24 + j3 + 1, k24 + 16, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[6], i24 + j3 + 1, k24 + 19, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[7], i24 + j3 + 1, k24 + 22, 3, 3);*/
						}
					}
				}
			}
		}
		else
		{
			for (i = 0; i < mk; ++i) // mk个源点单元
			{
				it = ACAMatrix::RePID[mbegin + i]; // 源点单元编号
				i24 = 3 * i;
				for (j = 0; j < 1; ++j)
				{
					ik = 1 * it + j; // 第ik个源点
					j3 = 3 * j;

					for (k = 0; k < nk; ++k) // nk的场点单元
					{
						jt = ACAMatrix::RePID[nbegin + k]; // 场点单元编号
						k24 = 3 * k;
						Cal_Flag = knownSubMatrix(m_DSE, ik, jt);

						Sum_Flag += Cal_Flag;

						if (Cal_Flag)
						{
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[0], i24 + j3 + 1, k24 + 1, 3, 3);
							/*m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[1], i24 + j3 + 1, k24 + 4, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[2], i24 + j3 + 1, k24 + 7, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[3], i24 + j3 + 1, k24 + 10, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[4], i24 + j3 + 1, k24 + 13, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[5], i24 + j3 + 1, k24 + 16, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[6], i24 + j3 + 1, k24 + 19, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistCoef[7], i24 + j3 + 1, k24 + 22, 3, 3);*/
						}
					}
				}
			}
		}

		if (Sum_Flag == 0) // 零矩阵
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 0;
			m_ACAMatrix.m_ACAEntry[num].Entry.finish();
		}
		else
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 1;

			if (m_ACAMatrix.m_ACAEntry[num].flag == 1)
				m_ACAMatrix.m_ACAEntry[num].FullACA();
		}

	} // end cycle for all ACAEntrys
	char name_temp1[20] = "UnknownMatrix";
	char name_temp2[20] = "KnownMatrix";
	if (ACAMatrix::Flag_ReadWrite == 1) // 首先读入矩阵
	{
		if (flag == 0)
			m_ACAMatrix.WriteToFile(name_temp1, 0);
		else
			m_ACAMatrix.WriteToFile(name_temp2, 0);
	}

}

void CalACAEntryUij(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, double n, double dt)
{
	// 计算ACA矩阵列表，选择Fully ACA
	// 前提：1. ACA矩阵列表初始化完毕
	//       2. 如果ACA矩阵直接计算，则直接计算矩阵已经分配空间
	//          如果ACA矩阵使用ACA分解，则其中的矩阵和ACA向量均没有分配空间
	//       3. RePID数组初始化和计算完毕

	long num, i, i24, j, j3, k, k24, it, jt, mbegin, nbegin, ik, jk, mk, nk;

	Wmatrix sub(3, 3);

	int Cal_Flag, Sum_Flag;

	for (num = 0; num < m_ACAMatrix.EntryCount; ++num)
	{
		mk = m_ACAMatrix.m_ACAEntry[num].m;
		nk = m_ACAMatrix.m_ACAEntry[num].n;

		mbegin = m_ACAMatrix.m_ACAEntry[num].mbeginID;
		nbegin = m_ACAMatrix.m_ACAEntry[num].nbeginID;

		Sum_Flag = 0;

		// 矩阵初始化置零
		m_ACAMatrix.m_ACAEntry[num].Entry.initial(mk * 3, nk * 3, 0);//(mk * 24, nk * 24, 0);

		// 计算矩阵

		for (i = 0; i < mk; ++i) // mk个源点单元
		{
			it = ACAMatrix::RePID[mbegin + i]; // 源点单元实际编号
			i24 = 3 * i;//24 * i;
			for (j = 0; j < 1; ++j)
			{
				ik = 1 * it + j; // 第ik个源点编号
				j3 = 3 * j;

				for (k = 0; k < nk; ++k) // nk的场点单元
				{
					jt = ACAMatrix::RePID[nbegin + k]; // 场点单元使劲编号
					k24 = 3 * k;
					Cal_Flag = m_DSE[jt].IntDynaUijJudge(m_DSE[jt].m_nodelist[ik], n, dt);

					Sum_Flag += Cal_Flag;

					if (Cal_Flag)
					{
						m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[0], i24 + j3 + 1, k24 + 1, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[1], i24 + j3 + 1, k24 + 4, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[2], i24 + j3 + 1, k24 + 7, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[3], i24 + j3 + 1, k24 + 10, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[4], i24 + j3 + 1, k24 + 13, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[5], i24 + j3 + 1, k24 + 16, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[6], i24 + j3 + 1, k24 + 19, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistUij[7], i24 + j3 + 1, k24 + 22, 3, 3);
					}
				}
			}
		}

		// ACA矩阵处理

		if (Sum_Flag == 0)
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 0;
			m_ACAMatrix.m_ACAEntry[num].Entry.finish();
		}
		else
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 1;
			//m_ACAMatrix.m_ACAEntry[num].flag = 0;
			if (m_ACAMatrix.m_ACAEntry[num].flag) // ACA矩阵
			{
				m_ACAMatrix.m_ACAEntry[num].FullACA();
			} // end ACA
		}
	} // end cycle for all ACAEntrys

	long nn = n;

	if ((n - nn) > 0.5)
		nn++;
	else if ((n - nn) < -0.5)
		nn--;
	char name_temp3[5] = "G";
	ACAMatrix::Flag_ReadWrite = 1;
	if (ACAMatrix::Flag_ReadWrite == 1) // 写入矩阵
		m_ACAMatrix.WriteToFile(name_temp3, nn);

}

void CalACAEntryT1ij(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, double n, double dt)
{
	// 计算ACA矩阵列表，选择Fully ACA
	// 前提：1. ACA矩阵列表初始化完毕
	//       2. 如果ACA矩阵直接计算，则直接计算矩阵已经分配空间
	//          如果ACA矩阵使用ACA分解，则其中的矩阵和ACA向量均没有分配空间
	//       3. RePID数组初始化和计算完毕

	long num, i, i24, j, j3, k, k24, it, jt, mbegin, nbegin, ik, jk, mk, nk;
	Wmatrix sub(3, 3);

	int Cal_Flag, Sum_Flag;

	for (num = 0; num < m_ACAMatrix.EntryCount; ++num)
	{
		mk = m_ACAMatrix.m_ACAEntry[num].m;
		nk = m_ACAMatrix.m_ACAEntry[num].n;

		mbegin = m_ACAMatrix.m_ACAEntry[num].mbeginID;
		nbegin = m_ACAMatrix.m_ACAEntry[num].nbeginID;

		Sum_Flag = 0;

		// 矩阵初始化置零
		m_ACAMatrix.m_ACAEntry[num].Entry.initial(mk * 3, nk * 3, 0);//(mk * 24, nk * 24, 0);

		// 计算矩阵

		for (i = 0; i < mk; ++i) // mk个源点单元
		{
			it = ACAMatrix::RePID[mbegin + i]; // 源点单元实际编号
			i24 = 3 * i;//i24 = 24 * i;
			for (j = 0; j < 1; ++j)
			{
				ik = 1* it + j;//8 * it + j; // 第ik个源点编号
				j3 = 3 * j;

				for (k = 0; k < nk; ++k) // nk的场点单元
				{
					jt = ACAMatrix::RePID[nbegin + k]; // 场点单元使劲编号
					k24 = 3 * k; //24 * k;
					Cal_Flag = m_DSE[jt].IntDynaTijJudge(1, m_DSE[jt].m_nodelist[ik], n, dt);

					Sum_Flag += Cal_Flag;

					if (Cal_Flag)
					{
						m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[0], i24 + j3 + 1, k24 + 1, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[1], i24 + j3 + 1, k24 + 4, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[2], i24 + j3 + 1, k24 + 7, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[3], i24 + j3 + 1, k24 + 10, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[4], i24 + j3 + 1, k24 + 13, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[5], i24 + j3 + 1, k24 + 16, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[6], i24 + j3 + 1, k24 + 19, 3, 3);
						//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[7], i24 + j3 + 1, k24 + 22, 3, 3);
					}
				}
			}
		}

		// ACA矩阵处理

		if (Sum_Flag == 0)
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 0;
			m_ACAMatrix.m_ACAEntry[num].Entry.finish();
		}
		else
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 1;

			if (m_ACAMatrix.m_ACAEntry[num].flag) // ACA矩阵
			{
				m_ACAMatrix.m_ACAEntry[num].FullACA();
			} // end ACA
		}
	} // end cycle for all ACAEntrys

	long nn = n;

	if ((n - nn) > 0.5)
		nn++;
	else if ((n - nn) < -0.5)
		nn--;
	char name_temp4[5] = "T1";
	if (ACAMatrix::Flag_ReadWrite == 1) // 写入矩阵
		m_ACAMatrix.WriteToFile(name_temp4, nn);
}

void CalACAEntryT2ij(DSquareElement* m_DSE, ACAMatrix& m_ACAMatrix, double n, double dt)
{
	// 计算ACA矩阵列表，选择Fully ACA
	// 前提：1. ACA矩阵列表初始化完毕
	//       2. 如果ACA矩阵直接计算，则直接计算矩阵已经分配空间
	//          如果ACA矩阵使用ACA分解，则其中的矩阵和ACA向量均没有分配空间
	//       3. RePID数组初始化和计算完毕

	long num, i, i24, j, j3, k, k24, it, jt, mbegin, nbegin, ik, jk, mk, nk, l;

	int ipos;

	Wmatrix sub(3, 3);

	int Cal_Flag, Sum_Flag;

	for (num = 0; num < m_ACAMatrix.EntryCount; ++num)
	{
		mk = m_ACAMatrix.m_ACAEntry[num].m;
		nk = m_ACAMatrix.m_ACAEntry[num].n;

		mbegin = m_ACAMatrix.m_ACAEntry[num].mbeginID;
		nbegin = m_ACAMatrix.m_ACAEntry[num].nbeginID;

		Sum_Flag = 0;

		// 矩阵初始化置零
		m_ACAMatrix.m_ACAEntry[num].Entry.initial(mk * 3, nk *3, 0);//(mk * 24, nk * 24, 0);

		// 计算矩阵

		for (i = 0; i < mk; ++i) // mk个源点单元
		{
			it = ACAMatrix::RePID[mbegin + i]; // 源点单元实际编号
			i24 = 3 * i;//24 * i;
			for (j = 0; j < 1; ++j)//j < 8
			{
				ik = 1 * it + j; //8 * it + j; // 第ik个源点编号
				j3 = 3 * j;

				for (k = 0; k < nk; ++k) // nk的场点单元
				{
					jt = ACAMatrix::RePID[nbegin + k]; // 场点单元使劲编号
					k24 = 3 * k;//24 * k;

					if (n > 1.01) // 无奇异积分
					{
						Cal_Flag = m_DSE[jt].IntDynaTijJudge(2, m_DSE[jt].m_nodelist[ik], n, dt);

						Sum_Flag += Cal_Flag;

						if (Cal_Flag)
						{
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[0], i24 + j3 + 1, k24 + 1, 3, 3);
						/*	m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[1], i24 + j3 + 1, k24 + 4, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[2], i24 + j3 + 1, k24 + 7, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[3], i24 + j3 + 1, k24 + 10, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[4], i24 + j3 + 1, k24 + 13, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[5], i24 + j3 + 1, k24 + 16, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[6], i24 + j3 + 1, k24 + 19, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[7], i24 + j3 + 1, k24 + 22, 3, 3);*/
						}
					}
					else // 奇异积分
					{
						if (m_DSE[jt].IsIn(ik, ipos))
						{
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][0], i24 + j3 + 1, k24 + 1, 3, 3);
	/*						m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][1], i24 + j3 + 1, k24 + 4, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][2], i24 + j3 + 1, k24 + 7, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][3], i24 + j3 + 1, k24 + 10, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][4], i24 + j3 + 1, k24 + 13, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][5], i24 + j3 + 1, k24 + 16, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][6], i24 + j3 + 1, k24 + 19, 3, 3);
							m_ACAMatrix.m_ACAEntry[num].Entry.assign(m_DSE[jt].m_OnElementT2ij[ipos][7], i24 + j3 + 1, k24 + 22, 3, 3);*/

							Sum_Flag++;
						}
						else
						{
							Cal_Flag = m_DSE[jt].IntDynaTijJudge(2, m_DSE[jt].m_nodelist[ik], n, dt);

							Sum_Flag += Cal_Flag;

							if (Cal_Flag)
							{
								m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[0], i24 + j3 + 1, k24 + 1, 3, 3);
								//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[1], i24 + j3 + 1, k24 + 4, 3, 3);
								//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[2], i24 + j3 + 1, k24 + 7, 3, 3);
								//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[3], i24 + j3 + 1, k24 + 10, 3, 3);
								//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[4], i24 + j3 + 1, k24 + 13, 3, 3);
								//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[5], i24 + j3 + 1, k24 + 16, 3, 3);
								//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[6], i24 + j3 + 1, k24 + 19, 3, 3);
								//m_ACAMatrix.m_ACAEntry[num].Entry.assign(DSquareElement::m_AssistTij[7], i24 + j3 + 1, k24 + 22, 3, 3);
							}
						}
					}
				}
			}
		}

		// ACA矩阵处理

		if (Sum_Flag == 0)
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 0;
			m_ACAMatrix.m_ACAEntry[num].Entry.finish();
		}
		else
		{
			m_ACAMatrix.m_ACAEntry[num].flag_zero = 1;

			if (m_ACAMatrix.m_ACAEntry[num].flag) // ACA矩阵
			{
				m_ACAMatrix.m_ACAEntry[num].FullACA();
			} // end ACA
		}
	} // end cycle for all ACAEntrys

	long nn = n;

	if ((n - nn) > 0.5)
		nn++;
	else if ((n - nn) < -0.5)
		nn--;
	char name_temp5[5] = "T2";
	if (ACAMatrix::Flag_ReadWrite == 1) // 写入矩阵
		m_ACAMatrix.WriteToFile(name_temp5, nn);
}


void ReleaseUV(ACAMatrix& m_ACAMatrix)
{

	long i;

	for (i = 0; i < m_ACAMatrix.EntryCount; ++i)
	{
		if (m_ACAMatrix.m_ACAEntry[i].u)
		{
			delete[] m_ACAMatrix.m_ACAEntry[i].u;
			m_ACAMatrix.m_ACAEntry[i].u = 0;
		}
		if (m_ACAMatrix.m_ACAEntry[i].v)
		{
			delete[] m_ACAMatrix.m_ACAEntry[i].v;
			m_ACAMatrix.m_ACAEntry[i].v = 0;
		}

		m_ACAMatrix.m_ACAEntry[i].Entry.finish();
	}
}

//计算AX
int ACA_BuildAX(ACAMatrix& m_ACAMatrix, Wvector& X, Wvector& AX)
{

	long i, j, k, m, n, ik, jk, mbegin, nbegin, mmax, nmax;

	Wvector TX, TAX;
	double sub;

	// 判断ACA列表矩阵的最大值
	mmax = ACAEntry::MAX_M;
	nmax = ACAEntry::MAX_N;

	// 初始化临时向量
	TX.initial(nmax * 3);
	TAX.initial(mmax * 3);

	// 初始化AX
	AX = 0.0;

	for (k = 0; k < m_ACAMatrix.EntryCount; ++k)
	{
		m = m_ACAMatrix.m_ACAEntry[k].m;
		n = m_ACAMatrix.m_ACAEntry[k].n;
		mbegin = m_ACAMatrix.m_ACAEntry[k].mbeginID;
		nbegin = m_ACAMatrix.m_ACAEntry[k].nbeginID;

		// 计算临时向量TX
		for (i = 0; i < n; ++i)
		{
			ik = 3 * ACAMatrix::RePID[nbegin + i];
			jk = 3 * i;
			for (j = 0; j < 3; ++j)
				TX.b[jk + j] = X.b[ik + j];
		}

		if (m_ACAMatrix.m_ACAEntry[k].flag_zero == 0) // 零矩阵
		{
			;
		}
		else if (m_ACAMatrix.m_ACAEntry[k].flag == 0 || (m_ACAMatrix.m_ACAEntry[k].flag == 1 && m_ACAMatrix.m_ACAEntry[k].flag_overflow == 0)) // 直接计算
		{
			mat_vec_product(m_ACAMatrix.m_ACAEntry[k].Entry, TX, TAX);
		}
		else if (m_ACAMatrix.m_ACAEntry[k].flag == 1 && m_ACAMatrix.m_ACAEntry[k].flag_overflow == 1) // ACA计算
		{
			TAX = 0.0;
			for (i = 0; i < m_ACAMatrix.m_ACAEntry[k].num; ++i)
			{
				sub = d_dot(m_ACAMatrix.m_ACAEntry[k].v[i], TX, m_ACAMatrix.m_ACAEntry[k].v[i].n);
				for (j = 0; j < m_ACAMatrix.m_ACAEntry[k].u[i].n; ++j)
					TAX.b[j] += sub * m_ACAMatrix.m_ACAEntry[k].u[i].b[j];
			}
		}

		// 利用TAX给AX赋值
		if (m_ACAMatrix.m_ACAEntry[k].flag_zero == 1) // 非零矩阵
		{
			for (i = 0; i < m; ++i)
			{
				ik = 3 * ACAMatrix::RePID[mbegin + i];
				jk = 3 * i;
				for (j = 0; j < 3; ++j)
					AX.b[ik + j] += TAX.b[jk + j];
			}
		}
	}

	// 释放临时向量
	TX.finish();
	TAX.finish();

	return 1;
}

//计算AX
int ACA_BuildAX_ReadWrite(ACAMatrix& m_ACAMatrix, Wvector& X, Wvector& AX)
{

	long i, j, k, m, n, ik, jk, mbegin, nbegin, mmax, nmax;

	Wvector TX, TAX;
	double sub;

	// 读入矩阵
	if (ACAMatrix::Flag_ReadWrite == 1)
		m_ACAMatrix.ReadFromFile();

	// 判断ACA列表矩阵的最大值
	mmax = ACAEntry::MAX_M;
	nmax = ACAEntry::MAX_N;

	// 初始化临时向量
	TX.initial(nmax * 3);//(nmax * 24);
	TAX.initial(mmax * 3);//(nmax * 24);

	// 初始化AX
	AX = 0.0;

	for (k = 0; k < m_ACAMatrix.EntryCount; ++k)
	{
		m = m_ACAMatrix.m_ACAEntry[k].m;
		n = m_ACAMatrix.m_ACAEntry[k].n;
		mbegin = m_ACAMatrix.m_ACAEntry[k].mbeginID;
		nbegin = m_ACAMatrix.m_ACAEntry[k].nbeginID;

		// 计算临时向量TX
		for (i = 0; i < n; ++i)
		{
			ik = 3 * ACAMatrix::RePID[nbegin + i];
			jk = 3 * i;
			for (j = 0; j < 3; ++j)
				TX.b[jk + j] = X.b[ik + j];
		}

		if (m_ACAMatrix.m_ACAEntry[k].flag_zero == 0) // 零矩阵
		{
			;
		}
		else if (m_ACAMatrix.m_ACAEntry[k].flag == 0 || (m_ACAMatrix.m_ACAEntry[k].flag == 1 && m_ACAMatrix.m_ACAEntry[k].flag_overflow == 0)) // 直接计算
		{
			mat_vec_product(m_ACAMatrix.m_ACAEntry[k].Entry, TX, TAX);
		}
		else if (m_ACAMatrix.m_ACAEntry[k].flag == 1 && m_ACAMatrix.m_ACAEntry[k].flag_overflow == 1) // ACA计算
		{
			TAX = 0.0;
			for (i = 0; i < m_ACAMatrix.m_ACAEntry[k].num; ++i)
			{
				sub = d_dot(m_ACAMatrix.m_ACAEntry[k].v[i], TX, m_ACAMatrix.m_ACAEntry[k].v[i].n);
				for (j = 0; j < m_ACAMatrix.m_ACAEntry[k].u[i].n; ++j)
					TAX.b[j] += sub * m_ACAMatrix.m_ACAEntry[k].u[i].b[j];
			}
		}

		// 利用TAX给AX赋值
		if (m_ACAMatrix.m_ACAEntry[k].flag_zero == 1) // 非零矩阵
		{
			for (i = 0; i < m; ++i)
			{
				ik = 3 * ACAMatrix::RePID[mbegin + i];
				jk = 3 * i;
				for (j = 0; j < 3; ++j)
					AX.b[ik + j] += TAX.b[jk + j];
			}
		}
	}

	// 释放临时向量
	TX.finish();
	TAX.finish();

	// 写入矩阵
	if (ACAMatrix::Flag_ReadWrite == 1)
		m_ACAMatrix.WriteToFile();

	return 1;
}

int ACA_GMRES(ACAMatrix& m_ACAMatrix, Wvector& b, Wvector& x0, Wvector& x, ACAPreConditioner& m_Pre, long restart, double tol, long maxit, int& flag, int* iter)
{
	//A:             系数矩阵
	//b:             右手向量
	//restart:       内循环迭代次数
	//tol:           迭代允许误差
	//maxit:         外循环次数
	//M:             预处理矩阵
	//x0:            初始向量
	//x:             结果向量
	//flag:          结果是否正确0，1，3
	//iter:          内外迭代次数    iter[0]:inner      iter[1]:outer
	long n = b.n;
	long i, j, k, l;                       //used for loop
	double tolb;                             //tol*norm(b)
	double n2b;                              //b的2范数
	Wvector r(n);                    //r=b-A*x
	Wvector vh(n);                   //相当于r0
	Wvector u(n);
	Wvector u2(n);
	Wvector q(restart);
	double normr;                            //r的2范数
	int stag;                           //是否中断
	double phibar;
	double rt, c, s, temp;

	Wmatrix V(n, restart + 1);             //Arnoldi 向量
	Wvector h(restart + 1);               //上Hessenberg矩阵的一列 st A*V = V*H
	Wmatrix QT(restart + 1, restart + 1);    //旋转正交矩阵 st QT*H = R
	Wmatrix R(restart, restart);         //旋转后上三角矩阵 st H = Q*R
	Wvector f(restart);                 //y = R\f => x = x0 + V*y
	Wmatrix W(n, restart);               //W = V*inv(R)

	n2b = vector_norm2(b, n);

	if (n2b <= 0.00000000001)
	{
		x = 0.0;
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	x = x0;
	flag = 1;
	tolb = tol * n2b;

	// A * x = r
	ACA_BuildAX(m_ACAMatrix, x, r);

	for (i = 1; i <= n; ++i)
	{
		r[i] = b[i] - r[i];
	}//end for ii

	normr = vector_norm2(r, n);

	if (normr <= tolb)
	{
		flag = 0;
		iter[0] = 0;
		iter[1] = 0;
		return 1;
	}//end if

	stag = 0;

	for (i = 1; i <= maxit; ++i)  //外循环
	{

		QT = 0.0;
		R = 0.0;

		// 预处理
		ACA_MX(m_Pre, r, vh);

		j = 0;
		h[1] = vector_norm2(vh, n);
		for (j = 1; j <= n; ++j)  //V(:,1) = vh / h(1)
			V(j, 1) = vh[j] / h[1];

		QT(1, 1) = 1.0;

		phibar = h[1];

		for (j = 1; j <= restart; ++j)        //内循环
		{
			printf("当前步：外循环：%ld  内循环：%ld  残差: %lf\n", i, j, normr / n2b);

			for (k = 1; k <= n; ++k)
				u[k] = V(k, j);

			// A * u = u2
			ACA_BuildAX(m_ACAMatrix, u, u2);

			// 预处理
			ACA_MX(m_Pre, u2, u);

			for (k = 1; k <= j; ++k)
			{
				h[k] = 0.0;
				for (l = 1; l <= n; ++l)//h(k) = V(:,k)' * u;
					h[k] += V(l, k)*u[l];
				for (l = 1; l <= n; ++l)//u = u - h(k) * V(:,k)
					u[l] -= h[k] * V(l, k);
			}

			h[j + 1] = vector_norm2(u, n);

			for (k = 1; k <= n; ++k)//V(:,j+1) = u / h(j+1);
				V(k, j + 1) = u[k] / h[j + 1];

			for (k = 1; k <= j; ++k)
			{
				R(k, j) = 0.0;
				for (l = 1; l <= j; ++l)
					R(k, j) += QT(k, l)*h[l];
			}//R(1:j,j) = QT(1:j,1:j) * h(1:j)

			rt = R(j, j);

			if (h[j + 1] == 0.0)
			{
				c = 1.0;
				s = 0.0;
			}
			else if (h[j + 1] * h[j + 1] > rt*rt)
			{
				temp = rt / h[j + 1];
				s = 1.0 / sqrt(1.0 + temp * temp);
				c = -temp * s;
			}
			else
			{
				temp = h[j + 1] / rt;
				c = 1.0 / sqrt(1.0 + temp * temp);
				s = -temp * c;
			}

			R(j, j) = c * rt - s * h[j + 1];

			for (k = 1; k <= j; ++k)
				q[k] = QT(j, k);

			for (k = 1; k <= j; ++k)
				QT(j, k) = c * q[k];

			for (k = 1; k <= j; ++k)
				QT(j + 1, k) = s * q[k];

			QT(j, j + 1) = -s;
			QT(j + 1, j + 1) = c;
			f[j] = c * phibar;
			phibar = s * phibar;

			if (j < restart)
			{
				if (f[j] == 0.0)
					stag = 1;
				if (j == 1)
				{
					for (k = 1; k <= n; ++k)
						W(k, j) = V(k, j) / R(j, j);
				}
				else
				{
					for (k = 1; k <= n; ++k)
					{
						W(k, j) = V(k, j) / R(j, j);
						for (l = 1; l <= j - 1; ++l)
							W(k, j) -= W(k, l)*R(l, j) / R(j, j);
					}
				}
				for (k = 1; k <= n; ++k)
					x[k] += f[j] * W(k, j);
			}
			else
			{
				inv_mat(R, j);
				mat_vec_product(R, f, q);
				mat_vec_product(V, q, x);

				x += x0;
			}

			// A * x = vh
			ACA_BuildAX(m_ACAMatrix, x, vh);

			for (k = 1; k <= n; ++k)
				vh[k] = b[k] - vh[k];

			normr = vector_norm2(vh, n);

			if (normr <= tolb)
			{
				// TEST
				printf("stop here!\n");

				flag = 0;
				iter[0] = i;
				iter[1] = j;
				return 1;
			}

			if (stag == 1)
			{
				flag = 3;
				iter[0] = i;
				iter[1] = j;
				return -1;
			}

		}//end for j

		if (flag == 1)
		{
			x0 = x;

			// A * x0 = r
			ACA_BuildAX(m_ACAMatrix, x0, r);

			for (k = 1; k <= n; ++k)
				r[k] = b[k] - r[k];
			iter[0] = i;
			iter[1] = j;
		}

	}//end for i


	return 1;
}



void ACAMatrixCount(Tree* &m_tree, long& sum)
{
	if (m_tree == 0)
		return;

	int i;

	TreeNeighbor* tempt;

	if (m_tree->Flag) //叶子节点
	{
		tempt = m_tree->m_Neibor;
		while (tempt)
		{
			sum++;
			tempt = tempt->next;
		}
	}

	tempt = m_tree->m_Interaction;
	while (tempt)
	{
		sum++;
		tempt = tempt->next;
	}

	tempt = 0;

	for (i = 0; i < CHILDNUMBER; i++)
		ACAMatrixCount(m_tree->m_Children[i], sum);
}

int DynaACASolver(DSquareElement* m_DSE, Point* m_NodeList, Point* m_EleCenterList, BoundaryValue* bd, long NodeNum, long EleNum, long MaxIterations, long MaxLeafPointCount, long MaxLeafPointCountPre, double Error, double ACAErrorValue, double ACAErrorStop, long NStep, double MaxLength, int Flag_ReadWrite)
{
	long i, j;

	ACAMatrix* MT;
	ACAMatrix* MG;

	// 计算存储量
	double* ACAMatrixMem_MG1 = new double[NStep + 1];
	double* ClaMatrixMem_MG1 = new double[NStep + 1];
	double* ACAMatrixMem_MT1 = new double[NStep + 1];
	double* ClaMatrixMem_MT1 = new double[NStep + 1];
	double* ACAMatrixMem_MT2 = new double[NStep + 1];
	double* ClaMatrixMem_MT2 = new double[NStep + 1];

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// 预处理类
	ACAPreConditioner m_Pre;

	int flag;                                      //结果信号
	int iter[2];

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	ACAMatrix::Flag_ReadWrite = Flag_ReadWrite;

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

	//生成预处理矩阵
	ACAGMRESPreConditioner(m_DSE, m_Pre, m_NodeList, NodeNum, MaxLeafPointCountPre);

	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	// ACA矩阵静态初始化
	ACAEntry::StaticInit(EleNum, ACAErrorValue, ACAErrorStop);

	MT = new ACAMatrix[NStep + 1];
	MG = new ACAMatrix[NStep + 1];

	// 初始化第一个ACA矩阵和重新编号数组
	Setup_ACAEntry_RePID(MT[0], m_EleCenterList, EleNum, MaxLeafPointCount);

	// 初始化其余ACA矩阵
	for (i = 1; i < NStep; ++i)
	{
		if (i <= MaxN)
			GetInfoFromACAEntry(MT[0], MT[i]);
	}
	for (i = 0; i <= NStep; ++i)
	{
		if (i <= MaxN)
			GetInfoFromACAEntry(MT[0], MG[i]);
	}

	// 2. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 3. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("ACA Solver: Step 1.\n");

	CalACAEntry(m_DSE, MT[0], 0);
	CalACAEntry(m_DSE, MG[0], 1);

	x0 = 0.0;

	// TEST
	printf("ACA Solver: Step 1. Calculate T2ij1 and G2ij1...\n");

	CalACAEntryT2ij(m_DSE, MT[1], 1.0, dt);

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	// mat_vec_product(MG[0],bd[1].lu,B);
	ACA_BuildAX_ReadWrite(MG[0], bd[1].lu, B);

	printf("ACA Solver: Step 1. Calculate Result...\n");

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中

	if (Flag_ReadWrite)
		MT[0].ReadFromFile();

	ACA_GMRES(MT[0], B, x0, bd[1].lt, m_Pre, MaxIterations, Error, 1, flag, iter);

	if (Flag_ReadWrite)
		MT[0].WriteToFile();

	// TEST
	printf("ACA Solver: Step 1. Convert...\n");

	// 将 n = 1时刻的未知量/已知量转化为位移和面力
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);





	// 4. 计算 n = 2 ~ NStep的结果

	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("ACA Solver: Step %ld.\n", i);

		if (i <= MaxN)
		{
			// 使用MG[i]暂存MT_1(i-1)和MG_1(i-1)，并分别累加到原来存贮MT_2(i-1)和MG_2(i-1)的MT(i-1)和MG(i-1)
			// GetMatrixDynaU1ij(m_DSE,MG[i],NodeNum,EleNum,(i-1)*1.0,dt);
			CalACAEntryUij(m_DSE, MG[i - 1], (i - 1)*1.0, dt);
			MG[i - 1].Memory(ACAMatrixMem_MG1[i - 1], ClaMatrixMem_MG1[i - 1]);

			// GetMatrixDynaT1ij(m_DSE,MG[i],NodeNum,EleNum,(i-1)*1.0,dt);
			CalACAEntryT1ij(m_DSE, MG[i], (i - 1)*1.0, dt);
			MG[i].Memory(ACAMatrixMem_MT1[i - 1], ClaMatrixMem_MT1[i - 1]);
			MT[i - 1] += MG[i];
			MT[i - 1].Memory(ACAMatrixMem_MT2[i - 1], ClaMatrixMem_MT2[i - 1]);
			ReleaseUV(MG[i]);

			// 计算MT_2(i)(如果i < NStep)，存在MT(i)中
			if (i < NStep)
				CalACAEntryT2ij(m_DSE, MT[i], i*1.0, dt);
		}

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		ACA_BuildAX_ReadWrite(MG[0], bd[i].lu, B);



		// j = 1..i-1
		for (j = 1; j <= i - 1; ++j)
		{
			if (j <= MaxN)
			{
				// CB1 = [MG_1(j) + MG_2(j)] * bd[i-j].lt，这里的bd[i-j].lt明确指的是i-j时刻的面力
				// mat_vec_product(MG[j],bd[i-j].lt,CB1);
				ACA_BuildAX_ReadWrite(MG[j], bd[i - j].lt, CB1);
				B += CB1;

				// CB1 = [MT_1(j) + MT_2(j)] * bd[i-j].lu，这里的bd[i-j].lu明确指的是i-j时刻的位移
				// mat_vec_product(MT[j],bd[i-j].lu,CB1);
				ACA_BuildAX_ReadWrite(MT[j], bd[i - j].lu, CB1);
				B -= CB1;
			}
		}

		if (Flag_ReadWrite)
			MT[0].ReadFromFile();

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		ACA_GMRES(MT[0], B, x0, bd[i].lt, m_Pre, MaxIterations, Error, 1, flag, iter);

		if (Flag_ReadWrite)
			MT[0].WriteToFile();

		// 将 n = i时刻的未知量/已知量转化为位移和面力
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}

	//输出存贮量

	FILE* fp;
	fopen_s(&fp, "ACAMemCompare.txt", "w");
	fprintf(fp, "STEP  ACA_MG_to_CLA_MG  ACA_MT1_to_CLA_MT1  ACA_MT12_to_CLA_MT12\n");
	for (i = 1; i < MaxN; ++i)
		fprintf(fp, "%5ld  %15.10e  %15.10e  %15.10e\n", i, ACAMatrixMem_MG1[i] / ClaMatrixMem_MG1[i], ACAMatrixMem_MT1[i] / ClaMatrixMem_MT1[i], ACAMatrixMem_MT2[i] / ClaMatrixMem_MT2[i]);
	fclose(fp);


	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;

	delete[] ACAMatrixMem_MG1;
	ACAMatrixMem_MG1 = 0;
	delete[] ClaMatrixMem_MG1;
	ClaMatrixMem_MG1 = 0;
	delete[] ACAMatrixMem_MT1;
	ACAMatrixMem_MT1 = 0;
	delete[] ClaMatrixMem_MT1;
	ClaMatrixMem_MT1 = 0;
	delete[] ACAMatrixMem_MT2;
	ACAMatrixMem_MT2 = 0;
	delete[] ClaMatrixMem_MT2;
	ClaMatrixMem_MT2 = 0;

	ACAMatrix::StaticRelease();
	ACAEntry::StaticRelease();

	m_Pre.Finish();

	// End Counting...
	m_Counter.EndCount();

	FILE* pf;
	fopen_s(&pf, "ACATime.txt", "w");
	fprintf(pf, "求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);


	return 1;
}

int DynaACASolverNew(DSquareElement* m_DSE, Point* m_NodeList, Point* m_EleCenterList, BoundaryValue* bd, long NodeNum, long EleNum, long MaxIterations, long MaxLeafPointCount, long MaxLeafPointCountPre, double Error, double ACAErrorValue, double ACAErrorStop, long NStep, double MaxLength, int Flag_ReadWrite)
{
	long i, j;

	ACAMatrix* MT;
	ACAMatrix* MG;

	// 计算存储量
	double* ACAMatrixMem_MG1 = new double[NStep + 1];
	double* ClaMatrixMem_MG1 = new double[NStep + 1];
	double* ACAMatrixMem_MT1 = new double[NStep + 1];
	double* ClaMatrixMem_MT1 = new double[NStep + 1];
	double* ACAMatrixMem_MT2 = new double[NStep + 1];
	double* ClaMatrixMem_MT2 = new double[NStep + 1];

	double dt = DSquareElement::m_DMat.Dt;

	Wvector B(3 * NodeNum), CB1(3 * NodeNum), CB2(3 * NodeNum);

	Wvector x0(3 * NodeNum);

	// 预处理类
	ACAPreConditioner m_Pre;

	int flag;                                      //结果信号
	int iter[2];

	Counter m_Counter;

	// Start Counting...
	m_Counter.StartCount();

	ACAMatrix::Flag_ReadWrite = Flag_ReadWrite;

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

	//生成预处理矩阵
	ACAGMRESPreConditioner(m_DSE, m_Pre, m_NodeList, NodeNum, MaxLeafPointCountPre);

	// 1. 初始化矩阵系列
	// MT1: 0 ~ NStep - 1     MT2: 1 ~ NStep - 1, 因为没有MT2(0)，MT2(NStep)也由于初始位移为0没有用上
	// GT1: 0 ~ NStep - 1     MG2: 1 ~ NStep
	// MT(0)代表MT1[0],也代表未知量矩阵； MG(0)代表MG1[0]，也代表已知量矩阵
	// MT(1~NStep-1)代表MT1 + MT2，MG(1~NStep-1)代表MG1 + MG2
	// MG(NStep)代表MG2[NStep]

	// ACA矩阵静态初始化
	ACAEntry::StaticInit(EleNum, ACAErrorValue, ACAErrorStop);

	MT = new ACAMatrix[NStep + 1];
	MG = new ACAMatrix[NStep + 1];

	// 初始化第一个ACA矩阵和重新编号数组
	Setup_ACAEntry_RePID(MT[0], m_EleCenterList, EleNum, MaxLeafPointCount);

	// 初始化其余ACA矩阵
	for (i = 1; i < NStep; ++i)
	{
		if (i <= MaxN)
			GetInfoFromACAEntry(MT[0], MT[i]);
	}
	for (i = 0; i <= NStep; ++i)
	{
		if (i <= MaxN)
			GetInfoFromACAEntry(MT[0], MG[i]);
	}

	// 2. 计算矩阵
	i = 50;
	CalACAEntryUij(m_DSE, MG[i], i * 1.0, dt);
	






	CalACAEntry(m_DSE, MT[0], 0);
	CalACAEntry(m_DSE, MG[0], 1);







	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th ACA Matrices.\n", i);

		CalACAEntryT1ij(m_DSE, MT[i], i*1.0, dt);
		MT[i].Memory(ACAMatrixMem_MT1[i], ClaMatrixMem_MT1[i]);

		CalACAEntryT2ij(m_DSE, MG[MaxN], i*1.0, dt);
		MT[i] += MG[MaxN];
		MT[i].Memory(ACAMatrixMem_MT2[i], ClaMatrixMem_MT2[i]);
		ReleaseUV(MG[MaxN]);

		CalACAEntryUij(m_DSE, MG[i], i*1.0, dt);
		MG[i].Memory(ACAMatrixMem_MG1[i], ClaMatrixMem_MG1[i]);
	}

	// 3. 将0时刻的已知量和未知量调整顺序为局部坐标下的位移和面力，其中位移为0

	bd[0].lt = 0.0;
	bd[0].Convert(m_DSE, EleNum);
	bd[0].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 4. 计算n = 1时刻的系数矩阵并运算得到 n = 1 时刻的结果

	// TEST
	printf("ACA Solver: Step 1.\n");

	x0 = 0.0;

	// B = MG_1(0) * bd[1].lu，这里的MG_1(0)实际上是已知量矩阵，bd[1].lu是1时刻的已知向量
	// mat_vec_product(MG[0],bd[1].lu,B);
	ACA_BuildAX_ReadWrite(MG[0], bd[1].lu, B);

	// 计算 n = 1时刻的未知向量，存在bd[1].lt中

	if (Flag_ReadWrite)
		MT[0].ReadFromFile();

	ACA_GMRES(MT[0], B, x0, bd[1].lt, m_Pre, MaxIterations, Error, 1, flag, iter);

	if (Flag_ReadWrite)
		MT[0].WriteToFile();

	// 将 n = 1时刻的未知量/已知量转化为位移和面力
	bd[1].Convert(m_DSE, EleNum);
	bd[1].AllLocToAbs(NodeNum, DSquareElement::m_transmat);


	// 5. 计算 n = 2 ~ NStep的结果
	m_Counter.EndCount();
	double DifTime1 = m_Counter.TimeDiff();

	m_Counter.StartCount();



	for (i = 2; i <= NStep; ++i)
	{
		// TEST
		printf("ACA Solver: Step %ld.\n", i);

		// B = MG_1(0) * bd[i].lu，这里的MG_1(0)实际上是已知量矩阵，bd[i].lu是i时刻的已知向量
		ACA_BuildAX_ReadWrite(MG[0], bd[i].lu, B);

		// j = 1..i-1
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

				// CB1 = [MG_1(j) + MG_2(j)] * bd[i-j].lt，这里的bd[i-j].lt明确指的是i-j时刻的面力
				// mat_vec_product(MG[j],bd[i-j].lt,CB1);
				ACA_BuildAX_ReadWrite(MG[j], CB1, CB2);
				B += CB2;

				// [T]
				CB1 = 0.0;
				if ((i - j + 1) > 0 && (i - j + 1) < i)
					CB1 += bd[i - j + 1].lu;
				if ((i - j) > 0 && (i - j) < i)
					CB1.plus(bd[i - j].lu, 2.0);
				if ((i - j - 1) > 0 && (i - j - 1) < i)
					CB1 += bd[i - j - 1].lu;
				// CB1 = [MT_1(j) + MT_2(j)] * bd[i-j].lu，这里的bd[i-j].lu明确指的是i-j时刻的位移
				// mat_vec_product(MT[j],bd[i-j].lu,CB1);
				ACA_BuildAX_ReadWrite(MT[j], CB1, CB2);
				B -= CB2;
			}
		}

		if (Flag_ReadWrite)
			MT[0].ReadFromFile();

		// 计算 n = i时刻的未知向量，存在bd[i].lt中
		ACA_GMRES(MT[0], B, x0, bd[i].lt, m_Pre, MaxIterations, Error, 1, flag, iter);

		if (Flag_ReadWrite)
			MT[0].WriteToFile();

		// 将 n = i时刻的未知量/已知量转化为位移和面力
		bd[i].Convert(m_DSE, EleNum);
		bd[i].AllLocToAbs(NodeNum, DSquareElement::m_transmat);
	}

	//输出存贮量

	FILE* fp;
	char* CurrPath = _getcwd(NULL, 0);//返回当前路径
	string StrCurrPath(CurrPath);
	string FileName, Path;
	unsigned found = StrCurrPath.find_last_of("/\\");
	FileName = StrCurrPath.substr(found + 1);//返回当前文件夹名称
	Path = StrCurrPath.substr(0, found + 1);//返回上级路径
	string OutFile = "output\\";
	string TempPath = Path + OutFile;
	string Temp = "ACAMemCompare.txt";
	Path = TempPath + Temp;
	const char* OutPath = Path.c_str();
	fopen_s(&fp, OutPath, "w");


	fprintf(fp, "STEP  ACA_MG_to_CLA_MG  ACA_MT1_to_CLA_MT1  ACA_MT12_to_CLA_MT12\n");
	for (i = 1; i < MaxN; ++i)
		fprintf(fp, "%5ld  %15.10e  %15.10e  %15.10e\n", i, ACAMatrixMem_MG1[i] / ClaMatrixMem_MG1[i], ACAMatrixMem_MT1[i] / ClaMatrixMem_MT1[i], ACAMatrixMem_MT2[i] / ClaMatrixMem_MT2[i]);
	fclose(fp);


	// 释放空间
	showMemoryInfo();

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;

	delete[] ACAMatrixMem_MG1;
	ACAMatrixMem_MG1 = 0;
	delete[] ClaMatrixMem_MG1;
	ClaMatrixMem_MG1 = 0;
	delete[] ACAMatrixMem_MT1;
	ACAMatrixMem_MT1 = 0;
	delete[] ClaMatrixMem_MT1;
	ClaMatrixMem_MT1 = 0;
	delete[] ACAMatrixMem_MT2;
	ACAMatrixMem_MT2 = 0;
	delete[] ClaMatrixMem_MT2;
	ClaMatrixMem_MT2 = 0;

	ACAMatrix::StaticRelease();
	ACAEntry::StaticRelease();

	m_Pre.Finish();

	// End Counting...
	m_Counter.EndCount();


	Temp = "ACATime.txt";
	Path = TempPath + Temp;
	OutPath = Path.c_str();

	FILE* pf;

	fopen_s(&pf, OutPath, "w");
	fprintf(pf, "时间:\n");
	fprintf(pf, "组装时间:%lf\n", DifTime1);
	fprintf(pf, "迭代求解时间:%lf\n", m_Counter.TimeDiff());
	fclose(pf);



	return 1;
}

int DynaACASolverMem(DSquareElement* m_DSE, Point* m_NodeList, Point* m_EleCenterList, BoundaryValue* bd, long NodeNum, long EleNum, long MaxIterations, long MaxLeafPointCount, long MaxLeafPointCountPre, double Error, double ACAErrorValue, double ACAErrorStop, long NStep, double MaxLength, int Flag_ReadWrite)
{
	long i, j;

	ACAMatrix* MT;
	ACAMatrix* MG;

	// 计算存储量
	double* ACAMatrixMem_MG1 = new double[NStep + 1];
	double* ClaMatrixMem_MG1 = new double[NStep + 1];
	double* ACAMatrixMem_MT1 = new double[NStep + 1];
	double* ClaMatrixMem_MT1 = new double[NStep + 1];
	double* ACAMatrixMem_MT2 = new double[NStep + 1];
	double* ClaMatrixMem_MT2 = new double[NStep + 1];

	double dt = DSquareElement::m_DMat.Dt;

	ACAMatrix::Flag_ReadWrite = Flag_ReadWrite;

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

	// ACA矩阵静态初始化
	ACAEntry::StaticInit(EleNum, ACAErrorValue, ACAErrorStop);

	MT = new ACAMatrix[NStep + 1];
	MG = new ACAMatrix[NStep + 1];

	// 初始化第一个ACA矩阵和重新编号数组
	Setup_ACAEntry_RePID(MT[0], m_EleCenterList, EleNum, MaxLeafPointCount);

	// 初始化其余ACA矩阵
	for (i = 1; i < NStep; ++i)
	{
		if (i <= MaxN)
			GetInfoFromACAEntry(MT[0], MT[i]);
	}
	for (i = 0; i <= NStep; ++i)
	{
		if (i <= MaxN)
			GetInfoFromACAEntry(MT[0], MG[i]);
	}

	// 2. 计算矩阵

	for (i = 1; i <= MaxN; ++i)
	{
		printf("Calculate %ld-th ACA Matrices.\n", i);

		CalACAEntryT1ij(m_DSE, MT[i], i*1.0, dt);
		MT[i].Memory(ACAMatrixMem_MT1[i], ClaMatrixMem_MT1[i]);
		ReleaseUV(MT[i]);

		CalACAEntryUij(m_DSE, MG[i], i*1.0, dt);
		MG[i].Memory(ACAMatrixMem_MG1[i], ClaMatrixMem_MG1[i]);
		ReleaseUV(MG[i]);
	}


	//输出存贮量

	FILE* fp;
	fopen_s(&fp, "ACAMemCompare.txt", "w");
	fprintf(fp, "STEP  ACA_MG_to_CLA_MG  ACA_MT1_to_CLA_MT1  ACA_MT12_to_CLA_MT12\n");
	for (i = 1; i < MaxN; ++i)
		fprintf(fp, "%5ld  %15.10e  %15.10e  %15.10e\n", i, ACAMatrixMem_MG1[i] / ClaMatrixMem_MG1[i], ACAMatrixMem_MT1[i] / ClaMatrixMem_MT1[i], ACAMatrixMem_MT2[i] / ClaMatrixMem_MT2[i]);
	fclose(fp);


	// 释放空间

	delete[] MT;
	MT = 0;
	delete[] MG;
	MG = 0;

	delete[] ACAMatrixMem_MG1;
	ACAMatrixMem_MG1 = 0;
	delete[] ClaMatrixMem_MG1;
	ClaMatrixMem_MG1 = 0;
	delete[] ACAMatrixMem_MT1;
	ACAMatrixMem_MT1 = 0;
	delete[] ClaMatrixMem_MT1;
	ClaMatrixMem_MT1 = 0;
	delete[] ACAMatrixMem_MT2;
	ACAMatrixMem_MT2 = 0;
	delete[] ClaMatrixMem_MT2;
	ClaMatrixMem_MT2 = 0;

	ACAMatrix::StaticRelease();
	ACAEntry::StaticRelease();

	return 1;
}