#include "mat_vec.h"
#include "stdio.h"
#include "math.h"
#include"vector"
#include "DBEM.h"
using namespace std;
//------------------------------------------------------------------------------

//Wmatrix
Wmatrix::Wmatrix()
{
	m = 1; n = 1;
	a = 0;
}

Wmatrix::Wmatrix(long mm, long nn)
{
	m = mm;
	n = nn;
	long i;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];
}

Wmatrix::Wmatrix(long mm, long nn, int flag)
{
	long i, j;
	m = mm;
	n = nn;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];

	if (flag == 0)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				a[i][j] = 0.0;
			}
	}
	else if (flag == 1 && m == n)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				if (i == j)
					a[i][j] = 1.0;
				else
					a[i][j] = 0.0;
			}
	}
}

Wmatrix::Wmatrix(const Wmatrix& aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	m = aa.m;
	n = aa.n;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
		{
			a[i][j] = aa.a[i][j];
		}
}

Wmatrix::~Wmatrix()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
}

Wmatrix& Wmatrix::operator=(const Wmatrix& aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa.a[i][j];
	return *this;
}

Wmatrix& Wmatrix::operator=(double* aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa[j * m + i];

	return *this;
}

Wmatrix& Wmatrix::operator=(const double aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa;
	return *this;
}

void Wmatrix::initial(long mm, long nn)
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	m = mm;
	n = nn;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];
}

void Wmatrix::initial(long mm, long nn, int flag)
{
	long i, j;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	m = mm;
	n = nn;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];

	if (flag == 0)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				a[i][j] = 0.0;
			}
	}
	else if (flag == 1 && m == n)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				if (i == j)
					a[i][j] = 1.0;
				else
					a[i][j] = 0.0;
			}
	}
}

void Wmatrix::finish()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
}

void Wmatrix::assign(Wmatrix& sub, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 1; i <= mm; ++i)
		for (j = 1; j <= nn; ++j)
			a[br + i - 2][bc + j - 2] = sub(i, j);
}

void Wmatrix::assign(double* aa, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 0; i < mm; ++i)
		for (j = 0; j < nn; ++j)
			a[br + i - 1][bc + j - 1] = aa[j * mm + i];
}

void Wmatrix::Add(Wmatrix& sub, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 1; i <= mm; ++i)
		for (j = 1; j <= nn; ++j)
			a[br + i - 2][bc + j - 2] += sub(i, j);
}

void Wmatrix::Trans()
{
	//Only avaiable for square matrix
	if (m != n)
		return;
	long i, j;
	double temp;
	for (i = 1; i < m; ++i)
	{
		for (j = i + 1; j <= n; ++j)
		{
			temp = a[i - 1][j - 1];
			a[i - 1][j - 1] = a[j - 1][i - 1];
			a[j - 1][i - 1] = temp;

		}
	}
}

void Wmatrix::Minus()
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = -a[i][j];
}

void Wmatrix::Add(Wmatrix& aa)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += aa.a[i][j];
}

Wmatrix& Wmatrix::operator+=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += v;

	return *this;
}

Wmatrix& Wmatrix::operator-=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] -= v;

	return *this;
}

Wmatrix& Wmatrix::operator+=(Wmatrix& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += v.a[i][j];

	return *this;
}

Wmatrix& Wmatrix::operator-=(Wmatrix& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] -= v.a[i][j];

	return *this;
}

Wmatrix& Wmatrix::operator*=(Wmatrix& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] *= v.a[i][j];

	return *this;
}

Wmatrix& Wmatrix::operator*=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] *= v;

	return *this;
}

Wmatrix& Wmatrix::operator/=(Wmatrix& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] /= v.a[i][j];

	return *this;
}

Wmatrix& Wmatrix::operator/=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] /= v;

	return *this;
}

//Wvector
Wvector::Wvector()
{
	n = 1;
	b = 0;
}
Wvector::Wvector(long nn)
{
	n = nn;
	b = new double[n];
}

Wvector::Wvector(long nn, int flag)
{
	n = nn;
	b = new double[n];
	if (flag == 0)
		for (long i = 0; i < n; ++i)
			b[i] = 0.0;
}

void Wvector::initial(long nn)
{
	n = nn;
	b = new double[n];
}

void Wvector::finish()
{
	if (b)
	{
		delete[] b;
		b = 0;
	}
}

Wvector::Wvector(const Wvector& bb)
{
	n = bb.n;
	b = new double[n];
	long i;
	for (i = 0; i < n; ++i)
		b[i] = bb.b[i];
}

Wvector::~Wvector()
{
	if (b)
	{
		delete[] b;
		b = 0;
	}
}

Wvector& Wvector::operator=(const Wvector& bb)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] = bb.b[i];
	return *this;
}

Wvector& Wvector::operator=(const double bb)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] = bb;
	return *this;
}

void Wvector::assign(Wvector& sub, long br, long n)
{
	for (long i = 0; i < n; ++i)
		b[br + i - 1] = sub[i + 1];
}

Wvector& Wvector::operator+=(double v)
{
	long i;

	//#pragma omp parallel for num_threads(thread_num)

	for (i = 0; i < n; ++i)
		b[i] += v;

	return *this;
}

Wvector& Wvector::operator-=(double v)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] -= v;

	return *this;
}

Wvector& Wvector::operator+=(Wvector& v)
{
	
	//#pragma omp parallel for num_threads(thread_num)
	for (long i = 0; i < n; ++i)
		b[i] += v.b[i];

	return *this;
}

Wvector& Wvector::operator-=(Wvector& v)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] -= v.b[i];

	return *this;
}

int Wvector::plus(Wvector& v, double c)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] += v.b[i] * c;

	return 1;
}

int Wvector::minus(Wvector& v, double c)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] -= v.b[i] * c;

	return 1;
}

Wvector& Wvector::operator*=(Wvector& v)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] *= v.b[i];

	return *this;
}

Wvector& Wvector::operator*=(double v)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] *= v;

	return *this;
}

Wvector& Wvector::operator/=(Wvector& v)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] /= v.b[i];

	return *this;
}

Wvector& Wvector::operator/=(double v)
{
	long i;
	for (i = 0; i < n; ++i)
		b[i] /= v;

	return *this;
}

int Wvector::OutputData(char* Name, long ID)
{
	long i;
	char O_Name[100];
	double V_norm = 0.0;
	sprintf_s(O_Name, "%s_%ld", Name, ID);

	FILE* fp;
	fopen_s(&fp, O_Name, "w");
	for (i = 0; i < n; ++i)
	{
		fprintf(fp, "%5ld %12.8e\n", i + 1, b[i]);
		V_norm += b[i] * b[i];
	}
	fprintf(fp, "\n\n%12.8e\n", sqrt(V_norm));
	fclose(fp);

	return 1;
}

int Wvector::OutputData(char* Name, long ID, double* Value)
{
	long i;
	char O_Name[100];
	double V_norm = 0.0;
	sprintf_s(O_Name, "%s_%ld", Name, ID);

	FILE* fp;
	fopen_s(&fp, O_Name, "w");
	for (i = 0; i < n; ++i)
	{
		fprintf(fp, "%5ld %12.8e\n", i + 1, Value[i]);
		V_norm += Value[i] * Value[i];
	}
	fprintf(fp, "\n\n%12.8e\n", sqrt(V_norm));
	fclose(fp);

	return 1;
}


//WmatrixCSR
WmatrixCSR::WmatrixCSR()
{
	m = 1; n = 1;
	a = 0;
}

WmatrixCSR::WmatrixCSR(long mm, long nn)
{
	m = mm;
	n = nn;
	long i;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];
}

WmatrixCSR::WmatrixCSR(long mm, long nn, int flag)
{
	long i, j;
	m = mm;
	n = nn;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];

	if (flag == 0)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				a[i][j] = 0.0;
			}
	}
	else if (flag == 1 && m == n)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				if (i == j)
					a[i][j] = 1.0;
				else
					a[i][j] = 0.0;
			}
	}
}

WmatrixCSR::WmatrixCSR(const WmatrixCSR& aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	m = aa.m;
	n = aa.n;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
		{
			a[i][j] = aa.a[i][j];
		}
}

WmatrixCSR::~WmatrixCSR()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	if (col)
	{
		for (i = 0; i < m; ++i)
		{
			if (col[i])
			{
				delete[] col[i];
				col[i] = 0;
			}
		}
		delete[] col;
		col = 0;
	}
	if (row)
	{
		delete[] row;
		row = 0;
	}
	m = n = num = 0;
}

WmatrixCSR& WmatrixCSR::operator=(const WmatrixCSR& aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa.a[i][j];
	return *this;
}

WmatrixCSR& WmatrixCSR::operator=(double* aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa[j * m + i];

	return *this;
}

WmatrixCSR& WmatrixCSR::operator=(const double aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa;
	return *this;
}



void WmatrixCSR::initial(long mm, long nn)
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	if (col)
	{
		for (i = 0; i < m; ++i)
		{
			if (col[i])
			{
				delete[] col[i];
				col[i] = 0;
			}
		}
		delete[] col;
		col = 0;
	}
	if (row)
	{
		delete[] row;
		row = 0;
	}
	m = mm;
	n = nn;
	row = new int[m];
	num = 0;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];

	col = new int* [m];
	for (i = 0; i < m; ++i)
		col[i] = new int[n];
}

void WmatrixCSR::sum()
{
	num = 0;
	for (int i = 0; i < n; ++i) {
		num += row[i];
	}
}

void WmatrixCSR::reinitial()
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			col[i][j] = 0;
			a[i][j] = 0;
		}
		row[i] = 0;
	}
	num = 0;

}

void WmatrixCSR::initial(long mm, long nn, int flag)
{
	long i, j;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	m = mm;
	n = nn;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];

	if (flag == 0)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				a[i][j] = 0.0;
			}
	}
	else if (flag == 1 && m == n)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				if (i == j)
					a[i][j] = 1.0;
				else
					a[i][j] = 0.0;
			}
	}
}

void WmatrixCSR::finish()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	//if (col)
	//{
	//	for (i = 0; i < m; ++i)
	//	{
	//		if (col[i])
	//		{
	//			delete[] col[i];
	//			col[i] = 0;
	//		}
	//	}
	//	delete[] col;
	//	col = 0;
	//}
	//if (row)
	//{
	//	delete[] row;
	//	row = 0;
	//}
	for (i = 0; i < m; ++i)
	{
		delete[] col[i];
		col[i] = 0;
	}
	delete[] col;
	col = 0;
	delete[] row;
	row = 0;

}

void WmatrixCSR::assign(WmatrixCSR& sub, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 1; i <= mm; ++i)
		for (j = 1; j <= nn; ++j)
			a[br + i - 2][bc + j - 2] = sub(i, j);
}

void WmatrixCSR::assign(double* aa, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 0; i < mm; ++i) {
		for (j = 0; j < nn; ++j) {
			if (aa[j * mm + i] != 0.0) {
				a[br + i][row[br + i]] = aa[j * mm + i];
				col[br + i][row[br + i]] = bc + j;
				++row[br + i];
			}
		}
	}
}

void WmatrixCSR::Add(WmatrixCSR& sub, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 1; i <= mm; ++i)
		for (j = 1; j <= nn; ++j)
			a[br + i - 2][bc + j - 2] += sub(i, j);
}

void WmatrixCSR::Trans()
{
	//Only avaiable for square matrix
	if (m != n)
		return;
	long i, j;
	double temp;
	for (i = 1; i < m; ++i)
	{
		for (j = i + 1; j <= n; ++j)
		{
			temp = a[i - 1][j - 1];
			a[i - 1][j - 1] = a[j - 1][i - 1];
			a[j - 1][i - 1] = temp;

		}
	}
}

void WmatrixCSR::Minus()
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = -a[i][j];
}

void WmatrixCSR::Add(WmatrixCSR& aa)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += aa.a[i][j];
}

WmatrixCSR& WmatrixCSR::operator+=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += v;

	return *this;
}

WmatrixCSR& WmatrixCSR::operator-=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] -= v;

	return *this;
}

WmatrixCSR& WmatrixCSR::operator+=(WmatrixCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += v.a[i][j];

	return *this;
}

WmatrixCSR& WmatrixCSR::operator-=(WmatrixCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] -= v.a[i][j];

	return *this;
}

WmatrixCSR& WmatrixCSR::operator*=(WmatrixCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] *= v.a[i][j];

	return *this;
}

WmatrixCSR& WmatrixCSR::operator*=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] *= v;

	return *this;
}

WmatrixCSR& WmatrixCSR::operator/=(WmatrixCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] /= v.a[i][j];

	return *this;
}

WmatrixCSR& WmatrixCSR::operator/=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] /= v;

	return *this;
}

//WvectorCSR
void WvectorCSR::initial(long nn)
{
	n = nn;
	value = new double[n];
	col = new int[n];
}

void WvectorCSR::reinitial()
{
	for (int i = 0; i < n; ++i) {
		value[i] = 0;
		col[i] = 0;
	}
}

void WvectorCSR::finish()
{
	if (value)
	{
		delete[] value;
		value = 0;
	}
	if (col)
	{
		delete[] col;
		col = 0;
	}
}



CSRMat::CSRMat()
{
	n = 1;
	row_ptr = new int[1];
	col_index = new int[1];
	value = new double[1];
}

CSRMat::~CSRMat()
{
	if (row_ptr)
	{
		delete[] row_ptr;
		row_ptr = 0;
	}
	if (col_index)
	{
		delete[] col_index;
		col_index = 0;
	}
	if (value)
	{
		delete[] value;
		value = 0;
	}

}

//CSRMat
void CSRMat::initial(long nn, long nnum) {

	if (row_ptr)
	{
		delete[] row_ptr;
		row_ptr = 0;
	}
	if (col_index)
	{
		delete[] col_index;
		col_index = 0;
	}
	if (value)
	{
		delete[] value;
		value = 0;
	}

	n = nn;
	long num = nnum;
	long numrow = n + 1;
	row_ptr = new int[numrow];//矩阵每行第一列在value中对应的位置，最后一个数为value元素总数
	col_index = new int[num];//矩阵列序号  从0开始
	value = new double[num];//矩阵值 
}

void CSRMat::Matdelete()
{
	n = 0;
	//num = 0;
	delete[] row_ptr;
	row_ptr = 0;
	delete[] col_index;
	col_index = 0;
	delete[] value;
	value = 0;
}

//WmatrixCCSR
WmatrixCCSR::WmatrixCCSR()
{
	m = 1; n = 1;
	a = 0;
}

WmatrixCCSR::WmatrixCCSR(long mm, long nn)
{
	m = mm;
	n = nn;
	long i;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];
}

WmatrixCCSR::WmatrixCCSR(long mm, long nn, int flag)
{
	long i, j;
	m = mm;
	n = nn;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];

	if (flag == 0)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				a[i][j] = 0.0;
			}
	}
	else if (flag == 1 && m == n)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				if (i == j)
					a[i][j] = 1.0;
				else
					a[i][j] = 0.0;
			}
	}
}

WmatrixCCSR::WmatrixCCSR(const WmatrixCCSR& aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	m = aa.m;
	n = aa.n;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
		{
			a[i][j] = aa.a[i][j];
		}
}

WmatrixCCSR::~WmatrixCCSR()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	if (col)
	{
		for (i = 0; i < m; ++i)
		{
			if (col[i])
			{
				delete[] col[i];
				col[i] = 0;
			}
		}
		delete[] col;
		col = 0;
	}
	if (row)
	{
		delete[] row;
		row = 0;
	}
	m = 0;
	n = 0;
}

WmatrixCCSR& WmatrixCCSR::operator=(const WmatrixCCSR& aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n * 9; ++j)
			a[i][j] = aa.a[i][j];
	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator=(double* aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa[j * m + i];

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator=(const double aa)
{
	//前提：两个矩阵都初始化，而且维数一样
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = aa;
	return *this;
}



void WmatrixCCSR::initial(long mm, long nn)
{
	m = mm / 3;
	n = nn / 3;
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	if (col)
	{
		for (i = 0; i < m; ++i)
		{
			if (col[i])
			{
				delete[] col[i];
				col[i] = 0;
			}
		}
		delete[] col;
		col = 0;
	}
	if (row)
	{
		delete[] row;
		row = 0;
	}
	row = new int[m];
	num = 0;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n * 9];
	col = new int* [m];
	for (i = 0; i < m; ++i)
		col[i] = new int[n];
}

void WmatrixCCSR::sum()
{
	num = 0;
	for (int i = 0; i < m; ++i) {
		num += row[i];
	}
}

void WmatrixCCSR::reinitial()
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n * 9; ++j) {
			a[i][j] = 0;
		}
		for (int j = 0; j < n; ++j) {
			col[i][j] = 0;
		}
		row[i] = 0;
	}
	num = 0;

}

void WmatrixCCSR::initial(long mm, long nn, int flag)
{
	long i, j;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	m = mm;
	n = nn;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n];

	if (flag == 0)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				a[i][j] = 0.0;
			}
	}
	else if (flag == 1 && m == n)
	{
		for (i = 0; i < m; ++i)
			for (j = 0; j < n; ++j)
			{
				if (i == j)
					a[i][j] = 1.0;
				else
					a[i][j] = 0.0;
			}
	}
}

void WmatrixCCSR::finish()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	//if (col)
	//{
	//	for (i = 0; i < m; ++i)
	//	{
	//		if (col[i])
	//		{
	//			delete[] col[i];
	//			col[i] = 0;
	//		}
	//	}
	//	delete[] col;
	//	col = 0;
	//}
	//if (row)
	//{
	//	delete[] row;
	//	row = 0;
	//}
	for (i = 0; i < m; ++i)
	{
		delete[] col[i];
		col[i] = 0;
	}
	delete[] col;
	col = 0;
	delete[] row;
	row = 0;

}



void WmatrixCCSR::assign(double* aa, long br, long bc, long mm, long nn)
{
	long i, j, k, l;
	long num;
	for (i = 0; i < mm; ++i) {
		for (j = 0; j < nn; ++j) {
			if (aa[j * mm + i] != 0.0) {
				num = 0;
				for (k = 0; k < 3; k++) {
					for (l = 0; l < 3; l++) {
						a[br / 3][row[br / 3] * 9 + num] = aa[l * mm + k];
						num++;
					}
				}
				col[br / 3][row[br / 3]] = bc / 3;
				row[br / 3]++;
				return;

			}
		}
	}
	
	return;
	
}

void WmatrixCCSR::assignSym(double* aa, long br, long bc, long mm, long nn)
{
	long i, j, k, l;
	long num;
	for (i = 0; i < mm; ++i) {
		for (j = 0; j < nn; ++j) {
			if (aa[j * mm + i] != 0.0) {
				num = 0;
				for (k = 0; k < 3; k++) {
					for (l = k; l < 3; l++) {
						a[br / 3][row[br / 3] * 6 + num] = aa[l * mm + k];
						num++;
					}
				}
				col[br / 3][row[br / 3]] = bc / 3;
				row[br / 3]++;
				return;

			}
		}
	}
	return;
}

void WmatrixCCSR::Add(WmatrixCCSR& sub, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 1; i <= mm; ++i)
		for (j = 1; j <= nn; ++j)
			a[br + i - 2][bc + j - 2] += sub(i, j);
}

void WmatrixCCSR::Trans()
{
	//Only avaiable for square matrix
	if (m != n)
		return;
	long i, j;
	double temp;
	for (i = 1; i < m; ++i)
	{
		for (j = i + 1; j <= n; ++j)
		{
			temp = a[i - 1][j - 1];
			a[i - 1][j - 1] = a[j - 1][i - 1];
			a[j - 1][i - 1] = temp;

		}
	}
}

void WmatrixCCSR::Minus()
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] = -a[i][j];
}

void WmatrixCCSR::Add(WmatrixCCSR& aa)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += aa.a[i][j];
}

WmatrixCCSR& WmatrixCCSR::operator+=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += v;

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator-=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] -= v;

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator+=(WmatrixCCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] += v.a[i][j];

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator-=(WmatrixCCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] -= v.a[i][j];

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator*=(WmatrixCCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] *= v.a[i][j];

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator*=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] *= v;

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator/=(WmatrixCCSR& v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] /= v.a[i][j];

	return *this;
}

WmatrixCCSR& WmatrixCCSR::operator/=(double v)
{
	long i, j;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			a[i][j] /= v;

	return *this;
}

WvectorCCSR::WvectorCCSR()
{
	n = 1;
	value = 0;
	col = 0;

}

//WvectorCCSR
void WvectorCCSR::initial(long nn)
{
	n = nn / 3;
	value = new double[n * 9];
	col = new int[n];
}

void WvectorCCSR::reinitial()
{
	for (int i = 0; i < n * 9; ++i) {
		value[i] = 0;
	}
	for (int i = 0; i < n; ++i) {
		col[i] = 0;
	}
}

void WvectorCCSR::finish()
{
	if (value)
	{
		delete[] value;
		value = 0;
	}
	if (col)
	{
		delete[] col;
		col = 0;
	}
}

WvectorCCSR::~WvectorCCSR()
{
	if (value)
	{
		delete[] value;
		value = 0;
	}
	if (col)
	{
		delete[] col;
		col = 0;
	}
	n = 0;
}



CCSRMat::CCSRMat()
{
	n = 1;
	row_ptr = new int[1];
	col_index = new int[1];
	value = new double[1];
	row_ptr[0] = 1;
	col_index[0] = 1;
	value[0] = 1;

}

CCSRMat::~CCSRMat()
{
	if (row_ptr)
	{
		delete[] row_ptr;
		row_ptr = 0;
	}
	if (col_index)
	{
		delete[] col_index;
		col_index = 0;
	}
	if (value)
	{
		delete[] value;
		value = 0;
	}
}

//CCSRMat
void CCSRMat::initial(long nn, long nnum) {
	if (row_ptr)
	{
		delete[] row_ptr;
		row_ptr = 0;
	}
	if (col_index)
	{
		delete[] col_index;
		col_index = 0;
	}
	if (value)
	{
		delete[] value;
		value = 0;
	}
	n = nn;//细胞矩阵总行数为    原矩阵行/3
	long long num = nnum;//细胞矩阵总列号数量
	long numrow = n + 1;
	row_ptr = new int[numrow];//矩阵每细胞行第一列在value中对应的位置，最后一个数为value元素总数
	col_index = new int[num];//矩阵列序号  从0开始
	long long numvalue = num * 9;
	value = new double[numvalue];//矩阵值 默认列号数量 * 9
}

void CCSRMat::Matdelete()
{
	n = 0;
	//num = 0;
	delete[] row_ptr;
	row_ptr = 0;
	delete[] col_index;
	col_index = 0;
	delete[] value;
	value = 0;

}

long CCSRMat::NonZeroBlocks() const
{
	if (!row_ptr || n < 0)
		return 0;
	return row_ptr[n];
}

bool CCSRMat::WriteBinary(FILE* fp) const
{
	if (!fp || !row_ptr || !col_index || !value)
		return false;
	long nnzBlocks = NonZeroBlocks();
	if (fwrite(&n, sizeof(long), 1, fp) != 1)
		return false;
	if (fwrite(&nnzBlocks, sizeof(long), 1, fp) != 1)
		return false;
	if (fwrite(row_ptr, sizeof(int), n + 1, fp) != (size_t)(n + 1))
		return false;
	if (nnzBlocks > 0 && fwrite(col_index, sizeof(int), nnzBlocks, fp) != (size_t)nnzBlocks)
		return false;
	long long valueCount = (long long)nnzBlocks * 9;
	if (valueCount > 0 && fwrite(value, sizeof(double), (size_t)valueCount, fp) != (size_t)valueCount)
		return false;
	return true;
}

bool CCSRMat::ReadBinary(FILE* fp)
{
	if (!fp)
		return false;
	long nn = 0;
	long nnzBlocks = 0;
	if (fread(&nn, sizeof(long), 1, fp) != 1)
		return false;
	if (fread(&nnzBlocks, sizeof(long), 1, fp) != 1)
		return false;
	if (nn < 0 || nnzBlocks < 0)
		return false;
	initial(nn, nnzBlocks);
	if (fread(row_ptr, sizeof(int), n + 1, fp) != (size_t)(n + 1))
		return false;
	if (nnzBlocks > 0 && fread(col_index, sizeof(int), nnzBlocks, fp) != (size_t)nnzBlocks)
		return false;
	long long valueCount = (long long)nnzBlocks * 9;
	if (valueCount > 0 && fread(value, sizeof(double), (size_t)valueCount, fp) != (size_t)valueCount)
		return false;
	return row_ptr[n] == nnzBlocks;
}

//————————SymCCSR————————

//WmatrixSymCCSR
WmatrixSymCCSR::WmatrixSymCCSR()
{
	m = 1; n = 1;
	a = 0;
}

WmatrixSymCCSR::~WmatrixSymCCSR()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	if (col)
	{
		for (i = 0; i < m; ++i)
		{
			if (col[i])
			{
				delete[] col[i];
				col[i] = 0;
			}
		}
		delete[] col;
		col = 0;
	}
	if (row)
	{
		delete[] row;
		row = 0;
	}
	m = n = num = 0;
}



void WmatrixSymCCSR::initial(long mm, long nn)
{
	m = mm / 3;
	n = nn / 3;
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	if (col)
	{
		for (i = 0; i < m; ++i)
		{
			if (col[i])
			{
				delete[] col[i];
				col[i] = 0;
			}
		}
		delete[] col;
		col = 0;
	}
	if (row)
	{
		delete[] row;
		row = 0;
	}
	row = new int[m];
	num = 0;
	a = new double* [m];
	for (i = 0; i < m; ++i)
		a[i] = new double[n * 6];
	col = new int* [m];
	for (i = 0; i < m; ++i)
		col[i] = new int[n];
}

void WmatrixSymCCSR::sum()
{
	num = 0;
	for (int i = 0; i < n; ++i) {
		num += row[i];
	}
}

void WmatrixSymCCSR::reinitial()
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n * 6; ++j) {
			a[i][j] = 0;
		}
		for (int j = 0; j < n; ++j) {
			col[i][j] = 0;
		}
		row[i] = 0;
	}
	num = 0;

}

void WmatrixSymCCSR::finish()
{
	long i;
	if (a)
	{
		for (i = 0; i < m; ++i)
		{
			if (a[i])
			{
				delete[] a[i];
				a[i] = 0;
			}
		}
		delete[] a;
		a = 0;
	}
	//if (col)
	//{
	//	for (i = 0; i < m; ++i)
	//	{
	//		if (col[i])
	//		{
	//			delete[] col[i];
	//			col[i] = 0;
	//		}
	//	}
	//	delete[] col;
	//	col = 0;
	//}
	//if (row)
	//{
	//	delete[] row;
	//	row = 0;
	//}
	for (i = 0; i < m; ++i)
	{
		delete[] col[i];
		col[i] = 0;
	}
	delete[] col;
	col = 0;
	delete[] row;
	row = 0;

}

void WmatrixSymCCSR::assign(WmatrixCCSR& sub, long br, long bc, long mm, long nn)
{
	long i, j;
	for (i = 1; i <= mm; ++i)
		for (j = 1; j <= nn; ++j)
			a[br + i - 2][bc + j - 2] = sub(i, j);
}

void WmatrixSymCCSR::assign(double* aa, long br, long bc, long mm, long nn)
{
	long i, j, k, l;
	long num;
	for (i = 0; i < mm; ++i) {
		for (j = 0; j < nn; ++j) {
			if (aa[j * mm + i] != 0.0) {
				num = 0;
				for (k = 0; k < 3; k++) {
					for (l = k; l < 3; l++) {
						a[br / 3][row[br / 3] * 6 + num] = aa[l * mm + k];
						num++;
					}
				}
				col[br / 3][row[br / 3]] = bc / 3;
				row[br / 3]++;
				return;

			}
		}
	}
	return;
}


//WvectorSymCCSR
WvectorSymCCSR::WvectorSymCCSR()
{
	n = 0;

}

void WvectorSymCCSR::initial(long nn)
{
	n = nn / 3;
	value = new double[n * 6];
	col = new int[n];
}

void WvectorSymCCSR::reinitial()
{
	for (int i = 0; i < n * 6; ++i) {
		value[i] = 0;
	}
	for (int i = 0; i < n; ++i) {
		col[i] = 0;
	}
}

void WvectorSymCCSR::finish()
{
	if (value)
	{
		delete[] value;
		value = 0;
	}
	if (col)
	{
		delete[] col;
		col = 0;
	}
}

WvectorSymCCSR::~WvectorSymCCSR()
{
	n = 0;
	if (value)
	{
		delete[] value;
		value = 0;
	}
	if (col)
	{
		delete[] col;
		col = 0;
	}
}

//CCSRMat

SymCCSRMat::SymCCSRMat()
{
	n = 1;
	row_ptr = new int[1];
	col_index = new int[1];
	value = new double[1];
	row_ptr[0] = 1;
	col_index[0] = 1;
	value[0] = 1;

}

SymCCSRMat::~SymCCSRMat()
{
	n = 0;
	if (row_ptr)
	{
		delete[] row_ptr;
		row_ptr = 0;
	}
	if (col_index)
	{
		delete[] col_index;
		col_index = 0;
	}
	if (value)
	{
		delete[] value;
		value = 0;
	}
}


void SymCCSRMat::initial(long nn, long nnum) {
	if (row_ptr)
	{
		delete[] row_ptr;
		row_ptr = 0;
	}
	if (col_index)
	{
		delete[] col_index;
		col_index = 0;
	}
	if (value)
	{
		delete[] value;
		value = 0;
	}
	n = nn;//细胞矩阵总行数为    原矩阵行/3
	long num = nnum;//细胞矩阵总列号数量
	long numrow = n + 1;
	row_ptr = new int[numrow];//矩阵每细胞行第一列在value中对应的位置，最后一个数为value元素总数
	col_index = new int[num];//矩阵列序号  从0开始
	long long numvalue = num * 6;
	value = new double[numvalue];//矩阵值 默认列号数量 * 9
}

void SymCCSRMat::Matdelete()
{
	n = 0;
	//num = 0;
	delete[] row_ptr;
	row_ptr = 0;
	delete[] col_index;
	col_index = 0;
	delete[] value;
	value = 0;

}

long SymCCSRMat::NonZeroBlocks() const
{
	if (!row_ptr || n < 0)
		return 0;
	return row_ptr[n];
}

bool SymCCSRMat::WriteBinary(FILE* fp) const
{
	if (!fp || !row_ptr || !col_index || !value)
		return false;
	long nnzBlocks = NonZeroBlocks();
	if (fwrite(&n, sizeof(long), 1, fp) != 1)
		return false;
	if (fwrite(&nnzBlocks, sizeof(long), 1, fp) != 1)
		return false;
	if (fwrite(row_ptr, sizeof(int), n + 1, fp) != (size_t)(n + 1))
		return false;
	if (nnzBlocks > 0 && fwrite(col_index, sizeof(int), nnzBlocks, fp) != (size_t)nnzBlocks)
		return false;
	long long valueCount = (long long)nnzBlocks * 6;
	if (valueCount > 0 && fwrite(value, sizeof(double), (size_t)valueCount, fp) != (size_t)valueCount)
		return false;
	return true;
}

bool SymCCSRMat::ReadBinary(FILE* fp)
{
	if (!fp)
		return false;
	long nn = 0;
	long nnzBlocks = 0;
	if (fread(&nn, sizeof(long), 1, fp) != 1)
		return false;
	if (fread(&nnzBlocks, sizeof(long), 1, fp) != 1)
		return false;
	if (nn < 0 || nnzBlocks < 0)
		return false;
	initial(nn, nnzBlocks);
	if (fread(row_ptr, sizeof(int), n + 1, fp) != (size_t)(n + 1))
		return false;
	if (nnzBlocks > 0 && fread(col_index, sizeof(int), nnzBlocks, fp) != (size_t)nnzBlocks)
		return false;
	long long valueCount = (long long)nnzBlocks * 6;
	if (valueCount > 0 && fread(value, sizeof(double), (size_t)valueCount, fp) != (size_t)valueCount)
		return false;
	return row_ptr[n] == nnzBlocks;
}


//------------------------------------------------------------------------------

PthArg::PthArg()
{
	targs = 0;
	start = 0;
	end = 0;
	sum = 0;
	AddStrat = 0;
	AddEnd = 0;
	AddSum = 0;
}

PthArg& PthArg::operator=(PthArg& M)
{
	targs = M.targs;
	start = M.start;
	end = M.end;
	sum = M.sum;
	return *this;
	// TODO: 在此处插入 return 语句
}

PthArg::~PthArg()
{
	targs = 0;
	start = 0;
	end = 0;
	sum = 0;
	AddStrat = 0;
	AddEnd = 0;
	AddSum = 0;
}

PthConvol::PthConvol()
{
	PthStep = 0;
	targs = 0;
	MaxN = 0;

}

void PthConvol::initial(long nn)
{
	long i;
	PthStep = nn;
	MaxN = nn - 1;
	start = new long[PthStep];
	end = new long[PthStep];
	for (i = 0; i < PthStep; i++) {
		start[i] = 0;
		end[i] = 0;
	}
}

PthConvol::~PthConvol()
{
	long targs;
	long PthStep;
	delete[] start;
	delete[] end;
	start = 0;
	end = 0;

}


//————————————————————
void SparseConvert(WmatrixCCSR& A, CCSRMat& M)
{
	int i, j;
	int nnum = 0;
	int colnum = 0;
	//A.sum();
	//M.initial(A.n, A.num);//A.num为细胞矩阵个数
	//M.row_ptr[0] = 0;
	for (i = 0; i < A.n; ++i) {
		for (j = 0; j < A.row[i] * 9; ++j) {
			M.value[nnum] = A.a[i][j];
			++nnum;
		}
		for (j = 0; j < A.row[i]; j++) {
			M.col_index[colnum] = A.col[i][j];
			colnum++;
		}
		M.row_ptr[i + 1] = M.row_ptr[i] + A.row[i];
	}
	A.reinitial();
}

void SparseAddition(WmatrixCCSR& A, CCSRMat& M, WvectorCCSR& b)
{
	long nn = A.n;//细胞矩阵行数   原矩阵/3
	long i, j, k;
	long b_row;//b中非零元素个数
	int A_index, M_index;//计算M=A+M一行中的A和M的指示器
	int M_row_ptr;
	int M_rownum;//M中第i行非零元素个数
	for (i = 0; i < nn; ++i) {
		A_index = 0;
		M_index = 0;
		M_row_ptr = M.row_ptr[i];//第i行第一个非零元素在M.value中的位置；
		b_row = 0;
		M_rownum = M.row_ptr[i + 1] - M.row_ptr[i];
		b.reinitial();//b归0
		while (A_index < A.row[i] && M_index < M_rownum) {
			if (A.col[i][A_index] < M.col_index[M_row_ptr + M_index]) {//A<B
				for (k = 0; k < 9; k++) {
					b.value[b_row * 9 + k] = A.a[i][A_index * 9 + k];
				}

				b.col[b_row] = A.col[i][A_index];
				++b_row;
				++A_index;

			}
			else if (A.col[i][A_index] > M.col_index[M_row_ptr + M_index]) {//A>B
				for (k = 0; k < 9; k++) {
					b.value[b_row * 9 + k] = M.value[(M_row_ptr + M_index) * 9 + k];
				}
				b.col[b_row] = M.col_index[M_row_ptr + M_index];
				++b_row;
				++M_index;
			}
			else {//A=B
				for (k = 0; k < 9; k++) {
					b.value[b_row * 9 + k] = A.a[i][A_index * 9 + k] + M.value[(M_row_ptr + M_index) * 9 + k];
				}
				b.col[b_row] = A.col[i][A_index];
				++b_row;
				++M_index;
				++A_index;
			}
		}
		while (A_index < A.row[i]) {

			for (k = 0; k < 9; k++) {
				b.value[b_row * 9 + k] = A.a[i][A_index * 9 + k];
			}
			b.col[b_row] = A.col[i][A_index];
			++b_row;
			++A_index;
		}
		while (M_index < M_rownum) {
			for (k = 0; k < 9; k++) {
				b.value[b_row * 9 + k] = M.value[(M_row_ptr + M_index) * 9 + k];
			}
			b.col[b_row] = M.col_index[M_row_ptr + M_index];
			++b_row;
			++M_index;
		}
		A.row[i] = b_row;
		for (j = 0; j < b_row; ++j) {
			for (k = 0; k < 9; k++) {
				A.a[i][j * 9 + k] = b.value[j * 9 + k];
			}
			A.col[i][j] = b.col[j];
		}
	}

	M.Matdelete();
	MCCSRtemp.sum();
	MTS[StepNum].initial(MCCSRtemp.m, MCCSRtemp.num);
	MTS[StepNum].row_ptr[0] = 0;
	SparseConvert(A, M);
	A.reinitial();
	b.reinitial();
}

int SparseMul(CCSRMat& M, Wvector& x, Wvector& b)
{
	int i, j, k;
	int n = M.n;
	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

	for (i = 0; i < n; ++i) {
		b[i * 3 + 1] = 0.0;
		b[i * 3 + 2] = 0.0;
		b[i * 3 + 3] = 0.0;
		rownum = M.row_ptr[i + 1] - M.row_ptr[i];//某细胞行非零细胞矩阵个数
		rowstart = M.row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			btemp[0] = x[M.col_index[rowstart + j] * 3 + 1];
			btemp[1] = x[M.col_index[rowstart + j] * 3 + 2];
			btemp[2] = x[M.col_index[rowstart + j] * 3 + 3];
			b[i * 3 + 1] += M.value[rowstart * 9 + j * 9] * btemp[0] + M.value[rowstart * 9 + j * 9 + 1] * btemp[1] + M.value[rowstart * 9 + j * 9 + 2] * btemp[2];
			b[i * 3 + 2] += M.value[rowstart * 9 + j * 9 + 3] * btemp[0] + M.value[rowstart * 9 + j * 9 + 4] * btemp[1] + M.value[rowstart * 9 + j * 9 + 5] * btemp[2];
			b[i * 3 + 3] += M.value[rowstart * 9 + j * 9 + 6] * btemp[0] + M.value[rowstart * 9 + j * 9 + 7] * btemp[1] + M.value[rowstart * 9 + j * 9 + 8] * btemp[2];
		}
	}

	return 0;

}




//————————SymCCSR_Operaate————————————

void SparseConvert(WmatrixSymCCSR& A, SymCCSRMat& M)
{
	int i, j;
	int nnum = 0;
	int colnum = 0;
	//A.sum();
	//M.initial(A.n, A.num);//A.num为细胞矩阵个数
	//M.row_ptr[0] = 0;
	for (i = 0; i < A.n; ++i) {
		for (j = 0; j < A.row[i] * 6; ++j) {
			M.value[nnum] = A.a[i][j];
			++nnum;
		}
		for (j = 0; j < A.row[i]; j++) {
			M.col_index[colnum] = A.col[i][j];
			colnum++;
		}
		M.row_ptr[i + 1] = M.row_ptr[i] + A.row[i];
	}
	A.reinitial();

}

void SparseAddition(WmatrixSymCCSR& A, SymCCSRMat& M, WvectorSymCCSR& b)
{
	long nn = A.n;//细胞矩阵行数   原矩阵/3
	long i, j, k;
	long b_row;//b中非零元素个数
	int A_index, M_index;//计算M=A+M一行中的A和M的指示器
	int M_row_ptr;
	int M_rownum;//M中第i行非零元素个数
	for (i = 0; i < nn; ++i) {
		A_index = 0;
		M_index = 0;
		M_row_ptr = M.row_ptr[i];//第i行第一个非零元素在M.value中的位置；
		b_row = 0;
		M_rownum = M.row_ptr[i + 1] - M.row_ptr[i];
		b.reinitial();//b归0
		while (A_index < A.row[i] && M_index < M_rownum) {
			if (A.col[i][A_index] < M.col_index[M_row_ptr + M_index]) {//A<B
				for (k = 0; k < 6; k++) {
					b.value[b_row * 6 + k] = A.a[i][A_index * 6 + k];
				}

				b.col[b_row] = A.col[i][A_index];
				++b_row;
				++A_index;

			}
			else if (A.col[i][A_index] > M.col_index[M_row_ptr + M_index]) {//A>B
				for (k = 0; k < 6; k++) {
					b.value[b_row * 6 + k] = M.value[(M_row_ptr + M_index) * 6 + k];
				}
				b.col[b_row] = M.col_index[M_row_ptr + M_index];
				++b_row;
				++M_index;
			}
			else {//A=B
				for (k = 0; k < 6; k++) {
					b.value[b_row * 6 + k] = A.a[i][A_index * 6 + k] + M.value[(M_row_ptr + M_index) * 6 + k];
				}
				b.col[b_row] = A.col[i][A_index];
				++b_row;
				++M_index;
				++A_index;
			}
		}
		while (A_index < A.row[i]) {

			for (k = 0; k < 6; k++) {
				b.value[b_row * 6 + k] = A.a[i][A_index * 6 + k];
			}
			b.col[b_row] = A.col[i][A_index];//断点
			++b_row;
			++A_index;
		}
		while (M_index < M_rownum) {
			for (k = 0; k < 6; k++) {
				b.value[b_row * 6 + k] = M.value[(M_row_ptr + M_index) * 6 + k];
			}
			b.col[b_row] = M.col_index[M_row_ptr + M_index];
			++b_row;
			++M_index;
		}
		A.row[i] = b_row;
		for (j = 0; j < b_row; ++j) {
			for (k = 0; k < 6; k++) {
				A.a[i][j * 6 + k] = b.value[j * 6 + k];
			}
			A.col[i][j] = b.col[j];
		}
	}

	M.Matdelete();
	SparseConvert(A, M);
	A.reinitial();
	b.reinitial();

}

int SparseMul(SymCCSRMat& M, Wvector& x, Wvector& b)
{
	int i, j, k;
	int n = M.n;
	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

	for (i = 0; i < n; ++i) {
		b[i * 3 + 1] = 0.0;
		b[i * 3 + 2] = 0.0;
		b[i * 3 + 3] = 0.0;
		rownum = M.row_ptr[i + 1] - M.row_ptr[i];//某细胞行非零细胞矩阵个数
		rowstart = M.row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			btemp[0] = x[M.col_index[rowstart + j] * 3 + 1];
			btemp[1] = x[M.col_index[rowstart + j] * 3 + 2];
			btemp[2] = x[M.col_index[rowstart + j] * 3 + 3];
			b[i * 3 + 1] += M.value[rowstart * 6 + j * 6] * btemp[0] + M.value[rowstart * 6 + j * 6 + 1] * btemp[1] + M.value[rowstart * 6 + j * 6 + 2] * btemp[2];
			b[i * 3 + 2] += M.value[rowstart * 6 + j * 6 + 1] * btemp[0] + M.value[rowstart * 6 + j * 6 + 3] * btemp[1] + M.value[rowstart * 6 + j * 6 + 4] * btemp[2];
			b[i * 3 + 3] += M.value[rowstart * 6 + j * 6 + 2] * btemp[0] + M.value[rowstart * 6 + j * 6 + 4] * btemp[1] + M.value[rowstart * 6 + j * 6 + 5] * btemp[2];
		}
	}

	return 0;

}






void SparseConvert(WmatrixCSR& A, CSRMat& M)
{
	int i, j;
	int nnum = 0;
	A.sum();
	M.initial(A.n, A.num);
	M.row_ptr[0] = 0;
	for (i = 0; i < A.n; ++i) {
		for (j = 0; j < A.row[i]; ++j) {
			M.value[nnum] = A.a[i][j];
			M.col_index[nnum] = A.col[i][j];
			++nnum;
		}
		M.row_ptr[i + 1] = M.row_ptr[i] + A.row[i];
	}
	A.reinitial();
}

void SparseAddition(WmatrixCSR& A, CSRMat& M, WvectorCSR& b)
{
	long nn = A.n;
	long i, j;
	long b_row;//b中非零元素个数
	int A_index, M_index;//计算M=A+M一行中的A和M的指示器
	int M_row_ptr;
	int M_rownum;//M中第i行非零元素个数
	for (i = 0; i < nn; ++i) {
		A_index = 0;
		M_index = 0;
		M_row_ptr = M.row_ptr[i];//第i行第一个非零元素在M.value中的位置；
		b_row = 0;
		M_rownum = M.row_ptr[i + 1] - M.row_ptr[i];
		b.reinitial();//b归0
		while (A_index < A.row[i] && M_index < M_rownum) {
			if (A.col[i][A_index] < M.col_index[M_row_ptr + M_index]) {//A<B
				b.value[b_row] = A.a[i][A_index];
				b.col[b_row] = A.col[i][A_index];
				++b_row;
				++A_index;

			}
			else if (A.col[i][A_index] > M.col_index[M_row_ptr + M_index]) {//A>B
				b.value[b_row] = M.value[M_row_ptr + M_index];
				b.col[b_row] = M.col_index[M_row_ptr + M_index];
				++b_row;
				++M_index;
			}
			else {//A=B
				b.value[b_row] = M.value[M_row_ptr + M_index] + A.a[i][A_index];
				b.col[b_row] = A.col[i][A_index];
				++b_row;
				++M_index;
				++A_index;
			}
		}
		while (A_index < A.row[i]) {
			b.value[b_row] = A.a[i][A_index];
			b.col[b_row] = A.col[i][A_index];//断点
			++b_row;
			++A_index;
		}
		while (M_index < M_rownum) {
			b.value[b_row] = M.value[M_row_ptr + M_index];
			b.col[b_row] = M.col_index[M_row_ptr + M_index];
			++b_row;
			++M_index;
		}
		A.row[i] = b_row;
		for (j = 0; j < b_row; ++j) {
			A.a[i][j] = b.value[j];
			A.col[i][j] = b.col[j];
		}
	}

	M.Matdelete();
	SparseConvert(A, M);
	A.reinitial();
	b.reinitial();
}

int SparseMul(CSRMat& M, Wvector& x, Wvector& b)
{
	int i, j;
	int n = M.n;
	int rownum;
	for (i = 0; i < n; ++i) {
		b[i + 1] = 0.0;
		rownum = M.row_ptr[i + 1] - M.row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			b[i + 1] += M.value[M.row_ptr[i] + j] * x[M.col_index[M.row_ptr[i] + j] + 1];

		}

	}


	return 0;
}


