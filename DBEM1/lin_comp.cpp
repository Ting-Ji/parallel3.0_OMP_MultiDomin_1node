#include "lin_comp.h"
#include "math.h"

//-----------------------------------------------------

int gauss_mat(Wmatrix& A, Wvector& B, long n)
{
	long *js, l, k, i, j, is;
	double d, t;
	js = new long[n + 1];
	l = 1;
	for (k = 1; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n; i++)
			for (j = k; j <= n; j++)
			{
				t = A(i, j)*A(i, j);
				if (t > d) { d = t; js[k] = j; is = i; }
			}
		if (d + 1.0 == 1.0) l = 0;
		else
		{
			if (js[k] != k)
				for (i = 1; i <= n; i++)
				{
					t = A(i, k); A(i, k) = A(i, js[k]); A(i, js[k]) = t;
				}
			if (is != k)
			{
				for (j = k; j <= n; j++)
				{
					t = A(k, j); A(k, j) = A(is, j); A(is, j) = t;
				}
				t = B[k]; B[k] = B[is]; B[is] = t;
			}
		}
		if (l == 0)
		{
			delete[] js;
			return 0;
		}
		d = A(k, k);
		for (j = k + 1; j <= n; j++)
		{
			A(k, j) = A(k, j) / d;
		}
		B[k] = B[k] / d;
		for (i = k + 1; i <= n; i++)
		{
			for (j = k + 1; j <= n; j++)
			{
				A(i, j) = A(i, j) - A(i, k)*A(k, j);
			}
			B[i] = B[i] - A(i, k)*B[k];
		}
	}
	d = A(n, n);
	if (d + 1.0 == 1.0)
	{
		delete[] js;
		return 0;
	}
	B[n] = B[n] / d;
	for (i = n - 1; i >= 1; i--)
	{
		t = 0.0;
		for (j = i + 1; j <= n; j++)
			t = t + A(i, j)*B[j];
		B[i] = B[i] - t;
	}
	js[n] = n;
	for (k = n; k >= 1; k--)
		if (js[k] != k)
		{
			t = B[k]; B[k] = B[js[k]]; B[js[k]] = t;
		}
	delete[] js;

	return 1;
}

int inv_mat(Wmatrix& A, long n)
{
	double eps = 0.0000000001;
	long k, i, j, *is, *js;
	double t, d, temp;
	is = new long[n + 1];
	js = new long[n + 1];

	if (is == 0 || js == 0)
	{
		delete[] is;
		delete[] js;
		return -1;
	}

	for (k = 1; k <= n; k++)
	{
		d = 0.0;
		for (i = k; i <= n; i++)
		{
			for (j = k; j <= n; j++)
			{
				t = sqrt(A(i, j)*A(i, j));
				if (t > d)
				{
					d = t; is[k] = i; js[k] = j;
				}
			}//end for j
		}  //end for i

		if (d <= eps)
		{
			delete[] is;
			delete[] js;
			return -1;
		}

		for (j = 1; j <= n; j++)
		{
			i = is[k];
			temp = A(k, j);
			A(k, j) = A(i, j);
			A(i, j) = temp;
		}

		for (i = 1; i <= n; i++)
		{
			j = js[k];
			temp = A(i, k);
			A(i, k) = A(i, j);
			A(i, j) = temp;
		}

		A(k, k) = 1.0 / A(k, k);

		for (j = 1; j <= n; j++)
		{
			if (j != k)
				A(k, j) *= A(k, k);
		}

		for (i = 1; i <= n; i++)
		{
			if (i != k)
			{
				for (j = 1; j <= n; j++)
				{
					if (j != k)
						A(i, j) -= A(i, k)*A(k, j);
				}
			}
		}//end for i

		for (i = 1; i <= n; i++)
		{
			if (i != k)
				A(i, k) *= -1.0*A(k, k);
		}//end for i

	}//end for k

	for (k = n; k >= 1; k--)
	{
		for (j = 1; j <= n; j++)
		{
			i = js[k];
			temp = A(k, j);
			A(k, j) = A(i, j);
			A(i, j) = temp;
		}//end for j

		for (i = 1; i <= n; i++)
		{
			j = is[k];
			temp = A(i, k);
			A(i, k) = A(i, j);
			A(i, j) = temp;

		}//end for i

	}//end for k

	delete[] is;
	delete[] js;

	return 1;
}

double d_dot(Wvector& a, Wvector& b, long n)
{
	double res = 0.0;
	for (long i = 1; i <= n; ++i)
		res += a[i] * b[i];
	return res;
}

double vector_norm2(Wvector& a, long n)
{
	double res = 0.0;
	for (long i = 1; i <= n; ++i)
		res += a[i] * a[i];
	return sqrt(res);
}

int mat_vec_product(Wmatrix& A, Wvector& x, Wvector& b)
{
	long m = A.m;
	long n = A.n;
	long i, j;
	for (i = 1; i <= m; ++i)
	{
		b[i] = 0.0;
		for (j = 1; j <= n; ++j)
			b[i] += A(i, j)*x[j];
	}

	return 1;
}

int MMMul(Wmatrix& A, Wmatrix& B, Wmatrix& C)
{
	long m = A.m;
	long n = A.n;
	long l = B.n;

	long i, j, k;
	for (i = 1; i <= m; ++i)
	{
		for (j = 1; j <= l; ++j)
		{
			C(i, j) = 0.0;
			for (k = 1; k <= n; ++k)
				C(i, j) += A(i, k)*B(k, j);
		}
	}

	return 1;
}

//取矩阵的最大值和位置
int ArgMax(Wmatrix& A, long& i, long& j, double& value)
{
	long it, jt;
	double tvalue = 0.0;

	it = 1;
	jt = 1;
	value = 0.0;

	for (it = 1; it <= A.m; ++it)
	{
		for (jt = 1; jt <= A.n; ++jt)
		{
			if (tvalue <= fabs(A(it, jt)))
			{
				i = it - 1;
				j = jt - 1;
				tvalue = fabs(A(it, jt));
				value = A(it, jt);
			}
		}
	}

	return 1;
}

//取矩阵的最大值和位置
int ArgMax(Wvector& b, long& i, double& value)
{
	long it;
	double tvalue = 0.0;

	it = 1;
	value = 0.0;

	for (it = 1; it <= b.n; ++it)
	{
		if (tvalue <= fabs(b[it]))
		{
			i = it - 1;
			tvalue = fabs(b[it]);
			value = b[it];
		}
	}

	return 1;
}

//取矩阵的F范数
int normF(Wmatrix& A, double& value)
{
	long i, j;
	value = 0.0;
	for (i = 1; i <= A.m; ++i)
		for (j = 1; j <= A.n; ++j)
			value += A(i, j)*A(i, j);
	value = sqrt(value);

	return 1;
}

//取向量的F范数
int normF(Wvector& b, double& value)
{
	long i;
	value = 0.0;
	for (i = 1; i <= b.n; ++i)
		value += b[i] * b[i];
	value = sqrt(value);

	return 1;
}

//-----------------------------------------------------------
