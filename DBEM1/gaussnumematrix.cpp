#include "GaussNumeMatrix.h"
#include "stdio.h"

GaussNumeMatrix::GaussNumeMatrix()
{
	FILE* gaussp;
	fopen_s(&gaussp, "gaussposition", "r");
	FILE* gaussw;
	fopen_s(&gaussw, "gaussweight", "r");
	long i, j, m;
	double dg[1600];
	for (i = 0; i < 1600; ++i)
	{
		fscanf_s(gaussp, "%lf,", &dg[i]);
	}

	for (i = 0; i < GAUSSPOINT; ++i)
		pos[i] = dg[i + 40 * (GAUSSPOINT - 1)];
	for (i = 0; i < SINGAUSSPOINT; ++i)
		spos[i] = dg[i + 40 * (SINGAUSSPOINT - 1)];
	for (i = 0; i < GAUSSPOINTPW; ++i)
		pwpos[i] = dg[i + 40 * (GAUSSPOINTPW - 1)];
	for (i = 0; i < SINGAUSSPOINTPW; ++i)
		pwspos[i] = dg[i + 40 * (SINGAUSSPOINTPW - 1)];

	for (i = 0; i < 1600; ++i)
	{
		fscanf_s(gaussw, "%lf,", &dg[i]);
	}

	for (i = 0; i < GAUSSPOINT; ++i)
		b[i] = dg[i + 40 * (GAUSSPOINT - 1)];
	for (i = 0; i < SINGAUSSPOINT; ++i)
		sb[i] = dg[i + 40 * (SINGAUSSPOINT - 1)];
	for (i = 0; i < GAUSSPOINTPW; ++i)
		pwb[i] = dg[i + 40 * (GAUSSPOINTPW - 1)];
	for (i = 0; i < SINGAUSSPOINTPW; ++i)
		pwsb[i] = dg[i + 40 * (SINGAUSSPOINTPW - 1)];

	for (i = 0; i < GAUSSPOINT; ++i)
	{
		for (j = 0; j < GAUSSPOINT; ++j)
		{
			m = i * GAUSSPOINT + j;
			rp[m][0] = pos[j];
			rp[m][1] = pos[i];
			ra[m] = b[j] * b[i];
		}
	}

	for (i = 0; i < SINGAUSSPOINT; ++i)
	{
		for (j = 0; j < SINGAUSSPOINT; ++j)
		{
			m = i * SINGAUSSPOINT + j;
			sp[m][0] = spos[j];
			sp[m][1] = spos[i];
			sa[m] = sb[j] * sb[i];
		}
	}

	for (i = 0; i < GAUSSPOINTPW; ++i)
	{
		for (j = 0; j < GAUSSPOINTPW; ++j)
		{
			m = i * GAUSSPOINTPW + j;
			pwrp[m][0] = pwpos[j];
			pwrp[m][1] = pwpos[i];
			pwra[m] = pwb[j] * pwb[i];
		}
	}

	for (i = 0; i < SINGAUSSPOINTPW; ++i)
	{
		for (j = 0; j < SINGAUSSPOINTPW; ++j)
		{
			m = i * SINGAUSSPOINTPW + j;
			pwsp[m][0] = pwspos[j];
			pwsp[m][1] = pwspos[i];
			pwsa[m] = pwsb[j] * pwsb[i];
		}
	}

	fclose(gaussp);
	fclose(gaussw);
}