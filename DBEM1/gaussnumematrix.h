#ifndef _GAUSSNUMEMATRIX_
#define _GAUSSNUMEMATRIX_

#include "Precompiler.h"
//---------------------------------------------------------------
class GaussNumeMatrix
{
public:
	double ra[GAUSSPOINT2];
	double rp[GAUSSPOINT2][2];
	double b[GAUSSPOINT];
	double pos[GAUSSPOINT];
	//singular
	double sa[SINGAUSSPOINT2];
	double sp[SINGAUSSPOINT2][2];
	double sb[SINGAUSSPOINT];
	double spos[SINGAUSSPOINT];

	// 分片积分用
	//far field integral
	double pwra[GAUSSPOINTPW2];
	double pwrp[GAUSSPOINTPW2][2];
	double pwb[GAUSSPOINTPW];
	double pwpos[GAUSSPOINTPW];

	double pwsa[SINGAUSSPOINTPW2];
	double pwsp[SINGAUSSPOINTPW2][2];
	double pwsb[SINGAUSSPOINTPW];
	double pwspos[SINGAUSSPOINTPW];

	GaussNumeMatrix();
};
//---------------------------------------------------------------
#endif