#ifndef _PRECOMPILER_
#define _PRECOMPILER_
#include "Dimension.h"
#include "math.h"
//------------------------------------------------------------
const double pi=3.141592653589793239;
const double ONE_OVER_4PI = 0.0795774715459476679;

const int DIRECTION=(DIMENSION==2?4:6);
const int CHILDNUMBER=(DIMENSION==2?4:8);
const double ERRORMAC = 1.0e-20;

const int GAUSSPOINT = 4;
const int GAUSSPOINT2 = GAUSSPOINT * GAUSSPOINT;
const int SINGAUSSPOINT = 4;
const int SINGAUSSPOINT2 = SINGAUSSPOINT * SINGAUSSPOINT;

// 分片积分定义常数
const int GAUSSPOINTPW = 3;                                    // 用于分片积分
const int GAUSSPOINTPW2 = GAUSSPOINTPW * GAUSSPOINTPW;         // GAUSSPOINTPW^2
const int NUMPW = 12;                                          // 用于分片数量
const int NUMPW2 = NUMPW * NUMPW;                              // NUMPW^2

const int SINGAUSSPOINTPW = 3;                                          // 用于分片积分
const int SINGAUSSPOINTPW2 = SINGAUSSPOINTPW * SINGAUSSPOINTPW;         // GAUSSPOINTPW^2
const int SINNUMPW = 20;                                                // 用于分片数量
const int SINNUMPW2 = SINNUMPW * SINNUMPW;                              // NUMPW^2

//------------------------------------------------------------
#endif
