#include "Counter.h"
#include "stdio.h"
//---------------------------------------------------------------
double Counter::RegularUnknownTime = 0.0;
double Counter::SingularUnknownTime = 0.0;
double Counter::NearSingularUnknownTime = 0.0;
double Counter::RegularKnownTime = 0.0;
double Counter::SingularKnownTime = 0.0;
double Counter::NearSingularKnownTime = 0.0;

long Counter::RegularUnknownNum = 0;
long Counter::SingularUnknownNum = 0;
long Counter::NearSingularUnknownNum = 0;
long Counter::RegularKnownNum = 0;
long Counter::SingularKnownNum = 0;
long Counter::NearSingularKnownNum = 0;

double Counter::CollectionTime = 0.0;
double Counter::FMMTime = 0.0;
double Counter::DirectTime = 0.0;

void Counter::PrintInfo()
{
	printf("*************Count Result******************\n");
	//	printf("RegularUnknown:      times:%ld   cost:%lf  cost per time:%lf\n",Counter::RegularUnknownNum,Counter::RegularUnknownTime,Counter::RegularUnknownNum==0?0.0:(Counter::RegularUnknownTime/Counter::RegularUnknownNum));
	//	printf("SingularUnknown:     times:%ld   cost:%lf  cost per time:%lf\n",Counter::SingularUnknownNum,Counter::SingularUnknownTime,Counter::SingularUnknownNum==0?0.0:(Counter::SingularUnknownTime/Counter::SingularUnknownNum));
	//	printf("NearSingularUnknown: times:%ld   cost:%lf  cost per time:%lf\n",Counter::NearSingularUnknownNum,Counter::NearSingularUnknownTime,Counter::NearSingularUnknownNum==0?0.0:(Counter::NearSingularUnknownTime/Counter::NearSingularUnknownNum));
	//	printf("RegularTij:      times:%ld   cost:%lf  cost per time:%lf\n",Counter::RegularKnownNum,Counter::RegularKnownTime,Counter::RegularKnownNum==0?0.0:(Counter::RegularKnownTime/Counter::RegularKnownNum));
	//	printf("SingularTij:     times:%ld   cost:%lf  cost per time:%lf\n",Counter::SingularKnownNum,Counter::SingularKnownTime,Counter::SingularKnownNum==0?0.0:(Counter::SingularKnownTime/Counter::SingularKnownNum));
	//	printf("NearSingularTij: times:%ld   cost:%lf  cost per time:%lf\n",Counter::NearSingularKnownNum,Counter::NearSingularKnownTime,Counter::NearSingularKnownNum==0?0.0:(Counter::NearSingularKnownTime/Counter::NearSingularKnownNum));
	printf("Collection cost:%lf\n", Counter::CollectionTime);
	printf("Fmm        cost:%lf\n", Counter::FMMTime);
	printf("Direct     cost:%lf\n", Counter::DirectTime);
	printf("*******************************************\n");
}
//---------------------------------------------------------------