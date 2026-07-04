#ifndef _COUNTER_
#define _COUNTER_
#include <sys/timeb.h>
#include <sys/types.h>
//-----------------------------------------------------------------------
class Counter
{
public:
	static double RegularUnknownTime;
	static double SingularUnknownTime;
	static double NearSingularUnknownTime;
	static double RegularKnownTime;
	static double SingularKnownTime;
	static double NearSingularKnownTime;

	static long RegularUnknownNum;
	static long SingularUnknownNum;
	static long NearSingularUnknownNum;
	static long RegularKnownNum;
	static long SingularKnownNum;
	static long NearSingularKnownNum;

	static double CollectionTime;
	static double FMMTime;
	static double DirectTime;

	timeb BeginTime;
	timeb EndTime;
	Counter() { ; }
	void StartCount() { ftime(&BeginTime); }
	void EndCount() { ftime(&EndTime); }
	double TimeDiff() { return (0.001*(EndTime.millitm - BeginTime.millitm) + EndTime.time - BeginTime.time); }
	static void PrintInfo();
};
//-----------------------------------------------------------------------
#endif