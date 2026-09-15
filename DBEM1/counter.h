#ifndef _COUNTER_
#define _COUNTER_
#include <chrono>
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

	std::chrono::steady_clock::time_point BeginTime;
	std::chrono::steady_clock::time_point EndTime;
	Counter() { ; }
	void StartCount() { BeginTime = std::chrono::steady_clock::now(); }
	void EndCount() { EndTime = std::chrono::steady_clock::now(); }
	double TimeDiff() { return std::chrono::duration<double>(EndTime - BeginTime).count(); }
	static void PrintInfo();
};
//-----------------------------------------------------------------------
#endif