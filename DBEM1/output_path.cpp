#include "output_path.h"

#include <direct.h>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <map>
#include <share.h>

// Only the serial coordinator writes timing records.
struct DBEMTimingLog {
    FILE* file;
    std::map<std::string, std::pair<long, double> > totals;
    std::map<std::string, double> pendingMerges;
    DBEMTimingLog() : file(0) {}
    ~DBEMTimingLog() { if (file) fclose(file); }
};
static DBEMTimingLog g_timing;

double DBEMWallTime()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

void DBEMWriteTiming(const char* stage, double seconds, long step)
{
    if (!g_timing.file) return;
    std::pair<long, double>& total = g_timing.totals[stage];
    ++total.first;
    total.second += seconds;
    // Batch frequent merge records until the next enclosing stage is reported.
    if (strstr(stage, "multidomain.assembly.merge.") == stage) {
        g_timing.pendingMerges[stage] += seconds;
        return;
    }
    for (const auto& pending : g_timing.pendingMerges) {
        const auto& cumulative = g_timing.totals[pending.first];
        fprintf(g_timing.file, "%s\t%ld\t%.9f\t%ld\t%.9f\n",
            pending.first.c_str(), step, pending.second, cumulative.first, cumulative.second);
    }
    g_timing.pendingMerges.clear();
    fprintf(g_timing.file, "%s\t%ld\t%.9f\t%ld\t%.9f\n", stage, step, seconds, total.first, total.second);
    fflush(g_timing.file);
}

static std::string g_outputDirectory = "output";

static std::string DBEMBaseName(const char* value)
{
	if (!value || !value[0])
		return "case";
	const char* start = value;
	for (const char* p = value; *p; ++p)
	{
		if (*p == '/' || *p == '\\')
			start = p + 1;
	}
	std::string name(start);
	size_t dot = name.find_last_of('.');
	if (dot != std::string::npos)
		name = name.substr(0, dot);
	if (name.empty())
		name = "case";
	return name;
}

static std::string DBEMSanitizeName(const std::string& name)
{
	std::string result;
	for (size_t i = 0; i < name.size(); ++i)
	{
		unsigned char ch = (unsigned char)name[i];
		if (std::isalnum(ch) || ch == '_' || ch == '-')
			result.push_back((char)ch);
		else
			result.push_back('_');
	}
	if (result.empty())
		result = "case";
	return result;
}

void DBEMInitializeOutputDirectory(const char* modelName, long NStep)
{
	char stepText[64];
	sprintf_s(stepText, "%ld", NStep);
	std::string caseName = DBEMSanitizeName(DBEMBaseName(modelName)) + "_" + stepText;
	_mkdir("output");
	g_outputDirectory = std::string("output\\") + caseName;
	_mkdir(g_outputDirectory.c_str());
    if (g_timing.file) fclose(g_timing.file);
    g_timing.file = 0;
    g_timing.totals.clear();
    g_timing.pendingMerges.clear();
    const std::string timingPath = g_outputDirectory + "\\timing.txt";
    g_timing.file = _fsopen(timingPath.c_str(), "w", _SH_DENYNO);
    if (!g_timing.file)
        fprintf(stderr, "Cannot write timing file: %s\n", timingPath.c_str());
    else
    {
        fprintf(g_timing.file, "# Wall-clock seconds; step=-1 means not applicable. Last record per stage gives cumulative totals.\n");
        fprintf(g_timing.file, "# Parent stages include children; do not sum all stages. Merge intervals exclude logging I/O.\n");
        fprintf(g_timing.file, "# Merge seconds are batched since the previous stage report; calls remain cumulative merge invocations.\n");
        fprintf(g_timing.file, "stage\tstep\tseconds\tcalls\tcumulative_seconds\n");
        fflush(g_timing.file);
    }

}

const std::string& DBEMOutputDirectory()
{
	_mkdir("output");
	_mkdir(g_outputDirectory.c_str());
	return g_outputDirectory;
}

std::string DBEMOutputPath(const char* filename)
{
	const char* leaf = filename ? filename : "";
	if (_strnicmp(leaf, "output\\", 7) == 0 || _strnicmp(leaf, "output/", 7) == 0)
		leaf += 7;
	while (*leaf == '\\' || *leaf == '/')
		++leaf;
	return DBEMOutputDirectory() + "\\" + leaf;
}
