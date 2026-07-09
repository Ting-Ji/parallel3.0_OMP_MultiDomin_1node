#include "output_path.h"

#include <direct.h>
#include <cctype>
#include <cstdio>
#include <cstring>

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
