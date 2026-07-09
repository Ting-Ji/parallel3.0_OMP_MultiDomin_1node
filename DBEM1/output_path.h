#pragma once

#include <string>

void DBEMInitializeOutputDirectory(const char* modelName, long NStep);
const std::string& DBEMOutputDirectory();
std::string DBEMOutputPath(const char* filename);
