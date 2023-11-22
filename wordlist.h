#pragma once

#include <string>
#include <vector>

using MasterWordlist = std::vector<std::string>;

MasterWordlist readMasterWordlistFromFile(const std::string &filename,
                                          int minimumScore);
