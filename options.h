#pragma once

#include <string>
#include <vector>

struct Options {
  struct Seed {
    int row;
    int col;
    char direction;
    std::string word;
  };

  int minimumScore = 0;
  std::string wordlistPath;
  std::string gridfilePath;
  std::string outputfilePath;
  std::vector<Seed> seeds;
};

Options parseCommandLineOptions(int argc, char **argv);
