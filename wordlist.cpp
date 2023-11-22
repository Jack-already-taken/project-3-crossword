#include "wordlist.h"
#include <fstream>
#include <iostream>

MasterWordlist readMasterWordlistFromFile(const std::string &filename,
                                          int minimumScore) {
  MasterWordlist retval;

  std::ifstream is(filename);
  std::string line;
  while (std::getline(is, line)) {
    const size_t i = line.find(';');
    if (i == std::string::npos) {
      std::cerr << "skipping bad line: " << line << std::endl;
      continue;
    }

    std::string word = line.substr(0, i);
    const int score = std::stoi(line.substr(i + 1));
    if (score >= minimumScore) {
      retval.emplace_back(word);
    }
  }

  return retval;
}
