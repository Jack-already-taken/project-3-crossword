#include "options.h"

static Options::Seed parseSeed(const std::string &arg) {
  std::vector<std::string> v;
  size_t i = 0;
  for (size_t j = arg.find(','); j != std::string::npos; j = arg.find(',', i)) {
    v.push_back(arg.substr(i, j - i));
    i = j + 1;
  }
  v.push_back(arg.substr(i));

  // TODO: check correct number of arguments
  // TODO: check these are integers, use proper C++ to parse them (maybe
  // istringstream?)
  int row = std::stoi(v[0]);
  int col = std::stoi(v[1]);
  char direction = v[2][0];
  std::string word = v[3];

  return {row, col, direction, word};
}

Options parseCommandLineOptions(int argc, char **argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    const std::string optarg(i + 1 < argc ? argv[i + 1] : "");
    if (arg == "-m") {
      options.minimumScore = std::stoi(optarg);
      ++i;
    } else if (arg == "-g") {
      options.gridfilePath = optarg;
      ++i;
    } else if (arg == "-w") {
      options.wordlistPath = optarg;
      ++i;
    } else if (arg == "-s") {
      options.seeds.push_back(parseSeed(optarg));
      ++i;
    } else if (arg == "-o") {
      options.outputfilePath = optarg;
    } else {
      // report bad argument
    }
  }
  return options;
}
