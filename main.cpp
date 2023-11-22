#ifndef _WIN32
#include <time.h> // https://stackoverflow.com/questions/459691/best-timing-method-in-c
#endif
#include "grid.h"
#include "options.h"
#include "wordlist.h"
#include <cassert>
#include <fstream>
#include <map>
#include <vector>
#include <memory>
#include <cilk/cilk.h>

// #define FANCY_ALLOCATION
// uncomment this to use the non-thread-safe fancy allocator

/* everything in this project is addressed as (row, col) */

/***** DATATYPES *****/

// we can use "const char*" as our Word type because we know the word list
// read from the dictionary will be preserved, so these raw pointers will
// remain valid. This eliminates a potential source of a lot of overhead
// (string copying) in some cases, and at a minimum makes it easier to profile
// the code and understand where the "important" performance is going.

using Word = const char *; // using Word = std::string;

using Wordlist = std::vector<Word>;

using SearchSpace = std::vector<std::shared_ptr<const Wordlist>>;

// <entry number, character position in that entry>
using Constraint = std::pair<size_t, size_t>;
using Constraints = std::vector<Constraint>;

struct Entry {
  Entry(int r, int c, char d, const Constraints &constraints_)
      : row(r), col(c), direction(d), constraints(constraints_) {}
  int row;
  int col;
  char direction;
  Constraints constraints;
};

using Entries = std::vector<Entry>;

/***** CORE LOGIC *****/

Entries constructEntriesForGrid(const Grid &grid) {
  const int n_rows(grid.getNumRows());
  const int n_cols(grid.getNumCols());

  // the first entry is a dummy entry of zero length
  const Entry nullEntry(-1, -1, ' ', Constraints());

  Entries entries(1, nullEntry);

  std::vector<std::vector<std::pair<size_t, size_t>>> entryMap(n_rows);
  for (int row = 0; row < n_rows; ++row) {
    entryMap[row].resize(n_cols);
  }

  for (int row = 0; row < n_rows; ++row) {
    for (int col = 0; col < n_cols; ++col) {
      if (grid(row, col) == '#')
        continue;

      const Extent across = grid.findExtent(row, col, 'a');
      if (across.min == col) {
        const int len = across.max - across.min + 1;
        const size_t entry_no = entries.size();
        entries.push_back(Entry(row, col, 'a', Constraints(len)));
        for (int x = across.min; x <= across.max; ++x) {
          entryMap[row][x].first = entry_no;
        }
      }

      const Extent down = grid.findExtent(row, col, 'd');
      if (down.min == row) {
        const int len = down.max - down.min + 1;
        const size_t entry_no = entries.size();
        entries.push_back(Entry(row, col, 'd', Constraints(len)));
        for (int y = down.min; y <= down.max; ++y) {
          entryMap[y][col].second = entry_no;
        }
      }
    }
  }

  // every across clue is checked by a down clue, which means all entries in
  // the EntryMap should be non-zero
  for (int row = 0; row < n_rows; ++row) {
    for (int col = 0; col < n_cols; ++col) {
      if (grid(row, col) == '#')
        continue;
      assert(entryMap[row][col].first != 0 && entryMap[row][col].second != 0);
    }
  }

  // fill in the constraints for each entry, we skip the dummy entry which is
  // just there so 0 is a useful sentinel
  for (size_t i = 1; i != entries.size(); ++i) {
    Entry &entry(entries[i]);
    for (size_t j = 0; j != entry.constraints.size(); ++j) {
      if (entry.direction == 'a') {
        const size_t otherEntryNo = entryMap[entry.row][entry.col + j].second;
        const size_t otherEntryLetter = entry.row - entries[otherEntryNo].row;
        entry.constraints[j] = Constraint(otherEntryNo, otherEntryLetter);
      } else {
        const size_t otherEntryNo = entryMap[entry.row + j][entry.col].first;
        const size_t otherEntryLetter = entry.col - entries[otherEntryNo].col;
        entry.constraints[j] = Constraint(otherEntryNo, otherEntryLetter);
      }
    }
  }

  return entries;
}

#ifdef FANCY_ALLOCATION
// this is a really clunky allocator meant to help insure that we aren't
// needlessly allocating memory repeatedly.  When a Wordlist is no longer
// referenced anywhere it gets put on the freelist where it can be reused by
// later allocations.  Note that this isn't thread safe in any way!

static std::vector<std::unique_ptr<Wordlist>> freelist;
static int allocCount = 0; // counts new Wordlist allocations
static int reuseCount = 0; // counts number of reused Wordlists

void releaseWordlist(Wordlist *p) {
  freelist.emplace_back(std::unique_ptr<Wordlist>(p));
}

std::shared_ptr<Wordlist> allocateWordlist() {
  std::unique_ptr<Wordlist> p;
  if (freelist.empty()) {
    ++allocCount;
    p = std::make_unique<Wordlist>();
  } else {
    ++reuseCount;
    std::swap(p, freelist.back());
    freelist.pop_back();
    p->resize(0);
  }

  return std::shared_ptr<Wordlist>(p.release(), releaseWordlist);
}
#else
std::shared_ptr<Wordlist> allocateWordlist() {
  return std::make_shared<Wordlist>();
}
#endif

SearchSpace
initializeSearchSpaceForEntries(const Entries &entries, const Grid &grid,
                                const MasterWordlist &masterWordlist) {
  // starts with an empty wordlist for the initial dummy entry
  SearchSpace searchSpace(1, allocateWordlist());

  for (size_t i = 1; i != entries.size(); ++i) {
    const Entry &entry(entries[i]);
    const size_t entryLen(entry.constraints.size());

    using GridConstraint = std::pair<char, size_t>;
    std::vector<GridConstraint> gridConstraints;
    for (int j = 0; j != entryLen; ++j) {
      const char gridLetter = (entry.direction == 'a')
                                  ? grid(entry.row, entry.col + j)
                                  : grid(entry.row + j, entry.col);
      if (gridLetter != '.') {
        gridConstraints.push_back(GridConstraint(gridLetter, j));
      }
    }

    auto thisWordlist = allocateWordlist();
    for (const auto &word :
         masterWordlist) { // loop through all words in the master wordlist ...
      if (word.size() ==
          entryLen) { // ... and only consider those of the right length
        bool okay = true;
        for (const auto &gc : gridConstraints) { // loop through all constraints
          if (word[gc.second] != gc.first) {
            okay = false;
            break;
          }
        }
        if (okay) {
          thisWordlist->push_back(word.c_str());
        }
      }
    }
    searchSpace.push_back(thisWordlist);
  }

  return searchSpace;
}

std::shared_ptr<const Wordlist> constrainWordlist(const Wordlist &srcList,
                                                  char letter, size_t pos) {
  auto dstList = allocateWordlist();

  for (const auto &word : srcList) {
    if (word[pos] == letter) {
      dstList->push_back(word);
    }
  }

  return dstList;
}

using SolveOrder = std::vector<size_t>;

SolveOrder pickSolveOrder(const Entries &entries) {
  // constrainedLetters: mapping between Entry (actually, Entry's
  // index into `entries`) and the constrained letters for that Entry
  std::map<size_t, size_t> constrainedLetters;
  for (size_t i = 1; i != entries.size(); ++i) {
    constrainedLetters[i] = 0;
  }

  SolveOrder solveOrder;
  while (!constrainedLetters.empty()) {

    size_t maxCount = 0;
    size_t maxEntry = constrainedLetters.begin()->first;
    for (const auto &i : constrainedLetters) {
      if (maxCount < i.second) {
        maxEntry = i.first;
        maxCount = i.second;
      }
    }

    solveOrder.push_back(maxEntry);
    constrainedLetters.erase(maxEntry);

    // for every other unsolved entry that maxEntry crosses, increment its
    // constrained letter count (unknown) letter count by 1
    for (const auto &constraint : entries[maxEntry].constraints) {
      auto j = constrainedLetters.find(constraint.first);
      if (j != constrainedLetters.end()) {
        j->second += 1;
      }
    }
  }

  return solveOrder;
}

class PossibleFitComputer {

public:
  PossibleFitComputer(const Constraints &, const SearchSpace &);

  bool operator()(const Word &word) const;

private:
  std::vector<int> letterMasks;
};

PossibleFitComputer::PossibleFitComputer(const Constraints &constraints,
                                         const SearchSpace &searchSpace)
    : letterMasks(constraints.size(), 0) {
  for (size_t i = 0; i != constraints.size(); ++i) {
    const size_t entryNo(constraints[i].first);
    const size_t entryLetter(constraints[i].second);
    int mask = 0;
    for (const Word &word : *searchSpace[entryNo]) {
      mask |= 1 << (word[entryLetter] - 'a');
    }
    letterMasks[i] = mask;
  }
}

bool PossibleFitComputer::operator()(const Word &word) const {
  const size_t n(letterMasks.size());
  for (size_t i = 0; i != n; ++i) {
    if (((1 << (word[i] - 'a')) & letterMasks[i]) == 0)
      return false;
  }
  return true;
}

class Solver {

public:
  Solver(const Grid &initialGrid, const Entries &, const SolveOrder &);

  const Wordlist &getSolution() const { return solution; }
  bool solve(const SearchSpace &, size_t i = 0);

private:
  std::unique_ptr<SearchSpace>
  computeNewSearchSpaceGivenGuess(const SearchSpace &, const Entry &,
                                  const Word &guess) const;
  void reportProgress() const;

  const Grid &initialGrid;
  const Entries entries;
  const SolveOrder solveOrder;
  Wordlist solution;
  const clock_t startTime;
  const int printInterval;
  clock_t printTime;
  int numCalls;
};

Solver::Solver(const Grid &initialGrid_, const Entries &entries_,
               const SolveOrder &solveOrder_)
    : initialGrid(initialGrid_), entries(entries_), solveOrder(solveOrder_),
      solution(entries_.size()), startTime(clock()),
      printInterval(5 * CLOCKS_PER_SEC), numCalls(0) {
  assert(solveOrder.size() + 1 == entries.size());
  printTime = startTime + printInterval;
}

bool Solver::solve(const SearchSpace &searchSpace, size_t i) {
  // invariant upon calling Solver::solve: the current state is valid
  if (i == solveOrder.size())
    return true;

#if 1
  ++numCalls;
  if (clock() > printTime) {
    printTime += printInterval;
    reportProgress();
  }
#endif

  const size_t entryNo(solveOrder[i]);
  const Entry &entry(entries[entryNo]);

  PossibleFitComputer possibleFit(entry.constraints, searchSpace);

  // for (const auto &word : *searchSpace[entryNo])
  // the above line is the fancy C++ way to write the below two lines
  for (int w = 0; w < searchSpace[entryNo]->size(); w++) {
    Word const &word = (*(searchSpace[entryNo]))[w];

    // don't allow duplicate entries in the solution
    if (std::find(solution.begin(), solution.end(), word) != solution.end())
      continue;

    if (!possibleFit(word))
      continue;

    // build a pruned search space, each of the entries that are intersected
    // by word need to have their dictionaries culled to reflect the new
    // constraints we've added
    std::unique_ptr<SearchSpace> searchSpaceCopy(
        computeNewSearchSpaceGivenGuess(searchSpace, entry, word));

    // returns null if any of the updated dictionaries are empty
    if (!searchSpaceCopy)
      continue;

    solution[entryNo] = word;

    if (solve(*searchSpaceCopy, i + 1))
      return true;

    // word didn't work, so clear it out so we don't confuse anybody into
    // thinking it's currently taken by entryNo
    solution[entryNo] = "";
  }

  return false;
}

void Solver::reportProgress() const {
  Grid grid(initialGrid);
  for (size_t i = 1; i != entries.size(); ++i) {
    const Entry &entry(entries[i]);
    if (solution[i]) {
      grid.setWord(solution[i], entry.row, entry.col, entry.direction);
    }
  }

  const double elapsed(double(clock() - startTime) / CLOCKS_PER_SEC);

  std::cout << "Current status (" << numCalls << " calls, " << elapsed
            << " seconds): " << std::endl
            << grid;
}

std::unique_ptr<SearchSpace>
Solver::computeNewSearchSpaceGivenGuess(const SearchSpace &currentSearchSpace,
                                        const Entry &entry,
                                        const Word &guess) const {
  std::unique_ptr<SearchSpace> retval(
      std::make_unique<SearchSpace>(currentSearchSpace));

  const size_t wordLen(entry.constraints.size());
  for (size_t j = 0; j != wordLen; ++j) {
    const Constraint &constraint(entry.constraints[j]);
    auto srcList = currentSearchSpace[constraint.first];
    auto dstList = constrainWordlist(*srcList, guess[j], constraint.second);
    if (dstList->empty()) {
      // if we got rid of all possible words for a given entry clearly
      // this word can't be correct, as we won't be able to solve the
      // remaining entries
      return nullptr;
    }

    (*retval)[constraint.first] = dstList;
  }

  return retval;
}

/* this print routine is for debugging only */
void printOutSeedsForGrid(const Grid &grid) {
  std::cout << "seeds for grid:";

  const Entries entries = constructEntriesForGrid(grid);

  for (size_t i = 1; i != entries.size(); ++i) {
    const Entry &entry(entries[i]);

    std::string word;
    bool okay = true;
    for (int j = 0; j != entry.constraints.size(); ++j) {
      const char letter = entry.direction == 'a'
                              ? grid(entry.row, entry.col + j)
                              : grid(entry.row + j, entry.col);
      if (letter == '.') {
        okay = false;
      } else {
        word += letter;
      }
    }

    if (okay) {
      std::cout << " -s " << entry.row << "," << entry.col << ","
                << entry.direction << "," << word;
    }
  }
  std::cout << std::endl;
}

Grid fillGrid(const Grid &initialGrid, MasterWordlist &masterWordlist) {
  const Entries entries = constructEntriesForGrid(initialGrid);
  const SearchSpace searchSpace =
      initializeSearchSpaceForEntries(entries, initialGrid, masterWordlist);

  Solver solver(initialGrid, entries, pickSolveOrder(entries));

  Grid grid(initialGrid);
  if (solver.solve(searchSpace)) {
    const Wordlist &solution(solver.getSolution());
    for (size_t i = 1; i != entries.size(); ++i) {
      const Entry &entry(entries[i]);
      const std::string &word(solution[i]);

      assert(
          grid.checkWord(word.c_str(), entry.row, entry.col, entry.direction));
      grid.setWord(word.c_str(), entry.row, entry.col, entry.direction);
    }
  } else {
    grid.invalidate();
  }

  return grid;
}

int main(int argc, char **argv) {
  /* begin: parse command line options, initialize data structures */
  const Options options = parseCommandLineOptions(argc, argv);

  MasterWordlist masterWordlist;
  if (!options.wordlistPath.empty()) {
    masterWordlist =
        readMasterWordlistFromFile(options.wordlistPath, options.minimumScore);
  }

  if (options.gridfilePath.empty()) {
    std::cerr << "Must specify gridfile with -g <gridfile>" << std::endl;
    return 1;
  }

  Grid initialGrid = readGridfile(options.gridfilePath);

  for (const auto &seed : options.seeds) {
    std::cout << "Seed: " << seed.row << " " << seed.col << " "
              << seed.direction << " " << seed.word << std::endl;
    initialGrid.setWord(seed.word.c_str(), seed.row, seed.col, seed.direction);

    if (std::find(masterWordlist.begin(), masterWordlist.end(), seed.word) ==
        masterWordlist.end()) {
      std::cout << "  Adding \"" << seed.word << "\" to word list" << std::endl;
      masterWordlist.push_back(seed.word);
    }
  }
  /* end: parse command line options, initialize data structures */

  std::cout << "Input grid: " << initialGrid;
  const clock_t startTime = clock();

  /* this call fills the grid. */
  Grid grid = fillGrid(initialGrid, masterWordlist);

  /* report the results */
  const clock_t finishTime = clock();
  const double elapsed = double(finishTime - startTime) / CLOCKS_PER_SEC;
  std::cout << "Elapsed Time: " << elapsed << " seconds" << std::endl;
#ifdef FANCY_ALLOCATION
  std::cout << "allocCount=" << allocCount << ", reuseCount=" << reuseCount
            << std::endl;
#endif

  if (grid.isValid() && grid.isFilled()) {
    std::cout << "Solution found." << std::endl << grid << std::endl;
  } else {
    std::cout << "No solution found." << std::endl;
  }
  if (!options.outputfilePath.empty()) {
    grid.outputToFile(options.outputfilePath);
  }
  return 0;
}
