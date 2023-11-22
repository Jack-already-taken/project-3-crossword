#pragma once
#include "wordlist.h"
#include <iostream>

/* Extent represents a range [min, max]
 *
 * Intent is that extent reflects [min, max], i.e, min represents the
 * first letter in the word and max represents the last, so loop from
 * min <= index <= max */
class Extent {
public:
  Extent() : min(-1), max(-1) {}
  Extent(int min_, int max_) : min(min_), max(max_) {}
  int min;
  int max;
};

class Grid {
public:
  Grid(int rows_, int cols_);

  int getNumRows() const { return rows; }
  int getNumCols() const { return cols; }

  char &operator()(int row, int col) { return cells[row * cols + col]; }
  char operator()(int row, int col) const { return cells[row * cols + col]; }

  // both set and reallySet set a particular grid cell to a letter
  // - set checks to make sure we're not overwriting a different letter
  // - reallySet does not do that check
  void reallySetCell(char c, int row, int col) { (*this)(row, col) = c; }
  void setCell(char c, int row, int col);

  // write an entire word
  void setWord(const char *word, int row, int col, char direction);

  // checkWord says if it is possible to place a word in the current grid
  // at the specified location (it checks to see if it conflicts with
  // the letters that are already in the grid)
  bool checkWord(const char *word, int row, int col, char direction) const;

  // given a cell in the grid and a direction, how far does the word that
  // crosses that cell extend?
  Extent findExtent(int row, int col, char direction) const;

  // is every cell in this grid filled?
  bool isFilled() const;

  bool isValid() const { return valid; }
  void invalidate() { valid = false; }

  void outputToFile(const std::string &filename) const;

  friend std::ostream &operator<<(std::ostream &os, const Grid &grid);

private:
  bool checkForCellConflict(char c, int row, int col) const;
  int cols, rows;
  /* a grid is "invalid" (valid == false) if either:
   * - there is a full word in the grid that is not in the wordlist
   * - the program has concluded that the grid cannot be completed (there
   *   are no solutions possible)
   */
  bool valid;
  std::vector<char> cells; // actually stores the grid data
};

// allocates and initializes grid, from an input file
// also sets cols and rows
Grid readGridfile(const std::string &gridfile);
