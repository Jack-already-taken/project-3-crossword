#include "grid.h"
#ifndef _WIN32
#include <unistd.h>
#endif
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// empty cell: stored as '.'
// can be read from gridfile as '.' or ' '
// black cell: stored as '#'

Grid readGridfile(const std::string &gridfile) {
  std::ifstream fp(gridfile.c_str());
  if (!fp.is_open()) {
    std::cerr << "Error: " << gridfile << " cannot be opened\n";
    exit(1);
  }

  int rows, cols;
  fp >> rows >> cols;
  if (!fp || (cols <= 0) || (rows <= 0)) {
    std::cerr << "Grid read of " << gridfile << " failed on cols/rows"
              << std::endl;
    exit(1);
  }

  Grid grid(rows, cols);

  char c;
  std::string restofline; // ignored
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      fp >> c;
      if (c == ' ') {
        c = '.';
      }
      grid.reallySetCell(c, row, col);
    }
    std::getline(fp, restofline);
    if (!fp) {
      std::cerr << "Grid read of " << gridfile << " failed on row " << row
                << std::endl;
      exit(1);
    }
  }

  return grid;
}

Grid::Grid(int rows_, int cols_)
    : cols(cols_), rows(rows_), valid(true), cells(cols_ * rows_, ' ') {}

void Grid::setWord(const char *word, int row, int col, char direction) {
  for (int j = 0; j < strlen(word); j++) {
    switch (direction) {
    case 'a':
      setCell(word[j], row, col + j);
      break;
    case 'd':
      setCell(word[j], row + j, col);
      break;
    default:
      std::cerr << "Grid::setWord: invalid direction " << direction
                << std::endl;
      break;
    }
  }
}

bool Grid::checkWord(const char *word, int row, int col, char direction) const {
  // if I try to place word in (r,c) in a direction, will it work?
  if ((direction != 'a') && (direction != 'd')) {
    std::cerr << "Grid::checkWord: invalid direction " << direction
              << std::endl;
    exit(1);
  }
  for (int j = 0; j < strlen(word); j++) {
    switch (direction) {
    case 'a':
      if (checkForCellConflict(word[j], row, col + j)) {
        return false;
      }
      break;
    case 'd':
      if (checkForCellConflict(word[j], row + j, col)) {
        return false;
      }
      break;
    }
  }
  return true;
}

bool Grid::checkForCellConflict(char c, int row, int col) const {
  char existing = (*this)(row, col);
  return (existing == '#') || ((existing != '.') && (existing != c));
}

void Grid::setCell(char c, int row, int col) {
  // checks against what's already there
  char existing = (*this)(row, col);
  if (checkForCellConflict(existing, row, col)) {
    std::cerr << "Grid::set: probable error: writing " << c << " onto "
              << existing << " at (" << row << ", " << col << ")" << std::endl;
  }
  reallySetCell(c, row, col);
}

bool Grid::isFilled() const {
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      // if it's not # or a letter, the grid is not filled yet
      char entry = (*this)(r, c);
      if (entry == '#') {
        continue;
      }
      if ((entry < 'a') || (entry > 'z')) {
        return false;
      }
    }
  }
  return true;
}

Extent Grid::findExtent(int row, int col, char direction) const {
  Extent extent;
  switch (direction) {
  case 'a':
    /* row is constant */
    for (extent.min = col;
         /* keep walking left as long as the next cell is still on the grid
            and haven't hit a # yet */
         ((extent.min - 1) >= 0) && ((*this)(row, extent.min - 1) != '#');
         (extent.min)--)
      ;
    for (extent.max = col; ((extent.max + 1) <= cols - 1) &&
                           ((*this)(row, extent.max + 1) != '#');
         (extent.max)++)
      ;
    break;
  case 'd':
    /* col is constant */
    for (extent.min = row;
         ((extent.min - 1) >= 0) && ((*this)(extent.min - 1, col) != '#');
         (extent.min)--)
      ;
    for (extent.max = row; ((extent.max + 1) <= rows - 1) &&
                           ((*this)(extent.max + 1, col) != '#');
         (extent.max)++)
      ;
  }
  return extent;
}

void Grid::outputToFile(const std::string &outputfilePath) const {
  std::ofstream out(outputfilePath);
  if (out.is_open()) {
    out << *this;
    out.close();
  } else {
    std::cerr << "Error opening/creating output file " << outputfilePath << "."
              << std::endl;
  }
}

std::ostream &operator<<(std::ostream &os, const Grid &grid) {
  os << "Grid is (" << grid.rows << ", " << grid.cols << ") and "
     << (grid.valid ? "VALID" : "INVALID") << std::endl;
  for (int r = 0; r < grid.rows; r++) {
    for (int c = 0; c < grid.cols; c++) {
      os << grid(r, c);
    }
    os << std::endl;
  }
  return os;
}
