#ifndef QTAROPTS_H_
#define QTAROPTS_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

using namespace std;

struct Option {
  string mir, tar, out, mod;
  int len, thd, bub;
  float aln, mfe;
  bool strict, quiet, help;
};

struct Sequence {
  string nam, seq;
};

#endif
