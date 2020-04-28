#ifndef QTARCORE_H_
#define QTARCORE_H_

#include <pthread.h>
#include <string.h>
#include <ctype.h>
#include "qtar_opts.h"
#include "qtar_subs.h"
#include "sw/ssw_cpp.h"

int _target();

struct Locate {
  string nam;
  int beg, end;
  bool operator<(const Locate& ob) const {
    return (beg < ob.beg) || (beg == ob.beg && nam < ob.nam);
  }
};

static map<int, vector<Locate> > loc;
static pthread_mutex_t mutex;
static string rec;
extern Option opt;
extern Sequence seq;
extern ofstream ofs;
extern map<string, string> ref;
extern map<string, string> anchor;
extern map<string, float> energy;
extern map<char, char> base;

#endif

