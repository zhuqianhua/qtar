#ifndef QTARSUBS_H_
#define QTARSUBS_H_

#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include "qtar_opts.h"

void _usage();
void _getopt(int argc, char* argv[]);
int _fasta(string &str);
void _anchor();
string _revc(string &str, int flag);
void _info(string des);

static map<string, int> dups;
static int num;
extern Option opt;
extern map<string, string> ref;
extern map<string, string> anchor;
extern map<string, float> energy;
extern map<char, char> base;

#endif

