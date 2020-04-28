#include "qtar_opts.h"

Option opt = {"", "", "./qtar_output.xls", "",
              6, 6, 6, 
	          10.0, -8,
             };

ofstream ofs;
Sequence seq = {"", ""};
map<string, string> ref;
map<string, string> anchor;
map<string, float> energy;
map<char, char> base;
