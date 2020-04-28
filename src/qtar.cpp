#include "qtar.h"

int main(int argc, char* argv[])
{
  _getopt(argc, argv);
  _anchor(); 
  if (opt.quiet != true)
    _info("target and output ...");
  ofs.open(opt.out.c_str()); 
  ofs << "miRNA\tTarget\tPosition\tSeed type\tMFE\tAlign score" << endl;
  ifstream ifs(opt.mir.c_str());
  string str;
  while (!ifs.eof())
  {
    getline(ifs, str);
    int flag = _fasta(str);
    if (flag == 0) {
      if (seq.nam.compare("") != 0)
        _target();
      seq.nam = str; 
      seq.seq = "";
    } else {
      seq.seq += str;
    }
  }
  ifs.close();
  _target();
  ofs.close();
  if (opt.quiet != true)
    _info("done ...");
  return 0;
}

