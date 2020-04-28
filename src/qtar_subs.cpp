#include "qtar_subs.h"

void _usage()
{
  printf ("\nProgram: qtar (Quickly Target for smRNA)\n");
  printf ("Contact: Qianhua ZHU <zhuqianhua@bgi.com>\n\n");
  printf ("Usage  : qtar [options]\n\n");
  printf ("Options: -s, --smrna    *<s>  small rna sequence in fasta format\n");
  printf ("         -t, --target   *<s>  target sequence in fasta format\n");
  printf ("         -o, --output    <s>  output, default ./qtar_output.xls\n");
  printf ("         -m, --mode      <s>  mode of target, a for animal and p for plant\n");
  printf ("                               a same as: -l 6 -b 6 -a 15 -e -12\n");
  printf ("                               p same as: -S -l 7 -b 2 -a 15 -e -9\n");
  printf ("         -S, --strict    <b>  only involve pairing of the miRNA seed\n");
  printf ("         -l, --length    <i>  length of anchor, default 6\n");
  printf ("         -p, --process   <i>  number of threads, default 6\n");
  printf ("         -b, --bubble    <i>  length of bubbles, default 6\n");
  printf ("         -a, --aln       <f>  cutoff of align score, default 10.0\n");
  printf ("         -e, --mfe       <f>  cutoff of minimal free energy, default -8.0\n");
  printf ("         -q, --quiet     <b>  print nothing except serious errors\n");
  printf ("         -h, --help      <b>  print this information\n\n");
  exit(0);
}

void _getopt(int argc, char* argv[])
{
  char const * shortOpt = "s:t:o:m:Sl:p:b:a:e:qh";
  struct option longOpt[] = {
    {"smrna", 1, NULL, 's'},
    {"target", 1, NULL, 't'},
    {"output", 1, NULL, 'o'},
    {"mode", 1, NULL, 'm'},
    {"strict", 0, NULL, 'S'},
    {"length", 1, NULL, 'l'},
    {"process", 1, NULL, 'p'},
    {"bubble", 1, NULL, 'b'},
    {"aln", 1, NULL, 'a'},
    {"mfe", 1, NULL, 'e'},
    {"quiet", 0, NULL, 'q'},
    {"help", 0, NULL, 'h'},
    {NULL, 0, NULL, 0},
  };
  int nextOpt;
  while ((nextOpt = getopt_long(argc, argv, shortOpt, longOpt, NULL)) != -1) {
    switch (nextOpt) {
      case 's':
        opt.mir = optarg;
        break;
      case 't':
        opt.tar = optarg;
        break;
      case 'o':
        opt.out = optarg;
        break;
      case 'm':
        opt.mod = optarg;
        break;
      case 'S':
        opt.strict = true;
        break;
      case 'l':
        opt.len = atoi(optarg);
        break;
      case 'p':
        opt.thd = atoi(optarg);
        break;
      case 'b':
        opt.bub = atoi(optarg);
        break;
      case 'a':
        opt.aln = atof(optarg);
        break;
      case 'e':
        opt.mfe = atof(optarg);
        break;
      case 'q':
        opt.quiet = true;
        break;
      case 'h':
        opt.help = true;
        break;
    }
  }
  if (opt.help == true or opt.mir.compare("") == 0 or opt.tar.compare("") == 0)
    _usage();
  if (opt.len < 0)
    opt.len = 5;
  if (opt.mod.compare("a") == 0) {
    opt.len = 6;
    opt.aln = 15;
    opt.mfe = -12;
    opt.bub = 6;
  } else if (opt.mod.compare("p") == 0) {
    opt.len = 7;
    opt.aln = 15;
    opt.mfe = -9;
    opt.bub = 2;
    opt.strict = true;
  }
}

void _info(string des)
{
  string str = "[%Y-%m-%d %X] " + des;
  time_t t = time(0);
  char tmp[64];
  strftime (tmp, sizeof(tmp), str.c_str(), localtime(&t));
  puts(tmp);
}

void _index(string &nam, string &seq)
{
  // duplicates in the reference
  if (ref.find(nam) != ref.end()) {
    if (dups.find(nam) != dups.end())
      dups[nam] += 1;
    else
      dups[nam] = 1;
    stringstream ss;
    ss << dups[nam];
    nam = nam + "_" + ss.str();
  }
  for (int i = 0; i <= seq.length()-opt.len; i ++) {
    stringstream ss;
    ss << i;
    string tac = seq.substr(i, opt.len);
    string tps = nam + '%' + ss.str();
    if (anchor.find(tac) == anchor.end())
      continue;
    if (anchor[tac].compare("") == 0)
      anchor[tac] = tps;
    else 
      anchor[tac] += ' ' + tps;
  }
  ref[nam] = seq;
  seq = "";
}

int _fasta(string &str) 
{
  int flag = 0;
  if (str.find_first_of('>') != string::npos) {
    int end = str.length();
    if (str.find_first_of(' ') != string::npos)
      end = str.find_first_of(' ');
    str = str.substr(str.find_first_of('>') + 1, end - str.find_first_of('>') - 1);
    return 0;
  } else {
    transform(str.begin(), str.end(), str.begin(), ::toupper);
    return 1;
  }
}

void _mfe()
{
  // reverse and complementary
  base['A'] = 'T'; base['T'] = 'A'; base['C'] = 'G'; base['G'] = 'C'; base['U'] = 'A';
  // Initialization: Freier et al., Improved free-energy parameters for predictions of RNA duplex stability
  // Table 2
  energy["AATT"] = -0.9; energy["ATAT"] = -0.9; energy["TATA"] = -1.1; energy["CATG"] = -1.8;
  energy["CTAG"] = -1.7; energy["GATC"] = -2.3; energy["GTAC"] = -2.1; energy["CGCG"] = -2.0;
  energy["GCGC"] = -3.4; energy["GGCC"] = -2.9;
  // Table 3
  energy["AAXT"] = -0.8; energy["ACXT"] = -0.5; energy["AGXT"] = -0.8; energy["ATXT"] = -0.6;
  energy["CAXG"] = -1.7; energy["CCXG"] = -0.8; energy["CGXG"] = -1.7; energy["CTXG"] = -1.2;
  energy["GAXC"] = -1.1; energy["GCXC"] = -0.4; energy["GGXC"] = -1.3; energy["GTXC"] = -0.6;
  energy["TAXA"] = -0.7; energy["TCXA"] = -0.1; energy["TGXA"] = -0.7; energy["TTXA"] = -0.1;
  energy["AATX"] = -0.3; energy["CATX"] = -0.3; energy["GATX"] = -0.4; energy["TATX"] = -0.2;
  energy["ACGX"] = -0.5; energy["CCGX"] = -0.2; energy["GCGX"] = -0.2; energy["TCGX"] = -0.1;
  energy["AGCX"] = -0.2; energy["CGCX"] = -0.3; energy["GGCX"] = -0.0; energy["TGCX"] = -0.0;
  energy["ATAX"] = -0.3; energy["CTAX"] = -0.2; energy["GTAX"] = -0.2; energy["TTAX"] = -0.2;
  // Table 4
  energy["GAAC"] = -1.1; energy["GACC"] = -1.3; energy["GAGC"] = -1.3; energy["GATC"] = -2.3;
  energy["GCAC"] = -1.1; energy["GCCC"] = -0.6; energy["GCGC"] = -3.4; energy["GCTC"] = -0.5;
  energy["GGAC"] = -1.6; energy["GGCC"] = -2.9; energy["GGGC"] = -1.4; energy["GGTC"] = -1.4;
  energy["GTAC"] = -2.1; energy["GTCC"] = -0.8; energy["GTGC"] = -2.3; energy["GTTC"] = -0.7;
  energy["CAAG"] = -1.9; energy["CACG"] = -2.0; energy["CAGG"] = -1.9; energy["CATG"] = -1.8;
  energy["CCAG"] = -1.0; energy["CCCG"] = -1.1; energy["CCGG"] = -2.9; energy["CCTG"] = -0.8;
  energy["CGAG"] = -1.9; energy["CGCG"] = -2.0; energy["CGGG"] = -1.9; energy["CGTG"] = -1.6;
  energy["CTAG"] = -1.7; energy["CTCG"] = -1.5; energy["CTGG"] = -1.9; energy["CTTG"] = -1.2;
  energy["AAAT"] = -0.8; energy["AACT"] = -1.0; energy["AAGT"] = -1.0; energy["AATT"] = -0.9;
  energy["ACAT"] = -0.7; energy["ACCT"] = -0.7; energy["ACGT"] = -2.1; energy["ACTT"] = -0.7;
  energy["AGAT"] = -0.8; energy["AGCT"] = -1.7; energy["AGGT"] = -1.0; energy["AGTT"] = -0.9;
  energy["ATAT"] = -0.9; energy["ATCT"] = -0.8; energy["ATGT"] = -0.9; energy["ATTT"] = -0.8;
  energy["TAAA"] = -1.0; energy["TACA"] = -0.8; energy["TAGA"] = -1.1; energy["TATA"] = -1.1;
  energy["TCAA"] = -0.7; energy["TCCA"] = -0.6; energy["TCGA"] = -2.3; energy["TCTA"] = -0.5;
  energy["TGAA"] = -1.1; energy["TGCA"] = -1.8; energy["TGGA"] = -1.2; energy["TGTA"] = -0.9;
  energy["TTAA"] = -0.9; energy["TTCA"] = -0.6; energy["TTGA"] = -1.0; energy["TTTA"] = -0.5;
  // Table 5
  energy["AGTT"] = -0.5; energy["ATGT"] = -0.7;
  energy["CGTG"] = -1.5; energy["CTGG"] = -1.5;
  energy["GGTC"] = -1.3; energy["GTGC"] = -1.9;
  energy["TGTA"] = -0.7; energy["TTGA"] = -0.5;
  energy["GGTT"] = -0.5; energy["GTGT"] = -0.5;
  energy["TGTG"] = -0.6; energy["TTGG"] = -0.5;
  // Table 6
  energy["XX"] = 0.8; energy["XXX"] = 1.3; energy["XXXX"] = 1.7; energy["XXXXX"] = 2.1;
  energy["XXXXXX"] = 2.5; energy["XXXXXXX"] = 2.6; energy["XXXXXXXX"] = 2.8;
  energy["XXXXXXXXX"] = 3.1;
}

string _revc(string &str, int flag)
{
  string rev(str);
  reverse(rev.begin(), rev.end());
  if (flag == 1) {
    for (int i = 0; i < str.size(); i++)
      if (base.find(rev[i]) != base.end())
        rev[i] = base[rev[i]];
  }
  return rev;
}

void _anc(int typ)
{
  string inf;
  if (typ == 0)
    inf = opt.mir;
  else 
    inf = opt.tar;
  ifstream ifs(inf.c_str());
  string nam(""), seq(""), str("");
  while (!ifs.eof())
  {
    getline(ifs, str);
    int flag = _fasta(str);
    if (flag == 0) {
      if (nam.compare("") != 0 && seq.length() > opt.len) {
        if (typ == 0) {
          string tpm = _revc(seq, 1);
          for (int i = 0; i <= tpm.length()-opt.len; i ++)
            anchor[tpm.substr(i, opt.len)] = "";
        } else { 
          num ++;
          if (opt.quiet != true && num % 10000 == 0) {
            stringstream ss;
            ss << num;
            _info("anchor " + ss.str() + " targets ...");
          }
          _index(nam, seq);
        }
      }
      nam = str;
    } else {
      seq += str;
    }
  }
  ifs.close();
  if (seq.length() > opt.len)
    if (typ == 0)
      for (int i = 0; i <= seq.length()-opt.len; i ++)
        anchor[seq.substr(i, opt.len)] = "";
    else {
      if (opt.quiet != true) {
        stringstream ss;
        ss << num;
        _info("anchor " + ss.str() + " targets ...");
      }
      _index(nam, seq);
    }
} 

void _anchor()
{
  if (opt.quiet != true)
    _info("build anchor sets ...");
  _anc(0);
  _anc(1);
  dups.clear();
  if (opt.quiet != true)
    _info("initialize mfe score ...");
  _mfe();
}
