#include "qtar_core.h"

vector<string> _split(string str, char delim)
{
  vector<string> vs;
  string temp;
  for (int i = 0; i < str.length(); ++ i) {
    if (str[i] != delim)
      temp += str[i];
    else {
      if (temp.compare("") != 0)
        vs.push_back(temp);
      temp = "";
    }
  }
  if (temp.compare("") != 0)
    vs.push_back(temp);
  return vs;
}

void _locate()
{
  loc.clear();
  map<Locate, char> pis;
  Locate tag;
  for (int i = 0; i <= rec.length()-opt.len; i ++) {
    string loc = rec.substr(i, opt.len);
    if (anchor.find(loc) == anchor.end())
      continue;
    vector<string> vanc = _split(anchor[loc], ' ');
    for (vector<string>::iterator iv = vanc.begin(); iv != vanc.end(); iv++) {
      vector<string> vtar = _split(*iv, '%');
      tag.beg = atoi(vtar[1].c_str())-i-opt.bub;
      tag.end = atoi(vtar[1].c_str())+seq.seq.length()-i-1+opt.bub;
      if (tag.beg < 0)
        tag.beg = 0;
      if (tag.end > ref[vtar[0]].length()-1)
        tag.end = ref[vtar[0]].length()-1;
      tag.nam = vtar[0];
      pis[tag] = ' '; 
    }
  }
  int i = 0; 
  for (map<Locate, char>::iterator im = pis.begin(); im != pis.end(); im++) {
    loc[i % opt.thd].push_back(im->first);  
    i ++;
  }
  pis.clear();
}

vector<int> _mfe(string &seq) {
  int beg = 0, len = 0;
  vector<int> vs;
  string tmp;
  while (1)
  {
    if (beg >= seq.size() - 1)
      break;
    if (len == 0) {
      if (seq[beg] == 'X')
        len ++;
      else {
        vs.push_back(beg);
        vs.push_back(2);
        beg += 2;
      }
    } else {
      if (seq[beg+len] == 'X')
        len ++;
      else {
        if (len > 1) {
          vs.push_back(beg);
          vs.push_back(len);
          beg += len;
        }
        len = 0;
        if (beg >= seq.size() - 1)
          break;
        vs.push_back(beg);
        vs.push_back(2);
        beg += 2;
      }
    }
  }
  return vs;
}

void *_score(void *thread) {
  int THD = *(int *) thread;
  if (loc.find(THD) == loc.end())
    return 0;
  for (vector<Locate>::iterator iv = loc[THD].begin(); iv != loc[THD].end(); iv++) {
    string tar = ref[iv->nam].substr(iv->beg, iv->end - iv->beg + 1);
    int32_t mask = strlen(tar.c_str()) / 2;
    mask = mask < opt.len ? opt.len : mask;
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment align;
    string sqt = _revc(seq.seq, 0);
    aligner.Align(rec.c_str(), tar.c_str(), tar.size(), filter, &align, mask); 
    /* Threshold 0, clip
    if (opt.strict == true) {
      int five = align.cigar_string.find_first_of('S'), third = align.cigar_string.find_last_of('S');
      if (five != string::npos && third != string::npos && five != third)
        continue;
    } */
    // Threshold 1, alignment score
    if (align.sw_score < opt.aln)
      continue;
    // Decode the alignment
    string dgt, tpt, tpm;
    int pst = align.ref_begin, psm = 0, clp = 0;
    int five = 0, third = 0, bub = 0;
    for (int i = 0; i < align.cigar_string.size(); i++) {
      if (isdigit(align.cigar_string[i]))
        dgt += align.cigar_string[i];
      else {
        string tmp(atoi(dgt.c_str()), 'X');
        switch(align.cigar_string[i]) 
        {
          case 'S':
            clp = atoi(dgt.c_str());   
            psm += atoi(dgt.c_str());
            if (i == align.cigar_string.size() - 1)
              third = atoi(dgt.c_str());
            else
              five = atoi(dgt.c_str());
            break;
          case '=':
            tpt += tar.substr(pst, atoi(dgt.c_str()));
            tpm += sqt.substr(psm, atoi(dgt.c_str()));
            pst += atoi(dgt.c_str());
            psm += atoi(dgt.c_str());		
            break;
          case 'I':
            tpt += tmp;
            tpm += sqt.substr(psm, atoi(dgt.c_str()));
            psm += atoi(dgt.c_str());
            bub += atoi(dgt.c_str());			
            break;
          case 'D':
            tpt += tar.substr(pst, atoi(dgt.c_str()));
            tpm += tmp;
            pst += atoi(dgt.c_str());
            bub += atoi(dgt.c_str());			
            break;
          case 'X':
            tpt += tar.substr(pst, atoi(dgt.c_str()));
            tpm += sqt.substr(psm, atoi(dgt.c_str()));			
            pst += atoi(dgt.c_str());
            psm += atoi(dgt.c_str());			
            bub += atoi(dgt.c_str());
            break;
          default:
            cerr << "Unkonwn align type: " << align.cigar_string[i] << endl;
            exit(0);
        }
        dgt = "";
      }
    }
    // Threshold 1, align details
    if (five + third > opt.bub || bub > opt.bub)
      continue;
    // Caculate the minimal free eneregy
    vector<int> pos = _mfe(tpm); 
    float mfe = 0.0;
    string fet(""), fem("");
    if (pos.size() >= 1) {
      for (int i = 0; i < pos.size() - 1; i += 2) {
        fet = tpt.substr(pos[i], pos[i+1]);
        fem = tpm.substr(pos[i], pos[i+1]);
        reverse(fem.begin(), fem.end());
        if (energy.find(fet+fem) != energy.end())
          mfe += energy[fet+fem];
        else if (energy.find(fet) != energy.end())
          mfe += energy[fet];
        else if (energy.find(fem) != energy.end())
          mfe += energy[fem];
      }
    }
    // Threshold 2, minimal free energy
    if (mfe > opt.mfe)
      continue;
    // Determine the align type of miRNA seed: 6-mer, 7-mer, 8-mer, [6-mer, 7-mer, 8-mer]_gu, none
    string typ("none");
    int m = 0, g = 0;
    if (clp <= 2 && tpm.length() > 8 - clp) {
      for (int i = 0; i < 8 - clp; i ++) {
        string bst = tpt.substr(tpt.length()-1-i, 1), bsm = tpm.substr(tpm.length()-1-i, 1);
        if (bst[0] == base[bsm[0]])
          m ++;
        else if ((bst.compare("G") == 0 && bsm.compare("T") == 0) || ((bst.compare("T") == 0 && bsm.compare("G") == 0)))
          g ++;
      }
      if (m + g >= 6) {
        stringstream ss;
        ss << m+g;
        typ = ss.str() + "-mer";
        if (g > 0)
          typ += "_gu";
      }
    } 
    // Threshold 3, align type of miRNA seed
    if (opt.strict == true && typ.compare("none") == 0)
      continue;
    pthread_mutex_lock(&mutex);
    int beg = align.ref_begin + iv->beg + 1;
    ofs << seq.nam << "\t" << iv->nam << "\t" << beg << "\t" << typ << "\t" << mfe << "\t" << align.sw_score << endl;
    pthread_mutex_unlock(&mutex);
  }
}

int _target()
{
  if (seq.seq.length() < opt.len)
    return 0;
  rec = _revc(seq.seq, 1); 
  _locate();
  int THREAD[opt.thd];
  pthread_t thread[opt.thd];
  for (int i = 0; i < opt.thd; i++) {
    THREAD[i] = i;
    int res = pthread_create(&thread[i], NULL, _score, (void *) &THREAD[i]);
    if (res != 0) {
      printf("Thread create %d failed\n", i);
      exit(0);
    }
  }
  while (1)
  {
    int mak = 0;
    void *status;
    for (int i = 0; i < opt.thd; i++)
      if (pthread_join(thread[i], &status) != 0)
        mak ++;
    if (mak == 0)
      break;
  }
  return 0;
}

