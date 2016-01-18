#include "parmanager.h"

#include <sstream>

using namespace std;

string ParManager::WTagFile(const DataClass& dc){
  const string prefix("/home/vitaly/B0toD0pipi/libTatami/params/");
  stringstream out;
  out.str(""); out << prefix << "wtag_svd" << dc.Svd();
  if(dc.MC()) out << "_mc";
  out << ".txt";
  return out.str();
}

string ParManager::BkgParFile(const DataClass& dc){
  const string prefix("/home/vitaly/B0toD0pipi/libTatami/params/");
  return prefix + string("def_bkg.txt");
}

string ParManager::SigParFile(const DataClass& dc){
  const string prefix("/home/vitaly/B0toD0pipi/libTatami/params/");
  stringstream out;
  out.str(""); out << prefix << "sig_";
  out << (dc.Bp() ? "bp" : "b0") << "_svd" << dc.Svd();
  if(dc.MC()) out << "_mc";
  out << ".txt";
  return out.str();
}

