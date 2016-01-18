#include "icpvevent.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

ICPVEvt::ICPVEvt(const string &fname){
  ReadStructure(fname);
}

int ICPVEvt::ReadStructure(const string &fname){
  ifstream ifile(fname.c_str(),ofstream::in);
  if(!ifile.is_open()){
    cout << "ICPVEvt::ReadStructure: can't open file " << fname << endl;
    return -1;
  } else{
    cout << "Getting event structure from file " << fname << endl;
  }
  IVars.clear(); DVars.clear();
  string line;
  istringstream out;
  string name,type,dscr;
  while(getline(ifile,line)){
    out.str(line); out >> name; out >> type; out >> dscr;
         if(type == string("D")) DVars.push_back(DVar(0,name));
    else if(type == string("I")) IVars.push_back(IVar(0,name));
    else cout << "ICPVEvt::ReadStructure: wrong var type " << type << endl;
  }
  return 0;
}

void ICPVEvt::PrintStructure(void) const{
  cout << "Integer variables (" << IVars.size() << "):" << endl;
  for(unsigned i=0; i<IVars.size(); i++) cout << " " << IVars[i].Name();
  cout << endl;
  cout << "Real variables (" << DVars.size() << "):" << endl;
  for(unsigned i=0; i<DVars.size(); i++) cout << " " << DVars[i].Name();
  cout << endl;
  return;
}

int ICPVEvt::GetIVar(const int i) const{
  if(i<0 || i>=(int)IVars.size()){
    cout << "ICPVEvt::GetIVar: wrong var index " << i << endl;
    return 0;
  }
  return IVars[i].Val();
}

double ICPVEvt::GetDVar(const int i) const{
  if(i<0 || i>=(int)DVars.size()){
    cout << "ICPVEvt::GetDVar: wrong var index " << i << endl;
    return 0;
  }
  return DVars[i].Val();
}

int ICPVEvt::GetIVar(const string& name) const{
  return GetIVar(FindIVar(name));
}

double ICPVEvt::GetDVar(const string& name) const{
  return GetDVar(FindDVar(name));
}

int ICPVEvt::FindIVar(const string& name) const{
  const int index = FindIndex(IVars,name);
  if(index == -1) cout << "ICPVEvt::FindIVar: no variable with name " << name << endl;
  return index;
}

int ICPVEvt::FindDVar(const string& name) const{
  const int index = FindIndex(DVars,name);
  if(index == -1) cout << "ICPVEvt::FindDVar: no variable with name " << name << endl;
  return index;
}

int ICPVEvt::FindIndex(const vector<ICPVVar<auto> >& v, const string& name) const{
  for(unsigned i=0; i<v.size(); i++){ if(name == v[i].Name()) return i;}
  return -1;
}
