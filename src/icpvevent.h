#ifndef ICPVEVENT_H
#define ICPVEVENT_H

#include <vector>
#include <string>

template<class T> class ICPVVar{
public:
  ICPVVar(const T& v, const std::string& n): val(v),name(n) {}
  T Val(void)            const {return val;}
  std::string Name(void) const {return name;}

  ICPVVar& operator=(const ICPVVar& ovar){
    this->val  = ovar.val;
    this->name = ovar.name;
    return *this;
  }
private:
  T val;
  std::string name;
};

typedef ICPVVar<double> DVar;
typedef ICPVVar<int>    IVar;

/// \brief Class for representation of tuple
/// tuple structure is initialized with text config file

class ICPVEvt{
public:
  ICPVEvt(const std::string& fname);

  ICPVEvt& operator=(ICPVEvt& othevt);

  int GetIVar(const int i) const;
  int GetIVar(const std::string& name) const;

  double GetDVar(const int i) const;
  double GetDVar(const std::string& name) const;

  int FindIVar(const std::string& name) const;
  int FindDVar(const std::string& name) const;

  void PrintStructure(void) const;
private:
  int FindIndex(const std::vector<ICPVVar<auto> > &v, const std::string& name) const;
  int ReadStructure(const std::string& fname);
  std::vector<IVar> IVars;
  std::vector<DVar> DVars;
};

#endif // ICPVEVENT_H
