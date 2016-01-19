#ifndef ICPVEVENT_H
#define ICPVEVENT_H

#include <vector>
#include <string>

///
/// \brief The ICPVVar class
///

template<class T> class ICPVVar{
public:
  ICPVVar(const T& v, const std::string& n): val(v),name(n) {}
  void SetVal(const T& x) {val = x;}
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

typedef ICPVVar<double> dvar;
typedef ICPVVar<int>    ivar;

///
/// \brief The ICPVEvt class. Class for representation of tuple
/// tuple structure is initialized with text config file
///

class ICPVEvt{
public:
  ICPVEvt(const std::string& fname);
  ICPVEvt(const ICPVEvt& ev);
  ICPVEvt(const std::vector<ivar>& IV,const std::vector<dvar>& DV);

  ICPVEvt& operator=(const ICPVEvt& vt);

  void Set(const std::vector<ivar>& IV,const std::vector<dvar>& DV);
  void SetIVar(const int i,const int val);
  void SetIVar(const std::string& name,const int val);
  void SetDVar(const int i, const double& val);
  void SetDVar(const std::string& name, const double& val);

  std::string IName(const int i) const;
  int IVar(const int i) const;
  int IVar(const std::string& name) const;

  std::string DName(const int i) const;
  double DVar(const int i) const;
  double DVar(const std::string& name) const;

  int FindIVar(const std::string& name) const;
  int FindDVar(const std::string& name) const;

  const std::vector<ivar>& IVars(void) const { return m_IVars;}
  const std::vector<dvar>& DVars(void) const { return m_DVars;}

  void PrintStructure(void) const;
private:
  int FindIndex(const std::vector<ICPVVar<auto> > &v, const std::string& name) const;
  int ReadStructure(const std::string& fname);

  std::vector<ivar> m_IVars;
  std::vector<dvar> m_DVars;
};

#endif // ICPVEVENT_H
