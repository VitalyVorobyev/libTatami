#ifndef PARMANAGER_H
#define PARMANAGER_H

#include <string>

class DataClass{
public:
  DataClass(const int svd, const bool mc, const bool bp): m_svd(svd),m_mc(mc),m_bp(bp) {}
  DataClass(const int svd, const bool mc): DataClass(svd,mc,false) {}
  DataClass(const int svd): DataClass(svd,false,false) {}
  DataClass(void): DataClass(0,false,false) {}

  void SetSVD(const int v) {m_svd = v;}
  void SetMC(const bool v) {m_mc  = v;}
  void SetBp(const bool v) {m_bp  = v;}

  int Svd(void) const {return m_svd;}
  bool MC(void) const {return m_mc;}
  bool Bp(void) const {return m_bp;}

private:
  int m_svd;
  bool m_mc;
  bool m_bp;
};

class ParManager{
public:
  ParManager() {}

  static std::string WTagFile(const DataClass& dc);
  static std::string BkgParFile(const DataClass& dc);
  static std::string SigParFile(const DataClass& dc);
};

#endif // PARMANAGER_H
