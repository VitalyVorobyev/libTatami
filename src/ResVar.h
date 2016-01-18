#ifndef RESVAR_H
#define RESVAR_H

#include "icpvevent.h"

///
/// \brief The RdetVar class. Event-dependent variables for vertex resilution
///

class RdetVar{
public:
  RdetVar(void);
  RdetVar(const RdetVar& var);

  RdetVar& operator=(const RdetVar& var);
  int ReadVars(const ICPVEvt &evt, const bool type);

  void Set_ntrk(const int v)      {m_ntrk = v;}
  void Set_sz(const double& v)    {m_sz = v;}
  void Set_chisq(const double& v) {m_chisq = v;}
  void Set_ndf(const int v)       {m_ndf = v;}

  int ntrk(void)     const {return m_ntrk;}
  double sz(void)    const {return m_sz;}
  double chisq(void) const {return m_chisq;}
  int ndf(void)      const {return m_ndf;}

  static const bool RecSide;
  static const bool AscSide;

private:
  int m_ntrk;
  double m_sz;
  double m_chisq;
  int m_ndf;
};

#endif
