#ifndef ABSPDF_H
#define ABSPDF_H

#include "cnvl.h"
#include "icpvevent.h"

///
/// \brief The AbsPdf class. Abstract class for PDF based on cnvl library
///

class AbsPdf : public cnvl{
public:
  AbsPdf(): cnvl(), m_ll(-70), m_ul(70) {}
  /// Virtual method for calculation of PDF
  virtual double operator()(const ICPVEvt& ext) = 0;
  /// Virtual method for calculation of PDF
  virtual double operator()(const double& x) = 0;

  /// Lower limit of dt vatiable
  double ll(void) const {return m_ll;}
  /// Upper limit of dt vatiable
  double ul(void) const {return m_ul;}

  /// Set symmetrical range for dt variable
  void SetRange(const double& v){m_ll = v, m_ul = -v;}
  /// Set range for dt variable
  void SetRange(const double& min,const double& max){m_ll = max, m_ul = min;}

protected:
  double m_ll;
  double m_ul;
};

#endif // ABSPDF_H
