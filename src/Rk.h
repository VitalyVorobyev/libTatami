#ifndef RK_H
#define RK_H

#include "absicpvpdf.h"
#include "rkparam.h"

///
/// \brief The RkPdf class. Describes time resolution part related to kinematic approximation dz -> dt.
///

class RkPdf: public AbsICPVPdf{
public:
  RkPdf(): AbsICPVPdf() {}
  /// Calculate PDF
  double operator()(const ICPVEvt& evt);
  /// Calculate PDF
  double operator()(const double& x) const;

private:
  double Pdf(const double& x) const;
  double norm_EfRk(void) const;
  double norm_AfRk(void) const;
  double norm_MfRk(void) const;
  double EfRk(const double& x) const;
  double AfRk(const double& x) const;
  double MfRk(const double& x) const;

  RkPar m_pars;  
};

#endif // RK_H
