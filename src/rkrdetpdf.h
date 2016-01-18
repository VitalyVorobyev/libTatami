#ifndef RKRDETPDF_H
#define RKRDETPDF_H

#include "absicpvpdf.h"
#include "ResConst.h"
#include "ResVar.h"

///
/// \brief The RkRdetPdf class describes time resolution w/o non-primary vertices effect
///

class RkRdetPdf: public AbsICPVPdf{
public:
  RkRdetPdf(void): AbsICPV() {}
  double operator()(const ICPVEvt& evt);
  double operator()(const double& x);

  double      EfRkRdet_fullrec(const ICPVEvt& evt);
  double      AfRkRdet_fullrec(const ICPVEvt& evt);
  double      MfRkRdet_fullrec(const ICPVEvt& evt);
  double norm_EfRkRdet_fullrec(const ICPVEvt& evt);
  double norm_AfRkRdet_fullrec(const ICPVEvt& evt);
  double norm_MfRkRdet_fullrec(const ICPVEvt& evt);
private:
  double      EfRkRdet_fullrec(void) const;
  double      AfRkRdet_fullrec(void) const;
  double      MfRkRdet_fullrec(void) const;
  double norm_EfRkRdet_fullrec(void) const;
  double norm_AfRkRdet_fullrec(void) const;
  double norm_MfRkRdet_fullrec(void) const;

  double      EfRkRdet_full_sup(const double& mu, const double& sigma) const;
  double      AfRkRdet_full_sup(const double& mu, const double& sigma) const;
  double      MfRkRdet_full_sup(const double& mu, const double& sigma) const;
  double norm_EfRkRdet_full_sup(const double& mu, const double& sigma) const;
  double norm_AfRkRdet_full_sup(const double& mu, const double& sigma) const;
  double norm_MfRkRdet_full_sup(const double& mu, const double& sigma) const;

  double NoNPPdf(const ICPVEvt& evt, const bool otlr, const bool no_interf);
};

#endif // RKRDETPDF_H
