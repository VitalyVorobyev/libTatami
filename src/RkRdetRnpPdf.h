#ifndef RKRDETRNPRDF_H
#define RKRDETRNPRDF_H

#include "conv_coef.h"
#include "ResConst.h"
#include "ResVar.h"
#include "parmanager.h"
#include "wtag.h"
#include "icpvevent.h"
#include "ttools.h"
#include "rrecpars.h"
#include "rascrnppars.h"
#include "rkparam.h"
#include "absicpvpdf.h"

///
/// \brief The RkRdetRnpPdf class. Provides full parameterization of signal dt distribution at Belle
///

class RkRdetRnpPdf : public AbsICPVPdf{
public:
  RkRdetRnpPdf(const DataClass &dc);

  /// Calculate PDF
  double operator()(const ICPVEvt& ext);
  /// Calculate PDF
  double operator()(const double& x);

  /// Unnormalized convolution Exp x Gauss
  double      EfRkRdetRnp_fullrec(const ICPVEvt& evt);
  /// Unnormalized convolution Exp*cos x Gauss
  double      MfRkRdetRnp_fullrec(const ICPVEvt& evt);
  /// Unnormalized convolution Exp*sin x Gauss
  double      AfRkRdetRnp_fullrec(const ICPVEvt& evt);
  /// Normalization for convolution Exp x Gauss
  double norm_EfRkRdetRnp_fullrec(const ICPVEvt& evt);
  /// Normalization for convolution Exp*cos x Gauss
  double norm_MfRkRdetRnp_fullrec(const ICPVEvt& evt);
  /// Normalization for convolution Exp*sin x Gauss
  double norm_AfRkRdetRnp_fullrec(const ICPVEvt& evt);

private:
  int ReadVars(const ICPVEvt& evt);
  int ReadAndCalc(const ICPVEvt& evt);

  double PdfAB(  const ICPVEvt& evt, const bool otlr = true, const bool no_interf = false);
  double PdfAB(const bool otlr = true, const bool no_interf = false) const;

  double      EfRkRdetRnp_fullrec(void) const;
  double      MfRkRdetRnp_fullrec(void) const;
  double      AfRkRdetRnp_fullrec(void) const;
  double norm_EfRkRdetRnp_fullrec(void) const;
  double norm_MfRkRdetRnp_fullrec(void) const;
  double norm_AfRkRdetRnp_fullrec(void) const;

  double      EfRkRdetRnp_full_sup(const double& mu, const double& sigma) const;
  double      MfRkRdetRnp_full_sup(const double& mu, const double& sigma) const;
  double      AfRkRdetRnp_full_sup(const double& mu, const double& sigma) const;
  double norm_EfRkRdetRnp_full_sup(const double& mu, const double& sigma) const;
  double norm_MfRkRdetRnp_full_sup(const double& mu, const double& sigma) const;
  double norm_AfRkRdetRnp_full_sup(const double& mu, const double& sigma) const;
//  double AddOutlier( const double& x, const double& Lin, const double& nLi = 1.0, const double& alpha = 1.0);
//  double Add_Outlier(const double& x, const double& Lin, const double& nLi, const double& alpha);
private:
  int flavor;
  int keeptagl;
  int m_svd;
  double cexp, amix;
  WTag m_wtag;

  double dt;
  RdetVar m_rvar;
  RdetVar m_avar;
  double costhBcms;

  RascRnpPars m_apar;
  RrecPars    m_rpar;
  RkPar       m_kpar;
  ResConst    m_cnst;

  conv_coef coco;
};

#endif // RKRDETRNPRDF_H
