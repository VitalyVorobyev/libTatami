#ifndef RbkgPdf_H
#define RbkgPdf_H

#include "abspdf.h"
#include "icpvevent.h"

///
/// \brief The BkgPDFParSet class. Parametrs for background dt distribution
///

class BkgPDFParSet{
public:
  BkgPDFParSet();
  BkgPDFParSet(const BkgPDFParSet& x);

  BkgPDFParSet& operator=(const BkgPDFParSet& x);

  void SetTau(const double& x)         {m_tau         = x; return;}
  void Set_f_tail_mlt(const double& x) {m_f_tail_mlt  = x; return;}
  void Set_S_main_mlt(const double& x) {m_S_main_mlt  = x; return;}
  void Set_S_tail_mlt(const double& x) {m_S_tail_mlt  = x; return;}
  void Set_f_tail_sgl(const double& x) {m_f_tail_sgl  = x; return;}
  void Set_S_main_sgl(const double& x) {m_S_main_sgl  = x; return;}
  void Set_S_tail_sgl(const double& x) {m_S_tail_sgl  = x; return;}
  void Set_f_delta_mlt(const double& x){m_f_delta_mlt = x; return;}
  void Set_f_delta_sgl(const double& x){m_f_delta_sgl = x; return;}
  void Set_mu_delta(const double& x)   {m_mu_delta    = x; return;}
  void Set_mu(const double& x)         {m_mu          = x; return;}
  void Set_f_otlr(const double& x)     {m_f_otlr      = x; return;}
  void Set_s_otlr(const double& x)     {m_s_otlr      = x; return;}

  double tau(void)         const {return m_tau;}
  double f_tail_mlt(void)  const {return m_f_tail_mlt;}
  double S_main_mlt(void)  const {return m_S_main_mlt;}
  double S_tail_mlt(void)  const {return m_S_tail_mlt;}
  double f_tail_sgl(void)  const {return m_f_tail_sgl;}
  double S_main_sgl(void)  const {return m_S_main_sgl;}
  double S_tail_sgl(void)  const {return m_S_tail_sgl;}
  double f_delta_mlt(void) const {return m_f_delta_mlt;}
  double f_delta_sgl(void) const {return m_f_delta_sgl;}
  double mu_delta(void)    const {return m_mu_delta;}
  double mu(void)          const {return m_mu;}
  double f_otlr(void)      const {return m_f_otlr;}
  double s_otlr(void)      const {return m_s_otlr;}

  int GetParametersFromFile(const std::string& fname);

private:
  double m_S_main_mlt;
  double m_S_tail_mlt;
  double m_f_tail_mlt;
  double m_f_delta_mlt;

  double m_S_main_sgl;
  double m_S_tail_sgl;
  double m_f_tail_sgl;
  double m_f_delta_sgl;

  double m_mu;
  double m_mu_delta;
  double m_tau;

  double m_f_otlr;
  double m_s_otlr;
};

///
/// \brief The RbkgPdf class. Parameterization of background dt distribution
///

class RbkgPdf : public AbsPdf{
public:
  RbkgPdf(const std::string &fname);
  RbkgPdf(void);
  double operator()(const double& x);
  double operator()(const ICPVEvt& evt);
  BkgPDFParSet& Pars(void) {return m_pars;}

  void SetSigma(const double& x) {m_sigma = x;}
  void SetNDF(const int& x) {m_ndf = x;}
  void SetScale(const double& x){m_scale = x;}
  void SetShift(const double& x){m_shift = x;}
  void SetScaleShift(const double& scale,const double& shift){m_scale = scale; m_shift = shift;}
  double GetScale(void) const { return m_scale;}
  double GetShift(void) const { return m_shift;}

private:
  double Pdf(const double &x, const double &s, const int ndf);
  void DefaultInit(const int mode, const int svd, const bool mc, const bool tune);
  double AddOutlier(const double& x, const double Lin,const double& nLi);
  static const double cm2ps;

  BkgPDFParSet m_pars;
  double m_scale;
  double m_shift;

  /// Event-dependent estimation of uncertainty
  double m_sigma;
  /// Number of degrees of freedom
  int m_ndf;
};

#endif // RbkgPdf_H
