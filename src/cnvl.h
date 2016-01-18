#ifndef CNVL_H
#define CNVL_H

#include <iostream>
#include <cmath>
#include <float.h>

#include <complex>

typedef double (*nXi_t)(const double& t, const double& xd);

class cnvl{
public:
  cnvl() {}
  ~cnvl() {}//delete roomath;}

//  int SetKKCS(const double& _K, const double& _Kb, const double& _C, const double& _S);
//  int SetFlvXi(const int _flv, const int _xi);
//  int SetXi(const int _xi);
//  int SetSinCos(const double& sin, const double& cos);

//  double GetSin2Beta(void) const {return sin2beta;}
//  double GetCos2Beta(void) const {return cos2beta;}

protected:
  double recexp(const double& re, const double& im) const;
  double imcexp(const double& re, const double& im) const;
  double rewerf(const double& re, const double& im) const;
  double imwerf(const double& re, const double& im) const;

  double Ep(const double& t, const double& tau) const;
  double En(const double& t, const double& tau) const;
  double Ef(const double& t, const double& tau) const;
  double Enp(const double& t, const double& tau_n, const double& tau_p) const;
  double xEp(const double& t, const double& tau) const;
  double xEn(const double& t, const double& tau) const;
  double xEf(const double& t, const double& tau) const;
  static double nMp(const double& t, const double& xd);
  double Mp(const double& t, const double& tau, const double& dm) const;
  static double nMn(const double& t, const double& xd);
  double Mn(const double& t, const double& tau, const double& dm) const;
  double nMf(const double& t, const double& xd) const;
  double Mf(const double& t, const double& tau, const double& dm) const;
  static double nAp(const double& t, const double& xd);
  double Ap(const double& t, const double& tau, const double& dm) const;
  static double nAn(const double& t, const double& xd);
  double An(const double& t, const double& tau, const double& dm) const;
  double nAf(const double& t, const double& xd) const;
  double Af(const double& t, const double& tau, const double& dm) const;

  double norm_nEp(const double& _ll, const double& _ul, const double& o = 0) const;
  double norm_nEn(const double& _ll, const double& _ul, const double& o = 0) const;
  double norm_nEf(const double& _ll, const double& _ul, const double& o = 0) const;
  double norm_Ep(const double& _ll, const double& _ul, const double& tau, const double& o = 0) const;
  double norm_En(const double& _ll, const double& _ul, const double& tau, const double& o = 0) const;
  double norm_Ef(const double& _ll, const double& _ul, const double& tau, const double& o = 0) const;
  double norm_Ap(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0) const;
  double norm_Ax_sup(const double& x, const double& tau, const double& dm) const;
  double norm_Mp(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0) const;
  double norm_Mn(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0) const;
  double norm_Mx_sup(const double& x, const double& tau, const double& dm) const;
  double norm_An(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0) const;

  double nEp_conv_gauss(const double& t, const double& m, const double& s) const;
  double nEn_conv_gauss(const double& t, const double& m, const double& s) const;
  double Ep_conv_gauss(const double& t, const double& tau, const double& m, const double& s) const;
  double En_conv_gauss(const double& t, const double& tau, const double& m, const double& s) const;
  double Ef_conv_gauss(const double& t, const double& tau, const double& m, const double& s) const;
  double Enp_conv_gauss(const double& t, const double& tau_n, const double& tau_p, const double& m, const double& s) const;
  double approx_exp2erfc(const double& x) const;

  double nMp_conv_gauss(const double& t, const double& xd, const double& m, const double& s) const;
  double nMn_conv_gauss(const double& t, const double& xd, const double& m, const double& s) const;
  double nAp_conv_gauss(const double& t, const double& xd, const double& m, const double& s) const;
  double nAn_conv_gauss(const double& t, const double& xd, const double& m, const double& s) const;
  double Mp_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s) const;
  double Mn_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s) const;
  double Mf_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s) const;
  double Ap_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s) const;
  double An_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s) const;
  double Af_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s) const;

  double _IM(const double& x, const double& m, const double& s, const double& beta, const double& gamma, const double& b) const;
  double _IA(const double& x, const double& m, const double& s, const double& beta, const double& gamma, const double& b) const;

  double gaussian(const double& x, const double& m, const double& s) const;
  double norm_gaussian_w_cutoff(const double& cutoff, const double& m, const double& s) const;
  double norm_gaussian(const double& ll, const double& ul, const double& m, const double& s) const;
  double DiracDelta(const double& x) const;

  double xXi_conv_gauss_by_int(nXi_t p_func, const double& t, const double& xd, const double& mu, const double& sigma) const;
  double xEp_conv_gauss(const double& t, const double& tau, const double& m, const double& s) const;
  double xEn_conv_gauss(const double& t, const double& tau, const double& m, const double& s) const;
  double xEf_conv_gauss(const double& t, const double& tau, const double& m, const double& s) const;

  double norm_neg_sup(const double& m, const double& s) const;
  double norm_nEp_conv_gauss_sub(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.) const;
  double norm_nEn_conv_gauss_sub(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.) const;
  double norm_nEp_conv_gauss(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.) const;
  double norm_nEn_conv_gauss(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.) const;
  double norm_Ep_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.) const;
  double norm_En_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.) const;
  double norm_nEf_conv_gauss(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.) const;
  double norm_Ef_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.) const;
  double norm_Enp_conv_gauss(const double& _ll, const double& _ul,const double& tau_n, const double& tau_p,const double& m, const double& s, const double& o = 0.) const;
  double norm_nag_sup(const double& m, const double& s, const double& xd) const;
  double norm_nmg_sup(const double& m, const double& s, const double& xd) const;
  double norm_nAn_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.) const;
  double norm_nAp_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.) const;
  double norm_nAf_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.) const;
  double norm_nMn_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.) const;
  double norm_nMp_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.) const;
  double norm_nMf_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.) const;
  double norm_An_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.) const;
  double norm_Ap_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.) const;
  double norm_Af_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.) const;
  double norm_Mn_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.) const;
  double norm_Mp_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.) const;
  double norm_Mf_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.) const;

  double int_polyexp2(const double& _ll, const double& _ul,const double& alpha, const double& beta, const double& gamma,const double& a) const;
  double int_polyexp_erfc(const double& _ll, const double& _ul,const double& alpha, const double& beta, const double& gamma,const double& a) const;
  double norm_xEp_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.) const;
  double norm_xEn_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.) const;
  double norm_xEf_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.) const;

  double norm_xEp(const double& _ll, const double& _ul, const double& tau, const double& o = 0.) const;
  double norm_xEn(const double& _ll, const double& _ul, const double& tau, const double& o = 0.) const;
  double norm_xEf(const double& _ll, const double& _ul, const double& tau, const double& o = 0.) const;

protected:
  // ** Interface parameters for inherits classes **
//  double K,Kb,C,S;
//  int flv, xi;
//  double sin2beta;
//  double cos2beta;

  // ** Constants **
  static const double inv_sqrt2;
  static const double inv_sqrt_pi;
  static const double sqrt_pi;
  static const double inv_sqrt_2pi;
  static const double sqrt2;
};

#endif // CNVL_H
