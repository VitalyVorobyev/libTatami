/** Copyright 2016 Vitaly Vorobyev
 ** @file cnvl.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_CNVL_H_
#define INCLUDE_CNVL_H_

#include <iostream>
#include <cmath>
#include <cfloat>

namespace libTatami {

typedef double (*nXi_t)(const double& t, const double& xd);

///
/// \brief The cnvl class
///
class cnvl {
 public:
    ///
    /// \brief cnvl
    ///
    cnvl() {}
    ~cnvl() {}

 protected:
    ///
    /// \brief recexp. Real part of exp(re + i*im)
    /// \param re. Real part of argument
    /// \param im. Imag part of argument
    /// \return
    ///
    double recexp(const double& re, const double& im) const;
    ///
    /// \brief imcexp. Imag part of exp(re + i*im)
    /// \param re. Real part of the argument
    /// \param im. Imag part of the argument
    /// \return
    ///
    double imcexp(const double& re, const double& im) const;
    ///
    /// \brief rewerf. Real part of exp(-z^2) erfc(-iz),
    /// where z = re + i*im
    /// \param re. Real part of the argument
    /// \param im. Imag part of the argument
    /// \return
    ///
    double rewerf(const double& re, const double& im) const;
    ///
    /// \brief imwerf. Imag part of exp(-z^2) erfc(-iz),
    /// where z = re + i*im
    /// \param re. Real part of the argument
    /// \param im. Imag part of the argument
    /// \return
    ///
    double imwerf(const double& re, const double& im) const;
    ///
    /// \brief Ep. Normalized exponential distribution for $t > 0$
    /// and $|\tau|$.
    /// \param t. Argument of the exponential distribution
    /// \param $\tau$. Parameter of the exponential distribution
    /// \return Exponential PDF for t >= 0. Returns zero if t < 0
    ///
    double Ep(const double& t, const double& tau) const;
    ///
    /// \brief En. Normalized exponential distribution for t < 0
    /// and |tau|.
    /// \param t. Argument of the exponential distribution
    /// \param tau. Parameter of the exponential distribution
    /// \return Exponential PDF for t < 0. Returns zero if t >= 0
    ///
    double En(const double& t, const double& tau) const;
    ///
    /// \brief Ef. Symmetrized normalized exponential distribution for |tau|
    /// \param t. Argument of the exponential distribution
    /// \param tau. Parameter of the exponential distribution
    /// \return PDF
    ///
    double Ef(const double& t, const double& tau) const;
    ///
    /// \brief Enp. Bifurcated normalized exponential distribution.
    /// The parameter |tau_n| is used for negative argument and
    /// the parameter |tau_p| is used for positive argument.
    /// \param t. Argument of the exponential distribution
    /// \param tau_n. Parameter for negative arguments
    /// \param tau_p. Parameter for positive arguments
    /// \return PDF
    ///
    double Enp(const double& t, const double& tau_n,
               const double& tau_p) const;
    ///
    /// \brief xEp. t * Ep(t, tau) / |tau|
    /// \param t. Argument
    /// \param tau. Parameter
    /// \return t * Ep(t, tau) / |tau|
    ///
    double xEp(const double& t, const double& tau) const;
    ///
    /// \brief xEn. -t * En(t, tau) / |tau|
    /// \param t. Argument
    /// \param tau. Parameter
    /// \return
    ///
    double xEn(const double& t, const double& tau) const;
    ///
    /// \brief xEf. |t| * Ef(t, tau) / |tau|
    /// \param t. Argument
    /// \param tau. Parameter
    /// \return
    ///
    double xEf(const double& t, const double& tau) const;
    ///
    /// \brief nMp. e^{-|t|} * cos(xd * t)
    /// \param t
    /// \param xd
    /// \return Returns 0 if t < 0
    ///
    static double nMp(const double& t, const double& xd);
    ///
    /// \brief Mp. nMp(t/|tau|, dm * |tau|), t >= 0
    /// \param t
    /// \param tau
    /// \param dm
    /// \return Returns 0 if t < 0
    ///
    double Mp(const double& t, const double& tau, const double& dm) const;
    ///
    /// \brief nMn. e^{-|t|} * cos(xd * t), t < 0
    /// \param t
    /// \param xd
    /// \return Returns 0 if t >= 0
    ///
    static double nMn(const double& t, const double& xd);
    ///
    /// \brief Mn. nMn(t/|tau|, dm * |tau|), t < 0
    /// \param t
    /// \param tau
    /// \param dm
    /// \return Returns 0 if t >= 0
    ///
    double Mn(const double& t, const double& tau, const double& dm) const;
    ///
    /// \brief nMf. e^(-|t|) * cos(xd * t)
    /// \param t
    /// \param xd
    /// \return
    ///
    double nMf(const double& t, const double& xd) const;
    ///
    /// \brief Mf. nMf(t/|tau|, dm * |tau|)
    /// \param t
    /// \param tau
    /// \param dm
    /// \return
    ///
    double Mf(const double& t, const double& tau, const double& dm) const;
    ///
    /// \brief nAp. e^{-|t|} * sin(xd * t), t >= 0
    /// \param t
    /// \param xd
    /// \return Returns 0 if t < 0
    ///
    static double nAp(const double& t, const double& xd);
    ///
    /// \brief Ap. nAp(t/|tau|, dm * |tau|)
    /// \param t
    /// \param tau
    /// \param dm
    /// \return Returns 0 if t < 0
    ///
    double Ap(const double& t, const double& tau, const double& dm) const;
    ///
    /// \brief nAn. e^{-|t|} * sin(xd * t), t < 0
    /// \param t
    /// \param xd
    /// \return Returns 0 if t >= 0
    ///
    static double nAn(const double& t, const double& xd);
    ///
    /// \brief An. nAn(t/|tau|, dm * |tau|)
    /// \param t
    /// \param tau
    /// \param dm
    /// \return Returns 0 if t >= 0
    ///
    double An(const double& t, const double& tau, const double& dm) const;
    ///
    /// \brief nAf. e^(-|t|) * sin(xd * t)
    /// \param t
    /// \param xd
    /// \return
    ///
    double nAf(const double& t, const double& xd) const;
    ///
    /// \brief Af. nAf(t/|tau|, dm * |tau|)
    /// \param t. Argument
    /// \param tau. The lifetime parameter
    /// \param dm. The mass difference parameter
    /// \return
    ///
    double Af(const double& t, const double& tau, const double& dm) const;
    ///
    /// \brief norm_nEp. Normalization coefficient for the nEp function.
    /// Normalization interval is defined as (ll, ul), where
    /// ll = max(_ll - o, 0) and ul = max(_ul - o, 0)
    /// \param _ll. Lower limit
    /// \param _ul. Upper limit
    /// \param o. Offset
    /// \return e^{-ll} - e^{-ul}
    ///
    double norm_nEp(const double& _ll, const double& _ul,
                    const double& o = 0) const;
    ///
    /// \brief norm_Ap. Normalization coefficient for the Ap function.
    /// Normalization interval is defined as (ll, ul), where
    /// ll = min(_ll, _ul) - o and
    /// ul = max(_ll, _ul) - o
    /// \param _ll.
    /// \param _ul
    /// \param tau
    /// \param dm
    /// \param o. Offset.
    /// \return
    ///
    double norm_Ap(const double& _ll, const double& _ul, const double& tau,
                   const double& dm, const double& o = 0) const;
    ///
    /// \brief norm_nEn. Normalization coefficient for the nEn function.
    /// Normalization interval is defined as (ll, ul), where
    /// ll = min(_ll - o, 0) and ul = min(_ul - o, 0)
    /// \param _ll. Lower limit
    /// \param _ul. Upper limit
    /// \param o. Offset
    /// \return e^{ll} - e^{ul}
    ///
    double norm_nEn(const double& _ll, const double& _ul,
                    const double& o = 0) const;
    ///
    /// \brief norm_nEf. Normalization coefficient for the nEf function.
    /// \param _ll. Lower limit
    /// \param _ul. Upper limit
    /// \param o. Offset
    /// \return 0.5 * (norm_nEn() + norm_nEp)
    ///
    double norm_nEf(const double& _ll, const double& _ul,
                    const double& o = 0) const;
    ///
    /// \brief norm_Ep. Normalization coefficient for the Ep function.
    /// \param _ll. Lower limit
    /// \param _ul. Upper limit
    /// \param tau. Lifetime
    /// \param o. Offset
    /// \return norm_nEp(ll, ul, no), where ll = _ll / |tau|,
    /// ul = _ul / |tau| and no = o / |tau|
    ///
    double norm_Ep(const double& _ll, const double& _ul,
                   const double& tau, const double& o = 0) const;
    ///
    /// \brief norm_En. Normalization coefficient for the En function.
    /// \param _ll. Lower limit
    /// \param _ul. Upper limit
    /// \param tau. Lifetime
    /// \param o. Offset
    /// \return norm_nEn(ll, ul, no), where ll = _ll / |tau|,
    /// ul = _ul / |tau| and no = o / |tau|
    ///
    double norm_En(const double& _ll, const double& _ul,
                   const double& tau, const double& o = 0) const;
    ///
    /// \brief norm_Ef. Normalization coefficient for the Ef function.
    /// \param _ll. Lower limit
    /// \param _ul. Upper limit
    /// \param tau. Lifetime
    /// \param o. Offset
    /// \return norm_nEf(ll, ul, no), where ll = _ll / |tau|,
    /// ul = _ul / |tau| and no = o / |tau|
    ///
    double norm_Ef(const double& _ll, const double& _ul,
                   const double& tau, const double& o = 0) const;
    ///
    /// \brief norm_Ax_sup. Supplementary routine for normalization
    /// of the Ap and An functions
    /// \param x
    /// \param tau
    /// \param dm
    /// \return
    ///
    double norm_Ax_sup(const double& x, const double& tau,
                       const double& dm) const;
    ///
    /// \brief norm_Mp
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param dm
    /// \param o
    /// \return
    ///
    double norm_Mp(const double& _ll, const double& _ul, const double& tau,
                   const double& dm, const double& o = 0) const;
    ///
    /// \brief norm_Mn
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param dm
    /// \param o
    /// \return
    ///
    double norm_Mn(const double& _ll, const double& _ul, const double& tau,
                   const double& dm, const double& o = 0) const;
    ///
    /// \brief norm_Mx_sup
    /// \param x
    /// \param tau
    /// \param dm
    /// \return
    ///
    double norm_Mx_sup(const double& x, const double& tau,
                       const double& dm) const;
    ///
    /// \brief norm_An
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param dm
    /// \param o
    /// \return
    ///
    double norm_An(const double& _ll, const double& _ul, const double& tau,
                   const double& dm, const double& o = 0) const;
    ///
    /// \brief nEp_conv_gauss
    /// \param t
    /// \param m
    /// \param s
    /// \return
    ///
    double nEp_conv_gauss(const double& t, const double& m,
                          const double& s) const;
    ///
    /// \brief nEn_conv_gauss
    /// \param t
    /// \param m
    /// \param s
    /// \return
    ///
    double nEn_conv_gauss(const double& t, const double& m,
                          const double& s) const;
    ///
    /// \brief Ep_conv_gauss
    /// \param t
    /// \param tau
    /// \param m
    /// \param s
    /// \return
    ///
    double Ep_conv_gauss(const double& t, const double& tau,
                         const double& m, const double& s) const;
    ///
    /// \brief En_conv_gauss
    /// \param t
    /// \param tau
    /// \param m
    /// \param s
    /// \return
    ///
    double En_conv_gauss(const double& t, const double& tau,
                         const double& m, const double& s) const;
    ///
    /// \brief Ef_conv_gauss
    /// \param t
    /// \param tau
    /// \param m
    /// \param s
    /// \return
    ///
    double Ef_conv_gauss(const double& t, const double& tau,
                         const double& m, const double& s) const;
    ///
    /// \brief Enp_conv_gauss
    /// \param t
    /// \param tau_n
    /// \param tau_p
    /// \param m
    /// \param s
    /// \return
    ///
    double Enp_conv_gauss(const double& t, const double& tau_n,
                          const double& tau_p, const double& m,
                          const double& s) const;
    ///
    /// \brief approx_exp2erfc. Tailor expansion of ???
    /// \param x. Argument
    /// \return (1/x - 1/(2*x^3)) / sqrt(pi)
    ///
    double approx_exp2erfc(const double& x) const;
    ///
    /// \brief nMp_conv_gauss
    /// \param t
    /// \param xd
    /// \param m
    /// \param s
    /// \return
    ///
    double nMp_conv_gauss(const double& t, const double& xd, const double& m,
                          const double& s) const;
    ///
    /// \brief nMn_conv_gauss
    /// \param t
    /// \param xd
    /// \param m
    /// \param s
    /// \return
    ///
    double nMn_conv_gauss(const double& t, const double& xd, const double& m,
                          const double& s) const;
    ///
    /// \brief nAp_conv_gauss
    /// \param t
    /// \param xd
    /// \param m
    /// \param s
    /// \return
    ///
    double nAp_conv_gauss(const double& t, const double& xd, const double& m,
                          const double& s) const;
    ///
    /// \brief nAn_conv_gauss
    /// \param t
    /// \param xd
    /// \param m
    /// \param s
    /// \return
    ///
    double nAn_conv_gauss(const double& t, const double& xd, const double& m,
                          const double& s) const;
    ///
    /// \brief Mp_conv_gauss
    /// \param t
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \return
    ///
    double Mp_conv_gauss(const double& t, const double& tau, const double& dm,
                         const double& m, const double& s) const;
    ///
    /// \brief Mn_conv_gauss
    /// \param t
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \return
    ///
    double Mn_conv_gauss(const double& t, const double& tau, const double& dm,
                         const double& m, const double& s) const;
    ///
    /// \brief Mf_conv_gauss
    /// \param t
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \return
    ///
    double Mf_conv_gauss(const double& t, const double& tau, const double& dm,
                         const double& m, const double& s) const;
    ///
    /// \brief Ap_conv_gauss
    /// \param t
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \return
    ///
    double Ap_conv_gauss(const double& t, const double& tau, const double& dm,
                         const double& m, const double& s) const;
    ///
    /// \brief An_conv_gauss
    /// \param t
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \return
    ///
    double An_conv_gauss(const double& t, const double& tau, const double& dm,
                         const double& m, const double& s) const;
    ///
    /// \brief Af_conv_gauss
    /// \param t
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \return
    ///
    double Af_conv_gauss(const double& t, const double& tau, const double& dm,
                         const double& m, const double& s) const;
    ///
    /// \brief _IM
    /// \param x
    /// \param m
    /// \param s
    /// \param beta
    /// \param gamma
    /// \param b
    /// \return
    ///
    double _IM(const double& x, const double& m, const double& s,
               const double& beta, const double& gamma, const double& b) const;
    ///
    /// \brief _IA
    /// \param x
    /// \param m
    /// \param s
    /// \param beta
    /// \param gamma
    /// \param b
    /// \return
    ///
    double _IA(const double& x, const double& m, const double& s,
               const double& beta, const double& gamma, const double& b) const;
    ///
    /// \brief gaussian
    /// \param x
    /// \param m
    /// \param s
    /// \return
    ///
    double gaussian(const double& x, const double& m, const double& s) const;
    ///
    /// \brief norm_gaussian_w_cutoff
    /// \param cutoff
    /// \param m
    /// \param s
    /// \return
    ///
    double norm_gaussian_w_cutoff(const double& cutoff, const double& m,
                                  const double& s) const;
    ///
    /// \brief norm_gaussian
    /// \param ll
    /// \param ul
    /// \param m
    /// \param s
    /// \return
    ///
    double norm_gaussian(const double& ll, const double& ul, const double& m,
                         const double& s) const;
    ///
    /// \brief DiracDelta
    /// \param x
    /// \return
    ///
    double DiracDelta(const double& x) const;
    ///
    /// \brief xXi_conv_gauss_by_int
    /// \param p_func
    /// \param t
    /// \param xd
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double xXi_conv_gauss_by_int(nXi_t p_func, const double& t,
                                 const double& xd, const double& mu,
                                 const double& sigma) const;
    ///
    /// \brief xEp_conv_gauss
    /// \param t
    /// \param tau
    /// \param m
    /// \param s
    /// \return
    ///
    double xEp_conv_gauss(const double& t, const double& tau, const double& m,
                          const double& s) const;
    ///
    /// \brief xEn_conv_gauss
    /// \param t
    /// \param tau
    /// \param m
    /// \param s
    /// \return
    ///
    double xEn_conv_gauss(const double& t, const double& tau, const double& m,
                          const double& s) const;
    ///
    /// \brief xEf_conv_gauss
    /// \param t
    /// \param tau
    /// \param m
    /// \param s
    /// \return
    ///
    double xEf_conv_gauss(const double& t, const double& tau, const double& m,
                          const double& s) const;
    ///
    /// \brief norm_neg_sup
    /// \param m
    /// \param s
    /// \return
    ///
    double norm_neg_sup(const double& m, const double& s) const;
    ///
    /// \brief norm_nEp_conv_gauss_sub
    /// \param _ll
    /// \param _ul
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nEp_conv_gauss_sub(const double& _ll, const double& _ul,
                                   const double& m, const double& s,
                                   const double& o = 0.) const;
    ///
    /// \brief norm_nEn_conv_gauss_sub
    /// \param _ll
    /// \param _ul
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nEn_conv_gauss_sub(const double& _ll, const double& _ul,
                                   const double& m, const double& s,
                                   const double& o = 0.) const;
    ///
    /// \brief norm_nEp_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nEp_conv_gauss(const double& _ll, const double& _ul,
                               const double& m, const double& s,
                               const double& o = 0.) const;
    ///
    /// \brief norm_nEn_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nEn_conv_gauss(const double& _ll, const double& _ul,
                               const double& m, const double& s,
                               const double& o = 0.) const;
    ///
    /// \brief norm_Ep_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Ep_conv_gauss(const double& _ll, const double& _ul,
                              const double& tau, const double& m,
                              const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_En_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_En_conv_gauss(const double& _ll, const double& _ul,
                              const double& tau, const double& m,
                              const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_nEf_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nEf_conv_gauss(const double& _ll, const double& _ul,
                               const double& m, const double& s,
                               const double& o = 0.) const;
    ///
    /// \brief norm_Ef_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Ef_conv_gauss(const double& _ll, const double& _ul,
                              const double& tau, const double& m,
                              const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_Enp_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param tau_n
    /// \param tau_p
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Enp_conv_gauss(const double& _ll, const double& _ul,
                               const double& tau_n, const double& tau_p,
                               const double& m, const double& s,
                               const double& o = 0.) const;
    ///
    /// \brief norm_nag_sup
    /// \param m
    /// \param s
    /// \param xd
    /// \return
    ///
    double norm_nag_sup(const double& m, const double& s,
                        const double& xd) const;
    ///
    /// \brief norm_nmg_sup
    /// \param m
    /// \param s
    /// \param xd
    /// \return
    ///
    double norm_nmg_sup(const double& m, const double& s,
                        const double& xd) const;
    ///
    /// \brief norm_nAn_conv_gauss
    /// \param ll
    /// \param ul
    /// \param xd
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nAn_conv_gauss(const double& ll, const double& ul,
                               const double& xd, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_nAp_conv_gauss
    /// \param ll
    /// \param ul
    /// \param xd
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nAp_conv_gauss(const double& ll, const double& ul,
                               const double& xd, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_nAf_conv_gauss
    /// \param ll
    /// \param ul
    /// \param xd
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nAf_conv_gauss(const double& ll, const double& ul,
                               const double& xd, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_nMn_conv_gauss
    /// \param ll
    /// \param ul
    /// \param xd
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nMn_conv_gauss(const double& ll, const double& ul,
                               const double& xd, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_nMp_conv_gauss
    /// \param ll
    /// \param ul
    /// \param xd
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nMp_conv_gauss(const double& ll, const double& ul,
                               const double& xd, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_nMf_conv_gauss
    /// \param ll
    /// \param ul
    /// \param xd
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_nMf_conv_gauss(const double& ll, const double& ul,
                               const double& xd, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_An_conv_gauss
    /// \param ll
    /// \param ul
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_An_conv_gauss(const double& ll, const double& ul,
                              const double& tau, const double& dm,
                              const double& m, const double& s,
                              const double& o = 0.) const;
    ///
    /// \brief norm_Ap_conv_gauss
    /// \param ll
    /// \param ul
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Ap_conv_gauss(const double& ll, const double& ul,
                              const double& tau, const double& dm,
                              const double& m, const double& s,
                              const double& o = 0.) const;
    ///
    /// \brief norm_Af_conv_gauss
    /// \param ll
    /// \param ul
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Af_conv_gauss(const double& ll, const double& ul,
                              const double& tau, const double& dm,
                              const double& m, const double& s,
                              const double& o = 0.) const;
    ///
    /// \brief norm_Mn_conv_gauss
    /// \param ll
    /// \param ul
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Mn_conv_gauss(const double& ll, const double& ul,
                              const double& tau, const double& dm,
                              const double& m, const double& s,
                              const double& o = 0.) const;
    ///
    /// \brief norm_Mp_conv_gauss
    /// \param ll
    /// \param ul
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Mp_conv_gauss(const double& ll, const double& ul,
                              const double& tau, const double& dm,
                              const double& m, const double& s,
                              const double& o = 0.) const;
    ///
    /// \brief norm_Mf_conv_gauss
    /// \param ll
    /// \param ul
    /// \param tau
    /// \param dm
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_Mf_conv_gauss(const double& ll, const double& ul,
                              const double& tau, const double& dm,
                              const double& m, const double& s,
                              const double& o = 0.) const;
    ///
    /// \brief int_polyexp2
    /// \param _ll
    /// \param _ul
    /// \param alpha
    /// \param beta
    /// \param gamma
    /// \param a
    /// \return
    ///
    double int_polyexp2(const double& _ll, const double& _ul,
                        const double& alpha, const double& beta,
                        const double& gamma, const double& a) const;
    ///
    /// \brief int_polyexp_erfc
    /// \param _ll
    /// \param _ul
    /// \param alpha
    /// \param beta
    /// \param gamma
    /// \param a
    /// \return
    ///
    double int_polyexp_erfc(const double& _ll, const double& _ul,
                            const double& alpha, const double& beta,
                            const double& gamma, const double& a) const;
    ///
    /// \brief norm_xEp_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_xEp_conv_gauss(const double& _ll, const double& _ul,
                               const double& tau, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_xEn_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_xEn_conv_gauss(const double& _ll, const double& _ul,
                               const double& tau, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_xEf_conv_gauss
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param m
    /// \param s
    /// \param o
    /// \return
    ///
    double norm_xEf_conv_gauss(const double& _ll, const double& _ul,
                               const double& tau, const double& m,
                               const double& s, const double& o = 0.) const;
    ///
    /// \brief norm_xEp
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param o
    /// \return
    ///
    double norm_xEp(const double& _ll, const double& _ul,
                    const double& tau, const double& o = 0.) const;
    ///
    /// \brief norm_xEn
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param o
    /// \return
    ///
    double norm_xEn(const double& _ll, const double& _ul,
                    const double& tau, const double& o = 0.) const;
    ///
    /// \brief norm_xEf
    /// \param _ll
    /// \param _ul
    /// \param tau
    /// \param o
    /// \return
    ///
    double norm_xEf(const double& _ll, const double& _ul,
                    const double& tau, const double& o = 0.) const;

 protected:
  // ** Constants ** //
    static const double inv_sqrt2;
    static const double inv_sqrt_pi;
    static const double sqrt_pi;
    static const double inv_sqrt_2pi;
    static const double sqrt2;
};

}  // namespace libTatami

#endif  // INCLUDE_CNVL_H_
