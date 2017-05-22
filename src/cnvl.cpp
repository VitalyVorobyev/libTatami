/** Copyright 2016 Vitaly Vorobyev
 ** @file cnvl.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include <complex>
#include <cmath>

#include "../include/cnvl.h"
#include "../include/typedefs.h"
#include "mylibs/Faddeeva/Faddeeva.h"

typedef std::complex<double> compld;

using std::fabs;
using std::exp;
using std::cos;
using std::erfc;
using std::log;

namespace libTatami {

const double cnvl::inv_sqrt2 =
        0.707106781186547461715008466853760182857513427734375;
const double cnvl::inv_sqrt_pi =
        0.56418958354775627928034964497783221304416656494140625;
const double cnvl::sqrt_pi =
        1.772453850905515881919427556567825376987457275390625;
const double cnvl::inv_sqrt_2pi =
        0.398942280401432702863218082711682654917240142822265625;
const double cnvl::sqrt2 =
        1.4142135623730951454746218587388284504413604736328125;

double cnvl::xXi_conv_gauss_by_int(nXi_t p_func,
                                 const double& t, const double& xd,
                                 const double& mu, const double& sigma) const {
    double area = 0.0;
    if (sigma == 0.0) return (*p_func)(t, xd);
    /* integral (-20sigma -- +20 sigma)*/
    const double t_ini = -20.0 * sigma + mu;
    const double t_end = +20.0 * sigma + mu;
    const int ndiv = 400;
    const double t_del = (t_end - t_ini) / ndiv;

    for (int i = 0; i < ndiv; ++i) {
        double t_run = t_ini + t_del * i;
        const double err = gaussian(t_run, t-mu, sigma);
        area += err * (*p_func)(t_run, xd);
    }
    area -= 0.5 * gaussian(t_ini, t-mu, sigma) * (*p_func)(t_ini, xd);
    area -= 0.5 * gaussian(t_end, t-mu, sigma) * (*p_func)(t_end, xd);
    area *= t_del;
    return area;
}

double cnvl::recexp(const double& re, const double& im) const {
    compld z(re, im);
    return exp(z).real();
}

double cnvl::imcexp(const double& re, const double& im) const {
    compld z(re, im);
    return exp(z).imag();
}

double cnvl::rewerf(const double& re, const double& im) const {
    compld z(re, im);
    return Faddeeva::w(z).real();
}

double cnvl::imwerf(const double& re, const double& im) const {
    compld z(re, im);
    return Faddeeva::w(z).imag();
}

double cnvl::Ep(const double& t, const double& tau) const {
    if (t < 0.) return 0.;
    const double inv_atau = 1. / fabs(tau);
    return inv_atau * exp(-fabs(t) * inv_atau);
}

double cnvl::En(const double& t, const double& tau) const {
    if (t >= 0.) return 0.;
    const double inv_atau = 1. / fabs(tau);
    return inv_atau*exp(-fabs(t) * inv_atau);
}

double cnvl::Ef(const double& t, const double& tau) const {
    const double inv_atau = 1. / fabs(tau);
    return 0.5 * inv_atau * exp(-fabs(t) * inv_atau);
}

double cnvl::Enp(const double& t, const double& tau_n,
                 const double& tau_p) const {
    const double norm = 1. / (fabs(tau_n) + fabs(tau_p));
    double inv_atau = (t >= 0) ? 1. / fabs(tau_p) : 1. / fabs(tau_n);
    return norm * exp(-fabs(t) * inv_atau);
}

double cnvl::xEp(const double& t, const double& tau) const {
    return t * Ep(t, tau) / fabs(tau);
}

double cnvl::xEn(const double& t, const double& tau) const {
    return -t * En(t, tau) / fabs(tau);
}

double cnvl::xEf(const double& t, const double& tau) const {
    return fabs(t) * Ef(t, tau) / fabs(tau);
}

double cnvl::nMp(const double& t, const double& xd) {
    if (t < 0.0) return 0.;
    return exp(-fabs(t)) * cos(xd * t);
}

double cnvl::Mp(const double& t, const double& tau, const double& dm) const {
    if (t < 0.0) return 0.0;
    return nMp(t / fabs(tau), dm * fabs(tau));
}

double cnvl::nMn(const double& t, const double& xd) {
    if (t >= 0.0) return 0.0;
    return exp(-fabs(t)) * cos(xd * t);
}

double cnvl::Mn(const double& t, const double& tau, const double& dm) const {
    if (t >= 0.) return 0.0;
    return nMn(t / fabs(tau), dm * fabs(tau));
}

double cnvl::nMf(const double& t, const double& xd) const {
    return exp(-fabs(t)) * cos(xd * t);
}

double cnvl::Mf(const double& t, const double& tau, const double& dm) const {
    return nMf(t / fabs(tau), dm * fabs(tau));
}

double cnvl::nAp(const double& t, const double& xd) {
    if (t < 0.0) return 0.0;
    return exp(-fabs(t)) * sin(xd * t);
}

double cnvl::Ap(const double& t, const double& tau, const double& dm) const {
    if (t < 0.0) return 0.0;
    return nAp(t / fabs(tau), dm * fabs(tau));
}

double cnvl::nAn(const double& t, const double& xd) {
    if (t >= 0.0) return 0.;
    return exp(-fabs(t)) * sin(xd * t);
}

double cnvl::An(const double& t, const double& tau, const double& dm) const {
    if (t >= 0.0) return 0.0;
    return nAn(t / fabs(tau), dm * fabs(tau));
}

double cnvl::nAf(const double& t, const double& xd) const {
    return exp(-fabs(t)) * sin(xd * t);
}

double cnvl::Af(const double& t, const double& tau, const double& dm) const {
    return nAf(t / fabs(tau), dm * fabs(tau));
}

double cnvl::norm_nEp(const double& _ll, const double& _ul,
                      const double& o) const {
    const double nul = (_ul - o >= 0.0) ? (_ul - o) : 0.0;
    const double nll = (_ll - o >= 0.0) ? (_ll - o) : 0.0;
    if (nul == nll) return 0;
    return exp(-nll) - exp(-nul);
}

double cnvl::norm_nEn(const double& _ll, const double& _ul,
                      const double& o) const {
    const double nul = (_ul - o < 0.0) ? (_ul - o) : 0.0;
    const double nll = (_ll - o < 0.0) ? (_ll - o) : 0.0;
    if (nul == nll) return 0;
    return exp(nul) - exp(nll);
}

double cnvl::norm_nEf(const double& _ll, const double& _ul,
                      const double& o) const {
    return 0.5 * (norm_nEn(_ll, _ul, o) + norm_nEp(_ll, _ul, o));
}

double cnvl::norm_Ep(const double& _ll, const double& _ul,
                     const double& tau, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    const double nll = _ll * inv_atau;
    const double nul = _ul * inv_atau;
    const double no  = o * inv_atau;
    return norm_nEp(nll, nul, no);
}

double cnvl::norm_En(const double& _ll, const double& _ul,
                     const double& tau, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    const double nll = _ll * inv_atau;
    const double nul = _ul * inv_atau;
    const double no  = o * inv_atau;
    return norm_nEn(nll, nul, no);
}

double cnvl::norm_Ef(const double& _ll, const double& _ul,
                     const double& tau, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    const double nll = _ll * inv_atau;
    const double nul = _ul * inv_atau;
    const double no  = o * inv_atau;
    return norm_nEf(nll, nul, no);
}

double cnvl::norm_Ap(const double& _ll, const double& _ul, const double& tau,
                     const double& dm, const double& o) const {
    double nll = (_ll < _ul) ? _ll : _ul;
    double nul = (_ll < _ul) ? _ul : _ll;
    nll -= o;
    nul -= o;
    if (nul <= 0.0) return 0.0;
    double f = norm_Ax_sup(nul, tau, dm);
    f -= norm_Ax_sup((nll > 0.0) ? nll : 0.0, tau, dm);
    return (_ll < _ul) ? f : -f;
}

double cnvl::norm_Ax_sup(const double& x, const double& tau,
                         const double& dm) const {
    const double dmt   = dm * x;
    const double dmtau = dm * tau;
    double f = -sin(dmt) - dmtau * cos(dmt);
    f *= tau * exp(-x / tau) / (1 + dmtau * dmtau);
    return f;
}

double cnvl::norm_Mp(const double& _ll, const double& _ul, const double& tau,
                     const double& dm, const double& o) const {
    double nll = (_ll < _ul) ? _ll : _ul;
    double nul = (_ll < _ul) ? _ul : _ll;
    nll -= o;
    nul -= o;
    if (nul <= 0.0) return 0.0;
    double f = norm_Mx_sup(nul, tau, dm);
    f -= norm_Mx_sup((nll > 0.0) ? nll : 0.0, tau, dm);
    return (_ll < _ul) ? f : -f;
}

double cnvl::norm_Mn(const double& _ll, const double& _ul, const double& tau,
                     const double& dm, const double&o) const {
    double nll = (_ll < _ul) ? _ll : _ul;
    double nul = (_ll < _ul) ? _ul : _ll;
    nll -= o;
    nul -= o;
    if (nll > 0.0) return 0.0;
    double f = norm_Mx_sup(-nll, tau, dm);
    f -= norm_Mx_sup((nul < 0.0) ? -nul : 0.0, tau, dm);
    return (_ll < _ul) ? f : -f;
}

double cnvl::norm_Mx_sup(const double& x, const double& tau,
                         const double& dm) const {
    const double dmt   = dm * x;
    const double dmtau = dm * tau;
    double f = -cos(dmt) - dmtau * sin(dmt);
    f *= tau  *exp(-x / tau) / (1 + dmtau * dmtau);
    return f;
}

double cnvl::norm_An(const double& _ll, const double& _ul, const double& tau,
                     const double& dm, const double& o) const {
    double nll = (_ll < _ul) ? _ll : _ul;
    double nul = (_ll < _ul) ? _ul : _ll;
    nll -= o;
    nul -= o;
    if (nll > 0.0) return 0.0;
    double f = norm_Ax_sup((nul < 0.0) ? -nul : 0.0, tau, dm);
    f -= norm_Ax_sup(-nll, tau, dm);
    return (_ll  <_ul) ? f : -f;
}

double cnvl::nEp_conv_gauss(const double& t, const double& m,
                            const double& s) const {
    static const double Tc = DBL_MAX_10_EXP * log(10.0);
    if (s == 0.0) return Ep(t - m, 1.0);
    double inv_s = 1.0 / fabs(s);
    double dt = -t + m;
    double ex = 0.5 * s * s + dt;
    if (ex < Tc) return 0.5 * exp(ex) * erfc(inv_sqrt2 * inv_s * (s * s + dt));
    const double gf = exp(-0.5 * dt * dt * inv_s * inv_s);
    const double x  = inv_sqrt2*(s + dt * inv_s);
    return 0.5 * approx_exp2erfc(x) * gf;
}

double cnvl::nEn_conv_gauss(const double& t, const double& m,
                            const double& s) const {
    if (s == 0.0) return En(t - m, 1.0);
    static const double Tc = DBL_MAX_10_EXP * log(10.0);
    double inv_s = 1.0 / fabs(s);
    double dt = t - m;
    double ex = 0.5 * s * s + dt;
    if (ex < Tc) return 0.5 * exp(ex) * erfc(inv_sqrt2 * inv_s * (s * s + dt));
    const double gf = exp(-0.5 * dt * dt * inv_s * inv_s);
    const double x  = inv_sqrt2 * (s + dt * inv_s);
    return 0.5 * approx_exp2erfc(x) * gf;
}

double cnvl::Ep_conv_gauss(const double& t, const double& tau,
                           const double& m, const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    return inv_atau * nEp_conv_gauss(nt, nm, ns);
}

double cnvl::En_conv_gauss(const double& t, const double& tau,
                           const double& m, const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    return inv_atau * nEn_conv_gauss(nt, nm, ns);
}

double cnvl::Ef_conv_gauss(const double& t, const double& tau,
                           const double& m, const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    double f = nEn_conv_gauss(nt, nm, ns) + nEp_conv_gauss(nt, nm, ns);
    return f * 0.5 * inv_atau;
}

double cnvl::Enp_conv_gauss(const double& t, const double& tau_n,
                            const double& tau_p, const double& m,
                            const double& s) const {
    const double inv_atau_n = 1. / fabs(tau_n);
    const double nt_n = t * inv_atau_n;
    const double nm_n = m * inv_atau_n;
    const double ns_n = fabs(s) * inv_atau_n;

    const double inv_atau_p = 1. / fabs(tau_p);
    const double nt_p = t * inv_atau_p;
    const double nm_p = m * inv_atau_p;
    const double ns_p = fabs(s) * inv_atau_p;

    double f = nEn_conv_gauss(nt_n, nm_n, ns_n) +
               nEp_conv_gauss(nt_p, nm_p, ns_p);
    return f * 1. / (fabs(tau_n) + fabs(tau_p));
}

double cnvl::approx_exp2erfc(const double& x) const {
    const double inv_x = 1. / x;
    return inv_sqrt_pi * (inv_x - 0.5 * inv_x * inv_x * inv_x);
}

double cnvl::nMp_conv_gauss(const double& t, const double& xd,
                            const double& m, const double& s) const {
    if (!finite(s)) return 0.0;
    if (s == 0.0)   return Mp(t - m, 1.0, xd);
    const double inv_s = 1. / fabs(s);
    if (!finite(inv_s)) return Mp(t - m, 1.0, xd);
    const double inv_s2 = 1. / (s * s);
    double f = _IM(t, m, s, 0.5 * inv_s2, 1 - (t - m) * inv_s2, xd);
    if (!finite(f)) return xXi_conv_gauss_by_int(nMp, t, xd, m, s);
    return f;
}

double cnvl::nMn_conv_gauss(const double& t, const double& xd,
                            const double& m, const double& s) const {
    if (!finite(s)) return 0.0;
    if (s == 0.0)   return Mn(t - m, 1., xd);
    const double inv_s = 1. / fabs(s);
    if (!finite(inv_s)) return Mn(t - m, 1., xd);
    const double inv_s2 = 1. / (s * s);
    double f = _IM(t, m, s, 0.5 * inv_s2, 1 + (t - m) * inv_s2, xd);
    if (!finite(f)) return xXi_conv_gauss_by_int(nMn, t, xd, m, s);
    return f;
}

double cnvl::nAp_conv_gauss(const double& t, const double& xd,
                            const double& m, const double& s) const {
    if (!finite(s)) return 0.;
    if (s == 0.)    return Ap(t - m, 1., xd);
    const double inv_s = 1. / fabs(s);
    if (!finite(inv_s)) return Ap(t - m, 1., xd);
    const double inv_s2 = 1. / s / s;
    double f = _IA(t, m, s, 0.5 * inv_s2, 1 - (t - m) * inv_s2, xd);
    if (!finite(f)) return xXi_conv_gauss_by_int(nAp, t, xd, m, s);
    return f;
}

double cnvl::nAn_conv_gauss(const double& t, const double& xd,
                            const double& m, const double& s) const {
    if (!finite(s)) return 0.0;
    if (s == 0.0)   return An(t - m, 1., xd);
    const double inv_s = 1. / fabs(s);
    if (!finite(inv_s)) return An(t - m, 1.0, xd);
    const double inv_s2 = 1. / s / s;
    double f = -1. * _IA(t, m, s, 0.5 * inv_s2, 1 + (t - m) * inv_s2, xd);
    if (!finite(f)) return xXi_conv_gauss_by_int(nAn, t, xd, m, s);
    return f;
}

double cnvl::Mp_conv_gauss(const double& t, const double& tau,
                           const double& dm, const double& m,
                           const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    return nMp_conv_gauss(t * inv_atau, dm * fabs(tau),
                          m * inv_atau, fabs(s) * inv_atau);
}

double cnvl::Mn_conv_gauss(const double& t, const double& tau,
                           const double& dm, const double& m,
                           const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    return nMn_conv_gauss(t * inv_atau, dm * fabs(tau),
                          m * inv_atau, fabs(s) * inv_atau);
}

double cnvl::Mf_conv_gauss(const double& t, const double& tau,
                           const double& dm, const double& m,
                           const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    const double xd = dm * fabs(tau);
    return nMn_conv_gauss(nt, xd, nm, ns) + nMp_conv_gauss(nt, xd, nm, ns);
}

double cnvl::Ap_conv_gauss(const double& t, const double& tau,
                           const double& dm, const double& m,
                           const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    return nAp_conv_gauss(t * inv_atau, dm * fabs(tau),
                          m * inv_atau, fabs(s) * inv_atau);
}

double cnvl::An_conv_gauss(const double& t, const double& tau,
                           const double& dm, const double& m,
                           const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    return nAn_conv_gauss(t * inv_atau, dm * fabs(tau),
                          m * inv_atau, fabs(s) * inv_atau);
}

double cnvl::Af_conv_gauss(const double& t, const double& tau,
                           const double& dm, const double& m,
                           const double& s) const {
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    const double xd = dm * fabs(tau);
    return nAn_conv_gauss(nt, xd, nm, ns) + nAp_conv_gauss(nt, xd, nm, ns);
}

double cnvl::_IM(const double& x, const double& m, const double& s,
                 const double& beta, const double& gamma,
                 const double& b) const {
#if defined(USE_CERNLIB_WERF)&&(IX_TREAT_NEGATIVE_GAMMA == 0)
    const double inv_2sqrtbeta = 0.5 / sqrt(beta);
    double f = rewerf(inv_2sqrtbeta * b, inv_2sqrtbeta * gamma);
    returnf * gaussian(x, m, s) * sqrt_pi * inv_2sqrtbeta;
#else /* USE_CERNLIB_WERF */
    const double inv_2sqrtbeta = 0.5 / sqrt(beta);
    double f = rewerf(inv_2sqrtbeta * b, inv_2sqrtbeta * fabs(gamma));
    f *= gaussian(x, m, s) * sqrt_pi * inv_2sqrtbeta;
    if (/* gamma */ inv_2sqrtbeta * gamma >=0) {
        return f;
    }
    f *= -1.0;
    const double inv_beta = 1. / beta;
    const double dx = x - m;
    const double inv_s = 1. / fabs(s);
    const double rex = 0.25 * (gamma * gamma - b * b) * inv_beta -
                  0.5 * dx * dx * inv_s * inv_s;
    double imx = 0.5 * fabs(gamma) * b * inv_beta;
    return f + inv_sqrt_2pi * inv_s * sqrt(M_PI * inv_beta) * recexp(rex, imx);
#endif
}

double cnvl::_IA(const double& x, const double& m, const double& s,
                 const double& beta, const double& gamma,
                 const double& b) const {
#if defined(USE_CERNLIB_WERF)&&(IX_TREAT_NEGATIVE_GAMMA == 0)
    const double inv_2sqrtbeta = 0.5 / sqrt(beta);
    double f = imwerf(inv_2sqrtbeta * b, inv_2sqrtbeta * gamma);
    return f * gaussian(x, m, s) * sqrt_pi * inv_2sqrtbeta;
#else /* USE_CERNLIB_WERF */
    const double inv_2sqrtbeta = 0.5 / sqrt(beta);
    double f = imwerf(inv_2sqrtbeta * b, inv_2sqrtbeta * fabs(gamma));
    f *= gaussian(x, m, s) * sqrt_pi * inv_2sqrtbeta;
    if (/*  gamma */inv_2sqrtbeta * gamma  >= 0) {
        return f;
    }
    const double inv_beta = 1. / beta;
    const double dx = x - m;
    const double inv_s = 1. / fabs(s);
    const double rex = 0.25 * (gamma * gamma - b * b) * inv_beta -
            0.5 * dx * dx * inv_s * inv_s;
    const double imx = 0.5 * fabs(gamma) * b * inv_beta;
    return f + inv_sqrt_2pi * inv_s * sqrt(M_PI * inv_beta) * imcexp(rex, imx);
#endif  /* USE_CERNLIB_WERF */
}

double cnvl::gaussian(const double& x, const double& m,
                      const double& s) const {
    if (!finite(s)) return 0;
    if (s == 0.)    return DiracDelta(x - m);
    double inv_s = 1. / fabs(s);
    if (!finite(inv_s)) return DiracDelta(x - m);
    const double dx = x - m;
    return inv_sqrt_2pi * inv_s * exp(-0.5 * dx * dx * inv_s * inv_s);
}

double cnvl::norm_gaussian_w_cutoff(const double& cutoff, const double& m,
                                    const double& s) const {
    double a_s = fabs(s);
    if (s == 0.) return 1;
    double inv_s = 1. / a_s;
    double x1 = (cutoff + m) * inv_sqrt2 * inv_s;
    double x2 = (cutoff - m) * inv_sqrt2 * inv_s;
    return 1. - 0.5 * erfc(x1) - 0.5 * erfc(x2);
}

double cnvl::norm_gaussian(const double& ll, const double& ul,
                           const double& m, const double& s) const {
    double a_s = fabs(s);
    if (s == 0.0) return 1;
    double inv_s = 1. / a_s;
    double x1 = (-ll + m) * inv_sqrt2 * inv_s;
    double x2 =  (ul - m) * inv_sqrt2 * inv_s;
    return 1. - 0.5 * (erfc(x1) + erfc(x2));
}

double cnvl::DiracDelta(const double& x) const {
    if (x == 0.) return FLT_MAX;
    return 0.;
}

double cnvl::xEn_conv_gauss(const double& t, const double& tau,
                            const double& m, const double& s) const {
    if (s == 0.) return xEn(t - m, tau);
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    const double ns2 = ns * ns;

    return (-ns2 - nt + nm) * En_conv_gauss(t, tau, m, s) +
            ns2 * gaussian(t, m, s);
}

double cnvl::xEf_conv_gauss(const double& t, const double& tau,
                            const double& m, const double& s) const {
    if (s == 0.) return xEn(t - m, tau);
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    const double ns2 = ns * ns;

    return 0.5 * (-ns2 + nt-nm) * Ep_conv_gauss(t, tau, m, s) +
           0.5 * (-ns2 - nt+nm) * En_conv_gauss(t, tau, m, s) +
           ns2*gaussian(t, m, s);
}

double cnvl::xEp_conv_gauss(const double& t, const double& tau,
                            const double& m, const double& s) const {
    if (s == 0.) return xEp(t - m, tau);
    const double inv_atau = 1. / fabs(tau);
    const double nt = t * inv_atau;
    const double nm = m * inv_atau;
    const double ns = fabs(s) * inv_atau;
    const double ns2 = ns * ns;

    return (-ns2 + nt - nm) * Ep_conv_gauss(t, tau, m, s) +
            ns2 * gaussian(t, m, s);
}

double cnvl::norm_neg_sup(const double& m, const double& s) const {
    const double Tc = DBL_MAX_10_EXP*log(10.0);

    const double as = fabs(s);
    const double ex = 0.5 * as * as - m;
    const double x = inv_sqrt2 * (as - m / as);
    if (ex < Tc) {
        double f = exp(ex);
        double g = erfc(x);
        if (g > 0.0 && finite(f)) {
            return f * g;
        }
    }
    double f = exp(-0.5 * m * m / (s * s));
    double g = approx_exp2erfc(x);
    return f * g;
}

double cnvl::norm_nEp_conv_gauss_sub(const double& _ll, const double& _ul,
                                     const double& m, const double& s,
                                     const double& o) const {
    const double f1 = erfc(inv_sqrt2 * (-_ll + o + m) / s) -
                 norm_neg_sup(_ll - o - m, s);
    const double f2 = erfc(inv_sqrt2 *  (_ul - o - m) / s) +
                 norm_neg_sup(_ul - o - m, s);
    return 0.5 * (f1 + f2);
}

double cnvl::norm_nEn_conv_gauss_sub(const double& _ll, const double& _ul,
                                     const double& m, const double& s,
                                     const double& o) const {
    const double f1 = erfc(inv_sqrt2 * (-_ll + o + m) / s) +
            norm_neg_sup(-_ll + o + m, s);
    const double f2 = erfc(inv_sqrt2 *  (_ul - o - m) / s)
            - norm_neg_sup(-_ul + o + m, s);
    return 0.5 * (f1 + f2);
}

double cnvl::norm_nEp_conv_gauss(const double& _ll, const double& _ul,
                                 const double& m, const double& s,
                                 const double& o) const {
    return 1. - norm_nEp_conv_gauss_sub(_ll, _ul, m, s, o);
}

double cnvl::norm_nEn_conv_gauss(const double& _ll, const double& _ul,
                                 const double& m, const double& s,
                                 const double& o) const {
    return 1. - norm_nEn_conv_gauss_sub(_ll, _ul, m, s, o);
}

double cnvl::norm_Ep_conv_gauss(const double& _ll, const double& _ul,
                                const double& tau, const double& m,
                                const double& s, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return 1. - norm_nEp_conv_gauss_sub(_ll * inv_atau, _ul * inv_atau,
                                        m * inv_atau, s * inv_atau,
                                        o * inv_atau);
}

double cnvl::norm_En_conv_gauss(const double& _ll, const double& _ul,
                                const double& tau, const double& m,
                                const double& s, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return 1. - norm_nEn_conv_gauss_sub(_ll * inv_atau, _ul * inv_atau,
                                        m * inv_atau, s * inv_atau,
                                        o * inv_atau);
}

double cnvl::norm_nEf_conv_gauss(const double& _ll, const double& _ul,
                                 const double& m, const double& s,
                                 const double& o) const {
  return 0.5 * (norm_nEn_conv_gauss(_ll, _ul, m, s, o) +
                norm_nEp_conv_gauss(_ll, _ul, m, s, o));
}

double cnvl::norm_Ef_conv_gauss(const double& _ll, const double& _ul,
                                const double& tau, const double& m,
                                const double& s, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    const double nll = _ll * inv_atau;
    const double nul = _ul * inv_atau;
    const double nm  = m * inv_atau;
    const double ns  = s * inv_atau;
    const double no  = (o == 0.) ? 0. : o * inv_atau;
    return 1. - 0.5 * (norm_nEp_conv_gauss_sub(nll, nul, nm, ns, no) +
                       norm_nEn_conv_gauss_sub(nll, nul, nm, ns, no));
}

double cnvl::norm_Enp_conv_gauss(const double& _ll, const double& _ul,
                                 const double& tau_n, const double& tau_p,
                                 const double& m, const double& s,
                                 const double& o) const {
    const double inv_atau_n = 1. / fabs(tau_n);
    const double nll_n = _ll * inv_atau_n;
    const double nul_n = _ul * inv_atau_n;
    const double nm_n  = m * inv_atau_n;
    const double ns_n  = s * inv_atau_n;
    const double no_n  = o * inv_atau_n;

    const double inv_atau_p = 1. / fabs(tau_p);
    const double nll_p = _ll * inv_atau_p;
    const double nul_p = _ul * inv_atau_p;
    const double nm_p  = m * inv_atau_p;
    const double ns_p  = s * inv_atau_p;
    const double no_p  = o * inv_atau_p;

    const double fn = fabs(tau_n) / (fabs(tau_n) + fabs(tau_p));
    const double fp = fabs(tau_p) / (fabs(tau_n) + fabs(tau_p));
    double f = fn * norm_nEn_conv_gauss_sub(nll_n, nul_n, nm_n, ns_n, no_n) +
               fp * norm_nEp_conv_gauss_sub(nll_p, nul_p, nm_p, ns_p, no_p);
    return 1. - f;
}

double cnvl::norm_nag_sup(const double& m, const double& s,
                          const double& xd) const {
    const double as = fabs(s);
    const double inv_as = 1. / as;
    const double minvassq = 0.5 * m * m * inv_as * inv_as;

    const double rex = inv_sqrt2 * xd * as;
    const double imx = inv_sqrt2 * (as + m * inv_as);

    if (imx >= 0) {
        return -(imwerf(rex, imx) + xd * rewerf(rex, imx)) * exp(-minvassq);
    }

    const double rexx = -rex * rex + imx * imx - minvassq;
    const double imxx = -2.0 * rex * imx;

    return exp(-minvassq) * (imwerf(-rex, -imx) + xd * rewerf(-rex, -imx)) -
            2. * (imcexp(rexx, imxx) + xd * recexp(rexx, imxx));
}

double cnvl::int_polyexp2(const double& _ll, const double& _ul,
                          const double& alpha, const double& beta,
                          const double& gamma, const double& a) const {
    const double inv_2a = 0.5 / a;
    const double sqrt_a = sqrt(a);

    return -inv_2a * (alpha * _ul + beta) * exp(-a * _ul * _ul) +
            inv_2a * (alpha * _ll + beta) * exp(-a * _ll * _ll) +
            (alpha * inv_2a + gamma) * sqrt_pi / sqrt_a *
            (1. - 0.5 * erfc(sqrt_a * _ul) - 0.5 * erfc(sqrt_a * _ll));
}

double cnvl::int_polyexp_erfc(const double& _ll, const double& _ul,
                              const double& alpha, const double& beta,
                              const double& gamma, const double& a) const {
    const double inv_a = 1. / a;
    double f = 0.0;

    const double a1 = -2. * alpha * inv_a + beta;
    const double a0 =  2. * alpha * alpha * inv_a * inv_a +
            beta * inv_a + a * gamma;

    f +=  inv_a * exp(a * _ul) *
            (alpha * _ul * _ul + a1 * _ul + a0) * erfc(_ul);
    f += -inv_a * exp(a * _ll) *
            (alpha * _ll * _ll + a1 * _ll + a0) * erfc(_ll);

    const double aa1 = a1 + a * alpha;
    const double aa0 = a0 + 0.5 * a * a1 + 0.25 * a * a * alpha;

    f += exp(0.25 * a * a) *
            int_polyexp2(_ll - 0.5 * a, _ul - 0.5 * a, alpha, aa1, aa0, 1.);
    return f;
}

double cnvl::norm_xEp_conv_gauss(const double& _ll, const double& _ul,
                                 const double& tau, const double& m,
                                 const double& s, const double& o) const {
    if (s == 0.) return norm_xEp(_ll - m, _ul - m, tau, o);
    const double inv_atau = 1. / fabs(tau);
    const double inv_s = 1. / fabs(s);

    const double ns = fabs(s) * inv_atau;
    const double ns2 = ns * ns;

    const double nll = (-(_ll - m - o) - ns * s) * inv_sqrt2 * inv_s;
    const double nul = (-(_ul - m - o) - ns * s) * inv_sqrt2 * inv_s;

    return ns2 * exp(1.5 * ns2) *
           int_polyexp_erfc(nll, nul, 0., 1., 0., ns * sqrt2) +
           ns2 * norm_gaussian(_ll, _ul, m + o, s);
}

double cnvl::norm_xEn_conv_gauss(const double& _ll, const double& _ul,
                                 const double& tau, const double& m,
                                 const double& s, const double& o) const {
    if (s == 0.) return norm_xEn(_ll - m, _ul - m, tau, o);
    const double inv_atau = 1. / fabs(tau);
    const double inv_s = 1. / fabs(s);

    const double ns = fabs(s) * inv_atau;
    const double ns2 = ns * ns;

    const double nll = ((_ll - m - o) - ns * s) * inv_sqrt2 * inv_s;
    const double nul = ((_ul - m - o) - ns * s) * inv_sqrt2 * inv_s;

    return ns2 * exp(1.5 * ns2) *
           int_polyexp_erfc(nll, nul, 0., 1., 0., ns * sqrt2) +
           ns2 * norm_gaussian(_ll, _ul, m + o, s);
}

double cnvl::norm_xEf_conv_gauss(const double& _ll, const double& _ul,
                                 const double& tau, const double& m,
                                 const double& s, const double& o) const {
    return 0.5 * (norm_xEp_conv_gauss(_ll, _ul, tau, m, s, o) +
                  norm_xEn_conv_gauss(_ll, _ul, tau, m, s, o));
}

double cnvl::norm_xEp(const double& _ll, const double& _ul,
                      const double& tau, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    const double nul = (_ul-o >= 0.) ? (_ul - o) * inv_atau : 0.;
    const double nll = (_ll-o >= 0.) ? (_ll - o) * inv_atau : 0.;
    if (nul == nll) return 0;
    const double vu = exp(-nul) * (-1. - nul);
    const double vl = exp(-nll) * (-1. - nll);
    return vu - vl;
}

double cnvl::norm_xEn(const double& _ll, const double& _ul,
                      const double& tau, const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    const double nul = (_ul - o < 0.) ? (_ul - o) * inv_atau : 0.;
    const double nll = (_ll - o < 0.) ? (_ll - o) * inv_atau : 0.;
    if (nul == nll) return 0;
    const double vu = exp(nul) * (1. - nul);
    const double vl = exp(nll) * (1. - nll);
    return vu - vl;
}

double cnvl::norm_xEf(const double& _ll, const double& _ul,
                      const double& tau, const double& o) const {
    return 0.5 * (norm_xEn(_ll, _ul, tau, o) + norm_xEp(_ll, _ul, tau, o));
}

double cnvl::norm_nmg_sup(const double& m, const double& s,
                          const double& xd) const {
    const double as = fabs(s);
    const double inv_as = 1. / as;

    double rex = inv_sqrt2 * xd * as;
    double imx = inv_sqrt2 * (as + m * inv_as);
    if (imx >= 0) {
        return -(rewerf(rex, imx) - xd * imwerf(rex, imx)) *
                 exp(-0.5 * m * m * inv_as * inv_as);
    }
    double rexx = -rex * rex + imx * imx - 0.5 * m * m * inv_as * inv_as;
    double imxx = -2. * rex * imx;

    return exp(-0.5 * m * m * inv_as * inv_as) *
           (rewerf(-rex, -imx) - xd * imwerf(-rex, -imx)) -
           2. * (recexp(rexx, imxx) - xd * imcexp(rexx, imxx));
}

double cnvl::norm_nAn_conv_gauss(const double& ll, const double& ul,
                                 const double& xd, const double& m,
                                 const double& s, const double& o) const {
    const double inv_as = 1. / fabs(s);
    return (-xd + 0.5 * (xd * erfc(inv_sqrt2 * inv_as * (-ll + o + m)) -
                         norm_nag_sup(ll - o - m, s, xd)) +
                  0.5 * (xd * erfc(inv_sqrt2 * inv_as *  (ul - o - m)) +
                         norm_nag_sup(ul - o - m, s, xd)) ) /
            (1. + xd * xd);
}

double cnvl::norm_nAp_conv_gauss(const double& ll, const double& ul,
                                 const double& xd, const double& m,
                                 const double& s, const double& o) const {
    const double inv_as = 1. / fabs(s);
    return (xd - 0.5 * (xd * erfc(inv_sqrt2 * inv_as * (-ll + o + m)) +
                        norm_nag_sup(-ll + o + m, s, xd)) -
                 0.5 * (xd * erfc(inv_sqrt2 * inv_as *  (ul - o - m)) -
                        norm_nag_sup(-ul + o + m, s, xd)) ) /
            (1. + xd * xd);
}

double cnvl::norm_nAf_conv_gauss(const double& ll, const double& ul,
                                 const double& xd, const double& m,
                                 const double& s, const double& o) const {
    return norm_nAn_conv_gauss(ll, ul, xd, m, s, o) +
           norm_nAp_conv_gauss(ll, ul, xd, m, s, o);
}

double cnvl::norm_nMn_conv_gauss(const double& ll, const double& ul,
                                 const double& xd, const double& m,
                                 const double& s, const double& o) const {
    const double inv_as = 1. / fabs(s);
    return (1. - 0.5 * (erfc(inv_sqrt2 * inv_as * (-ll + o + m)) -
                        norm_nmg_sup(ll - o - m, s, xd)) -
                 0.5 * (erfc(inv_sqrt2 * inv_as *  (ul - o - m)) +
                        norm_nmg_sup(ul - o - m, s, xd)) ) /
            (1. + xd * xd);
}

double cnvl::norm_nMp_conv_gauss(const double& ll, const double& ul,
                                 const double& xd, const double& m,
                                 const double& s, const double& o) const {
    const double inv_as = 1. / fabs(s);
    return (1. - 0.5 * (erfc(inv_sqrt2 * inv_as * (-ll + o + m)) +
                        norm_nmg_sup(-ll + o + m, s, xd)) -
                 0.5 * (erfc(inv_sqrt2 * inv_as *  (ul - o - m)) -
                        norm_nmg_sup(-ul + o + m, s, xd)) ) /
            (1. + xd * xd);
}

double cnvl::norm_nMf_conv_gauss(const double& ll, const double& ul,
                                 const double& xd, const double& m,
                                 const double& s, const double& o) const {
    return norm_nMn_conv_gauss(ll, ul, xd, m, s, o) +
           norm_nMp_conv_gauss(ll, ul, xd, m, s, o);
}

double cnvl::norm_An_conv_gauss(const double& ll, const double& ul,
                                const double& tau, const double& dm,
                                const double& m, const double& s,
                                const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return tau * norm_nAn_conv_gauss(ll * inv_atau, ul * inv_atau,
                                     dm * fabs(tau), m * inv_atau,
                                     s  * inv_atau,  o * inv_atau);
}

double cnvl::norm_Ap_conv_gauss(const double& ll, const double& ul,
                                const double& tau, const double& dm,
                                const double& m, const double& s,
                                const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return tau * norm_nAp_conv_gauss(ll * inv_atau, ul * inv_atau,
                                     dm * fabs(tau), m * inv_atau,
                                     s  * inv_atau,  o * inv_atau);
}

double cnvl::norm_Af_conv_gauss(const double& ll, const double& ul,
                                const double& tau, const double& dm,
                                const double& m, const double& s,
                                const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return tau * norm_nAf_conv_gauss(ll * inv_atau, ul * inv_atau,
                                     dm * fabs(tau), m * inv_atau,
                                     s  * inv_atau,  o * inv_atau);
}

double cnvl::norm_Mn_conv_gauss(const double& ll, const double& ul,
                                const double& tau, const double& dm,
                                const double& m, const double& s,
                                const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return tau * norm_nMn_conv_gauss(ll * inv_atau, ul * inv_atau,
                                     dm * fabs(tau), m * inv_atau,
                                     s  * inv_atau,  o * inv_atau);
}

double cnvl::norm_Mp_conv_gauss(const double& ll, const double& ul,
                                const double& tau, const double& dm,
                                const double& m, const double& s,
                                const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return tau * norm_nMp_conv_gauss(ll * inv_atau, ul * inv_atau,
                                     dm * fabs(tau), m * inv_atau,
                                     s  * inv_atau,  o * inv_atau);
}

double cnvl::norm_Mf_conv_gauss(const double& ll, const double& ul,
                                const double& tau, const double& dm,
                                const double& m, const double& s,
                                const double& o) const {
    const double inv_atau = 1. / fabs(tau);
    return tau * norm_nMf_conv_gauss(ll * inv_atau, ul * inv_atau,
                                     dm * fabs(tau), m * inv_atau,
                                     s  * inv_atau,  o * inv_atau);
}

}  // namespace libTatami
