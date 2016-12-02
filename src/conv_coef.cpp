/** Copyright 2016 Vitaly Vorobyev
 ** @file conv_coef.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/conv_coef.h"

using std::abs;

namespace libTatami {

const int conv_coef::E = 1;
const int conv_coef::A = 2;
const int conv_coef::M = 3;

const fEnp& conv_coef::Ecoefs(const RascRnpPars& pRasc, const RkPar& pRk) {
    if (state == E) return f;
    state = E; f.Clear();
    cdouble fe = 1. - pRasc.fd();
    cdouble nfn = (1. - pRasc.fp()) * fe;
    cdouble nfp = pRasc.fp() * fe;
    f.fEn = pRk.fact_n_e() * pRasc.fd();
    f.fEp = pRk.fact_p_e() * pRasc.fd();
    add_EnEn_coef(&f.fEn, &f.fEn_np, &f.fxEn,
                  pRk.ntau_n(), pRasc.tau_np_n(), pRk.fact_n_e() * nfn);
    add_EnEp_coef(&f.fEn, &f.fEp_np,
                  pRk.ntau_n(), pRasc.tau_np_p(), pRk.fact_n_e() * nfp);
    add_EpEn_coef(&f.fEp, &f.fEn_np,
                  pRk.ntau_p(), pRasc.tau_np_n(), pRk.fact_p_e() * nfn);
    add_EpEp_coef(&f.fEp, &f.fEp_np, &f.fxEp,
                  pRk.ntau_p(), pRasc.tau_np_p(), pRk.fact_p_e() * nfp);
    return f;
}

const fEnp& conv_coef::Acoefs(const RascRnpPars& pRasc, const RkPar& pRk) {
    if (state == A) return f;
    state = A; f.Clear();
    cdouble fe = 1.0 - pRasc.fd();
    cdouble nfn = (1.0 - pRasc.fp()) * fe;
    cdouble nfp = pRasc.fp() * fe;
    cdouble w_mn_n = -pRk.fact_am() * pRk.ndmtau();
    cdouble w_mn_p = -pRk.fact_am() * pRk.ndmtau();
    f.fAn = pRk.fact_am() * pRasc.fd();
    f.fAp = pRk.fact_am() * pRasc.fd();
    f.fMn = w_mn_n * pRasc.fd();
    f.fMp = w_mn_p * pRasc.fd();
    add_AnEn_coef(&f.fAn, &f.fMn, &f.fEn_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_n(), pRk.fact_am() * nfn);
    add_AnEp_coef(&f.fAn, &f.fMn, &f.fEp_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_p(), pRk.fact_am() * nfp);
    add_MnEn_coef(&f.fMn, &f.fAn, &f.fEn_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_n(), w_mn_n * nfn);
    add_MnEp_coef(&f.fMn, &f.fAn, &f.fEp_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_p(), w_mn_n * nfp);
    add_ApEn_coef(&f.fAp, &f.fMp, &f.fEn_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_n(), pRk.fact_am() * nfn);
    add_ApEp_coef(&f.fAp, &f.fMp, &f.fEp_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_p(), pRk.fact_am() * nfp);
    add_MpEn_coef(&f.fMp, &f.fAp, &f.fEn_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_n(), w_mn_p * nfn);
    add_MpEp_coef(&f.fMp, &f.fAp, &f.fEp_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_p(), w_mn_p * nfp);
    return f;
}

const fEnp& conv_coef::Mcoefs(const RascRnpPars& pRasc, const RkPar& pRk) {
    if (state == M) return f;
    state = M; f.Clear();

    cdouble fe = 1. - pRasc.fd();
    cdouble nfn = (1. - pRasc.fp()) * fe;
    cdouble nfp = pRasc.fp() * fe;
    cdouble w_an_n = pRk.fact_am() * pRk.ndmtau();
    cdouble w_an_p = pRk.fact_am() * pRk.ndmtau();
    f.fMn = pRk.fact_am() * pRasc.fd();
    f.fMp = pRk.fact_am() * pRasc.fd();
    f.fAn = w_an_n * pRasc.fd();
    f.fAp = w_an_p * pRasc.fd();
    add_AnEn_coef(&f.fAn, &f.fMn, &f.fEn_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_n(), w_an_n * nfn);
    add_AnEp_coef(&f.fAn, &f.fMn, &f.fEp_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_p(), w_an_n * nfp);
    add_MnEn_coef(&f.fMn, &f.fAn, &f.fEn_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_n(), pRk.fact_am() * nfn);
    add_MnEp_coef(&f.fMn, &f.fAn, &f.fEp_np, pRk.ntau_n(),
                  pRk.ndm_n(), pRasc.tau_np_p(), pRk.fact_am() * nfp);
    add_ApEn_coef(&f.fAp, &f.fMp, &f.fEn_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_n(), w_an_p * nfn);
    add_ApEp_coef(&f.fAp, &f.fMp, &f.fEp_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_p(), w_an_p * nfp);
    add_MpEn_coef(&f.fMp, &f.fAp, &f.fEn_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_n(), pRk.fact_am() * nfn);
    add_MpEp_coef(&f.fMp, &f.fAp, &f.fEp_np, pRk.ntau_p(),
                  pRk.ndm_p(), pRasc.tau_np_p(), pRk.fact_am() * nfp);
    return f;
}

void conv_coef::add_EpEn_coef(double* fEp1, double* fEn2,
                    cdouble& tau1, cdouble& tau2, cdouble& weight) const {
    cdouble inv_tausum = 1. / (tau1 + tau2);
    *fEp1 += weight * tau1 * inv_tausum;
    *fEn2 += weight * tau2 * inv_tausum;
}

void conv_coef::add_EnEp_coef(double* fEn1, double* fEp2,
                    cdouble& tau1, cdouble& tau2, cdouble& weight) const {
    cdouble inv_tausum = 1. / (tau1 + tau2);
    *fEn1 += weight * tau1 * inv_tausum;
    *fEp2 += weight * tau2 * inv_tausum;
}

void conv_coef::add_EnEn_coef(double* fEn1, double* fEn2, double* fxEn1,
                    cdouble& tau1, cdouble& tau2, cdouble& weight) const {
    if (tau1 == tau2  /* ||!finite(1.0/(tau1-tau2)) */) {
        *fxEn1 += weight;
    } else {
        cdouble inv_tausub = 1. / (tau1 - tau2);
        *fEn1 +=  weight * tau1 * inv_tausub;
        *fEn2 += -weight * tau2 * inv_tausub;
    }
}

void conv_coef::add_EpEp_coef(double* fEp1, double* fEp2, double* fxEp1,
                    cdouble& tau1, cdouble& tau2, cdouble& weight) const {
    if (tau1 == tau2  /* ||!finite(1.0/(tau1-tau2)) */) {
        *fxEp1 += weight;
    } else {
        cdouble inv_tausub = 1. / (tau1 - tau2);
        *fEp1 +=  weight * tau1 * inv_tausub;
        *fEp2 += -weight * tau2 * inv_tausub;
    }
}

void conv_coef::add_ApEp_coef(double* fAp1, double* fMp1, double* fEp2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1 - inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fAp1 += -weight * A * inv_at2 * inv_tau;
    *fMp1 += -weight * A * inv_at2 * dm;
    *fEp2 +=  weight * A * dm;
}

void conv_coef::add_AnEn_coef(double* fAn1, double* fMn1, double* fEn2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1 - inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fAn1 += -weight * A * inv_at2 * inv_tau;
    *fMn1 += +weight * A * inv_at2 * dm;
    *fEn2 += -weight * A * dm;
}

void conv_coef::add_ApEn_coef(double* fAp1, double* fMp1, double* fEn2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1 + inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fAp1 += weight * A * inv_at2 * inv_tau;
    *fMp1 += weight * A * inv_at2 * dm;
    *fEn2 += weight * A * dm;
}

void conv_coef::add_AnEp_coef(double* fAn1, double* fMn1, double* fEp2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1 + inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fAn1 +=  weight * A * inv_at2 * inv_tau;
    *fMn1 += -weight * A * inv_at2 * dm;
    *fEp2 += -weight * A * dm;
}

void conv_coef::add_MpEp_coef(double* fMp1, double* fAp1, double* fEp2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1 - inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fMp1 += -weight * A * inv_at2 * inv_tau;
    *fAp1 +=  weight * A * inv_at2 * dm;
    *fEp2 +=  weight * A * inv_tau;
}

void conv_coef::add_MnEn_coef(double* fMn1, double* fAn1, double* fEn2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1  -inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fMn1 += -weight * A * inv_at2 * inv_tau;
    *fAn1 += -weight * A * inv_at2 * dm;
    *fEn2 +=  weight * A * inv_tau;
}

void conv_coef::add_MpEn_coef(double* fMp1, double* fAp1, double* fEn2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1 + inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fMp1 +=  weight * A * inv_at2 * inv_tau;
    *fAp1 += -weight * A * inv_at2 * dm;
    *fEn2 +=  weight * A * inv_tau;
}

void conv_coef::add_MnEp_coef(double* fMn1, double* fAn1, double* fEp2,
          cdouble& tau1, cdouble& dm, cdouble& tau2, cdouble& weight) const {
    cdouble inv_at1 = 1. / abs(tau1);
    cdouble inv_at2 = 1. / abs(tau2);
    cdouble inv_tau = inv_at1 + inv_at2;
    cdouble A = 1. / (inv_tau * inv_tau + dm * dm);
    *fMn1 +=  weight * A * inv_at2 * inv_tau;
    *fAn1 +=  weight * A * inv_at2 * dm;
    *fEp2 +=  weight * A * inv_tau;
}

}  // namespace libTatami
