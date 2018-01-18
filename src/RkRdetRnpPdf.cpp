/** Copyright 2016 Vitaly Vorobyev
 ** @file RkRdetRnpPdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "RkRdetRnpPdf.h"

#include "icpvevent.h"
#include "parmanager.h"
#include "ttools.h"

using std::fabs;
using std::cout;
using std::endl;

namespace libTatami {

#ifndef DTRES_EXTERAM_THRE
#define DTRES_EXTERAM_THRE -FLT_MAX
#endif

RkRdetRnpPdf::RkRdetRnpPdf(const DataClass& dc) :
    AbsICPVPdf(),
    m_wtag(ParManager::WTagFile(dc)),
    m_cnst(ParManager::SigParFile(dc)) {}

int RkRdetRnpPdf::ReadVars(const ICPVEvt &evt) const {
    dt  = evt.FindDVar("dt");
    m_rvar.ReadVars(evt, RdetVar::RecSide);
    m_rvar.ReadVars(evt, RdetVar::RecSide);
    keeptagl = evt.FindIVar("keeptagl");
    costhBcms = evt.FindDVar("costhBcms");
    return 0;
}

int RkRdetRnpPdf::ReadAndCalc(const ICPVEvt& evt) const {
    ReadVars(evt);
    m_rpar.Calculate(m_cnst, m_rvar);
    m_apar.Calculate(m_cnst, m_avar, keeptagl);
    m_apar.swap_rnp_param();
    m_kpar.SetAkCk(costhBcms, 0.5*10.58, tau(), dm());

    // precalculate means and sigmas //
    mu_mm = m_rpar.mu_main_rec() + m_apar.mu_main_asc();
    s_mm = TTools::sum_sigma(m_rpar.Smain_rec(), m_apar.Smain_asc());
    if (m_apar.ftail_asc() > 0.) {
        mu_mt = m_rpar.mu_main_rec() + m_apar.mu_tail_asc();
        s_mt = TTools::sum_sigma(m_rpar.Smain_rec(), m_apar.Stail_asc());
        if (m_rpar.ftail_rec() > 0.) {  // single track track (Asc) &&
                                        // single track vertex (rec)
            mu_tm = m_rpar.mu_tail_rec() + m_apar.mu_main_asc();
            s_tm = TTools::sum_sigma(m_rpar.Stail_rec(), m_apar.Smain_asc());
            mu_tt = m_rpar.mu_tail_rec() + m_apar.mu_tail_asc();
            s_tt = TTools::sum_sigma(m_rpar.Stail_rec(), m_apar.Stail_asc());
        }
    } else if (m_rpar.ftail_rec() > 0.) {  // multiple track track (Asc) &&
                                           // single track vertex (rec)
        mu_tm = m_rpar.mu_tail_rec() + m_apar.mu_main_asc();
        s_tm = TTools::sum_sigma(m_rpar.Stail_rec(), m_apar.Smain_asc());
    }
    return 0;
}

double RkRdetRnpPdf::Xsum4(double Xmm, double Xmt, double Xtm,
                           double Xtt) const {
    return (1. - m_rpar.ftail_rec() ) *
           (1. - m_apar.ftail_asc() ) * Xmm +
                 m_rpar.ftail_rec() *
           (1. - m_apar.ftail_asc() ) * Xtm +
           (1. - m_rpar.ftail_rec() ) *
                 m_apar.ftail_asc() * Xmt +
                 m_rpar.ftail_rec() *
                 m_apar.ftail_asc() * Xtt;
}

double RkRdetRnpPdf::Xsum2mmmt(double Xmm, double Xmt) const {
    return (1. - m_apar.ftail_asc()) * Xmm + m_apar.ftail_asc() * Xmt;
}

double RkRdetRnpPdf::Xsum2mmtm(double Xmm, double Xtm) const {
    return (1. - m_rpar.ftail_rec()) * Xmm + m_rpar.ftail_rec() * Xtm;
}

double RkRdetRnpPdf::EfRkRdetRnp_fullrec(const ICPVEvt& evt) const {
    ReadAndCalc(evt);
    return EfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::EfRkRdetRnp_fullrec(void) const {
    const double Li_mm = EfRkRdetRnp_full_sup(mu_mm, s_mm);
    if (m_apar.ftail_asc() > 0.) {
        const double Li_mt = EfRkRdetRnp_full_sup(mu_mt, s_mt);
        if (m_rpar.ftail_rec() > 0.) {  // single track track (Asc) &&
                                        // single track vertex (rec)
            const double Li_tm = EfRkRdetRnp_full_sup(mu_tm, s_tm);
            const double Li_tt = EfRkRdetRnp_full_sup(mu_tt, s_tt);
            return Xsum4(Li_mm, Li_mt, Li_tm, Li_tt);
        } else {  /* single track track (Asc) && multiple track vertex (rec)*/
            return Xsum2mmmt(Li_mm, Li_mt);
        }
    } else if (m_rpar.ftail_rec() > 0.) {  // multiple track track (Asc) &&
                                           //   single track vertex (rec)
        const double Li_tm = EfRkRdetRnp_full_sup(mu_tm, s_tm);
        return Xsum2mmtm(Li_mm, Li_tm);
    }
    /* multiple track track (Asc) && multiple track vertex (rec)*/
    return Li_mm;
}

double RkRdetRnpPdf::EfRkRdetRnp_full_sup(double mu, double sigma) const {
    const fEnp& f = coco.Ecoefs(m_apar, m_kpar);
    double Li = 0.;
    if (f.fEn != 0.)
        Li += f.fEn    *  En_conv_gauss(dt, m_kpar.ntau_n(),   mu, sigma);
    if (f.fEp != 0.)
        Li += f.fEp    *  Ep_conv_gauss(dt, m_kpar.ntau_p(),   mu, sigma);
    if (f.fEn_np != 0.)
        Li += f.fEn_np *  En_conv_gauss(dt, m_apar.tau_np_n(), mu, sigma);
    if (f.fEp_np != 0.)
        Li += f.fEp_np *  Ep_conv_gauss(dt, m_apar.tau_np_p(), mu, sigma);
    if (f.fxEn != 0.)
        Li += f.fxEn   * xEn_conv_gauss(dt, m_kpar.ntau_n(),   mu, sigma);
    if (f.fxEp != 0.)
        Li += f.fxEp   * xEp_conv_gauss(dt, m_kpar.ntau_p(),   mu, sigma);
    return Li;
}

double RkRdetRnpPdf::AfRkRdetRnp_fullrec(const ICPVEvt& evt) const {
    ReadAndCalc(evt);
    return AfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::AfRkRdetRnp_fullrec(void) const {
    const double Li_mm = AfRkRdetRnp_full_sup(mu_mm, s_mm);
    if (m_apar.ftail_asc() > 0.) {
        const double Li_mt = AfRkRdetRnp_full_sup(mu_mt, s_mt);
        if (m_rpar.ftail_rec() > 0.) {  // single track track (Asc) &&
                                        // single track vertex (rec)
            const double Li_tm = AfRkRdetRnp_full_sup(mu_tm, s_tm);
            const double Li_tt = AfRkRdetRnp_full_sup(mu_tt, s_tt);
            return Xsum4(Li_mm, Li_mt, Li_tm, Li_tt);
        } else {  /* single track track (Asc) && multiple track vertex (rec)*/
            return Xsum2mmmt(Li_mm, Li_mt);
        }
    } else if (m_rpar.ftail_rec() > 0.) {  // multiple track track (Asc) &&
                                           // single track vertex (rec)
        const double Li_tm = AfRkRdetRnp_full_sup(mu_tm, s_tm);
        return Xsum2mmtm(Li_mm, Li_tm);
      }
    /* multiple track track (Asc) && multiple track vertex (rec)*/
    return Li_mm;
}

double RkRdetRnpPdf::AfRkRdetRnp_full_sup(double mu, double sigma) const {
    const fEnp& f = coco.Acoefs(m_apar, m_kpar);
    double Li = 0.;
    if (f.fAn != 0.)    Li += f.fAn * An_conv_gauss(dt, m_kpar.ntau_n(),
                                                    m_kpar.ndm_n(), mu, sigma);
    if (f.fAp != 0.)    Li += f.fAp * Ap_conv_gauss(dt, m_kpar.ntau_p(),
                                                    m_kpar.ndm_p(), mu, sigma);
    if (f.fMn != 0.)    Li += f.fMn * Mn_conv_gauss(dt, m_kpar.ntau_n(),
                                                    m_kpar.ndm_n(), mu, sigma);
    if (f.fMp != 0.)    Li += f.fMp * Mp_conv_gauss(dt, m_kpar.ntau_p(),
                                                    m_kpar.ndm_p(), mu, sigma);
    if (f.fEn_np != 0.) Li += f.fEn_np * En_conv_gauss(dt,
                                                 m_apar.tau_np_n(), mu, sigma);
    if (f.fEp_np != 0.) Li += f.fEp_np * Ep_conv_gauss(dt,
                                                 m_apar.tau_np_p(), mu, sigma);
    return Li;
}

double RkRdetRnpPdf::MfRkRdetRnp_fullrec(const ICPVEvt& evt) const {
    ReadAndCalc(evt);
    return MfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::MfRkRdetRnp_fullrec(void) const {
    const double Li_mm = MfRkRdetRnp_full_sup(mu_mm, s_mm);
    if (m_apar.ftail_asc() > 0.) {
        const double Li_mt = MfRkRdetRnp_full_sup(mu_mt, s_mt);
        if (m_rpar.ftail_rec() > 0.) {  // single track track (Asc) &&
                                        // single track vertex (rec)
            const double Li_tm = MfRkRdetRnp_full_sup(mu_tm, s_tm);
            const double Li_tt = MfRkRdetRnp_full_sup(mu_tt, s_tt);
            return Xsum4(Li_mm, Li_mt, Li_tm, Li_tt);
        } else {  /* single track track (Asc) && multiple track vertex (rec)*/
            return Xsum2mmmt(Li_mm, Li_mt);
        }
    } else if (m_rpar.ftail_rec() > 0.) {  // multiple track track (Asc) &&
                                           //   single track vertex (rec)
        const double Li_tm = MfRkRdetRnp_full_sup(mu_tm, s_tm);
        return Xsum2mmtm(Li_mm, Li_tm);
    }
    /* multiple track track (Asc) && multiple track vertex (rec)*/
    return Li_mm;
}

double RkRdetRnpPdf::MfRkRdetRnp_full_sup(double mu, double sigma) const {
    const fEnp& f = coco.Mcoefs(m_apar, m_kpar);
    double Li = 0.;
    if (f.fAn != 0.)    Li += f.fAn *     An_conv_gauss(dt,
                                   m_kpar.ntau_n(), m_kpar.ndm_n(), mu, sigma);
    if (f.fAp != 0.)    Li += f.fAp *     Ap_conv_gauss(dt,
                                   m_kpar.ntau_p(), m_kpar.ndm_p(), mu, sigma);
    if (f.fMn != 0.)    Li += f.fMn *     Mn_conv_gauss(dt,
                                   m_kpar.ntau_n(), m_kpar.ndm_n(), mu, sigma);
    if (f.fMp != 0.)    Li += f.fMp *     Mp_conv_gauss(dt,
                                   m_kpar.ntau_p(), m_kpar.ndm_p(), mu, sigma);
    if (f.fEn_k != 0.)  Li += f.fEn_k *   En_conv_gauss(dt,
                                                   -m_kpar.cktau(), mu, sigma);
    if (f.fEp_k != 0.)  Li += f.fEp_k *   Ep_conv_gauss(dt,
                                                    m_kpar.cktau(), mu, sigma);
    if (f.fxEn_k != 0.) Li += f.fxEn_k * xEn_conv_gauss(dt,
                                                   -m_kpar.cktau(), mu, sigma);
    if (f.fxEp_k != 0.) Li += f.fxEp_k * xEp_conv_gauss(dt,
                                                    m_kpar.cktau(), mu, sigma);
    if (f.fEn_np != 0.) Li += f.fEn_np *  En_conv_gauss(dt,
                                                 m_apar.tau_np_n(), mu, sigma);
    if (f.fEp_np != 0.) Li += f.fEp_np *  Ep_conv_gauss(dt,
                                                 m_apar.tau_np_p(), mu, sigma);
    return Li;
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_fullrec(const ICPVEvt& evt) const {
    ReadAndCalc(evt);
    return norm_EfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_fullrec(void) const {
    const double Li_mm = norm_EfRkRdetRnp_full_sup(mu_mm, s_mm);
    if (m_apar.ftail_asc() > 0.) {
        const double Li_mt = norm_EfRkRdetRnp_full_sup(mu_mt, s_mt);
        if (m_rpar.ftail_rec() > 0.) {  // single track track (Asc) &&
                                        // single track vertex (rec)
            const double Li_tm = norm_EfRkRdetRnp_full_sup(mu_tm, s_tm);
            const double Li_tt = norm_EfRkRdetRnp_full_sup(mu_tt, s_tt);
            return Xsum4(Li_mm, Li_mt, Li_tm, Li_tt);
        } else { /* single track track (Asc) && multiple track vertex (rec) */
            return Xsum2mmmt(Li_mm, Li_mt);
        }
    } else if (m_rpar.ftail_rec() > 0.) {  // multiple track track (Asc) &&
                                           // single track vertex (rec)
        const double Li_tm = norm_EfRkRdetRnp_full_sup(mu_tm, s_tm);
        return Xsum2mmtm(Li_mm, Li_tm);
  }
  /* multiple track track (Asc) && multiple track vertex (rec) */
  return Li_mm;
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_full_sup(double mu, double sigma) const {
    const fEnp& f = coco.Ecoefs(m_apar, m_kpar);
    double Li = 0.;
    if (f.fEn != 0.)    Li += f.fEn *    norm_En_conv_gauss(ll(), ul(),
                                           m_kpar.ntau_n(), mu, sigma);
    if (f.fEp != 0.)    Li += f.fEp *    norm_Ep_conv_gauss(ll(), ul(),
                                           m_kpar.ntau_p(), mu, sigma);
    if (f.fEn_np != 0.) Li += f.fEn_np * norm_En_conv_gauss(ll(), ul(),
                                         m_apar.tau_np_n(), mu, sigma);
    if (f.fEp_np != 0.) Li += f.fEp_np * norm_Ep_conv_gauss(ll(), ul(),
                                         m_apar.tau_np_p(), mu, sigma);
    if (f.fxEn != 0.)   Li += f.fxEn *  norm_xEn_conv_gauss(ll(), ul(),
                                           m_kpar.ntau_n(), mu, sigma);
    if (f.fxEp != 0.)   Li += f.fxEp *  norm_xEp_conv_gauss(ll(), ul(),
                                           m_kpar.ntau_p(), mu, sigma);
    return Li;
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_fullrec(const ICPVEvt& evt) const {
    ReadAndCalc(evt);
    return norm_MfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_fullrec(void) const {
    const double Li_mm = norm_MfRkRdetRnp_full_sup(mu_mm, s_mm);
    if (m_apar.ftail_asc() > 0.) {
        const double Li_mt = norm_MfRkRdetRnp_full_sup(mu_mt, s_mt);
        if (m_rpar.ftail_rec() > 0.) {  // single track track (Asc) &&
                                        // single track vertex (rec)
            const double Li_tm = norm_MfRkRdetRnp_full_sup(mu_tm, s_tm);
            const double Li_tt = norm_MfRkRdetRnp_full_sup(mu_tt, s_tt);
            return Xsum4(Li_mm, Li_mt, Li_tm, Li_tt);
        } else { /* single track track (Asc) && multiple track vertex (rec) */
            return Xsum2mmmt(Li_mm, Li_mt);
        }
    } else if (m_rpar.ftail_rec() > 0.) {  // multiple track track (Asc) &&
                                           // single track vertex (rec)
        const double Li_tm = norm_MfRkRdetRnp_full_sup(mu_tm, s_tm);
        return Xsum2mmtm(Li_mm, Li_tm);
    }
    /* multiple track track (Asc) && multiple track vertex (rec) */
    return Li_mm;
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_full_sup(double mu, double sigma) const {
    const fEnp& f = coco.Mcoefs(m_apar, m_kpar);
    double Li = 0.;
    if (f.fAn != 0.)    Li += f.fAn *    norm_An_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_n(), m_kpar.ndm_n(), mu, sigma);
    if (f.fAp != 0.)    Li += f.fAp *    norm_Ap_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_p(), m_kpar.ndm_p(), mu, sigma);
    if (f.fMn != 0.)    Li += f.fMn *    norm_Mn_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_n(), m_kpar.ndm_n(), mu, sigma);
    if (f.fMp != 0.)    Li += f.fMp *    norm_Mp_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_p(), m_kpar.ndm_p(), mu, sigma);
    if (f.fEn_k != 0.)  Li += f.fEn_k *  norm_En_conv_gauss(ll(), ul(),
                                         -m_kpar.cktau(), mu, sigma);
    if (f.fEp_k != 0.)  Li += f.fEp_k *  norm_Ep_conv_gauss(ll(), ul(),
                                          m_kpar.cktau(), mu, sigma);
    if (f.fxEn_k != 0.) Li += f.fxEn_k * norm_xEn_conv_gauss(ll(), ul(),
                                         -m_kpar.cktau(), mu, sigma);
    if (f.fxEp_k != 0.) Li += f.fxEp_k * norm_xEp_conv_gauss(ll(), ul(),
                                          m_kpar.cktau(), mu, sigma);
    if (f.fEn_np != 0.) Li += f.fEn_np * norm_En_conv_gauss(ll(), ul(),
                                          m_apar.tau_np_n(), mu, sigma);
    if (f.fEp_np != 0.) Li += f.fEp_np * norm_Ep_conv_gauss(ll(), ul(),
                                          m_apar.tau_np_p(), mu, sigma);
    return Li;
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_fullrec(const ICPVEvt& evt) const {
    ReadAndCalc(evt);
    return norm_AfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_fullrec(void) const {
    const double Li_mm = norm_AfRkRdetRnp_full_sup(mu_mm, s_mm);
    if (m_apar.ftail_asc() > 0.) {
        const double Li_mt = norm_AfRkRdetRnp_full_sup(mu_mt, s_mt);
        if (m_rpar.ftail_rec() > 0.) {  // single track track (Asc) &&
                                        // single track vertex (rec)
            const double Li_tm = norm_AfRkRdetRnp_full_sup(mu_tm, s_tm);
            const double Li_tt = norm_AfRkRdetRnp_full_sup(mu_tt, s_tt);
            return Xsum4(Li_mm, Li_mt, Li_tm, Li_tt);
        } else { /* single track track (Asc) && multiple track vertex (rec) */
            return Xsum2mmmt(Li_mm, Li_mt);
        }
    } else if (m_rpar.ftail_rec() > 0.) {  // multiple track track (Asc) &&
                                           // single track vertex (rec)
        const double Li_tm = norm_AfRkRdetRnp_full_sup(mu_tm, s_tm);
        return Xsum2mmtm(Li_mm, Li_tm);
    }
    /* multiple track track (Asc) && multiple track vertex (rec) */
    return Li_mm;
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_full_sup(double mu, double sigma) const {
    const fEnp& f = coco.Acoefs(m_apar, m_kpar);
    double Li = 0.;
    if (f.fAn != 0.) Li += f.fAn * norm_An_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_n(), m_kpar.ndm_n(), mu, sigma);
    if (f.fAp != 0.) Li += f.fAp * norm_Ap_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_p(), m_kpar.ndm_p(), mu, sigma);
    if (f.fMn != 0.) Li += f.fMn * norm_Mn_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_n(), m_kpar.ndm_n(), mu, sigma);
    if (f.fMp != 0.) Li += f.fMp * norm_Mp_conv_gauss(ll(), ul(),
                                   m_kpar.ntau_p(), m_kpar.ndm_p(), mu, sigma);
    if (f.fEn_np != 0.) Li += f.fEn_np * norm_En_conv_gauss(ll(), ul(),
                                         m_apar.tau_np_n(), mu, sigma);
    if (f.fEp_np != 0.) Li += f.fEp_np * norm_Ep_conv_gauss(ll(), ul(),
                                         m_apar.tau_np_p(), mu, sigma);
    return Li;
}

double RkRdetRnpPdf::PdfAB(const ICPVEvt& evt, bool otlr, bool no_interf) const {
    ReadAndCalc(evt);
    return PdfAB(otlr, no_interf);
}

double RkRdetRnpPdf::PdfAB(bool otlr, bool no_interf) const {
    const double Ef = EfRkRdetRnp_fullrec();
    const double norm_Ef = norm_EfRkRdetRnp_fullrec();
    if (!no_interf) {
        const double Mf = 0.5 / tau() * MfRkRdetRnp_fullrec();
        const double Af = 0.5 / tau() * AfRkRdetRnp_fullrec();
        const double norm_Mf = 0.5 / tau() *norm_MfRkRdetRnp_fullrec();
        const double pdf = Ef * cexp + amix * (C() * Mf - S() * Af);
        const double pdf_norm = norm_Ef * cexp + amix * C() * norm_Mf;
        if (pdf <= 0 || pdf_norm <= 0) {
            cout << "PdfAB. pdf: " << pdf << ", norm: " << pdf_norm
                 << ", cexp: " << cexp << ", amix: " << amix
                 << endl;
            return -fabs(pdf / pdf_norm);
        }
        return pdf / pdf_norm;
//    if(otlr) return AddOutlier(dt,pdf,pdf_norm);
//    else     return pdf/pdf_norm;
    } else {
        return Ef / norm_Ef;
//    if(otlr) return AddOutlier(dt,Ef,norm_Ef);
//    else     return Ef/norm_Ef;
    }
    return -999;
}

// double RkRdetRnpPdf::AddOutlier(double x, double Lin,
// double nLi, double alpha){
//  return Add_Outlier(x,Lin,nLi,alpha);
// }

// double RkRdetRnpPdf::Add_Outlier(double x, double Lin,
//                                  double nLi, double alpha){
//  const double fol = ((m_svd == 2 || ntrk_rec>1) && ntrk_asc>1) ? f_ol_mul :
// f_ol_sgl;
//  const double m = 0.0;
//  const double Lol = gaussian(x,m,sigma_ol);
//  const double nLol = norm_gaussian(ll(),ul(),m,sigma_ol);
//  const double Li = (1.0-fol)*Lin/nLi + fol*alpha*Lol/nLol;
//  return Li;
// }

double RkRdetRnpPdf::operator() (const ICPVEvt& evt) const {
    ReadAndCalc(evt);
    return PdfAB();
}

double RkRdetRnpPdf::operator() (double x) const {
    dt = x;
    return PdfAB();
}

}  // namespace libTatami
