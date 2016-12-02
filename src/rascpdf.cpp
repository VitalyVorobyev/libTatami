/** Copyright 2016 Vitaly Vorobyev
 ** @file rascpdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/rascpdf.h"

namespace libTatami {

int RascPdf::ReadVars(const ICPVEvt& evt) {
    m_vars.ReadVars(evt, RdetVar::AscSide);
    dz       = evt.DVar("dz_asc");
    keeptagl = evt.IVar("keeptagl");
    return 0;
}

double RascPdf::operator()(const ICPVEvt& evt) {
    return PdfRascRnp(evt);
}

double RascPdf::operator()(const double& x) {
    dz = x;
    return PdfRascRnp();
}

double RascPdf::RascRnp(const ICPVEvt &evt) {
    ReadVars(evt); m_pars.Calculate(m_cnst, m_vars, keeptagl);
    return RascRnp();
}

double RascPdf::RascRnp(void) const {
    cdouble Li_md = gaussian(dz, m_pars.mu_main_asc(), m_pars.Smain_asc());
    cdouble Li_mp = Ep_conv_gauss(dz, m_pars.tau_np_p(),
                                  m_pars.mu_main_asc(), m_pars.Smain_asc());
    cdouble Li_mn = En_conv_gauss(dz, m_pars.tau_np_n(),
                                  m_pars.mu_main_asc(), m_pars.Smain_asc());
    cdouble Li_me = m_pars.fp()*Li_mp + (1. - m_pars.fp()) * Li_mn;
    cdouble Li_mt = m_pars.fd()*Li_md + (1. - m_pars.fd()) * Li_me;
    if (m_pars.ftail_asc() == 0.) {  /* Mainly multiple track case */
        return Li_mt;
    }
    cdouble Li_td = gaussian(dz, m_pars.mu_tail_asc(), m_pars.Stail_asc());
    cdouble Li_tp = Ep_conv_gauss(dz, m_pars.tau_np_p_tl(),
                                  m_pars.mu_tail_asc(), m_pars.Stail_asc());
    cdouble Li_tn = En_conv_gauss(dz, m_pars.tau_np_n_tl(),
                                  m_pars.mu_tail_asc(), m_pars.Stail_asc());
    cdouble Li_te = m_pars.fp() * Li_tp + (1. - m_pars.fp()) * Li_tn;
    cdouble Li_tt = m_pars.fd() * Li_td + (1. - m_pars.fd()) * Li_te;
    double Li  = (1. - m_pars.ftail_asc()) * Li_mt +
                       m_pars.ftail_asc()  * Li_tt;
    return Li;
}

double RascPdf::norm_RascRnp(const ICPVEvt& evt) {
      ReadVars(evt); m_pars.Calculate(m_cnst, m_vars, keeptagl);
    return norm_RascRnp();
}

double RascPdf::norm_RascRnp(void) const {
    cdouble Li_md = norm_gaussian(m_ll, m_ul, m_pars.mu_main_asc(),
                                  m_pars.Smain_asc());
    cdouble Li_mp = norm_Ep_conv_gauss(m_ll, m_ul, m_pars.tau_np_p(),
                                     m_pars.mu_main_asc(), m_pars.Smain_asc());
    cdouble Li_mn = norm_En_conv_gauss(m_ll, m_ul, m_pars.tau_np_n(),
                                     m_pars.mu_main_asc(), m_pars.Smain_asc());
    cdouble Li_me = m_pars.fp() * Li_mp + (1. - m_pars.fp()) * Li_mn;
    cdouble Li_mt = m_pars.fd() * Li_md + (1. - m_pars.fd()) * Li_me;
    if (m_pars.ftail_asc() == 0.) { /* Mainly multiple track case */
        return Li_mt;
    }
    cdouble Li_td = norm_gaussian(m_ll, m_ul, m_pars.mu_tail_asc(),
                                  m_pars.Stail_asc());
    cdouble Li_tp = norm_Ep_conv_gauss(m_ll, m_ul, m_pars.tau_np_p_tl(),
                                     m_pars.mu_tail_asc(), m_pars.Stail_asc());
    cdouble Li_tn = norm_En_conv_gauss(m_ll, m_ul, m_pars.tau_np_n_tl(),
                                     m_pars.mu_tail_asc(), m_pars.Stail_asc());
    cdouble Li_te = m_pars.fp() * Li_tp + (1. - m_pars.fp()) * Li_tn;
    cdouble Li_tt = m_pars.fd() * Li_td + (1. - m_pars.fd()) * Li_te;
    double Li  = (1. - m_pars.ftail_asc()) * Li_mt +
                       m_pars.ftail_asc()  * Li_tt;
    return Li;
}

double RascPdf::PdfRascRnp(const ICPVEvt& evt) {
    ReadVars(evt); m_pars.Calculate(m_cnst, m_vars, keeptagl);
    return PdfRascRnp();
}

double RascPdf::PdfRascRnp(void) const {
    cdouble      pdf =      RascRnp();
    cdouble norm_pdf = norm_RascRnp();
    return pdf / norm_pdf;
//    return AddOutlier(dz_asc,pdf,norm_pdf);
}

}  // namespace libTatami
