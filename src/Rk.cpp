/** Copyright 2016 Vitaly Vorobyev
 ** @file Rk.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/Rk.h"

namespace libTatami {

double RkPdf::norm_EfRk(void) const {
    return 0.5 * ((1. - m_pars.r_ckak()) *
                  norm_En(m_ll, 0., m_tau * (m_pars.ak() - m_pars.ck())) +
                  (1. + m_pars.r_ckak()) *
                  norm_Ep(0., m_ul, m_tau * (m_pars.ak() + m_pars.ck())));
}

double RkPdf::norm_AfRk(void) const {
    return m_pars.fact_am() *
                          (norm_An(m_ll, 0., m_pars.ntau_n(), m_pars.ndm_n()) -
         m_pars.ndmtau() * norm_Mn(m_ll, 0., m_pars.ntau_n(), m_pars.ndm_n()) +
                           norm_Ap(0., m_ul, m_pars.ntau_p(), m_pars.ndm_p()) -
         m_pars.ndmtau() * norm_Mp(0., m_ul, m_pars.ntau_p(), m_pars.ndm_p()) );
}

double RkPdf::norm_MfRk(void) const {
    return m_pars.fact_am() *
                          (norm_Mn(m_ll, 0., m_pars.ntau_n(), m_pars.ndm_n()) +
         m_pars.ndmtau() * norm_An(m_ll, 0., m_pars.ntau_n(), m_pars.ndm_n()) +
                           norm_Mp(0., m_ul, m_pars.ntau_p(), m_pars.ndm_p()) +
         m_pars.ndmtau() * norm_Ap(0., m_ul, m_pars.ntau_p(), m_pars.ndm_p()) );
}

double RkPdf::EfRk(const double& x) const {
    return (x < 0. ? m_pars.fact_n_e() *
                     En(x, m_tau * (m_pars.ak() - m_pars.ck())) :
                     m_pars.fact_p_e() *
                     Ep(x, m_tau * (m_pars.ak() + m_pars.ck())) );
}

double RkPdf::AfRk(const double& x) const {
    return (x < 0. ?
                m_pars.fact_am() * (An(x, m_pars.ntau_n(), m_pars.ndm_n()) -
                m_pars.ndmtau()  *  Mn(x, m_pars.ntau_n(), m_pars.ndm_n())) :
                m_pars.fact_am() * (Ap(x, m_pars.ntau_p(), m_pars.ndm_p()) -
                m_pars.ndmtau()  *  Mp(x, m_pars.ntau_p(), m_pars.ndm_p())) );
}

double RkPdf::MfRk(const double& x) const {
    return (x < 0. ?
                m_pars.fact_am() * (Mn(x, m_pars.ntau_n(), m_pars.ndm_n()) +
                m_pars.ndmtau()  *  An(x, m_pars.ntau_n(), m_pars.ndm_n())) :
                m_pars.fact_am() * (Mp(x, m_pars.ntau_p(), m_pars.ndm_p()) +
                m_pars.ndmtau()  *  Ap(x, m_pars.ntau_p(), m_pars.ndm_p())) );
}

double RkPdf::Pdf(const double& x) const {
    double pdf      =      EfRk(x);
    double pdf_norm = norm_EfRk();
    if (m_c != 0) {
        pdf      += -0.5 / m_tau * m_c *      MfRk(x);
        pdf_norm += -0.5 / m_tau * m_c * norm_MfRk();
    }
    if (m_s != 0) pdf += 0.5 / m_tau * m_s * AfRk(x);
    return pdf / pdf_norm;
}

double RkPdf::operator()(const double& x) const {
    return Pdf(x);
}

double RkPdf::operator()(const ICPVEvt& evt) {
    const double dz        = evt.FindDVar("dz");
    const double costhBcms = evt.FindDVar("costhBcms");
    m_pars.SetAkCk(costhBcms, 0.5*10.58, m_tau, m_dm);
    return Pdf(dz);
}

}  // namespace libTatami
