/** Copyright 2016 Vitaly Vorobyev
 ** @file Rk.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "Rk.h"

#include "icpvevent.h"

namespace libTatami {

double RkPdf::norm_EfRk(void) const {
    return 0.5 * ((1. - m_pars.r_ckak()) *
                  norm_En(ll(), 0., tau() * (m_pars.ak() - m_pars.ck())) +
                  (1. + m_pars.r_ckak()) *
                  norm_Ep(0., ul(), tau() * (m_pars.ak() + m_pars.ck())));
}

double RkPdf::norm_AfRk(void) const {
    return m_pars.fact_am() *
                          (norm_An(ll(), 0., m_pars.ntau_n(), m_pars.ndm_n()) -
         m_pars.ndmtau() * norm_Mn(ll(), 0., m_pars.ntau_n(), m_pars.ndm_n()) +
                           norm_Ap(0., ul(), m_pars.ntau_p(), m_pars.ndm_p()) -
         m_pars.ndmtau() * norm_Mp(0., ul(), m_pars.ntau_p(), m_pars.ndm_p()) );
}

double RkPdf::norm_MfRk(void) const {
    return m_pars.fact_am() *
                          (norm_Mn(ll(), 0., m_pars.ntau_n(), m_pars.ndm_n()) +
         m_pars.ndmtau() * norm_An(ll(), 0., m_pars.ntau_n(), m_pars.ndm_n()) +
                           norm_Mp(0., ul(), m_pars.ntau_p(), m_pars.ndm_p()) +
         m_pars.ndmtau() * norm_Ap(0., ul(), m_pars.ntau_p(), m_pars.ndm_p()) );
}

double RkPdf::EfRk(double x) const {
    return (x < 0. ? m_pars.fact_n_e() *
                     En(x, tau() * (m_pars.ak() - m_pars.ck())) :
                     m_pars.fact_p_e() *
                     Ep(x, tau() * (m_pars.ak() + m_pars.ck())) );
}

double RkPdf::AfRk(double x) const {
    return (x < 0. ?
                m_pars.fact_am() * (An(x, m_pars.ntau_n(), m_pars.ndm_n()) -
                m_pars.ndmtau()  *  Mn(x, m_pars.ntau_n(), m_pars.ndm_n())) :
                m_pars.fact_am() * (Ap(x, m_pars.ntau_p(), m_pars.ndm_p()) -
                m_pars.ndmtau()  *  Mp(x, m_pars.ntau_p(), m_pars.ndm_p())) );
}

double RkPdf::MfRk(double x) const {
    return (x < 0. ?
                m_pars.fact_am() * (Mn(x, m_pars.ntau_n(), m_pars.ndm_n()) +
                m_pars.ndmtau()  *  An(x, m_pars.ntau_n(), m_pars.ndm_n())) :
                m_pars.fact_am() * (Mp(x, m_pars.ntau_p(), m_pars.ndm_p()) +
                m_pars.ndmtau()  *  Ap(x, m_pars.ntau_p(), m_pars.ndm_p())) );
}

double RkPdf::Pdf(double x) const {
    double pdf      =      EfRk(x);
    double pdf_norm = norm_EfRk();
    if (C() != 0) {
        pdf      += -0.5 / tau() * C() *      MfRk(x);
        pdf_norm += -0.5 / tau() * C() * norm_MfRk();
    }
    if (S() != 0) pdf += 0.5 / tau() * S() * AfRk(x);
    return pdf / pdf_norm;
}

double RkPdf::operator()(double x) const {
    return Pdf(x);
}

double RkPdf::operator()(const ICPVEvt& evt) const {
    const double dz = evt.FindDVar("dz");
    const double costhBcms = evt.FindDVar("costhBcms");
    m_pars.SetAkCk(costhBcms, 0.5*10.58, tau(), dm());
    return Pdf(dz);
}

}  // namespace libTatami
