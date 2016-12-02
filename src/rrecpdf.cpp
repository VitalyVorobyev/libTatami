/** Copyright 2016 Vitaly Vorobyev
 ** @file rrecpdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/rrecpdf.h"
#include "../include/ttools.h"

namespace libTatami {

double RrecPdf::Rrec(const ICPVEvt& evt) {
    m_vars.ReadVars(evt, RdetVar::RecSide);
    m_pars.Calculate(m_cnst, m_vars);
    return Rrec();
}

double RrecPdf::Rrec(void) const {
    cdouble Li_mn = gaussian(dz, m_pars.mu_main_rec(), m_pars.Smain_rec());
    if (m_pars.ftail_rec() > 0.) { /* Mainly single track case */
        cdouble Li_tl = gaussian(dz, m_pars.mu_tail_rec(), m_pars.Stail_rec());
        return (1. - m_pars.ftail_rec()) * Li_mn + m_pars.ftail_rec() * Li_tl;
    }
    return Li_mn;
}

double RrecPdf::norm_Rrec(const ICPVEvt &evt) {
    m_vars.ReadVars(evt, RdetVar::RecSide);
    m_pars.Calculate(m_cnst, m_vars);
    return norm_Rrec();
}

double RrecPdf::norm_Rrec(void) const {
    cdouble Li_mn = norm_gaussian(m_ll, m_ul, m_pars.mu_main_rec(),
                                              m_pars.Smain_rec());
    if (m_pars.ftail_rec() > 0.) {  // Mainly multiple track case
        cdouble Li_tl = norm_gaussian(m_ll, m_ul, m_pars.mu_tail_rec(),
                                                  m_pars.Stail_rec());
        return (1. - m_pars.ftail_rec()) * Li_mn + m_pars.ftail_rec() * Li_tl;
    }
    return Li_mn;
}

double RrecPdf::operator() (cdouble& x) {
    dz = x;
    return Pdf();
}

double RrecPdf::Pdf(const ICPVEvt& evt) {
    m_vars.ReadVars(evt, RdetVar::RecSide);
    m_pars.Calculate(m_cnst, m_vars);
    dz = evt.FindDVar("dz");
    return Pdf();
}

double RrecPdf::Pdf(void) const {
    return Rrec() / norm_Rrec();
}

double RrecPdf::operator() (const ICPVEvt& evt) {
    return Pdf(evt);
}

}  // namespace libTatami
