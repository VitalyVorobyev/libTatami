/** Copyright 2016 Vitaly Vorobyev
 ** @file rrecpdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "rrecpdf.h"
#include "ttools.h"
#include "parmanager.h"
#include "icpvevent.h"

namespace libTatami {

RrecPdf::RrecPdf(const DataClass &dc) :
    AbsPdf(), m_cnst(ParManager::SigParFile(dc)) {}

double RrecPdf::Rrec(const ICPVEvt& evt) {
    m_vars.ReadVars(evt, RdetVar::RecSide);
    m_pars.Calculate(m_cnst, m_vars);
    return Rrec();
}

double RrecPdf::Rrec(void) const {
    const double Li_mn = gaussian(dz, m_pars.mu_main_rec(),
                                  m_pars.Smain_rec());
    if (m_pars.ftail_rec() > 0.) { /* Mainly single track case */
        const double Li_tl = gaussian(dz, m_pars.mu_tail_rec(),
                                      m_pars.Stail_rec());
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
    const double Li_mn = norm_gaussian(ll(), ul(), m_pars.mu_main_rec(),
                                       m_pars.Smain_rec());
    if (m_pars.ftail_rec() > 0.) {  // Mainly multiple track case
        const double Li_tl = norm_gaussian(ll(), ul(), m_pars.mu_tail_rec(),
                                           m_pars.Stail_rec());
        return (1. - m_pars.ftail_rec()) * Li_mn + m_pars.ftail_rec() * Li_tl;
    }
    return Li_mn;
}

double RrecPdf::operator() (double x) const {
    dz = x;
    return Pdf();
}

double RrecPdf::Pdf(const ICPVEvt& evt) const {
    m_vars.ReadVars(evt, RdetVar::RecSide);
    m_pars.Calculate(m_cnst, m_vars);
    dz = evt.FindDVar("dz");
    return Pdf();
}

double RrecPdf::Pdf(void) const {
    return Rrec() / norm_Rrec();
}

double RrecPdf::operator() (const ICPVEvt& evt) const {
    return Pdf(evt);
}

}  // namespace libTatami
