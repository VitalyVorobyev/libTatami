/** Copyright 2016 Vitaly Vorobyev
 ** @file ResVar.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/ResVar.h"

namespace libTatami {

const bool RdetVar::RecSide = true;
const bool RdetVar::AscSide = false;

RdetVar::RdetVar(void) :
    m_ntrk(1), m_sz(1.), m_chisq(1.), m_ndf(0)
{}

RdetVar::RdetVar(const RdetVar& var) {
    *this = var;
}

RdetVar& RdetVar::operator=(const RdetVar& var) {
    m_ntrk  = var.m_ntrk;
    m_sz    = var.m_sz;
    m_chisq = var.m_chisq;
    m_ndf   = var.m_ndf;
    return *this;
}

int RdetVar::ReadVars(const ICPVEvt &evt, const bool type) {
    if (type == RecSide) {
        m_ntrk  = evt.FindIVar("ntrk_rec");
        m_sz    = evt.FindDVar("sz_rec");
        m_chisq = evt.FindDVar("chisq_z_rec");
        m_ndf   = evt.FindIVar("ndf_z_rec");
    } else {
        m_ntrk  = evt.FindIVar("ntrk_asc");
        m_sz    = evt.FindDVar("sz_asc");
        m_chisq = evt.FindDVar("chisq_z_asc");
        m_ndf   = evt.FindIVar("ndf_z_asc");
    }
    return 0;
}

}  // namespace libTatami
