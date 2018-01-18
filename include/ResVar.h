/** Copyright 2016 Vitaly Vorobyev
 ** @file ResVar.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 */

#pragma once

#include <cstdint>

namespace libTatami {

class ICPVEvt;

///
/// \brief The RdetVar class. Event-dependent variables for vertex resilution
///
class RdetVar {
    ///
    /// \brief m_ntrk
    ///
    int m_ntrk;
    ///
    /// \brief m_sz
    ///
    double m_sz;
    ///
    /// \brief m_chisq
    ///
    double m_chisq;
    ///
    /// \brief m_ndf
    ///
    int m_ndf;

 public:
    ///
    /// \brief RdetVar
    ///
    RdetVar(void);
    ///
    /// \brief ReadVars
    /// \param evt
    /// \param type
    /// \return
    ///
    int ReadVars(const ICPVEvt &evt, const bool type);
    ///
    /// \brief Set_ntrk
    /// \param v
    ///
    void Set_ntrk(uint16_t v) { m_ntrk = v;}
    ///
    /// \brief Set_sz
    /// \param v
    ///
    void Set_sz(double v) { m_sz = v;}
    ///
    /// \brief Set_chisq
    /// \param v
    ///
    void Set_chisq(double v) { m_chisq = v;}
    ///
    /// \brief Set_ndf
    /// \param v
    ///
    void Set_ndf(uint16_t v) { m_ndf = v;}
    ///
    /// \brief ntrk
    /// \return
    ///
    int ntrk(void) const { return m_ntrk;}
    ///
    /// \brief sz
    /// \return
    ///
    double sz(void) const { return m_sz;}
    ///
    /// \brief chisq
    /// \return
    ///
    double chisq(void) const { return m_chisq;}
    ///
    /// \brief ndf
    /// \return
    ///
    int ndf(void) const {return m_ndf;}
    ///
    /// \brief RecSide
    ///
    static const bool RecSide;
    ///
    /// \brief AscSide
    ///
    static const bool AscSide;
};

}  // namespace libTatami
