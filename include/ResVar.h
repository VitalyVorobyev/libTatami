/** Copyright 2016 Vitaly Vorobyev
 ** @file ResVar.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 */

/** Copyright 2016 Vitaly Vorobyev
 * @file ResVar.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#ifndef INCLUDE_RESVAR_H_
#define INCLUDE_RESVAR_H_

#include "./icpvevent.h"

namespace libTatami {

///
/// \brief The RdetVar class. Event-dependent variables for vertex resilution
///
class RdetVar {
 public:
    ///
    /// \brief RdetVar
    ///
    RdetVar(void);
    ///
    /// \brief RdetVar
    /// \param var
    ///
    RdetVar(const RdetVar& var);
    ///
    /// \brief operator =
    /// \param var
    /// \return
    ///
    RdetVar& operator= (const RdetVar& var);
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
    void Set_ntrk(const int v) { m_ntrk = v;}
    ///
    /// \brief Set_sz
    /// \param v
    ///
    void Set_sz(const double& v) { m_sz = v;}
    ///
    /// \brief Set_chisq
    /// \param v
    ///
    void Set_chisq(const double& v) { m_chisq = v;}
    ///
    /// \brief Set_ndf
    /// \param v
    ///
    void Set_ndf(const int v) { m_ndf = v;}
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

 private:
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
};

}  // namespace libTatami

#endif  // INCLUDE_RESVAR_H_
