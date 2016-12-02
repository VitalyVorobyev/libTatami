/** Copyright 2016 Vitaly Vorobyev
 * @file mevt.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#ifndef INCLUDE_MEVT_H_
#define INCLUDE_MEVT_H_

#include "./typedefs.h"

namespace libTatami {

///
/// \brief The MEvt class
///
class MEvt {
 public:
    ///
    /// \brief MEvt
    /// \param fname
    ///
    explicit MEvt(const str& fname);
    ///
    /// \brief IVal
    /// \param name
    /// \return
    ///
    int IVal(const str& name) const { return m_ivars.at(name); }
    ///
    /// \brief DVal
    /// \param name
    /// \return
    ///
    double DVal(const str& name) { return m_dvars[name]; }
    ///
    /// \brief IVars
    /// \return
    ///
    const imap& IVars(void) const { return m_ivars;}
    ///
    /// \brief DVars
    /// \return
    ///
    const dmap& DVars(void) const { return m_dvars;}
    ///
    /// \brief PrintStructure
    ///
    void PrintStructure(void) const;

 private:
    ///
    /// \brief ReadStructure
    /// \param fname
    /// \return
    ///
    int ReadStructure(const str& fname);
    ///
    /// \brief m_ivars
    ///
    imap m_ivars;
    ///
    /// \brief m_dvars
    ///
    dmap m_dvars;
};

}  // namespace libTatami

#endif  // INCLUDE_MEVT_H_
