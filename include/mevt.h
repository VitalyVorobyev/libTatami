/** Copyright 2016 Vitaly Vorobyev
 * @file mevt.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#pragma once

#include <map>

namespace libTatami {

///
/// \brief The MEvt class
/// \todo Not used class yet
///
class MEvt {
    // Aliases
    using imap = std::map<std::string, int>;
    using dmap = std::map<std::string, double>;
    using ivar = std::pair<std::string, int>;
    using dvar = std::pair<std::string, double>;
    ///
    /// \brief ReadStructure
    /// \param fname
    /// \return
    ///
    int ReadStructure(const std::string& fname);
    ///
    /// \brief m_ivars
    ///
    imap m_ivars;
    ///
    /// \brief m_dvars
    ///
    dmap m_dvars;

 public:
    ///
    /// \brief MEvt
    /// \param fname
    ///
    explicit MEvt(const std::string& fname);
    ///
    /// \brief IVal
    /// \param name
    /// \return
    ///
    int IVal(const std::string& name) const { return m_ivars.at(name); }
    ///
    /// \brief DVal
    /// \param name
    /// \return
    ///
    double DVal(const std::string& name) { return m_dvars[name]; }
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
};

}  // namespace libTatami
