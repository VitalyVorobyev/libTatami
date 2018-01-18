/** Copyright 2016 Vitaly Vorobyev
 ** @file icpvevent.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#pragma once

#include <string>

#include "icpvvar.h"

namespace libTatami {

///
/// \brief The ICPVEvt class. Class for representation of tuple
/// tuple structure is initialized with text config file
///
class ICPVEvt {
    ///
    /// \brief ReadStructure
    /// \param fname
    /// \return
    ///
    int ReadStructure(const std::string& fname);
    ///
    /// \brief m_IVars
    ///
    ivarvec m_IVars;
    ///
    /// \brief m_DVars
    ///
    dvarvec m_DVars;

 public:
    ///
    /// \brief ICPVEvt
    /// \param fname
    ///
    explicit ICPVEvt(const std::string& fname);
    ///
    /// \brief ICPVEvt
    /// \param ev
    ///
    ICPVEvt(const ICPVEvt& ev);
    ///
    /// \brief ICPVEvt
    /// \param IV
    /// \param DV
    ///
    ICPVEvt(const ivarvec& IV, const dvarvec& DV);
    ///
    /// \brief operator =
    /// \param vt
    /// \return
    ///
    ICPVEvt& operator=(const ICPVEvt& vt);
    ///
    /// \brief Set
    /// \param IV
    /// \param DV
    ///
    void Set(const ivarvec& IV, const dvarvec& DV);
    ///
    /// \brief SetIVar
    /// \param i
    /// \param val
    ///
    void SetIVar(uint32_t i, uint32_t val);
    ///
    /// \brief SetIVar
    /// \param name
    /// \param val
    ///
    void SetIVar(const std::string& name, uint32_t val);
    ///
    /// \brief SetDVar
    /// \param i
    /// \param val
    ///
    void SetDVar(uint32_t i, double val);
    ///
    /// \brief SetDVar
    /// \param name
    /// \param val
    ///
    void SetDVar(const std::string& name, double val);
    ///
    /// \brief IName
    /// \param i
    /// \return
    ///
    std::string IName(uint32_t i) const;
    ///
    /// \brief IVar
    /// \param i
    /// \return
    ///
    int IVar(uint32_t i) const;
    ///
    /// \brief IVar
    /// \param name
    /// \return
    ///
    int IVar(const std::string& name) const;
    ///
    /// \brief DName
    /// \param i
    /// \return
    ///
    std::string DName(uint32_t i) const;
    ///
    /// \brief DVar
    /// \param i
    /// \return
    ///
    double DVar(uint32_t i) const;
    ///
    /// \brief DVar
    /// \param name
    /// \return
    ///
    double DVar(const std::string& name) const;
    ///
    /// \brief FindIVar
    /// \param name
    /// \return
    ///
    int FindIVar(const std::string& name) const;
    ///
    /// \brief FindDVar
    /// \param name
    /// \return
    ///
    int FindDVar(const std::string& name) const;
    ///
    /// \brief IVars
    /// \return
    ///
    const ivarvec& IVars(void) const { return m_IVars;}
    ///
    /// \brief DVars
    /// \return
    ///
    const dvarvec& DVars(void) const { return m_DVars;}
    ///
    /// \brief PrintStructure
    ///
    void PrintStructure(void) const;
};

}  // namespace libTatami
