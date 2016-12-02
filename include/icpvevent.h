/** Copyright 2016 Vitaly Vorobyev
 ** @file icpvevent.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_ICPVEVENT_H_
#define INCLUDE_ICPVEVENT_H_

#include <string>

#include "./icpvvar.h"

namespace libTatami {

typedef std::string str;

///
/// \brief The ICPVEvt class. Class for representation of tuple
/// tuple structure is initialized with text config file
///
class ICPVEvt {
 public:
    ///
    /// \brief ICPVEvt
    /// \param fname
    ///
    explicit ICPVEvt(const str& fname);
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
    void SetIVar(const unsigned i, const int val);
    ///
    /// \brief SetIVar
    /// \param name
    /// \param val
    ///
    void SetIVar(const str& name, const int val);
    ///
    /// \brief SetDVar
    /// \param i
    /// \param val
    ///
    void SetDVar(const unsigned i, const double& val);
    ///
    /// \brief SetDVar
    /// \param name
    /// \param val
    ///
    void SetDVar(const str& name, const double& val);
    ///
    /// \brief IName
    /// \param i
    /// \return
    ///
    str IName(const unsigned i) const;
    ///
    /// \brief IVar
    /// \param i
    /// \return
    ///
    int IVar(const unsigned i) const;
    ///
    /// \brief IVar
    /// \param name
    /// \return
    ///
    int IVar(const str& name) const;
    ///
    /// \brief DName
    /// \param i
    /// \return
    ///
    str DName(const unsigned i) const;
    ///
    /// \brief DVar
    /// \param i
    /// \return
    ///
    double DVar(const unsigned i) const;
    ///
    /// \brief DVar
    /// \param name
    /// \return
    ///
    double DVar(const str& name) const;
    ///
    /// \brief FindIVar
    /// \param name
    /// \return
    ///
    int FindIVar(const str& name) const;
    ///
    /// \brief FindDVar
    /// \param name
    /// \return
    ///
    int FindDVar(const str& name) const;
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

 private:
    ///
    /// \brief ReadStructure
    /// \param fname
    /// \return
    ///
    int ReadStructure(const str& fname);
    ///
    /// \brief m_IVars
    ///
    ivarvec m_IVars;
    ///
    /// \brief m_DVars
    ///
    dvarvec m_DVars;
};

}  // namespace libTatami

#endif  // INCLUDE_ICPVEVENT_H_
