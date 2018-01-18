/** Copyright 2016 Vitaly Vorobyev
 ** @file icpvvar.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#pragma once

#include <string>
#include <vector>

namespace libTatami {

///
/// \brief The ICPVVar class
///
template <class T> class ICPVVar {
 public:
    ///
    /// \brief ICPVVar
    /// \param v
    /// \param n
    ///
    ICPVVar(const T& v, const std::string& n): val(v), name(n) {}
    ///
    /// \brief SetVal
    /// \param x
    ///
    void SetVal(const T& x) {val = x;}
    ///
    /// \brief Val
    /// \return
    ///
    T Val(void) const {return val;}
    ///
    /// \brief Name
    /// \return
    ///
    std::string Name(void) const {return name;}
    ///
    /// \brief operator =
    /// \param ovar
    /// \return
    ///
    ICPVVar& operator=(const ICPVVar& ovar) {
        this->val  = ovar.val;
        this->name = ovar.name;
        return *this;
    }
    ///
    /// \brief FindIndex
    /// \param v
    /// \param name
    /// \return
    ///
    static int FindIndex(const std::vector<ICPVVar<auto> >& v,
                         const std::string& name) {
        for (unsigned i = 0; i < v.size(); i++) {
            if (name == v[i].Name()) return i;
        }
        return -1;
    }

 private:
    ///
    /// \brief val
    ///
    T val;
    ///
    /// \brief name
    ///
    std::string name;
};

// Aliases
using dvar = ICPVVar<double>;
using ivar = ICPVVar<int>;
using dvarvec = std::vector<dvar>;
using ivarvec = std::vector<ivar> ;

}  // namespace libTatami

