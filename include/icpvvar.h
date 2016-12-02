/** Copyright 2016 Vitaly Vorobyev
 ** @file icpvvar.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_ICPVVAR_H_
#define INCLUDE_ICPVVAR_H_

#include <string>
#include <vector>
// #include <algorithm>
// #include <map>

namespace libTatami {

typedef std::string str;

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
    ICPVVar(const T& v, const str& n): val(v), name(n) {}
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
    str Name(void) const {return name;}
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
                         const str& name) {
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
    str name;
};

///
/// \brief dvar
///
typedef ICPVVar<double> dvar;

///
/// \brief ivar
///
typedef ICPVVar<int> ivar;

///
/// \brief dvarvec
///
typedef std::vector<dvar> dvarvec;
// typedef std::map<std::string, double> dvarvec;

///
/// \brief ivarvec
///
typedef std::vector<ivar> ivarvec;
// typedef std::map<std::string, int> ivarvec;

}  // namespace libTatami

#endif  // INCLUDE_ICPVVAR_H_
