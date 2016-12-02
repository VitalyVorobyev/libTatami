/** Copyright 2016 Vitaly Vorobyev
 **
 ** @file typedefs.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_TYPEDEFS_H_
#define INCLUDE_TYPEDEFS_H_

#include <vector>
#include <string>
#include <map>
// #include <utility>

namespace libTatami {

typedef std::string str;

typedef const str cstr;
typedef const double cdouble;
typedef const int cint;

typedef std::vector<double> vectd;
typedef std::vector<int>    vecti;

typedef std::map<str, int>    imap;
typedef std::map<str, double> dmap;

// typedef std::pair<str, int>    ivar;
// typedef std::pair<str, double> dvar;

}  // namespace libTatami

#endif  // INCLUDE_TYPEDEFS_H_
