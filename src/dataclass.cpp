/** Copyright 2016 Vitaly Vorobyev
 ** @file dataclass.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "dataclass.h"

namespace libTatami {

DataClass::DataClass(uint16_t svd, bool mc, bool bp) :
    m_svd(svd), m_mc(mc), m_bp(bp) {}

DataClass::DataClass(uint16_t svd, bool mc) :
    DataClass(svd, mc, false) {}

DataClass::DataClass(uint16_t svd) :
    DataClass(svd, false, false) {}

DataClass::DataClass(void) :
    DataClass(0, false, false) {}

}  // namespace libTatami
