/** Copyright 2016 Vitaly Vorobyev
 ** @file fenp.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "fenp.h"

namespace libTatami {

void fEnp::Clear() {
    fEn = 0.;
    fEp = 0.;
    fEn_np = 0.;
    fEp_np = 0.;
    fxEn = 0.;
    fxEp = 0.;
    fAn = 0.;
    fAp = 0.;
    fMn = 0.;
    fMp = 0.;
    fxEn_k = 0.;
    fxEp_k = 0.;
    fEn_k = 0.;
    fEp_k = 0.;
}

}  // namespace libTatami
