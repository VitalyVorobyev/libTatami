/** Copyright 2016 Vitaly Vorobyev
 ** @file ttools.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "ttools.h"

#include <cmath>

namespace libTatami {

const double TTools::cm2ps  = 78.48566945838871754705;
const double TTools::beta   = 0.39114064034485169394;  // 0.425/sqrt(1.+0.425^2)
const double TTools::mbzero = 5.2796;

double TTools::sum_sigma(double s1, double s2) {
  return std::sqrt(s1 * s1 + s2 * s2);
}

void TTools::constraint(double& x, double ll, double ul) {
    if (x < ll) {
        x = ll;
    } else if (x > ul) {
        x = ul;
    }
}

void TTools::constraint(double& x, double ll) {
    if (x < ll) x = ll;
}

}  // namespace libTatami
