#include "ttools.h"
#include <cmath>

const double TTools::cm2ps  = 78.48566945838871754705;
const double TTools::beta   = 0.39114064034485169394; // 0.425/sqrt(1.0+0.425^2)
const double TTools::mbzero = 5.2796;

double TTools::sum_sigma(const double& s1, const double& s2){
  return sqrt(s1*s1+s2*s2);
}

void TTools::constraint(double& x, const double& ll, const double& ul){
  if(x<ll){
    x = ll;
  } else if(x>ul){
    x = ul;
  }
  return;
}

void TTools::constraint(double& x, const double& ll){
  if(x<ll) x = ll;
  return;
}
