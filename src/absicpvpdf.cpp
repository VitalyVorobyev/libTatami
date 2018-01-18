/** Copyright 2016 Vitaly Vorobyev
 ** @file parmanager.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "absicpvpdf.h"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

constexpr bool dump = false;

namespace libTatami {

AbsICPVPdf::AbsICPVPdf(double tau, double dm, double c, double s) :
    m_tau(tau), m_dm(dm), m_c(c), m_s(s), m_tag(1), m_wtag(0.) {}

void AbsICPVPdf::print_params(void) const {
    cout << "   tau: " << m_tau << endl
         << "    dm: " << m_dm << endl
         << "  wtag: " << m_wtag << endl;
}

void AbsICPVPdf::SetCS(double c, double s) const {
    m_c = c; m_s = s;
    if (dump)
        cout << "C: " << m_c << ", S: " << m_s << endl;
}

void AbsICPVPdf::SetCS(const std::pair<double, double>& x) const {
    m_c = x.first; m_s = x.second;
    if (dump)
        cout << "SetCS: C: " << m_c << ", S: " << m_s << endl;
}

double AbsICPVPdf::alpha() const {
    return (1. - 2. * m_wtag) / (1. + std::pow(m_dm * m_tau, 2));
}

}  // namespace libTatami
