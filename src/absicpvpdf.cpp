/** Copyright 2016 Vitaly Vorobyev
 ** @file parmanager.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/absicpvpdf.h"

#include <iostream>

using std::cout;
using std::endl;

namespace libTatami {

AbsICPVPdf::AbsICPVPdf(double tau, double dm, double c, double s) :
    m_tau(tau), m_dm(dm), m_c(c), m_s(s), m_tag(1), m_wtag(0.) {}

void AbsICPVPdf::print_params(void) const {
    cout << "   tau: " << m_tau << endl
         << "    dm: " << m_dm << endl
         << "  wtag: " << m_wtag << endl;
}

}  // namespace libTatami
