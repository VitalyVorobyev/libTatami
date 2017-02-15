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

AbsICPVPdf::AbsICPVPdf(const double& c, const double& s,
                       const double& tau, const double& dm) :
    m_c(c), m_s(s), m_tau(tau), m_dm(dm)
{}

AbsICPVPdf::AbsICPVPdf(void) :
    AbsICPVPdf(1., 0., 1.520, 0.505)
{}

AbsICPVPdf::AbsICPVPdf(const double& tau, const double& dm) :
    AbsICPVPdf(1., 0., tau, dm)
{}

void AbsICPVPdf::print_params(void) const {
    cout << "  tau: " << m_tau << endl;
    cout << "  dm: " << m_dm << endl;
}

}  // namespace libTatami
