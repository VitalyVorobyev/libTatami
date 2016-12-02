/** Copyright 2016 Vitaly Vorobyev
 ** @file RbkgPdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include <fstream>

#include "../include/RbkgPdf.h"
#include "../include/typedefs.h"

using std::cout;
using std::endl;

namespace libTatami {

cdouble RbkgPdf::cm2ps = 78.48566945838871754705;

RbkgPdf::RbkgPdf(cstr& fname) :
    AbsPdf(), m_scale(1), m_shift(0),
    m_sigma(1.), m_ndf(1) {
    m_pars.GetParametersFromFile(fname);
}

RbkgPdf::RbkgPdf(void) :
    RbkgPdf("../params/def_bkg.txt")
{}

double RbkgPdf::operator()(cdouble& x) {
    return Pdf(x, m_sigma, m_ndf);
}

double RbkgPdf::operator()(const ICPVEvt& evt) {
    cdouble dt = evt.FindDVar("dt");
    cdouble s  = evt.FindDVar("s");
    const int ndf = evt.FindIVar("ndf");
    return Pdf(dt, s, ndf);
}

double RbkgPdf::Pdf(const double &x, const double &s, const int ndf) {
    double smain = ndf ? s * m_pars.S_main_mlt() * cm2ps :
                         s * m_pars.S_main_sgl() * cm2ps;
    double stail = ndf ? smain * m_pars.S_tail_mlt() :
                         smain * m_pars.S_tail_sgl();
    cdouble f_tail  = ndf ? m_pars.f_tail_mlt() :
                            m_pars.f_tail_sgl();
    cdouble f_delta = ndf ? m_pars.f_delta_mlt() :
                            m_pars.f_delta_sgl();
    cdouble tau = m_pars.tau();
    cdouble mu = m_pars.mu() + m_shift;
    cdouble mu_delta = m_pars.mu_delta() + m_shift;
    smain *= m_scale;
    stail *= m_scale;

    cdouble pdf_l = f_tail  * Enp_conv_gauss(x, tau, tau, mu, stail) +
               (1 - f_tail) * Enp_conv_gauss(x, tau, tau, mu, smain);
    cdouble int_pdf_l =
               f_tail  * norm_Enp_conv_gauss(m_ll, m_ul, tau, tau, mu, stail) +
          (1 - f_tail) * norm_Enp_conv_gauss(m_ll, m_ul, tau, tau, mu, smain);


    cdouble pdf_d = f_tail  * gaussian(x, mu_delta, stail) +
               (1 - f_tail) * gaussian(x, mu_delta, smain);
    cdouble int_pdf_d = f_tail  * norm_gaussian(m_ll, m_ul, mu_delta, stail) +
                   (1 - f_tail) * norm_gaussian(m_ll, m_ul, mu_delta, smain);

    cdouble pdf = f_delta  * pdf_d +
             (1 - f_delta) * pdf_l;
    cdouble int_pdf = f_delta  * int_pdf_d +
                 (1 - f_delta) * int_pdf_l;
    if (pdf >= 0 && int_pdf >= 0) {
        if (m_pars.f_otlr() > 0.0001) return AddOutlier(x, pdf, int_pdf);
        else
            return pdf / int_pdf;
    } else {
        cout << "RbkgPdf::Pdf: " << pdf << ", norm = " << int_pdf << endl;
        return 0;
    }
}

double RbkgPdf::AddOutlier(const double& x, const double Lin,
                           const double& nLi) {
    cdouble Lol  = gaussian(x, 0., m_pars.s_otlr());
    cdouble nLol = norm_gaussian(m_ll, m_ul, 0., m_pars.s_otlr());
    cdouble Li   = (1. - m_pars.f_otlr() ) * Lin / nLi +
                         m_pars.f_otlr()   * Lol / nLol;
    return Li;
}

}  // namespace libTatami
