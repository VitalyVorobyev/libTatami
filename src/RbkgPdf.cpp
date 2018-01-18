/** Copyright 2016 Vitaly Vorobyev
 ** @file RbkgPdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "RbkgPdf.h"

#include <fstream>

#include "icpvevent.h"

using std::cout;
using std::endl;

namespace libTatami {

const double RbkgPdf::cm2ps = 78.48566945838871754705;

RbkgPdf::RbkgPdf(const std::string& fname) :
    AbsPdf(), m_pars(fname), m_scale(1), m_shift(0),
    m_sigma(1.), m_ndf(1) {}

RbkgPdf::RbkgPdf(void) : RbkgPdf("../params/def_bkg.txt") {}

double RbkgPdf::operator()(double x) const {
    return Pdf(x, m_sigma, m_ndf);
}

double RbkgPdf::operator()(const ICPVEvt& evt) const {
    const double dt = evt.FindDVar("dt");
    const double s  = evt.FindDVar("s");
    int ndf = evt.FindIVar("ndf");
    return Pdf(dt, s, ndf);
}

double RbkgPdf::Pdf(double x, double s, int ndf) const {
    double smain = ndf ? s * m_pars.S_main_mlt() * cm2ps :
                         s * m_pars.S_main_sgl() * cm2ps;
    double stail = ndf ? smain * m_pars.S_tail_mlt() :
                         smain * m_pars.S_tail_sgl();
    const double f_tail  = ndf ? m_pars.f_tail_mlt() :
                            m_pars.f_tail_sgl();
    const double f_delta = ndf ? m_pars.f_delta_mlt() :
                            m_pars.f_delta_sgl();
    const double tau = m_pars.tau();
    const double mu = m_pars.mu() + m_shift;
    const double mu_delta = m_pars.mu_delta() + m_shift;
    smain *= m_scale;
    stail *= m_scale;

    const double pdf_l = f_tail  * Enp_conv_gauss(x, tau, tau, mu, stail) +
               (1 - f_tail) * Enp_conv_gauss(x, tau, tau, mu, smain);
    const double int_pdf_l =
               f_tail  * norm_Enp_conv_gauss(ll(), ul(), tau, tau, mu, stail) +
          (1 - f_tail) * norm_Enp_conv_gauss(ll(), ul(), tau, tau, mu, smain);


    const double pdf_d = f_tail  * gaussian(x, mu_delta, stail) +
               (1 - f_tail) * gaussian(x, mu_delta, smain);
    const double int_pdf_d = f_tail  * norm_gaussian(ll(), ul(), mu_delta, stail) +
                   (1 - f_tail) * norm_gaussian(ll(), ul(), mu_delta, smain);

    const double pdf = f_delta  * pdf_d +
             (1 - f_delta) * pdf_l;
    const double int_pdf = f_delta  * int_pdf_d +
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

double RbkgPdf::AddOutlier(double x, double Lin, double nLi) const {
    const double Lol  = gaussian(x, 0., m_pars.s_otlr());
    const double nLol = norm_gaussian(ll(), ul(), 0., m_pars.s_otlr());
    const double Li   = (1. - m_pars.f_otlr() ) * Lin / nLi +
                         m_pars.f_otlr()   * Lol / nLol;
    return Li;
}

}  // namespace libTatami
