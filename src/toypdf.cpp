/** Copyright 2016 Vitaly Vorobyev
 ** @file toypdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include <iostream>
#include <cmath>

#include "../include/toypdf.h"
#include "../include/typedefs.h"

using std::cout;
using std::endl;
using std::pow;

namespace libTatami {

ToyPdf::ToyPdf(const ToyPdf& opdf)
  : ToyPdf(opdf.ResMean(), opdf.ResWidth(), opdf.Fbkg()) {
    this->m_ll = opdf.m_ll;
    this->m_ul = opdf.m_ul;
}

double ToyPdf::operator() (cdouble& dt) {
    cdouble pdf_s = pdfSig(dt, m_w);
    cdouble pdf_b = m_fbkg > 0 ? pdfBkg(dt, m_w) : 0;
    return (1. - m_fbkg) * pdf_s + m_fbkg * pdf_b;
}

double ToyPdf::operator() (const ICPVEvt& evt) {
    cdouble dt = evt.DVar("dt");
    return (*this)(dt);
}

double ToyPdf::operator() (cdouble& dt, cdouble& fbkg, cdouble& scale) {
    cdouble pdf_s =            pdfSig(dt, m_w * scale);
    cdouble pdf_b = fbkg > 0 ? pdfBkg(dt, m_w * scale) : 0;
    return (1. - fbkg) * pdf_s + fbkg * pdf_b;
}

double ToyPdf::operator() (cdouble& dt, cdouble& c, cdouble& s,
                           cdouble& fbkg, cdouble& scale) {
    m_c = c; m_s = s;
    return this->operator()(dt, fbkg, scale);
}

double ToyPdf::pdfSig(cdouble& dt, cdouble& wid) {
    double pdf = 0; double norm_pdf = 0;
    if (wid > 0) {
        pdf = Ef_conv_gauss(dt, m_tau, m_m, wid) + 0.5 / m_tau * (
                    m_c * Mf_conv_gauss(dt, m_tau, m_dm, m_m, wid) +
                    m_s * Af_conv_gauss(dt, m_tau, m_dm, m_m, wid) );
        norm_pdf = norm_Ef_conv_gauss(m_ll, m_ul, m_tau, m_m, wid) +
                0.5 / m_tau *
                m_c * norm_Mf_conv_gauss(m_ll, m_ul, m_tau, m_dm, m_m, wid);
        pdf /= norm_pdf;
    } else {
        pdf = exp(-fabs(dt) / m_tau) * (1. + m_c * cos(m_dm * dt) +
                                             m_s * sin(m_dm * dt));
        norm_pdf = 2 * m_tau * (1. + 1. / (1. + pow(m_tau * m_dm, 2)));
        pdf /= norm_pdf;
    }
    if (pdf <= 0 || norm_pdf <= 0 || std::isnan(pdf) || pdf > 1000) {
        cout << "ToyPdf::pdfSig: bad pdf: " << pdf << " norm: " <<
                norm_pdf << " dt: " << dt << endl;
        cout << "  tau " << m_tau << ", mean " << m_m << ", w " <<
                wid << ", pdf " << pdf << ", norm_pdf " << norm_pdf << endl;
    }
    return pdf;
}

double ToyPdf::pdfBkg(cdouble& dt, cdouble& wid) {
    if (wid <= 0) return 0;
    return gaussian(dt, m_m, wid) / norm_gaussian(m_ll, m_ul, m_m, wid);
}

}  // namespace libTatami
