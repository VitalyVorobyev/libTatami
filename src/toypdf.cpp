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
using std::cerr;
using std::endl;
using std::pow;
using std::isnan;

namespace libTatami {

ToyPdf::ToyPdf(const double& m, const double& w, const double& fb,
               const double& wtag, const double& ll, const double& ul) :
  AbsICPVPdf(), m_m(m), m_w(w), m_fbkg(fb) {
    SetWTag(wtag);
    SetRange(ll, ul);
    print_params();
}

ToyPdf::ToyPdf(const ToyPdf& opdf) :
    ToyPdf(opdf.ResMean(), opdf.ResWidth(), opdf.Fbkg()) {
    this->m_ll = opdf.m_ll;
    this->m_ul = opdf.m_ul;
}

double ToyPdf::operator() (const double& dt) {
    const double pdf_s = pdfSig(dt, m_w);
    const double pdf_b = m_fbkg > 0 ? pdfBkg(dt, m_w) : 0;
    const double pdf = (1. - m_fbkg) * pdf_s + m_fbkg * pdf_b;
    if (pdf < 0) {
        cerr << "ToyPdf::operator(): Negative ToyPdf: fbkg " << m_fbkg
             << ", pdfs " << pdf_s
             << ", pdfb " << pdf_b
             << endl;
    }
    return pdf;
}

double ToyPdf::operator() (const double& dt, const int tag) {
    if ((tag != 1) && (tag != -1)) {
        cerr << "Wrong tag value " << tag << endl;
        return -1.;
    }
    SetTag(tag);
    return (*this)(dt);
}

double ToyPdf::operator() (const double& dt, const int tag,
                           const double& c, const double& s) {
    m_c = c; m_s = s;
    SetTag(tag);
    return (*this)(dt);
}

double ToyPdf::operator() (const ICPVEvt& evt) {
    const double dt = evt.DVar("dt");
    return (*this)(dt);
}

double ToyPdf::operator() (const double& dt, const double& fbkg,
                           const double& scale) {
    const double pdf_s =            pdfSig(dt, m_w * scale);
    const double pdf_b = fbkg > 0 ? pdfBkg(dt, m_w * scale) : 0;
    return (1. - fbkg) * pdf_s + fbkg * pdf_b;
}

double ToyPdf::operator() (const double& dt, const double& c, const double& s,
                           const double& fbkg, const double& scale) {
    m_c = c; m_s = s;
    return this->operator()(dt, fbkg, scale);
}

double ToyPdf::pdfSig(const double& dt, const double& wid) {
    double pdf = 0; double norm_pdf = 0;
    if (wid > 0) {
        pdf = Ef_conv_gauss(dt, m_tau, m_m, wid) + m_tag * 0.5 * (1. - 2.*m_wtag) / m_tau * (
              m_c * Mf_conv_gauss(dt, m_tau, m_dm, m_m, wid) +
              m_s * Af_conv_gauss(dt, m_tau, m_dm, m_m, wid));
        if (pdf < 0 || isnan(pdf)) {
//            return 0.;
            cerr << "pdf(" << dt << ") = " << Ef_conv_gauss(dt, m_tau, m_m, wid) << " + " << 0.5 / m_tau << " * ("
                 << "c " << m_c << " * " << Mf_conv_gauss(dt, m_tau, m_dm, m_m, wid) << ") + "
                 << "(s " << m_s << " * " << Af_conv_gauss(dt, m_tau, m_dm, m_m, wid) << ")"
                 << endl;
        }
        norm_pdf = norm_Ef_conv_gauss(m_ll, m_ul, m_tau, m_m, wid) +
                   m_tag * 0.5 / m_tau *
                   m_c * norm_Mf_conv_gauss(m_ll, m_ul, m_tau, m_dm, m_m, wid);
        if (norm_pdf <= 0) {
            cerr << "Negative norm: " << norm_pdf << endl;
            return 0.;
        }
        pdf /= norm_pdf;
    } else {
        pdf = exp(-fabs(dt) / m_tau) * (1. + m_tag * (1. - 2.*m_wtag) *
                                        (m_c * cos(m_dm * dt) + m_s * sin(m_dm * dt)));
        norm_pdf = 2 * m_tau * (1. + m_tag * m_c * (1. - 2.*m_wtag) / (1. + pow(m_tau * m_dm, 2)));
        pdf /= norm_pdf;
    }
    return pdf;
}

double ToyPdf::pdfBkg(const double& dt, const double& wid) {
    if (wid <= 0) return 0;
    return gaussian(dt, m_m, wid) / norm_gaussian(m_ll, m_ul, m_m, wid);
}

void ToyPdf::print_params(void) const {
    cout << "ToyPdf init parameters:" << endl;
    AbsICPVPdf::print_params();
    cout << "  mean: " << m_m << endl
         << "  sigm: " << m_w << endl
         << "  fbkg: " << m_fbkg << endl
         << "  range: " << ll() << " " << ul() << endl;
}

}  // namespace libTatami
