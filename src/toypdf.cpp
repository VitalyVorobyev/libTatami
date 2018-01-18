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

#include "toypdf.h"
#include "typedefs.h"
#include "icpvevent.h"

using std::cout;
using std::cerr;
using std::endl;
using std::pow;
using std::isnan;

namespace libTatami {

ToyPdf::ToyPdf(double m, double w, double fb,
               double wtag, double ll, double ul) :
    AbsICPVPdf(), m_m(m), m_w(w), m_fbkg(fb) {
    SetWTag(wtag);
    SetRange(ll, ul);
}

double ToyPdf::operator() (double dt) const {
    const double pdf_s = pdfSig(dt, m_w);
    const double pdf_b = m_fbkg > 0 ? pdfBkg(dt, m_w) : 0;
    const double pdf = (1. - m_fbkg) * pdf_s + m_fbkg * pdf_b;
    if (pdf < 0) {
        cerr << "ToyPdf::operator(): Negative ToyPdf: fbkg: " << m_fbkg
             << ", pdfs: " << pdf_s
             << ", pdfb: " << pdf_b
             << ", dt: " << dt
             << endl;
    }
    return pdf;
}

double ToyPdf::operator() (double dt, int32_t tag) {
    if ((tag != 1) && (tag != -1)) {
        cerr << "Wrong tag value " << tag << endl;
        return -1.;
    }
    SetTag(tag);
    return (*this)(dt);
}

double ToyPdf::operator() (double dt, int32_t tag, double c, double s) const {
    SetCS(c, s);
    SetTag(tag);
    return (*this)(dt);
}

double ToyPdf::operator() (const ICPVEvt& evt) const {
    const double dt = evt.DVar("dt");
    return (*this)(dt);
}

double ToyPdf::operator() (double dt, double fbkg, double scale) const {
    const double pdf_s =            pdfSig(dt, m_w * scale);
    const double pdf_b = fbkg > 0 ? pdfBkg(dt, m_w * scale) : 0;
    return (1. - fbkg) * pdf_s + fbkg * pdf_b;
}

double ToyPdf::operator() (double dt, double c, double s,
                           double fbkg, double scale) const {
    SetCS(c, s);
    return this->operator()(dt, fbkg, scale);
}

double ToyPdf::pdfSig(double dt, double wid) const {
    double pdf = 0; double norm_pdf = 0;
    if (wid > 0) {
        pdf = Ef_conv_gauss(dt, tau(), m_m, wid) +
                tag() * 0.5 * (1. - 2.*wtag()) / tau() * (
              C() * Mf_conv_gauss(dt, tau(), dm(), m_m, wid) +
              S() * Af_conv_gauss(dt, tau(), dm(), m_m, wid));
        if (pdf < 0 || isnan(pdf)) {
            cerr << "pdf(" << dt << ") = "
                 << Ef_conv_gauss(dt, tau(), m_m, wid)
                 << " + " << 0.5 / tau() << " * ("
                 << "c " << C() << " * "
                 << Mf_conv_gauss(dt, tau(), dm(), m_m, wid) << ") + "
                 << "(s " << S() << " * "
                 << Af_conv_gauss(dt, tau(), dm(), m_m, wid) << ")"
                 << endl;
        }
        norm_pdf = norm_Ef_conv_gauss(ll(), ul(), tau(), m_m, wid) +
                   tag() * 0.5 / tau() *
                   C() * norm_Mf_conv_gauss(ll(), ul(), tau(), dm(), m_m, wid);
        if (norm_pdf <= 0) {
            cerr << "Negative norm: " << norm_pdf << endl;
            return 0.;
        }
        pdf /= norm_pdf;
    } else {
        pdf = exp(-fabs(dt) / tau()) * (1. + tag() * (1. - 2.*wtag()) *
                 (C() * cos(dm() * dt) + S() * sin(dm() * dt)));
        norm_pdf = 2 * tau() * (1. + tag() * C() * (1. - 2.*wtag()) /
                               (1. + pow(tau() * dm(), 2)));
        pdf /= norm_pdf;
        if (pdf < 0) {
            static int counter = 0;
            if (counter++ < 10)
                cerr << "ToyPdf::pdfSig < 0: dt " << dt
                     << ", c " << C() << ", s " << S() << endl;
            return 0.;
        }
    }
    return pdf;
}

double ToyPdf::pdfBkg(double dt, double wid) const {
    if (wid <= 0) return 0;
    return gaussian(dt, m_m, wid) / norm_gaussian(ll(), ul(), m_m, wid);
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
