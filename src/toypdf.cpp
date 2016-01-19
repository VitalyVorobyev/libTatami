#include "toypdf.h"

//#include <iostream>

using namespace std;

ToyPdf::ToyPdf(const ToyPdf& opdf)
  : ToyPdf()
{
  this->m_ll = opdf.m_ll;
  this->m_ul = opdf.m_ul;
}

double ToyPdf::operator() (const double& dt){
  const double pdf_s = pdfSig(dt,m_w);
  const double pdf_b = m_fbkg>0 ? pdfBkg(dt,m_w) : 0;
//  cout << "dt " << dt << ", pdf_s " << pdf_s << ", pdf_b " << pdf_b << ", fb " << m_fbkg << endl;
  return (1.-m_fbkg)*pdf_s + m_fbkg*pdf_b;
}

double ToyPdf::operator() (const double& dt, const double& fbkg, const double& scale){
  const double pdf_s = pdfSig(dt,m_w*scale);
  const double pdf_b = pdfBkg(dt,m_w*scale);
  return (1.-fbkg)*pdf_s + fbkg*pdf_b;
}

double ToyPdf::operator() (const double& dt, const double& c, const double& s, const double& fbkg, const double& scale){
  m_c = c; m_s = s;
  return this->operator()(dt,fbkg,scale);
}

double ToyPdf::pdfSig(const double& dt, const double& wid){
  const double pdf = Ef_conv_gauss(dt,m_tau,m_m,wid) + 0.5*(1.-2.*m_wrtag)/m_tau*(m_c*Mf_conv_gauss(dt,m_tau,m_dm,m_m,wid) + m_s*Af_conv_gauss(dt,m_tau,m_dm,m_m,wid));
  const double norm_pdf = norm_Ef_conv_gauss(m_ll,m_ul,m_tau,m_m,wid) + 0.5*(1.-2.*m_wrtag)/m_tau*m_c*norm_Mf_conv_gauss(m_ll,m_ul,m_tau,m_dm,m_m,wid);
//  cout << "tau " << m_tau << ", m " << m_m << ", w " << wid << ", pdf " << pdf << ", norm_pdf " << norm_pdf << endl;
  return pdf/norm_pdf;
}

double ToyPdf::pdfBkg(const double& dt, const double& wid){
  return gaussian(dt,m_m,wid)/norm_gaussian(m_ll,m_ul,m_m,wid);
}
