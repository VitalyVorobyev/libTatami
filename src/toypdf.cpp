#include "toypdf.h"

ToyPdf::ToyPdf(const ToyPdf& opdf)
  : ToyPdf()
{
  this->m_ll = opdf.m_ll;
  this->m_ul = opdf.m_ul;
}

double ToyPdf::operator() (const double& dt){
  return (1.-m_fbkg)*pdfSig(dt,m_w) + m_fbkg*pdfBkg(dt,m_w);
}

double ToyPdf::operator() (const double& dt, const double& fbkg, const double& scale){
  return (1.-fbkg)*pdfSig(dt,m_w*scale) + fbkg*pdfBkg(dt,m_w*scale);
}

double ToyPdf::operator() (const double& dt, const double& c, const double& s, const double& fbkg, const double& scale){
  m_c = c; m_s = s;
  return this->operator()(dt,fbkg,scale);
}

double ToyPdf::pdfSig(const double& dt, const double& wid){
  const double pdf = Ef_conv_gauss(dt,m_tau,m_m,wid) + 0.5*(1.-2.*m_wrtag)/m_tau*(m_c*Mf_conv_gauss(dt,m_tau,m_dm,m_m,wid) + m_s*Af_conv_gauss(dt,m_tau,m_dm,m_m,wid));
  const double norm_pdf = norm_Ef_conv_gauss(m_ll,m_ul,m_tau,m_m,wid) + 0.5*(1.-2.*m_wrtag)/m_tau*m_c*norm_Mf_conv_gauss(m_ll,m_ul,m_tau,m_dm,m_m,wid);
  return pdf/norm_pdf;
}

double ToyPdf::pdfBkg(const double& dt, const double& wid){
  return gaussian(dt,m_m,wid)/norm_gaussian(m_ll,m_ul,m_m,wid);
}
