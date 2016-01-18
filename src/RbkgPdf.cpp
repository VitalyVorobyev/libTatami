#include "RbkgPdf.h"
#include <fstream>

using namespace std;

const double RbkgPdf::cm2ps = 78.48566945838871754705;

BkgPDFParSet::BkgPDFParSet():
 m_S_main_mlt(1.), m_S_tail_mlt(3.), m_f_tail_mlt(0.2), m_f_delta_mlt(0.7),
 m_S_main_sgl(1.), m_S_tail_sgl(3.), m_f_tail_sgl(0.2), m_f_delta_sgl(0.6),
 m_mu(0.), m_mu_delta(0.), m_tau(1.3), m_f_otlr(0.), m_s_otlr(30.)
{
}

BkgPDFParSet::BkgPDFParSet(const BkgPDFParSet& x){
  *this = x;
}

BkgPDFParSet& BkgPDFParSet::operator=(const BkgPDFParSet& x){
  m_S_main_mlt = x.m_S_main_mlt;
  m_S_tail_mlt = x.m_S_tail_mlt;
  m_f_tail_mlt = x.m_f_tail_mlt;
  m_f_delta_mlt= x.m_f_delta_mlt;
  m_S_main_sgl = x.m_S_main_sgl;
  m_S_tail_sgl = x.m_S_tail_sgl;
  m_f_tail_sgl = x.m_f_tail_sgl;
  m_f_delta_sgl= x.m_f_delta_sgl;
  m_mu         = x.m_mu;
  m_mu_delta   = x.m_mu_delta;
  m_tau        = x.m_tau;
  m_f_otlr     = x.m_f_otlr;
  m_s_otlr     = x.m_s_otlr;
  return *this;
}

RbkgPdf::RbkgPdf(const string& fname):
  AbsPdf(), m_scale(1), m_shift(0),
  m_sigma(1.), m_ndf(1)
{
  m_pars.GetParametersFromFile(fname);
}

RbkgPdf::RbkgPdf(void):
  RbkgPdf("../params/def_bkg.txt")
{
}

int BkgPDFParSet::GetParametersFromFile(const string& fname){
  ifstream ifile(fname.c_str(),ofstream::in);
  if(!ifile.is_open()){
    cout << "Can't open file " << fname << endl;
    return -1;
  } else{
    cout << "Getting background description from file " << fname << endl;
  }
  string line,name;
  double val,err;
  char namech[15];
  int counter = 0;
  for(int i=0; i<13; i++){
    getline(ifile,line);
    sscanf(line.c_str(),"%s = %lf +- %lf",namech,&val,&err);
    name = string(namech);
    cout << name << " " << val << " +- " << err << endl;
    if(name == string("tau"))        { m_tau         = val; counter++; continue;}
    if(name == string("mu"))         { m_mu          = val; counter++; continue;}
    if(name == string("mu_delta"))   { m_mu_delta    = val; counter++; continue;}
    if(name == string("f_delta_mlt")){ m_f_delta_mlt = val; counter++; continue;}
    if(name == string("f_tail_mlt")) { m_f_tail_mlt  = val; counter++; continue;}
    if(name == string("S_main_mlt")) { m_S_main_mlt  = val; counter++; continue;}
    if(name == string("S_tail_mlt")) { m_S_tail_mlt  = val; counter++; continue;}
    if(name == string("f_delta_sgl")){ m_f_delta_sgl = val; counter++; continue;}
    if(name == string("f_tail_sgl")) { m_f_tail_sgl  = val; counter++; continue;}
    if(name == string("S_main_sgl")) { m_S_main_sgl  = val; counter++; continue;}
    if(name == string("S_tail_sgl")) { m_S_tail_sgl  = val; counter++; continue;}
    if(name == string("f_otlr"))     { m_f_otlr      = val; counter++; continue;}
    if(name == string("s_otlr"))     { m_s_otlr      = val; counter++; continue;}
  }
  return counter;
}

double RbkgPdf::operator ()(const double& x){
  return Pdf(x,m_sigma,m_ndf);
}

double RbkgPdf::operator ()(const ICPVEvt& evt){
  const double dt = evt.FindDVar("dt");
  const double s  = evt.FindDVar("s");
  const int ndf   = evt.FindIVar("ndf");
  return Pdf(dt,s,ndf);
}

double RbkgPdf::Pdf(const double &x, const double &s, const int ndf){
  double sigma_main    = ndf ? s*m_pars.S_main_mlt()*cm2ps    : s*m_pars.S_main_sgl()*cm2ps;
  double sigma_tail    = ndf ? sigma_main*m_pars.S_tail_mlt() : sigma_main*m_pars.S_tail_sgl();
  const double f_tail  = ndf ? m_pars.f_tail_mlt()            : m_pars.f_tail_sgl();
  const double f_delta = ndf ? m_pars.f_delta_mlt()           : m_pars.f_delta_sgl();
  const double tau      = m_pars.tau();
  const double mu       = m_pars.mu()       + m_shift;
  const double mu_delta = m_pars.mu_delta() + m_shift;
  sigma_main *= m_scale;
  sigma_tail *= m_scale;

  const double pdf_l = (1-f_tail)*Enp_conv_gauss(x,tau,tau,mu,sigma_main)+f_tail*Enp_conv_gauss(x,tau,tau,mu,sigma_tail);
  const double int_pdf_l = (1-f_tail)*norm_Enp_conv_gauss(m_ll,m_ul,tau,tau,mu,sigma_main)+f_tail*norm_Enp_conv_gauss(m_ll,m_ul,tau,tau,mu,sigma_tail);

  const double pdf_d = (1-f_tail)*gaussian(x,mu_delta,sigma_main)+f_tail*gaussian(x,mu_delta,sigma_tail);
  const double int_pdf_d = (1-f_tail)*norm_gaussian(m_ll,m_ul,mu_delta,sigma_main)+f_tail*norm_gaussian(m_ll,m_ul,mu_delta,sigma_tail);

  const double pdf = f_delta*pdf_d + (1-f_delta)*pdf_l;
  const double int_pdf = f_delta*int_pdf_d + (1-f_delta)*int_pdf_l;
  if(pdf>=0 && int_pdf>=0){
    if(m_pars.f_otlr()>0.0001){ return AddOutlier(x,pdf,int_pdf);}
    else{                     return pdf/int_pdf;}
  } else{
    cout << "RbkgPdf::Pdf: " << pdf << ", norm = " << int_pdf << endl;
    return 0;
  }
}

double RbkgPdf::AddOutlier(const double& x, const double Lin,const double& nLi){
  const double Lol  = gaussian(x,0.,m_pars.s_otlr());
  const double nLol = norm_gaussian(m_ll,m_ul,0.,m_pars.s_otlr());
  const double Li   = (1.0-m_pars.f_otlr())*Lin/nLi + m_pars.f_otlr()*Lol/nLol;
  return Li;
}
