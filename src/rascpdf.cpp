#include "rascpdf.h"

int RascPdf::ReadVars(const ICPVEvt& evt){
  m_vars.ReadVars(evt,RdetVar::AscSide);
  dz       = evt.GetDVar("dz_asc");
  keeptagl = evt.GetIVar("keeptagl");
  return 0;
}

double RascPdf::operator()(const ICPVEvt& evt){
  return PdfRascRnp(evt);
}

double RascPdf::operator()(const double& x){
  dz = x;
  return PdfRascRnp();
}

double RascPdf::RascRnp(const ICPVEvt &evt){
  ReadVars(evt); m_pars.Calculate(m_cnst,m_vars,keeptagl);
  return RascRnp();
}

double RascPdf::RascRnp(void) const{
  const double Li_md = gaussian(dz, m_pars.mu_main_asc(), m_pars.Smain_asc());
  const double Li_mp = Ep_conv_gauss(dz, m_pars.tau_np_p(), m_pars.mu_main_asc(), m_pars.Smain_asc());
  const double Li_mn = En_conv_gauss(dz, m_pars.tau_np_n(), m_pars.mu_main_asc(), m_pars.Smain_asc());
  const double Li_me = m_pars.fp()*Li_mp + (1.0 - m_pars.fp())*Li_mn;
  const double Li_mt = m_pars.fd()*Li_md + (1.0 - m_pars.fd())*Li_me;
  if(m_pars.ftail_asc()==0.0){ /* Mainly multiple track case */
    return Li_mt;
  }
  const double Li_td = gaussian(dz, m_pars.mu_tail_asc(), m_pars.Stail_asc());
  const double Li_tp = Ep_conv_gauss(dz,m_pars.tau_np_p_tl(), m_pars.mu_tail_asc(), m_pars.Stail_asc());
  const double Li_tn = En_conv_gauss(dz,m_pars.tau_np_n_tl(), m_pars.mu_tail_asc(), m_pars.Stail_asc());
  const double Li_te = m_pars.fp()*Li_tp + (1.0 - m_pars.fp())*Li_tn;
  const double Li_tt = m_pars.fd()*Li_td + (1.0 - m_pars.fd())*Li_te;
  double Li  = (1.0 - m_pars.ftail_asc())*Li_mt + m_pars.ftail_asc()*Li_tt;
  return Li;
}

double RascPdf::norm_RascRnp(const ICPVEvt& evt){
  ReadVars(evt); m_pars.Calculate(m_cnst,m_vars,keeptagl);
  return norm_RascRnp();
}

double RascPdf::norm_RascRnp(void) const{
  const double Li_md = norm_gaussian(m_ll, m_ul, m_pars.mu_main_asc(), m_pars.Smain_asc());
  const double Li_mp = norm_Ep_conv_gauss(m_ll, m_ul, m_pars.tau_np_p(), m_pars.mu_main_asc(), m_pars.Smain_asc());
  const double Li_mn = norm_En_conv_gauss(m_ll, m_ul, m_pars.tau_np_n(), m_pars.mu_main_asc(), m_pars.Smain_asc());
  const double Li_me = m_pars.fp()*Li_mp + (1.0 - m_pars.fp()) * Li_mn;
  const double Li_mt = m_pars.fd()*Li_md + (1.0 - m_pars.fd()) * Li_me;
  if(m_pars.ftail_asc()==0.0){ /* Mainly multiple track case */
    return Li_mt;
  }
  const double Li_td = norm_gaussian(m_ll, m_ul, m_pars.mu_tail_asc(), m_pars.Stail_asc());
  const double Li_tp = norm_Ep_conv_gauss(m_ll, m_ul, m_pars.tau_np_p_tl(), m_pars.mu_tail_asc(), m_pars.Stail_asc());
  const double Li_tn = norm_En_conv_gauss(m_ll, m_ul, m_pars.tau_np_n_tl(), m_pars.mu_tail_asc(), m_pars.Stail_asc());
  const double Li_te = m_pars.fp()*Li_tp + (1.0 - m_pars.fp()) * Li_tn;
  const double Li_tt = m_pars.fd()*Li_td + (1.0 - m_pars.fd()) * Li_te;
  double Li  = (1.0 - m_pars.ftail_asc()) * Li_mt + m_pars.ftail_asc()*Li_tt;
  return Li;
}

double RascPdf::PdfRascRnp(const ICPVEvt& evt){
  ReadVars(evt); m_pars.Calculate(m_cnst,m_vars,keeptagl);
  return PdfRascRnp();
}

double RascPdf::PdfRascRnp(void) const{
  const double      pdf =      RascRnp();
  const double norm_pdf = norm_RascRnp();
  return pdf/norm_pdf;
//  return AddOutlier(dz_asc,pdf,norm_pdf);
}

