#include "rrecpdf.h"
#include "ttools.h"

double RrecPdf::Rrec(const ICPVEvt& evt){
  m_vars.ReadVars(evt,RdetVar::RecSide);
  m_pars.Calculate(m_cnst,m_vars);
  return Rrec();
}

double RrecPdf::Rrec(void) const{
  const double Li_mn = gaussian(dz,m_pars.mu_main_rec(),m_pars.Smain_rec());
  if(m_pars.ftail_rec()>0.0){ /* Mainly single track case */
    const double Li_tl = gaussian(dz,m_pars.mu_tail_rec(),m_pars.Stail_rec());
    const double Li    = (1.0-m_pars.ftail_rec())*Li_mn + m_pars.ftail_rec()*Li_tl;
    return Li;
  }
  return Li_mn;
}

double RrecPdf::norm_Rrec(const ICPVEvt &evt){
  m_vars.ReadVars(evt,RdetVar::RecSide);
  m_pars.Calculate(m_cnst,m_vars);
  return norm_Rrec();
}

double RrecPdf::norm_Rrec(void) const{
  const double Li_mn = norm_gaussian(m_ll,m_ul,m_pars.mu_main_rec(),m_pars.Smain_rec());
  if(m_pars.ftail_rec()>0.0){ /* Mainly multiple track case */
    const double Li_tl = norm_gaussian(m_ll,m_ul,m_pars.mu_tail_rec(),m_pars.Stail_rec());
    const double Li    = (1.0-m_pars.ftail_rec())*Li_mn + m_pars.ftail_rec()*Li_tl;
    return Li;
  }
  return Li_mn;
}

double RrecPdf::operator()(const double& x){
  dz = x;
  return Pdf();
}

double RrecPdf::Pdf(const ICPVEvt& evt){
  m_vars.ReadVars(evt,RdetVar::RecSide);
  m_pars.Calculate(m_cnst,m_vars);
  dz = evt.FindDVar("dz");
  return Pdf();
}

double RrecPdf::Pdf(void) const{
  const double      pdf =      Rrec();
  const double norm_pdf = norm_Rrec();
  return pdf/norm_pdf;
//  return AddOutlier(dz,pdf,norm_pdf);
}

double RrecPdf::operator()(const ICPVEvt& evt){
  return Pdf(evt);
}

//double RrecPdf::AddOutlier(const double& x, const double& pdf, const double& norm){
//  x += 1;
//  return pdf/norm;
//}

