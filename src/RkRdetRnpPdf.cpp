#include "RkRdetRnpPdf.h"

using namespace std;

#ifndef DTRES_EXTERAM_THRE
#define DTRES_EXTERAM_THRE -FLT_MAX
#endif

//void RkRdetRnpPdf::SetFlv(const int x){
//  flv = x;
//  amix = flv;
//  cexp = 1;
//}

//void RkRdetRnpPdf::SetTag(const double& x){
//  if(fabs(x)>1.0001){
//    cout << x << " --- tag>1 !!!" << endl;
//    return;
//  }
//  flv  = x>0 ? 1 : -1;
//  amix = flv*m_wtag.Delut(x);
//  cexp = 1;//-flv*dw[i];
//  return;
//}

RkRdetRnpPdf::RkRdetRnpPdf(const DataClass& dc):
  AbsICPVPdf(), 
  m_wtag(ParManager::WTagFile(dc)),
  m_cnst(ParManager::SigParFile(dc))
{
}

int RkRdetRnpPdf::ReadVars(const ICPVEvt &evt){
  dt        = evt.FindDVar("dt");
  m_rvar.ReadVars(evt,RdetVar::RecSide);
  m_rvar.ReadVars(evt,RdetVar::RecSide);
  keeptagl  = evt.FindIVar("keeptagl");
  costhBcms = evt.FindDVar("costhBcms");
  return 0;
}

int RkRdetRnpPdf::ReadAndCalc(const ICPVEvt& evt){
  ReadVars(evt);
  m_rpar.Calculate(m_cnst,m_rvar);
  m_apar.Calculate(m_cnst,m_avar,keeptagl);
  m_apar.swap_rnp_param();
  m_kpar.SetAkCk(costhBcms,0.5*10.58,m_tau,m_dm);
  return 0;
}

double RkRdetRnpPdf::EfRkRdetRnp_fullrec(const ICPVEvt& evt){
  ReadAndCalc(evt);
  return EfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::EfRkRdetRnp_fullrec(void) const{
  const double Li_mm = EfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Smain_asc()));
  if(m_apar.ftail_asc()>0.0){
    const double Li_mt = EfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Stail_asc()));
    if(m_rpar.ftail_rec()>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = EfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
      const double Li_tt = EfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Stail_asc()));
      const double Li = ((1.0-m_rpar.ftail_rec())*(1.0-m_apar.ftail_asc())*Li_mm + m_rpar.ftail_rec()*(1.0-m_apar.ftail_asc())*Li_tm
                        +(1.0-m_rpar.ftail_rec())*m_apar.ftail_asc()*Li_mt + m_rpar.ftail_rec()*m_apar.ftail_asc()*Li_tt);
      return Li;
    } else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-m_apar.ftail_asc())*Li_mm + m_apar.ftail_asc()*Li_mt;
      return Li;
    }
  } else if(m_rpar.ftail_rec()>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = EfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
    const double Li = (1.0-m_rpar.ftail_rec())*Li_mm + m_rpar.ftail_rec()*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::EfRkRdetRnp_full_sup(const double& mu, const double& sigma) const{
  const double fd       = m_apar.fd();
  const double fp       = m_apar.fp();
  const double ntau_n   = m_kpar.ntau_n();
  const double ntau_p   = m_kpar.ntau_p();
  const double tau_np_n = m_apar.tau_np_n();
  const double tau_np_p = m_apar.tau_np_p();
  const double fact_n_e = m_kpar.fact_n_e();
  const double fact_p_e = m_kpar.fact_p_e();

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  double fEn = fact_n_e*fd;
  double fEp = fact_p_e*fd;
  double fEn_np = 0.0, fEp_np = 0.0, fxEn = 0.0, fxEp = 0.0;
  coco.add_EnEn_coef(fEn, fEn_np, fxEn,  ntau_n, tau_np_n, fact_n_e*nfn);
  coco.add_EnEp_coef(fEn, fEp_np,        ntau_n, tau_np_p, fact_n_e*nfp);
  coco.add_EpEn_coef(fEp, fEn_np,        ntau_p, tau_np_n, fact_p_e*nfn);
  coco.add_EpEp_coef(fEp, fEp_np, fxEp,  ntau_p, tau_np_p, fact_p_e*nfp);
  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  En_conv_gauss(dt, ntau_n,   mu, sigma);
  if(fEp!=0.0)    Li += fEp    *  Ep_conv_gauss(dt, ntau_p,   mu, sigma);
  if(fEn_np!=0.0) Li += fEn_np *  En_conv_gauss(dt, tau_np_n, mu, sigma);
  if(fEp_np!=0.0) Li += fEp_np *  Ep_conv_gauss(dt, tau_np_p, mu, sigma);
  if(fxEn!=0.0)   Li += fxEn   * xEn_conv_gauss(dt, ntau_n,   mu, sigma);
  if(fxEp!=0.0)   Li += fxEp   * xEp_conv_gauss(dt, ntau_p,   mu, sigma);
  return Li;
}

double RkRdetRnpPdf::AfRkRdetRnp_fullrec(const ICPVEvt& evt){
  ReadAndCalc(evt);
  return AfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::AfRkRdetRnp_fullrec(void) const{
  const double Li_mm = AfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Smain_asc()));
  if(m_apar.ftail_asc()>0.0){
    const double Li_mt = AfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Stail_asc()));
    if(m_rpar.ftail_rec()>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = AfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
      const double Li_tt = AfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Stail_asc()));
      const double Li = ((1.0-m_rpar.ftail_rec())*(1.0-m_apar.ftail_asc())*Li_mm + m_rpar.ftail_rec()*(1.0-m_apar.ftail_asc())*Li_tm
                        +(1.0-m_rpar.ftail_rec())*m_apar.ftail_asc()*Li_mt + m_rpar.ftail_rec()*m_apar.ftail_asc()*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-m_apar.ftail_asc())*Li_mm + m_apar.ftail_asc()*Li_mt;
      return Li;
    }
  }else if(m_rpar.ftail_rec()>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = AfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
    const double Li = (1.0-m_rpar.ftail_rec())*Li_mm + m_rpar.ftail_rec()*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::AfRkRdetRnp_full_sup(const double& mu, const double& sigma) const{
  const double fd       = m_apar.fd();
  const double fp       = m_apar.fp();
  const double ntau_n   = m_kpar.ntau_n();
  const double ntau_p   = m_kpar.ntau_p();
  const double tau_np_n = m_apar.tau_np_n();
  const double tau_np_p = m_apar.tau_np_p();
  const double fact_am  = m_kpar.fact_am();
  const double ndmtau   = m_kpar.ndmtau();
  const double ndm_n    = m_kpar.ndm_n();
  const double ndm_p    = m_kpar.ndm_p();

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_mn_n = -fact_am*ndmtau;
  const double w_mn_p = -fact_am*ndmtau;
  double fAn = fact_am*fd;
  double fAp = fact_am*fd;
  double fMn = w_mn_n*fd;
  double fMp = w_mn_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0;
  coco.add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_am*nfn);
  coco.add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_am*nfp);
  coco.add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, w_mn_n*nfn);
  coco.add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, w_mn_n*nfp);
  coco.add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_am*nfn);
  coco.add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_am*nfp);
  coco.add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, w_mn_p*nfn);
  coco.add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, w_mn_p*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * An_conv_gauss(dt, ntau_n, ndm_n, mu, sigma);
  if(fAp!=0.0)    Li += fAp    * Ap_conv_gauss(dt, ntau_p, ndm_p, mu, sigma);
  if(fMn!=0.0)    Li += fMn    * Mn_conv_gauss(dt, ntau_n, ndm_n, mu, sigma);
  if(fMp!=0.0)    Li += fMp    * Mp_conv_gauss(dt, ntau_p, ndm_p, mu, sigma);
  if(fEn_np!=0.0) Li += fEn_np * En_conv_gauss(dt, tau_np_n,      mu, sigma);
  if(fEp_np!=0.0) Li += fEp_np * Ep_conv_gauss(dt, tau_np_p,      mu, sigma);
  return Li;
}

double RkRdetRnpPdf::MfRkRdetRnp_fullrec(const ICPVEvt& evt){
  ReadAndCalc(evt);
  return MfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::MfRkRdetRnp_fullrec(void) const{
  const double Li_mm = MfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Smain_asc()));
  if(m_apar.ftail_asc()>0.0){
    const double Li_mt = MfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Stail_asc()));
    if(m_rpar.ftail_rec()>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = MfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
      const double Li_tt = MfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Stail_asc()));
      const double Li = ((1.0-m_rpar.ftail_rec())*(1.0-m_apar.ftail_asc())*Li_mm + m_rpar.ftail_rec()*(1.0-m_apar.ftail_asc())*Li_tm
                        +(1.0-m_rpar.ftail_rec())*m_apar.ftail_asc()*Li_mt + m_rpar.ftail_rec()*m_apar.ftail_asc()*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-m_apar.ftail_asc())*Li_mm + m_apar.ftail_asc()*Li_mt;
      return Li;
    }
  }else if(m_rpar.ftail_rec()>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = MfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
    const double Li = (1.0-m_rpar.ftail_rec())*Li_mm + m_rpar.ftail_rec()*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::MfRkRdetRnp_full_sup(const double& mu, const double& sigma) const{
  const double fd       = m_apar.fd();
  const double fp       = m_apar.fp();
  const double ntau_n   = m_kpar.ntau_n();
  const double ntau_p   = m_kpar.ntau_p();
  const double tau_np_n = m_apar.tau_np_n();
  const double tau_np_p = m_apar.tau_np_p();
  const double fact_am  = m_kpar.fact_am();
  const double ndmtau   = m_kpar.ndmtau();
  const double ndm_n    = m_kpar.ndm_n();
  const double ndm_p    = m_kpar.ndm_p();
  const double cktau    = m_kpar.cktau();

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_an_n = fact_am*ndmtau;
  const double w_an_p = fact_am*ndmtau;
  double fMn = fact_am*fd;
  double fMp = fact_am*fd;
  double fAn = w_an_n*fd;
  double fAp = w_an_p*fd;
  double fEn_k = 0.0, fEp_k = 0.0, fxEn_k = 0.0, fxEp_k = 0.0;
  double fEn_np = 0.0, fEp_np = 0.0;
  coco.add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, w_an_n*nfn);
  coco.add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, w_an_n*nfp);
  coco.add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_am*nfn);
  coco.add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_am*nfp);
  coco.add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, w_an_p*nfn);
  coco.add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, w_an_p*nfp);
  coco.add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_am*nfn);
  coco.add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_am*nfp);

  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    *  An_conv_gauss(dt, ntau_n, ndm_n, mu, sigma);
  if(fAp!=0.0)    Li += fAp    *  Ap_conv_gauss(dt, ntau_p, ndm_p, mu, sigma);
  if(fMn!=0.0)    Li += fMn    *  Mn_conv_gauss(dt, ntau_n, ndm_n, mu, sigma);
  if(fMp!=0.0)    Li += fMp    *  Mp_conv_gauss(dt, ntau_p, ndm_p, mu, sigma);
  if(fEn_k!=0.0)  Li += fEn_k  *  En_conv_gauss(dt,-cktau,         mu, sigma);
  if(fEp_k!=0.0)  Li += fEp_k  *  Ep_conv_gauss(dt, cktau,         mu, sigma);
  if(fxEn_k!=0.0) Li += fxEn_k * xEn_conv_gauss(dt,-cktau,         mu, sigma);
  if(fxEp_k!=0.0) Li += fxEp_k * xEp_conv_gauss(dt, cktau,         mu, sigma);
  if(fEn_np!=0.0) Li += fEn_np *  En_conv_gauss(dt, tau_np_n,      mu, sigma);
  if(fEp_np!=0.0) Li += fEp_np *  Ep_conv_gauss(dt, tau_np_p,      mu, sigma);
  return Li;
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_fullrec(const ICPVEvt& evt){
  ReadAndCalc(evt);
  return norm_EfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_fullrec(void) const{
  const double Li_mm = norm_EfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Smain_asc()));
  if(m_apar.ftail_asc()>0.0){
    const double Li_mt= norm_EfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Stail_asc()));
    if(m_rpar.ftail_rec()>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = norm_EfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
      const double Li_tt = norm_EfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Stail_asc()));
      const double Li = ((1.0-m_rpar.ftail_rec())*(1.0-m_apar.ftail_asc())*Li_mm + m_rpar.ftail_rec()*(1.0-m_apar.ftail_asc())*Li_tm
                        +(1.0-m_rpar.ftail_rec())*     m_apar.ftail_asc() *Li_mt + m_rpar.ftail_rec()*     m_apar.ftail_asc() *Li_tt);
      return Li;
    } else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-m_apar.ftail_asc())*Li_mm + m_apar.ftail_asc()*Li_mt;
      return Li;
    }
  } else if(m_rpar.ftail_rec()>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = norm_EfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
    const double Li = (1.0-m_rpar.ftail_rec())*Li_mm + m_rpar.ftail_rec()*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_full_sup(const double& mu, const double& sigma) const{
  const double fd       = m_apar.fd();
  const double fp       = m_apar.fp();
  const double ntau_n   = m_kpar.ntau_n();
  const double ntau_p   = m_kpar.ntau_p();
  const double tau_np_n = m_apar.tau_np_n();
  const double tau_np_p = m_apar.tau_np_p();
  const double fact_n_e = m_kpar.fact_n_e();
  const double fact_p_e = m_kpar.fact_p_e();

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  double fEn = fact_n_e*fd;
  double fEp = fact_p_e*fd;
  double fEn_np = 0.0, fEp_np = 0.0, fxEn = 0.0, fxEp = 0.0;
  coco.add_EnEn_coef(fEn,fEn_np,fxEn,ntau_n, tau_np_n, fact_n_e*nfn);
  coco.add_EnEp_coef(fEn,fEp_np,     ntau_n, tau_np_p, fact_n_e*nfp);
  coco.add_EpEn_coef(fEp,fEn_np,     ntau_p, tau_np_n, fact_p_e*nfn);
  coco.add_EpEp_coef(fEp,fEp_np,fxEp,ntau_p, tau_np_p, fact_p_e*nfp);
  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  norm_En_conv_gauss(m_ll, m_ul, ntau_n,   mu, sigma);
  if(fEp!=0.0)    Li += fEp    *  norm_Ep_conv_gauss(m_ll, m_ul, ntau_p,   mu, sigma);
  if(fEn_np!=0.0) Li += fEn_np *  norm_En_conv_gauss(m_ll, m_ul, tau_np_n, mu, sigma);
  if(fEp_np!=0.0) Li += fEp_np *  norm_Ep_conv_gauss(m_ll, m_ul, tau_np_p, mu, sigma);
  if(fxEn!=0.0)   Li += fxEn   * norm_xEn_conv_gauss(m_ll, m_ul, ntau_n,   mu, sigma);
  if(fxEp!=0.0)   Li += fxEp   * norm_xEp_conv_gauss(m_ll, m_ul, ntau_p,   mu, sigma);
  return Li;
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_fullrec(const ICPVEvt& evt){
  ReadAndCalc(evt);
  return norm_MfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_fullrec(void) const{
  const double Li_mm = norm_MfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Smain_asc()));
  if(m_apar.ftail_asc()>0.0){
    const double Li_mt = norm_MfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Stail_asc()));
    if(m_rpar.ftail_rec()>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = norm_MfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
      const double Li_tt = norm_MfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Stail_asc()));
      const double Li = ((1.0-m_rpar.ftail_rec())*(1.0-m_apar.ftail_asc())*Li_mm + m_rpar.ftail_rec()*(1.0-m_apar.ftail_asc())*Li_tm
                        +(1.0-m_rpar.ftail_rec())*m_apar.ftail_asc()*Li_mt + m_rpar.ftail_rec()*m_apar.ftail_asc()*Li_tt);
      return Li;
    } else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-m_apar.ftail_asc())*Li_mm + m_apar.ftail_asc()*Li_mt;
      return Li;
    }
  } else if(m_rpar.ftail_rec()>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = norm_MfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
    const double Li = (1.0-m_rpar.ftail_rec())*Li_mm + m_rpar.ftail_rec()*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_full_sup(const double& mu, const double& sigma) const{
  const double fd       = m_apar.fd();
  const double fp       = m_apar.fp();
  const double ntau_n   = m_kpar.ntau_n();
  const double ntau_p   = m_kpar.ntau_p();
  const double tau_np_n = m_apar.tau_np_n();
  const double tau_np_p = m_apar.tau_np_p();
  const double fact_am  = m_kpar.fact_am();
  const double ndmtau   = m_kpar.ndmtau();
  const double ndm_n    = m_kpar.ndm_n();
  const double ndm_p    = m_kpar.ndm_p();
  const double cktau    = m_kpar.cktau();

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_an_n = fact_am*ndmtau;
  const double w_an_p = fact_am*ndmtau;
  double fMn = fact_am*fd;
  double fMp = fact_am*fd;
  double fAn = w_an_n*fd;
  double fAp = w_an_p*fd;
  double fEn_k = 0.0, fEp_k = 0.0, fxEn_k = 0.0, fxEp_k = 0.0;
  double fEn_np = 0.0, fEp_np = 0.0;

  coco.add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, w_an_n*nfn);
  coco.add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, w_an_n*nfp);
  coco.add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_am*nfn);
  coco.add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_am*nfp);
  coco.add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, w_an_p*nfn);
  coco.add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, w_an_p*nfp);
  coco.add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_am*nfn);
  coco.add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_am*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    *  norm_An_conv_gauss(m_ll, m_ul, ntau_n, ndm_n, mu, sigma);
  if(fAp!=0.0)    Li += fAp    *  norm_Ap_conv_gauss(m_ll, m_ul, ntau_p, ndm_p, mu, sigma);
  if(fMn!=0.0)    Li += fMn    *  norm_Mn_conv_gauss(m_ll, m_ul, ntau_n, ndm_n, mu, sigma);
  if(fMp!=0.0)    Li += fMp    *  norm_Mp_conv_gauss(m_ll, m_ul, ntau_p, ndm_p, mu, sigma);
  if(fEn_k!=0.0)  Li += fEn_k  *  norm_En_conv_gauss(m_ll, m_ul,-cktau,         mu, sigma);
  if(fEp_k!=0.0)  Li += fEp_k  *  norm_Ep_conv_gauss(m_ll, m_ul, cktau,         mu, sigma);
  if(fxEn_k!=0.0) Li += fxEn_k * norm_xEn_conv_gauss(m_ll, m_ul,-cktau,         mu, sigma);
  if(fxEp_k!=0.0) Li += fxEp_k * norm_xEp_conv_gauss(m_ll, m_ul, cktau,         mu, sigma);
  if(fEn_np!=0.0) Li += fEn_np *  norm_En_conv_gauss(m_ll, m_ul, tau_np_n,      mu, sigma);
  if(fEp_np!=0.0) Li += fEp_np *  norm_Ep_conv_gauss(m_ll, m_ul, tau_np_p,      mu, sigma);
  return Li;
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_fullrec(const ICPVEvt& evt){
  ReadAndCalc(evt); 
  return norm_AfRkRdetRnp_fullrec();
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_fullrec(void) const{
  const double Li_mm = norm_AfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Smain_asc()));
  if(m_apar.ftail_asc()>0.0){
    const double Li_mt = norm_AfRkRdetRnp_full_sup(m_rpar.mu_main_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Smain_rec(),m_apar.Stail_asc()));
    if(m_rpar.ftail_rec()>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = norm_AfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
      const double Li_tt = norm_AfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_tail_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Stail_asc()));
      const double Li = ((1.0-m_rpar.ftail_rec())*(1.0-m_apar.ftail_asc())*Li_mm + m_rpar.ftail_rec()*(1.0-m_apar.ftail_asc())*Li_tm
                        +(1.0-m_rpar.ftail_rec())*m_apar.ftail_asc()*Li_mt + m_rpar.ftail_rec()*m_apar.ftail_asc()*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      return (1.0-m_apar.ftail_asc())*Li_mm + m_apar.ftail_asc()*Li_mt;
    }
  }else if(m_rpar.ftail_rec()>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = norm_AfRkRdetRnp_full_sup(m_rpar.mu_tail_rec()+m_apar.mu_main_asc(),TTools::sum_sigma(m_rpar.Stail_rec(),m_apar.Smain_asc()));
    const double Li = (1.0-m_rpar.ftail_rec())*Li_mm + m_rpar.ftail_rec()*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_full_sup(const double& mu, const double& sigma) const{
  const double fd       = m_apar.fd();
  const double fp       = m_apar.fp();
  const double ntau_n   = m_kpar.ntau_n();
  const double ntau_p   = m_kpar.ntau_p();
  const double tau_np_n = m_apar.tau_np_n();
  const double tau_np_p = m_apar.tau_np_p();
  const double fact_am = m_kpar.fact_am();
  const double ndmtau   = m_kpar.ndmtau();
  const double ndm_n    = m_kpar.ndm_n();
  const double ndm_p    = m_kpar.ndm_p();

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_mn_n = -fact_am*ndmtau;
  const double w_mn_p = -fact_am*ndmtau;
  double fAn = fact_am*fd;
  double fAp = fact_am*fd;
  double fMn = w_mn_n*fd;
  double fMp = w_mn_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0;
  coco.add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_am*nfn);
  coco.add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_am*nfp);
  coco.add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, w_mn_n*nfn);
  coco.add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, w_mn_n*nfp);
  coco.add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_am*nfn);
  coco.add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_am*nfp);
  coco.add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, w_mn_p*nfn);
  coco.add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, w_mn_p*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * norm_An_conv_gauss(m_ll, m_ul, ntau_n, ndm_n, mu, sigma);
  if(fAp!=0.0)    Li += fAp    * norm_Ap_conv_gauss(m_ll, m_ul, ntau_p, ndm_p, mu, sigma);
  if(fMn!=0.0)    Li += fMn    * norm_Mn_conv_gauss(m_ll, m_ul, ntau_n, ndm_n, mu, sigma);
  if(fMp!=0.0)    Li += fMp    * norm_Mp_conv_gauss(m_ll, m_ul, ntau_p, ndm_p, mu, sigma);
  if(fEn_np!=0.0) Li += fEn_np * norm_En_conv_gauss(m_ll, m_ul, tau_np_n,      mu, sigma);
  if(fEp_np!=0.0) Li += fEp_np * norm_Ep_conv_gauss(m_ll, m_ul, tau_np_p,      mu, sigma);
  return Li;
}

double RkRdetRnpPdf::PdfAB(const ICPVEvt& evt, const bool otlr, const bool no_interf){
  ReadAndCalc(evt);
  return PdfAB(otlr,no_interf);
}

double RkRdetRnpPdf::PdfAB(const bool otlr, const bool no_interf) const{
  const double      Ef =      EfRkRdetRnp_fullrec();
  const double norm_Ef = norm_EfRkRdetRnp_fullrec();
  if(!no_interf){
    const double      Mf = 0.5/m_tau*     MfRkRdetRnp_fullrec();
    const double      Af = 0.5/m_tau*     AfRkRdetRnp_fullrec();
    const double norm_Mf = 0.5/m_tau*norm_MfRkRdetRnp_fullrec();
    const double pdf =           Ef*cexp + amix*(m_c*Mf - m_s*Af);
    const double pdf_norm = norm_Ef*cexp + amix*m_c*norm_Mf;
    if(pdf<=0 || pdf_norm<=0){
      cout << "PdfAB. pdf: " << pdf << ", norm: " << pdf_norm;
      cout << ", cexp: " << cexp << ", amix: " << amix;// << ", flv: " << flv;
      cout << endl;
      return -fabs(pdf/pdf_norm);
    }
    return pdf/pdf_norm;
//    if(otlr) return AddOutlier(dt,pdf,pdf_norm);
//    else     return pdf/pdf_norm;
  } else{
    return Ef/norm_Ef;
//    if(otlr) return AddOutlier(dt,Ef,norm_Ef);
//    else     return Ef/norm_Ef;
  }
  return -999;
}

//double RkRdetRnpPdf::Pdf(const ICPVEvt& evt, const bool otlr, const bool no_interf){
//  ReadVars(evt);
//  double m_pdf = 0l;
//  const double      Ef =      EfRkRdetRnp_fullrec();
//  const double norm_Ef = norm_EfRkRdetRnp_fullrec();
//  if(!no_interf){
//    const double      Mf = 0.5/m_tau*     MfRkRdetRnp_fullrec();
//    const double      Af = 0.5/m_tau*     AfRkRdetRnp_fullrec();
//    const double norm_Mf = 0.5/m_tau*norm_MfRkRdetRnp_fullrec();
//    const double pdf      = (K+Kb)*Ef*cexp      + amix*((K-Kb)*Mf - 2.*xi*sqrt(K*Kb)*Af*(S*cos2beta-C*sin2beta));
//    const double pdf_norm = (K+Kb)*norm_Ef*cexp + amix*(K-Kb)*norm_Mf;
//    if(pdf<=0. || pdf_norm<=0.){
//      cout << "Pdf1: pdf = " << pdf << ", norm = " << pdf_norm << ", cexp: " << cexp << ", amix: " << amix << endl;
//      return pdf/pdf_norm;
//    }
//    if(otlr) m_pdf = AddOutlier(dt,pdf,pdf_norm);
//    else     m_pdf = pdf/pdf_norm;
//    return m_pdf;
//  } else{
//    if(Ef<=0. || norm_Ef<=0.){
//      cout << "Pdf: pdf = " << Ef << ", norm = " << norm_Ef << endl;
//      return Ef/norm_Ef;
//    }
//    if(otlr) m_pdf = AddOutlier(dt,Ef,norm_Ef);
//    else     m_pdf = Ef/norm_Ef;
//    return m_pdf;
//  }
//  return -999;
//}

//double RkRdetRnpPdf::AddOutlier(const double& x, const double& Lin, const double& nLi, const double& alpha){
//  return Add_Outlier(x,Lin,nLi,alpha);
//}

//double RkRdetRnpPdf::Add_Outlier(const double& x, const double& Lin, const double& nLi, const double& alpha){
//  const double fol = ((m_svd == 2 || ntrk_rec>1) && ntrk_asc>1) ? f_ol_mul : f_ol_sgl;
//  const double m = 0.0;
//  const double Lol = gaussian(x,m,sigma_ol);
//  const double nLol = norm_gaussian(m_ll,m_ul,m,sigma_ol);
//  const double Li = (1.0-fol)*Lin/nLi + fol*alpha*Lol/nLol;
//  return Li;
//}

double RkRdetRnpPdf::operator()(const ICPVEvt& evt){
  ReadAndCalc(evt);
  return PdfAB();
}

double RkRdetRnpPdf::operator()(const double& x){
  dt = x;
  return PdfAB();
}

