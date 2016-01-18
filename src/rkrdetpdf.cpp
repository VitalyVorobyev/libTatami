#include "rkrdetpdf.h"

double RkRdetRnpPdf::EfRkRdet_fullrec(const ICPVEvt& evt){
  ReadVars(evt);
  return EfRkRdet_fullrec();
}

double RkRdetRnpPdf::EfRkRdet_fullrec(void){
//  Set_fact_Ef();
  calc_vtxparam_rec(); Rrec_param();
  calc_vtxparam_asc(); Rasc_param();

  const double Li_mm = EfRkRdet_full_sup(mu_main_rec+mu_main_asc,sum_sigma(Smain_rec,Smain_asc));
  if(ftail_asc>0.0){
    const double Li_mt = EfRkRdet_full_sup(mu_main_rec+mu_tail_asc,sum_sigma(Smain_rec,Stail_asc));
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = EfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
      const double Li_tt = EfRkRdet_full_sup(mu_tail_rec+mu_tail_asc,sum_sigma(Stail_rec,Stail_asc));
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                        +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = EfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::EfRkRdet_full_sup(const double& mu, const double& sigma){
  return fact_n_e*En_conv_gauss(dt,ntau_n,mu,sigma)
        +fact_p_e*Ep_conv_gauss(dt,ntau_p,mu,sigma);
}

double RkRdetRnpPdf::AfRkRdet_fullrec(const ICPVEvt& evt){
  ReadVars(evt);
  return AfRkRdet_fullrec();
}

double RkRdetRnpPdf::AfRkRdet_fullrec(void){
//  Set_fact_AfMf();
  calc_vtxparam_rec(); Rrec_param();
  calc_vtxparam_asc(); Rasc_param();

  const double Li_mm = AfRkRdet_full_sup(mu_main_rec+mu_main_asc,sum_sigma(Smain_rec,Smain_asc));
  if(ftail_asc>0.0){
    const double Li_mt = AfRkRdet_full_sup(mu_main_rec+mu_tail_asc,sum_sigma(Smain_rec,Stail_asc));
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = AfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
      const double Li_tt = AfRkRdet_full_sup(mu_tail_rec+mu_tail_asc,sum_sigma(Stail_rec,Stail_asc));
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                        +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = AfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::AfRkRdet_full_sup(const double& mu, const double& sigma){
  return fact_am*(An_conv_gauss(dt,ntau_n,ndm_n,mu,sigma)-ndmtau*Mn_conv_gauss(dt,ntau_n,ndm_n,mu,sigma))
       + fact_am*(Ap_conv_gauss(dt,ntau_p,ndm_p,mu,sigma)-ndmtau*Mp_conv_gauss(dt,ntau_p,ndm_p,mu,sigma));
}

double RkRdetRnpPdf::MfRkRdet_fullrec(const ICPVEvt& evt){
  ReadVars(evt);
  return MfRkRdet_fullrec();
}

double RkRdetRnpPdf::MfRkRdet_fullrec(void){
//  Set_fact_AfMf();
  calc_vtxparam_rec(); Rrec_param();
  calc_vtxparam_asc(); Rasc_param();

  const double Li_mm = MfRkRdet_full_sup(mu_main_rec+mu_main_asc,sum_sigma(Smain_rec,Smain_asc));
  if(ftail_asc>0.0){
    const double Li_mt = MfRkRdet_full_sup(mu_main_rec+mu_tail_asc,sum_sigma(Smain_rec,Stail_asc));
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = MfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
      const double Li_tt = MfRkRdet_full_sup(mu_tail_rec+mu_tail_asc,sum_sigma(Stail_rec,Stail_asc));
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                        +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = MfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::MfRkRdet_full_sup(const double& mu, const double& sigma){
  return fact_am*(Mn_conv_gauss(dt,ntau_n,ndm_n,mu,sigma)+ndmtau*An_conv_gauss(dt,ntau_n,ndm_n,mu,sigma))
       + fact_am*(Mp_conv_gauss(dt,ntau_p,ndm_p,mu,sigma)+ndmtau*Ap_conv_gauss(dt,ntau_p,ndm_p,mu,sigma));
}

double RkRdetRnpPdf::norm_EfRkRdet_fullrec(const ICPVEvt& evt){
  ReadVars(evt);
  return norm_EfRkRdet_fullrec();
}

double RkRdetRnpPdf::norm_EfRkRdet_fullrec(void){
//  Set_fact_Ef();
  calc_vtxparam_rec(); Rrec_param();
  calc_vtxparam_asc(); Rasc_param();

  const double Li_mm = norm_EfRkRdet_full_sup(mu_main_rec+mu_main_asc,sum_sigma(Smain_rec,Smain_asc));
  if(ftail_asc>0.0){
    const double Li_mt = norm_EfRkRdet_full_sup(mu_main_rec+mu_tail_asc,sum_sigma(Smain_rec,Stail_asc));
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = norm_EfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
      const double Li_tt = norm_EfRkRdet_full_sup(mu_tail_rec+mu_tail_asc,sum_sigma(Stail_rec,Stail_asc));
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                        +(1.0-ftail_rec)*ftail_asc*Li_mt +       ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  } else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = norm_EfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_EfRkRdet_full_sup(const double& mu, const double& sigma){
  return fact_n_e*norm_En_conv_gauss(m_ll,m_ul,ntau_n,mu,sigma)
       + fact_p_e*norm_Ep_conv_gauss(m_ll,m_ul,ntau_p,mu,sigma);
}

double RkRdetRnpPdf::norm_AfRkRdet_fullrec(const ICPVEvt& evt){
  ReadVars(evt);
  return norm_AfRkRdet_fullrec();
}

double RkRdetRnpPdf::norm_AfRkRdet_fullrec(void){
//  Set_fact_AfMf();
  calc_vtxparam_rec(); Rrec_param();
  calc_vtxparam_asc(); Rasc_param();

  const double Li_mm = norm_AfRkRdet_full_sup(mu_main_rec+mu_main_asc,sum_sigma(Smain_rec,Smain_asc));
  if(ftail_asc>0.0){
    const double Li_mt    = norm_AfRkRdet_full_sup(mu_main_rec+mu_tail_asc,sum_sigma(Smain_rec,Stail_asc));
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = norm_AfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
      const double Li_tt = norm_AfRkRdet_full_sup(mu_tail_rec+mu_tail_asc,sum_sigma(Stail_rec,Stail_asc));
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                        +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = norm_AfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_AfRkRdet_full_sup(const double& mu, const double& sigma){
  return fact_am*(norm_An_conv_gauss(m_ll,m_ul,ntau_n,ndm_n,mu,sigma)-ndmtau*norm_Mn_conv_gauss(m_ll,m_ul,ntau_n,ndm_n,mu,sigma))
       + fact_am*(norm_Ap_conv_gauss(m_ll,m_ul,ntau_p,ndm_p,mu,sigma)-ndmtau*norm_Mp_conv_gauss(m_ll,m_ul,ntau_p,ndm_p,mu,sigma));
}

double RkRdetRnpPdf::norm_MfRkRdet_fullrec(const ICPVEvt& evt){
  ReadVars(evt);
  return norm_MfRkRdet_fullrec();
}

double RkRdetRnpPdf::norm_MfRkRdet_fullrec(void){
  calc_vtxparam_rec(); Rrec_param();
  calc_vtxparam_asc(); Rasc_param();

  const double Li_mm = norm_MfRkRdet_full_sup(mu_main_rec+mu_main_asc,sum_sigma(Smain_rec,Smain_asc));
  if(ftail_asc>0.0){
    const double Li_mt = norm_MfRkRdet_full_sup(mu_main_rec+mu_tail_asc,sum_sigma(Smain_rec,Stail_asc));
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double Li_tm = norm_MfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec,Smain_asc));
      const double Li_tt = norm_MfRkRdet_full_sup(mu_tail_rec+mu_tail_asc,sum_sigma(Stail_rec,Stail_asc));
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                        +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double Li_tm = norm_MfRkRdet_full_sup(mu_tail_rec+mu_main_asc,sum_sigma(Stail_rec, Smain_asc));
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_MfRkRdet_full_sup(const double& mu, const double& sigma){
  return fact_am*(norm_Mn_conv_gauss(m_ll,m_ul,ntau_n,ndm_n,mu,sigma)+ndmtau*norm_An_conv_gauss(m_ll,m_ul,ntau_n,ndm_n,mu,sigma))
       + fact_am*(norm_Mp_conv_gauss(m_ll,m_ul,ntau_p,ndm_p,mu,sigma)+ndmtau*norm_Ap_conv_gauss(m_ll,m_ul,ntau_p,ndm_p,mu,sigma));
}

double RkRdetRnpPdf::NoNPPdf(const ICPVEvt& evt, const bool otlr, const bool no_interf){
  ReadVars(evt);
  const double      Ef =      EfRkRdet_fullrec();
  const double norm_Ef = norm_EfRkRdet_fullrec();
  if(!no_interf){
    const double      Mf = 0.5/m_tau*MfRkRdet_fullrec();
    const double      Af = 0.5/m_tau*AfRkRdet_fullrec();
    const double norm_Mf = 0.5/m_tau*norm_MfRkRdet_fullrec();
    const double pdf      = (K+Kb)*Ef      + flv*(K-Kb)*Mf - 2.*flv*xi*sqrt(K*Kb)*Af*(S*cos2beta-C*sin2beta);
    const double pdf_norm = (K+Kb)*norm_Ef + flv*(K-Kb)*norm_Mf;
    if(pdf<0 || pdf_norm<=0){
      cout << "Pdf: pdf = " << pdf << ", norm = " << pdf_norm << endl;
      return 0;
    }
    if(otlr) return AddOutlier(dt,pdf,pdf_norm);
    else     return pdf/pdf_norm;
  } else{
    if(otlr) return AddOutlier(dt,Ef,norm_Ef);
    else     return Ef/norm_Ef;
  }
  return -999;
}

