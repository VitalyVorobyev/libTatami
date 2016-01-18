#include "rascrnppars.h"
#include "ttools.h"

#include <cfloat>

#ifndef DTRES_EXTERAM_THRE
#define DTRES_EXTERAM_THRE -FLT_MAX
//#define DTRES_EXTERAM_THRE -999
#endif /* DTRES_EXTERAM_THRE */

RascRnpPars::RascRnpPars(void):
  m_Smain_asc(1.), m_Stail_asc(1.), m_ftail_asc(0.),
  m_mu_main_asc(0.), m_mu_tail_asc(0.), m_xi_asc(0.), m_st_asc(1.),
  m_fd(0.), m_fp(0.), m_tau_np_p(1.), m_tau_np_n(1.), m_tau_np_p_tl(1.), m_tau_np_n_tl(1.)
{
}

void RascRnpPars::calc_vtxparam_asc(const RdetVar& var){
  m_xi_asc = (var.ndf() > 0 ? var.chisq()/var.ndf() : 1);
  m_st_asc = var.sz()*TTools::cm2ps;
  return;
}

void RascRnpPars::Rasc_param(const ResConst& cnst,const RdetVar& var){
  if(var.ntrk()>1){ /*  multiple track vertex  */
    m_ftail_asc = 0;//ftl_asc_mlt[0] + ftl_asc_mlt[1]*xi_asc;
    m_Smain_asc = (cnst.Sasc(0) + cnst.Sasc(1)*m_xi_asc)*m_st_asc;
    m_Stail_asc = m_Smain_asc;
    TTools::constraint(m_ftail_asc, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    TTools::constraint(m_Smain_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }else{ /*  single track vertex */
    m_ftail_asc = cnst.ftl_asc();
    m_Smain_asc = cnst.Smn_asc()*m_st_asc;
    m_Stail_asc = cnst.Stl_asc()*m_st_asc;
    TTools::constraint(m_ftail_asc,0.0,1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    TTools::constraint(m_Smain_asc, __POSTIVEPARAM_MIN__);
    TTools::constraint(m_Stail_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
//  m_mu_main_asc = 0.0;
//  m_mu_tail_asc = 0.0;
  return;
}

void RascRnpPars::Rnp_paramOld(const ResConst& cnst,const RdetVar& var, const int keeptagl){
  if(var.ntrk()>1){ /*  multiple track vertex  */
    m_fd = cnst.fd_np_mlt(keeptagl);
    m_fp = cnst.fp_np_mlt();
    m_tau_np_p = cnst.tau_np_p_mlt(0) + cnst.tau_np_p_mlt(1)*m_Smain_asc;
    m_tau_np_n = cnst.tau_np_n_mlt(0) + cnst.tau_np_n_mlt(1)*m_Smain_asc;
    m_tau_np_p_tl = 0.0;
    m_tau_np_n_tl = 0.0;
  }else{ /*  single track vertex */
    m_fd = keeptagl ? 1. : cnst.fd_np_sgl();
    m_fp = cnst.fp_np_sgl();
    m_tau_np_p    = cnst.tau_np_p_sgl(0) + cnst.tau_np_p_sgl(1)*m_Smain_asc;
    m_tau_np_n    = cnst.tau_np_n_sgl(0) + cnst.tau_np_n_sgl(1)*m_Smain_asc;
    m_tau_np_p_tl = cnst.tau_np_p_sgl(0) + cnst.tau_np_p_sgl(1)*m_Stail_asc;
    m_tau_np_n_tl = cnst.tau_np_n_sgl(0) + cnst.tau_np_n_sgl(1)*m_Stail_asc;
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    TTools::constraint(m_tau_np_p_tl, __POSTIVEPARAM_MIN__);
    TTools::constraint(m_tau_np_n_tl, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  TTools::constraint(m_tau_np_p,  __POSTIVEPARAM_MIN__);
  TTools::constraint(m_tau_np_n,  __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  TTools::constraint(m_fp,0.,1.);
  TTools::constraint(m_fd,0.,1.);
  return;
}

void RascRnpPars::Rnp_param(const ResConst& cnst,const RdetVar& var, const int keeptagl){
  //if((param->fn_np_mlt)[0]>DTRES_EXTERAM_THRE){
  if(cnst.fn_np_mlt() >= -0.00001) Rnp_param_10(cnst,var,keeptagl);
  else                           Rnp_param_03(cnst,var,keeptagl);
  return;
}

void RascRnpPars::Rnp_param_10(const ResConst& cnst,const RdetVar& var, const int keeptagl){
  if(var.ntrk()>1){ /*  multiple track vertex  */
    double st_asc_forf = (m_st_asc > cnst.rnp_kink_st() ? cnst.rnp_kink_st() : m_st_asc);
    double xi_asc_forf = (m_xi_asc > cnst.rnp_kink_xi() ? cnst.rnp_kink_xi() : m_xi_asc);

    double f_neg = cnst.fn_np_mlt();
    double f_delta = cnst.fd_np_mlt(keeptagl) + cnst.fd_np_st_mlt()*st_asc_forf + cnst.fd_np_xi_mlt()*xi_asc_forf + cnst.fd_np_stxi_mlt()*st_asc_forf*xi_asc_forf;

    TTools::constraint(f_delta, 0.0, 1.0);
    TTools::constraint(f_neg, 0.0, 1.0);

    m_fd = (1. - f_neg) * f_delta;
    m_fp = (m_fd < 1. ? (1. - f_neg)*(1. - f_delta)/(1. - m_fd) : 0.5);
    m_tau_np_p = cnst.Snp_global()*(cnst.tau_np_p_mlt(0) + cnst.tau_np_p_mlt(1)*m_st_asc + cnst.tau_np_p_xi_mlt()*m_xi_asc + cnst.tau_np_p_stxi_mlt()*m_st_asc*m_xi_asc);
    m_tau_np_n = cnst.Snp_global()*(cnst.tau_np_n_mlt(0) + cnst.tau_np_n_mlt(1)*m_st_asc + cnst.tau_np_n_xi_mlt()*m_xi_asc + cnst.tau_np_n_stxi_mlt()*m_st_asc*m_xi_asc);
  }else{ /*  single track vertex */
    m_fd = keeptagl ? 1. : cnst.fd_np_sgl();
    m_fp = cnst.fp_np_sgl();
    m_tau_np_p = cnst.Snp_global()*(cnst.tau_np_p_sgl(0) + cnst.tau_np_p_sgl(1)*m_st_asc);
    m_tau_np_n = cnst.Snp_global()*(cnst.tau_np_n_sgl(0) + cnst.tau_np_n_sgl(1)*m_st_asc);
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  TTools::constraint(m_tau_np_p,  __POSTIVEPARAM_MIN__);
  TTools::constraint(m_tau_np_n,  __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  TTools::constraint(m_fp, 0.0, 1.0);
  m_tau_np_p_tl = m_tau_np_p;
  m_tau_np_n_tl = m_tau_np_n;
  return;
}

void RascRnpPars::Rnp_param_03(const ResConst& cnst,const RdetVar& var, const int keeptagl){
  if(var.ntrk()>1){ /*  multiple track vertex  */
    const double Smain_asc_in = (cnst.Snp() <= DTRES_EXTERAM_THRE) ? m_Smain_asc : (1.0+cnst.Snp()*m_xi_asc)*m_st_asc;
    m_fd = cnst.fd_np_mlt(keeptagl);
    m_fp = cnst.fp_np_mlt();
    m_tau_np_p = cnst.Snp_global()*(cnst.tau_np_p_mlt(0) + cnst.tau_np_p_mlt(1)*Smain_asc_in);
    m_tau_np_n = cnst.Snp_global()*(cnst.tau_np_n_mlt(0) + cnst.tau_np_n_mlt(1)*Smain_asc_in);
    const double Stail_asc_in = (cnst.Snp() <= DTRES_EXTERAM_THRE) ? m_Stail_asc : (1.0+cnst.Snp()*m_xi_asc)*m_st_asc;
    m_tau_np_p_tl = cnst.Snp_global()*(cnst.tau_np_p_mlt(0) + cnst.tau_np_p_mlt(1)*Stail_asc_in);
    m_tau_np_n_tl = cnst.Snp_global()*(cnst.tau_np_n_mlt(0) + cnst.tau_np_n_mlt(1)*Stail_asc_in);
  }else{ /*  single track vertex */
    const double Smain_asc_in = (cnst.Snp() <= DTRES_EXTERAM_THRE) ? m_Smain_asc : m_st_asc;
    const double Stail_asc_in = (cnst.Snp() <= DTRES_EXTERAM_THRE) ? m_Stail_asc : m_st_asc;
    m_fd = keeptagl ? 1. : cnst.fd_np_sgl();
    m_fp = cnst.fp_np_sgl();
    m_tau_np_p    = cnst.Snp_global()*(cnst.tau_np_p_sgl(0) + cnst.tau_np_p_sgl(1)*Smain_asc_in);
    m_tau_np_n    = cnst.Snp_global()*(cnst.tau_np_n_sgl(0) + cnst.tau_np_n_sgl(1)*Smain_asc_in);
    m_tau_np_p_tl = cnst.Snp_global()*(cnst.tau_np_p_sgl(0) + cnst.tau_np_p_sgl(1)*Stail_asc_in);
    m_tau_np_n_tl = cnst.Snp_global()*(cnst.tau_np_n_sgl(0) + cnst.tau_np_n_sgl(1)*Stail_asc_in);
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  TTools::constraint(m_tau_np_p,  __POSTIVEPARAM_MIN__);
  TTools::constraint(m_tau_np_n,  __POSTIVEPARAM_MIN__);
  TTools::constraint(m_tau_np_p_tl, __POSTIVEPARAM_MIN__);
  TTools::constraint(m_tau_np_n_tl, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  TTools::constraint(m_fp,0.0,1.0);
  TTools::constraint(m_fd,0.0,1.0);
  return;
}

void RascRnpPars::swap_rnp_param(void){
  m_fp = 1.0 - m_fp;
  const double taunp_n_tmp    = m_tau_np_n;
  m_tau_np_n = m_tau_np_p;
  m_tau_np_p = taunp_n_tmp;
  const double taunp_n_tl_tmp = m_tau_np_n_tl;
  m_tau_np_n_tl = m_tau_np_p_tl;
  m_tau_np_p_tl = taunp_n_tl_tmp;
  return;
}

