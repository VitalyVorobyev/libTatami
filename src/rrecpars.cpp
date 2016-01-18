#include "rrecpars.h"
#include "ttools.h"

RrecPars::RrecPars(void):
  m_Smain_rec(1.), m_Stail_rec(1.), m_ftail_rec(0.),
  m_mu_main_rec(0.), m_mu_tail_rec(0.)
{
}

void RrecPars::Rrec_param(const ResConst& cnst,const RdetVar& var){
  if(var.ntrk()>1){ /*  multiple track vertex  */
    m_ftail_rec = 0;//ftl_rec_mlt[0] + ftl_rec_mlt[1]*xi_rec;
    m_Smain_rec = (cnst.Srec(0) + cnst.Srec(1)*m_xi_rec)*m_st_rec;
    m_Stail_rec = m_Smain_rec;
    TTools::constraint(m_ftail_rec, 0.0, 1.0);
//#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
//    constraint(m_Smain_rec, __POSTIVEPARAM_MIN__);
//#endif /* __FORCE_PARAM_CONSTRAINT__ */
  } else{ /*  single track vertex */
    m_ftail_rec = cnst.ftl_rec();
    m_Smain_rec = cnst.Smn_rec()*m_st_rec;
    m_Stail_rec = cnst.Stl_rec()*m_st_rec;
    TTools::constraint(m_ftail_rec, 0.0, 1.0);
//#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
//    constraint(m_Smain_rec, __POSTIVEPARAM_MIN__);
//    constraint(m_Stail_rec, __POSTIVEPARAM_MIN__);
//#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
  m_mu_main_rec = 0.0;
  m_mu_tail_rec = 0.0;
  return;
}

void RrecPars::calc_vtxparam_rec(const RdetVar& var){
  m_xi_rec = (var.ndf() > 0 ? var.chisq()/var.ndf() : 1);
  m_st_rec = var.sz()*TTools::cm2ps;
  return;
}

