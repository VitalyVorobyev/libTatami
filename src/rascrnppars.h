#ifndef RASCRNPPARS_H
#define RASCRNPPARS_H

#include "ResVar.h"
#include "ResConst.h"

///
/// \brief The RascRnpPars class calculates and keeps auxiliary parameters for asc side vertex resolution
///

class RascRnpPars{
public:
  RascRnpPars(void);
/// Calculate parameters
  void Calculate(const ResConst& cnst,const RdetVar& var, const int keeptagl){
    calc_vtxparam_asc(var); Rasc_param(cnst,var);
    Rnp_param(cnst,var,keeptagl);
  }
  void swap_rnp_param(void);

  double Smain_asc()   const {return m_Smain_asc;}
  double Stail_asc()   const {return m_Stail_asc;}
  double ftail_asc()   const {return m_ftail_asc;}
  double mu_main_asc() const {return m_mu_main_asc;}
  double mu_tail_asc() const {return m_mu_tail_asc;}
  double xi_asc()      const {return m_xi_asc;}
  double st_asc()      const {return m_st_asc;}
  // Rnp //
  double fd()          const {return m_fd;}
  double fp()          const {return m_fp;}
  double tau_np_p()    const {return m_tau_np_p;}
  double tau_np_n()    const {return m_tau_np_n;}
  double tau_np_p_tl() const {return m_tau_np_p_tl;}
  double tau_np_n_tl() const {return m_tau_np_n_tl;}
private:
  void calc_vtxparam_asc(const RdetVar& var);
  void Rasc_param(  const ResConst& pars,const RdetVar& var);
  void Rnp_paramOld(const ResConst& cnst,const RdetVar& var, const int keeptagl);
  void Rnp_param(   const ResConst& cnst,const RdetVar& var, const int keeptagl);
  void Rnp_param_03(const ResConst& cnst,const RdetVar& var, const int keeptagl);
  void Rnp_param_10(const ResConst& cnst,const RdetVar& var, const int keeptagl);

  // Rasc //
  double m_Smain_asc;
  double m_Stail_asc;
  double m_ftail_asc;
  double m_mu_main_asc;
  double m_mu_tail_asc;
  double m_xi_asc;
  double m_st_asc;

  // Rnp //
  double m_fd;
  double m_fp;
  double m_tau_np_p;
  double m_tau_np_n;
  double m_tau_np_p_tl;
  double m_tau_np_n_tl;
};

#endif
