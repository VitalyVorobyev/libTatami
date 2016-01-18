#ifndef RRECPARS_H
#define RRECPARS_H

#include "ResVar.h"
#include "ResConst.h"

///
/// \brief The RrecPars class calculates and keeps auxiliary parameters for rec side vertex resolution
///

class RrecPars{
public:
  RrecPars(void);
  void Calculate(const ResConst& cnst,const RdetVar& var){
    calc_vtxparam_rec(var); Rrec_param(cnst,var);
  }

  double Smain_rec()   const {return m_Smain_rec;}
  double Stail_rec()   const {return m_Stail_rec;}
  double ftail_rec()   const {return m_ftail_rec;}
  double mu_main_rec() const {return m_mu_main_rec;}
  double mu_tail_rec() const {return m_mu_tail_rec;}
private:
  void calc_vtxparam_rec(const RdetVar& var);
  void Rrec_param(const ResConst& cnst,const RdetVar& var);

  double m_Smain_rec;
  double m_Stail_rec;
  double m_ftail_rec;
  double m_mu_main_rec;
  double m_mu_tail_rec;
  double m_xi_rec;
  double m_st_rec;
};

#endif // RRECPARS_H
