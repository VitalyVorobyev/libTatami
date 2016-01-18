#include "rkparam.h"
#include "ttools.h"

#include <cmath>
#include <iostream>

using namespace std;

RkPar::RkPar(void):
  m_ak(1.), m_ck(0.), m_r_ckak(0.),
  m_ndm_n(1.), m_ndm_p(1.), m_ndmtau(1.),
  m_cktau(1.), m_ntau_n(1.), m_ntau_p(1.),
  m_fact_n_e(0.), m_fact_p_e(0.), m_fact_am(0.)
{
}

int RkPar::SetAkCk(const double& costh, const double& ecm, const double& tau, const double& dm){
  m_ak = ecm/TTools::mbzero;
  const double pcm = sqrt(ecm*ecm-TTools::mbzero*TTools::mbzero);
  m_ck = pcm*costh/(TTools::beta*TTools::mbzero);
  if(m_ak==0.0 || m_ak==m_ck){
    cout << "RkPar::SetAkCk: Invalid ak = " << m_ak << ", ck = " << m_ck << ", where they should be ak!=0.0&&ak!=ck." << endl;
    return -1;
  }// else{
//    cout << "SetAkCk: ak = " << ak << ", ck = " << ck << endl;
//  }

  m_r_ckak = m_ck/m_ak;
  m_ndmtau = m_r_ckak*dm*tau;
  m_ndm_n  = dm/(m_ak-m_ck);
  m_ndm_p  = dm/(m_ak+m_ck);
  m_ntau_n = tau*(m_ak-m_ck);
  m_ntau_p = tau*(m_ak+m_ck);
  m_fact_n_e = 0.5*(1-m_r_ckak);
  m_fact_p_e = 0.5*(1+m_r_ckak);
  const double inv_ak = 1.0/m_ak;
  const double fact = 1.0/(1.0+m_ndmtau*m_ndmtau);
  m_fact_am = inv_ak*fact;
  return 0;
}

