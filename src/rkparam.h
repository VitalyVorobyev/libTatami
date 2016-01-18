#ifndef RKPARAM_H
#define RKPARAM_H

class RkPar{
public:
  RkPar(void);

  int SetAkCk(const double& costh, const double& ecm, const double& tau, const double& dm);
  double ak(void)       const {return m_ak;}
  double ck(void)       const {return m_ck;}
  double r_ckak(void)   const {return m_r_ckak;}
  double ndm_n(void)    const {return m_ndm_n;}
  double ndm_p(void)    const {return m_ndm_p;}
  double ndmtau(void)   const {return m_ndmtau;}
  double cktau(void)    const {return m_cktau;}
  double ntau_n(void)   const {return m_ntau_n;}
  double ntau_p(void)   const {return m_ntau_p;}
  double fact_n_e(void) const {return m_fact_n_e;}
  double fact_p_e(void) const {return m_fact_p_e;}
  double fact_am(void)  const {return m_fact_am;}

private:
  double m_ak;
  double m_ck;
  double m_r_ckak;
  double m_ndm_n;
  double m_ndm_p;
  double m_ndmtau;
  double m_cktau;
  double m_ntau_n;
  double m_ntau_p;

  double m_fact_n_e;
  double m_fact_p_e;
  double m_fact_am;
};



#endif

