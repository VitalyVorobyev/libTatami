#ifndef RESCONST_H
#define RESCONST_H

#include <iostream>
#include <string>

class ResConst{
public:
  ResConst(const std::string& fname);
  int LoadParameters(const std::string& fname);

//  void Set_ftl_rec_mlt(const double& x1, const double& x2) {ftl_rec_mlt[0]    = x1; ftl_rec_mlt[1]  = x2; return;}
  void Set_Srec(const double& x1, const double& x2)        {m_Srec[0]           = x1; m_Srec[1]         = x2; return;}
//  void Set_Stl_rec_mlt(const double& x)                    {Stl_rec_mlt       = x;                        return;}
  void Set_ftl_rec(const double& x)                        {m_ftl_rec           = x;                        return;}
  void Set_Smn_rec(const double& x)                        {m_Smn_rec           = x;                        return;}
  void Set_Stl_rec(const double& x)                        {m_Stl_rec           = x;                        return;}

//  void Set_ftl_asc_mlt(const double& x1, const double& x2) {ftl_asc_mlt[0]    = x1; ftl_asc_mlt[1]  = x2; return;}
  void Set_Sasc(const double& x1, const double& x2)        {m_Sasc[0]           = x1; m_Sasc[1]         = x2; return;}
//  void Set_Stl_asc_mlt(const double& x)                    {Stl_asc_mlt       = x;                        return;}
  void Set_ftl_asc(const double& x)                        {m_ftl_asc           = x;                        return;}
  void Set_Smn_asc(const double& x)                        {m_Smn_asc           = x;                        return;}
  void Set_Stl_asc(const double& x)                        {m_Stl_asc           = x;                        return;}

//  double Get_ftl_rec_mlt(const int i) const {return ftl_rec_mlt[i];}
  double Srec(const int i)        const {return m_Srec[i];}
//  double Get_Stl_rec_mlt(void)        const {return Stl_rec_mlt;}
  double ftl_rec(void)            const {return m_ftl_rec;}
  double Smn_rec(void)            const {return m_Smn_rec;}
  double Stl_rec(void)            const {return m_Stl_rec;}

//  double Get_ftl_asc_mlt(const int i) const {return ftl_asc_mlt[i];}
  double Sasc(const int i)        const {return m_Sasc[i];}
//  double Get_Stl_asc_mlt(void)        const {return Stl_asc_mlt;}
  double ftl_asc(void)            const {return m_ftl_asc;}
  double Smn_asc(void)            const {return m_Smn_asc;}
  double Stl_asc(void)            const {return m_Stl_asc ;}

  // * NP tracks *
  void Set_fd_np_mlt(const double& x1, const double& x2)   {m_fd_np_mlt[0]      = x1; m_fd_np_mlt[1]    = x2; return;}
  void Set_fp_np_mlt(const double& x)                      {m_fp_np_mlt         = x;                        return;}
  void Set_tau_np_p_mlt(const double& x1, const double& x2){m_tau_np_p_mlt[0]   = x1; m_tau_np_p_mlt[1] = x2; return;}
  void Set_tau_np_n_mlt(const double& x1, const double& x2){m_tau_np_n_mlt[0]   = x1; m_tau_np_n_mlt[1] = x2; return;}
  void Set_fd_np_sgl(const double& x)                      {m_fd_np_sgl         = x;                        return;}
  void Set_fp_np_sgl(const double& x)                      {m_fp_np_sgl         = x;                        return;}
  void Set_tau_np_p_sgl(const double& x1, const double& x2){m_tau_np_p_sgl[0]   = x1; m_tau_np_p_sgl[1] = x2; return;}
  void Set_tau_np_n_sgl(const double& x1, const double& x2){m_tau_np_n_sgl[0]   = x1; m_tau_np_n_sgl[1] = x2; return;}
  void Set_fn_np_mlt(const double& x)                      {m_fn_np_mlt         = x;                        return;}
  void Set_rnp_kink_st(const double& x)                    {m_rnp_kink_st       = x;                        return;}
  void Set_rnp_kink_xi(const double& x)                    {m_rnp_kink_xi       = x;                        return;}
  void Set_fd_np_st_mlt(const double& x)                   {m_fd_np_st_mlt      = x;                        return;}
  void Set_fd_np_xi_mlt(const double& x)                   {m_fd_np_xi_mlt      = x;                        return;}
  void Set_fd_np_stxi_mlt(const double& x)                 {m_fd_np_stxi_mlt    = x;                        return;}
  void Set_Snp_global(const double& x)                     {m_Snp_global        = x;                        return;}
  void Set_tau_np_p_xi_mlt(const double& x)                {m_tau_np_p_xi_mlt   = x;                        return;}
  void Set_tau_np_p_stxi_mlt(const double& x)              {m_tau_np_p_stxi_mlt = x;                        return;}
  void Set_tau_np_n_xi_mlt(const double& x)                {m_tau_np_n_xi_mlt   = x;                        return;}
  void Set_tau_np_n_stxi_mlt(const double& x)              {m_tau_np_n_stxi_mlt = x;                        return;}
  void Set_Snp(const double& x)                            {m_Snp               = x;                        return;}

  double fd_np_mlt(const int i)    const {return m_fd_np_mlt[i];}
  double fp_np_mlt(void)           const {return m_fp_np_mlt;}
  double tau_np_p_mlt(const int i) const {return m_tau_np_p_mlt[i];}
  double tau_np_n_mlt(const int i) const {return m_tau_np_n_mlt[i];}
  double fd_np_sgl(void)           const {return m_fd_np_sgl;}
  double fp_np_sgl(void)           const {return m_fp_np_sgl;}
  double tau_np_p_sgl(const int i) const {return m_tau_np_p_sgl[i];}
  double tau_np_n_sgl(const int i) const {return m_tau_np_n_sgl[i];}
  double fn_np_mlt(void)           const {return m_fn_np_mlt;}
  double rnp_kink_st(void)         const {return m_rnp_kink_st;}
  double rnp_kink_xi(void)         const {return m_rnp_kink_xi;}
  double fd_np_st_mlt(void)        const {return m_fd_np_st_mlt;}
  double fd_np_xi_mlt(void)        const {return m_fd_np_xi_mlt;}
  double fd_np_stxi_mlt(void)      const {return m_fd_np_stxi_mlt;}
  double Snp_global(void)          const {return m_Snp_global;}
  double tau_np_p_xi_mlt(void)     const {return m_tau_np_p_xi_mlt;}
  double tau_np_p_stxi_mlt(void)   const {return m_tau_np_p_stxi_mlt;}
  double tau_np_n_xi_mlt(void)     const {return m_tau_np_n_xi_mlt;}
  double tau_np_n_stxi_mlt(void)   const {return m_tau_np_n_stxi_mlt;}
  double Snp(void)                 const {return m_Snp;}

// Otlr
  void Set_f_ol_sgl(const double& x) {m_f_ol_sgl = x; return;}
  void Set_f_ol_mlt(const double& x) {m_f_ol_mul = x; return;}
  void Set_sigma_ol(const double& x) {m_sigma_ol = x; return;}

  double f_ol_sgl(void) const {return m_f_ol_sgl;}
  double f_ol_mlt(void) const {return m_f_ol_mul;}
  double sigma_ol(void) const {return m_sigma_ol;}

protected:
  // ** Resolution parameters **
  // * Rec side *
//  double ftl_rec_mlt[2];
  double m_Srec[2];
//  double Stl_rec_mlt;
  double m_ftl_rec;
  double m_Smn_rec;
  double m_Stl_rec;

  // * Tag side *
//  double ftl_asc_mlt[2];
  double m_Sasc[2];
//  double Stl_asc_mlt;
  double m_ftl_asc;
  double m_Smn_asc;
  double m_Stl_asc;

  // * NP tracks *
  double m_fd_np_mlt[2];
  double m_fp_np_mlt;
  double m_tau_np_p_mlt[2];
  double m_tau_np_n_mlt[2];
  double m_fd_np_sgl;
  double m_fp_np_sgl;
  double m_tau_np_p_sgl[2];
  double m_tau_np_n_sgl[2];
  double m_fn_np_mlt;
  double m_rnp_kink_st;
  double m_rnp_kink_xi;
  double m_fd_np_st_mlt;
  double m_fd_np_xi_mlt;
  double m_fd_np_stxi_mlt;
  double m_Snp_global;
  double m_tau_np_p_xi_mlt;
  double m_tau_np_p_stxi_mlt;
  double m_tau_np_n_xi_mlt;
  double m_tau_np_n_stxi_mlt;
  double m_Snp;

  // * Outlayer *
  double m_sigma_ol;
  double m_f_ol_sgl;
  double m_f_ol_mul;
};

#endif
