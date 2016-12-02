/** Copyright 2016 Vitaly Vorobyev
 * @file ResConst.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#ifndef INCLUDE_RESCONST_H_
#define INCLUDE_RESCONST_H_

#include <iostream>
#include <string>

#include "./typedefs.h"

namespace libTatami {

///
/// \brief The ResConst class
///
class ResConst {
 public:
    ///
    /// \brief ResConst
    /// \param fname
    ///
    explicit ResConst(const str& fname);
    ///
    /// \brief LoadParameters
    /// \param fname
    /// \return
    ///
    int LoadParameters(const str& fname);
    ///
    /// \brief Set_Srec
    /// \param x1
    /// \param x2
    ///
    void Set_Srec(const double& x1, const double& x2) {
        m_Srec[0] = x1; m_Srec[1] = x2;
    }
    ///
    /// \brief Set_ftl_rec
    /// \param x
    ///
    void Set_ftl_rec(const double& x) {m_ftl_rec = x;}
    ///
    /// \brief Set_Smn_rec
    /// \param x
    ///
    void Set_Smn_rec(const double& x) {m_Smn_rec = x;}
    ///
    /// \brief Set_Stl_rec
    /// \param x
    ///
    void Set_Stl_rec(const double& x) {m_Stl_rec = x;}
    ///
    /// \brief Set_Sasc
    /// \param x1
    /// \param x2
    ///
    void Set_Sasc(const double& x1, const double& x2) {
        m_Sasc[0] = x1; m_Sasc[1] = x2;
    }
    ///
    /// \brief Set_ftl_asc
    /// \param x
    ///
    void Set_ftl_asc(const double& x) {m_ftl_asc = x;}
    ///
    /// \brief Set_Smn_asc
    /// \param x
    ///
    void Set_Smn_asc(const double& x) {m_Smn_asc = x;}
    ///
    /// \brief Set_Stl_asc
    /// \param x
    ///
    void Set_Stl_asc(const double& x) {m_Stl_asc = x;}
    ///
    /// \brief Srec
    /// \param i
    /// \return
    ///
    double Srec(cint i) const {return m_Srec[i];}
    ///
    /// \brief ftl_rec
    /// \return
    ///
    double ftl_rec(void) const {return m_ftl_rec;}
    ///
    /// \brief Smn_rec
    /// \return
    ///
    double Smn_rec(void) const {return m_Smn_rec;}
    ///
    /// \brief Stl_rec
    /// \return
    ///
    double Stl_rec(void) const {return m_Stl_rec;}
    ///
    /// \brief Sasc
    /// \param i
    /// \return
    ///
    double Sasc(cint i) const {return m_Sasc[i];}
    ///
    /// \brief ftl_asc
    /// \return
    ///
    double ftl_asc(void) const {return m_ftl_asc;}
    ///
    /// \brief Smn_asc
    /// \return
    ///
    double Smn_asc(void) const {return m_Smn_asc;}
    ///
    /// \brief Stl_asc
    /// \return
    ///
    double Stl_asc(void) const {return m_Stl_asc ;}

  // * NP tracks *
    ///
    /// \brief Set_fd_np_mlt
    /// \param x1
    /// \param x2
    ///
    void Set_fd_np_mlt(const double& x1, const double& x2) {
        m_fd_np_mlt[0] = x1; m_fd_np_mlt[1] = x2;
    }
    ///
    /// \brief Set_fp_np_mlt
    /// \param x
    ///
    void Set_fp_np_mlt(const double& x) {m_fp_np_mlt = x;}
    ///
    /// \brief Set_tau_np_p_mlt
    /// \param x1
    /// \param x2
    ///
    void Set_tau_np_p_mlt(const double& x1, const double& x2) {
        m_tau_np_p_mlt[0] = x1; m_tau_np_p_mlt[1] = x2;
    }
    ///
    /// \brief Set_tau_np_n_mlt
    /// \param x1
    /// \param x2
    ///
    void Set_tau_np_n_mlt(const double& x1, const double& x2) {
        m_tau_np_n_mlt[0] = x1; m_tau_np_n_mlt[1] = x2;
    }
    ///
    /// \brief Set_fd_np_sgl
    /// \param x
    ///
    void Set_fd_np_sgl(const double& x) {m_fd_np_sgl = x;}
    ///
    /// \brief Set_fp_np_sgl
    /// \param x
    ///
    void Set_fp_np_sgl(const double& x) {m_fp_np_sgl = x;}
    ///
    /// \brief Set_tau_np_p_sgl
    /// \param x1
    /// \param x2
    ///
    void Set_tau_np_p_sgl(const double& x1, const double& x2) {
        m_tau_np_p_sgl[0] = x1; m_tau_np_p_sgl[1] = x2;
    }
    ///
    /// \brief Set_tau_np_n_sgl
    /// \param x1
    /// \param x2
    ///
    void Set_tau_np_n_sgl(const double& x1, const double& x2) {
        m_tau_np_n_sgl[0] = x1; m_tau_np_n_sgl[1] = x2;
    }
    ///
    /// \brief Set_fn_np_mlt
    /// \param x
    ///
    void Set_fn_np_mlt(const double& x) {m_fn_np_mlt = x;}
    ///
    /// \brief Set_rnp_kink_st
    /// \param x
    ///
    void Set_rnp_kink_st(const double& x) {m_rnp_kink_st = x;}
    ///
    /// \brief Set_rnp_kink_xi
    /// \param x
    ///
    void Set_rnp_kink_xi(const double& x) {m_rnp_kink_xi = x;}
    ///
    /// \brief Set_fd_np_st_mlt
    /// \param x
    ///
    void Set_fd_np_st_mlt(const double& x) {m_fd_np_st_mlt = x;}
    ///
    /// \brief Set_fd_np_xi_mlt
    /// \param x
    ///
    void Set_fd_np_xi_mlt(const double& x) {m_fd_np_xi_mlt = x;}
    ///
    /// \brief Set_fd_np_stxi_mlt
    /// \param x
    ///
    void Set_fd_np_stxi_mlt(const double& x) {m_fd_np_stxi_mlt = x;}
    ///
    /// \brief Set_Snp_global
    /// \param x
    ///
    void Set_Snp_global(const double& x) {m_Snp_global = x;}
    ///
    /// \brief Set_tau_np_p_xi_mlt
    /// \param x
    ///
    void Set_tau_np_p_xi_mlt(const double& x) {m_tau_np_p_xi_mlt = x;}
    ///
    /// \brief Set_tau_np_p_stxi_mlt
    /// \param x
    ///
    void Set_tau_np_p_stxi_mlt(const double& x) {m_tau_np_p_stxi_mlt = x;}
    ///
    /// \brief Set_tau_np_n_xi_mlt
    /// \param x
    ///
    void Set_tau_np_n_xi_mlt(const double& x) {m_tau_np_n_xi_mlt = x;}
    ///
    /// \brief Set_tau_np_n_stxi_mlt
    /// \param x
    ///
    void Set_tau_np_n_stxi_mlt(const double& x) {m_tau_np_n_stxi_mlt = x;}
    ///
    /// \brief Set_Snp
    /// \param x
    ///
    void Set_Snp(const double& x) {m_Snp = x;}
    ///
    /// \brief fd_np_mlt
    /// \param i
    /// \return
    ///
    double fd_np_mlt(cint i) const {return m_fd_np_mlt[i];}
    ///
    /// \brief fp_np_mlt
    /// \return
    ///
    double fp_np_mlt(void) const {return m_fp_np_mlt;}
    ///
    /// \brief tau_np_p_mlt
    /// \param i
    /// \return
    ///
    double tau_np_p_mlt(cint i) const {return m_tau_np_p_mlt[i];}
    ///
    /// \brief tau_np_n_mlt
    /// \param i
    /// \return
    ///
    double tau_np_n_mlt(cint i) const {return m_tau_np_n_mlt[i];}
    ///
    /// \brief fd_np_sgl
    /// \return
    ///
    double fd_np_sgl(void) const {return m_fd_np_sgl;}
    ///
    /// \brief fp_np_sgl
    /// \return
    ///
    double fp_np_sgl(void) const {return m_fp_np_sgl;}
    ///
    /// \brief tau_np_p_sgl
    /// \param i
    /// \return
    ///
    double tau_np_p_sgl(cint i) const {return m_tau_np_p_sgl[i];}
    ///
    /// \brief tau_np_n_sgl
    /// \param i
    /// \return
    ///
    double tau_np_n_sgl(cint i) const {return m_tau_np_n_sgl[i];}
    ///
    /// \brief fn_np_mlt
    /// \return
    ///
    double fn_np_mlt(void) const {return m_fn_np_mlt;}
    ///
    /// \brief rnp_kink_st
    /// \return
    ///
    double rnp_kink_st(void) const {return m_rnp_kink_st;}
    ///
    /// \brief rnp_kink_xi
    /// \return
    ///
    double rnp_kink_xi(void) const {return m_rnp_kink_xi;}
    ///
    /// \brief fd_np_st_mlt
    /// \return
    ///
    double fd_np_st_mlt(void) const {return m_fd_np_st_mlt;}
    ///
    /// \brief fd_np_xi_mlt
    /// \return
    ///
    double fd_np_xi_mlt(void) const {return m_fd_np_xi_mlt;}
    ///
    /// \brief fd_np_stxi_mlt
    /// \return
    ///
    double fd_np_stxi_mlt(void) const {return m_fd_np_stxi_mlt;}
    ///
    /// \brief Snp_global
    /// \return
    ///
    double Snp_global(void) const {return m_Snp_global;}
    ///
    /// \brief tau_np_p_xi_mlt
    /// \return
    ///
    double tau_np_p_xi_mlt(void) const {return m_tau_np_p_xi_mlt;}
    ///
    /// \brief tau_np_p_stxi_mlt
    /// \return
    ///
    double tau_np_p_stxi_mlt(void) const {return m_tau_np_p_stxi_mlt;}
    ///
    /// \brief tau_np_n_xi_mlt
    /// \return
    ///
    double tau_np_n_xi_mlt(void) const {return m_tau_np_n_xi_mlt;}
    ///
    /// \brief tau_np_n_stxi_mlt
    /// \return
    ///
    double tau_np_n_stxi_mlt(void) const {return m_tau_np_n_stxi_mlt;}
    ///
    /// \brief Snp
    /// \return
    ///
    double Snp(void) const {return m_Snp;}

  // Outlier
    ///
    /// \brief Set_f_ol_sgl
    /// \param x
    ///
    void Set_f_ol_sgl(const double& x) {m_f_ol_sgl = x; return;}
    ///
    /// \brief Set_f_ol_mlt
    /// \param x
    ///
    void Set_f_ol_mlt(const double& x) {m_f_ol_mul = x; return;}
    ///
    /// \brief Set_sigma_ol
    /// \param x
    ///
    void Set_sigma_ol(const double& x) {m_sigma_ol = x; return;}
    ///
    /// \brief f_ol_sgl
    /// \return
    ///
    double f_ol_sgl(void) const {return m_f_ol_sgl;}
    ///
    /// \brief f_ol_mlt
    /// \return
    ///
    double f_ol_mlt(void) const {return m_f_ol_mul;}
    ///
    /// \brief sigma_ol
    /// \return
    ///
    double sigma_ol(void) const {return m_sigma_ol;}

 protected:
  // ** Resolution parameters **
  // * Rec side *
    ///
    /// \brief m_Srec
    ///
    double m_Srec[2];
    ///
    /// \brief m_ftl_rec
    ///
    double m_ftl_rec;
    ///
    /// \brief m_Smn_rec
    ///
    double m_Smn_rec;
    ///
    /// \brief m_Stl_rec
    ///
    double m_Stl_rec;

  // * Tag side *
    ///
    /// \brief m_Sasc
    ///
    double m_Sasc[2];
    ///
    /// \brief m_ftl_asc
    ///
    double m_ftl_asc;
    ///
    /// \brief m_Smn_asc
    ///
    double m_Smn_asc;
    ///
    /// \brief m_Stl_asc
    ///
    double m_Stl_asc;

  // * NP tracks *
    ///
    /// \brief m_fd_np_mlt
    ///
    double m_fd_np_mlt[2];
    ///
    /// \brief m_fp_np_mlt
    ///
    double m_fp_np_mlt;
    ///
    /// \brief m_tau_np_p_mlt
    ///
    double m_tau_np_p_mlt[2];
    ///
    /// \brief m_tau_np_n_mlt
    ///
    double m_tau_np_n_mlt[2];
    ///
    /// \brief m_fd_np_sgl
    ///
    double m_fd_np_sgl;
    ///
    /// \brief m_fp_np_sgl
    ///
    double m_fp_np_sgl;
    ///
    /// \brief m_tau_np_p_sgl
    ///
    double m_tau_np_p_sgl[2];
    ///
    /// \brief m_tau_np_n_sgl
    ///
    double m_tau_np_n_sgl[2];
    ///
    /// \brief m_fn_np_mlt
    ///
    double m_fn_np_mlt;
    ///
    /// \brief m_rnp_kink_st
    ///
    double m_rnp_kink_st;
    ///
    /// \brief m_rnp_kink_xi
    ///
    double m_rnp_kink_xi;
    ///
    /// \brief m_fd_np_st_mlt
    ///
    double m_fd_np_st_mlt;
    ///
    /// \brief m_fd_np_xi_mlt
    ///
    double m_fd_np_xi_mlt;
    ///
    /// \brief m_fd_np_stxi_mlt
    ///
    double m_fd_np_stxi_mlt;
    ///
    /// \brief m_Snp_global
    ///
    double m_Snp_global;
    ///
    /// \brief m_tau_np_p_xi_mlt
    ///
    double m_tau_np_p_xi_mlt;
    ///
    /// \brief m_tau_np_p_stxi_mlt
    ///
    double m_tau_np_p_stxi_mlt;
    ///
    /// \brief m_tau_np_n_xi_mlt
    ///
    double m_tau_np_n_xi_mlt;
    ///
    /// \brief m_tau_np_n_stxi_mlt
    ///
    double m_tau_np_n_stxi_mlt;
    ///
    /// \brief m_Snp
    ///
    double m_Snp;

  // * Outlayer *
    ///
    /// \brief m_sigma_ol
    ///
    double m_sigma_ol;
    ///
    /// \brief m_f_ol_sgl
    ///
    double m_f_ol_sgl;
    ///
    /// \brief m_f_ol_mul
    ///
    double m_f_ol_mul;
};

}  // namespace libTatami

#endif  // INCLUDE_RESCONST_H_
