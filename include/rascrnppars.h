/** Copyright 2016 Vitaly Vorobyev
 ** @file rascrnppars.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#pragma once

namespace libTatami {

class ResConst;
class RdetVar;

///
/// \brief The RascRnpPars class calculates and keeps auxiliary parameters
/// for asc side vertex resolution
///
class RascRnpPars {
    ///
    /// \brief calc_vtxparam_asc
    /// \param var
    ///
    void calc_vtxparam_asc(const RdetVar& var);
    ///
    /// \brief Rasc_param
    /// \param pars
    /// \param var
    ///
    void Rasc_param(const ResConst& pars, const RdetVar& var);
    ///
    /// \brief Rnp_paramOld
    /// \param cnst
    /// \param var
    /// \param keeptagl
    ///
    void Rnp_paramOld(const ResConst& cnst, const RdetVar& var, int keeptagl);
    ///
    /// \brief Rnp_param
    /// \param cnst
    /// \param var
    /// \param keeptagl
    ///
    void Rnp_param(const ResConst& cnst, const RdetVar& var, int keeptagl);
    ///
    /// \brief Rnp_param_03
    /// \param cnst
    /// \param var
    /// \param keeptagl
    ///
    void Rnp_param_03(const ResConst& cnst, const RdetVar& var, int keeptagl);
    ///
    /// \brief Rnp_param_10
    /// \param cnst
    /// \param var
    /// \param keeptagl
    ///
    void Rnp_param_10(const ResConst& cnst, const RdetVar& var, int keeptagl);

  // Rasc //
    ///
    /// \brief m_Smain_asc
    ///
    double m_Smain_asc;
    ///
    /// \brief m_Stail_asc
    ///
    double m_Stail_asc;
    ///
    /// \brief m_ftail_asc
    ///
    double m_ftail_asc;
    ///
    /// \brief m_mu_main_asc
    ///
    double m_mu_main_asc;
    ///
    /// \brief m_mu_tail_asc
    ///
    double m_mu_tail_asc;
    ///
    /// \brief m_xi_asc
    ///
    double m_xi_asc;
    ///
    /// \brief m_st_asc
    ///
    double m_st_asc;

  // Rnp //
    ///
    /// \brief m_fd
    ///
    double m_fd;
    ///
    /// \brief m_fp
    ///
    double m_fp;
    ///
    /// \brief m_tau_np_p
    ///
    double m_tau_np_p;
    ///
    /// \brief m_tau_np_n
    ///
    double m_tau_np_n;
    ///
    /// \brief m_tau_np_p_tl
    ///
    double m_tau_np_p_tl;
    ///
    /// \brief m_tau_np_n_tl
    ///
    double m_tau_np_n_tl;

 public:
    ///
    /// \brief RascRnpPars
    ///
    RascRnpPars(void);
    ///
    /// \brief Calculate. Calculate parameters
    /// \param cnst
    /// \param var
    /// \param keeptagl
    ///
    void Calculate(const ResConst& cnst, const RdetVar& var, int keeptagl) {
        calc_vtxparam_asc(var);
        Rasc_param(cnst, var);
        Rnp_param(cnst, var, keeptagl);
    }
    ///
    /// \brief swap_rnp_param
    ///
    void swap_rnp_param(void);
    ///
    /// \brief Smain_asc
    /// \return
    ///
    double Smain_asc() const {return m_Smain_asc;}
    ///
    /// \brief Stail_asc
    /// \return
    ///
    double Stail_asc() const {return m_Stail_asc;}
    ///
    /// \brief ftail_asc
    /// \return
    ///
    double ftail_asc() const {return m_ftail_asc;}
    ///
    /// \brief mu_main_asc
    /// \return
    ///
    double mu_main_asc() const {return m_mu_main_asc;}
    ///
    /// \brief mu_tail_asc
    /// \return
    ///
    double mu_tail_asc() const {return m_mu_tail_asc;}
    ///
    /// \brief xi_asc
    /// \return
    ///
    double xi_asc() const {return m_xi_asc;}
    ///
    /// \brief st_asc
    /// \return
    ///
    double st_asc() const {return m_st_asc;}
  // Rnp //
    ///
    /// \brief fd
    /// \return
    ///
    double fd() const {return m_fd;}
    ///
    /// \brief fp
    /// \return
    ///
    double fp() const {return m_fp;}
    ///
    /// \brief tau_np_p
    /// \return
    ///
    double tau_np_p() const {return m_tau_np_p;}
    ///
    /// \brief tau_np_n
    /// \return
    ///
    double tau_np_n() const {return m_tau_np_n;}
    ///
    /// \brief tau_np_p_tl
    /// \return
    ///
    double tau_np_p_tl() const {return m_tau_np_p_tl;}
    ///
    /// \brief tau_np_n_tl
    /// \return
    ///
    double tau_np_n_tl() const {return m_tau_np_n_tl;}
};

}  // namespace libTatami
