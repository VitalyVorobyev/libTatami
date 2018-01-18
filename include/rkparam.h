/** Copyright 2016 Vitaly Vorobyev
 ** @file rkparam.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **
 **/

#pragma once

namespace libTatami {

class RkPar {
    ///
    /// \brief m_ak
    ///
    double m_ak;
    ///
    /// \brief m_ck
    ///
    double m_ck;
    ///
    /// \brief m_r_ckak
    ///
    double m_r_ckak;
    ///
    /// \brief m_ndm_n
    ///
    double m_ndm_n;
    ///
    /// \brief m_ndm_p
    ///
    double m_ndm_p;
    ///
    /// \brief m_ndmtau
    ///
    double m_ndmtau;
    ///
    /// \brief m_cktau
    ///
    double m_cktau;
    ///
    /// \brief m_ntau_n
    ///
    double m_ntau_n;
    ///
    /// \brief m_ntau_p
    ///
    double m_ntau_p;
    ///
    /// \brief m_fact_n_e
    ///
    double m_fact_n_e;
    ///
    /// \brief m_fact_p_e
    ///
    double m_fact_p_e;
    ///
    /// \brief m_fact_am
    ///
    double m_fact_am;

 public:
    ///
    /// \brief RkPar
    ///
    RkPar(void);

    ///
    /// \brief SetAkCk
    /// \param costh
    /// \param ecm
    /// \param tau
    /// \param dm
    /// \return
    ///
    int SetAkCk(double costh, double ecm, double tau, double dm);
    ///
    /// \brief ak
    /// \return
    ///
    double ak(void) const {return m_ak;}
    ///
    /// \brief ck
    /// \return
    ///
    double ck(void) const {return m_ck;}
    ///
    /// \brief r_ckak
    /// \return
    ///
    double r_ckak(void) const {return m_r_ckak;}
    ///
    /// \brief ndm_n
    /// \return
    ///
    double ndm_n(void) const {return m_ndm_n;}
    ///
    /// \brief ndm_p
    /// \return
    ///
    double ndm_p(void) const {return m_ndm_p;}
    ///
    /// \brief ndmtau
    /// \return
    ///
    double ndmtau(void) const {return m_ndmtau;}
    ///
    /// \brief cktau
    /// \return
    ///
    double cktau(void) const {return m_cktau;}
    ///
    /// \brief ntau_n
    /// \return
    ///
    double ntau_n(void) const {return m_ntau_n;}
    ///
    /// \brief ntau_p
    /// \return
    ///
    double ntau_p(void) const {return m_ntau_p;}
    ///
    /// \brief fact_n_e
    /// \return
    ///
    double fact_n_e(void) const {return m_fact_n_e;}
    ///
    /// \brief fact_p_e
    /// \return
    ///
    double fact_p_e(void) const {return m_fact_p_e;}
    ///
    /// \brief fact_am
    /// \return
    ///
    double fact_am(void) const {return m_fact_am;}
};

}  // namespace libTatami
