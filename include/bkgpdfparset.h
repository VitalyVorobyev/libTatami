/** Copyright 2016 Vitaly Vorobyev
 * @file absicpvpdf.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#pragma once

#include <string>

namespace libTatami {

///
/// \brief The BkgPDFParSet class. Parametrs for background dt distribution
/// \todo implement parameters with std::map<std::string, double>
///
class BkgPDFParSet{
    ///
    /// \brief m_S_main_mlt
    ///
    double m_S_main_mlt;
    ///
    /// \brief m_S_tail_mlt
    ///
    double m_S_tail_mlt;
    ///
    /// \brief m_f_tail_mlt
    ///
    double m_f_tail_mlt;
    ///
    /// \brief m_f_delta_mlt
    ///
    double m_f_delta_mlt;
    ///
    /// \brief m_S_main_sgl
    ///
    double m_S_main_sgl;
    ///
    /// \brief m_S_tail_sgl
    ///
    double m_S_tail_sgl;
    ///
    /// \brief m_f_tail_sgl
    ///
    double m_f_tail_sgl;
    ///
    /// \brief m_f_delta_sgl
    ///
    double m_f_delta_sgl;
    ///
    /// \brief m_mu
    ///
    double m_mu;
    ///
    /// \brief m_mu_delta
    ///
    double m_mu_delta;
    ///
    /// \brief m_tau
    ///
    double m_tau;
    ///
    /// \brief m_f_otlr
    ///
    double m_f_otlr;
    ///
    /// \brief m_s_otlr
    ///
    double m_s_otlr;

 public:
    ///
    /// \brief BkgPDFParSet
    ///
    BkgPDFParSet();
    /**
     * @brief BkgPDFParSet
     * @param fname
     */
    explicit BkgPDFParSet(const std::string fname);
    ///
    /// \brief SetTau
    /// \param x
    ///
    void SetTau(double& x) {m_tau = x;}
    ///
    /// \brief Set_f_tail_mlt
    /// \param x
    ///
    void Set_f_tail_mlt(double& x) {m_f_tail_mlt = x;}
    ///
    /// \brief Set_S_main_mlt
    /// \param x
    ///
    void Set_S_main_mlt(double& x) {m_S_main_mlt = x;}
    ///
    /// \brief Set_S_tail_mlt
    /// \param x
    ///
    void Set_S_tail_mlt(double& x) {m_S_tail_mlt = x;}
    ///
    /// \brief Set_f_tail_sgl
    /// \param x
    ///
    void Set_f_tail_sgl(double& x) {m_f_tail_sgl = x;}
    ///
    /// \brief Set_S_main_sgl
    /// \param x
    ///
    void Set_S_main_sgl(double& x) {m_S_main_sgl = x;}
    ///
    /// \brief Set_S_tail_sgl
    /// \param x
    ///
    void Set_S_tail_sgl(double& x) {m_S_tail_sgl = x;}
    ///
    /// \brief Set_f_delta_mlt
    /// \param x
    ///
    void Set_f_delta_mlt(double& x) {m_f_delta_mlt = x;}
    ///
    /// \brief Set_f_delta_sgl
    /// \param x
    ///
    void Set_f_delta_sgl(double& x) {m_f_delta_sgl = x;}
    ///
    /// \brief Set_mu_delta
    /// \param x
    ///
    void Set_mu_delta(double& x) {m_mu_delta = x;}
    ///
    /// \brief Set_mu
    /// \param x
    ///
    void Set_mu(double& x) {m_mu = x;}
    ///
    /// \brief Set_f_otlr
    /// \param x
    ///
    void Set_f_otlr(double& x) {m_f_otlr = x;}
    ///
    /// \brief Set_s_otlr
    /// \param x
    ///
    void Set_s_otlr(double& x) {m_s_otlr = x;}
    ///
    /// \brief tau
    /// \return
    ///
    double tau(void) const {return m_tau;}
    ///
    /// \brief f_tail_mlt
    /// \return
    ///
    double f_tail_mlt(void) const {return m_f_tail_mlt;}
    ///
    /// \brief S_main_mlt
    /// \return
    ///
    double S_main_mlt(void) const {return m_S_main_mlt;}
    ///
    /// \brief S_tail_mlt
    /// \return
    ///
    double S_tail_mlt(void) const {return m_S_tail_mlt;}
    ///
    /// \brief f_tail_sgl
    /// \return
    ///
    double f_tail_sgl(void) const {return m_f_tail_sgl;}
    ///
    /// \brief S_main_sgl
    /// \return
    ///
    double S_main_sgl(void) const {return m_S_main_sgl;}
    ///
    /// \brief S_tail_sgl
    /// \return
    ///
    double S_tail_sgl(void) const {return m_S_tail_sgl;}
    ///
    /// \brief f_delta_mlt
    /// \return
    ///
    double f_delta_mlt(void) const {return m_f_delta_mlt;}
    ///
    /// \brief f_delta_sgl
    /// \return
    ///
    double f_delta_sgl(void) const {return m_f_delta_sgl;}
    ///
    /// \brief mu_delta
    /// \return
    ///
    double mu_delta(void) const {return m_mu_delta;}
    ///
    /// \brief mu
    /// \return
    ///
    double mu(void) const {return m_mu;}
    ///
    /// \brief f_otlr
    /// \return
    ///
    double f_otlr(void) const {return m_f_otlr;}
    ///
    /// \brief s_otlr
    /// \return
    ///
    double s_otlr(void) const {return m_s_otlr;}
    ///
    /// \brief GetParametersFromFile
    /// \param fname
    /// \return
    ///
    int GetParametersFromFile(const std::string& fname);
};

}  // namespace libTatami
