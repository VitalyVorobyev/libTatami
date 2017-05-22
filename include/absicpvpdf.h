/** Copyright 2016 Vitaly Vorobyev
 ** @file absicpvpdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_ABSICPVPDF_H_
#define INCLUDE_ABSICPVPDF_H_

#include "./abspdf.h"

namespace libTatami {

///
/// \brief The AbsICPVPdf class. Abstract class for PDF
/// describing CP-violating dt distribution
///
class AbsICPVPdf : public AbsPdf {
 public:
    ///
    /// \brief AbsICPVPdf
    /// \param c
    /// \param s
    /// \param tau
    /// \param dm
    ///
    AbsICPVPdf(const double& tau=1.520, const double& dm=0.505,
               const double& c=0, const double& s=0);
    ///
    /// \brief SetC. Set coefficient near cos(dt)
    /// \param v
    ///
    void SetC(const double& v) {m_c = v;}
    ///
    /// \brief SetS. Set coefficient near sin(dt)
    /// \param v
    ///
    void SetS(const double& v) {m_s = v;}
    /**
     * @brief SetTag
     * @param x
     */
    void SetTag(const int x) {m_tag = x;}
    ///
    /// \brief SetTauDm. Set lifetime and mass difference
    /// \param tau
    /// \param dm
    ///
    void SetTauDm(const double& tau, const double& dm) {
        SetTau(tau); SetDm(dm);
    }
    ///
    /// \brief SetTau. Set lifetime
    /// \param tau
    ///
    void SetTau(const double& tau) {m_tau = tau;}
    ///
    /// \brief SetDm. Set mass difference
    /// \param dm
    ///
    void SetDm(const double& dm) {m_dm = dm;}
    ///
    /// \brief C. Coefficient near cos(dt)
    /// \return
    ///
    double C(void) const {return m_c;}
    ///
    /// \brief S. Coefficient near sin(dt)
    /// \return
    ///
    double S(void) const {return m_s;}
    ///
    /// \brief tau. Lifetime
    /// \return
    ///
    double tau(void) const {return m_tau;}
    ///
    /// \brief dm. Mass difference
    /// \return
    ///
    double dm(void) const {return m_dm;}
    ///
    /// \brief print_params
    ///
    void print_params(void) const;

 protected:
    ///
    /// \brief Lifetime
    ///
    double m_tau;
    ///
    /// \brief Mass difference
    ///
    double m_dm;
    ///
    /// \brief Coefficient near cos(dm*dt)
    ///
    double m_c;
    ///
    /// \brief Coefficient near sin(dm*dt)
    ///
    double m_s;
    /**
     * @brief m_tag
     */
    int m_tag;
};

}  // namespace libTatami

#endif  // INCLUDE_ABSICPVPDF_H_
