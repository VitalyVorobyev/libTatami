/** Copyright 2016 Vitaly Vorobyev
 ** @file absicpvpdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#pragma once

#include <utility>

#include "abspdf.h"

namespace libTatami {

///
/// \brief The AbsICPVPdf class. Abstract class for PDF
/// describing CP-violating dt distribution
/// \todo make const operator(c,s,...) instead of const SetC() and const SetS()
///
class AbsICPVPdf : public AbsPdf {
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
    mutable double m_c;
    ///
    /// \brief Coefficient near sin(dm*dt)
    ///
    mutable double m_s;
    /**
     * @brief m_tag
     */
    mutable int m_tag;
    /**
     * @brief m_wtag. Wrong raggin probability
     */
    double m_wtag;

 public:
    ///
    /// \brief AbsICPVPdf
    /// \param c
    /// \param s
    /// \param tau
    /// \param dm
    ///
    AbsICPVPdf(double tau=1.520, double dm=0.505, double c=0, double s=0);
    ///
    /// \brief SetC. Set coefficient near cos(dt)
    /// \param v
    ///
    void SetC(double v) const {m_c = v;}
    ///
    /// \brief SetS. Set coefficient near sin(dt)
    /// \param v
    ///
    void SetS(double v) const {m_s = v;}
    /**
     * @brief SetTag
     * @param x
     */
    void SetTag(const int x) const {m_tag = x;}
    void SetCS(double c, double s) const;
    void SetCS(const std::pair<double, double>& x) const;
    /**
     * @brief SetWTag. Set wrong taggin probability
     * @param x
     */
    void SetWTag(double x) {m_wtag = x;}
    ///
    /// \brief SetTauDm. Set lifetime and mass difference
    /// \param tau
    /// \param dm
    ///
    void SetTauDm(double tau, double dm) {
        SetTau(tau); SetDm(dm);
    }
    ///
    /// \brief SetTau. Set lifetime
    /// \param tau
    ///
    void SetTau(double tau) {m_tau = tau;}
    ///
    /// \brief SetDm. Set mass difference
    /// \param dm
    ///
    void SetDm(double dm) {m_dm = dm;}
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
    int tag() const {return m_tag;}

    double wtag(void) const {return m_wtag;}
    ///
    /// \brief print_params
    ///
    void print_params(void) const;
    /** Returns (1 - 2 * wtag) / (1 + dm**2 * tau**2 ) */
    double alpha() const;
};

}  // namespace libTatami
