/** Copyright 2016 Vitaly Vorobyev
 ** @file rbkgpdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#pragma once

#include "abspdf.h"
#include "bkgpdfparset.h"

namespace libTatami {

///
/// \brief The RbkgPdf class. Parameterization of background dt distribution
///
class RbkgPdf : public AbsPdf {
    ///
    /// \brief Pdf
    /// \param x
    /// \param s
    /// \param ndf
    /// \return
    ///
    double Pdf(double x, double s, int ndf) const;
    ///
    /// \brief DefaultInit
    /// \param mode
    /// \param svd
    /// \param mc
    /// \param tune
    ///
//    void DefaultInit(int mode, int svd, bool mc, bool tune);
    ///
    /// \brief AddOutlier
    /// \param x
    /// \param Lin
    /// \param nLi
    /// \return
    ///
    double AddOutlier(double x, double Lin, double nLi) const;
    ///
    /// \brief cm2ps
    ///
    static const double cm2ps;
    ///
    /// \brief m_pars
    ///
    const BkgPDFParSet m_pars;
    ///
    /// \brief m_scale
    ///
    mutable double m_scale;
    ///
    /// \brief m_shift
    ///
    double m_shift;
    ///
    /// \brief m_sigma. Event-dependent estimation of uncertainty
    ///
    mutable double m_sigma;
    ///
    /// \brief m_ndf. Number of degrees of freedom
    ///
    mutable int m_ndf;

 public:
    ///
    /// \brief RbkgPdf
    /// \param fname
    ///
    explicit RbkgPdf(const std::string& fname);
    ///
    /// \brief RbkgPdf
    ///
    RbkgPdf(void);
    ///
    /// \brief operator ()
    /// \param x
    /// \return
    ///
    double operator()(double x) const override final;
    ///
    /// \brief operator ()
    /// \param evt
    /// \return
    ///
    double operator()(const ICPVEvt& evt) const override final;
    ///
    /// \brief Pars
    /// \return
    ///
    const BkgPDFParSet& Pars(void) const {return m_pars;}
    ///
    /// \brief SetSigma
    /// \param x
    ///
    void SetSigma(double x) {m_sigma = x;}
    ///
    /// \brief SetNDF
    /// \param x
    ///
    void SetNDF(int x) {m_ndf = x;}
    ///
    /// \brief SetScale
    /// \param x
    ///
    void SetScale(double x) {m_scale = x;}
    ///
    /// \brief SetShift
    /// \param x
    ///
    void SetShift(double x) {m_shift = x;}
    ///
    /// \brief SetScaleShift
    /// \param scale
    /// \param shift
    ///
    void SetScaleShift(double scale, double shift) {
        m_scale = scale; m_shift = shift;
    }
    ///
    /// \brief GetScale
    /// \return
    ///
    double GetScale(void) const { return m_scale;}
    ///
    /// \brief GetShift
    /// \return
    ///
    double GetShift(void) const { return m_shift;}
};

}  // namespace libTatami
