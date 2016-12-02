/** Copyright 2016 Vitaly Vorobyev
 ** @file rbkgpdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_RBKGPDF_H_
#define INCLUDE_RBKGPDF_H_

#include "./abspdf.h"
#include "./icpvevent.h"
#include "./bkgpdfparset.h"

namespace libTatami {

///
/// \brief The RbkgPdf class. Parameterization of background dt distribution
///
class RbkgPdf : public AbsPdf {
 public:
    ///
    /// \brief RbkgPdf
    /// \param fname
    ///
    explicit RbkgPdf(const str &fname);
    ///
    /// \brief RbkgPdf
    ///
    RbkgPdf(void);
    ///
    /// \brief operator ()
    /// \param x
    /// \return
    ///
    double operator()(const double& x);
    ///
    /// \brief operator ()
    /// \param evt
    /// \return
    ///
    double operator()(const ICPVEvt& evt);
    ///
    /// \brief Pars
    /// \return
    ///
    BkgPDFParSet& Pars(void) {return m_pars;}
    ///
    /// \brief SetSigma
    /// \param x
    ///
    void SetSigma(const double& x) {m_sigma = x;}
    ///
    /// \brief SetNDF
    /// \param x
    ///
    void SetNDF(const int& x) {m_ndf = x;}
    ///
    /// \brief SetScale
    /// \param x
    ///
    void SetScale(const double& x) {m_scale = x;}
    ///
    /// \brief SetShift
    /// \param x
    ///
    void SetShift(const double& x) {m_shift = x;}
    ///
    /// \brief SetScaleShift
    /// \param scale
    /// \param shift
    ///
    void SetScaleShift(const double& scale, const double& shift) {
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

 private:
    ///
    /// \brief Pdf
    /// \param x
    /// \param s
    /// \param ndf
    /// \return
    ///
    double Pdf(const double &x, const double &s, const int ndf);
    ///
    /// \brief DefaultInit
    /// \param mode
    /// \param svd
    /// \param mc
    /// \param tune
    ///
    void DefaultInit(const int mode, const int svd,
                     const bool mc, const bool tune);
    ///
    /// \brief AddOutlier
    /// \param x
    /// \param Lin
    /// \param nLi
    /// \return
    ///
    double AddOutlier(const double& x, const double Lin, const double& nLi);
    ///
    /// \brief cm2ps
    ///
    static const double cm2ps;
    ///
    /// \brief m_pars
    ///
    BkgPDFParSet m_pars;
    ///
    /// \brief m_scale
    ///
    double m_scale;
    ///
    /// \brief m_shift
    ///
    double m_shift;
    ///
    /// \brief m_sigma. Event-dependent estimation of uncertainty
    ///
    double m_sigma;
    ///
    /// \brief m_ndf. Number of degrees of freedom
    ///
    int m_ndf;
};

}  // namespace libTatami

#endif  // INCLUDE_RBKGPDF_H_
