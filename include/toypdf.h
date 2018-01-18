/** Copyright 2016-2018 Vitaly Vorobyev
 ** @file toypdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **
 **/

#pragma once

#include "./absicpvpdf.h"

#include <cstdint>

namespace libTatami {

///
/// \brief The ToyPdf class. Simple CP-violating PDF with Gaussian resolution
/// function. Intended for toy studies.
///
class ToyPdf : public AbsICPVPdf {
    ///
    /// \brief pdfSig
    /// \param dt
    /// \param wid
    /// \return
    ///
    double pdfSig(double dt, double wid) const;
    ///
    /// \brief pdfBkg
    /// \param dt
    /// \param wid
    /// \return
    ///
    double pdfBkg(double dt, double wid) const;
    ///
    /// \brief m_m. Mean of resolution Gaussian
    ///
    double m_m;
    ///
    /// \brief m_w. Width of resolution Gaussian
    ///
    double m_w;
    ///
    /// \brief m_fbkg. Backgroud fraction
    ///
    double m_fbkg;

 public:
    ///
    /// \brief ToyPdf
    /// \param m
    /// \param w
    /// \param fb
    ///
    ToyPdf(double m, double w, double fb=0, double wtag=0,
           double ll=-10, double ul=10);
    ///
    /// \brief ToyPdf
    ///
    ToyPdf(void): ToyPdf(0., 1.0, 0.) {}
    ///
    /// \brief operator ()
    /// \param evt
    /// \return
    ///
    double operator() (const ICPVEvt& evt) const override final;
    /**
     * @brief operator ()
     * @param dt. Time difference between signal and tagging B0 meson decays
     * @param tag. B0 tag
     * @return PDF value
     */
    double operator() (double dt, int32_t tag);
    /**
     * @brief operator ()
     * @param dt. Time difference between signal and tagging B0 meson decays
     * @param tag. B0 tag
     * @param c. Coeff near cos
     * @param s. Coeff near sin
     * @return PDF value
     */
    double operator() (double dt, int32_t tag, double c, double s) const;
    ///
    /// \brief operator () Calculate PDF
    /// \param dt
    /// \return
    ///
    double operator() (double dt) const override final;
    ///
    /// \brief operator (). Calculate PDF
    /// \param dt
    /// \param fbkg
    /// \param scale
    /// \return
    ///
    double operator() (double dt, double fbkg, double scale = 1) const;
    ///
    /// \brief operator (). Calculate PDF
    /// \param dt
    /// \param c
    /// \param s
    /// \param fbkg
    /// \param scale
    /// \return
    ///
    double operator() (double dt, double c, double s,
                       double fbkg, double scale = 1) const;
    ///
    /// \brief SetResMean. Set mean of resolution Gaussian
    /// \param v
    ///
    void SetResMean(double v)  {m_m = v;}
    ///
    /// \brief SetResWidth. Set width of resolution Gaussian
    /// \param v
    ///
    void SetResWidth(double v) {m_w = v;}
    ///
    /// \brief SetFbkg. Set backgroung fraction
    /// \param v
    ///
    void SetFbkg(double v) {m_fbkg = v;}
    ///
    /// \brief ResMean. Mean of resolution Gaussian
    /// \return
    ///
    double ResMean(void)  const {return m_m;}
    ///
    /// \brief ResWidth. Width of resolution Gaussian
    /// \return
    ///
    double ResWidth(void) const {return m_w;}
    ///
    /// \brief Fbkg. Backgroung fraction
    /// \return
    ///
    double Fbkg(void) const {return m_fbkg;}
    ///
    /// \brief print_params
    ///
    void print_params(void) const;
};

}  // namespace libTatami
