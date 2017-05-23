/** Copyright 2016 Vitaly Vorobyev
 ** @file toypdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **
 **/

#ifndef INCLUDE_TOYPDF_H_
#define INCLUDE_TOYPDF_H_

#include "./absicpvpdf.h"

namespace libTatami {

///
/// \brief The ToyPdf class. Simple CP-violating PDF with Gaussian resolution
/// function. Intended for toy studies.
///
class ToyPdf : public AbsICPVPdf {
 public:
    ///
    /// \brief ToyPdf
    /// \param m
    /// \param w
    /// \param fb
    ///
    ToyPdf(const double& m, const double& w, const double& fb=0,
           const double &wtag=0, const double& ll=-10, const double& ul=10);
    ///
    /// \brief ToyPdf
    ///
    ToyPdf(void): ToyPdf(0., 1.0, 0.) {}
    ///
    /// \brief ToyPdf
    /// \param opdf
    ///
    ToyPdf(const ToyPdf& opdf);
    ~ToyPdf() {}
    ///
    /// \brief operator ()
    /// \param evt
    /// \return
    ///
    double operator() (const ICPVEvt& evt);
    /**
     * @brief operator ()
     * @param dt. Time difference between signal and tagging B0 meson decays
     * @param tag. B0 tag
     * @return PDF value
     */
    double operator() (const double& dt, const int tag);
    /**
     * @brief operator ()
     * @param dt. Time difference between signal and tagging B0 meson decays
     * @param tag. B0 tag
     * @param c. Coeff near cos
     * @param s. Coeff near sin
     * @return PDF value
     */
    double operator() (const double& dt, const int tag,
                       const double& c, const double& s);
    ///
    /// \brief operator () Calculate PDF
    /// \param dt
    /// \return
    ///
    double operator() (const double& dt);
    ///
    /// \brief operator (). Calculate PDF
    /// \param dt
    /// \param fbkg
    /// \param scale
    /// \return
    ///
    double operator() (const double& dt, const double& fbkg,
                       const double& scale = 1);
    ///
    /// \brief operator (). Calculate PDF
    /// \param dt
    /// \param c
    /// \param s
    /// \param fbkg
    /// \param scale
    /// \return
    ///
    double operator() (const double& dt, const double& c, const double& s,
                       const double& fbkg, const double& scale = 1);
    ///
    /// \brief SetResMean. Set mean of resolution Gaussian
    /// \param v
    ///
    void SetResMean(const double& v)  {m_m = v; return;}
    ///
    /// \brief SetResWidth. Set width of resolution Gaussian
    /// \param v
    ///
    void SetResWidth(const double& v) {m_w = v; return;}
    ///
    /// \brief SetFbkg. Set backgroung fraction
    /// \param v
    ///
    void SetFbkg(const double& v)     {m_fbkg = v; return;}
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

 private:
    ///
    /// \brief pdfSig
    /// \param dt
    /// \param wid
    /// \return
    ///
    double pdfSig(const double& dt, const double& wid);
    ///
    /// \brief pdfBkg
    /// \param dt
    /// \param wid
    /// \return
    ///
    double pdfBkg(const double& dt, const double& wid);
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
};

}  // namespace libTatami

#endif  // INCLUDE_TOYPDF_H_
