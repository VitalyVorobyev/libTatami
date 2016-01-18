#ifndef TOYPDF_H
#define TOYPDF_H

#include "absicpvpdf.h"

///
/// \brief The ToyPdf class. Simple CP-violating PDF with Gaussian resolution function. Intended for toy studies.
///

class ToyPdf : public AbsICPVPdf{
public:
  ToyPdf(): AbsICPVPdf(), m_m(0.), m_w(1.), m_fbkg(0.), m_wrtag(0.) {}
  ToyPdf(const ToyPdf& opdf);

  /// Calculate PDF
  double operator() (const double& dt);
  /// Calculate PDF
  double operator() (const double& dt, const double& fbkg, const double& scale = 1);
  /// Calculate PDF
  double operator() (const double& dt, const double& c, const double& s, const double& fbkg, const double& scale = 1);

  /// Set mean of resolution Gaussian
  void SetResMean(const double& v)  {m_m = v; return;}
  /// Set width of resolution Gaussian
  void SetResWidth(const double& v) {m_w = v; return;}
  /// Set backgroung fraction
  void SetFbkg(const double& v)     {m_fbkg = v; return;}
  /// Set wrong tagging probability
  void SetWrTag(const double& v)    {m_wrtag = v; return;}

  /// Mean of resolution Gaussian
  double ResMean(void)  const {return m_m;}
  /// Width of resolution Gaussian
  double ResWidth(void) const {return m_w;}
  /// Backgroung fraction
  double Fbkg(void)     const {return m_fbkg;}
  /// Wrong tagging probability
  double WrTag(void)    const {return m_wrtag;}

private:
  double pdfSig(const double& dt, const double& wid);
  double pdfBkg(const double& dt, const double& wid);

  /// Mean of resolution Gaussian
  double m_m;
  /// Width of resolution Gaussian
  double m_w;
  /// Backgroud fraction
  double m_fbkg;
  /// Wrong tag probability
  double m_wrtag;
};

#endif // TOYPDF_H
