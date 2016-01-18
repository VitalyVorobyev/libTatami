#ifndef ABSICPVPDF_H
#define ABSICPVPDF_H

#include "abspdf.h"

///
/// \brief The AbsICPVPdf class. Abstract class for PDF describing CP-violating dt distribution
///

class AbsICPVPdf : public AbsPdf{
public:
  AbsICPVPdf() {}

  /// Set coefficient near cos(dt)
  void SetC(const double& v) {m_c = v; return;}
  /// Set coefficient near sin(dt)
  void SetS(const double& v) {m_s = v; return;}
  /// Set lifetime and mass difference
  void SetTauDm(const double& tau, const double& dm) {SetTau(tau); SetDm(dm);}
  /// Set lifetime
  void SetTau(const double& tau) {m_tau = tau;}
  /// Set mass difference
  void SetDm(const double& dm) {m_dm = dm;}

  /// Coefficient near cos(dt)
  double C(void) const {return m_c;}
  /// Coefficient near sin(dt)
  double S(void) const {return m_s;}
  /// Lifetime
  double tau(void) const {return m_tau;}
  /// Mass difference
  double dm(void) const {return m_dm;}

protected:
  /// Coefficient near cos(dt)
  double m_c;
  /// Coefficient near sin(dt)
  double m_s;
  /// Lifetime
  double m_tau;
  /// Mass difference
  double m_dm;
};

#endif // ABSICPVPDF_H
