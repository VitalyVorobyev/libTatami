#ifndef RRECPDF_H
#define RRECPDF_H

#include "abspdf.h"
#include "ResConst.h"
#include "ResVar.h"
#include "rrecpars.h"
#include "parmanager.h"

///
/// \brief The RrecPdf class describes rec side vertex resolution.
///

class RrecPdf: public AbsPdf{
public:
  RrecPdf(const DataClass &dc): AbsPdf(), m_cnst(ParManager::SigParFile(dc)) {}
/// Normalized PDF for rec side resolution
  double operator()(const ICPVEvt& evt);
/// Normalized PDF for rec side resolution
  double operator()(const double& x);
/// Non-normalized PDF for rec side resolution
  double      Rrec(const ICPVEvt& evt);
/// Normalization for rec side resolution
  double norm_Rrec(const ICPVEvt& evt);

private:
  double Pdf(const ICPVEvt& evt);
  double Pdf(void) const;

/// Non-normalized PDF for rec side resolution (const)
  double      Rrec(void) const;
/// Normalization for rec side resolution (const)
  double norm_Rrec(void) const;
//  double AddOutlier(const double& x, const double& pdf, const double& norm);

  /// Auxiliary parameters
  RrecPars m_pars;
  /// Resolution constants
  ResConst m_cnst;
  /// Event-dependent variables
  RdetVar m_vars;
  double dz;
};

#endif // RRECPDF_H
