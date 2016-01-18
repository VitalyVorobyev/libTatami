#ifndef RASCPDF_H
#define RASCPDF_H

#include "abspdf.h"
#include "ResConst.h"
#include "ResVar.h"
#include "rascrnppars.h"
#include "parmanager.h"

///
/// \brief The RascPdf class describes asc vertex resolution
///

class RascPdf: public AbsPdf{
public:
  RascPdf(const DataClass &dc): AbsPdf(), m_cnst(ParManager::SigParFile(dc)) {}
  double operator()(const ICPVEvt& evt);
  double operator()(const double& x);

  double      RascRnp(const ICPVEvt& evt);
  double norm_RascRnp(const ICPVEvt& evt);

private:
  int ReadVars(const ICPVEvt& evt);
  double PdfRascRnp(const ICPVEvt& evt);
  double PdfRascRnp(void) const;
  double      RascRnp(void) const;
  double norm_RascRnp(void) const;

  RascRnpPars m_pars;
  ResConst    m_cnst;
  RdetVar     m_vars;
  int keeptagl;
  double dz;
};

#endif // RASCPDF_H

