/** Copyright 2016 Vitaly Vorobyev
 ** @file rrecpdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **
 **/

#ifndef INCLUDE_RRECPDF_H_
#define INCLUDE_RRECPDF_H_

#include "./abspdf.h"
#include "./ResConst.h"
#include "./ResVar.h"
#include "./rrecpars.h"
#include "./parmanager.h"

namespace libTatami {

///
/// \brief The RrecPdf class describes rec side vertex resolution.
///
class RrecPdf: public AbsPdf {
 public:
    ///
    /// \brief RrecPdf
    /// \param dc
    ///
    explicit RrecPdf(const DataClass &dc) :
        AbsPdf(), m_cnst(ParManager::SigParFile(dc)) {}
    ///
    /// \brief operator (). Normalized PDF for rec side resolution
    /// \param evt
    /// \return
    ///
    double operator()(const ICPVEvt& evt);
    ///
    /// \brief operator (). Normalized PDF for rec side resolution
    /// \param x
    /// \return
    ///
    double operator()(const double& x);
    ///
    /// \brief Rrec. Non-normalized PDF for rec side resolution
    /// \param evt
    /// \return
    ///
    double Rrec(const ICPVEvt& evt);
    ///
    /// \brief norm_Rrec. Normalization for rec side resolution
    /// \param evt
    /// \return
    ///
    double norm_Rrec(const ICPVEvt& evt);

 private:
    ///
    /// \brief Pdf
    /// \param evt
    /// \return
    ///
    double Pdf(const ICPVEvt& evt);
    ///
    /// \brief Pdf
    /// \return
    ///
    double Pdf(void) const;
    ///
    /// \brief Rrec. Non-normalized PDF for rec side resolution (const)
    /// \return
    ///
    double Rrec(void) const;
    ///
    /// \brief norm_Rrec. Normalization for rec side resolution (const)
    /// \return
    ///
    double norm_Rrec(void) const;
    ///
    /// \brief m_pars. Auxiliary parameters
    ///
    RrecPars m_pars;
    ///
    /// \brief m_cnst. Resolution constants
    ///
    ResConst m_cnst;
    ///
    /// \brief m_vars. Event-dependent variables
    ///
    RdetVar m_vars;
    ///
    /// \brief dz
    ///
    double dz;
};

}  // namespace libTatami

#endif  // INCLUDE_RRECPDF_H_
