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

#pragma once

#include "abspdf.h"
#include "ResConst.h"
#include "rrecpars.h"
#include "ResVar.h"

namespace libTatami {

class ICPVEvt;
class DataClass;

///
/// \brief The RrecPdf class describes rec side vertex resolution.
///
class RrecPdf: public AbsPdf {
    ///
    /// \brief Pdf
    /// \param evt
    /// \return
    ///
    double Pdf(const ICPVEvt& evt) const;
    ///ll
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
    mutable RrecPars m_pars;
    ///
    /// \brief m_cnst. Resolution constants
    ///
    ResConst m_cnst;
    ///
    /// \brief m_vars. Event-dependent variables
    ///
    mutable RdetVar m_vars;
    ///
    /// \brief dz
    ///
    mutable double dz;

 public:
    ///
    /// \brief RrecPdf
    /// \param dc
    ///
    explicit RrecPdf(const DataClass &dc);
    ///
    /// \brief operator (). Normalized PDF for rec side resolution
    /// \param evt
    /// \return
    ///
    double operator()(const ICPVEvt& evt) const override final;
    ///
    /// \brief operator (). Normalized PDF for rec side resolution
    /// \param x
    /// \return
    ///
    double operator()(double x) const override final;
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
};

}  // namespace libTatami
