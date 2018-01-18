/** Copyright 2016 Vitaly Vorobyev
 * @file rascpdf.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#pragma once

#include "abspdf.h"
#include "ResConst.h"
#include "ResVar.h"
#include "rascrnppars.h"

namespace libTatami {

class DataClass;

///
/// \brief The RascPdf class describes asc vertex resolution
///
class RascPdf: public AbsPdf {
    ///
    /// \brief ReadVars
    /// \param evt
    /// \return
    ///
    int ReadVars(const ICPVEvt& evt) const;
    ///
    /// \brief PdfRascRnp
    /// \param evt
    /// \return
    ///
    double PdfRascRnp(const ICPVEvt& evt) const;
    ///
    /// \brief PdfRascRnp
    /// \return
    ///
    double PdfRascRnp(void) const;
    ///
    /// \brief RascRnp
    /// \return
    ///
    double RascRnp(void) const;
    ///
    /// \brief norm_RascRnp
    /// \return
    ///
    double norm_RascRnp(void) const;
    ///
    /// \brief m_pars
    ///
    mutable RascRnpPars m_pars;
    ///
    /// \brief m_cnst
    ///
    const ResConst m_cnst;
    ///
    /// \brief m_vars
    ///
    mutable RdetVar m_vars;
    ///
    /// \brief keeptagl
    ///
    mutable int keeptagl;
    ///
    /// \brief dz
    ///
    mutable double dz;

 public:
    ///
    /// \brief RascPdf
    /// \param dc
    ///
    explicit RascPdf(const DataClass &dc);
    ///
    /// \brief operator ()
    /// \param evt
    /// \return
    ///
    double operator()(const ICPVEvt& evt) const override final;
    ///
    /// \brief operator ()
    /// \param x
    /// \return
    ///
    double operator()(double x) const override final;
    ///
    /// \brief RascRnp
    /// \param evt
    /// \return
    ///
    double RascRnp(const ICPVEvt& evt) const;
    ///
    /// \brief norm_RascRnp
    /// \param evt
    /// \return
    ///
    double norm_RascRnp(const ICPVEvt& evt) const;
};

}  // namespace libTatami
