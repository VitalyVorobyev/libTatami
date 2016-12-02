/** Copyright 2016 Vitaly Vorobyev
 * @file rascpdf.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#ifndef INCLUDE_RASCPDF_H_
#define INCLUDE_RASCPDF_H_

#include "./abspdf.h"
#include "./ResConst.h"
#include "./ResVar.h"
#include "./rascrnppars.h"
#include "./parmanager.h"

namespace libTatami {

///
/// \brief The RascPdf class describes asc vertex resolution
///
class RascPdf: public AbsPdf {
 public:
    ///
    /// \brief RascPdf
    /// \param dc
    ///
    explicit RascPdf(const DataClass &dc) :
        AbsPdf(), m_cnst(ParManager::SigParFile(dc)) {}
    ///
    /// \brief operator ()
    /// \param evt
    /// \return
    ///
    double operator()(const ICPVEvt& evt);
    ///
    /// \brief operator ()
    /// \param x
    /// \return
    ///
    double operator()(const double& x);
    ///
    /// \brief RascRnp
    /// \param evt
    /// \return
    ///
    double RascRnp(const ICPVEvt& evt);
    ///
    /// \brief norm_RascRnp
    /// \param evt
    /// \return
    ///
    double norm_RascRnp(const ICPVEvt& evt);

 private:
    ///
    /// \brief ReadVars
    /// \param evt
    /// \return
    ///
    int ReadVars(const ICPVEvt& evt);
    ///
    /// \brief PdfRascRnp
    /// \param evt
    /// \return
    ///
    double PdfRascRnp(const ICPVEvt& evt);
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
    RascRnpPars m_pars;
    ///
    /// \brief m_cnst
    ///
    ResConst    m_cnst;
    ///
    /// \brief m_vars
    ///
    RdetVar     m_vars;
    ///
    /// \brief keeptagl
    ///
    int keeptagl;
    ///
    /// \brief dz
    ///
    double dz;
};

}  // namespace libTatami

#endif  // INCLUDE_RASCPDF_H_
