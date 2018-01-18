/** Copyright 2016 Vitaly Vorobyev
 ** @file rk.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 */

#pragma once

#include "absicpvpdf.h"
#include "rkparam.h"

namespace libTatami {

///
/// \brief The RkPdf class. Describes time resolution part related to
/// kinematic approximation dz -> dt.
///
class RkPdf: public AbsICPVPdf {
    ///
    /// \brief Pdf
    /// \param x
    /// \return
    ///
    double Pdf(double x) const;
    ///
    /// \brief norm_EfRk
    /// \return
    ///
    double norm_EfRk(void) const;
    ///
    /// \brief norm_AfRk
    /// \return
    ///
    double norm_AfRk(void) const;
    ///
    /// \brief norm_MfRk
    /// \return
    ///
    double norm_MfRk(void) const;
    ///
    /// \brief EfRk
    /// \param x
    /// \return
    ///
    double EfRk(double x) const;
    ///
    /// \brief AfRk
    /// \param x
    /// \return
    ///
    double AfRk(double x) const;
    ///
    /// \brief MfRk
    /// \param x
    /// \return
    ///
    double MfRk(double x) const;
    ///
    /// \brief m_pars
    ///
    mutable RkPar m_pars;

 public:
    ///
    /// \brief RkPdf
    ///
    RkPdf(): AbsICPVPdf() {}
    ///
    /// \brief operator (). Calculate PDF
    /// \param evt
    /// \return
    ///
    double operator()(const ICPVEvt& evt) const override final;
    ///
    /// \brief operator (). Calculate PDF
    /// \param x
    /// \return
    ///
    double operator()(double x) const override final;
};

}  // namespace libTatami
