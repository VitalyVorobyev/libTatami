/** Copyright 2016 Vitaly Vorobyev
 ** @file rk.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 */

#ifndef INCLUDE_RK_H_
#define INCLUDE_RK_H_

#include "./absicpvpdf.h"
#include "./rkparam.h"

namespace libTatami {

///
/// \brief The RkPdf class. Describes time resolution part related to
/// kinematic approximation dz -> dt.
///
class RkPdf: public AbsICPVPdf {
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
    double operator()(const ICPVEvt& evt);
    ///
    /// \brief operator (). Calculate PDF
    /// \param x
    /// \return
    ///
    double operator()(const double& x) const;

 private:
    ///
    /// \brief Pdf
    /// \param x
    /// \return
    ///
    double Pdf(const double& x) const;
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
    double EfRk(const double& x) const;
    ///
    /// \brief AfRk
    /// \param x
    /// \return
    ///
    double AfRk(const double& x) const;
    ///
    /// \brief MfRk
    /// \param x
    /// \return
    ///
    double MfRk(const double& x) const;
    ///
    /// \brief m_pars
    ///
    RkPar m_pars;
};

}  // namespace libTatami

#endif  // INCLUDE_RK_H_
