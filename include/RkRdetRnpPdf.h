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

#include "absicpvpdf.h"
#include "wtag.h"
#include "ResVar.h"
#include "conv_coef.h"
#include "rascrnppars.h"
#include "rkparam.h"
#include "rrecpars.h"
#include "ResConst.h"

namespace libTatami {

// Forward declarations
class ICPVEvt;
class DataClass;

///
/// \brief The RkRdetRnpPdf class. Provides full parameterization of
/// signal dt distribution at Belle
/// \todo !!!! cexp and amix not initialized !!!!
///
class RkRdetRnpPdf : public AbsICPVPdf {
    ///
    /// \brief ReadVars
    /// \param evt
    /// \return
    ///
    int ReadVars(const ICPVEvt& evt) const;
    ///
    /// \brief ReadAndCalc
    /// \param evt
    /// \return
    ///
    int ReadAndCalc(const ICPVEvt& evt) const;
    ///
    /// \brief PdfAB
    /// \param evt
    /// \param otlr
    /// \param no_interf
    /// \return
    ///
    double PdfAB(const ICPVEvt& evt, bool otlr = true,
                 bool no_interf = false) const;
    ///
    /// \brief PdfAB
    /// \param otlr
    /// \param no_interf
    /// \return
    ///
    double PdfAB(bool otlr = true, bool no_interf = false) const;
    ///
    /// \brief EfRkRdetRnp_fullrec
    /// \return
    ///
    double EfRkRdetRnp_fullrec(void) const;
    ///
    /// \brief MfRkRdetRnp_fullrec
    /// \return
    ///
    double MfRkRdetRnp_fullrec(void) const;
    ///
    /// \brief AfRkRdetRnp_fullrec
    /// \return
    ///
    double AfRkRdetRnp_fullrec(void) const;
    ///
    /// \brief norm_EfRkRdetRnp_fullrec
    /// \return
    ///
    double norm_EfRkRdetRnp_fullrec(void) const;
    ///
    /// \brief norm_MfRkRdetRnp_fullrec
    /// \return
    ///
    double norm_MfRkRdetRnp_fullrec(void) const;
    ///
    /// \brief norm_AfRkRdetRnp_fullrec
    /// \return
    ///
    double norm_AfRkRdetRnp_fullrec(void) const;
    ///
    /// \brief EfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double EfRkRdetRnp_full_sup(double mu, double sigma) const;
    ///
    /// \brief MfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double MfRkRdetRnp_full_sup(double mu, double sigma) const;
    ///
    /// \brief AfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double AfRkRdetRnp_full_sup(double mu, double sigma) const;
    ///
    /// \brief norm_EfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double norm_EfRkRdetRnp_full_sup(double mu, double sigma) const;
    ///
    /// \brief norm_MfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double norm_MfRkRdetRnp_full_sup(double mu, double sigma) const;
    ///
    /// \brief norm_AfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double norm_AfRkRdetRnp_full_sup(double mu, double sigma) const;
    ///
    /// \brief Xsum4
    /// \param Xmm
    /// \param Xmt
    /// \param Xtm
    /// \param Xtt
    /// \return
    ///
    double Xsum4(double Xmm, double Xmt, double Xtm, double Xtt) const;
    ///
    /// \brief Xsum2mmmt
    /// \param Xmm
    /// \param Xmt
    /// \return
    ///
    double Xsum2mmmt(double Xmm, double Xmt) const;
    ///
    /// \brief Xsum2mmtm
    /// \param Xmm
    /// \param Xtm
    /// \return
    ///
    double Xsum2mmtm(double Xmm, double Xtm) const;
    ///
    /// \brief flavor
    ///
//    int flavor;
    ///
    /// \brief keeptagl
    ///
    mutable int keeptagl;
    ///
    /// \brief m_svd
    ///
//    int m_svd;
    ///
    /// \brief cexp
    ///
    mutable double cexp;
    ///
    /// \brief amix
    ///
    double amix;
    ///
    /// \brief dt
    ///
    mutable double dt;
    ///
    /// \brief m_rvar
    ///
    mutable RdetVar m_rvar;
    ///
    /// \brief m_avar
    ///
    mutable RdetVar m_avar;
    ///
    /// \brief costhBcms
    ///
    mutable double costhBcms;
    ///
    /// \brief m_apar
    ///
    mutable RascRnpPars m_apar;
    ///
    /// \brief m_rpar
    ///
    mutable RrecPars m_rpar;
    ///
    /// \brief m_kpar
    ///
    mutable RkPar m_kpar;
    ///
    /// \brief m_wtag
    ///
    const WTag m_wtag;
    ///
    /// \brief m_cnst
    ///
    const ResConst m_cnst;
    ///
    /// \brief coco
    ///
    const conv_coef coco;
    ///
    /// \brief mu_mm
    ///
    mutable double mu_mm;
    ///
    /// \brief mu_mt
    ///
    mutable double mu_mt;
    ///
    /// \brief mu_tm
    ///
    mutable double mu_tm;
    ///
    /// \brief mu_tt
    ///
    mutable double mu_tt;
    ///
    /// \brief s_mm
    ///
    mutable double s_mm;
    ///
    /// \brief s_mt
    ///
    mutable double s_mt;
    ///
    /// \brief s_tm
    ///
    mutable double s_tm;
    ///
    /// \brief s_tt
    ///
    mutable double s_tt;

 public:
    explicit RkRdetRnpPdf(const DataClass &dc);
    ///
    /// \brief operator (). Calculate PDF
    /// \param ext
    /// \return
    ///
    double operator()(const ICPVEvt& ext) const override final;
    ///
    /// \brief operator (). Calculate PDF
    /// \param x
    /// \return
    ///
    double operator()(double x) const override final;
    ///
    /// \brief EfRkRdetRnp_fullrec. Unnormalized convolution Exp x Gauss
    /// \param evt
    /// \return
    ///
    double EfRkRdetRnp_fullrec(const ICPVEvt& evt) const;
    ///
    /// \brief MfRkRdetRnp_fullrec. Unnormalized convolution Exp*cos x Gauss
    /// \param evt
    /// \return
    ///
    double MfRkRdetRnp_fullrec(const ICPVEvt& evt) const;
    ///
    /// \brief AfRkRdetRnp_fullrec. Unnormalized convolution Exp*sin x Gauss
    /// \param evt
    /// \return
    ///
    double AfRkRdetRnp_fullrec(const ICPVEvt& evt) const;
    ///
    /// \brief norm_EfRkRdetRnp_fullrec. Normalization for
    /// convolution Exp x Gauss
    /// \param evt
    /// \return
    ///
    double norm_EfRkRdetRnp_fullrec(const ICPVEvt& evt) const;
    ///
    /// \brief norm_MfRkRdetRnp_fullrec. Normalization for
    /// convolution Exp*cos x Gauss
    /// \param evt
    /// \return
    ///
    double norm_MfRkRdetRnp_fullrec(const ICPVEvt& evt) const;
    ///
    /// \brief norm_AfRkRdetRnp_fullrec. Normalization for
    /// convolution Exp*sin x Gauss
    /// \param evt
    /// \return
    ///
    double norm_AfRkRdetRnp_fullrec(const ICPVEvt& evt) const;
};

}  // namespace libTatami
