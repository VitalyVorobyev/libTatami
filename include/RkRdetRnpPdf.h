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

#ifndef INCLUDE_RKRDETRNPPDF_H_
#define INCLUDE_RKRDETRNPPDF_H_

#include "./absicpvpdf.h"
#include "./conv_coef.h"
#include "./icpvevent.h"
#include "./parmanager.h"
#include "./rascrnppars.h"
#include "./rkparam.h"
#include "./rrecpars.h"
#include "./ResConst.h"
#include "./ResVar.h"
#include "./ttools.h"
#include "./wtag.h"

namespace libTatami {

///
/// \brief The RkRdetRnpPdf class. Provides full parameterization of
/// signal dt distribution at Belle
///
class RkRdetRnpPdf : public AbsICPVPdf {
 public:
    explicit RkRdetRnpPdf(const DataClass &dc);
    ///
    /// \brief operator (). Calculate PDF
    /// \param ext
    /// \return
    ///
    double operator()(const ICPVEvt& ext);
    ///
    /// \brief operator (). Calculate PDF
    /// \param x
    /// \return
    ///
    double operator()(const double& x);
    ///
    /// \brief EfRkRdetRnp_fullrec. Unnormalized convolution Exp x Gauss
    /// \param evt
    /// \return
    ///
    double EfRkRdetRnp_fullrec(const ICPVEvt& evt);
    ///
    /// \brief MfRkRdetRnp_fullrec. Unnormalized convolution Exp*cos x Gauss
    /// \param evt
    /// \return
    ///
    double MfRkRdetRnp_fullrec(const ICPVEvt& evt);
    ///
    /// \brief AfRkRdetRnp_fullrec. Unnormalized convolution Exp*sin x Gauss
    /// \param evt
    /// \return
    ///
    double AfRkRdetRnp_fullrec(const ICPVEvt& evt);
    ///
    /// \brief norm_EfRkRdetRnp_fullrec. Normalization for
    /// convolution Exp x Gauss
    /// \param evt
    /// \return
    ///
    double norm_EfRkRdetRnp_fullrec(const ICPVEvt& evt);
    ///
    /// \brief norm_MfRkRdetRnp_fullrec. Normalization for
    /// convolution Exp*cos x Gauss
    /// \param evt
    /// \return
    ///
    double norm_MfRkRdetRnp_fullrec(const ICPVEvt& evt);
    ///
    /// \brief norm_AfRkRdetRnp_fullrec. Normalization for
    /// convolution Exp*sin x Gauss
    /// \param evt
    /// \return
    ///
    double norm_AfRkRdetRnp_fullrec(const ICPVEvt& evt);

 private:
    ///
    /// \brief ReadVars
    /// \param evt
    /// \return
    ///
    int ReadVars(const ICPVEvt& evt);
    ///
    /// \brief ReadAndCalc
    /// \param evt
    /// \return
    ///
    int ReadAndCalc(const ICPVEvt& evt);
    ///
    /// \brief PdfAB
    /// \param evt
    /// \param otlr
    /// \param no_interf
    /// \return
    ///
    double PdfAB(const ICPVEvt& evt, const bool otlr = true,
                 const bool no_interf = false);
    ///
    /// \brief PdfAB
    /// \param otlr
    /// \param no_interf
    /// \return
    ///
    double PdfAB(const bool otlr = true, const bool no_interf = false);
    ///
    /// \brief EfRkRdetRnp_fullrec
    /// \return
    ///
    double EfRkRdetRnp_fullrec(void);
    ///
    /// \brief MfRkRdetRnp_fullrec
    /// \return
    ///
    double MfRkRdetRnp_fullrec(void);
    ///
    /// \brief AfRkRdetRnp_fullrec
    /// \return
    ///
    double AfRkRdetRnp_fullrec(void);
    ///
    /// \brief norm_EfRkRdetRnp_fullrec
    /// \return
    ///
    double norm_EfRkRdetRnp_fullrec(void);
    ///
    /// \brief norm_MfRkRdetRnp_fullrec
    /// \return
    ///
    double norm_MfRkRdetRnp_fullrec(void);
    ///
    /// \brief norm_AfRkRdetRnp_fullrec
    /// \return
    ///
    double norm_AfRkRdetRnp_fullrec(void);
    ///
    /// \brief EfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double EfRkRdetRnp_full_sup(const double& mu, const double& sigma);
    ///
    /// \brief MfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double MfRkRdetRnp_full_sup(const double& mu, const double& sigma);
    ///
    /// \brief AfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double AfRkRdetRnp_full_sup(const double& mu, const double& sigma);
    ///
    /// \brief norm_EfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double norm_EfRkRdetRnp_full_sup(const double& mu, const double& sigma);
    ///
    /// \brief norm_MfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double norm_MfRkRdetRnp_full_sup(const double& mu, const double& sigma);
    ///
    /// \brief norm_AfRkRdetRnp_full_sup
    /// \param mu
    /// \param sigma
    /// \return
    ///
    double norm_AfRkRdetRnp_full_sup(const double& mu, const double& sigma);
//    double AddOutlier( const double& x, const double& Lin,
//    const double& nLi = 1.0, const double& alpha = 1.0);
//    double Add_Outlier(const double& x, const double& Lin,
//    const double& nLi, const double& alpha);
    ///
    /// \brief Xsum4
    /// \param Xmm
    /// \param Xmt
    /// \param Xtm
    /// \param Xtt
    /// \return
    ///
    double Xsum4(const double& Xmm, const double& Xmt,
                 const double& Xtm, const double& Xtt) const;
    ///
    /// \brief Xsum2mmmt
    /// \param Xmm
    /// \param Xmt
    /// \return
    ///
    double Xsum2mmmt(const double& Xmm, const double& Xmt) const;
    ///
    /// \brief Xsum2mmtm
    /// \param Xmm
    /// \param Xtm
    /// \return
    ///
    double Xsum2mmtm(const double& Xmm, const double& Xtm) const;
    ///
    /// \brief flavor
    ///
    int flavor;
    ///
    /// \brief keeptagl
    ///
    int keeptagl;
    ///
    /// \brief m_svd
    ///
    int m_svd;
    ///
    /// \brief cexp
    ///
    double cexp;
    ///
    /// \brief amix
    ///
    double amix;
    ///
    /// \brief m_wtag
    ///
    WTag m_wtag;
    ///
    /// \brief dt
    ///
    double dt;
    ///
    /// \brief m_rvar
    ///
    RdetVar m_rvar;
    ///
    /// \brief m_avar
    ///
    RdetVar m_avar;
    ///
    /// \brief costhBcms
    ///
    double costhBcms;
    ///
    /// \brief m_apar
    ///
    RascRnpPars m_apar;
    ///
    /// \brief m_rpar
    ///
    RrecPars m_rpar;
    ///
    /// \brief m_kpar
    ///
    RkPar m_kpar;
    ///
    /// \brief m_cnst
    ///
    ResConst m_cnst;
    ///
    /// \brief coco
    ///
    conv_coef coco;
    ///
    /// \brief mu_mm
    ///
    double mu_mm;
    ///
    /// \brief mu_mt
    ///
    double mu_mt;
    ///
    /// \brief mu_tm
    ///
    double mu_tm;
    ///
    /// \brief mu_tt
    ///
    double mu_tt;
    ///
    /// \brief s_mm
    ///
    double s_mm;
    ///
    /// \brief s_mt
    ///
    double s_mt;
    ///
    /// \brief s_tm
    ///
    double s_tm;
    ///
    /// \brief s_tt
    ///
    double s_tt;
};

}  // namespace libTatami

#endif  // INCLUDE_RKRDETRNPPDF_H_
