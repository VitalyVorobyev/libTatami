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

#ifndef INCLUDE_RRECPARS_H_
#define INCLUDE_RRECPARS_H_

#include "./ResVar.h"
#include "./ResConst.h"

namespace libTatami {

///
/// \brief The RrecPars class calculates and keeps auxiliary parameters
/// for rec side vertex resolution
///
class RrecPars {
 public:
    ///
    /// \brief RrecPars
    ///
    RrecPars(void);
    ///
    /// \brief Calculate
    /// \param cnst
    /// \param var
    ///
    void Calculate(const ResConst& cnst, const RdetVar& var) {
        calc_vtxparam_rec(var); Rrec_param(cnst, var);
    }
    ///
    /// \brief Smain_rec
    /// \return
    ///
    double Smain_rec()   const {return m_Smain_rec;}
    ///
    /// \brief Stail_rec
    /// \return
    ///
    double Stail_rec()   const {return m_Stail_rec;}
    ///
    /// \brief ftail_rec
    /// \return
    ///
    double ftail_rec()   const {return m_ftail_rec;}
    ///
    /// \brief mu_main_rec
    /// \return
    ///
    double mu_main_rec() const {return m_mu_main_rec;}
    ///
    /// \brief mu_tail_rec
    /// \return
    ///
    double mu_tail_rec() const {return m_mu_tail_rec;}

 private:
    ///
    /// \brief calc_vtxparam_rec
    /// \param var
    ///
    void calc_vtxparam_rec(const RdetVar& var);
    ///
    /// \brief Rrec_param
    /// \param cnst
    /// \param var
    ///
    void Rrec_param(const ResConst& cnst, const RdetVar& var);
    ///
    /// \brief m_Smain_rec
    ///
    double m_Smain_rec;
    ///
    /// \brief m_Stail_rec
    ///
    double m_Stail_rec;
    ///
    /// \brief m_ftail_rec
    ///
    double m_ftail_rec;
    ///
    /// \brief m_mu_main_rec
    ///
    double m_mu_main_rec;
    ///
    /// \brief m_mu_tail_rec
    ///
    double m_mu_tail_rec;
    ///
    /// \brief m_xi_rec
    ///
    double m_xi_rec;
    ///
    /// \brief m_st_rec
    ///
    double m_st_rec;
};

}  // namespace libTatami

#endif  // INCLUDE_RRECPARS_H_
