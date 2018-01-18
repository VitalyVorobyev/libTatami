/** Copyright 2016 Vitaly Vorobyev
 * @file conv_coef.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#pragma once

#include <iostream>

#include "fenp.h"

namespace libTatami {

class RascRnpPars;
class RkPar;

///
/// \brief The conv_coef class.
/// \todo Simplify methodss arguments: use dedicated structures
///
class conv_coef {
    ///
    /// \brief add_EpEn_coef
    /// \param fEp1
    /// \param fEn2
    /// \param tau1
    /// \param tau2
    /// \param weight
    ///
    void add_EpEn_coef(double& fEp1, double& fEn2,
                       double tau1, double tau2, double weight = 1) const;
    ///
    /// \brief add_EnEp_coef
    /// \param fEn1
    /// \param fEp2
    /// \param tau1
    /// \param tau2
    /// \param weight
    ///
    void add_EnEp_coef(double& fEn1, double& fEp2,
                       double tau1, double tau2, double weight = 1) const;
    ///
    /// \brief add_EnEn_coef
    /// \param fEn1
    /// \param fEn2
    /// \param fxEn1
    /// \param tau1
    /// \param tau2
    /// \param weight
    ///
    void add_EnEn_coef(double& fEn1, double& fEn2, double& fxEn1,
                       double tau1, double tau2, double weight = 1) const;
    ///
    /// \brief add_EpEp_coef
    /// \param fEp1
    /// \param fEp2
    /// \param fxEp1
    /// \param tau1
    /// \param tau2
    /// \param weight
    ///
    void add_EpEp_coef(double& fEp1, double& fEp2, double& fxEp1,
                       double tau1, double tau2, double weight = 1) const;
    ///
    /// \brief add_ApEp_coef
    /// \param fAp1
    /// \param fMp1
    /// \param fEp2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_ApEp_coef(double& fAp1, double& fMp1, double& fEp2,
                       double tau1, double dm, double tau2,
                       double weight = 1.) const;
    ///
    /// \brief add_AnEn_coef
    /// \param fAn1
    /// \param fMn1
    /// \param fEn2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_AnEn_coef(double& fAn1, double& fMn1, double& fEn2,
                       double tau1, double dm,
                       double tau2, double weight = 1.) const;
    ///
    /// \brief add_ApEn_coef
    /// \param fAp1
    /// \param fMp1
    /// \param fEn2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_ApEn_coef(double& fAp1, double& fMp1, double& fEn2,
                       double tau1, double dm,
                       double tau2, double weight = 1.) const;
    ///
    /// \brief add_AnEp_coef
    /// \param fAn1
    /// \param fMn1
    /// \param fEp2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_AnEp_coef(double& fAn1, double& fMn1, double& fEp2,
                       double tau1, double dm,
                       double tau2, double weight = 1.) const;
    ///
    /// \brief add_MpEp_coef
    /// \param fMp1
    /// \param fAp1
    /// \param fEp2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_MpEp_coef(double& fMp1, double& fAp1, double& fEp2,
                       double tau1, double dm,
                       double tau2, double weight = 1.) const;
    ///
    /// \brief add_MnEn_coef
    /// \param fMn1
    /// \param fAn1
    /// \param fEn2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_MnEn_coef(double& fMn1, double& fAn1, double& fEn2,
                       double tau1, double dm,
                       double tau2, double weight = 1.) const;
    ///
    /// \brief add_MpEn_coef
    /// \param fMp1
    /// \param fAp1
    /// \param fEn2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_MpEn_coef(double& fMp1, double& fAp1, double& fEn2,
                       double tau1, double dm,
                       double tau2, double weight = 1.) const;
    ///
    /// \brief add_MnEp_coef
    /// \param fMn1
    /// \param fAn1
    /// \param fEp2
    /// \param tau1
    /// \param dm
    /// \param tau2
    /// \param weight
    ///
    void add_MnEp_coef(double& fMn1, double& fAn1, double& fEp2,
                       double tau1, double dm, double tau2,
                       double weight = 1.) const;
    ///
    /// \brief f
    ///
    mutable fEnp f;
    enum class State : int {A = 1, E = 2, M = 3, Null = 0};
    ///
    /// \brief state
    ///
    mutable State state;

 public:
    ///
    /// \brief conv_coef
    ///
    conv_coef(void) : state(State::Null) {}
    ///
    /// \brief Clear
    ///
    void Clear() const {state = State::Null;}
    ///
    /// \brief Ecoefs
    /// \param pRasc
    /// \param pRk
    /// \return
    ///
    const fEnp& Ecoefs(const RascRnpPars& pRasc, const RkPar& pRk) const;
    ///
    /// \brief Acoefs
    /// \param pRasc
    /// \param pRk
    /// \return
    ///
    const fEnp& Acoefs(const RascRnpPars& pRasc, const RkPar& pRk) const;
    ///
    /// \brief Mcoefs
    /// \param pRasc
    /// \param pRk
    /// \return
    ///
    const fEnp& Mcoefs(const RascRnpPars& pRasc, const RkPar& pRk) const;
};

}  // namespace libTatami
