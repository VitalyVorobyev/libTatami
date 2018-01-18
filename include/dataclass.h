/** Copyright 2016 Vitaly Vorobyev
 ** @file dataclass.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#pragma once

#include <cstdint>

namespace libTatami {

///
/// \brief The DataClass class
///
class DataClass {
    ///
    /// \brief m_svd
    ///
    uint16_t m_svd;
    ///
    /// \brief m_mc
    ///
    bool m_mc;
    ///
    /// \brief m_bp
    ///
    bool m_bp;

 public:
    ///
    /// \brief DataClass
    /// \param svd
    /// \param mc
    /// \param bp
    ///
    DataClass(uint16_t svd, bool mc, bool bp);
    ///
    /// \brief DataClass
    /// \param svd
    /// \param mc
    ///
    DataClass(uint16_t svd, bool mc);
    ///
    /// \brief DataClass
    /// \param svd
    ///
    explicit DataClass(uint16_t svd);
    ///
    /// \brief DataClass
    ///
    DataClass(void);
    ///
    /// \brief SetSVD
    /// \param v
    ///
    void SetSVD(uint16_t v) {m_svd = v;}
    ///
    /// \brief SetMC
    /// \param v
    ///
    void SetMC(bool v) {m_mc  = v;}
    ///
    /// \brief SetBp
    /// \param v
    ///
    void SetBp(bool v) {m_bp  = v;}
    ///
    /// \brief Svd
    /// \return
    ///
    int Svd(void) const {return m_svd;}
    ///
    /// \brief MC
    /// \return
    ///
    bool MC(void) const {return m_mc;}
    ///
    /// \brief Bp
    /// \return
    ///
    bool Bp(void) const {return m_bp;}
};

}  // namespace libTatami

