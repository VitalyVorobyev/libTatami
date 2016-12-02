/** Copyright 2016 Vitaly Vorobyev
 ** @file dataclass.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_DATACLASS_H_
#define INCLUDE_DATACLASS_H_

namespace libTatami {

///
/// \brief The DataClass class
///
class DataClass {
 public:
    ///
    /// \brief DataClass
    /// \param svd
    /// \param mc
    /// \param bp
    ///
    DataClass(const int svd, const bool mc, const bool bp);
    ///
    /// \brief DataClass
    /// \param svd
    /// \param mc
    ///
    DataClass(const int svd, const bool mc);
    ///
    /// \brief DataClass
    /// \param svd
    ///
    explicit DataClass(const int svd);
    ///
    /// \brief DataClass
    ///
    DataClass(void);
    ///
    /// \brief SetSVD
    /// \param v
    ///
    void SetSVD(const int v) {m_svd = v;}
    ///
    /// \brief SetMC
    /// \param v
    ///
    void SetMC(const bool v) {m_mc  = v;}
    ///
    /// \brief SetBp
    /// \param v
    ///
    void SetBp(const bool v) {m_bp  = v;}
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

 private:
    ///
    /// \brief m_svd
    ///
    int m_svd;
    ///
    /// \brief m_mc
    ///
    bool m_mc;
    ///
    /// \brief m_bp
    ///
    bool m_bp;
};

}  // namespace libTatami

#endif  // INCLUDE_DATACLASS_H_
