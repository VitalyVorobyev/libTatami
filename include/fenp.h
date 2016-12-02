/** Copyright 2016 Vitaly Vorobyev
 ** @file fenp.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef INCLUDE_FENP_H_
#define INCLUDE_FENP_H_

namespace libTatami {

///
/// \brief The fEnp class
///
class fEnp {
 public:
    ///
    /// \brief fEnp
    ///
    fEnp() {
        Clear();
    }
    ///
    /// \brief Clear
    ///
    void Clear();
    ///
    /// \brief fEn
    ///
    double fEn;
    ///
    /// \brief fEp
    ///
    double fEp;
    ///
    /// \brief fEn_np
    ///
    double fEn_np;
    ///
    /// \brief fEp_np
    ///
    double fEp_np;
    ///
    /// \brief fxEn
    ///
    double fxEn;
    ///
    /// \brief fxEp
    ///
    double fxEp;
    ///
    /// \brief fAn
    ///
    double fAn;
    ///
    /// \brief fAp
    ///
    double fAp;
    ///
    /// \brief fMn
    ///
    double fMn;
    ///
    /// \brief fMp
    ///
    double fMp;
    ///
    /// \brief fxEn_k
    ///
    double fxEn_k;
    ///
    /// \brief fxEp_k
    ///
    double fxEp_k;
    ///
    /// \brief fEn_k
    ///
    double fEn_k;
    ///
    /// \brief fEp_k
    ///
    double fEp_k;
};

}  // namespace libTatami

#endif  // INCLUDE_FENP_H_
