/** Copyright 2016 Vitaly Vorobyev
 ** @file wtag.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **
 **/

#ifndef INCLUDE_WTAG_H_
#define INCLUDE_WTAG_H_

#include <vector>
#include <string>

#include "./typedefs.h"

namespace libTatami {

///
/// \brief The WTag class. Class for calculation of delution factor
/// in tagging quality bins
///
class WTag {
 public:
    ///
    /// \brief WTag
    /// \param fname
    ///
    explicit WTag(const str &fname);
    ///
    /// \brief LoadParameters. Read wrong tagging probability coefficients
    /// from file
    /// \param fname
    /// \return
    ///
    int LoadParameters(const str& fname);
    ///
    /// \brief Delut. Delution factor for tagging quality bin
    /// \param bin
    /// \return
    ///
    double Delut(const unsigned bin) const;
    ///
    /// \brief Delut. Delution factor for tag quality variable
    /// \param q
    /// \return
    ///
    double Delut(const double& q) const;
    ///
    /// \brief WrProb. Wrong tagging probability
    /// \param q
    /// \return
    ///
    double WrProb(const double& q) const;
    ///
    /// \brief WrProb. Wrong tagging probability
    /// \param bin
    /// \return
    ///
    double WrProb(const unsigned bin) const;
    ///
    /// \brief GetBin. Wrong tagging bin
    /// \param q
    /// \return
    ///
    unsigned GetBin(const double& q) const;

 private:
    bool checkBin(const unsigned b) const;
    vectd m_w;
    vectd m_dw;
    static const vectd m_wbins;
};

}  // namespace libTatami

#endif  // INCLUDE_WTAG_H_
