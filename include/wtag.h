/** Copyright 2016-2018 Vitaly Vorobyev
 ** @file wtag.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **
 **/

#pragma once

#include <cstdint>
#include <vector>
#include <string>

namespace libTatami {

///
/// \brief The WTag class. Class for calculation of delution factor
/// in tagging quality bins
///
class WTag {
    bool checkBin(uint16_t b) const;
    std::vector<double> m_w;
    std::vector<double> m_dw;
    static const std::vector<double> m_wbins;

 public:
    ///
    /// \brief WTag
    /// \param fname
    ///
    explicit WTag(const std::string &fname);
    ///
    /// \brief LoadParameters. Read wrong tagging probability coefficients
    /// from file
    /// \param fname
    /// \return
    ///
    int LoadParameters(const std::string& fname);
    ///
    /// \brief Delut. Delution factor for tagging quality bin
    /// \param bin
    /// \return
    ///
    double Delut(uint16_t bin) const;
    ///
    /// \brief Delut. Delution factor for tag quality variable
    /// \param q
    /// \return
    ///
    double Delut(double q) const;
    ///
    /// \brief WrProb. Wrong tagging probability
    /// \param q
    /// \return
    ///
    double WrProb(double q) const;
    ///
    /// \brief WrProb. Wrong tagging probability
    /// \param bin
    /// \return
    ///
    double WrProb(uint16_t bin) const;
    ///
    /// \brief GetBin. Wrong tagging bin
    /// \param q
    /// \return
    ///
    uint16_t GetBin(double q) const;
};

}  // namespace libTatami
