/** Copyright 2016 Vitaly Vorobyev
 * @file parmanager.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#pragma once

#include <string>

namespace libTatami {

class DataClass;

///
/// \brief The ParManager class
/// \todo move prefix to a text config file
///
class ParManager {
 public:
    ///
    /// \brief ParManager
    ///
    ParManager() {}
    ///
    /// \brief WTagFile
    /// \param dc
    /// \return
    ///
    static std::string WTagFile(const DataClass& dc);
    ///
    /// \brief BkgParFile
    /// \param dc
    /// \return
    ///
    static std::string BkgParFile(const DataClass& dc);
    ///
    /// \brief SigParFile
    /// \param dc
    /// \return
    ///
    static std::string SigParFile(const DataClass& dc);
    ///
    /// \brief prefix
    ///
    static std::string prefix;
};

}  // namespace libTatami
