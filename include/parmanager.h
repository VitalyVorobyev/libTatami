/** Copyright 2016 Vitaly Vorobyev
 * @file parmanager.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#ifndef INCLUDE_PARMANAGER_H_
#define INCLUDE_PARMANAGER_H_

#include <string>

#include "./dataclass.h"

namespace libTatami {

typedef std::string str;

///
/// \brief The ParManager class
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
    static str WTagFile(const DataClass& dc);
    ///
    /// \brief BkgParFile
    /// \param dc
    /// \return
    ///
    static str BkgParFile(const DataClass& dc);
    ///
    /// \brief SigParFile
    /// \param dc
    /// \return
    ///
    static str SigParFile(const DataClass& dc);
    ///
    /// \brief prefix
    ///
    static std::string prefix;
};

}  // namespace libTatami

#endif  // INCLUDE_PARMANAGER_H_
