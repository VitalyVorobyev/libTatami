/** Copyright 2016 Vitaly Vorobyev
 ** @file parmanager.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/parmanager.h"

using std::stringstream;
using std::to_string;

namespace libTatami {

str ParManager::prefix("/home/vitaly/B0toD0pipi/libTatami/params/");

str ParManager::WTagFile(const DataClass& dc) {
    str fname = prefix + "wtag_svd" + to_string(dc.Svd());
    if (dc.MC()) fname += "_mc";
    return fname + ".txt";
}

str ParManager::BkgParFile(const DataClass& dc) {
    return prefix + str("def_bkg.txt");
}

str ParManager::SigParFile(const DataClass& dc) {
    str fname = prefix + "sig_" + (dc.Bp() ? "bp" : "b0") +
                "_svd" + to_string(dc.Svd());
    if (dc.MC()) fname += "_mc";
    return fname + ".txt";
}

}  // namespace libTatami
