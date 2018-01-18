/** Copyright 2016 Vitaly Vorobyev
 ** @file parmanager.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "parmanager.h"

#include "dataclass.h"

using std::stringstream;
using std::to_string;
using std::string;

namespace libTatami {

string ParManager::prefix("/home/vitaly/B0toD0pipi/libTatami/params/");

string ParManager::WTagFile(const DataClass& dc) {
    string fname = prefix + "wtag_svd" + to_string(dc.Svd());
    if (dc.MC()) fname += "_mc";
    return fname + ".txt";
}

string ParManager::BkgParFile(const DataClass& dc) {
    return prefix + string("def_bkg.txt");
}

string ParManager::SigParFile(const DataClass& dc) {
    string fname = prefix + "sig_" + (dc.Bp() ? "bp" : "b0") +
                "_svd" + to_string(dc.Svd());
    if (dc.MC()) fname += "_mc";
    return fname + ".txt";
}

}  // namespace libTatami
