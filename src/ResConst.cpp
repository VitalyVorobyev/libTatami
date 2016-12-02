/** Copyright 2016 Vitaly Vorobyev
 ** @file ResConst.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/ResConst.h"

#include <fstream>

using std::cout;
using std::endl;
using std::getline;
using std::sscanf;
using std::ifstream;

namespace libTatami {

ResConst::ResConst(const str& fname) {
    LoadParameters(fname);
}

int ResConst::LoadParameters(const str& fname) {
    ifstream ifile(fname.c_str(), ifstream::in);
    if (!ifile.is_open()) {
        cout << "ResConst::LoadParameters: can't open file " << fname << endl;
        return -1;
    } else {
        cout << "Getting dt resolution from file " << fname << endl;
    }
    str line, name;
    double val, errp, errn;
    char namech[15];
    int counter = 0;
    for (int i = 0; i < 13; i++) {
        getline(ifile, line);
        if (4 != sscanf(line.c_str(), "%s = %lf %lf %lf",
                        namech, &val, &errp, &errn)) {
            if (2 != sscanf(line.c_str(), "%s = %lf", namech, &val)) {
                cout << "ResConst::LoadParameters: can't parse string " <<
                        line << endl;
                continue;
            }
        }
        name = str(namech);
        cout << name << " " << val << endl;
        if (name == str("Srec0")) {
            m_Srec[0]           = val; counter++; continue;}
        if (name == str("Srec1")) {
            m_Srec[1]           = val; counter++; continue;}
        if (name == str("Sasc0")) {
            m_Sasc[0]           = val; counter++; continue;}
        if (name == str("Sasc1")) {
            m_Sasc[1]           = val; counter++; continue;}
        if (name == str("Snp_global")) {
            m_Snp_global        = val; counter++; continue;}
        if (name == str("Smn_rec")) {
            m_Smn_rec           = val; counter++; continue;}
        if (name == str("Stl_rec")) {
            m_Stl_rec           = val; counter++; continue;}
        if (name == str("ftl_rec")) {
            m_ftl_rec           = val; counter++; continue;}
        if (name == str("Smn_asc")) {
            m_Smn_asc           = val; counter++; continue;}
        if (name == str("Stl_asc")) {
            m_Stl_asc           = val; counter++; continue;}
        if (name == str("ftl_asc")) {
            m_ftl_asc           = val; counter++; continue;}

        if (name == str("fd_np_sgl")) {
            m_fd_np_sgl         = val; counter++; continue;}
        if (name == str("tau_np_p_sgl0")) {
            m_tau_np_p_sgl[0]   = val; counter++; continue;}
        if (name == str("tau_np_p_sgl1")) {
            m_tau_np_p_sgl[1]   = val; counter++; continue;}
        if (name == str("tau_np_n_sgl0")) {
            m_tau_np_n_sgl[0]   = val; counter++; continue;}
        if (name == str("tau_np_n_sgl1")) {
            m_tau_np_n_sgl[1]   = val; counter++; continue;}
        if (name == str("fd_np_mlt0")) {
            m_fd_np_mlt[0]      = val; counter++; continue;}
        if (name == str("fd_np_mlt1")) {
            m_fd_np_mlt[1]      = val; counter++; continue;}
        if (name == str("fd_np_st_mlt")) {
            m_fd_np_st_mlt      = val; counter++; continue;}
        if (name == str("fd_np_xi_mlt")) {
            m_fd_np_xi_mlt      = val; counter++; continue;}
        if (name == str("fd_np_stxi_mlt")) {
            m_fd_np_stxi_mlt    = val; counter++; continue;}
        if (name == str("fn_np_mlt")) {
            m_fn_np_mlt         = val; counter++; continue;}
        if (name == str("tau_np_p_mlt0")) {
            m_tau_np_p_mlt[0]   = val; counter++; continue;}
        if (name == str("tau_np_p_mlt1")) {
            m_tau_np_p_mlt[1]   = val; counter++; continue;}
        if (name == str("tau_np_p_xi_mlt")) {
            m_tau_np_p_xi_mlt   = val; counter++; continue;}
        if (name == str("tau_np_p_stxi_mlt")) {
            m_tau_np_p_stxi_mlt = val; counter++; continue;}
        if (name == str("tau_np_n_mlt0")) {
            m_tau_np_n_mlt[0]   = val; counter++; continue;}
        if (name == str("tau_np_n_mlt1")) {
            m_tau_np_n_mlt[1]   = val; counter++; continue;}
        if (name == str("tau_np_n_xi_mlt")) {
            m_tau_np_n_xi_mlt   = val; counter++; continue;}
        if (name == str("tau_np_n_stxi_mlt")) {
            m_tau_np_n_stxi_mlt = val; counter++; continue;}

        if (name == str("sigma_ol")) {
            m_sigma_ol          = val; counter++; continue;}
        if (name == str("f_ol_sgl")) {
            m_f_ol_sgl          = val; counter++; continue;}
        if (name == str("f_ol_mul")) {
            m_f_ol_mul          = val; counter++; continue;}

        if (name == str("rnp_kink_xi")) {
            m_rnp_kink_xi       = val; counter++; continue;}
        if (name == str("rnp_kink_st")) {
            m_rnp_kink_st       = val; counter++; continue;}
    }
    return counter;
}

}  // namespace libTatami
