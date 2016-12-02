/** Copyright 2016 Vitaly Vorobyev
 ** @file bkgpdfparset.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include <iostream>
#include <fstream>

#include "../include/bkgpdfparset.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::getline;
using std::sscanf;

namespace libTatami {

BkgPDFParSet::BkgPDFParSet():
    m_S_main_mlt(1.), m_S_tail_mlt(3.), m_f_tail_mlt(0.2), m_f_delta_mlt(0.7),
    m_S_main_sgl(1.), m_S_tail_sgl(3.), m_f_tail_sgl(0.2), m_f_delta_sgl(0.6),
    m_mu(0.), m_mu_delta(0.), m_tau(1.3), m_f_otlr(0.), m_s_otlr(30.)
{}

int BkgPDFParSet::GetParametersFromFile(const str& fname) {
    ifstream ifile(fname.c_str(), ifstream::in);
    if (!ifile.is_open()) {
        cout << "Can't open file " << fname << endl;
        return -1;
    } else {
        cout << "Getting background description from file " << fname << endl;
    }
    str line, name;
    double val, err;
    char namech[15];
    int counter = 0;
    for (int i = 0; i < 13; i++) {
        getline(ifile, line);
        sscanf(line.c_str(), "%s = %lf +- %lf", namech, &val, &err);
        name = str(namech);
        cout << name << " " << val << " +- " << err << endl;
        if (name == "tau")         { m_tau         = val; counter++; continue;}
        if (name == "mu")          { m_mu          = val; counter++; continue;}
        if (name == "mu_delta")    { m_mu_delta    = val; counter++; continue;}
        if (name == "f_delta_mlt") { m_f_delta_mlt = val; counter++; continue;}
        if (name == "f_tail_mlt")  { m_f_tail_mlt  = val; counter++; continue;}
        if (name == "S_main_mlt")  { m_S_main_mlt  = val; counter++; continue;}
        if (name == "S_tail_mlt")  { m_S_tail_mlt  = val; counter++; continue;}
        if (name == "f_delta_sgl") { m_f_delta_sgl = val; counter++; continue;}
        if (name == "f_tail_sgl")  { m_f_tail_sgl  = val; counter++; continue;}
        if (name == "S_main_sgl")  { m_S_main_sgl  = val; counter++; continue;}
        if (name == "S_tail_sgl")  { m_S_tail_sgl  = val; counter++; continue;}
        if (name == "f_otlr")      { m_f_otlr      = val; counter++; continue;}
        if (name == "s_otlr")      { m_s_otlr      = val; counter++; continue;}
    }
    return counter;
}

}  // namespace libTatami
