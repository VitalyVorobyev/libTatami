/** Copyright 2016 Vitaly Vorobyev
 ** @file wtag.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "wtag.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using std::cout;
using std::endl;
using std::ifstream;
using std::string;

namespace libTatami {

const std::vector<double> WTag::m_wbins =
    {0.0, 0.1, 0.25, 0.5, 0.625, 0.75, 0.875, 1.01};

WTag::WTag(const string &fname) {
    LoadParameters(fname);
}

int WTag::LoadParameters(const string& fname) {
    ifstream ifile(fname.c_str(), ifstream::in);
    if (!ifile.is_open()) {
        cout << "WTag::LoadParameters: can't open file " << fname << endl;
        return -1;
    } else {
        cout << "Getting wtag from file " << fname << endl;
    }
    m_w.clear();
    m_w.push_back(0.5);
    string line;
    double val, errp, errn;
    while (!ifile.eof()) {
        getline(ifile, line);
        if (3 != sscanf(line.c_str(), "%lf %lf %lf", &val, &errp, &errn)) {
            if (1 != sscanf(line.c_str(), "%lf", &val)) {
                cout << "WTag::LoadParameters: can't parse string " <<
                        line << endl;
                continue;
            }
        }
        m_w.push_back(val);
        cout << " " << val << endl;
    }
    return m_w.size();
}

double WTag::Delut(uint16_t bin) const {
    if (!checkBin(bin)) {
        cout << "WTag::Delut: wrong bin " << bin << endl;
        return 0;
    }
    return 1. - 2. * m_w[bin];
}

double WTag::Delut(double q) const {
    return Delut(GetBin(q));
}

double WTag::WrProb(double q) const {
    return WrProb(GetBin(q));
}

double WTag::WrProb(uint16_t bin) const {
    if (!checkBin(bin)) {
        cout << "WTag::WrProb: wrong bin " << bin << endl;
        return 0;
    }
    return m_w[bin];
}

uint16_t WTag::GetBin(double q) const {
    if (q <= 0 || q > 1.) return 0;
    return upper_bound(m_wbins.begin(), m_wbins.end(), q) - m_wbins.begin() - 1;
}

bool WTag::checkBin(uint16_t b) const {
    return (b < 0 || b >= m_w.size());
}

}  // namespace libTatami

// const double wbins[8] =
// {0.0,0.1,0.25,0.5,0.625,0.75,0.875,1.01};
// const double w_mc_svd1[7] =
// {0.5,0.420827,0.300296,0.219317,0.154636,0.0916131,0.0228891};
// const double w_data_svd1[7] =
// {0.5,0.418852,0.329879,0.233898,0.170608,0.099791, 0.0228501};
// const double w_mc_svd2[7] =
// {0.5,0.412222,0.307838,0.212765,0.149933,0.0913264,0.0218754};
// const double w_data_svd2[7] =
// {0.5,0.418826,0.319303,0.222948,0.163191,0.104085, 0.0251454};

// const double dw_mc_svd1[7] =
// {0.0, 0.0583019, 0.00573998,-0.0392635, 0.00474508,-0.0118737, -0.00585326};
// const double dw_data_svd1[7] =
// {0.0, 0.0569661, 0.0126192, -0.0147724,-0.000550289,0.00887704, 0.00465683};
// const double dw_mc_svd2[7] =
// {0.0, 0.00408778,0.010326,  -0.00479522,0.00151989, 0.0143633,  0.00189979};
// const double dw_data_svd2[7] =
// {0.0,-0.00877001,0.0103515, -0.0109253,-0.0186365,  0.00168037,-0.0036441};

// const double w_data_svd1_posi[7] =
// {0.,7.235697e-03,7.129388e-03,7.417778e-03,
//  6.885875e-03,6.761047e-03,4.336734e-03};
// const double w_data_svd1_nega[7] =
// {0.,6.001569e-03,6.430566e-03,7.693083e-03,
//  6.416449e-03,8.807757e-03,4.587614e-03};
// const double w_data_svd2_posi[7] =
// {0.,4.152612e-03,3.243236e-03,3.721417e-03,
//  3.315138e-03,3.180302e-03,2.175087e-03};
// const double w_data_svd2_nega[7] =
// {0.,3.577812e-03,2.803811e-03,3.486607e-03,
//  4.241595e-03,3.696399e-03,3.077622e-03};

// const double dw_data_svd1_posi[7] =
// {0.,3.961548e-03,3.543129e-03,4.129422e-03,
//  4.169570e-03,3.998982e-03,2.433324e-03};
// const double dw_data_svd1_nega[7] =
// {0.,3.927049e-03,3.698619e-03,4.179366e-03,
//  4.602366e-03,3.914627e-03,2.360543e-03};

// double w[7],dw[7];
