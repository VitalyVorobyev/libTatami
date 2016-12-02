/** Copyright 2016 Vitaly Vorobyev
 ** @file toypdfgen.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../include/toypdfgen.h"

using std::cout;
using std::endl;
using std::random_device;

namespace libTatami {

ToyPdfGen::ToyPdfGen(AbsPdf *pdf):
  m_seed(0), m_maxtries(1e7), m_pdf(pdf) {
    init();
}

int ToyPdfGen::Generate(const int N, vectd* vec, const bool silent) {
    const double maj = (*m_pdf)(0);
    if (!silent) cout << "Majorant: " << maj << endl;
    vec->clear();
    unidist unifMaj(0., maj);
    rndmeng ren;

    int64_t tries = 0;
    int64_t Evtn = 0;
    double dt, xi;
    if (!silent) cout << "Generating " << N << " events..." << endl;
    while (Evtn < N && tries++ < m_maxtries) {
        dt = (*unif)(re);
        xi = unifMaj(ren);
        if (xi < (*m_pdf)(dt)) {
          vec->push_back(dt);
          if (!(++Evtn % 1000)) cout << Evtn << " events" << endl;
        }
    }
    const double Eff = static_cast<double>(Evtn) / tries;
    if (tries == m_maxtries) {
        cout << "ToyPdfGen::Generate: tries limit exceed!" << endl;
        cout << "  majorant:         " << maj        << endl;
        cout << "  efficiency:       " << Eff        << endl;
        cout << "  events generated: " << Evtn       << endl;
        cout << "  tries limit: "      << m_maxtries << endl;
        return -1;
    }
    if (!silent) cout << "Done!. Efficiency " << Eff << endl;
    return 0;
}

void ToyPdfGen::init(void) {
    unif = new unidist(m_pdf->ul(), m_pdf->ll());
    SetSeed(m_seed);
}

void ToyPdfGen::SetSeed(const unsigned seed) {
  m_seed = seed;
  if (!seed) re.seed(random_device {}());
  else       re.seed(m_seed);
}

}  // namespace libTatami
