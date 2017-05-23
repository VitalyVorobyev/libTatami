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

typedef std::vector<double> vectd;

using std::cout;
using std::endl;
using std::cerr;
using std::random_device;

namespace libTatami {

ToyPdfGen::ToyPdfGen(AbsPdf *pdf):
  m_seed(0), m_maxtries(1e7), m_pdf(pdf) {
    init();
}

double ToyPdfGen::FindMaj(const uint64_t N) {
    double maj = 0;
    for (uint64_t i = 0; i < N; i++) {
        double dt = (*unif)(re);
        double pdf = (*m_pdf)(dt);
        if (pdf > maj) maj = pdf;
    }
    return maj;
}

int ToyPdfGen::Generate(const uint64_t N, vectd* vec, const bool silent) {
    double maj = 1.1 * FindMaj(N / 5);
    if (!silent) cout << "Majorant: " << maj << endl;
    vec->clear();
    unidist unifMaj(0., maj);
    rndmeng ren;

    uint64_t tries = 0;
    uint64_t Evtn = 0;
    if (!silent) cout << "Generating " << N << " events..." << endl;
    while (Evtn < N && tries++ < m_maxtries) {
        double dt = (*unif)(re);
        double xi = unifMaj(ren);
        double pdf = (*m_pdf)(dt);
        if (pdf > maj) {
            cerr << "Update maj: " << maj << " -> " << 1.1*pdf << endl;
            maj = 1.1*pdf;
        }
        if (xi < pdf) {
          vec->push_back(dt);
          Evtn++;
          if (!silent && !(Evtn % 10000)) cout << Evtn << " events" << endl;
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
//    if (!silent) cout << "Done!. Efficiency " << Eff << endl;
    return 0;
}

void ToyPdfGen::init(void) {
    unif = new unidist(m_pdf->ul(), m_pdf->ll());
    SetSeed(m_seed);
}

void ToyPdfGen::SetSeed(const uint32_t seed) {
  m_seed = seed;
  if (!seed) re.seed(random_device {}());
  else       re.seed(m_seed);
}

}  // namespace libTatami
