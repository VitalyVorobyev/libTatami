/** Copyright 2016 Vitaly Vorobyev
 ** @file toypdfgen.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "toypdfgen.h"

#include <algorithm>

#include "abspdf.h"

using std::cout;
using std::endl;
using std::cerr;
using std::random_device;
using std::max;
using std::vector;

namespace libTatami {

ToyPdfGen::ToyPdfGen(const AbsPdf& pdf):
  m_seed(0), m_maxtries(1e7), m_pdf(pdf) {
    init();
}

double ToyPdfGen::FindMaj(uint64_t N) {
    double maj = 0;
    for (uint64_t i = 0; i < N; i++) {
        double dt = (*unif)(re);
        double pdf = m_pdf(dt);
        if (pdf > maj) maj = pdf;
    }
    return maj;
}

vector<double> ToyPdfGen::Generate(uint64_t N, bool silent) {
    vector<double> vec;
    double maj = 1.1 * FindMaj(max(N / 5, 5000ul));
    if (!silent) cout << "Majorant: " << maj << endl;
    std::uniform_real_distribution<double> unifMaj(0., maj);
    std::default_random_engine ren;

    uint64_t tries = 0;
    uint64_t Evtn = 0;
    if (!silent) cout << "Generating " << N << " events..." << endl;
    while (Evtn < N && tries++ < m_maxtries) {
        double dt = (*unif)(re);
        double xi = unifMaj(ren);
        double pdf = m_pdf(dt);
        if (pdf > maj) {
            cerr << "Update maj: " << maj << " -> " << 1.1 * pdf << endl;
            maj = 1.1 * pdf;
        }
        if (xi < pdf) {
          vec.emplace_back(dt);
          Evtn++;
          if (!silent && !(Evtn % 10000)) cout << Evtn << " events" << endl;
        }
    }
    const double Eff = static_cast<double>(Evtn) / tries;
    if (tries == m_maxtries) {
        cout << "ToyPdfGen::Generate: tries limit exceed!" << endl
             << "  majorant:         " << maj        << endl
             << "  efficiency:       " << Eff        << endl
             << "  events generated: " << Evtn       << endl
             << "  tries limit: "      << m_maxtries << endl;
    }
    if (!silent) cout << "Done!. Efficiency " << Eff << endl;
    return std::move(vec);
}

void ToyPdfGen::init(void) {
    unif = std::make_unique<std::uniform_real_distribution<double>>(
                m_pdf.ul(), m_pdf.ll());
    SetSeed(m_seed);
}

void ToyPdfGen::SetSeed(uint32_t seed) {
  m_seed = seed;
  if (!seed) re.seed(random_device {}());
  else       re.seed(m_seed);
}

}  // namespace libTatami
