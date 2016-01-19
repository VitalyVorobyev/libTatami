#include "toypdfgen.h"

using namespace std;

ToyPdfGen::ToyPdfGen(AbsPdf *pdf):
  m_seed(0), m_maxtries(1e7), m_pdf(pdf)
{
  init();
}

int ToyPdfGen::Generate(const int N,vector<double>& vec){
  const double maj = (*m_pdf)(0);
  cout << "Majorant: " << maj << endl;
  vec.clear();
  uniform_real_distribution<double> unifMaj(0.,maj);
  default_random_engine ren;

  long tries = 0;
  int Evtn = 0;
  double dt,xi;
  cout << "Generating " << N << " events..." << endl;
  while(Evtn < N && tries < m_maxtries){
    dt = (*unif)(re);
    xi = unifMaj(ren);
    if(xi<(*m_pdf)(dt)){
      vec.push_back(dt);
      if(!(++Evtn%1000)) cout << Evtn << " events" << endl;
    }
    tries++;
  }
  const double Eff = (double)Evtn/tries;
  if(tries == m_maxtries){
    cout << "ToyPdfGen::Generate: tries limit exceed!";
    cout << "  efficiency: "       << Eff << endl;
    cout << "  events generated: " << Evtn << endl;
    cout << "  tries limit: "      << m_maxtries << endl;
    return -1;
  }
  cout << "Done!. Efficiency " << Eff << endl;
  return 0;
}

void ToyPdfGen::init(void){
  unif = new std::uniform_real_distribution<double>(m_pdf->ul(),m_pdf->ll());
  SetSeed(m_seed);
}

void ToyPdfGen::SetSeed(const unsigned seed){
  m_seed = seed;
  re.seed(m_seed);
}
