#ifndef TOYPDFGEN_H
#define TOYPDFGEN_H

#include "abspdf.h"

#include <vector>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>

///
/// \brief The ToyPdfGen class. Generator for toy studies.
/// Works with any class derived from AbsPdf
///

class ToyPdfGen{
public:
  /// Object of class derivative from AbsPdf should be provided to constructor
  ToyPdfGen(AbsPdf* pdf);
  /// Generate N events and write vector of dt variables in vec
  int Generate(const int N,std::vector<double>& vec);
  /// Initialize random numbers generator
  void SetSeed(const unsigned seed);

private:
  void init(void);
  std::uniform_real_distribution<double>* unif;
  std::default_random_engine re;
  unsigned m_seed;
  long m_maxtries;
  AbsPdf* m_pdf;
};

#endif // TOYPDFGEN_H
