/**
 * @file toypdfgen.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 * Copyright 2016 Vitaly Vorobyev
 */

#ifndef INCLUDE_TOYPDFGEN_H_
#define INCLUDE_TOYPDFGEN_H_

#include <vector>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
// #include <chrono>

#include "./abspdf.h"

namespace libTatami {

#include <cstdint>

typedef std::uniform_real_distribution<double> unidist;
typedef std::default_random_engine             rndmeng;
typedef std::vector<double>                    vectd;

///
/// \brief The ToyPdfGen class. Generator for toy studies.
/// Works with any class derived from AbsPdf
///
class ToyPdfGen {
 public:
    ///
    /// \brief ToyPdfGen. Object of class derivative from AbsPdf
    ///  should be provided to constructor
    /// \param pdf
    ///
    explicit ToyPdfGen(AbsPdf* pdf);
    ///
    /// \brief Generate. Generate N events and write vector of dt
    ///  variables in vec
    /// \param N
    /// \param vec
    /// \param silent
    /// \return
    ///
    int Generate(const int N, vectd* vec, const bool silent = false);
    ///
    /// \brief SetSeed. Initialize random numbers generator
    /// \param seed
    ///
    void SetSeed(const unsigned seed);

 private:
    ///
    /// \brief init
    ///
    void init(void);
    ///
    /// \brief unif
    ///
    unidist* unif;
    ///
    /// \brief re
    ///
    rndmeng  re;
    ///
    /// \brief m_seed
    ///
    unsigned m_seed;
    ///
    /// \brief m_maxtries
    ///
    int64_t m_maxtries;
    ///
    /// \brief m_pdf
    ///
    AbsPdf* m_pdf;
};

}  // namespace libTatami

#endif  // INCLUDE_TOYPDFGEN_H_
