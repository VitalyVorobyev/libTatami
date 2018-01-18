/** Copyright 2016 Vitaly Vorobyev
 *
 * @file toypdfgen.h
 *
 * @brief This message displayed in Doxygen Files index
 *
 * @author Vitaly Vorobyev
 * Contact: vit.vorobiev@gmail.com
 *
 */

#pragma once

#include <vector>
#include <random>       // std::default_random_engine
#include <cstdint>
#include <memory>

namespace libTatami {

class AbsPdf;

///
/// \brief The ToyPdfGen class. Generator for toy studies.
/// Works with any class derived from AbsPdf
///
class ToyPdfGen {
    /**
     * @brief FindMaj
     * @param N. How many tries of tries
     * @return
     */
    double FindMaj(uint64_t N);
    ///
    /// \brief init
    ///
    void init(void);
    ///
    /// \brief unif
    ///
    std::unique_ptr<std::uniform_real_distribution<double>> unif;
    ///
    /// \brief re
    ///
    std::default_random_engine  re;
    ///
    /// \brief m_seed
    ///
    uint32_t m_seed;
    ///
    /// \brief m_maxtries
    ///
    uint64_t m_maxtries;
    ///
    /// \brief m_pdf
    ///
    const AbsPdf& m_pdf;

 public:
    ///
    /// \brief ToyPdfGen. Object of class derivative from AbsPdf
    ///  should be provided to constructor
    /// \param pdf
    ///
    explicit ToyPdfGen(const AbsPdf& pdf);
    ///
    /// \brief Generate. Generate N events and write vector of dt
    ///  variables in vec
    /// \param N
    /// \param vec
    /// \param silent
    /// \return
    ///
    std::vector<double> Generate(uint64_t N, bool silent = false);
    ///
    /// \brief SetSeed. Initialize random numbers generator
    /// \param seed
    ///
    void SetSeed(uint32_t seed);
};

}  // namespace libTatami
