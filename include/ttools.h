/** Copyright 2016 Vitaly Vorobyev
 ** @file ttools.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **
 **/

#ifndef INCLUDE_TTOOLS_H_
#define INCLUDE_TTOOLS_H_

namespace libTatami {

///
/// /brief The TTools class. Constants and auxiliary function
///
class TTools {
 public:
    ///
    /// \brief TTools
    ///
    TTools(void) {}
    /// sqrt of sum of squares
    /// \brief sum_sigma
    /// \param s1
    /// \param s2
    /// \return
    ///
    static double sum_sigma(const double& s1, const double& s2);
    ///
    /// \brief constraint. set double sided limits on variable
    /// \param x
    /// \param ll
    /// \param ul
    ///
    static void constraint(double* x, const double& ll, const double& ul);
    ///
    /// \brief constraint. set left side limit on variable
    /// \param x
    /// \param ll
    ///
    static void constraint(double* x, const double& ll);
    ///
    /// \brief cm2ps. coefficient to transform cm to ps
    ///
    static const double cm2ps;
    ///
    /// \brief beta of Y(4S) at Belle
    ///
    static const double beta;
    ///
    /// \brief mbzero. Mass of B0 meson (GeV)
    ///
    static const double mbzero;
};

}  // namespace libTatami

#endif  // INCLUDE_TTOOLS_H_
