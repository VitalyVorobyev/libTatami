/** Copyright 2016 Vitaly Vorobyev
 ** @file abspdf.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "abspdf.h"

#include <cmath>

using std::fabs;

constexpr double eps = 0.000001;

namespace libTatami {

AbsPdf::AbsPdf(void) : cnvl(), m_ll(-70), m_ul(70) {}

std::vector<double> AbsPdf::tabulate(uint32_t ndots,
                                     double lo, double hi) const {
    if (fabs(lo) < eps) lo = m_ll;
    if (fabs(hi) < eps) hi = m_ul;
    if (lo > hi) {
        lo = m_ll;
        hi = m_ul;
    }
    if (!ndots) ndots = 1000;
    const double dt = (hi - lo) / ndots;
    std::vector<double> table(ndots);
    auto t = lo;
    for (auto tit = table.begin(); tit != table.end(); tit++, t += dt)
        *tit = (*this)(t);
    return std::move(table);
}

}  // namespace libTatami
