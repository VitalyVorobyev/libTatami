/** Copyright 2016 Vitaly Vorobyev
 ** @file abspdf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#pragma once

#include "cnvl.h"

namespace libTatami {

class ICPVEvt;

///
/// \brief The AbsPdf class. Abstract class for PDF based on cnvl library
///
class AbsPdf : public cnvl {
    ///
    /// \brief m_ll. lower dt limit
    ///
    double m_ll;
    ///
    /// \brief m_ul. upper dt limit
    ///
    double m_ul;

 public:
    ///
    /// \brief AbsPdf
    ///
    AbsPdf();
    ///
    /// \brief operator (). Virtual method for calculation of PDF
    /// \param ext
    /// \return
    ///
    virtual double operator() (const ICPVEvt& ext) const = 0;
    ///
    /// \brief operator (). Virtual method for calculation of PDF
    /// \param x
    /// \return
    ///
    virtual double operator() (double x) const = 0;
    ///
    /// \brief ll. Lower limit of dt vatiable
    /// \return
    ///
    double ll(void) const {return m_ll;}
    ///
    /// \brief ul. Upper limit of dt vatiable
    /// \return
    ///
    double ul(void) const {return m_ul;}
    ///
    /// \brief SetRange. Set symmetrical range for dt variable
    /// \param v
    ///
    void SetRange(double v) { m_ul = v, m_ll = -v;}
    ///
    /// \brief SetRange. Set range for dt variable
    /// \param min
    /// \param max
    ///
    void SetRange(double min, double max) {
        m_ll = min, m_ul = max;
    }
};

}  // namespace libTatami
