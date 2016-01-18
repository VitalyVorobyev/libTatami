#ifndef WTAG_H
#define WTAG_H

#include <vector>
#include <string>

///
/// \brief The WTag class. Class for calculation of delution factor in tagging quality bins
///

class WTag{
public:
  WTag(const std::string &fname);
  /// Read wrong tagging probability coefficients from file
  int LoadParameters(const std::string& fname);

  /// Delution factor for tagging quality bin
  double Delut(const int bin);
  /// Delution factor for tag quality variable
  double Delut(const double& q);

private:
  int GetBin(const double& q);
  std::vector<double> m_w;
  std::vector<double> m_dw;
  static const double m_wbins[8];
};

#endif // WTAG_H
