#ifndef TTOOLS_H
#define TTOOLS_H

///
/// /brief The TTools class. Proides some constants and auxiliary function
///

class TTools{
public:
  TTools(void) {}

  /// sqrt of sum of squares
  static double sum_sigma(const double& s1, const double& s2);
  /// set double sided limits on variable
  static void constraint(double& x, const double& ll, const double& ul);
  /// set left side limit on variable
  static void constraint(double& x, const double& ll);

  /// coefficient to transform cm to ps
  static const double cm2ps;
  /// beta of Y(4S) at Belle
  static const double beta;
  /// mass of B0 meson (GeV)
  static const double mbzero;
};

#endif // TTOOLS_H
