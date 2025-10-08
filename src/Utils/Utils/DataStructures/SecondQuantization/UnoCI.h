/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_UNOCI_H
#define UTILS_UNOCI_H

#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <Eigen/Core>
#include <tuple>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

template<Reference restrictedness>
class CasSpecifier;

namespace UnoCI {

struct ThresholdFilter {
  ThresholdFilter(std::pair<double, double> thresholds = {0.02, 1.98}) : thresh(std::move(thresholds)) {
  }
  std::pair<double, double> thresh;
  std::pair<std::vector<int>, long int> operator()(const Eigen::VectorXd& occupations) const;
};

struct CasSizeFilter {
  CasSizeFilter(int nOrbitalsInCas = 2) : casSize(nOrbitalsInCas) {
    assert(casSize >= 2);
  }
  std::pair<std::vector<int>, long int> operator()(const Eigen::VectorXd& occupations) const;
  int casSize;
};
/**
 * @brief Solves the orthogonal eigenvalue problem (D^a+D^b)C = nC.
 */
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> solveOrthogonal(const Eigen::MatrixXd& alphaDensity,
                                                             const Eigen::MatrixXd& betaDensity);
/**
 * @brief Solves the generalized eigenvalue problem (D^a+D^b)SC = nC.
 */
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> solve(const Eigen::MatrixXd& alphaDensity,
                                                   const Eigen::MatrixXd& betaDensity, const Eigen::MatrixXd& overlap);
/**
 * @brief Obtains a CAS with the Unrestricted Natural Orbital method, together with the NO occupations.
 * @returns std::tuple<CasSpecifier<Reference::Restricted>, Eigen::VectorXd> the CAS and the fractionak occupations
 * vector. If no overlap is given, an orthogonal AO basis is assumed. This is true for NDDO methods, for instance. The
 * CAS is then restricted, as only one set of natural orbitals is obtained. The UNOs can then be employed to decide
 * which orbitals are the ones that are the most active to get references determinants.
 */
template<typename FilterType = ThresholdFilter>
std::tuple<CasSpecifier<Reference::Restricted>, Eigen::VectorXd>
getCas(const Eigen::MatrixXd& alphaDensity, const Eigen::MatrixXd& betaDensity,
       const Eigen::MatrixXd& overlap = Eigen::MatrixXd(0, 0), FilterType filter = {});

/**
 * @brief Depending on the thresholds, decides which orbitals to retain for the CAS and which not.
 * It makes no difference which threshold is the high one and which the low one, it is automatically checked.
 */
std::pair<std::vector<int>, long int> filter(const Eigen::VectorXd& occupations, std::pair<double, double> thresholds);

} // namespace UnoCI
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_UNOCI_H
