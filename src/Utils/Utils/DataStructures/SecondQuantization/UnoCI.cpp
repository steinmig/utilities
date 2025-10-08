/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "UnoCI.h"
#include "CasSpecifier.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <deque>

namespace Scine {
namespace Utils {
namespace SecondQuantization {
namespace UnoCI {

std::tuple<Eigen::MatrixXd, Eigen::VectorXd> solveOrthogonal(const Eigen::MatrixXd& alphaDensity,
                                                             const Eigen::MatrixXd& betaDensity) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(alphaDensity + betaDensity);
  return {es.eigenvectors().rowwise().reverse(), es.eigenvalues().reverse()};
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd> solve(const Eigen::MatrixXd& alphaDensity,
                                                   const Eigen::MatrixXd& betaDensity, const Eigen::MatrixXd& overlap) {
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(alphaDensity + betaDensity, overlap);
  return {es.eigenvectors().rowwise().reverse(), es.eigenvalues().reverse()};
}

template<typename FilterType>
std::tuple<CasSpecifier<Reference::Restricted>, Eigen::VectorXd> getCas(const Eigen::MatrixXd& alphaDensity,
                                                                        const Eigen::MatrixXd& betaDensity,
                                                                        const Eigen::MatrixXd& overlap, FilterType filter) {
  Eigen::MatrixXd coefficients;
  Eigen::VectorXd occupations;
  int casElectrons{};
  std::vector<int> casIndices;

  std::tie(coefficients, occupations) =
      overlap.size() == 0 ? solveOrthogonal(alphaDensity, betaDensity) : solve(alphaDensity, betaDensity, overlap);
  auto totElectrons = std::lrint(occupations.sum());

  std::tie(casIndices, casElectrons) = filter(occupations);

  auto mos = MolecularOrbitals::createFromRestrictedCoefficients(coefficients);
  CasGenerator generator(mos);
  return {generator.generateRestricted(totElectrons, casIndices), occupations};
}

auto ThresholdFilter::operator()(const Eigen::VectorXd& occupations) const -> std::pair<std::vector<int>, long int> {
  std::vector<int> casIndices;
  double lowThresh = std::min(thresh.first, thresh.second) - 1e-12;
  double highThresh = std::max(thresh.first, thresh.second) + 1e-12;
  double casElectrons = 0;
  for (int i = 0; i < occupations.size(); ++i) {
    if (occupations(i) >= lowThresh && occupations(i) <= highThresh) {
      casIndices.push_back(i);
      casElectrons += occupations(i);
    }
  }
  return {casIndices, std::lrint(casElectrons)};
}

auto CasSizeFilter::operator()(const Eigen::VectorXd& occupations) const -> std::pair<std::vector<int>, long int> {
  if (casSize < 2) {
    throw std::runtime_error("Cas size of < 2 orbitals is too small!");
  }
  std::vector<int> casIndices;
  std::vector<int> occIndices;
  std::deque<int> virIndices;
  std::vector<int> somoIndices;
  for (int i = 0; i < occupations.size(); ++i) {
    if (occupations(i) == 1) {
      somoIndices.push_back(i);
    }
    else if (occupations(i) > 1) {
      occIndices.push_back(i);
    }
    else {
      virIndices.push_back(i);
    }
  }
  assert(somoIndices.size() % 2 == 0);

  double nElectronsInCas = 0.0;
  int orbitalsToChoose = std::min(casSize, int(occupations.size()));
  while (orbitalsToChoose > 0) {
    if (int(somoIndices.size()) > orbitalsToChoose) {
      casIndices.push_back(somoIndices[somoIndices.size() / 2]);
      nElectronsInCas += occupations(somoIndices[somoIndices.size() / 2]);
      somoIndices.erase(somoIndices.begin() + somoIndices.size() / 2);
      --orbitalsToChoose;
    }
    else {
      casIndices.push_back(occIndices.back());
      --orbitalsToChoose;
      casIndices.push_back(virIndices.front());
      --orbitalsToChoose;
      nElectronsInCas += occupations(occIndices.back()) + occupations(virIndices.front());

      occIndices.pop_back();
      virIndices.pop_front();
    }
  }
  std::sort(casIndices.begin(), casIndices.end());
  return {casIndices, std::lrint(nElectronsInCas)};
}

template std::tuple<CasSpecifier<Reference::Restricted>, Eigen::VectorXd>
getCas<ThresholdFilter>(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, ThresholdFilter);
template std::tuple<CasSpecifier<Reference::Restricted>, Eigen::VectorXd>
getCas<CasSizeFilter>(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, CasSizeFilter);
} // namespace UnoCI
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
