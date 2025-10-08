/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SPIN_SQUARED_EVALUATOR_H
#define UTILS_SPIN_SQUARED_EVALUATOR_H

#include "ElectronicDeterminant.h"
#include "HermitianHamiltonian.h"
#include <boost/functional/hash.hpp>
#include <unordered_map>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

namespace SpinSquaredEvaluation {
using Wavefunction = std::unordered_map<ElectronicDeterminant, double, ElectronicDeterminantHash>;
enum class DeterminantRepresentation { OrbitalAlphaOrbitalBeta, OnvAlphaOnvBeta };

inline auto applySMinus(const Wavefunction& ciWavefunction, DeterminantRepresentation rep) -> Wavefunction {
  Wavefunction result;
  for (const auto& det : ciWavefunction) {
    for (int p : det.first.getAllOccupied(SpinComponent::Alpha)) {
      if (!det.first.isOccupied(SpinComponent::Beta, p)) {
        auto tmp = det.first;
        tmp.applyExcitation(CrossExcitation(p, p, {SpinComponent::Alpha, SpinComponent::Beta}));
        result[tmp] += rep == DeterminantRepresentation::OrbitalAlphaOrbitalBeta ? det.second : -det.second;
      }
    }
  }
  return result;
}

inline auto applySPlus(const Wavefunction& ciWavefunction) -> Wavefunction {
  Wavefunction result;
  for (const auto& det : ciWavefunction) {
    for (int q : det.first.getAllOccupied(SpinComponent::Beta)) {
      if (!det.first.isOccupied(SpinComponent::Alpha, q)) {
        auto tmp = det.first;
        tmp.applyExcitation(CrossExcitation(q, q, {SpinComponent::Beta, SpinComponent::Alpha}));
        result[tmp] += det.second;
      }
    }
  }
  return result;
}
inline auto getMz(const Wavefunction& ciWavefunction) -> double {
  const auto nElectrons = ciWavefunction.begin()->first.countElectrons();
  return 0.5 * (nElectrons.first - nElectrons.second);
}

inline auto overlap(const Wavefunction& ciWavefunction, const Wavefunction& ciWavefunction2) -> double {
  const auto& smallerWF = ciWavefunction.size() < ciWavefunction2.size() ? ciWavefunction : ciWavefunction2;
  const auto& biggerWF = ciWavefunction.size() < ciWavefunction2.size() ? ciWavefunction2 : ciWavefunction;

  double overlap = 0.0;
  for (const auto& det : smallerWF) {
    auto it = biggerWF.find(det.first);
    if (it != biggerWF.end()) {
      overlap += det.second * it->second;
    }
  }
  return overlap;
}

inline auto evaluate(const Wavefunction& ciWavefunction,
                     DeterminantRepresentation rep = DeterminantRepresentation::OrbitalAlphaOrbitalBeta) -> double {
  double spinSPlusSMinus = overlap(ciWavefunction, applySPlus(applySMinus(ciWavefunction, rep)));
  auto mz = getMz(ciWavefunction);
  return spinSPlusSMinus + mz * (mz - 1);
}

} // namespace SpinSquaredEvaluation
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
#endif // UTILS_SPIN_SQUARED_EVALUATOR_H
