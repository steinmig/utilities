/**
 * @file IdealGasTranslator.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IdealGasTranslator.h"
#include <Utils/Geometry/Utilities/Properties.h>

namespace Scine {
namespace Utils {

IdealGasTranslator::IdealGasTranslator(const AtomCollection& atoms) {
  moleculeMass_ = 0.0;
  for (const auto& atomMass : Geometry::Properties::getMasses(atoms.getElements())) {
    moleculeMass_ += atomMass;
  }
  moleculeMass_ *= 1e-3 / Constants::avogadroNumber / Constants::electronRestMass;
}
IdealGasTranslator::IdealGasTranslator(double molarMass)
  : moleculeMass_(molarMass * 1e-3 / Constants::avogadroNumber / Constants::electronRestMass) {
}

IdealGasTranslator::~IdealGasTranslator() = default;

Eigen::VectorXd IdealGasTranslator::densityOfState(const Eigen::VectorXd&) {
  // If we use the particle-in-a-box model, this expression will be pressure/volume dependent.
  throw std::runtime_error(
      "The density and the sum of states are not implemented for the ideal gas translation model at the moment.");
}
Eigen::VectorXd IdealGasTranslator::sumOfStates(const Eigen::VectorXd&) {
  // If we use the particle-in-a-box model, this expression will be pressure/volume dependent.
  throw std::runtime_error(
      "The density and the sum of states are not implemented for the ideal gas translation model at the moment.");
}
double IdealGasTranslator::partitionFunction(const ThermodynamicReferenceState& state) {
  if (state.temperature < 1e-6) {
    return 1;
  }
  if (state.pressure < 1e-6) {
    throw std::runtime_error("The classical partitioning function for an ideal gas translator at zero pressure is"
                             "not defined! Please use a pressure of more than 1e-6 Pa.");
  }
  constexpr double au_per_pascal =
      Constants::hartree_per_joule / (Constants::bohr_per_meter * Constants::bohr_per_meter * Constants::bohr_per_meter);
  // q = V * (m k T /(2pi))^3/2 ; V = k T / p (atomic units: h/(2pi) = 1)
  const double pressureAtomicUnits = state.pressure * au_per_pascal;
  double volume = kB * state.temperature / pressureAtomicUnits;
  return volume * std::pow(moleculeMass_ * kB * state.temperature / (2.0 * Constants::pi), 1.5);
}
double IdealGasTranslator::entropy(const ThermodynamicReferenceState& state) {
  return kB * (5.0 / 2.0 + std::log(this->partitionFunction(state)));
}
double IdealGasTranslator::enthalpy(const ThermodynamicReferenceState& state) {
  return 2.5 * kB * state.temperature;
}
double IdealGasTranslator::heatCapacityCp(const ThermodynamicReferenceState&) {
  return 2.5 * kB;
}
double IdealGasTranslator::heatCapacityCv(const ThermodynamicReferenceState& state) {
  return heatCapacityCp(state) * 3.0 / 5.0;
}
} // namespace Utils
} // namespace Scine