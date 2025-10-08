/**
 * @file ElectronicDegreeOfFreedom.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElectronicDegreeOfFreedom.h"
#include <utility>

namespace Scine {
namespace Utils {
ElectronicDegreeOfFreedom::ElectronicDegreeOfFreedom(unsigned int spin, double electronicEnergy, Eigen::VectorXd energyLevels)
  : spin_(spin), electronicEnergy_(electronicEnergy), energyLevels_(std::move(energyLevels)) {
  if (energyLevels_.size() == 0) {
    energyLevels_ = Eigen::VectorXd::Constant(1, electronicEnergy_);
  }
  energyLevels_.array() -= electronicEnergy_;
}
Eigen::VectorXd ElectronicDegreeOfFreedom::densityOfState(const Eigen::VectorXd& energies) {
  return this->calculateQuantumDensityOfStates(energies, 0, energyLevels_.size() - 1);
}
Eigen::VectorXd ElectronicDegreeOfFreedom::sumOfStates(const Eigen::VectorXd& energies) {
  return this->calculateQuantumSumOfStates(energies, 0, energyLevels_.size() - 1);
}
double ElectronicDegreeOfFreedom::partitionFunction(const ThermodynamicReferenceState& state) {
  return this->calculateQuantumPartitioningFunction(state.temperature, 0, energyLevels_.size() - 1);
}
double ElectronicDegreeOfFreedom::entropy(const ThermodynamicReferenceState& state) {
  return this->calculateQuantumEntropy(state.temperature, 0, energyLevels_.size() - 1);
}
double ElectronicDegreeOfFreedom::enthalpy(const ThermodynamicReferenceState& state) {
  return this->calculateQuantumEnthalpy(state.temperature, 0, energyLevels_.size() - 1) + electronicEnergy_;
}
double ElectronicDegreeOfFreedom::heatCapacityCp(const ThermodynamicReferenceState&) {
  if (energyLevels_.size() > 1) {
    throw std::runtime_error(
        "The heat capacity cannot be calculated if excited states are close in energy to the ground state");
  }
  // TODO: This is only true if the excitation energies are well separated from the ground state. Implement the general
  // expression!
  return 0.0;
}
double ElectronicDegreeOfFreedom::heatCapacityCv(const ThermodynamicReferenceState&) {
  if (energyLevels_.size() > 1) {
    throw std::runtime_error(
        "The heat capacity cannot be calculated if excited states are close in energy to the ground state");
  }
  // TODO: This is only true if the excitation energies are well separated from the ground state. Implement the general
  // expression!
  return 0.0;
}
unsigned int ElectronicDegreeOfFreedom::degeneracy(unsigned int) {
  return spin_;
}
double ElectronicDegreeOfFreedom::energyLevel(unsigned int quantumNumber) {
  if (quantumNumber >= energyLevels_.size()) {
    return std::numeric_limits<double>::infinity();
  }
  return energyLevels_(quantumNumber);
}

ElectronicDegreeOfFreedom::~ElectronicDegreeOfFreedom() = default;

} // namespace Utils
} // namespace Scine