/**
 * @file DegreesOfFreedomCollection.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DegreesOfFreedomCollection.h"
#include <Utils/Properties/Thermochemistry/DegreeOfFreedom.h>
#include <Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h>

namespace Scine {
namespace Utils {

DegreesOfFreedomCollection::DegreesOfFreedomCollection() = default;
DegreesOfFreedomCollection::~DegreesOfFreedomCollection() = default;

double DegreesOfFreedomCollection::partitionFunction(const ThermodynamicReferenceState& state) {
  double q = 1.0;
  for (const auto& degreeOfFreedom : *this) {
    q *= degreeOfFreedom->partitionFunction(state);
  }
  return q;
}
double DegreesOfFreedomCollection::entropy(const ThermodynamicReferenceState& state) {
  double s = 0.0;
  for (const auto& degreeOfFreedom : *this) {
    s += degreeOfFreedom->entropy(state);
  }
  return s;
}
double DegreesOfFreedomCollection::enthalpy(const ThermodynamicReferenceState& state) {
  double h = 0.0;
  for (const auto& degreeOfFreedom : *this) {
    h += degreeOfFreedom->enthalpy(state);
  }
  return h;
}
double DegreesOfFreedomCollection::heatCapacityCp(const ThermodynamicReferenceState& state) {
  double cp = 0.0;
  for (const auto& degreeOfFreedom : *this) {
    cp += degreeOfFreedom->heatCapacityCp(state);
  }
  return cp;
}
double DegreesOfFreedomCollection::heatCapacityCv(const ThermodynamicReferenceState& state) {
  double cv = 0.0;
  for (const auto& degreeOfFreedom : *this) {
    cv += degreeOfFreedom->heatCapacityCv(state);
  }
  return cv;
}
Eigen::VectorXd DegreesOfFreedomCollection::densityOfState(const Eigen::VectorXd& energies) {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(energies.size());
  if (this->empty()) {
    result[0] = 1.0;
    return result;
  }
  result = this->at(0)->densityOfState(energies);
  for (unsigned int i = 1; i < this->size(); ++i) {
    result = convolution(result, this->at(i)->densityOfState(energies));
  }
  return result;
}
Eigen::VectorXd DegreesOfFreedomCollection::sumOfStates(const Eigen::VectorXd& energies) {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(energies.size());
  if (this->empty()) {
    result = Eigen::VectorXd::Ones(energies.size());
    return result;
  }
  result = this->at(0)->sumOfStates(energies);
  for (unsigned int i = 1; i < this->size(); ++i) {
    /*
     * The sum of states of two independent degrees of freedom is the convolution
     * of the sum of states for degree one and the density of states for the second
     * degree.
     *
     * N(E) = int_0^E dx omega(x)  | omega = density of states, N(E) = sum of states
     * omega(x) = int dy omega_a(y) omega_b(x - y)
     * --> N(E) = int_0^E dx int dy omega_a(y) omega_b(x - y)
     *          = int dy omega_a(y) int_0^E omega_b(x - y) dx
     *          = int dy omega_a(y) N_b(E - y)  | omega_b(z) = 0 for z < 0
     */
    result = convolution(result, this->at(i)->densityOfState(energies));
  }
  return result;
}
Eigen::VectorXd DegreesOfFreedomCollection::convolution(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
  if (a.size() != b.size()) {
    throw std::runtime_error(
        "The energy graining for both distributions must be identical to allow convolution of the distributions.");
  }
  /*
   * a(e) * b(e) = int a(x) b(e - x) dx
   *
   * For a density of states this can be interpreted as follows:
   * For any energy e you assign x to distribution a and the remaining energy
   * e - x to b.
   */
  Eigen::VectorXd result = Eigen::VectorXd::Zero(a.size());
  for (unsigned int e = 0; e < a.size(); ++e) {
    for (unsigned int x = 0; x <= e; ++x) {
      result(e) += a(x) * b(e - x);
    }
  }
  return result;
}
double DegreesOfFreedomCollection::helmholtzFreeEnergy(const ThermodynamicReferenceState& state) {
  return this->enthalpy(state) - state.temperature * this->entropy(state);
}

} // namespace Utils
} // namespace Scine
