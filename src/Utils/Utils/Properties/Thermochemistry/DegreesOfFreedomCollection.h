/**
 * @file DegreesOfFreedomCollection.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_DEGREESOFFREEDOMCOLLECTION_H
#define UTILSOS_DEGREESOFFREEDOMCOLLECTION_H

#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {
class DegreeOfFreedom;
class ThermodynamicReferenceState;

/**
 * @brief This class represents an arbitrary collection of degrees of freedom.
 */
class DegreesOfFreedomCollection : public std::vector<std::shared_ptr<DegreeOfFreedom>> {
 public:
  /**
   * @brief Default constructor. Degrees of freedom may be added by calling push_back.
   */
  DegreesOfFreedomCollection();
  /**
   * @brief Default destructor.
   */
  ~DegreesOfFreedomCollection();
  /**
   * @brief Calculate the density of states. If this container is empty, only one state is possible at zero energy.
   * @param energies The energy integration grid.
   * @return The density of states.
   */
  Eigen::VectorXd densityOfState(const Eigen::VectorXd& energies);
  /**
   * @brief Calculate the sum of states. If this container is empty, only one state is possible at zero energy.
   * @param energies The energy integration grid.
   * @return The sum of states.
   */
  Eigen::VectorXd sumOfStates(const Eigen::VectorXd& energies);
  /**
   * @brief Calculate the partition function of the system. If this container is empty, the partition function is 1.
   * @param state The reference state (temperature, pressure).
   * @return The partition function.
   */
  double partitionFunction(const ThermodynamicReferenceState& state);
  /**
   * @brief Calculate the entropy of the system.
   * @param state The reference state (temperature, pressure).
   * @return The entropy.
   */
  double entropy(const ThermodynamicReferenceState& state);
  /**
   * @brief Calculate the enthalpy of the system.
   * @param state The reference state (temperature, pressure).
   * @return The enthalpy.
   */
  double enthalpy(const ThermodynamicReferenceState& state);
  /**
   * @brief Calculate the heat capacity at constant pressure.
   * @param state The reference state (temperature, pressure).
   * @return The heat capacity.
   */
  double heatCapacityCp(const ThermodynamicReferenceState& state);
  /**
   * @brief Calculate the heat capacity at constant volume.
   * @param state The reference state (temperature, pressure).
   * @return The heat capacity.
   */
  double heatCapacityCv(const ThermodynamicReferenceState& state);
  /**
   * @brief Calculate the Helmholtz free energy.
   * @param state The reference state (temperature, pressure).
   * @return The Helmholtz free energy.
   */
  double helmholtzFreeEnergy(const ThermodynamicReferenceState& state);

 private:
  /**
   * @brief Calculate the convolution of two distributions \$f a(e) * b(e) = \int_0^e a(x) b(e - x) dx \$f.
   * The distributions must be provided as vectors where each entry corresponds to an energy on a grid. The energy
   * grids are assumed to be identical.
   * @param a The fist distribution.
   * @param b The second distribution.
   * @return The convolution.
   */
  static Eigen::VectorXd convolution(const Eigen::VectorXd& a, const Eigen::VectorXd& b);
};

} // namespace Utils
} // namespace Scine
#endif // UTILSOS_DEGREESOFFREEDOMCOLLECTION_H
