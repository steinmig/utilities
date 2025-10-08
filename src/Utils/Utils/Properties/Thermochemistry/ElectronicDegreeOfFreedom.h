/**
 * @file ElectronicDegreeOfFreedom.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_ELECTRONICDEGREEOFFREEDOM_H
#define UTILSOS_ELECTRONICDEGREEOFFREEDOM_H

#include "Utils/Properties/Thermochemistry/DegreeOfFreedom.h"

namespace Scine {
namespace Utils {

/**
 * @class
 * @brief This class represent the electronic degrees of freedom.
 */
class ElectronicDegreeOfFreedom : public DegreeOfFreedom {
 public:
  /**
   * @brief Constructor.
   * @param spin              The spin multiplicity of the electronic states.
   * @param electronicEnergy  The energy of the ground state.
   * @param energyLevels      The energy levels of ground and excited states. If not provided, only the ground state is
   * used.
   */
  ElectronicDegreeOfFreedom(unsigned int spin, double electronicEnergy,
                            Eigen::VectorXd energyLevels = Eigen::VectorXd::Zero(0));
  /**
   * @brief Default destructor.
   */
  ~ElectronicDegreeOfFreedom();

  /**
   * @brief Getter for the density of states at the given energies.
   * @param energies The energies.
   * @return The density of states.
   */
  virtual Eigen::VectorXd densityOfState(const Eigen::VectorXd& energies) override;
  /**
   * @brief Getter for the sum of states at the given energies.
   * @param energies The energies.
   * @return The sum of states.
   */
  virtual Eigen::VectorXd sumOfStates(const Eigen::VectorXd& energies) override;
  /**
   * @brief Getter for the partition function.
   * @param state The thermodynamic reference state.
   * @return The partition function (in atomic units).
   */
  virtual double partitionFunction(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the entropy.
   * @param state The thermodynamic reference state.
   * @return The entropy (in atomic units).
   */
  virtual double entropy(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the enthalpy.
   * @param state The thermodynamic reference state.
   * @return The enthalpy (in atomic units).
   */
  virtual double enthalpy(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the heat capacity at constant pressure.
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCp(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the heat capacity at constant volume.
   *
   * @note Only implemented for the case that only the ground state is available.
   *
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCv(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the degeneracy of the state with the given quantum number.
   *
   * @note Only implemented for the case that only the ground state is available.
   *
   * @param quantumNumber The quantum number.
   * @return The degeneracy.
   */
  virtual unsigned int degeneracy(unsigned int quantumNumber) override;
  /**
   * @brief Getter for the energy of the state with the given quantum number.
   *        May return inf if the quantum number exceeds the electronic states.
   * @param quantumNumber The quantum number.
   * @return The energy.
   */
  virtual double energyLevel(unsigned int quantumNumber) override;

 private:
  const unsigned int spin_;
  const double electronicEnergy_;
  Eigen::VectorXd energyLevels_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILSOS_ELECTRONICDEGREEOFFREEDOM_H
