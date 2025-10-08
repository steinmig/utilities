/**
 * @file IdealGasTranslator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_IDEALGASTRANSLATOR_H
#define UTILSOS_IDEALGASTRANSLATOR_H

#include "Utils/Properties/Thermochemistry/DegreeOfFreedom.h"
#include <memory>

namespace Scine {
namespace Utils {
class AtomCollection;

/**
 * @class
 * @brief This class represents the translational degrees of freedom of a molecule as an ideal gas.
 */
class IdealGasTranslator : public DegreeOfFreedom {
 public:
  /**
   * @brief Constructor
   * @param atoms The molecule's atoms.
   */
  IdealGasTranslator(const AtomCollection& atoms);
  /**
   * @brief Constructor.
   * @param molarMass The molecule's molar mass in g/mol.
   */
  IdealGasTranslator(double molarMass);
  /**
   * @brief Default destructor.
   */
  ~IdealGasTranslator();

  /**
   * @brief Getter for the density of states at the given energies.
   *
   * @note Not implemented yet since it depends on the box volume/pressure!
   *
   * @param energies The energies.
   * @return The density of states.
   */
  virtual Eigen::VectorXd densityOfState(const Eigen::VectorXd& energies) override;
  /**
   * @brief Getter for the sum of states at the given energies.
   *
   * @note Not implemented yet since it depends on the box volume/pressure!
   *
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
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCv(const ThermodynamicReferenceState& state) override;

 private:
  double moleculeMass_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILSOS_IDEALGASTRANSLATOR_H
