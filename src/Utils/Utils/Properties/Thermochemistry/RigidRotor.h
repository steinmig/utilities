/**
 * @file StaticRotor.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_RIGIDROTOR_H
#define UTILSOS_RIGIDROTOR_H

#include "Utils/Properties/Thermochemistry/DegreeOfFreedom.h"
#include "Utils/Typenames.h"
#include <memory>

namespace Scine {
namespace Utils {

namespace Geometry::Properties {
struct PrincipalMomentsOfInertia;
}
class AtomCollection;

/**
 * @class
 * @brief This class represent rotational degrees of freedom for the molecule as a rigid rotor.
 *
 * Classical statistics are available for linear and non-linear molecules. Quantum mechanical
 * statistics are only available for linear rotors.
 */
class RigidRotor : public DegreeOfFreedom {
 public:
  /**
   * @brief Constructor.
   * @param principalMomentsOfInertia  The principle moment of intertia.
   * @param linear                     If true, the molecule is assumed to be linear.
   * @param classicalMechanics         If true, classical mechanics are assumed. By default true.
   * @param symmetryNumber             The symmetry number, i.e., the factor of symmetry redundant microstates.
   */
  RigidRotor(std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> principalMomentsOfInertia, bool linear,
             bool classicalMechanics = true, unsigned int symmetryNumber = 1);
  /**
   * @brief Constructor
   * @param atoms              The molecule's atoms.
   * @param classicalMechanics If true, classical mechanics are assumed. By default true.
   * @param symmetryNumber     The symmetry number, i.e., the factor of symmetry redundant microstates.
   */
  explicit RigidRotor(const AtomCollection& atoms, bool classicalMechanics = true, unsigned int symmetryNumber = 1);
  /**
   * @brief Default destructor.
   */
  ~RigidRotor();

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
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCv(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the degeneracy of the state with the given quantum number.
   *        This function is only implemented for linear rotors.
   * @param quantumNumber The quantum number.
   * @return The degeneracy.
   */
  virtual unsigned int degeneracy(unsigned int quantumNumber) override;
  /**
   * @brief Getter for the energy of the state with the given quantum number.
   *        This function is only implemented for linear rotors.
   * @param quantumNumber The quantum number.
   * @return The energy.
   */
  virtual double energyLevel(unsigned int quantumNumber) override;
  /**
   * @brief Getter for the rotational constants.
   * @return The rotational constants.
   */
  Eigen::Vector3d getRotationalConstants();
  /**
   * @brief Getter for the symmetry number.
   * @return The symmetry number.
   */
  unsigned int getSymmetryNumber();
  /**
   * @brief Setter for the symmetry number.
   * @param symmetryNumber The symmetry number.
   */
  void setSymmetryNumber(unsigned int symmetryNumber);

 private:
  std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> principalMomentsOfInertia_;
  std::unique_ptr<Eigen::Vector3d> rotationalConstants_;
  const bool linear_;
  const bool classicalMechanics_;
  unsigned int symmetryNumber_;
  const bool singleAtom_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILSOS_RIGIDROTOR_H
