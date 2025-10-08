/**
 * @file MolecularDegreesOfFreedom.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_MOLECULARDEGREESOFFREEDOM_H
#define UTILSOS_MOLECULARDEGREESOFFREEDOM_H

#include <Utils/Properties/Thermochemistry/DegreesOfFreedomCollection.h>
#include <Utils/Typenames.h>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {
class DegreeOfFreedom;
class NormalModesContainer;
namespace Geometry::Properties {
struct PrincipalMomentsOfInertia;
}
class AtomCollection;
class ThermodynamicReferenceState;
class RigidRotor;
class HarmonicOscillator;
class ElectronicDegreeOfFreedom;
class IdealGasTranslator;

/**
 * @class
 * @brief This class represents all the degrees of freedom of a molecule. At the moment it relies on the rigid
 * rotor/particle in a box/ harmonic oscillator model.
 */
class MolecularDegreesOfFreedom : public std::enable_shared_from_this<MolecularDegreesOfFreedom>,
                                  public DegreesOfFreedomCollection {
 public:
  /**
   * @brief Constructor.
   * @param normalModesContainer            The molecule's vibrational modes.
   * @param principalMomentsOfInertia       The molecule's principal moment of inertia.
   * @param elements                        The elements.
   * @param spinMultiplicity                The spin multiplicity.
   * @param electronicEnergy                The ground state electronic energy.
   * @param symmetryNumber                  The molecule's symmetry number, i.e., the factor of symmetry redundant
   * rotational states.
   * @param classicalMechanicVibrations     If true, classical vibrations are used. By default false.
   * @param allElectronicEnergies           Electronic energies (including the ground state). If none are provided, only
   * the ground state is used.
   */
  MolecularDegreesOfFreedom(std::shared_ptr<NormalModesContainer> normalModesContainer,
                            std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> principalMomentsOfInertia,
                            const ElementTypeCollection& elements, int spinMultiplicity, double electronicEnergy,
                            unsigned int symmetryNumber = 1, bool classicalMechanicVibrations = false,
                            const Eigen::VectorXd& allElectronicEnergies = Eigen::VectorXd::Zero(0));
  /**
   * @brief Constructor
   * @param hessian                      The hessian.
   * @param atoms                        The atom collection.
   * @param spinMultiplicity             The spin multiplicity.
   * @param electronicEnergy             The ground state electronic energy.
   * @param symmetryNumber               The molecule's symmetry number, i.e., the factor of symmetry redundant
   * rotational states.
   * @param classicalMechanicVibrations  If true, classical vibrations are used. By default false.
   * @param allElectronicEnergies        Electronic energies (including the ground state). If none are provided, only
   * the ground state is used.
   */
  MolecularDegreesOfFreedom(const HessianMatrix& hessian, const AtomCollection& atoms, int spinMultiplicity,
                            double electronicEnergy, unsigned int symmetryNumber = 1, bool classicalMechanicVibrations = false,
                            const Eigen::VectorXd& allElectronicEnergies = Eigen::VectorXd::Zero(0));
  /**
   * @brief Getter for an energy range from zero to maxEnergy in steps of deltaE
   * @param deltaE     The energy increment.
   * @param maxEnergy  The max energy.
   * @return The energy list.
   */
  static Eigen::VectorXd getEnergyGraining(double deltaE, double maxEnergy);
  /**
   * @brief Getter for an energy range from startEnergy to endEnergy in steps of deltaE.
   * @param deltaE       The energy increment.
   * @param startEnergy  The minimum energy.
   * @param endEnergy    The max energy.
   * @return
   */
  static Eigen::VectorXd getEnergyGraining(double deltaE, double startEnergy, double endEnergy);

  /**
   * @brief Getter for the vibrational degrees of freedom.
   * @return The vibrational degrees of freedom.
   */
  std::shared_ptr<DegreeOfFreedom> getVibrationalDegreesOfFreedom();
  /**
   * @brief Getter for the rotational degrees of freedom.
   * @return The rotational degrees of freedom.
   */
  std::shared_ptr<DegreeOfFreedom> getRotationalDegreesOfFreedom();
  /**
   * @brief Getter for the translational degrees of freedom.
   * @return The translational degrees of freedom.
   */
  std::shared_ptr<DegreeOfFreedom> getTranslationalDegreesOfFreedom();
  /**
   * @brief Getter for the electronic degrees of freedom.
   * @return The electronic degrees of freedom.
   */
  std::shared_ptr<DegreeOfFreedom> getElectronicDegreesOfFreedom();
  /**
   * @brief Getter for the degrees of freedom used in RRKM theory.
   * @param includeRotation If true, the rotational degrees of freedom are included.
   * @return The degrees of freedom for RRKM theory.
   */
  std::shared_ptr<DegreesOfFreedomCollection> getRRKMDegreesOfFreedom(bool includeRotation = true);
  /**
   * @brief Getter for the zero point energy, i.e, energy in vacuum at zero kelvin.
   * @return The zero point energy.
   */
  double getVibrationalZeroPointEnergy();
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
  /**
   * @brief Getter for the normal modes container.
   * @return The normal modes container.
   */
  const NormalModesContainer& getNormalModesContainer();

 private:
  std::shared_ptr<HarmonicOscillator> vib_;
  std::shared_ptr<RigidRotor> rot_;
  std::shared_ptr<IdealGasTranslator> trans_;
  std::shared_ptr<ElectronicDegreeOfFreedom> elec_;
  std::shared_ptr<DegreesOfFreedomCollection> rrkmDegreesOfFreedomWithRotation_;
  std::shared_ptr<DegreesOfFreedomCollection> rrkmDegreesOfFreedomWithoutRotation_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILSOS_MOLECULARDEGREESOFFREEDOM_H
