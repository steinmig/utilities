/**
 * @file ThermochemistryContainer.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_THERMOCHEMISTRYCALCULATOR_H
#define UTILS_THERMOCHEMISTRYCALCULATOR_H

#include <Utils/DataStructures/PartialHessian.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry.h>
#include <Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h>

namespace Scine {
namespace Utils {
class MolecularDegreesOfFreedom;

/**
 * @brief Struct containing the single thermochemical properties of interest.
 */
struct ThermochemicalContainer {
  double entropy = 0, enthalpy = 0, heatCapacityP = 0, heatCapacityV = 0, gibbsFreeEnergy = 0,
         zeroPointVibrationalEnergy = 0, partitionFunction = 1.0;
  int symmetryNumber = 1;
  ThermochemicalContainer operator+(const ThermochemicalContainer& rhs) const {
    ThermochemicalContainer result = *this;
    result.symmetryNumber = this->symmetryNumber;
    result.entropy += rhs.entropy;
    result.enthalpy += rhs.enthalpy;
    result.heatCapacityP += rhs.heatCapacityP;
    result.heatCapacityV += rhs.heatCapacityV;
    result.gibbsFreeEnergy += rhs.gibbsFreeEnergy;
    result.zeroPointVibrationalEnergy += rhs.zeroPointVibrationalEnergy;
    result.partitionFunction *= rhs.partitionFunction;
    return result;
  }
  /** @brief Overloaded multiplication for unit conversion. */
  ThermochemicalContainer operator*(double factor) const {
    ThermochemicalContainer result = *this;
    result.entropy *= factor;
    result.enthalpy *= factor;
    result.heatCapacityP *= factor;
    result.heatCapacityV *= factor;
    result.gibbsFreeEnergy *= factor;
    result.zeroPointVibrationalEnergy *= factor;
    return result;
  }
  bool isApprox(const ThermochemicalContainer& rhs, double epsilon = 1e-12) const {
    return (std::fabs(this->enthalpy - rhs.enthalpy) < epsilon && std::fabs(this->entropy - rhs.entropy) < epsilon &&
            std::fabs(this->gibbsFreeEnergy - rhs.gibbsFreeEnergy) < epsilon &&
            std::fabs(this->zeroPointVibrationalEnergy - rhs.zeroPointVibrationalEnergy) < epsilon &&
            std::fabs(this->heatCapacityP - rhs.heatCapacityP) < epsilon &&
            std::fabs(this->heatCapacityV - rhs.heatCapacityV) < epsilon);
  }
};

/**
 * @brief Struct containing the vibrational, rotational, translational and the overall thermochemical components.
 */
struct ThermochemicalComponentsContainer {
  ThermochemicalContainer vibrationalComponent{}, rotationalComponent{}, translationalComponent{},
      electronicComponent{}, overall{};
  /** @brief Overloaded multiplication for unit conversion */
  ThermochemicalComponentsContainer operator*(double factor) const {
    ThermochemicalComponentsContainer result;
    result.vibrationalComponent = this->vibrationalComponent * factor;
    result.rotationalComponent = this->rotationalComponent * factor;
    result.translationalComponent = this->translationalComponent * factor;
    result.electronicComponent = this->electronicComponent * factor;
    result.overall = this->overall * factor;
    return result;
  }
  bool isApprox(const ThermochemicalComponentsContainer& rhs, double epsilon = 1e-12) const {
    return (this->overall.isApprox(rhs.overall, epsilon) &&
            this->electronicComponent.isApprox(rhs.electronicComponent, epsilon) &&
            this->rotationalComponent.isApprox(rhs.rotationalComponent, epsilon) &&
            this->translationalComponent.isApprox(rhs.translationalComponent, epsilon) &&
            this->vibrationalComponent.isApprox(rhs.vibrationalComponent, epsilon));
  }
};

/**
 * In NDDO semiempirical methods, the ZPVE is already included in the electronic energy
 * due to the way they are parametrized. In the thermochemical calculation it must then
 * not be included again.
 */
enum class ZPVEInclusion { alreadyIncluded, notIncluded };

/**
 * @class ThermochemistryCalculator @file ThermochemistryCalculator.h
 * @brief This class calculates and stores the most important thermochemical data.
 * The class calculates important thermochemical descriptors, as the Zero Point Vibrational Energy,
 * the standard enthalpy of formation, the system entropy, the heat capacity and the Gibbs' free energy
 * from the results of a vibrational analysis.
 *
 * By default a temperature of 298.15 K, a pressure of 1 atm, and a molecular symmetry number of one are assumed.
 */
class ThermochemistryCalculator {
 public:
  ThermochemistryCalculator(NormalModesContainer normalModesContainer,
                            Geometry::Properties::PrincipalMomentsOfInertia principalMomentsOfInertia,
                            ElementTypeCollection elements, int spinMultiplicity, double electronicEnergy);
  ThermochemistryCalculator(const HessianMatrix& hessian, const AtomCollection& atoms, int spinMultiplicity,
                            double electronicEnergy);
  ThermochemistryCalculator(const PartialHessian& hessian, const AtomCollection& atoms, int spinMultiplicity,
                            double electronicEnergy);
  ThermochemistryCalculator(const HessianMatrix& hessian, ElementTypeCollection elements,
                            const PositionCollection& positions, int spinMultiplicity, double electronicEnergy);
  ThermochemistryCalculator(const PartialHessian& hessian, ElementTypeCollection elements,
                            const PositionCollection& positions, int spinMultiplicity, double electronicEnergy);
  ~ThermochemistryCalculator() = default;

  /**
   * @brief Setter for the temperature at which the thermochemical calculation is performed.
   */
  void setTemperature(double temperature);
  /**
   * @brief Setter for the pressure at which the thermochemical calculation is performed.
   */
  void setPressure(double pressure);
  /**
   * In NDDO semiempirical methods, the ZPVE is already included in the electronic energy
   * due to the way they are parametrized. In the thermochemical calculation it must then
   * not be included again.
   * @param inclusion notIncluded is the standard way, alreadyIncluded is in case of the NDDO semiempirical methods.
   */
  void setZPVEInclusion(ZPVEInclusion inclusion);
  /**
   * @brief Runt the calculation.
   * @return The thermochemical information in a container object.
   */
  ThermochemicalComponentsContainer calculate();

  /**
   * @brief Sets the symmetry sigma factor.
   * @param sigma The symmetry factor related to the point group symmetry of the molecule.
   * Examples:
   * Cn,v/h : n
   * Dn,v/h : 2*n
   * C_inf,v : 1
   * C_inf,h : 2
   * S_n : n/2
   * T : 12
   * O : 24
   * I : 60
   */
  void setMolecularSymmetryNumber(int sigma);
  /**
   * @brief Getter for the uncerlying molecular degrees of freedom.
   * @return The degrees of freedom.
   */
  const MolecularDegreesOfFreedom& getMolecularDegreesOfFreedom() const;

 protected:
  ZPVEInclusion zpveIncluded{ZPVEInclusion::notIncluded};

 private:
  ThermodynamicReferenceState referenceState_ = ReferenceStates::StandardStateGas;
  std::shared_ptr<MolecularDegreesOfFreedom> molecularDegreesOfFreedom_;
  // Sets the symmetry number for diatomic molecules
  // TODO Automatically determine symmetry for all molecules
  static unsigned int calculateSigmaForDiatomicMolecule(const ElementTypeCollection& elements);
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_THERMOCHEMISTRYCALCULATOR_H
