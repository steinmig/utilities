/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ThermochemistryCalculator.h"
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/Properties/Thermochemistry/DegreeOfFreedom.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>

namespace Scine {
namespace Utils {

ThermochemistryCalculator::ThermochemistryCalculator(NormalModesContainer normalModesContainer,
                                                     Geometry::Properties::PrincipalMomentsOfInertia principalMomentsOfInertia,
                                                     ElementTypeCollection elements, int spinMultiplicity,
                                                     double electronicEnergy)
  : molecularDegreesOfFreedom_(std::make_shared<MolecularDegreesOfFreedom>(
        std::make_shared<NormalModesContainer>(normalModesContainer),
        std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(principalMomentsOfInertia), elements,
        spinMultiplicity, electronicEnergy, calculateSigmaForDiatomicMolecule(elements))) {
}

ThermochemistryCalculator::ThermochemistryCalculator(const HessianMatrix& hessian, const AtomCollection& atoms,
                                                     int spinMultiplicity, double electronicEnergy)
  : ThermochemistryCalculator(hessian, atoms.getElements(), atoms.getPositions(), spinMultiplicity, electronicEnergy) {
}

ThermochemistryCalculator::ThermochemistryCalculator(const PartialHessian& hessian, const AtomCollection& atoms,
                                                     int spinMultiplicity, double electronicEnergy)
  : ThermochemistryCalculator(hessian, atoms.getElements(), atoms.getPositions(), spinMultiplicity, electronicEnergy) {
}

ThermochemistryCalculator::ThermochemistryCalculator(const HessianMatrix& hessian, ElementTypeCollection elements,
                                                     const PositionCollection& positions, int spinMultiplicity,
                                                     double electronicEnergy)
  : molecularDegreesOfFreedom_(std::make_shared<MolecularDegreesOfFreedom>(hessian, AtomCollection(elements, positions),
                                                                           spinMultiplicity, electronicEnergy,
                                                                           calculateSigmaForDiatomicMolecule(elements))) {
}

ThermochemistryCalculator::ThermochemistryCalculator(const PartialHessian& hessian, ElementTypeCollection elements,
                                                     const PositionCollection& positions, int spinMultiplicity,
                                                     double electronicEnergy) {
  auto masses = Utils::Geometry::Properties::getMasses(elements);
  auto centerOfMass = Utils::Geometry::Properties::getCenterOfMass(positions, masses);
  auto principalMomentsOfInertia = Utils::Geometry::Properties::calculatePrincipalMoments(positions, masses, centerOfMass);
  auto normalModesContainer = Utils::NormalModeAnalysis::calculateNormalModes(hessian, elements, positions);
  molecularDegreesOfFreedom_ = std::make_shared<MolecularDegreesOfFreedom>(
      std::make_shared<NormalModesContainer>(normalModesContainer),
      std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(principalMomentsOfInertia), elements,
      spinMultiplicity, electronicEnergy, calculateSigmaForDiatomicMolecule(elements));
}

void ThermochemistryCalculator::setTemperature(double temperature) {
  if (temperature < 0.0) {
    throw std::runtime_error("A negative temperature was detected.");
  }
  referenceState_.temperature = temperature;
}

void ThermochemistryCalculator::setPressure(double pressure) {
  referenceState_.pressure = pressure;
}

void ThermochemistryCalculator::setZPVEInclusion(ZPVEInclusion inclusion) {
  zpveIncluded = inclusion;
}

ThermochemicalComponentsContainer ThermochemistryCalculator::calculate() {
  ThermochemicalComponentsContainer container;
  container.vibrationalComponent = molecularDegreesOfFreedom_->getVibrationalDegreesOfFreedom()->createThermochemicalContainer(
      referenceState_, zpveIncluded == ZPVEInclusion::alreadyIncluded);
  container.rotationalComponent =
      molecularDegreesOfFreedom_->getRotationalDegreesOfFreedom()->createThermochemicalContainer(referenceState_);
  container.translationalComponent =
      molecularDegreesOfFreedom_->getTranslationalDegreesOfFreedom()->createThermochemicalContainer(referenceState_);
  container.electronicComponent =
      molecularDegreesOfFreedom_->getElectronicDegreesOfFreedom()->createThermochemicalContainer(referenceState_);
  container.vibrationalComponent.zeroPointVibrationalEnergy = molecularDegreesOfFreedom_->getVibrationalZeroPointEnergy();
  container.overall = container.vibrationalComponent + container.rotationalComponent +
                      container.translationalComponent + container.electronicComponent;
  container.overall.symmetryNumber = int(molecularDegreesOfFreedom_->getSymmetryNumber());
  return container;
}

void ThermochemistryCalculator::setMolecularSymmetryNumber(int sigma) {
  molecularDegreesOfFreedom_->setSymmetryNumber(sigma);
}

unsigned int ThermochemistryCalculator::calculateSigmaForDiatomicMolecule(const ElementTypeCollection& elements) {
  if (elements.size() == 2 && elements[1] == elements[0]) {
    return 2;
  }
  return 1;
}
const MolecularDegreesOfFreedom& ThermochemistryCalculator::getMolecularDegreesOfFreedom() const {
  return *molecularDegreesOfFreedom_;
}

} // namespace Utils
} // namespace Scine
