/**
 * @file MolecularDegreesOfFreedom.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularDegreesOfFreedom.h"
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/Utilities/Properties.h>
#include <Utils/Properties/Thermochemistry/ElectronicDegreeOfFreedom.h>
#include <Utils/Properties/Thermochemistry/HarmonicOscillator.h>
#include <Utils/Properties/Thermochemistry/IdealGasTranslator.h>
#include <Utils/Properties/Thermochemistry/RigidRotor.h>
#include <iostream>

namespace Scine {
namespace Utils {

MolecularDegreesOfFreedom::MolecularDegreesOfFreedom(const HessianMatrix& hessian, const AtomCollection& atoms,
                                                     int spinMultiplicity, double electronicEnergy,
                                                     unsigned int symmetryNumber, bool classicalVibrations,
                                                     const Eigen::VectorXd& allElectronicEnergies) {
  elec_ = std::make_shared<ElectronicDegreeOfFreedom>(spinMultiplicity, electronicEnergy, allElectronicEnergies);
  trans_ = std::make_shared<IdealGasTranslator>(atoms);
  rot_ = std::make_shared<RigidRotor>(atoms, true, symmetryNumber);
  vib_ = std::make_shared<HarmonicOscillator>(hessian, atoms, classicalVibrations);
  this->clear();
  this->push_back(elec_);
  this->push_back(trans_);
  this->push_back(rot_);
  this->push_back(vib_);
}
MolecularDegreesOfFreedom::MolecularDegreesOfFreedom(
    std::shared_ptr<NormalModesContainer> normalModesContainer,
    std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> principalMomentsOfInertia,
    const ElementTypeCollection& elements, int spinMultiplicity, double electronicEnergy, unsigned int symmetryNumber,
    bool classicalVibrations, const Eigen::VectorXd& allElectronicEnergies) {
  double molarMass = 0.0;
  for (const auto& atomMass : Geometry::Properties::getMasses(elements)) {
    molarMass += atomMass;
  }
  bool isLinear = normalModesContainer->getWaveNumbers().size() == 3 * elements.size() - 5;
  elec_ = std::make_shared<ElectronicDegreeOfFreedom>(spinMultiplicity, electronicEnergy, allElectronicEnergies);
  trans_ = std::make_shared<IdealGasTranslator>(molarMass);
  rot_ = std::make_shared<RigidRotor>(principalMomentsOfInertia, isLinear, true, symmetryNumber);
  vib_ = std::make_shared<HarmonicOscillator>(normalModesContainer, classicalVibrations);
  this->clear();
  this->push_back(elec_);
  this->push_back(trans_);
  this->push_back(rot_);
  this->push_back(vib_);
}
std::shared_ptr<DegreeOfFreedom> MolecularDegreesOfFreedom::getVibrationalDegreesOfFreedom() {
  return vib_;
}
std::shared_ptr<DegreeOfFreedom> MolecularDegreesOfFreedom::getRotationalDegreesOfFreedom() {
  return rot_;
}
std::shared_ptr<DegreeOfFreedom> MolecularDegreesOfFreedom::getTranslationalDegreesOfFreedom() {
  return trans_;
}
std::shared_ptr<DegreeOfFreedom> MolecularDegreesOfFreedom::getElectronicDegreesOfFreedom() {
  return elec_;
}
double MolecularDegreesOfFreedom::getVibrationalZeroPointEnergy() {
  return this->getVibrationalDegreesOfFreedom()->enthalpy(ReferenceStates::VacuumZeroKelvin);
}
unsigned int MolecularDegreesOfFreedom::getSymmetryNumber() {
  return rot_->getSymmetryNumber();
}
void MolecularDegreesOfFreedom::setSymmetryNumber(unsigned int symmetryNumber) {
  rot_->setSymmetryNumber(symmetryNumber);
}
std::shared_ptr<DegreesOfFreedomCollection> MolecularDegreesOfFreedom::getRRKMDegreesOfFreedom(bool includeRotation) {
  auto degreesOfFreedomCollection = (includeRotation) ? rrkmDegreesOfFreedomWithRotation_ : rrkmDegreesOfFreedomWithoutRotation_;
  if (!degreesOfFreedomCollection) {
    degreesOfFreedomCollection = std::make_shared<DegreesOfFreedomCollection>();
    degreesOfFreedomCollection->push_back(vib_->getClassicalVariant());
    if (includeRotation) {
      degreesOfFreedomCollection->push_back(getRotationalDegreesOfFreedom());
      rrkmDegreesOfFreedomWithRotation_ = degreesOfFreedomCollection;
    }
    else {
      rrkmDegreesOfFreedomWithoutRotation_ = degreesOfFreedomCollection;
    }
  }
  return degreesOfFreedomCollection;
}
Eigen::VectorXd MolecularDegreesOfFreedom::getEnergyGraining(double deltaE, double maxEnergy) {
  return getEnergyGraining(deltaE, 0.0, maxEnergy);
}
Eigen::VectorXd MolecularDegreesOfFreedom::getEnergyGraining(double deltaE, double startEnergy, double endEnergy) {
  if (endEnergy < startEnergy) {
    throw std::runtime_error("The start energy must be smaller than the end energy.");
  }
  if (deltaE <= 0) {
    throw std::runtime_error("You cannot use zero or a negative value as an energy increment.");
  }
  const double difference = endEnergy - startEnergy;
  const unsigned int nGrains = std::ceil(difference / deltaE);
  Eigen::VectorXd energies = Eigen::VectorXd::Zero(nGrains);
  double currentEnergy = startEnergy;
  for (unsigned int i = 0; i < nGrains; ++i) {
    energies[i] = currentEnergy;
    currentEnergy += deltaE;
  }
  return energies;
}
const NormalModesContainer& MolecularDegreesOfFreedom::getNormalModesContainer() {
  return this->vib_->getNormalModesContainer();
}

} // namespace Utils
} // namespace Scine