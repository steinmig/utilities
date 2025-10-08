/**
 * @file StaticRotor.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RigidRotor.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/Utilities/Properties.h>
#include <iostream>

Scine::Utils::RigidRotor::~RigidRotor() = default;

Scine::Utils::RigidRotor::RigidRotor(std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> principalMomentsOfInertia,
                                     bool linear, bool classicalMechanics, unsigned int symmetryNumber)
  : principalMomentsOfInertia_(principalMomentsOfInertia),
    linear_(linear),
    classicalMechanics_(classicalMechanics),
    symmetryNumber_(symmetryNumber),
    singleAtom_(!bool((principalMomentsOfInertia->eigenvalues.array() > 0).count())) {
  if (!linear_ && !classicalMechanics_) {
    throw std::runtime_error("Quantum rotors are only available for linear molecules.");
  }
}
Scine::Utils::RigidRotor::RigidRotor(const AtomCollection& atoms, bool classicalMechanics, unsigned int symmetryNumber)
  : linear_(atoms.isLinear()),
    classicalMechanics_(classicalMechanics),
    symmetryNumber_(symmetryNumber),
    singleAtom_(atoms.size() == 1) {
  if (!linear_ && !classicalMechanics_) {
    throw std::runtime_error("Quantum rotors are only available for linear molecules.");
  }
  auto masses = Utils::Geometry::Properties::getMasses(atoms.getElements());
  auto centerOfMass = Utils::Geometry::Properties::getCenterOfMass(atoms.getPositions(), masses);
  principalMomentsOfInertia_ = std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(
      Utils::Geometry::Properties::calculatePrincipalMoments(atoms.getPositions(), masses, centerOfMass));
  if (atoms.size() == 2 && atoms.getElement(0) == atoms.getElement(1)) {
    symmetryNumber_ = 2;
  }
}
Eigen::VectorXd Scine::Utils::RigidRotor::densityOfState(const Eigen::VectorXd& energies) {
  if (energies.size() == 0) {
    return Eigen::VectorXd::Zero(0);
    ;
  }
  if (singleAtom_) {
    Eigen::VectorXd densityOfStates = Eigen::VectorXd::Zero(energies.size());
    densityOfStates[0] = 1.0;
    return densityOfStates;
  }
  assertNonNegativeEnergies(energies);
  if (linear_) {
    if (classicalMechanics_) {
      const double deltaE = (energies.size() == 1) ? 1.0 : energies[1] - energies[0];
      const double b = this->getRotationalConstants()[2];
      return Eigen::VectorXd::Constant(energies.size(), deltaE / b / symmetryNumber_);
    }
    else {
      return this->calculateQuantumDensityOfStates(energies) / symmetryNumber_;
    }
  }
  const double deltaE = (energies.size() == 1) ? 1.0 : energies[1] - energies[0];
  double prod = getRotationalConstants().prod();
  Eigen::VectorXd densityOfStates = (2.0 * (energies.array() / prod).sqrt() * deltaE / symmetryNumber_).matrix();
  return densityOfStates;
}
Eigen::VectorXd Scine::Utils::RigidRotor::sumOfStates(const Eigen::VectorXd& energies) {
  if (energies.size() == 0) {
    return Eigen::VectorXd::Zero(0);
    ;
  }
  assertNonNegativeEnergies(energies);
  if (singleAtom_) {
    Eigen::VectorXd sumOfStates = Eigen::VectorXd::Ones(energies.size());
    return sumOfStates;
  }
  if (linear_) {
    if (classicalMechanics_) {
      const double b = this->getRotationalConstants()[2];
      return (energies.array() / b / symmetryNumber_).matrix();
    }
    else {
      return this->calculateQuantumSumOfStates(energies) / symmetryNumber_;
    }
  }
  double prod = getRotationalConstants().prod();
  Eigen::VectorXd sumOfStates = (4.0 / 3.0 * energies.array() * (energies.array() / prod).sqrt() / symmetryNumber_).matrix();
  return sumOfStates;
}
Eigen::Vector3d Scine::Utils::RigidRotor::getRotationalConstants() {
  if (!rotationalConstants_) {
    /*
     * B = h^2 / (8 pi^2 c I)
     * with h / (2 pi) = 1
     * B = 1/(2 c I)
     */
    constexpr double rotationConstantPrefactor = 1.0 / (2 * Constants::electronRestMass_per_u);
    Eigen::Vector3d b =
        (this->principalMomentsOfInertia_->eigenvalues.array().inverse() * rotationConstantPrefactor).matrix();
    rotationalConstants_ = std::make_unique<Eigen::Vector3d>(b);
  }
  return *rotationalConstants_;
}
double Scine::Utils::RigidRotor::partitionFunction(const Scine::Utils::ThermodynamicReferenceState& state) {
  if (state.temperature < 1e-6 || singleAtom_) {
    return 1.0;
  }
  if (linear_) {
    if (classicalMechanics_) {
      const double b = getRotationalConstants()[2];
      return kB * state.temperature / b / symmetryNumber_;
    }
    else {
      return this->calculateQuantumPartitioningFunction(state.temperature) / symmetryNumber_;
    }
  }
  double prod = getRotationalConstants().prod();
  double kbTCubed = state.temperature * state.temperature * state.temperature * kB * kB * kB;
  return std::sqrt(Constants::pi * kbTCubed / prod) / symmetryNumber_;
}
unsigned int Scine::Utils::RigidRotor::degeneracy(unsigned int j) {
  if (!linear_) {
    throw std::runtime_error("Rigid rotor state enumeration is only implemented for linear rotors.");
  }
  return 2 * j + 1;
}
double Scine::Utils::RigidRotor::energyLevel(unsigned int j) {
  if (!linear_) {
    throw std::runtime_error("Rigid rotor state enumeration is only implemented for linear rotors.");
  }
  const double b = getRotationalConstants()[2];
  return b * double(j * (j + 1));
}
double Scine::Utils::RigidRotor::entropy(const Scine::Utils::ThermodynamicReferenceState& state) {
  if (singleAtom_) {
    return 0.0;
  }
  if (linear_ && !classicalMechanics_) {
    return calculateQuantumEntropy(state.temperature) - std::log(symmetryNumber_);
  }
  double q = partitionFunction(state);
  double directions = (linear_) ? 1.0 : 1.5;
  return kB * (std::log(q) + directions);
}
double Scine::Utils::RigidRotor::enthalpy(const Scine::Utils::ThermodynamicReferenceState& state) {
  if (singleAtom_) {
    return 0.0;
  }
  if (linear_ && !classicalMechanics_) {
    return calculateQuantumEnthalpy(state.temperature);
  }
  double prefactor = (linear_) ? 1.0 : 1.5;
  return prefactor * state.temperature * kB;
}
double Scine::Utils::RigidRotor::heatCapacityCp(const Scine::Utils::ThermodynamicReferenceState&) {
  if (singleAtom_) {
    return 0.0;
  }
  if (!classicalMechanics_) {
    throw std::runtime_error("There is no implementation for the heat capacity of a QM rotor available at the moment."
                             " Please use the classical implementation.");
  }
  return (linear_) ? kB : 1.5 * kB;
}
double Scine::Utils::RigidRotor::heatCapacityCv(const Scine::Utils::ThermodynamicReferenceState& state) {
  return heatCapacityCp(state) * 3.0 / 5.0;
}
unsigned int Scine::Utils::RigidRotor::getSymmetryNumber() {
  return symmetryNumber_;
}
void Scine::Utils::RigidRotor::setSymmetryNumber(unsigned int symmetryNumber) {
  symmetryNumber_ = symmetryNumber;
}
