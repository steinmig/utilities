/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/Utilities/Properties.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <Utils/Properties/Thermochemistry/RigidRotor.h>
#include <gmock/gmock.h>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {

class RigidRotorTest : public Test {
 public:
  Geometry::Properties::PrincipalMomentsOfInertia formaldehydePMI;
  Geometry::Properties::PrincipalMomentsOfInertia hfPMI;
  Geometry::Properties::PrincipalMomentsOfInertia arPMI;
  AtomCollection hfAtoms;

 protected:
  void SetUp() final {
    PositionCollection hfPositions(2, 3);
    hfPositions << 0.0000000000000, 0.0000000000000, 0.0000000000000, 0.9655884052935, 0.0000000000000, -0.0000001000000;
    hfPositions *= Constants::bohr_per_angstrom;

    hfAtoms = AtomCollection({ElementType::H, ElementType::F}, hfPositions);
    Eigen::Vector3d eigenValues(2.8969, 21.7672, 24.6640);
    // Convert from 1e-40 g cm^2 to amu*bohr^2
    eigenValues *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
    formaldehydePMI.eigenvalues = eigenValues;
    formaldehydePMI.eigenvectors = Eigen::Matrix3d::Random();

    Eigen::Vector3d eigenValuesHF(0.00000000, 1.4818, 1.4818);
    // Convert from 1e-40 g cm^2 to amu*bohr^2
    eigenValuesHF *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
    hfPMI.eigenvalues = eigenValuesHF;
    hfPMI.eigenvectors = Eigen::Matrix3d::Random();

    arPMI.eigenvalues = Eigen::Vector3d::Zero();
    arPMI.eigenvectors = Eigen::MatrixXd::Identity(3, 3);
  }
};

TEST_F(RigidRotorTest, DifferentConstructionsAreEqual) {
  std::vector<bool> classicalMechanics = {true, false};
  for (const auto& isClassicalMechanics : classicalMechanics) {
    RigidRotor rotorA(std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(hfPMI), true, isClassicalMechanics);
    RigidRotor rotorB(hfAtoms, isClassicalMechanics);

    auto referenceState = ReferenceStates::StandardStateGas;
    auto energies = MolecularDegreesOfFreedom::getEnergyGraining(1 * Constants::hartree_per_invCentimeter,
                                                                 60 * Constants::hartree_per_invCentimeter);

    double diff = (rotorA.densityOfState(energies) - rotorB.densityOfState(energies)).array().abs().maxCoeff();
    EXPECT_NEAR(0.0, diff, 1.0);
    diff = (rotorA.sumOfStates(energies) - rotorB.sumOfStates(energies)).array().abs().maxCoeff();
    EXPECT_NEAR(0.0, diff, 1e-3);
    // First excitation can be reached. So the sum of states should be just above 2.
    EXPECT_NEAR(rotorA.entropy(referenceState), rotorB.entropy(referenceState), 1e-6);
    EXPECT_NEAR(rotorA.enthalpy(referenceState), rotorB.enthalpy(referenceState), 1e-6);
    EXPECT_NEAR(rotorA.partitionFunction(referenceState), rotorB.partitionFunction(referenceState), 5e-1);
    if (isClassicalMechanics) {
      EXPECT_NEAR(rotorA.heatCapacityCv(referenceState), rotorB.heatCapacityCv(referenceState), 1e-6);
      EXPECT_NEAR(rotorA.heatCapacityCp(referenceState), rotorB.heatCapacityCp(referenceState), 1e-6);
    }
    else {
      EXPECT_THROW(rotorA.heatCapacityCv(referenceState), std::runtime_error);
      EXPECT_THROW(rotorA.heatCapacityCp(referenceState), std::runtime_error);
    }
    EXPECT_NEAR(rotorA.helmholtzFreeEnergy(referenceState), rotorB.helmholtzFreeEnergy(referenceState), 1e-6);
    EXPECT_THROW(rotorA.densityOfState(Eigen::VectorXd::Constant(5, -1.0)), std::runtime_error);
    EXPECT_THROW(rotorA.sumOfStates(Eigen::VectorXd::Constant(5, -1.0)), std::runtime_error);
  }
  EXPECT_THROW(RigidRotor(std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(hfPMI), false, false),
               std::runtime_error);
}

TEST_F(RigidRotorTest, linearThermodynamicFunctionsQMvsClassical) {
  RigidRotor qmRotor(hfAtoms, false);
  RigidRotor classicalRotor(hfAtoms, true);
  const auto state = ReferenceStates::StandardStateGas;

  // The classical and QM results should be similar at 298.15 K
  EXPECT_NEAR(qmRotor.partitionFunction(state), classicalRotor.partitionFunction(state),
              std::abs(classicalRotor.partitionFunction(state) * 1e-1));
  EXPECT_NEAR(qmRotor.entropy(state), classicalRotor.entropy(state), std::abs(classicalRotor.entropy(state) * 1e-4));
  EXPECT_NEAR(qmRotor.enthalpy(state), classicalRotor.enthalpy(state), std::abs(classicalRotor.enthalpy(state) * 1e-1));
  EXPECT_NEAR(qmRotor.helmholtzFreeEnergy(state), classicalRotor.helmholtzFreeEnergy(state),
              std::abs(classicalRotor.helmholtzFreeEnergy(state) * 1e-1));

  // Comparing the densities of state is not trivial since the QM result will be a step-wise function and
  // the classical result is continuous. However, we can at least check the sum of states and whether the
  // densities of state sum of to the sum of states.
  auto energies = MolecularDegreesOfFreedom::getEnergyGraining(1 * Constants::hartree_per_invCentimeter,
                                                               600 * Constants::hartree_per_invCentimeter);
  const auto qmDensityOfStates = qmRotor.densityOfState(energies);
  const auto classicalDensityOfStates = classicalRotor.densityOfState(energies);
  const auto qmSumOfStates = qmRotor.sumOfStates(energies);
  const auto classicalSumOfStates = classicalRotor.sumOfStates(energies);
  const double highEnergyDiff = (qmSumOfStates.tail(20) - classicalSumOfStates.tail(20)).array().abs().maxCoeff();
  const double orderOfMagnitude = qmSumOfStates[qmSumOfStates.size() - 1];
  EXPECT_NEAR(highEnergyDiff, 0.0, orderOfMagnitude / 5);
  const double trivialQMSum = qmDensityOfStates.sum();
  const double trivialClassicalSum = classicalDensityOfStates.sum();
  EXPECT_NEAR(trivialQMSum, trivialClassicalSum, orderOfMagnitude / 8);
  EXPECT_NEAR(trivialQMSum, orderOfMagnitude, 1e-12);
}
TEST_F(RigidRotorTest, linearMoleculeThermodynamicFunctions) {
  RigidRotor classicalRotor(hfAtoms, true);
  const ThermodynamicReferenceState state(298.00, 1e+5);
  // Reference data for the classical rotor.
  EXPECT_NEAR(classicalRotor.enthalpy(state) * Constants::kCalPerMol_per_hartree, 592.1875 / 1000, 1e-5);
  EXPECT_NEAR(classicalRotor.heatCapacityCp(state) * Constants::kCalPerMol_per_hartree, 1.9872 / 1000, 1e-5);
  EXPECT_NEAR(classicalRotor.entropy(state) * Constants::kCalPerMol_per_hartree, 6.7458 / 1000, 1e-5);
}

TEST_F(RigidRotorTest, nonLinearMoleculeThermodynamicFunctions) {
  RigidRotor rotor(std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(formaldehydePMI), false, true, 2);

  const ThermodynamicReferenceState state(298.00, 1e+5);
  EXPECT_NEAR(rotor.enthalpy(state) * Constants::kCalPerMol_per_hartree, 888.2813 / 1000, 1e-5);
  EXPECT_NEAR(rotor.heatCapacityCp(state) * Constants::kCalPerMol_per_hartree, 2.9808 / 1000, 1e-5);
  EXPECT_NEAR(rotor.entropy(state) * Constants::kCalPerMol_per_hartree, 16.0088 / 1000, 1e-5);

  auto energies = MolecularDegreesOfFreedom::getEnergyGraining(1 * Constants::hartree_per_invCentimeter,
                                                               600 * Constants::hartree_per_invCentimeter);
  const auto densityOfStates = rotor.densityOfState(energies);
  const auto sumOfStates = rotor.sumOfStates(energies);
  const double orderOfMagnitude = sumOfStates[sumOfStates.size() - 1];
  EXPECT_NEAR(densityOfStates.sum(), orderOfMagnitude, orderOfMagnitude * 5e-3);
}

TEST_F(RigidRotorTest, zeroTemperature) {
  RigidRotor rotor(std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(formaldehydePMI), false, true, 2);
  const ThermodynamicReferenceState state = ReferenceStates::VacuumZeroKelvin;
  EXPECT_NEAR(rotor.enthalpy(state), 0.0, 1e-5);
  EXPECT_NEAR(rotor.heatCapacityCp(state), 0.0, 1e-5);
  EXPECT_NEAR(rotor.entropy(state), 0.0, 1e-5);
}

TEST_F(RigidRotorTest, emptyEnergyRange) {
  RigidRotor rotor(std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>(formaldehydePMI), false, true, 2);
  Eigen::VectorXd energies = Eigen::VectorXd::Zero(0);
  EXPECT_EQ(energies.size(), rotor.densityOfState(energies).size());
  EXPECT_EQ(energies.size(), rotor.sumOfStates(energies).size());
}

} // namespace Utils
} // namespace Scine