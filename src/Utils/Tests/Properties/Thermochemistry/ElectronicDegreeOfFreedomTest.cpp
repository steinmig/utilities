/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/Utilities/Properties.h>
#include <Utils/Properties/Thermochemistry/ElectronicDegreeOfFreedom.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <gmock/gmock.h>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {

class ElectronicDegreeOfFreedomTest : public Test {
 public:
  const double groundStateEnergy = -230.544687652801;
};

TEST_F(ElectronicDegreeOfFreedomTest, onlyGroundStateThermodynamicFunctions) {
  ElectronicDegreeOfFreedom electronicDegreeOfFreedom(1, groundStateEnergy);
  const auto state = ReferenceStates::StandardStateLiquid;
  EXPECT_NEAR(electronicDegreeOfFreedom.enthalpy(state), groundStateEnergy, 1e-12);
  EXPECT_NEAR(electronicDegreeOfFreedom.helmholtzFreeEnergy(state), groundStateEnergy, 1e-12);
  EXPECT_NEAR(electronicDegreeOfFreedom.heatCapacityCp(state), 0.0, 1e-12);
  EXPECT_NEAR(electronicDegreeOfFreedom.heatCapacityCv(state), 0.0, 1e-12);
  EXPECT_NEAR(electronicDegreeOfFreedom.entropy(state), 0.0, 1e-12);

  const auto energies = MolecularDegreesOfFreedom::getEnergyGraining(1e-3, 0.1);
  EXPECT_NEAR(electronicDegreeOfFreedom.densityOfState(energies).sum(), 1.0, 1e-12);
  EXPECT_NEAR(electronicDegreeOfFreedom.sumOfStates(energies).sum(), energies.size(), 1e-12);
}

TEST_F(ElectronicDegreeOfFreedomTest, fakeSpectrum) {
  auto energyLevels = MolecularDegreesOfFreedom::getEnergyGraining(0.1, 1);
  energyLevels.array() += groundStateEnergy;
  ElectronicDegreeOfFreedom electronicDegreeOfFreedom(1, groundStateEnergy, energyLevels);
  const auto zeroKelvin = ReferenceStates::VacuumZeroKelvin;
  EXPECT_NEAR(electronicDegreeOfFreedom.enthalpy(zeroKelvin), groundStateEnergy, 1e-12);
  EXPECT_NEAR(electronicDegreeOfFreedom.helmholtzFreeEnergy(zeroKelvin), groundStateEnergy, 1e-12);
  EXPECT_NEAR(electronicDegreeOfFreedom.entropy(zeroKelvin), 0.0, 1e-12);

  // Check if the expected trends are fulfilled.
  const ThermodynamicReferenceState plasma(298150.0, 1e+5);
  EXPECT_GT(electronicDegreeOfFreedom.enthalpy(plasma), groundStateEnergy);
  EXPECT_LT(electronicDegreeOfFreedom.helmholtzFreeEnergy(plasma), groundStateEnergy);
  EXPECT_GT(electronicDegreeOfFreedom.entropy(plasma), 0.0);

  const auto energies = MolecularDegreesOfFreedom::getEnergyGraining(0.01, 11);
  const auto densityOfStates = electronicDegreeOfFreedom.densityOfState(energies);
  EXPECT_NEAR(densityOfStates.sum(), 10, 1e-12);
  EXPECT_NEAR(densityOfStates[10], 1, 1e-12);
  EXPECT_NEAR(densityOfStates[11], 0, 1e-12);
  EXPECT_NEAR(densityOfStates[9], 0, 1e-12);
  const auto sumOfStates = electronicDegreeOfFreedom.sumOfStates(energies);
  EXPECT_NEAR(sumOfStates[sumOfStates.size() - 1], 10, 1e-12);
  EXPECT_NEAR(sumOfStates[100], 10, 1e-12);
  EXPECT_THROW(electronicDegreeOfFreedom.densityOfState(Eigen::VectorXd::Constant(5, -1.0)), std::runtime_error);
}

TEST_F(ElectronicDegreeOfFreedomTest, emptyEnergyRange) {
  auto energyLevels = MolecularDegreesOfFreedom::getEnergyGraining(0.1, 1);
  energyLevels.array() += groundStateEnergy;
  ElectronicDegreeOfFreedom electronicDegreeOfFreedom(1, groundStateEnergy, energyLevels);
  Eigen::VectorXd energies = Eigen::VectorXd::Zero(0);
  EXPECT_EQ(energies.size(), electronicDegreeOfFreedom.densityOfState(energies).size());
  EXPECT_EQ(energies.size(), electronicDegreeOfFreedom.sumOfStates(energies).size());
}

TEST_F(ElectronicDegreeOfFreedomTest, spin) {
  const auto zeroKelvin = ReferenceStates::VacuumZeroKelvin;
  ElectronicDegreeOfFreedom electronicDegreeOfFreedom(2, groundStateEnergy);
  EXPECT_NEAR(electronicDegreeOfFreedom.entropy(zeroKelvin), kB * std::log(2), 1e-12);
}

} // namespace Utils
} // namespace Scine