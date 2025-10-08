/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/Utilities/Properties.h>
#include <Utils/Properties/Thermochemistry/IdealGasTranslator.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <gmock/gmock.h>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {

class IdealGasTranslatorTest : public Test {
 public:
  AtomCollection hfAtoms;

 protected:
  void SetUp() final {
    PositionCollection hfPositions(2, 3);
    hfPositions << 0.0000000000000, 0.0000000000000, 0.0000000000000, 0.9655884052935, 0.0000000000000, -0.0000001000000;
    hfPositions *= Constants::bohr_per_angstrom;

    hfAtoms = AtomCollection({ElementType::H, ElementType::F}, hfPositions);
  }
};

TEST_F(IdealGasTranslatorTest, DifferentConstructionsAreEqual) {
  double molarMass = 0.0;
  for (const auto& atomMass : Geometry::Properties::getMasses(hfAtoms.getElements())) {
    molarMass += atomMass;
  }
  IdealGasTranslator idealGasTranslatorA(molarMass);
  IdealGasTranslator idealGasTranslatorB(hfAtoms);
  const auto state = ReferenceStates::StandardStateGas;

  EXPECT_NEAR(idealGasTranslatorA.partitionFunction(state), idealGasTranslatorB.partitionFunction(state), 1e-9);
  EXPECT_NEAR(idealGasTranslatorA.enthalpy(state), idealGasTranslatorB.enthalpy(state), 1e-12);
  EXPECT_NEAR(idealGasTranslatorA.entropy(state), idealGasTranslatorB.entropy(state), 1e-12);
  EXPECT_NEAR(idealGasTranslatorA.heatCapacityCp(state), idealGasTranslatorB.heatCapacityCp(state), 1e-12);
  EXPECT_NEAR(idealGasTranslatorA.heatCapacityCv(state), idealGasTranslatorB.heatCapacityCv(state), 1e-12);
  EXPECT_NEAR(idealGasTranslatorA.helmholtzFreeEnergy(state), idealGasTranslatorB.helmholtzFreeEnergy(state), 1e-12);

  // Expect throws for density and sum of states.
  // This feature is not implemented yet.
  auto energies = MolecularDegreesOfFreedom::getEnergyGraining(1e-3, 1);
  EXPECT_THROW(idealGasTranslatorA.densityOfState(energies), std::runtime_error);
  EXPECT_THROW(idealGasTranslatorA.sumOfStates(energies), std::runtime_error);
}

TEST_F(IdealGasTranslatorTest, zeroTemperature) {
  IdealGasTranslator idealGasTranslator(hfAtoms);
  const ThermodynamicReferenceState state = ReferenceStates::VacuumZeroKelvin;
  EXPECT_NEAR(idealGasTranslator.enthalpy(state), 0.0, 1e-5);
  EXPECT_NEAR(idealGasTranslator.heatCapacityCp(state), 0.0, 1e-5);
  EXPECT_NEAR(idealGasTranslator.entropy(state), 0.0, 1e-5);
}

TEST_F(IdealGasTranslatorTest, zeroPressure) {
  IdealGasTranslator idealGasTranslator(hfAtoms);
  const ThermodynamicReferenceState state(298.15, 0.0);
  EXPECT_THROW(idealGasTranslator.entropy(state), std::runtime_error);
  EXPECT_NO_THROW(idealGasTranslator.enthalpy(state));
}

TEST_F(IdealGasTranslatorTest, standardStateChange) {
  const double standardStateShift = 7.9219 * Constants::hartree_per_kJPerMol;
  IdealGasTranslator idealGasTranslator(hfAtoms);
  const double freeEnergyGasPhase = idealGasTranslator.helmholtzFreeEnergy(ReferenceStates::StandardStateGas);
  const double freeEnergyLiquidPhase = idealGasTranslator.helmholtzFreeEnergy(ReferenceStates::StandardStateLiquid);
  EXPECT_NEAR(freeEnergyGasPhase + standardStateShift - freeEnergyLiquidPhase, 0.0, standardStateShift * 1e-2);
}

TEST_F(IdealGasTranslatorTest, referenceData) {
  double molarMass = 0.0;
  for (const auto& atomMass :
       Geometry::Properties::getMasses({ElementType::C, ElementType::O, ElementType::H, ElementType::H})) {
    molarMass += atomMass;
  }
  IdealGasTranslator idealGasTranslator(molarMass);
  const ThermodynamicReferenceState state(298.00, 101325.0);

  EXPECT_NEAR(idealGasTranslator.enthalpy(state) * Constants::kCalPerMol_per_hartree, 1480.4688 / 1000, 1e-5);
  EXPECT_NEAR(idealGasTranslator.heatCapacityCp(state) * Constants::kCalPerMol_per_hartree, 4.9680 / 1000, 1e-5);
  EXPECT_NEAR(idealGasTranslator.entropy(state) * Constants::kCalPerMol_per_hartree, 36.1295 / 1000, 1e-5);
}

} // namespace Utils
} // namespace Scine
