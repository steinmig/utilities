/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Constants.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/Utilities/Properties.h>
#include <Utils/Properties/Thermochemistry/HarmonicOscillator.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <gmock/gmock.h>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {

class HarmonicOscillatorTest : public Test {
 public:
  NormalModesContainer formaldehydeNormalModes;
  NormalModesContainer HFNormalModes;
  NormalModesContainer arNormalModes;
  Eigen::MatrixXd hfHessian = Eigen::MatrixXd::Zero(6, 6);
  AtomCollection hfAtoms;

 protected:
  void SetUp() final {
    NormalMode m1(1101.75, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes.add(m1);
    NormalMode m2(1157.94, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes.add(m2);
    NormalMode m3(1349.03, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes.add(m3);
    NormalMode m4(1791.24, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes.add(m4);
    NormalMode m5(2614.79, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes.add(m5);
    NormalMode m6(2664.54, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes.add(m6);
    Eigen::Vector3d eigenValues(2.8969, 21.7672, 24.6640);
    // Convert from 1e-40 g cm^2 to amu*bohr^2
    eigenValues *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
    NormalMode m1HF(3968.7, DisplacementCollection::Random(1, 3));
    HFNormalModes.add(m1HF);
    Eigen::Vector3d eigenValuesHF(0.00000000, 1.4818, 1.4818);
    // Convert from 1e-40 g cm^2 to amu*bohr^2
    eigenValuesHF *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
    // clang-format off
    // ref matrix lower triangular in MILLIDYNES/ANGSTROM/SQRT(MASS(I)*MASS(J))
    std::vector<double> refMatrixTriangular = {
       8.8120622747042,  0.0000000050605,   0.0047603389882,  -0.0000009019966,  -0.0000000000000,   0.0047603389883,
      -2.0296809143334, -0.0000000011656,   0.0000002077567,   0.4674960849783,  -0.0000000011656,  -0.0010964357188,
       0.0000000000000,  0.0000000002685,   0.0002525390079,   0.0000002077567,   0.0000000000000,  -0.0010964357188,
      -0.0000000478526, -0.0000000000000,   0.0002525390079};
    // clang-format on
    int count = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j <= i; ++j) {
        hfHessian(i, j) = refMatrixTriangular[count];
        hfHessian(j, i) = refMatrixTriangular[count];
        count++;
      }
    }
    // Back-scale the mass-weighted coordinates to Cartesian coordinates.
    auto masses = Geometry::Properties::getMasses({ElementType::H, ElementType::F});
    for (unsigned long i = 0; i < masses.size(); ++i) {
      hfHessian.middleRows(3 * i, 3) *= std::sqrt(masses[i]);
      hfHessian.middleCols(3 * i, 3) *= std::sqrt(masses[i]);
    }
    // unit conversion, original is milliDyn / angstrom
    hfHessian *= 1e-8;                                      // N / Angstrom = kg m / (s^2 angstrom)
    hfHessian *= Constants::meter_per_angstrom;             // J / angstrom^2
    hfHessian *= Constants::hartree_per_joule;              // hartree / angstrom^2
    hfHessian *= std::pow(Constants::angstrom_per_bohr, 2); // hartree / bohr^2

    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(3, 3);
    NormalMode m1ar(0.0, identity.row(0));
    NormalMode m2ar(0.0, identity.row(1));
    NormalMode m3ar(0.0, identity.row(2));
    arNormalModes.add(m1ar);
    arNormalModes.add(m2ar);
    arNormalModes.add(m3ar);

    PositionCollection hfPositions(2, 3);
    hfPositions << 0.0000000000000, 0.0000000000000, 0.0000000000000, 0.9655884052935, 0.0000000000000, -0.0000001000000;
    hfPositions *= Constants::bohr_per_angstrom;

    hfAtoms = AtomCollection({ElementType::H, ElementType::F}, hfPositions);
  }
};

TEST_F(HarmonicOscillatorTest, DifferentConstructionsAreEqual) {
  const Eigen::VectorXd energies = MolecularDegreesOfFreedom::getEnergyGraining(
      1 * Constants::hartree_per_invCentimeter, 6000 * Constants::hartree_per_invCentimeter);
  const auto referenceState = ReferenceStates::StandardStateGas;
  HarmonicOscillator harmonicOscillatorA(std::make_shared<NormalModesContainer>(HFNormalModes), true);
  HarmonicOscillator harmonicOscillatorB(hfHessian, hfAtoms, true);
  double diff =
      (harmonicOscillatorA.densityOfState(energies) - harmonicOscillatorB.densityOfState(energies)).array().abs().maxCoeff();
  EXPECT_NEAR(0.0, diff, 1.0);
  diff = (harmonicOscillatorA.sumOfStates(energies) - harmonicOscillatorB.sumOfStates(energies)).array().abs().maxCoeff();
  EXPECT_NEAR(0.0, diff, 1e-3);
  // First excitation can be reached. So the sum of states should be just above 2.
  EXPECT_NEAR(harmonicOscillatorA.sumOfStates(energies)[energies.size() - 1], 2.0, 0.1);
  EXPECT_NEAR(harmonicOscillatorA.entropy(referenceState), harmonicOscillatorB.entropy(referenceState), 1e-6);
  EXPECT_NEAR(harmonicOscillatorA.enthalpy(referenceState), harmonicOscillatorB.enthalpy(referenceState), 1e-6);
  EXPECT_NEAR(harmonicOscillatorA.partitionFunction(referenceState),
              harmonicOscillatorB.partitionFunction(referenceState), 1e-6);
  EXPECT_NEAR(harmonicOscillatorA.heatCapacityCv(referenceState), harmonicOscillatorB.heatCapacityCv(referenceState), 1e-6);
  EXPECT_NEAR(harmonicOscillatorA.heatCapacityCp(referenceState), harmonicOscillatorB.heatCapacityCp(referenceState), 1e-6);
  EXPECT_NEAR(harmonicOscillatorA.helmholtzFreeEnergy(referenceState),
              harmonicOscillatorB.helmholtzFreeEnergy(referenceState), 1e-6);
  EXPECT_THROW(harmonicOscillatorA.densityOfState(Eigen::VectorXd::Constant(5, -1.0)), std::runtime_error);
  EXPECT_THROW(harmonicOscillatorA.sumOfStates(Eigen::VectorXd::Constant(5, -1.0)), std::runtime_error);
}

TEST_F(HarmonicOscillatorTest, ThermodynamicFunctions) {
  HarmonicOscillator harmonicOscillator(std::make_shared<NormalModesContainer>(formaldehydeNormalModes));
  double zpe = harmonicOscillator.enthalpy(ReferenceStates::VacuumZeroKelvin);
  EXPECT_NEAR(zpe * Constants::kCalPerMol_per_hartree, 15.267, 1e-3);

  EXPECT_NEAR(harmonicOscillator.heatCapacityCp(ReferenceStates::StandardStateGas) * Constants::kCalPerMol_per_hartree,
              0.6650 / 1000, 1e-5);
  EXPECT_NEAR(harmonicOscillator.entropy(ReferenceStates::StandardStateGas) * Constants::kCalPerMol_per_hartree,
              0.1365 / 1000, 1e-6);
  EXPECT_NEAR((harmonicOscillator.enthalpy(ThermodynamicReferenceState(298.00, ReferenceStates::StandardStateGas.pressure)) - zpe) *
                  Constants::kCalPerMol_per_hartree,
              0.0345759, 1e-6);
}

TEST_F(HarmonicOscillatorTest, ZeroTemperature) {
  HarmonicOscillator harmonicOscillator(std::make_shared<NormalModesContainer>(formaldehydeNormalModes));
  EXPECT_NEAR(harmonicOscillator.entropy(ReferenceStates::VacuumZeroKelvin), 0.0, 1e-12);
  EXPECT_NEAR(harmonicOscillator.heatCapacityCp(ReferenceStates::VacuumZeroKelvin), 0.0, 1e-12);
  EXPECT_NEAR(harmonicOscillator.helmholtzFreeEnergy(ReferenceStates::VacuumZeroKelvin),
              harmonicOscillator.enthalpy(ReferenceStates::VacuumZeroKelvin), 1e-12);
}

TEST_F(HarmonicOscillatorTest, EmptyEnergyRange) {
  HarmonicOscillator harmonicOscillator(std::make_shared<NormalModesContainer>(formaldehydeNormalModes));
  Eigen::VectorXd energies = Eigen::VectorXd::Zero(0);
  EXPECT_EQ(energies.size(), harmonicOscillator.densityOfState(energies).size());
  EXPECT_EQ(energies.size(), harmonicOscillator.sumOfStates(energies).size());
}

TEST_F(HarmonicOscillatorTest, SingleAtom) {
  HarmonicOscillator harmonicOscillator(std::make_shared<NormalModesContainer>(arNormalModes));
  const Eigen::VectorXd energies = MolecularDegreesOfFreedom::getEnergyGraining(
      1 * Constants::hartree_per_invCentimeter, 100 * Constants::hartree_per_invCentimeter);
  EXPECT_NEAR(1.0, harmonicOscillator.densityOfState(energies).sum(), 1e-12);
  EXPECT_NEAR(energies.size(), harmonicOscillator.sumOfStates(energies).sum(), 1e-12);
  EXPECT_NEAR(1.0, harmonicOscillator.partitionFunction(ReferenceStates::StandardStateGas), 1e-12);
  EXPECT_NEAR(0.0, harmonicOscillator.enthalpy(ReferenceStates::StandardStateGas), 1e-12);
  EXPECT_NEAR(0.0, harmonicOscillator.entropy(ReferenceStates::StandardStateGas), 1e-12);
  EXPECT_NEAR(0.0, harmonicOscillator.heatCapacityCp(ReferenceStates::StandardStateGas), 1e-12);
  EXPECT_NEAR(0.0, harmonicOscillator.heatCapacityCv(ReferenceStates::StandardStateGas), 1e-12);
}

TEST_F(HarmonicOscillatorTest, ClassicalVsQuantum) {
  HarmonicOscillator quantumHarmonicOscillator(std::make_shared<NormalModesContainer>(formaldehydeNormalModes));
  HarmonicOscillator classicalHarmonicOscillator(std::make_shared<NormalModesContainer>(formaldehydeNormalModes), true);

  auto energies = MolecularDegreesOfFreedom::getEnergyGraining(1 * Constants::hartree_per_invCentimeter,
                                                               18000 * Constants::hartree_per_invCentimeter);
  auto quantumDensityOfStates = quantumHarmonicOscillator.densityOfState(energies);
  auto classicalDensityOfStates = classicalHarmonicOscillator.densityOfState(energies);
  auto quantumSumOfStates = quantumHarmonicOscillator.sumOfStates(energies);
  auto classicalSumOfStates = classicalHarmonicOscillator.sumOfStates(energies);

  EXPECT_NEAR(quantumDensityOfStates.sum(), quantumSumOfStates[quantumDensityOfStates.size() - 1], 1e-12);
  EXPECT_NEAR(classicalDensityOfStates.sum(), classicalSumOfStates[classicalSumOfStates.size() - 1], 1.0);
  double diff = (classicalSumOfStates - quantumSumOfStates).array().abs().maxCoeff();
  EXPECT_NEAR(diff, 0.0, 600);
  EXPECT_NEAR(quantumSumOfStates[quantumDensityOfStates.size() - 1], 9924, 1e-12);
  EXPECT_NEAR(classicalSumOfStates[quantumDensityOfStates.size() - 1], 10450, 1e-1);
}

} // namespace Utils
} // namespace Scine