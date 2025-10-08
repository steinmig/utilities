/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Constants.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/Utilities/Properties.h>
#include <Utils/Properties/Thermochemistry/DegreeOfFreedom.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h>
#include <gmock/gmock.h>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {

class MolecularDegreesOfFreedomTest : public Test {
 public:
  std::shared_ptr<NormalModesContainer> formaldehydeNormalModes;
  std::shared_ptr<NormalModesContainer> HFNormalModes;
  std::shared_ptr<NormalModesContainer> arNormalModes;
  std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> formaldehydePMI;
  std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> hfPMI;
  std::shared_ptr<Geometry::Properties::PrincipalMomentsOfInertia> arPMI;
  ElementTypeCollection formaldehydeElements;
  std::shared_ptr<AtomCollection> hfAtoms;
  std::shared_ptr<AtomCollection> arAtom;
  int formaldehydeMultiplicity = 1;
  int hfMultiplicity = 1;
  int arMultiplicity = 1;
  double arbitraryEnergy = 1.0;
  Eigen::MatrixXd hfHessian = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd arHessian = Eigen::MatrixXd::Zero(3, 3);

 protected:
  void SetUp() final {
    formaldehydeNormalModes = std::make_shared<NormalModesContainer>();
    HFNormalModes = std::make_shared<NormalModesContainer>();
    arNormalModes = std::make_shared<NormalModesContainer>();
    arPMI = std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>();
    hfPMI = std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>();
    formaldehydePMI = std::make_shared<Geometry::Properties::PrincipalMomentsOfInertia>();
    formaldehydeElements = {ElementType::C, ElementType::O, ElementType::H, ElementType::H};
    NormalMode m1(1101.75, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes->add(m1);
    NormalMode m2(1157.94, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes->add(m2);
    NormalMode m3(1349.03, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes->add(m3);
    NormalMode m4(1791.24, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes->add(m4);
    NormalMode m5(2614.79, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes->add(m5);
    NormalMode m6(2664.54, DisplacementCollection::Random(4, 3));
    formaldehydeNormalModes->add(m6);
    Eigen::Vector3d eigenValues(2.8969, 21.7672, 24.6640);
    // Convert from 1e-40 g cm^2 to amu*bohr^2
    eigenValues *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
    formaldehydePMI->eigenvalues = eigenValues;
    formaldehydePMI->eigenvectors = Eigen::Matrix3d::Random();

    ElementTypeCollection hfElements = {ElementType::H, ElementType::F};
    // clang-format off
    PositionCollection hfPositions(2, 3);
    hfPositions << 0.0000000000000,    0.0000000000000,    0.0000000000000,
                   0.9655884052935,    0.0000000000000,   -0.0000001000000;
    // clang-format on
    hfPositions *= Constants::bohr_per_angstrom;
    hfAtoms = std::make_shared<AtomCollection>(hfElements, hfPositions);
    NormalMode m1HF(3968.7, DisplacementCollection::Random(1, 3));
    HFNormalModes->add(m1HF);
    Eigen::Vector3d eigenValuesHF(0.00000000, 1.4818, 1.4818);
    // Convert from 1e-40 g cm^2 to amu*bohr^2
    eigenValuesHF *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
    hfPMI->eigenvalues = eigenValuesHF;
    hfPMI->eigenvectors = Eigen::Matrix3d::Random();
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
    auto masses = Geometry::Properties::getMasses(hfElements);
    for (unsigned long i = 0; i < masses.size(); ++i) {
      hfHessian.middleRows(3 * i, 3) *= std::sqrt(masses[i]);
      hfHessian.middleCols(3 * i, 3) *= std::sqrt(masses[i]);
    }
    // unit conversion, original is milliDyn / angstrom
    hfHessian *= 1e-8;                                      // N / Angstrom = kg m / (s^2 angstrom)
    hfHessian *= Constants::meter_per_angstrom;             // J / angstrom^2
    hfHessian *= Constants::hartree_per_joule;              // hartree / angstrom^2
    hfHessian *= std::pow(Constants::angstrom_per_bohr, 2); // hartree / bohr^2

    ElementTypeCollection arElements = {ElementType::Ar};
    // clang-format off
    PositionCollection arPositions(1, 3);
    arPositions << 0.0000000000000,    0.0000000000000,    0.0000000000000;
    // clang-format on
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(3, 3);
    NormalMode m1ar(0.0, identity.row(0));
    NormalMode m2ar(0.0, identity.row(1));
    NormalMode m3ar(0.0, identity.row(2));
    arNormalModes->add(m1ar);
    arNormalModes->add(m2ar);
    arNormalModes->add(m3ar);
    arPMI->eigenvalues = Eigen::Vector3d::Zero();
    arPMI->eigenvectors = Eigen::MatrixXd::Identity(3, 3);
  }
};

TEST_F(MolecularDegreesOfFreedomTest, DifferentConstructionsAreEqual) {
  MolecularDegreesOfFreedom molecularDegreesOfFreedomA(hfHessian, *hfAtoms, 1, arbitraryEnergy);
  MolecularDegreesOfFreedom molecularDegreesOfFreedomB(HFNormalModes, hfPMI, hfAtoms->getElements(), 1, arbitraryEnergy);

  const auto state = ReferenceStates::StandardStateGas;
  EXPECT_NEAR(molecularDegreesOfFreedomA.partitionFunction(state), molecularDegreesOfFreedomB.partitionFunction(state),
              molecularDegreesOfFreedomB.partitionFunction(state) * 1e-4);
  EXPECT_NEAR(molecularDegreesOfFreedomA.entropy(state), molecularDegreesOfFreedomB.entropy(state), 1e-6);
  EXPECT_NEAR(molecularDegreesOfFreedomA.enthalpy(state), molecularDegreesOfFreedomB.enthalpy(state), 1e-6);
  EXPECT_NEAR(molecularDegreesOfFreedomA.heatCapacityCp(state), molecularDegreesOfFreedomB.heatCapacityCp(state), 1e-9);
  EXPECT_NEAR(molecularDegreesOfFreedomA.heatCapacityCv(state), molecularDegreesOfFreedomB.heatCapacityCv(state), 1e-9);
  EXPECT_NEAR(molecularDegreesOfFreedomA.helmholtzFreeEnergy(state),
              molecularDegreesOfFreedomB.helmholtzFreeEnergy(state), 1e-6);
}

TEST_F(MolecularDegreesOfFreedomTest, thermodynamicFunctions) {
  MolecularDegreesOfFreedom molecularDegreesOfFreedom(formaldehydeNormalModes, formaldehydePMI, formaldehydeElements, 1,
                                                      arbitraryEnergy, 2);
  const ThermodynamicReferenceState state(298.00, ReferenceStates::StandardStateGas.pressure);
  const double zpe = molecularDegreesOfFreedom.getVibrationalZeroPointEnergy();
  // The reference data is without the vibrational ZPE.
  EXPECT_NEAR((molecularDegreesOfFreedom.enthalpy(state) - zpe) * Constants::kCalPerMol_per_hartree,
              2403.3261 / 1000 + arbitraryEnergy * Constants::kCalPerMol_per_hartree, 3e-5);
  EXPECT_NEAR(molecularDegreesOfFreedom.heatCapacityCp(state) * Constants::kCalPerMol_per_hartree, 8.6138 / 1000, 3e-5);
  EXPECT_NEAR(molecularDegreesOfFreedom.entropy(state) * Constants::kCalPerMol_per_hartree, 52.2772 / 1000, 3e-5);
  EXPECT_NEAR(molecularDegreesOfFreedom.getRotationalDegreesOfFreedom()->entropy(state), 2.55116e-05, 1e-9);
  EXPECT_NEAR(molecularDegreesOfFreedom.getTranslationalDegreesOfFreedom()->entropy(state), 5.75756e-05, 1e-9);
  EXPECT_NEAR(molecularDegreesOfFreedom.getVibrationalDegreesOfFreedom()->entropy(state), 2.17559e-07, 1e-9);
  EXPECT_NEAR(molecularDegreesOfFreedom.getElectronicDegreesOfFreedom()->entropy(state), 0.0, 1e-12);
}

} // namespace Utils
} // namespace Scine
