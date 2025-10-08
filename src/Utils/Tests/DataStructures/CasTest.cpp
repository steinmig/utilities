/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SecondQuantization/CasSpecifier.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <gmock/gmock.h>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

using namespace testing;

/**
 * @test ACasTest @file CasTest.cpp
 * Tests for the generation of CAS spaces
 */

class ACasTest : public Test {
 public:
 private:
};

TEST_F(ACasTest, CASParsingWorks) {
  std::string casString = "8, 9, 10, 11, 12, 13";
  auto casIndices = CasGenerator::parseCasString(casString);
  EXPECT_EQ(casIndices[0], 8);
  EXPECT_EQ(casIndices[1], 9);
  EXPECT_EQ(casIndices[2], 10);
  EXPECT_EQ(casIndices[3], 11);
  EXPECT_EQ(casIndices[4], 12);
  EXPECT_EQ(casIndices[5], 13);
}

TEST_F(ACasTest, CanGenerateCasRestricted) {
  int arbitraryTotalNumberElectrons = 20;
  int arbitraryDimension = 20;
  int numberCasOrbitals = 5;
  Eigen::MatrixXd randomCoefMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
  auto mos = Utils::MolecularOrbitals::createFromRestrictedCoefficients(randomCoefMatrix);
  auto generator = CasGenerator(mos);
  auto casSpec = generator.generateRestricted(arbitraryTotalNumberElectrons, numberCasOrbitals);

  std::vector<int> indices = {8, 9, 10, 11, 12, 13};
  auto casSpec2 = generator.generateRestricted(arbitraryTotalNumberElectrons, indices);

  EXPECT_EQ(casSpec.getActiveIndices().restricted.size(), 5);
  EXPECT_EQ(casSpec.getTotalNumberOfOrbitals(), 20);
  EXPECT_EQ(casSpec.getCasElectrons().restricted, 4);
  EXPECT_EQ(casSpec2.getActiveIndices().restricted.size(), 6);
  EXPECT_EQ(casSpec2.getCasElectrons().restricted, 4);

  auto activeMOs = casSpec.getActiveOrbitalsCoefficients();
  int index = 0;
  for (int activeOrbital : casSpec.getActiveIndices().restricted) {
    for (int row = 0; row < arbitraryDimension; ++row) {
      EXPECT_EQ(activeMOs.restrictedMatrix()(row, index), mos.restrictedMatrix()(row, activeOrbital));
    }
    ++index;
  }
  auto activeMOs2 = casSpec2.getActiveOrbitalsCoefficients();
  index = 0;
  for (int activeOrbital : casSpec2.getActiveIndices().restricted) {
    for (int row = 0; row < arbitraryDimension; ++row) {
      EXPECT_EQ(activeMOs2.restrictedMatrix()(row, index), mos.restrictedMatrix()(row, activeOrbital));
    }
    ++index;
  }
  EXPECT_EQ(casSpec.getActiveIndices().restricted[0], 8);
  EXPECT_EQ(casSpec.getActiveIndices().restricted[1], 9);
  EXPECT_EQ(casSpec.getActiveIndices().restricted[2], 10);
  EXPECT_EQ(casSpec.getActiveIndices().restricted[3], 11);
  EXPECT_EQ(casSpec.getActiveIndices().restricted[4], 12);

  EXPECT_EQ(casSpec.getInactiveIndices().restricted[0], 0);
  EXPECT_EQ(casSpec.getInactiveIndices().restricted[1], 1);
  EXPECT_EQ(casSpec.getInactiveIndices().restricted[2], 2);
  EXPECT_EQ(casSpec.getInactiveIndices().restricted[3], 3);
  EXPECT_EQ(casSpec.getInactiveIndices().restricted[4], 4);
  EXPECT_EQ(casSpec.getInactiveIndices().restricted[5], 5);
  EXPECT_EQ(casSpec.getInactiveIndices().restricted[6], 6);
  EXPECT_EQ(casSpec.getInactiveIndices().restricted[7], 7);
}

TEST_F(ACasTest, CanGenerateCasUnrestricted) {
  int arbitraryTotalNumberOfElectrons = 20;
  int spinMultiplicity = 5;
  int arbitraryDimension = 20;
  int numberCasOrbitals = 5;
  std::vector<int> alphaIndices = {9, 10, 11, 12, 13};
  std::vector<int> betaIndices = {6, 7, 8, 9, 10};
  Eigen::MatrixXd randomAlphaCoefMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
  Eigen::MatrixXd randomBetaCoefMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
  auto mos = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients(randomAlphaCoefMatrix, randomBetaCoefMatrix);
  auto generator = CasGenerator(mos);
  auto casSpec = generator.generateUnrestricted(arbitraryTotalNumberOfElectrons, spinMultiplicity, numberCasOrbitals);

  auto casSpec2 = generator.generateUnrestricted(arbitraryTotalNumberOfElectrons, spinMultiplicity, alphaIndices, betaIndices);

  EXPECT_EQ(casSpec.getActiveIndices().alpha.size(), 5);
  EXPECT_EQ(casSpec.getActiveIndices().beta.size(), 5);
  EXPECT_EQ(casSpec2.getActiveIndices().alpha.size(), 5);
  EXPECT_EQ(casSpec2.getActiveIndices().beta.size(), 5);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha.size(), 10);
  EXPECT_EQ(casSpec.getInactiveIndices().beta.size(), 6);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha.size(), 8);
  EXPECT_EQ(casSpec2.getInactiveIndices().beta.size(), 5);
  EXPECT_EQ(casSpec.getTotalNumberOfOrbitals(), 20);
  EXPECT_EQ(casSpec2.getTotalNumberOfOrbitals(), 20);
  EXPECT_EQ(casSpec.getCasElectrons().alpha, 2);
  EXPECT_EQ(casSpec.getCasElectrons().beta, 2);

  auto activeMOs = casSpec.getActiveOrbitalsCoefficients();
  int index = 0;
  for (int activeOrbital : casSpec.getActiveIndices().alpha) {
    for (int row = 0; row < arbitraryDimension; ++row) {
      EXPECT_EQ(activeMOs.alphaMatrix()(row, index), mos.alphaMatrix()(row, activeOrbital));
    }
    ++index;
  }
  index = 0;
  for (int activeOrbital : casSpec.getActiveIndices().beta) {
    for (int row = 0; row < arbitraryDimension; ++row) {
      EXPECT_EQ(activeMOs.betaMatrix()(row, index), mos.betaMatrix()(row, activeOrbital));
    }
    ++index;
  }

  auto activeMOs2 = casSpec2.getActiveOrbitalsCoefficients();
  index = 0;
  for (int activeOrbital : casSpec2.getActiveIndices().alpha) {
    for (int row = 0; row < arbitraryDimension; ++row) {
      EXPECT_EQ(activeMOs2.alphaMatrix()(row, index), mos.alphaMatrix()(row, activeOrbital));
    }
    ++index;
  }
  index = 0;
  for (int activeOrbital : casSpec2.getActiveIndices().beta) {
    for (int row = 0; row < arbitraryDimension; ++row) {
      EXPECT_EQ(activeMOs2.betaMatrix()(row, index), mos.betaMatrix()(row, activeOrbital));
    }
    ++index;
  }

  EXPECT_EQ(casSpec.getActiveIndices().alpha[0], 10);
  EXPECT_EQ(casSpec.getActiveIndices().alpha[1], 11);
  EXPECT_EQ(casSpec.getActiveIndices().alpha[2], 12);
  EXPECT_EQ(casSpec.getActiveIndices().alpha[3], 13);
  EXPECT_EQ(casSpec.getActiveIndices().alpha[4], 14);
  EXPECT_EQ(casSpec.getActiveIndices().beta[0], 6);
  EXPECT_EQ(casSpec.getActiveIndices().beta[1], 7);
  EXPECT_EQ(casSpec.getActiveIndices().beta[2], 8);
  EXPECT_EQ(casSpec.getActiveIndices().beta[3], 9);
  EXPECT_EQ(casSpec.getActiveIndices().beta[4], 10);
  EXPECT_EQ(casSpec2.getActiveIndices().alpha[0], 9);
  EXPECT_EQ(casSpec2.getActiveIndices().alpha[1], 10);
  EXPECT_EQ(casSpec2.getActiveIndices().alpha[2], 11);
  EXPECT_EQ(casSpec2.getActiveIndices().alpha[3], 12);
  EXPECT_EQ(casSpec2.getActiveIndices().alpha[4], 13);
  EXPECT_EQ(casSpec2.getActiveIndices().beta[0], 6);
  EXPECT_EQ(casSpec2.getActiveIndices().beta[1], 7);
  EXPECT_EQ(casSpec2.getActiveIndices().beta[2], 8);
  EXPECT_EQ(casSpec2.getActiveIndices().beta[3], 9);
  EXPECT_EQ(casSpec2.getActiveIndices().beta[4], 10);

  EXPECT_EQ(casSpec.getInactiveIndices().alpha[0], 0);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[1], 1);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[2], 2);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[3], 3);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[4], 4);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[5], 5);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[6], 6);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[7], 7);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[8], 8);
  EXPECT_EQ(casSpec.getInactiveIndices().alpha[9], 9);

  EXPECT_EQ(casSpec.getInactiveIndices().beta[0], 0);
  EXPECT_EQ(casSpec.getInactiveIndices().beta[1], 1);
  EXPECT_EQ(casSpec.getInactiveIndices().beta[2], 2);
  EXPECT_EQ(casSpec.getInactiveIndices().beta[3], 3);
  EXPECT_EQ(casSpec.getInactiveIndices().beta[4], 4);
  EXPECT_EQ(casSpec.getInactiveIndices().beta[5], 5);

  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[0], 0);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[1], 1);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[2], 2);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[3], 3);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[4], 4);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[5], 5);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[6], 6);
  EXPECT_EQ(casSpec2.getInactiveIndices().alpha[7], 7);

  EXPECT_EQ(casSpec2.getInactiveIndices().beta[0], 0);
  EXPECT_EQ(casSpec2.getInactiveIndices().beta[1], 1);
  EXPECT_EQ(casSpec2.getInactiveIndices().beta[2], 2);
  EXPECT_EQ(casSpec2.getInactiveIndices().beta[3], 3);
  EXPECT_EQ(casSpec2.getInactiveIndices().beta[4], 4);
}

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
