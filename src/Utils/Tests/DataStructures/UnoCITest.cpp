/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/SecondQuantization/CasSpecifier.h>
#include <Utils/DataStructures/SecondQuantization/UnoCI.h>
#include <Utils/Math/IterativeDiagonalizer/SubspaceOrthogonalizer.h>
#include <gmock/gmock.h>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {
namespace UnoCI {

using namespace testing;

/**
 * @test ACasTest @file CasTest.cpp
 * Tests for the generation of CAS spaces
 */

class AUnoCITest : public Test {
 public:
  Eigen::VectorXd restrictedOccupation;
  Eigen::VectorXd arbitraryOccupation;
  Eigen::MatrixXd arbitraryCoefficients;
  void SetUp() final {
    restrictedOccupation = Eigen::VectorXd::Zero(20);
    arbitraryOccupation = Eigen::VectorXd::Zero(20);
    std::fill(restrictedOccupation.data(), restrictedOccupation.data() + 10, 2);
    arbitraryOccupation << 2, 2, 2, 2, 2, 2, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0;
    arbitraryCoefficients = Eigen::MatrixXd::Random(20, 20);
    arbitraryCoefficients += Eigen::MatrixXd(arbitraryCoefficients.transpose());
    SubspaceOrthogonalizer::qrOrthogonalize(arbitraryCoefficients, 20);
  }

 private:
};

TEST_F(AUnoCITest, SolvingWorks) {
  Eigen::MatrixXd restrictedDensity =
      arbitraryCoefficients * restrictedOccupation.asDiagonal() * arbitraryCoefficients.transpose();
  Eigen::MatrixXd coefs;
  Eigen::VectorXd occs;
  std::tie(coefs, occs) = solveOrthogonal(restrictedDensity / 2, restrictedDensity / 2);
  for (int i = 0; i < 20; ++i) {
    EXPECT_NEAR(occs[i], restrictedOccupation[i], 1e-12);
  }

  Eigen::MatrixXd overlap = Eigen::MatrixXd::Identity(20, 20);
  std::tie(coefs, occs) = solve(restrictedDensity / 2, restrictedDensity / 2, overlap);
  for (int i = 0; i < 20; ++i) {
    EXPECT_NEAR(occs[i], restrictedOccupation[i], 1e-12);
  }
}

TEST_F(AUnoCITest, EmptyCASIfRestricted) {
  Eigen::MatrixXd restrictedDensity =
      arbitraryCoefficients * restrictedOccupation.asDiagonal() * arbitraryCoefficients.transpose();
  auto cas = std::get<0>(getCas(restrictedDensity / 2, restrictedDensity / 2));
  ASSERT_TRUE(cas.empty());
}

TEST_F(AUnoCITest, RecoversCorrectCAS) {
  Eigen::MatrixXd arbitraryDensity =
      arbitraryCoefficients * arbitraryOccupation.asDiagonal() * arbitraryCoefficients.transpose();
  auto cas = std::get<0>(getCas(arbitraryDensity / 2, arbitraryDensity / 2));
  auto expectedActiveIndices = std::vector<int>{6, 7, 8, 9, 10, 11, 12, 13};
  auto expectedCoreIndices = std::vector<int>{0, 1, 2, 3, 4, 5};
  ASSERT_EQ(cas.getCasElectrons().restricted, 8);
  ASSERT_EQ(cas.getTotalNumberOfOrbitals(), 20);
  ASSERT_EQ(cas.getActiveIndices().restricted.size(), 8);
  ASSERT_EQ(cas.getInactiveIndices().restricted.size(), 6);
  for (int i = 0; i < 8; ++i) {
    EXPECT_EQ(cas.getActiveIndices().restricted[i], expectedActiveIndices[i]);
  }
  for (int i = 0; i < 6; ++i) {
    EXPECT_EQ(cas.getInactiveIndices().restricted[i], expectedCoreIndices[i]);
  }
}

TEST_F(AUnoCITest, RecoversCorrectCASWithCasSize) {
  Eigen::MatrixXd arbitraryDensity =
      arbitraryCoefficients * arbitraryOccupation.asDiagonal() * arbitraryCoefficients.transpose();
  auto cas = std::get<0>(
      getCas(arbitraryDensity / 2, arbitraryDensity / 2, {0, 0}, Utils::SecondQuantization::UnoCI::CasSizeFilter{8}));
  auto expectedActiveIndices = std::vector<int>{6, 7, 8, 9, 10, 11, 12, 13};
  auto expectedCoreIndices = std::vector<int>{0, 1, 2, 3, 4, 5};
  ASSERT_EQ(cas.getCasElectrons().restricted, 8);
  ASSERT_EQ(cas.getTotalNumberOfOrbitals(), 20);
  ASSERT_EQ(cas.getActiveIndices().restricted.size(), 8);
  ASSERT_EQ(cas.getInactiveIndices().restricted.size(), 6);
  for (int i = 0; i < 8; ++i) {
    EXPECT_EQ(cas.getActiveIndices().restricted[i], expectedActiveIndices[i]);
  }
  for (int i = 0; i < 6; ++i) {
    EXPECT_EQ(cas.getInactiveIndices().restricted[i], expectedCoreIndices[i]);
  }
}
} // namespace UnoCI
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
