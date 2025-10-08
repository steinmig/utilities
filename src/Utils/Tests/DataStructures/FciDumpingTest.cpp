/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/SecondQuantization/FciDumping.h>
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

class AFciDumpingTest : public Test {
 public:
  std::stringstream fakeOutStream;
  int nMOs = 5;
  FciDumpData restrictedFciDumpData;
  FciDumpData unrestrictedFciDumpData;
  MOTypeTrait<Reference::Restricted>::Type restrictedEri;
  MOTypeTrait<Reference::Unrestricted>::Type unrestrictedEri;
  MOTypeTrait<Reference::Restricted>::OneElectronMatrixType restricted1el;
  MOTypeTrait<Reference::Unrestricted>::OneElectronMatrixType unrestricted1el;

  void SetUp() final {
    restrictedEri = Eigen::MatrixXd::Random(nMOs * nMOs, nMOs * nMOs);
    restrictedEri = (restrictedEri + restrictedEri.transpose()).eval();
    unrestrictedEri[{SpinComponent::Alpha, SpinComponent::Alpha}] = Eigen::MatrixXd::Random(nMOs * nMOs, nMOs * nMOs);
    unrestrictedEri[{SpinComponent::Alpha, SpinComponent::Beta}] = Eigen::MatrixXd::Random(nMOs * nMOs, nMOs * nMOs);
    unrestrictedEri[{SpinComponent::Beta, SpinComponent::Alpha}] = Eigen::MatrixXd::Random(nMOs * nMOs, nMOs * nMOs);
    unrestrictedEri[{SpinComponent::Beta, SpinComponent::Beta}] = Eigen::MatrixXd::Random(nMOs * nMOs, nMOs * nMOs);
    for (int row = 0; row < nMOs * nMOs; ++row) {
      for (int col = 0; col < nMOs * nMOs; ++col) {
        int i = row / nMOs;
        int j = row % nMOs;
        int k = col / nMOs;
        int l = col % nMOs;
        auto symIdx = getMappedIndex<ERISymmetry::eightfold>({{i, j, k, l}});
        unrestrictedEri.at({SpinComponent::Alpha, SpinComponent::Alpha})(row, col) = unrestrictedEri.at(
            {SpinComponent::Alpha, SpinComponent::Alpha})(symIdx[0] * nMOs + symIdx[1], symIdx[2] * nMOs + symIdx[3]);
        unrestrictedEri.at({SpinComponent::Alpha, SpinComponent::Beta})(row, col) = unrestrictedEri.at(
            {SpinComponent::Alpha, SpinComponent::Beta})(symIdx[0] * nMOs + symIdx[1], symIdx[2] * nMOs + symIdx[3]);
        unrestrictedEri.at({SpinComponent::Beta, SpinComponent::Alpha})(row, col) = unrestrictedEri.at(
            {SpinComponent::Beta, SpinComponent::Alpha})(symIdx[0] * nMOs + symIdx[1], symIdx[2] * nMOs + symIdx[3]);
        unrestrictedEri.at({SpinComponent::Beta, SpinComponent::Beta})(row, col) = unrestrictedEri.at(
            {SpinComponent::Beta, SpinComponent::Beta})(symIdx[0] * nMOs + symIdx[1], symIdx[2] * nMOs + symIdx[3]);
        restrictedEri(row, col) = restrictedEri(symIdx[0] * nMOs + symIdx[1], symIdx[2] * nMOs + symIdx[3]);
      }
    }
    restricted1el = Eigen::MatrixXd::Random(nMOs, nMOs);
    restricted1el = (restricted1el.transpose() + restricted1el).eval();
    unrestricted1el[{SpinComponent::Alpha}] = Eigen::MatrixXd::Random(nMOs, nMOs);
    unrestricted1el[{SpinComponent::Beta}] = Eigen::MatrixXd::Random(nMOs, nMOs);
    unrestricted1el.at({SpinComponent::Alpha}) =
        (unrestricted1el.at({SpinComponent::Alpha}).transpose() + unrestricted1el.at({SpinComponent::Alpha})).eval();
    unrestricted1el.at({SpinComponent::Beta}) =
        (unrestricted1el.at({SpinComponent::Beta}).transpose() + unrestricted1el.at({SpinComponent::Beta})).eval();

    restrictedFciDumpData.nOrbitals = nMOs;
    restrictedFciDumpData.nElectrons = 4;
    restrictedFciDumpData.coreEnergy = -42.42;
    restrictedFciDumpData.spinPolarization = 0;
    restrictedFciDumpData.unrestrictedReference = false;

    unrestrictedFciDumpData.nOrbitals = nMOs;
    unrestrictedFciDumpData.nElectrons = 4;
    unrestrictedFciDumpData.coreEnergy = -42.42;
    unrestrictedFciDumpData.spinPolarization = 0;
    unrestrictedFciDumpData.unrestrictedReference = true;
  }
  void TearDown() final {
    std::remove("fcidump.fci");
  }

 private:
};

TEST_F(AFciDumpingTest, CanWriteOneElectronIntegralsRestricted) {
  fakeOutStream << std::scientific << std::left << std::setprecision(15);
  FciDumping::writeOneElectronIntegrals<Utils::Reference::Restricted>(fakeOutStream, restricted1el);
  std::string buffer;
  double value = 0;
  int i = 0;
  int j = 0;
  while (std::getline(fakeOutStream, buffer)) {
    std::stringstream linebuffer(buffer);
    linebuffer >> value >> i >> j;
    EXPECT_NEAR(restricted1el(i - 1, j - 1), value, 1e-12);
  }
}

TEST_F(AFciDumpingTest, CanWriteOneElectronIntegralsUnrestricted) {
  fakeOutStream << std::scientific << std::left << std::setprecision(15);
  FciDumping::writeOneElectronIntegrals<Utils::Reference::Unrestricted>(fakeOutStream, unrestricted1el);
  std::string buffer;
  double value = 0;
  int i = 0;
  int j = 0;
  while (std::getline(fakeOutStream, buffer)) {
    std::stringstream linebuffer(buffer);
    linebuffer >> value >> i >> j;
    if ((i - 1) % 2 == 0) {
      EXPECT_NEAR(unrestricted1el.at(SpinComponent::Alpha)((i - 1) / 2, (j - 1) / 2), value, 1e-12);
    }
    if ((i - 1) % 2 == 1) {
      EXPECT_NEAR(unrestricted1el.at(SpinComponent::Beta)(i / 2 - 1, j / 2 - 1), value, 1e-12);
    }
  }
}

TEST_F(AFciDumpingTest, CanWriteTwoElectronIntegralsRestricted) {
  fakeOutStream << std::scientific << std::left << std::setprecision(15);
  FciDumping::writeTwoElectronIntegrals<Utils::Reference::Restricted>(fakeOutStream, restrictedEri);
  std::string buffer;
  double value = 0;
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  while (std::getline(fakeOutStream, buffer)) {
    std::stringstream linebuffer(buffer);
    linebuffer >> value >> i >> j >> k >> l;
    EXPECT_NEAR(restrictedEri((i - 1) * nMOs + (j - 1), (k - 1) * nMOs + l - 1), value, 1e-12);
  }
}

TEST_F(AFciDumpingTest, CanWriteTwoElectronIntegralsUnrestricted) {
  fakeOutStream << std::scientific << std::left << std::setprecision(15);
  FciDumping::writeTwoElectronIntegrals<Utils::Reference::Unrestricted>(fakeOutStream, unrestrictedEri);
  std::string buffer;
  double value = 0;
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int p = 0;
  int q = 0;
  int r = 0;
  int s = 0;
  SpinComponent spin1;
  SpinComponent spin2;
  while (std::getline(fakeOutStream, buffer)) {
    std::stringstream linebuffer(buffer);
    linebuffer >> value >> i >> j >> k >> l;
    if ((i - 1) % 2 == 0) {
      spin1 = SpinComponent::Alpha;
      p = (i - 1) / 2;
      q = (j - 1) / 2;
    }
    else if ((i - 1) % 2 == 1) {
      spin1 = SpinComponent::Beta;
      p = (i - 2) / 2;
      q = (j - 2) / 2;
    }
    if ((k - 1) % 2 == 0) {
      spin2 = SpinComponent::Alpha;
      r = (k - 1) / 2;
      s = (l - 1) / 2;
    }
    else if ((k - 1) % 2 == 1) {
      spin2 = SpinComponent::Beta;
      r = (k - 2) / 2;
      s = (l - 2) / 2;
    }
    EXPECT_NEAR(unrestrictedEri.at({spin1, spin2})(p * nMOs + q, r * nMOs + s), value, 1e-12);
  }
}

TEST_F(AFciDumpingTest, FciDumpCanWriteRestricted) {
  FciDumping::OneElectronMatrixType<Utils::Reference::Restricted> oneElectronMatrix;
  FciDumping::TwoElectronMatrixType<Utils::Reference::Restricted> twoElectronMatrix =
      Eigen::MatrixXd::Zero(nMOs * nMOs, nMOs * nMOs);
  FciDumpData data;

  FciDumping::write<Utils::Reference::Restricted>("fcidump.fci", restricted1el, restrictedEri, restrictedFciDumpData);
  std::tie(oneElectronMatrix, twoElectronMatrix, data) = FciDumping::read<Utils::Reference::Restricted>("fcidump.fci");

  ASSERT_EQ(restricted1el.rows(), oneElectronMatrix.rows());
  ASSERT_EQ(restricted1el.cols(), oneElectronMatrix.cols());
  ASSERT_EQ(restrictedEri.rows(), twoElectronMatrix.rows());
  ASSERT_EQ(restrictedEri.cols(), twoElectronMatrix.cols());
  for (int row = 0; row < restricted1el.rows(); ++row) {
    for (int col = 0; col < restricted1el.cols(); ++col) {
      EXPECT_NEAR(oneElectronMatrix(row, col), restricted1el(row, col), 1e-12);
    }
  }
  for (int row = 0; row < restrictedEri.rows(); ++row) {
    for (int col = 0; col < restrictedEri.cols(); ++col) {
      SCOPED_TRACE("p = " + std::to_string(row / nMOs) + " q=" + std::to_string(row % nMOs) +
                   " r=" + std::to_string(col / nMOs) + " s=" + std::to_string(col % nMOs));
      EXPECT_NEAR(twoElectronMatrix(row, col), restrictedEri(row, col), 1e-12);
    }
  }
}

TEST_F(AFciDumpingTest, FciDumpCanWriteUnrestricted) {
  FciDumping::OneElectronMatrixType<Utils::Reference::Unrestricted> oneElectronMatrix;
  FciDumping::TwoElectronMatrixType<Utils::Reference::Unrestricted> twoElectronMatrix;
  FciDumpData data;

  FciDumping::write<Utils::Reference::Unrestricted>("fcidump.fci", unrestricted1el, unrestrictedEri, unrestrictedFciDumpData);
  std::tie(oneElectronMatrix, twoElectronMatrix, data) =
      FciDumping::read<Utils::Reference::Unrestricted>("fcidump.fci");
  for (auto spin1 : {SpinComponent::Alpha, SpinComponent::Beta}) {
    ASSERT_EQ(unrestricted1el.at(spin1).rows(), oneElectronMatrix.at(spin1).rows());
    ASSERT_EQ(unrestricted1el.at(spin1).cols(), oneElectronMatrix.at(spin1).cols());
    for (int row = 0; row < restricted1el.rows(); ++row) {
      for (int col = 0; col < restricted1el.cols(); ++col) {
        EXPECT_NEAR(oneElectronMatrix.at(spin1)(row, col), unrestricted1el.at(spin1)(row, col), 1e-12);
      }
    }
    for (auto spin2 : {SpinComponent::Alpha, SpinComponent::Beta}) {
      ASSERT_EQ(unrestrictedEri.at({spin1, spin2}).rows(), twoElectronMatrix.at({spin1, spin2}).rows());
      ASSERT_EQ(unrestrictedEri.at({spin1, spin2}).cols(), twoElectronMatrix.at({spin1, spin2}).cols());
      for (int row = 0; row < restrictedEri.rows(); ++row) {
        for (int col = 0; col < restrictedEri.cols(); ++col) {
          EXPECT_NEAR(twoElectronMatrix.at({spin1, spin2})(row, col), unrestrictedEri.at({spin1, spin2})(row, col), 1e-12);
        }
      }
    }
  }
}
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
