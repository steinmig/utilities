/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/DataStructures/SecondQuantization/CICalculatorSettings.h"
#include "Utils/UniversalSettings/SettingsNames.h"
#include <Utils/CalculatorBasics/TestCalculator.h>
#include <Utils/DataStructures/SecondQuantization/CICalculator.h>
#include <Utils/DataStructures/SecondQuantization/CISigmaVectorEvaluator.h>
#include <Utils/DataStructures/SecondQuantization/CISolver.h>
#include <Utils/DataStructures/SecondQuantization/CasIndicesHandler.h>
#include <Utils/DataStructures/SecondQuantization/CasSpecifier.h>
#include <Utils/DataStructures/SecondQuantization/HermitianHamiltonian.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectPreconditionerEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectSigmaVectorEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/SubspaceOrthogonalizer.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

class SQTestCalculator : public Utils::CloneInterface<SQTestCalculator, Core::Calculator> {
 private:
  // PropertyList _requiredProperties{};
  // AtomCollection _structure;
  // Results _results;
  TestCalculator _calc;
  // std::shared_ptr<Settings> _settings;
  /**
   * @brief The default precision of the calculator.
   */
  double _precision = 13.0;
  /**
   * @brief Truncate a given value to the set precision.
   * @param value The value to be truncated
   * @return The truncated value.
   */
  double truncateOff(double value);

 public:
  SQTestCalculator() = default;
  ~SQTestCalculator() override = default;
  PropertyList possibleProperties() const final {
    return Property::Energy | Property::Gradients | Property::Hessian | Property::SQSpecifier;
  }

  PropertyList getRequiredProperties() const final {
    return _calc.getRequiredProperties();
  }

  std::string name() const final {
    return std::string("SQTestCalculator");
  }

  const Settings& settings() const final {
    return _calc.settings();
  }
  Settings& settings() final {
    return _calc.settings();
  }

  Utils::Results& results() final {
    return _calc.results();
  }
  const Utils::Results& results() const final {
    return _calc.results();
  }
  bool allowsPythonGILRelease() const final {
    return false;
  }
  bool supportsMethodFamily(const std::string& methodFamily) const override {
    return methodFamily == "TEST";
  }

  void setStructure(const AtomCollection& structure) final {
    _calc.setStructure(structure);
  }

  std::unique_ptr<Utils::AtomCollection> getStructure() const final {
    return _calc.getStructure();
  }

  void modifyPositions(PositionCollection newPositions) final {
    _calc.getStructure()->setPositions(newPositions);
  }

  void loadState(std::shared_ptr<Core::State> state /* state */) final {
  }

  std::shared_ptr<Core::State> getState() const final {
    return nullptr;
  }

  const PositionCollection& getPositions() const final {
    return _calc.getStructure()->getPositions();
  }

  void setRequiredProperties(const PropertyList& requiredProperties) final {
    _calc.setRequiredProperties(requiredProperties);
  }

  const Results& calculate(std::string dummy = "") final {
    const Results& results = _calc.calculate(dummy);
    if (getRequiredProperties().containsSubSet(Property::SQSpecifier)) {
      setSQSpecifier();
    }
    return results;
  }

  void setSQSpecifier() {
    // mo integrals from butadiene OM3
    const int casSize = 4;
    Eigen::MatrixXd oneElectronInts(casSize, casSize);
    oneElectronInts << -1.45070766, -0.02864755, -0.00469431, -0.00838718, -0.02864755, -1.45471637, 0.00968231,
        0.00159671, -0.00469431, 0.00968231, -0.90130713, 0.02252204, -0.00838718, 0.00159671, 0.02252204, -0.7521663;
    Eigen::MatrixXd twoElectronInts(casSize * casSize, casSize * casSize);
    twoElectronInts << 3.31751991e-01, 1.19451769e-02, -4.15526616e-03, 4.00263967e-03, 1.19451769e-02, 2.25599696e-01,
        -2.42485270e-03, -2.47617573e-03, -4.15526616e-03, -2.42485270e-03, 2.42113078e-01, -1.75001472e-02,
        4.00263967e-03, -2.47617573e-03, -1.75001472e-02, 2.34950947e-01, 1.19451769e-02, 1.07950393e-02, 3.72189910e-04,
        -9.94742513e-04, 1.07950393e-02, -7.16224593e-04, -3.96600320e-05, 6.39828696e-04, 3.72189910e-04,
        -3.96600320e-05, 2.19015755e-03, -9.14013377e-03, -9.94742513e-04, 6.39828696e-04, -9.14013377e-03,
        6.79268780e-04, -4.15526616e-03, 3.72189910e-04, 3.60094623e-03, -1.15858122e-03, 3.72189910e-04,
        4.40496050e-03, 7.99926052e-04, -6.64446133e-04, 3.60094623e-03, 7.99926052e-04, 2.50474908e-03, 5.64475869e-04,
        -1.15858122e-03, -6.64446133e-04, 5.64475869e-04, 3.88690201e-03, 4.00263967e-03, -9.94742513e-04,
        -1.15858122e-03, 2.42660052e-03, -9.94742513e-04, 2.51218698e-03, -1.38608690e-03, -4.77287791e-04,
        -1.15858122e-03, -1.38608690e-03, 5.42978223e-03, 9.14561879e-04, 2.42660052e-03, -4.77287791e-04,
        9.14561879e-04, 2.85169515e-03, 1.19451769e-02, 1.07950393e-02, 3.72189910e-04, -9.94742513e-04, 1.07950393e-02,
        -7.16224593e-04, -3.96600320e-05, 6.39828696e-04, 3.72189910e-04, -3.96600320e-05, 2.19015755e-03,
        -9.14013377e-03, -9.94742513e-04, 6.39828696e-04, -9.14013377e-03, 6.79268780e-04, 2.25599696e-01,
        -7.16224593e-04, 4.40496050e-03, 2.51218698e-03, -7.16224593e-04, 2.86672034e-01, -4.46041334e-03,
        2.36089780e-03, 4.40496050e-03, -4.46041334e-03, 2.81537043e-01, 1.89615263e-04, 2.51218698e-03, 2.36089780e-03,
        1.89615263e-04, 2.61950153e-01, -2.42485270e-03, -3.96600320e-05, 7.99926052e-04, -1.38608690e-03,
        -3.96600320e-05, -4.46041334e-03, 1.77913099e-02, -4.37728279e-04, 7.99926052e-04, 1.77913099e-02,
        -1.92229679e-02, 9.46282347e-04, -1.38608690e-03, -4.37728279e-04, 9.46282347e-04, -3.91936124e-03,
        -2.47617573e-03, 6.39828696e-04, -6.64446133e-04, -4.77287791e-04, 6.39828696e-04, 2.36089780e-03,
        -4.37728279e-04, 4.04945948e-03, -6.64446133e-04, -4.37728279e-04, -1.28286331e-03, -1.46808717e-03,
        -4.77287791e-04, 4.04945948e-03, -1.46808717e-03, -5.13424028e-03, -4.15526616e-03, 3.72189910e-04,
        3.60094623e-03, -1.15858122e-03, 3.72189910e-04, 4.40496050e-03, 7.99926052e-04, -6.64446133e-04, 3.60094623e-03,
        7.99926052e-04, 2.50474908e-03, 5.64475869e-04, -1.15858122e-03, -6.64446133e-04, 5.64475869e-04, 3.88690201e-03,
        -2.42485270e-03, -3.96600320e-05, 7.99926052e-04, -1.38608690e-03, -3.96600320e-05, -4.46041334e-03,
        1.77913099e-02, -4.37728279e-04, 7.99926052e-04, 1.77913099e-02, -1.92229679e-02, 9.46282347e-04,
        -1.38608690e-03, -4.37728279e-04, 9.46282347e-04, -3.91936124e-03, 2.42113078e-01, 2.19015755e-03,
        2.50474908e-03, 5.42978223e-03, 2.19015755e-03, 2.81537043e-01, -1.92229679e-02, -1.28286331e-03, 2.50474908e-03,
        -1.92229679e-02, 3.27147018e-01, -2.32944763e-03, 5.42978223e-03, -1.28286331e-03, -2.32944763e-03,
        2.67503360e-01, -1.75001472e-02, -9.14013377e-03, 5.64475869e-04, 9.14561879e-04, -9.14013377e-03,
        1.89615263e-04, 9.46282347e-04, -1.46808717e-03, 5.64475869e-04, 9.46282347e-04, -2.32944763e-03, 1.61277703e-02,
        9.14561879e-04, -1.46808717e-03, 1.61277703e-02, -1.17555198e-03, 4.00263967e-03, -9.94742513e-04,
        -1.15858122e-03, 2.42660052e-03, -9.94742513e-04, 2.51218698e-03, -1.38608690e-03, -4.77287791e-04,
        -1.15858122e-03, -1.38608690e-03, 5.42978223e-03, 9.14561879e-04, 2.42660052e-03, -4.77287791e-04,
        9.14561879e-04, 2.85169515e-03, -2.47617573e-03, 6.39828696e-04, -6.64446133e-04, -4.77287791e-04,
        6.39828696e-04, 2.36089780e-03, -4.37728279e-04, 4.04945948e-03, -6.64446133e-04, -4.37728279e-04,
        -1.28286331e-03, -1.46808717e-03, -4.77287791e-04, 4.04945948e-03, -1.46808717e-03, -5.13424028e-03,
        -1.75001472e-02, -9.14013377e-03, 5.64475869e-04, 9.14561879e-04, -9.14013377e-03, 1.89615263e-04,
        9.46282347e-04, -1.46808717e-03, 5.64475869e-04, 9.46282347e-04, -2.32944763e-03, 1.61277703e-02, 9.14561879e-04,
        -1.46808717e-03, 1.61277703e-02, -1.17555198e-03, 2.34950947e-01, 6.79268780e-04, 3.88690201e-03, 2.85169515e-03,
        6.79268780e-04, 2.61950153e-01, -3.91936124e-03, -5.13424028e-03, 3.88690201e-03, -3.91936124e-03,
        2.67503360e-01, -1.17555198e-03, 2.85169515e-03, -5.13424028e-03, -1.17555198e-03, 2.65841432e-01;

    SecondQuantization::MoIntegrals ints{oneElectronInts, twoElectronInts, {}};
    double coreEnergy = -17.29867357102363;

    Eigen::MatrixXd randomMOs = Eigen::MatrixXd::Random(32, 32);
    SubspaceOrthogonalizer::qrOrthogonalize(randomMOs, 32);
    auto mos = MolecularOrbitals::createFromRestrictedCoefficients(randomMOs);

    SecondQuantization::CasGenerator generator(mos);
    auto cas = generator.generateRestricted(32, 4);

    Eigen::VectorXd arbitraryOccs = Eigen::VectorXd::Zero(32);
    arbitraryOccs.head(16).array() = 2;

    _calc.results().set<Property::SQSpecifier>({ints, {}, cas, {}, arbitraryOccs, coreEnergy});
  }
};

using namespace testing;
class ACICalculation : public Test {
 public:
  std::shared_ptr<Core::Calculator> calc;
  std::shared_ptr<SecondQuantization::CICalculator> ci;
  Core::Log logger;

 private:
  void SetUp() final {
    ci = std::make_shared<SecondQuantization::CICalculator>();
    ci->settings().modifyBool(SettingsNames::directness, false);
    calc = std::make_shared<SQTestCalculator>();
    logger = Core::Log::silent();
    ci->setLog(logger);
  }
};

TEST_F(ACICalculation, CanInitializeCICalculator) {
  ASSERT_NO_THROW(ci->setReferenceCalculator(calc));
}

TEST_F(ACICalculation, ThrowsIfInvalidCalculatorSet) {
  ASSERT_THROW(ci->setReferenceCalculator(std::make_shared<TestCalculator>()), std::runtime_error);
}

TEST_F(ACICalculation, ThrowsInCalculationIfNoCalculatorSet) {
  ASSERT_THROW(ci->calculate(), std::runtime_error);
}

TEST_F(ACICalculation, CasIndicesAreCorrectlyGenerated) {
  Eigen::MatrixXd arbitraryMOs = Eigen::MatrixXd::Random(32, 32);
  SubspaceOrthogonalizer::qrOrthogonalize(arbitraryMOs, 32);
  auto mos = MolecularOrbitals::createFromRestrictedCoefficients(arbitraryMOs);
  CasGenerator generator(mos);

  // arbitrary cas space with non-contiguous orbitals
  auto arbitraryCas = generator.generateRestricted(32, {11, 12, 14, 15, 16, 17, 19, 20});
  CasIndicesHandler handler(arbitraryCas);

  auto idxAroundFermi = handler.getIndices(4);
  auto idxSpecified = handler.getIndices({14, 15, 16, 17});

  std::vector<int> idxExpected{2, 3, 4, 5};
  ASSERT_EQ(idxAroundFermi.size(), idxExpected.size());
  ASSERT_EQ(idxSpecified.size(), idxExpected.size());

  for (int i = 0; i < int(idxExpected.size()); ++i) {
    EXPECT_EQ(idxAroundFermi[i], idxExpected[i]);
    EXPECT_EQ(idxSpecified[i], idxExpected[i]);
  }

  // Test that magic numbers work:
  auto allIdx = handler.getIndices(DoubleExcitationSpecifiers::FullSpace);
  auto noIdx = handler.getIndices(DoubleExcitationSpecifiers::OnlySingleExcitations);
  auto noRefIdx = handler.getIndices(ReferenceSpecifiers::SingleReference);
  EXPECT_EQ(noIdx.size(), 0);
  EXPECT_EQ(noRefIdx.size(), 0);
  ASSERT_EQ(allIdx.size(), 8);
  for (int i = 0; i < int(allIdx.size()); ++i) {
    EXPECT_EQ(allIdx[i], i);
  }

  ASSERT_THROW(handler.getIndices({13, 14, 15, 16}), std::runtime_error);
}

TEST_F(ACICalculation, CanGenerateConnectedDeterminants) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  const auto& specifier = calc->results().get<Property::SQSpecifier>();
  HermitianHamiltonian ham(specifier);
  ElectronicDeterminant refDet("2200");
  auto dets = ham.generateAllConnected(refDet, {{SpinComponent::Alpha, {1}}, {SpinComponent::Beta, {1}}},
                                       {{SpinComponent::Alpha, {2}}, {SpinComponent::Beta, {2}}}, true);
  std::vector<ElectronicDeterminant> expectedDets{{"2200"}, {"2020"}, {"2ab0"}, {"2a0b"}, {"2ba0"},
                                                  {"2b0a"}, {"a2b0"}, {"a20b"}, {"b2a0"}, {"b20a"}};
  ASSERT_EQ(dets.size(), expectedDets.size());
  for (const auto& det : dets) {
    EXPECT_TRUE(std::find(expectedDets.begin(), expectedDets.end(), det) != expectedDets.end());
  }

  dets = ham.generateAllConnected(refDet, {{SpinComponent::Alpha, {0, 1}}, {SpinComponent::Beta, {0, 1}}},
                                  {{SpinComponent::Alpha, {2, 3}}, {SpinComponent::Beta, {2, 3}}}, true);
  expectedDets = ham.generateAllConnected(refDet, true);

  ASSERT_EQ(dets.size(), expectedDets.size());
  for (const auto& det : dets) {
    EXPECT_TRUE(std::find(expectedDets.begin(), expectedDets.end(), det) != expectedDets.end());
  }
}

TEST_F(ACICalculation, CanCalculateCIProblem) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  const auto& specifier = calc->results().get<Property::SQSpecifier>();
  CISolver solver(ci->settings(), specifier);
  solver.setIndices(MOIndicesForReferences{{0, 1, 2, 3}});
  auto result = solver.solve(logger);
  double firstExcitedState = (result.eigenValues(1) - result.eigenValues(0)) * Utils::Constants::ev_per_hartree;
  double secondExcitedState = (result.eigenValues(2) - result.eigenValues(0)) * Utils::Constants::ev_per_hartree;

  // triplet state, checked with DMRG
  EXPECT_NEAR(firstExcitedState, 14.9276, 1e-4);
  // singlet state, checked with DMRG
  EXPECT_NEAR(secondExcitedState, 15.256, 1e-3);
}

TEST_F(ACICalculation, CanCalculateCIProblemWithCalculator) {
  ci->setReferenceCalculator(calc);
  ci->settings().modifyString(SettingsNames::doublesIndicesOption, "");
  ci->settings().modifyBool(SettingsNames::casWithUnoCI, true);
  ci->settings().modifyInt(SettingsNames::casOrbitalsAroundFermiOption, 4);
  ci->settings().modifyInt(SettingsNames::referencesAroundFermiOption, 4);
  auto result = ci->calculate();
  auto res = result.get<Property::ExcitedStates>().unrestricted->eigenStates;
  double firstExcitedState = (res.eigenValues(1) - res.eigenValues(0)) * Utils::Constants::ev_per_hartree;
  double secondExcitedState = (res.eigenValues(2) - res.eigenValues(0)) * Utils::Constants::ev_per_hartree;

  EXPECT_EQ(result.get<Property::ExcitedStates>().transitionLabels.size(), res.eigenValues.size());

  // triplet state, checked with DMRG
  EXPECT_NEAR(firstExcitedState, 14.9276, 1e-4);
  // singlet state, checked with DMRG
  EXPECT_NEAR(secondExcitedState, 15.256, 1e-3);
}

TEST_F(ACICalculation, CanCalculateCIProblemWithCalculatorAndDavidson) {
  ci->setReferenceCalculator(calc);
  ci->settings().modifyString(SettingsNames::doublesIndicesOption, "");
  ci->settings().modifyBool(SettingsNames::casWithUnoCI, true);
  ci->settings().modifyInt(SettingsNames::casOrbitalsAroundFermiOption, 4);
  ci->settings().modifyInt(SettingsNames::referencesAroundFermiOption, 4);
  ci->settings().modifyInt(SettingsNames::numberOfEigenstates, 3);
  auto result = ci->calculate();
  auto res = result.get<Property::ExcitedStates>().unrestricted->eigenStates;
  double firstExcitedState = (res.eigenValues(1) - res.eigenValues(0)) * Utils::Constants::ev_per_hartree;
  double secondExcitedState = (res.eigenValues(2) - res.eigenValues(0)) * Utils::Constants::ev_per_hartree;
  Eigen::VectorXd multiplicities = result.get<Property::ExcitedStates>().unrestricted->multiplicities;

  EXPECT_EQ(result.get<Property::ExcitedStates>().transitionLabels.size(), 36);
  EXPECT_EQ(res.eigenValues.size(), 3);

  EXPECT_NEAR(multiplicities(0), 1, 1e-6);
  EXPECT_NEAR(multiplicities(1), 3, 1e-6);
  EXPECT_NEAR(multiplicities(2), 1, 1e-6);

  // triplet state, checked with DMRG
  EXPECT_NEAR(firstExcitedState, 14.9276, 1e-4);
  // singlet state, checked with DMRG
  EXPECT_NEAR(secondExcitedState, 15.256, 1e-3);
}

TEST_F(ACICalculation, CanCalculateCIProblemWithCalculatorAndDirectDavidson) {
  ci->setReferenceCalculator(calc);
  ci->settings().modifyString(SettingsNames::doublesIndicesOption, "");
  ci->settings().modifyBool(SettingsNames::casWithUnoCI, true);
  ci->settings().modifyInt(SettingsNames::casOrbitalsAroundFermiOption, 4);
  ci->settings().modifyInt(SettingsNames::referencesAroundFermiOption, 4);
  ci->settings().modifyInt(SettingsNames::numberOfEigenstates, 3);
  //  ci->settings().modifyDouble(convergence, 1e-9);
  auto resultWithMatrix = ci->calculate();
  auto res = resultWithMatrix.get<Property::ExcitedStates>().unrestricted->eigenStates;
  Eigen::VectorXd multiplicities = resultWithMatrix.get<Property::ExcitedStates>().unrestricted->multiplicities;

  ci->settings().modifyBool(SettingsNames::directness, true);

  auto resultDirect = ci->calculate();
  auto resDirect = resultWithMatrix.get<Property::ExcitedStates>().unrestricted->eigenStates;
  Eigen::VectorXd multiplicitiesDirect = resultDirect.get<Property::ExcitedStates>().unrestricted->multiplicities;

  EXPECT_EQ(resultWithMatrix.get<Property::ExcitedStates>().transitionLabels.size(), 36);
  EXPECT_EQ(resDirect.eigenValues.size(), 3);

  for (int state = 0; state < res.eigenValues.size(); ++state) {
    EXPECT_NEAR(res.eigenValues(state), resDirect.eigenValues(state), 1e-12);
    EXPECT_NEAR(multiplicities(state), multiplicitiesDirect(state), 1e-6);
    for (int det = 0; det < res.eigenVectors.rows(); ++det) {
      EXPECT_NEAR(std::abs(res.eigenVectors.col(state)(det)), std::abs(resDirect.eigenVectors.col(state)(det)), 1e-6);
    }
  }
}

TEST_F(ACICalculation, CanIndexOnv) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  const auto& specifier = calc->results().get<Property::SQSpecifier>();

  CISolver solver(ci->settings(), specifier);
  solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals);

  EXPECT_EQ(sve.onvToIndexMap(OccupationNumberVector{std::vector<int>{1, 1, 0, 0}}), 0);
  EXPECT_EQ(sve.onvToIndexMap(OccupationNumberVector{std::vector<int>{1, 0, 1, 0}}), 1);
  EXPECT_EQ(sve.onvToIndexMap(OccupationNumberVector{std::vector<int>{1, 0, 0, 1}}), 2);
  EXPECT_EQ(sve.onvToIndexMap(OccupationNumberVector{std::vector<int>{0, 1, 1, 0}}), 3);
  EXPECT_EQ(sve.onvToIndexMap(OccupationNumberVector{std::vector<int>{0, 1, 0, 1}}), 4);
  EXPECT_EQ(sve.onvToIndexMap(OccupationNumberVector{std::vector<int>{0, 0, 1, 1}}), 5);
}

TEST_F(ACICalculation, CanConstrucKroeneckerDelta) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  const auto& specifier = calc->results().get<Property::SQSpecifier>();

  CISolver solver(ci->settings(), specifier);
  solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals);
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(solver.basis());

  Eigen::MatrixXi kdMatrix;
  std::vector<std::vector<int>> kdDiagonal;
  std::tie(kdMatrix, kdDiagonal) = sve.constructKroeneckerDeltaMatrix(onvs);
  for (int onv1 = 0; onv1 < kdMatrix.rows(); ++onv1) {
    auto firstOnv = onvs[onv1];
    for (int onv2 = 0; onv2 < kdMatrix.rows(); ++onv2) {
      auto secondOnv = onvs[onv2];
      if (onv1 != onv2) {
        if (kdMatrix(onv1, onv2) != CISigmaVectorEvaluator::None) {
          int r = std::abs(kdMatrix(onv1, onv2)) / firstOnv.size();
          int s = std::abs(kdMatrix(onv1, onv2)) % firstOnv.size();
          auto expectedOnv = firstOnv;
          expectedOnv.flip(r);
          expectedOnv.flip(s);
          EXPECT_EQ(expectedOnv, secondOnv);
        }
      }
    }
  }

  ASSERT_EQ(kdDiagonal.size(), onvs.size());
  // diagonal elements
  for (int onv1 = 0; onv1 < kdMatrix.rows(); ++onv1) {
    auto firstOnv = onvs[onv1];
    for (int rs : kdDiagonal[onv1]) {
      int r = rs / (firstOnv.size() + 1);
      EXPECT_TRUE(firstOnv.isOccupied(r));
    }
  }
}

TEST_F(ACICalculation, CanCreateIndicesMap) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  const auto& specifier = calc->results().get<Property::SQSpecifier>();

  CISolver solver(ci->settings(), specifier);
  solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals);
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(solver.basis());

  // Test connection
  int idx = 0;
  for (const auto& onv : onvs) {
    for (int idxConnected : sve.connection()[idx]) {
      auto it = std::find(solver.basis().begin(), solver.basis().end(), ElectronicDeterminant({onv}, {onvs[idxConnected]}));
      ASSERT_TRUE(it != solver.basis().end());
      ASSERT_EQ(it->getOnv(SpinComponent::Alpha), onv);
      it = std::find(solver.basis().begin(), solver.basis().end(), ElectronicDeterminant({onvs[idxConnected]}, {onv}));
      ASSERT_TRUE(it != solver.basis().end());
      ASSERT_EQ(it->getOnv(SpinComponent::Beta), onv);
    }
    ++idx;
  }
  for (int detIdx = 0; detIdx < int(solver.basis().size()); ++detIdx) {
    const auto& det = solver.basis()[detIdx];
    const auto& onvAlpha = det.getOnv(SpinComponent::Alpha);
    const auto& onvBeta = det.getOnv(SpinComponent::Beta);

    const auto& connected = sve.connection()[sve.onvToIndexMap(onvAlpha)];
    EXPECT_TRUE(std::find(connected.begin(), connected.end(), sve.onvToIndexMap(onvBeta)) != connected.end());
  }

  // test alpha/betaIndices
  for (int onv1 = 0; onv1 < int(onvs.size()); ++onv1) {
    for (int onv2 = 0; onv2 < int(sve.connection()[onv1].size()); ++onv2) {
      auto itAlpha =
          std::find(solver.basis().begin(), solver.basis().end(), ElectronicDeterminant({onvs[onv1]}, {onvs[onv2]}));
      auto itBeta =
          std::find(solver.basis().begin(), solver.basis().end(), ElectronicDeterminant({onvs[onv2]}, {onvs[onv1]}));
      EXPECT_EQ(std::distance(solver.basis().begin(), itAlpha), sve.indices(SpinComponent::Alpha)[onv1][onv2]);
      EXPECT_EQ(std::distance(solver.basis().begin(), itBeta), sve.indices(SpinComponent::Beta)[onv1][onv2]);
    }
  }
  for (int detIdx = 0; detIdx < int(solver.basis().size()); ++detIdx) {
    const auto& det = solver.basis()[detIdx];
    const auto& onvAlpha = det.getOnv(SpinComponent::Alpha);
    const auto& onvBeta = det.getOnv(SpinComponent::Beta);

    const auto& connectedA = sve.connection()[sve.onvToIndexMap(onvAlpha)];
    const auto& connectedB = sve.connection()[sve.onvToIndexMap(onvBeta)];
    auto itA = std::find(connectedA.begin(), connectedA.end(), sve.onvToIndexMap(onvBeta));
    auto itB = std::find(connectedB.begin(), connectedB.end(), sve.onvToIndexMap(onvAlpha));

    ASSERT_TRUE(itA != connectedA.end());
    ASSERT_TRUE(itB != connectedB.end());
    EXPECT_EQ(sve.indices(SpinComponent::Alpha)[sve.onvToIndexMap(onvAlpha)][*itA], detIdx);
    EXPECT_EQ(sve.indices(SpinComponent::Beta)[sve.onvToIndexMap(onvBeta)][*itB], detIdx);
  }
}

TEST_F(ACICalculation, CanCalculateCIProblemWithCalculatorAndDirectDavidsonCIS) {
  ci->setReferenceCalculator(calc);
  ci->settings().modifyString(SettingsNames::doublesIndicesOption, "");
  ci->settings().modifyBool(SettingsNames::casWithUnoCI, true);
  ci->settings().modifyBool(SettingsNames::directness, false);
  ci->settings().modifyInt(SettingsNames::casOrbitalsAroundFermiOption, 4);
  ci->settings().modifyInt(SettingsNames::referencesAroundFermiOption, 0);
  ci->settings().modifyInt(SettingsNames::numberOfEigenstates, 10);
  auto result = ci->calculate();
  auto res = result.get<Property::ExcitedStates>().unrestricted->eigenStates;
  Eigen::VectorXd multiplicities = result.get<Property::ExcitedStates>().unrestricted->multiplicities;

  ci->settings().modifyString(SettingsNames::doublesIndicesOption, "");
  ci->settings().modifyBool(SettingsNames::casWithUnoCI, true);
  ci->settings().modifyBool(SettingsNames::directness, true);
  ci->settings().modifyInt(SettingsNames::casOrbitalsAroundFermiOption, 4);
  ci->settings().modifyInt(SettingsNames::referencesAroundFermiOption, 0);
  ci->settings().modifyInt(SettingsNames::numberOfEigenstates, 10);
  ci->referenceCalculation();
  result = ci->calculate();
  auto resDirect = result.get<Property::ExcitedStates>().unrestricted->eigenStates;
  Eigen::VectorXd multiplicitiesDirect = result.get<Property::ExcitedStates>().unrestricted->multiplicities;

  for (int state = 0; state < res.eigenValues.size(); ++state) {
    EXPECT_NEAR(res.eigenValues(state), resDirect.eigenValues(state), 1e-12);
    EXPECT_NEAR(multiplicities(state), multiplicitiesDirect(state), 1e-6);
    for (int det = 0; det < res.eigenVectors.rows(); ++det) {
      EXPECT_NEAR(std::abs(res.eigenVectors.col(state)(det)), std::abs(resDirect.eigenVectors.col(state)(det)), 1e-6);
    }
  }
}

TEST_F(ACICalculation, CanCalculateK) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  auto specifier = calc->results().get<Property::SQSpecifier>();

  CISolver solver(ci->settings(), specifier);
  //  solver.setIndices(MOIndicesForDoubles{{0, 1, 2, 3}});
  // solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals, specifier.coreEnergy);
  auto basis = solver.basis();
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(basis);

  Eigen::MatrixXd k =
      *CISigmaVectorEvaluator::constructK(specifier.integrals.oneBodyIntegrals, specifier.integrals.twoBodyIntegrals);

  Eigen::MatrixXd k2 = specifier.integrals.get(MoIntegrals::Type::OneBody);
  int nOrbs = onvs[0].size();
  for (int p = 0; p < nOrbs; ++p) {
    for (int q = 0; q < nOrbs; ++q) {
      for (int r = 0; r < nOrbs; ++r) {
        k2(p, q) -= 0.5 * specifier.integrals.get(MoIntegrals::Type::TwoBody)(p * nOrbs + r, r * nOrbs + q);
      }
    }
  }
  ASSERT_EQ(k.rows(), k2.rows());
  ASSERT_EQ(k.cols(), k2.cols());
  for (int row = 0; row < k2.rows(); ++row) {
    for (int col = 0; col < k2.cols(); ++col) {
      ASSERT_NEAR(k(row, col), k2(row, col), 1e-12);
    }
  }
}

TEST_F(ACICalculation, CanCalculatePrecontractedK) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  auto specifier = calc->results().get<Property::SQSpecifier>();

  specifier.integrals.twoBodyIntegrals = {};
  HermitianHamiltonian hamiltonian(specifier);

  CISolver solver(ci->settings(), specifier);
  //  solver.setIndices(MOIndicesForDoubles{{0, 1, 2, 3}});
  // solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals, specifier.coreEnergy);
  auto basis = solver.basis();
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(basis);

  Eigen::MatrixXd k = *CISigmaVectorEvaluator::constructK(specifier.integrals.get(MoIntegrals::Type::OneBody),
                                                          specifier.integrals.twoBodyIntegrals);
  Eigen::MatrixXd preconK = *sve.precontractK(k);

  Eigen::MatrixXd preconK2 = Eigen::MatrixXd::Zero(onvs.size(), onvs.size());
  for (int row = 0; row < preconK2.rows(); ++row) {
    const auto& onvI = onvs[row];
    for (int col = 0; col < preconK2.cols(); ++col) {
      const auto& onvJ = onvs[col];
      for (int p = 0; p < onvs[0].size(); ++p) {
        if (onvI.isOccupied(p)) {
          for (int q = 0; q < onvs[0].size(); ++q) {
            if (!onvI.isOccupied(q) || p == q) {
              auto onvTry = onvI;
              onvTry.flip(p);
              onvTry.flip(q);
              if (onvTry == onvJ) {
                preconK2(row, col) += CISigmaVectorEvaluator::excitationSign(onvI, p, q) * k(p, q);
              }
            }
          }
        }
      }
    }
  }

  ASSERT_EQ(preconK.rows(), preconK2.rows());
  ASSERT_EQ(preconK.cols(), preconK2.cols());

  for (int row = 0; row < preconK2.rows(); ++row) {
    for (int col = 0; col < preconK2.cols(); ++col) {
      ASSERT_NEAR(preconK(row, col), preconK2(row, col), 1e-12);
    }
  }
}

TEST_F(ACICalculation, DirectDiagonalizesOneBody) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  auto specifier = calc->results().get<Property::SQSpecifier>();

  specifier.integrals.twoBodyIntegrals = {};
  HermitianHamiltonian hamiltonian(specifier);

  CISolver solver(ci->settings(), specifier);
  //  solver.setIndices(MOIndicesForDoubles{{0, 1, 2, 3}});
  // solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  auto sve = std::make_shared<CISigmaVectorEvaluator>(solver.basis(), specifier.integrals, specifier.coreEnergy);
  auto basis = solver.basis();
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(basis);

  Eigen::MatrixXd ciMatrix = Eigen::MatrixXd::Zero(basis.size(), basis.size());
  for (int row = 0; row < int(basis.size()); ++row) {
    for (int col = 0; col < int(basis.size()); ++col) {
      ciMatrix(row, col) = hamiltonian.calculateMatrixElement(basis[row], basis[col]);
    }
  }

  auto indsve = std::make_shared<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(ciMatrix);

  auto diagonalize = [&](std::shared_ptr<SigmaVectorEvaluator> sigmaVectorEval) -> EigenContainer {
    NonOrthogonalDavidson davidson(9, ciMatrix.rows());
    davidson.setSigmaVectorEvaluator(std::move(sigmaVectorEval));
    davidson.setPreconditionerEvaluator(std::make_shared<IndirectPreconditionerEvaluator>(ciMatrix.diagonal()));
    return davidson.solve(logger);
  };

  auto ensve = diagonalize(sve);
  auto enindsve = diagonalize(indsve);
  for (int state = 0; state < 9; ++state) {
    EXPECT_NEAR(ensve.eigenValues(state), enindsve.eigenValues(state), 1e-12);
    for (int det = 0; det < ensve.eigenVectors.rows(); ++det) {
      EXPECT_NEAR(std::abs(ensve.eigenVectors.col(state)(det)), std::abs(enindsve.eigenVectors.col(state)(det)), 1e-6);
    }
  }
}

TEST_F(ACICalculation, CanCalculatePrecontractedG) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  auto specifier = calc->results().get<Property::SQSpecifier>();

  HermitianHamiltonian hamiltonian(specifier);

  CISolver solver(ci->settings(), specifier);
  //  solver.setIndices(MOIndicesForDoubles{{0, 1, 2, 3}});
  // solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals, specifier.coreEnergy);
  auto basis = solver.basis();
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(basis);

  Eigen::MatrixXd g = specifier.integrals.get(MoIntegrals::Type::TwoBody);

  Eigen::MatrixXd preconG = sve.precontractG(g);

  Eigen::MatrixXd preconG2 = Eigen::MatrixXd::Zero(onvs.size(), onvs.size());
  for (int row = 0; row < preconG2.rows(); ++row) {
    const auto& onvI = onvs[row];
    std::vector<std::tuple<OccupationNumberVector, int, int>> single;
    for (int p = 0; p < onvs[0].size(); ++p) {
      if (onvI.isOccupied(p)) {
        for (int q = 0; q < onvs[0].size(); ++q) {
          if (!onvI.isOccupied(q) || p == q) {
            auto onvSingle = onvI;
            onvSingle.flip(p);
            onvSingle.flip(q);
            single.emplace_back(onvSingle, p, q);
          }
        }
      }
    }
    for (const auto& onv2 : single) {
      const auto& onv_2 = std::get<0>(onv2);
      int p = std::get<1>(onv2);
      int q = std::get<2>(onv2);
      for (int r = 0; r < onvs[0].size(); ++r) {
        if (onv_2.isOccupied(r)) {
          for (int s = 0; s < onvs[0].size(); ++s) {
            if (!onv_2.isOccupied(s) || r == s) {
              auto onvJ = onv_2;
              onvJ.flip(r);
              onvJ.flip(s);
              try {
                int idxJ = sve.onvToIndexMap(onvJ);
                preconG2(row, idxJ) += CISigmaVectorEvaluator::excitationSign(onvI, p, q) *
                                       CISigmaVectorEvaluator::excitationSign(onv_2, r, s) *
                                       g(p * onvs[0].size() + q, r * onvs[0].size() + s);
              }
              catch (...) {
                // do nothing, I am misusing this as a if found then do
              }
            }
          }
        }
      }
    }
  }

  ASSERT_EQ(preconG.rows(), preconG2.rows());
  ASSERT_EQ(preconG.cols(), preconG2.cols());

  for (int row = 0; row < preconG2.rows(); ++row) {
    for (int col = 0; col < preconG2.cols(); ++col) {
      ASSERT_NEAR(preconG(row, col), 0.5 * preconG2(row, col), 1e-12);
    }
  }
}

TEST_F(ACICalculation, DirectCalculatesSigma1Correctly) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  auto specifier = calc->results().get<Property::SQSpecifier>();

  specifier.integrals.twoBodyIntegrals = {};
  HermitianHamiltonian hamiltonian(specifier);

  CISolver solver(ci->settings(), specifier);
  //  solver.setIndices(MOIndicesForDoubles{{0, 1, 2, 3}});
  // solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals, specifier.coreEnergy);
  auto basis = solver.basis();
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(basis);

  Eigen::MatrixXd ciMatrix = Eigen::MatrixXd::Zero(basis.size(), basis.size());
  for (int row = 0; row < int(basis.size()); ++row) {
    for (int col = 0; col < int(basis.size()); ++col) {
      ciMatrix(row, col) = hamiltonian.calculateMatrixElement(basis[row], basis[col]);
    }
  }

  IndirectSigmaVectorEvaluator<Eigen::MatrixXd> indsve(ciMatrix);

  // Be careful! Sigmas will be equal only for vectors with 1 entry.
  // Since the determinants in direct ci and in hamiltonian have a different
  // representation
  // direct ci ->    |onv_alpha> x |onv_beta>
  // hamiltonian ->  |\phi_{1\alpha}\phi_{1\beta}\phi_{2\alpha}\phi_{3\beta}...>
  // the sign is not the same. Since a sigma vector is not an observable,
  // the different sign give different sigma vectors in the contraction.
  //
  // If a vector with just one entry is given as a guess, then the difference
  // is at most the sign in the sigma vector.
  Eigen::MatrixXd guess = Eigen::MatrixXd::Identity(ciMatrix.rows(), 9);

  Eigen::MatrixXd indsigma = indsve.evaluate(guess);
  Eigen::MatrixXd dirsigma = sve.evaluate(guess);

  for (int state = 0; state < indsigma.cols(); ++state) {
    for (int det = 0; det < indsigma.rows(); ++det) {
      EXPECT_NEAR(std::abs(indsigma.col(state)(det)), std::abs(dirsigma.col(state)(det)), 1e-6);
    }
  }
}

TEST_F(ACICalculation, DirectDiagonalizesTwoBody) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  auto specifier = calc->results().get<Property::SQSpecifier>();

  specifier.integrals.oneBodyIntegrals = {};
  HermitianHamiltonian hamiltonian(specifier);

  CISolver solver(ci->settings(), specifier);
  //  solver.setIndices(MOIndicesForDoubles{{0, 1, 2, 3}});
  // solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  auto sve = std::make_shared<CISigmaVectorEvaluator>(solver.basis(), specifier.integrals, specifier.coreEnergy);
  auto basis = solver.basis();
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(basis);

  Eigen::MatrixXd ciMatrix = Eigen::MatrixXd::Zero(basis.size(), basis.size());
  for (int row = 0; row < int(basis.size()); ++row) {
    for (int col = 0; col < int(basis.size()); ++col) {
      ciMatrix(row, col) = hamiltonian.calculateMatrixElement(basis[row], basis[col]);
    }
  }

  auto indsve = std::make_shared<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(ciMatrix);

  auto diagonalize = [&](std::shared_ptr<SigmaVectorEvaluator> sigmaVectorEval) -> EigenContainer {
    NonOrthogonalDavidson davidson(9, ciMatrix.rows());
    davidson.setSigmaVectorEvaluator(std::move(sigmaVectorEval));
    davidson.setPreconditionerEvaluator(std::make_shared<IndirectPreconditionerEvaluator>(ciMatrix.diagonal()));
    return davidson.solve(logger);
  };

  auto ensve = diagonalize(sve);
  auto enindsve = diagonalize(indsve);
  for (int state = 0; state < 9; ++state) {
    EXPECT_NEAR(ensve.eigenValues(state), enindsve.eigenValues(state), 1e-12);
    for (int det = 0; det < ensve.eigenVectors.rows(); ++det) {
      ASSERT_NEAR(std::abs(ensve.eigenVectors.col(state)(det)), std::abs(enindsve.eigenVectors.col(state)(det)), 1e-6);
    }
  }
}
TEST_F(ACICalculation, DirectCalculatesSigma2Correctly) {
  ci->setReferenceCalculator(calc);
  ci->referenceCalculation();
  auto specifier = calc->results().get<Property::SQSpecifier>();

  specifier.integrals.oneBodyIntegrals = {};
  HermitianHamiltonian hamiltonian(specifier);

  CISolver solver(ci->settings(), specifier);
  //  solver.setIndices(MOIndicesForDoubles{{0, 1, 2, 3}});
  // solver.setIndices(MOIndicesForReferences{{1, 2}});
  solver.solve(logger);
  CISigmaVectorEvaluator sve(solver.basis(), specifier.integrals, specifier.coreEnergy);
  auto basis = solver.basis();
  auto onvs = CISigmaVectorEvaluator::getSortedOnvList(basis);

  Eigen::MatrixXd ciMatrix = Eigen::MatrixXd::Zero(basis.size(), basis.size());
  for (int row = 0; row < int(basis.size()); ++row) {
    for (int col = 0; col < int(basis.size()); ++col) {
      ciMatrix(row, col) = hamiltonian.calculateMatrixElement(basis[row], basis[col]);
    }
  }

  IndirectSigmaVectorEvaluator<Eigen::MatrixXd> indsve(ciMatrix);

  // Be careful! Sigmas will be equal only for vectors with 1 entry.
  // Since the determinants in direct ci and in hamiltonian have a different
  // representation
  // direct ci ->    |onv_alpha> x |onv_beta>
  // hamiltonian ->  |\phi_{1\alpha}\phi_{1\beta}\phi_{2\alpha}\phi_{3\beta}...>
  // the sign is not the same. Since a sigma vector is not an observable,
  // the different sign give different sigma vectors in the contraction.
  //
  // If a vector with just one entry is given as a guess, then the difference
  // is at most the sign in the sigma vector.
  Eigen::MatrixXd guess = Eigen::MatrixXd::Identity(ciMatrix.rows(), 9);

  Eigen::MatrixXd indsigma = indsve.evaluate(guess);
  Eigen::MatrixXd dirsigma = sve.evaluate(guess);

  for (int state = 0; state < indsigma.cols(); ++state) {
    for (int det = 0; det < indsigma.rows(); ++det) {
      EXPECT_NEAR(std::abs(indsigma.col(state)(det)), std::abs(dirsigma.col(state)(det)), 1e-6);
    }
  }
}

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
