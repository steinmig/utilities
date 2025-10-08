/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CICalculator.h"
#include "CICalculatorSettings.h"
#include "CISolver.h"
#include "CasIndicesHandler.h"
#include "SpinSquaredEvaluator.h"
#include "Utils/DataStructures/SecondQuantization/CasSpecifier.h"
#include "Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h"
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
namespace Scine {
namespace Utils {
namespace SecondQuantization {

CICalculator::CICalculator()
  : settings_(std::make_unique<CICalculatorSettings>()), results_(std::make_unique<Results>()) {
}

CICalculator::~CICalculator() = default;
void CICalculator::setReferenceCalculator(std::shared_ptr<Core::Calculator> referenceCalculator) {
  if (!referenceCalculator->possibleProperties().containsSubSet(Utils::Property::SQSpecifier)) {
    throw std::runtime_error("The reference calculator does not provide a Hamiltonian in second quantization.");
  }
  calc_ = std::move(referenceCalculator);
}

void CICalculator::referenceCalculation() {
  checkReferenceCalculator();
  applySettings();
  calc_->setRequiredProperties(Utils::Property::SQSpecifier | Utils::Property::DipoleMatrixAO);
  calc_->calculate("");
  MOIndicesForReferences refIndices;
  MOIndicesForDoubles doublesIndices;
  std::tie(refIndices, doublesIndices) = getIndices();
  if (!calc_->results().has<Property::SQSpecifier>()) {
    throw std::runtime_error("No SecondQuantization Specifier found in the reference calculation results.");
  }

  // Create CI solver and set indices.
  solver_ = std::make_unique<CISolver>(settings(), calc_->results().get<Property::SQSpecifier>());
  if (!refIndices.data.empty()) {
    solver_->setIndices(refIndices);
  }
  if (!doublesIndices.data.empty()) {
    solver_->setIndices(doublesIndices);
  }
  solver_->setReferenceFromDoubleExcitations(settings().getBool(SettingsNames::referenceFromDoubleExcitations));
}

Core::Calculator& CICalculator::getReferenceCalculator() {
  return *calc_;
}
const Core::Calculator& CICalculator::getReferenceCalculator() const {
  return *calc_;
}
const Results& CICalculator::calculate() {
  if (!calc_) {
    throw std::runtime_error("No reference calculator was set.");
  }
  // ALWAYS do a reference calculation, otherwise one could silently
  // always get the same stuff after changing settings without wanting this.
  referenceCalculation();
  eigenpairs_ = std::make_unique<EigenContainer>(solver_->solve(getLog()));

  Eigen::VectorXd multiplicity = calculateSpinMultiplicity(*eigenpairs_);
  Eigen::Matrix3Xd transitionDipoles = calculateTransitionDipoles(*eigenpairs_);
  // SpinAdaptedElectronicTransitionResult temp(std::make_shared<ElectronicTransitionResult>(ElectronicTransitionResult{
  //     *eigenpairs_, std::move(transitionDipoles), std::move(multiplicity)}),
  // std::shared_ptr<ElectronicTransitionResult>(),
  // std::shared_ptr<ElectronicTransitionResult>(),
  // getLabels(solver_->basis()), std::shared_ptr<ElectronicCIResult>());
  results_->set<Property::ExcitedStates>(SpinAdaptedElectronicTransitionResult(
      std::make_shared<ElectronicTransitionResult>(
          ElectronicTransitionResult{*eigenpairs_, std::move(transitionDipoles), std::move(multiplicity)}),
      std::shared_ptr<ElectronicTransitionResult>(), std::shared_ptr<ElectronicTransitionResult>(),
      getLabels(solver_->basis()), std::shared_ptr<ElectronicCIResult>()));

  return *results_;
}

auto CICalculator::getWavefunction(int state) const -> SpinSquaredEvaluation::Wavefunction {
  if (!eigenpairs_ || state > eigenpairs_->eigenValues.size()) {
    throw std::runtime_error("Required state has not been calculated.");
  }

  SpinSquaredEvaluation::Wavefunction wavefunction;
  for (unsigned i = 0; i < solver_->basis().size(); ++i) {
    wavefunction[solver_->basis()[i]] = eigenpairs_->eigenVectors.col(state)(i);
  }
  return wavefunction;
}

auto CICalculator::calculateSpinMultiplicity(const EigenContainer& results) const -> Eigen::VectorXd {
  Eigen::VectorXd multiplicity(results.eigenValues.size());

#pragma omp parallel for
  for (int state = 0; state < results.eigenValues.size(); ++state) {
    multiplicity(state) = std::sqrt(1 + 4.0 * SpinSquaredEvaluation::evaluate(getWavefunction(state)));
  }

  return multiplicity;
}

auto CICalculator::calculateTransitionDipoles(const EigenContainer& results) const -> Eigen::Matrix3Xd {
  if (!calc_->results().has<Utils::Property::DipoleMatrixAO>() || !settings().getBool("calculate_transition_dipole")) {
    return Eigen::Matrix3Xd::Zero(3, results.eigenValues.size());
  }
  Eigen::Matrix3Xd transitionDipole(3, results.eigenValues.size());

  const auto& dipoleMatAO = calc_->results().get<Utils::Property::DipoleMatrixAO>();
  const Eigen::MatrixXd cas =
      calc_->results().get<Property::SQSpecifier>().cas.getActiveOrbitalsCoefficients().restrictedMatrix();
  Eigen::MatrixXd dipMatMOx = cas.transpose() * dipoleMatAO[0] * cas;
  Eigen::MatrixXd dipMatMOy = cas.transpose() * dipoleMatAO[1] * cas;
  Eigen::MatrixXd dipMatMOz = cas.transpose() * dipoleMatAO[2] * cas;

  auto evaluateTransitionDipole = [this, &dipMatMOx, dipMatMOy,
                                   dipMatMOz](const SpinSquaredEvaluation::Wavefunction& groundState,
                                              const SpinSquaredEvaluation::Wavefunction& state) -> Eigen::Vector3d {
    SpinSquaredEvaluation::Wavefunction x, y, z;
    for (const auto& det : state) {
      for (SpinComponent spin : {SpinComponent::Alpha, SpinComponent::Beta}) {
        for (int p : det.first.getAllOccupied(spin)) {
          for (int q : det.first.getAllVirtual(spin)) {
            auto tmp = det.first;
            tmp.applyExcitation(SingleExcitation(p, q, spin));
            x[tmp] += solver_->getExcitationSign(det.first, SingleExcitation(p, q, spin)) * dipMatMOx(p, q) * det.second;
            y[tmp] += solver_->getExcitationSign(det.first, SingleExcitation(p, q, spin)) * dipMatMOy(p, q) * det.second;
            z[tmp] += solver_->getExcitationSign(det.first, SingleExcitation(p, q, spin)) * dipMatMOz(p, q) * det.second;
          }
        }
      }
    }
    return {SpinSquaredEvaluation::overlap(groundState, x), SpinSquaredEvaluation::overlap(groundState, y),
            SpinSquaredEvaluation::overlap(groundState, z)};
  };

  const auto groundState = getWavefunction(0);
#pragma omp parallel for
  for (int state = 0; state < results.eigenValues.size(); ++state) {
    transitionDipole.col(state) = evaluateTransitionDipole(groundState, getWavefunction(state));
  }

  return transitionDipole;
}

std::string CICalculator::name() const {
  return "CICalculator";
}

Settings& CICalculator::settings() {
  return *settings_;
}
const Settings& CICalculator::settings() const {
  return *settings_;
}
void CICalculator::applySettings() {
  for (const auto& keyValuePair : settings()) {
    if (calc_->settings().valueExists(keyValuePair.first) && keyValuePair.first != SettingsNames::unoCiLowThresholdOption) {
      calc_->settings().modifyValue(keyValuePair.first, keyValuePair.second);
    }
  }
}

Results& CICalculator::results() {
  return *results_;
}
const Results& CICalculator::results() const {
  return *results_;
}

void CICalculator::checkReferenceCalculator() {
  if (!calc_) {
    throw std::runtime_error("No reference calculation present.");
  }
  // TODO: Check settings
}
auto CICalculator::getIndices() const -> std::tuple<MOIndicesForReferences, MOIndicesForDoubles> {
  if (calc_->results().has<Property::SQSpecifier>()) {
    CasIndicesHandler handler(calc_->results().get<Property::SQSpecifier>().cas);

    auto doubleIndicesString = settings().getString(SettingsNames::doublesIndicesOption);
    std::vector<int> doubles;
    if (doubleIndicesString.empty()) {
      doubles = handler.getIndices(settings().getInt(SettingsNames::doublesAroundFermiOption));
    }
    else {
      doubles = handler.getIndices(CasGenerator::parseCasString(doubleIndicesString));
    }

    return {getReferenceIndices(handler), {doubles}};
  }
  throw std::runtime_error("No second quantization specifier in results.");
}

auto CICalculator::getReferenceIndices(const CasIndicesHandler& handler) const -> MOIndicesForReferences {
  auto referenceIndicesString = settings().getString(SettingsNames::referencesIndicesOption);
  std::vector<int> references;
  if (!referenceIndicesString.empty()) {
    references = handler.getIndices(CasGenerator::parseCasString(referenceIndicesString));
  }
  else if (settings().getInt(SettingsNames::referencesAroundFermiOption) != ReferenceSpecifiers::SingleReference) {
    references = handler.getIndices(settings().getInt(SettingsNames::referencesAroundFermiOption));
  }
  else if (settings().getBool(SettingsNames::casWithUnoCI)) {
    Eigen::VectorXd occs = calc_->results().get<Property::SQSpecifier>().occupations;
    auto activeIndices = calc_->results().get<Property::SQSpecifier>().cas.getActiveIndices().restricted;
    double threshold = settings().getDouble(SettingsNames::unoCiLowThresholdOption);
    for (int i = 0; i < int(occs.size()); ++i) {
      // numerical errors correction
      if (occs(i) > 2) {
        occs(i) = 2;
      }
      if (occs(i) < 0) {
        occs(i) = 0;
      }
      // add to reference set if UNO occupation between threshold and 2-threshold
      if (occs(i) > threshold && occs(i) < 2 - threshold) {
        references.push_back(std::distance(activeIndices.begin(), std::find(activeIndices.begin(), activeIndices.end(), i)));
      }
    }
  }
  return {references};
}

auto CICalculator::getLabels(const std::vector<ElectronicDeterminant>& basis) -> std::vector<std::string> {
  std::vector<std::string> labels;
  labels.reserve(basis.size());

  for (const auto& det : basis) {
    std::string label(det.size(SpinComponent::Alpha), '0');
    for (int i = 0; i < int(label.size()); ++i) {
      if (det.isOccupied(SpinComponent::Alpha, i) && det.isOccupied(SpinComponent::Beta, i)) {
        label[i] = '2';
      }
      else if (det.isOccupied(SpinComponent::Alpha, i)) {
        label[i] = 'a';
      }
      else if (det.isOccupied(SpinComponent::Beta, i)) {
        label[i] = 'b';
      }
    }
    labels.emplace_back(std::move(label));
  }

  return labels;
}

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
