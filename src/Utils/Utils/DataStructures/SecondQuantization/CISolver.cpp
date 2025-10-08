/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CISolver.h"
#include "CICalculatorSettings.h"
#include "CIDeterminantsCreator.h"
#include <Core/Log.h>
#include <Utils/DataStructures/SecondQuantization/CISigmaVectorEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Math/IterativeDiagonalizer/DiagonalizerSettings.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectPreconditionerEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectSigmaVectorEvaluator.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

namespace logging {
template<typename T>
void logTiming(Core::Log& log, const std::vector<T>& times) {
  auto time = [](auto t1, auto t0) -> double {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  };
  log.output << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(40) << ""
             << "CI TIMINGS [ms]" << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(5) << "" << std::setw(35) << "Diagonal Elements" << std::setw(35)
             << "Sorting of the Basis" << std::setw(35) << "CI Matrix Calculation" << Core::Log::nl;
  log.output << std::left << std::setw(5) << "" << std::setw(35) << time(times[1], times[0]) << std::setw(35)
             << time(times[2], times[1]) << std::setw(35) << time(times[3], times[2]) << Core::Log::endl;
  log.output << std::right;
}

void logActiveSpaceInformation(Core::Log& log, const ElectronicDeterminant& hfDet,
                               boost::optional<MOIndicesForReferences> refsMo,
                               boost::optional<MOIndicesForDoubles> doublesMo, bool refsFromDouble) {
  log.output << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(50) << ""
             << "ACTIVE SPACE INFORMATION" << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(5) << ""
             << "Hartree--Fock Determinant:" << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::setw(20) << "";
  log.output << [&](std::ostream& logger) { hfDet.print(logger); };
  log.output << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(5) << ""
             << "AS for generation of reference determinants:" << Core::Log::nl;
  log.output << Core::Log::nl;
  if (!refsMo) {
    log.output << std::setw(20) << ""
               << "Single reference";
  }
  else {
    log.output << std::setw(20) << "";
    log.output << [&](std::ostream& logger) {
      std::copy(refsMo->data.begin(), refsMo->data.end(), std::ostream_iterator<int>(logger, " "));
    };
    log.output << Core::Log::nl;
    if (refsFromDouble) {
      log.output << std::setw(20) << ""
                 << "References from symmetrical double excitations";
    }
  }
  log.output << Core::Log::nl;

  log.output << std::left << std::setw(5) << ""
             << "AS for generation of doubles:" << Core::Log::nl;
  log.output << Core::Log::nl;
  if (!doublesMo) {
    log.output << std::setw(20) << ""
               << "Only singly excited determinants";
  }
  else {
    log.output << std::setw(20) << "";
    log.output << [&](std::ostream& logger) {
      std::copy(doublesMo->data.begin(), doublesMo->data.end(), std::ostream_iterator<int>(logger, " "));
    };
  }
  log.output << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(5) << ""
             << "AS for generation of singles:" << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::setw(20) << ""
             << "Space spanned by Hartree--Fock determinant" << Core::Log::endl;
  log.output << std::right;
}

void logReferences(Core::Log& log, const std::vector<ElectronicDeterminant>& refs) {
  log.output << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(50) << ""
             << "REFERENCES" << Core::Log::nl;
  log.output << Core::Log::nl;
  for (const auto& ref : refs) {
    log.output << std::setw(20) << "";
    log.output << [&ref](std::ostream& logger) { ref.print(logger); };
    log.output << Core::Log::nl;
  }
}
void logBasis(Core::Log& log, const std::vector<ElectronicDeterminant>& basis) {
  log.output << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(50) << ""
             << "BASIS" << Core::Log::nl;
  log.output << Core::Log::nl;
  log.output << std::left << std::setw(5) << ""
             << "Dimension: " << basis.size() << Core::Log::nl;
}
} // namespace logging

CISolver::CISolver(const Settings& settings, const Specifier& specifier)
  : settings_(settings), integrals_(specifier.integrals), hamiltonian_(specifier), coreEnergy_(specifier.coreEnergy) {
  int nOrbs = specifier.cas.getActiveIndices().restricted.size();
  int nElectrons = specifier.cas.getCasElectrons().restricted;
  std::vector<int> onv(nOrbs, 0);
  std::fill(onv.begin(), onv.begin() + nElectrons / 2, 1);
  hfDet_ = ElectronicDeterminant(onv, onv);
}
void CISolver::setIndices(MOIndicesForReferences moIndicesReferences) {
  moIndicesReferences_ = std::move(moIndicesReferences);
}
void CISolver::setIndices(MOIndicesForDoubles moIndicesDoubles) {
  moIndicesDoubles_ = std::move(moIndicesDoubles);
}
void CISolver::setReferenceFromDoubleExcitations(bool referenceFromDoubles) {
  onlyDoubleExcitations_ = referenceFromDoubles;
}

auto CISolver::solve(Core::Log& log) -> EigenContainer {
  logging::logActiveSpaceInformation(log, hfDet_, moIndicesReferences_, moIndicesDoubles_, onlyDoubleExcitations_);
  //  DETERMINE REFS FOR MRCI OR JUST HF - DET
  std::vector<ElectronicDeterminant> refs = generateReferenceDeterminants();

  logging::logReferences(log, refs);
  // generate all connected determinants according to the options (CIS/CISD/MRCIS/MRCISD)
  if (basis_.empty()) {
    basis_ = generateCIBasis(refs);
  }

  logging::logBasis(log, basis_);

  auto t0 = std::chrono::system_clock::now();

  auto diagonal = getDiagonal();

  auto t1 = std::chrono::system_clock::now();
  // order
  std::sort(basis_.begin(), basis_.end(),
            [&](const auto& det1, const auto& det2) { return diagonal.at(det1) < diagonal.at(det2); });
  Eigen::VectorXd ciMatDiag(basis_.size());
  for (int i = 0; i < int(basis_.size()); ++i) {
    ciMatDiag(i) = diagonal.at(basis_[i]);
  }

  auto t2 = std::chrono::system_clock::now();

  if (settings_.getBool(SettingsNames::directness)) {
    auto sve = std::make_shared<CISigmaVectorEvaluator>(basis_, integrals_, coreEnergy_);
    auto t3 = std::chrono::system_clock::now();
    logging::logTiming(log, std::vector<decltype(t0)>{t0, t1, t2, t3});
    auto res = diagonalize(log, sve, ciMatDiag);
    convertSigns(res);
    return res;
  }
  // Fill ciMatrix
  Eigen::MatrixXd ciMat = Eigen::MatrixXd::Zero(basis_.size(), basis_.size());
  ciMat.diagonal() = ciMatDiag;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < int(basis_.size()); ++i) {
    for (int j = 0; j < i; ++j) {
      ciMat(i, j) = hamiltonian_.calculateMatrixElement(basis_[i], basis_[j]);
    }
  }

  auto t3 = std::chrono::system_clock::now();

  logging::logTiming(log, std::vector<decltype(t0)>{t0, t1, t2, t3});
  return diagonalize(
      log, std::make_shared<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(ciMat.selfadjointView<Eigen::Lower>()), ciMatDiag);
}

auto CISolver::generateCIBasis(const std::vector<ElectronicDeterminant>& references) const
    -> std::vector<ElectronicDeterminant> {
  std::vector<ElectronicDeterminant> basis = references;

  for (const auto& ref : references) {
    // generate list of double excitations from a reference
    auto lists = getOccupationLists(ref);
    auto connected = lists ? hamiltonian_.generateAllConnected(ref, (*lists)[0], (*lists)[1], true)
                           : hamiltonian_.generateSingleExcitations(ref, true);
    basis.reserve(connected.size() + basis.size());
    basis.insert(basis.end(), std::make_move_iterator(connected.begin()), std::make_move_iterator(connected.end()));
  }

  // Remove duplicates
  std::sort(basis.begin(), basis.end(), std::greater<>());
  auto last = std::unique(basis.begin(), basis.end());
  basis.erase(last, basis.end());

  return basis;
}

auto CISolver::getDiagonal() const -> std::unordered_map<ElectronicDeterminant, double, ElectronicDeterminantHash> {
  std::unordered_map<ElectronicDeterminant, double, ElectronicDeterminantHash> diagonal;
  // initialize elements
  for (const auto& det : basis_) {
    diagonal[det] = 0.0;
  }
#pragma omp parallel for
  for (int i = 0; i < int(basis_.size()); ++i) {
    double element = hamiltonian_.calculateMatrixElement({}, basis_[i].getAllOccupied(SpinComponent::Alpha),
                                                         basis_[i].getAllOccupied(SpinComponent::Beta));
    diagonal.at(basis_[i]) = element;
  }
  return diagonal;
}

const std::vector<ElectronicDeterminant>& CISolver::basis() const {
  return basis_;
}

std::vector<ElectronicDeterminant>& CISolver::basis() {
  return basis_;
}

auto CISolver::diagonalize(Core::Log& log, std::shared_ptr<SigmaVectorEvaluator> sve, const Eigen::VectorXd& diagonal)
    -> EigenContainer {
  int numberOfRoots = settings_.getInt(Utils::SettingsNames::numberOfEigenstates);
  int initialSubspaceDimension = settings_.getInt(Utils::SettingsNames::initialSubspaceDimension);
  checkAndCorrectNumberOfRoots(numberOfRoots, initialSubspaceDimension, diagonal.size());

  NonOrthogonalDavidson davidson(numberOfRoots, diagonal.size());

  davidson.settings().modifyInt(initialGuessDimensionOption, initialSubspaceDimension);
  davidson.settings().modifyDouble(residualNormToleranceOption, settings_.getDouble(convergence));
  davidson.settings().modifyString(gepAlgorithmForBalancedMethodOption,
                                   settings_.getString(gepAlgorithmForBalancedMethodOption));
  if (settings_.getInt(Utils::SettingsNames::maxDavidsonIterations) != 0) {
    davidson.settings().modifyInt(Utils::SettingsNames::maxDavidsonIterations,
                                  settings_.getInt(Utils::SettingsNames::maxDavidsonIterations));
  }
  davidson.setSigmaVectorEvaluator(std::move(sve));
  davidson.setPreconditionerEvaluator(std::make_shared<IndirectPreconditionerEvaluator>(diagonal));

  return davidson.solve(log);
}

void CISolver::convertSigns(EigenContainer& result) const {
  int idx = 0;
  for (const auto& det : basis_) {
    const auto& alphaOnv = det.getOccupation(SpinComponent::Alpha);
    int phase = 1;
    for (auto beta : det.getAllOccupied(SpinComponent::Beta)) {
      phase *= std::count(alphaOnv.begin() + beta + 1, alphaOnv.end(), true) % 2 == 0 ? 1 : -1;
    }
    result.eigenVectors.row(idx++).array() *= phase;
  }
}

auto CISolver::generateReferenceDeterminants() const -> std::vector<ElectronicDeterminant> {
  if (moIndicesReferences_) {
    CIDeterminantsCreator referenceCreator(hfDet_);
    return referenceCreator.generate(moIndicesReferences_->data, onlyDoubleExcitations_);
  }
  return {hfDet_};
}

auto CISolver::getOccupationLists(const ElectronicDeterminant& reference) const
    -> boost::optional<std::array<std::map<SpinComponent, std::vector<int>>, 2>> {
  // If no moIndicesDoubles_ specified, then only do CIS by returning no doubles indices
  if (!moIndicesDoubles_) {
    return {};
  }

  std::map<SpinComponent, std::vector<int>> occ = {{SpinComponent::Alpha, {}}, {SpinComponent::Beta, {}}};
  std::for_each(occ.begin(), occ.end(), [&](auto& element) { element.second.reserve(moIndicesDoubles_->data.size()); });
  auto vir = occ;

  for (int orbital : moIndicesDoubles_->data) {
    for (SpinComponent spin : {SpinComponent::Alpha, SpinComponent::Beta}) {
      if (reference.isOccupied(spin, orbital)) {
        occ.at(spin).push_back(orbital);
      }
      else {
        vir.at(spin).push_back(orbital);
      }
    }
  }
  return {{occ, vir}};
}

void CISolver::checkAndCorrectNumberOfRoots(int& numberOfEnergyLevels, int& initialSubspaceDimension, int nConfigurations) {
  // If 0 is given, then give all the energy levels
  if (numberOfEnergyLevels == 0 || numberOfEnergyLevels > nConfigurations) {
    numberOfEnergyLevels = nConfigurations;
  }
  if (initialSubspaceDimension == 0 || initialSubspaceDimension < numberOfEnergyLevels) {
    initialSubspaceDimension = numberOfEnergyLevels;
  }
  if (initialSubspaceDimension > nConfigurations) {
    initialSubspaceDimension = nConfigurations;
  }
}

auto CISolver::getExcitationSign(const ElectronicDeterminant& det, const SingleExcitation& excitation) const -> double {
  return hamiltonian_.calculateExcitationSign(det.getAllOccupied(SpinComponent::Alpha),
                                              det.getAllOccupied(SpinComponent::Beta), excitation);
}

CISolver::~CISolver() = default;
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
