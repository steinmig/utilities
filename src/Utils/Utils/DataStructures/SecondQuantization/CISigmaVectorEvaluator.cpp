/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CISigmaVectorEvaluator.h"
#include "ElectronicDeterminant.h"
#include "SQSpecifier.h"
#include <Eigen/Sparse>
#include <iomanip>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

constexpr const int CISigmaVectorEvaluator::None;

inline auto fillTwoBodyIntegrals(int nOrbs, const MoIntegrals& integrals) -> boost::optional<Eigen::MatrixXd> {
  if (integrals.twoBodyIntegrals) {
    if (integrals.get(MoIntegrals::Type::TwoBody).size() != std::pow(nOrbs, 4)) {
      throw std::runtime_error("Wrong dimension of two body tensor in CI sigma vector calculator.");
    }
    return integrals.twoBodyIntegrals;
  }
  return {};
}

CISigmaVectorEvaluator::CISigmaVectorEvaluator(const std::vector<ElectronicDeterminant>& basis,
                                               const MoIntegrals& integrals, double coreEnergy)
  : basis_(basis), coreEnergy_(coreEnergy) {
  onvs_ = getSortedOnvList(basis_);
  int index = 0;
  for (const auto& onv : onvs_) {
    indexMap_[onv] = index++;
  }

  constructIndicesMaps();
  std::tie(kroeneckerDeltaMatrix_, diagonalKDElements_) = constructKroeneckerDeltaMatrix(onvs_);

  g_ = fillTwoBodyIntegrals(onvs_[0].size(), integrals);
  precontractedK_ = precontractK(constructK(integrals.oneBodyIntegrals, g_));
  if (g_) {
    precontractedG_ = precontractG(*g_);
  }
}

auto CISigmaVectorEvaluator::evaluate(const Eigen::MatrixXd& guessVectors) const -> const Eigen::MatrixXd& {
  assert(guessVectors.rows() == Eigen::Index(basis_.size()));

  const int dimCol = guessVectors.cols();
  const int dimRow = guessVectors.rows();

  const int alreadyComputedSigmaVectors = sigmaVectors_.cols();
  const int vectorsToCompute = dimCol - alreadyComputedSigmaVectors;

  assert(vectorsToCompute > 0);
  Eigen::MatrixXd sigmaMatrix(dimRow, vectorsToCompute);

  for (int state = 0; state < vectorsToCompute; ++state) {
    int const colIndex = state + alreadyComputedSigmaVectors;
    const Eigen::VectorXd& guessVector = guessVectors.col(colIndex);

    Eigen::VectorXd sigma = coreEnergy_ * guessVector.array();

    if (precontractedK_) {
      formSigma1(sigma, guessVector);
    }
    if (g_) {
      formSigma2(sigma, guessVector);
    }
    sigmaMatrix.col(state) = sigma;
  }
  // Cache already calculated stuff
  sigmaVectors_.conservativeResize(dimRow, dimCol);
  sigmaVectors_.rightCols(vectorsToCompute) = sigmaMatrix;

  return sigmaVectors_;
}

void CISigmaVectorEvaluator::collapsed(int /*newSubspaceDimension*/) {
  sigmaVectors_ = Eigen::MatrixXd(0, 0);
}

auto CISigmaVectorEvaluator::getSortedOnvList(const std::vector<ElectronicDeterminant>& basis)
    -> std::vector<OccupationNumberVector> {
  std::vector<OccupationNumberVector> onvs;
  onvs.reserve(basis.size());

  // If ONV alpha =/= ONV beta, add also for spin component beta
  for (const auto& det : basis) {
    onvs.push_back(det.getOnv(SpinComponent::Alpha));
  }

  std::sort(onvs.begin(), onvs.end(),
            [](const auto& onv1, const auto& onv2) { return onv1.getOccupation() > onv2.getOccupation(); });
  auto last = std::unique(onvs.begin(), onvs.end());
  onvs.erase(last, onvs.end());

  return onvs;
}

void CISigmaVectorEvaluator::constructIndicesMaps() {
  detToOnvs_.resize(basis_.size());
  alphaIndices_.resize(onvs_.size());
  betaIndices_.resize(onvs_.size());
  onvConnection_.resize(onvs_.size());
#pragma omp parallel for schedule(dynamic)
  for (int onvIdx = 0; onvIdx < int(onvs_.size()); ++onvIdx) {
    const auto& alphaOnv = onvs_[onvIdx];
    std::vector<int> currentAlpha;
    std::vector<int> currentBeta;
    std::vector<int> connection;
    currentAlpha.reserve(onvs_.size());
    currentBeta.reserve(onvs_.size());
    connection.reserve(onvs_.size());

    int idx = 0;
    for (const auto& betaOnv : onvs_) {
      auto it = std::find(basis_.begin(), basis_.end(), ElectronicDeterminant{{alphaOnv}, {betaOnv}});
      auto it2 = std::find(basis_.begin(), basis_.end(), ElectronicDeterminant{{betaOnv}, {alphaOnv}});
      if (it != basis_.end()) {
        currentAlpha.push_back(std::distance(basis_.begin(), it));
        currentBeta.push_back(std::distance(basis_.begin(), it2));
        connection.push_back(idx);
      }
      idx++;
    }
    alphaIndices_[onvIdx] = std::move(currentAlpha);
    betaIndices_[onvIdx] = std::move(currentBeta);
    onvConnection_[onvIdx] = std::move(connection);
  }
  std::transform(
      basis_.begin(), basis_.end(), detToOnvs_.begin(), [this](const ElectronicDeterminant& det) -> std::pair<int, int> {
        return {indexMap_.at(det.getOnv(SpinComponent::Alpha)), indexMap_.at(det.getOnv(SpinComponent::Beta))};
      });
}

auto CISigmaVectorEvaluator::constructKroeneckerDeltaMatrix(const std::vector<OccupationNumberVector>& onvs) const
    -> std::tuple<Eigen::MatrixXi, std::vector<std::vector<int>>> {
  Eigen::MatrixXi kdMatrix = Eigen::MatrixXi::Constant(onvs.size(), onvs.size(), None);
  std::vector<std::vector<int>> diagonalKDIndices;
  diagonalKDIndices.reserve(onvs.size());

  for (long onvIdx = 0; onvIdx < long(onvs.size()); ++onvIdx) {
    const auto& onv = onvs[onvIdx];
    std::vector<int> indices;
    indices.reserve(onv.size());
    for (int occ = 0; occ < onv.size(); ++occ) {
      if (onv.isOccupied(occ)) {
        // diagonal elements
        indices.push_back(occ * (onv.size() + 1));

        for (int vir = 0; vir < onv.size(); ++vir) {
          if (!onv.isOccupied(vir)) {
            auto onv2 = onv;
            onv2.flip(occ);
            onv2.flip(vir);

            auto it = indexMap_.find(onv2);
            if (it != indexMap_.end()) {
              assert(kdMatrix(onvIdx, it->second) == None);
              kdMatrix(onvIdx, it->second) = excitationSign(onv, occ, vir) * (occ * onv.size() + vir);
            }
          }
        }
      }
    }
    diagonalKDIndices.emplace_back(std::move(indices));
  }

  return {kdMatrix, diagonalKDIndices};
}

void CISigmaVectorEvaluator::formSigma1(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector) const {
  contractSameSpinSigma(sigma, guessVector, *precontractedK_);
}

void CISigmaVectorEvaluator::formSigma2(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector) const {
  contractSameSpinSigma(sigma, guessVector, precontractedG_);
  formSigmaOppositeSpin(sigma, guessVector);
}

void CISigmaVectorEvaluator::formSigmaOppositeSpin(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector) const {
  int nOrbs = onvs_[0].size();

  auto dmat = constructD(guessVector);

#pragma omp parallel for default(none) shared(sigma, guessVector, dmat, nOrbs)
  for (int detIdx = 0; detIdx < int(basis_.size()); ++detIdx) {
    int ia = detToOnvs_[detIdx].first;
    int ib = detToOnvs_[detIdx].second;
    const Eigen::SparseMatrix<double>& dIaJbpq = dmat[ia];
    for (int k = 0; k < dIaJbpq.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(dIaJbpq, k); it; ++it) {
        int p = it.row() / nOrbs;
        int q = it.row() % nOrbs;
        int jb = k;
        double element = 0;
        int rsIdx = kroeneckerDeltaMatrix_(ib, jb);
        if (rsIdx != None) {
          element += (rsIdx > 0 ? 1.0 : -1.0) * (*g_)(p * nOrbs + q, std::abs(rsIdx));
        }
        else if (jb == ib) {
          for (int rsIdx : diagonalKDElements_[ib]) {
            assert(rsIdx >= 0);
            element += (*g_)(p * nOrbs + q, std::abs(rsIdx));
          }
        }
        sigma(detIdx) += it.value() * element;
      }
    }
  }
}

auto CISigmaVectorEvaluator::constructD(const Eigen::VectorXd& guessVector) const -> std::vector<Eigen::SparseMatrix<double>> {
  int nOrbs = onvs_[0].size();
  int nOcc = onvs_[0].countOccupied();

  std::vector<Eigen::SparseMatrix<double>> dMatrix(onvs_.size());

#pragma omp parallel for
  for (int onvIdx = 0; onvIdx < int(onvs_.size()); ++onvIdx) {
    const auto& onvI = onvs_[onvIdx];
    Eigen::SparseMatrix<double> dJbpq(nOrbs * nOrbs, onvs_.size());
    using T = Eigen::Triplet<double>;
    std::vector<T> tripletList;
    tripletList.reserve(onvs_.size() * nOcc * (nOrbs - nOcc));

    for (int p = 0; p < nOrbs; ++p) {
      if (onvI.isOccupied(p)) {
        for (int q = 0; q < nOrbs; ++q) {
          if (!onvI.isOccupied(q) || q == p) {
            const int sign = excitationSign(onvI, p, q);
            auto onvJa = onvI;
            onvJa.flip(p);
            onvJa.flip(q);
            for (int jb = 0; jb < int(onvs_.size()); ++jb) {
              double element = 0;
              // Construct D and G vectors for known p and q
              const auto& connectedToJb = onvConnection_[jb];
              const auto& betaIdxJb = betaIndices_[jb];
              for (int ja = 0; ja < int(connectedToJb.size()); ++ja) {
                if (onvJa == onvs_[connectedToJb[ja]]) {
                  element += sign * guessVector.data()[betaIdxJb[ja]];
                }
              }
              if (std::abs(element) > 1e-16) {
                tripletList.emplace_back(p * nOrbs + q, jb, element);
              }
            }
          }
        }
      }
    }
    dJbpq.setFromTriplets(tripletList.begin(), tripletList.end());
    dMatrix[onvIdx] = dJbpq;
  }
  return dMatrix;
}

void CISigmaVectorEvaluator::contractSameSpinSigma(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector,
                                                   const Eigen::MatrixXd& matrixToContract) const {
  auto contract = [&](int determinantIdx, int iSpin1, int iSpin2, const auto& connectedToSpin) {
    for (int j = 0; j < int(onvConnection_[iSpin1].size()); ++j) {
      sigma(determinantIdx) += matrixToContract(iSpin2, onvConnection_[iSpin1][j]) * guessVector(connectedToSpin[iSpin1][j]);
    }
  };

  for (int detIdx = 0; detIdx < int(basis_.size()); ++detIdx) {
    int ia = detToOnvs_[detIdx].first;
    int ib = detToOnvs_[detIdx].second;

    contract(detIdx, ia, ib, alphaIndices_);
    contract(detIdx, ib, ia, betaIndices_);
  }
}

auto CISigmaVectorEvaluator::precontractG(const Eigen::MatrixXd& twoElectronIntegrals) const -> Eigen::MatrixXd {
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(onvs_.size(), onvs_.size());

  const int nOrbs = onvs_[0].size();

  for (long i = 0; i < long(onvs_.size()); ++i) {
    const auto& onvI = onvs_[i];
    for (int occ = 0; occ < nOrbs; ++occ) {
      if (onvI.isOccupied(occ)) {
        for (int vir = 0; vir < nOrbs; ++vir) {
          if (!onvI.isOccupied(vir) || vir == occ) {
            auto onv2 = onvI;
            onv2.flip(occ);
            onv2.flip(vir);
            for (int occ2 = 0; occ2 < nOrbs; ++occ2) {
              if (onv2.isOccupied(occ2)) {
                for (int vir2 = 0; vir2 < nOrbs; ++vir2) {
                  if (!onv2.isOccupied(vir2) || vir2 == occ2) {
                    auto onvJ = onv2;
                    onvJ.flip(occ2);
                    onvJ.flip(vir2);

                    auto it = indexMap_.find(onvJ);
                    if (it != indexMap_.end()) {
                      G(i, it->second) += twoElectronIntegrals(occ * nOrbs + vir, occ2 * nOrbs + vir2) *
                                          excitationSign(onvI, occ, vir) * excitationSign(onv2, occ2, vir2);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return 0.5 * G;
}

auto CISigmaVectorEvaluator::constructK(const boost::optional<Eigen::MatrixXd>& oneBodyIntegral,
                                        const boost::optional<Eigen::MatrixXd>& twoBodyIntegral)
    -> boost::optional<Eigen::MatrixXd> {
  if (!oneBodyIntegral && !twoBodyIntegral) {
    throw std::runtime_error("No integrals provided to Sigma Vector evaluator for DirectCI.");
  }
  if (!twoBodyIntegral) {
    return oneBodyIntegral;
  }
  long nOrbs = oneBodyIntegral ? oneBodyIntegral->rows() : long(std::sqrt(twoBodyIntegral->rows()));
  Eigen::MatrixXd k = oneBodyIntegral ? *oneBodyIntegral : Eigen::MatrixXd::Zero(nOrbs, nOrbs);
  const Eigen::MatrixXd& g = *twoBodyIntegral;
  for (int p = 0; p < nOrbs; ++p) {
    for (int q = 0; q < nOrbs; ++q) {
      double twoBodyContributions = 0;
      for (int r = 0; r < nOrbs; ++r) {
        twoBodyContributions += g(p * nOrbs + r, r * nOrbs + q);
      }
      k(p, q) -= 0.5 * twoBodyContributions;
    }
  }

  return k;
}

auto CISigmaVectorEvaluator::precontractK(const boost::optional<Eigen::MatrixXd>& k) const -> boost::optional<Eigen::MatrixXd> {
  if (!k) {
    return {};
  }
  Eigen::MatrixXd precontractedK = Eigen::MatrixXd::Zero(onvs_.size(), onvs_.size());
  for (int i = 0; i < int(onvs_.size()); ++i) {
    for (int j = 0; j < int(onvs_.size()); ++j) {
      // off diagonal element
      int rsIdx = kroeneckerDeltaMatrix_(i, j);
      if (rsIdx != None) {
        precontractedK(i, j) += (rsIdx > 0 ? 1.0 : -1.0) * (k->data()[std::abs(rsIdx)]);
      }
      else if (i == j) {
        // Add all pq excitation where p and q are equal and an occupied orbital
        // in case of <K_alpha | E_rs | K_alpha>, where all occupied orbitals
        // connect the two ONV
        for (int rsIdx : diagonalKDElements_[i]) {
          precontractedK(i, i) += k->data()[rsIdx];
        }
      }
    }
  }

  return precontractedK;
}

auto CISigmaVectorEvaluator::indices(SpinComponent spin) const -> const std::vector<std::vector<int>>& {
  return spin == SpinComponent::Alpha ? alphaIndices_ : betaIndices_;
}
auto CISigmaVectorEvaluator::connection() const -> const std::vector<std::vector<int>>& {
  return onvConnection_;
}

auto CISigmaVectorEvaluator::onvToIndexMap(const OccupationNumberVector& onv) const -> int {
  return indexMap_.at(onv);
}

int CISigmaVectorEvaluator::excitationSign(const OccupationNumberVector& onv, int toDestroy, int toCreate) {
  const int min = std::min(toDestroy, toCreate);
  const int max = std::max(toDestroy, toCreate);
  const long nInterchanges =
      (toDestroy == toCreate) ? 0 : std::count(onv.getOccupation().begin() + min + 1, onv.getOccupation().begin() + max, true);
  return (nInterchanges % 2 == 0) ? 1 : -1;
}
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
