/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_CISIGMAVECTOREVALUATOR_H
#define UTILSOS_CISIGMAVECTOREVALUATOR_H

#include "ElectronicDeterminant.h"
#include "ExcitationHelpers.h"
#include "OccupationNumberVector.h"
#include <Utils/Math/IterativeDiagonalizer/SigmaVectorEvaluator.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <boost/optional.hpp>
#include <functional>
#include <map>
#include <unordered_map>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

class ElectronicDeterminant;
struct MoIntegrals;

/**
 * @class CISigmaVectorEvaluator @file SigmaVectorEvaluator.h
 * @brief Sigma vector evaluation for CI problems.
 * Implemented as indicated in chapter 11.8.3
 * Helgaker, Jorgensen, Olsen, Molecular Electronic Structure Theory
 *
 *
 * RIGHT NOW ONLY FOR CLOSED SHELLS!!!
 * S_z != 0, nAlpha == nBeta.
 */
class CISigmaVectorEvaluator final : public SigmaVectorEvaluator {
 public:
  static constexpr const int None = 0;

  CISigmaVectorEvaluator(const std::vector<ElectronicDeterminant>& basis, const MoIntegrals& integrals, double coreEnergy = 0.0);
  ~CISigmaVectorEvaluator() final = default;
  const Eigen::MatrixXd& evaluate(const Eigen::MatrixXd& guessVectors) const final;

  /**
   * @brief Allows for internal handling of new subspace dimension.
   * @param newSubspaceDimension The new guess vectors number after subspace collapse.
   */
  void collapsed(int /*newSubspaceDimension*/) final;

  /**
   *  @brief Functions to create the sigma vectors.
   *  @{
   */
  void formSigma1(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector) const;
  void formSigma2(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector) const;
  void contractSameSpinSigma(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector,
                             const Eigen::MatrixXd& matrixToContract) const;
  void formSigmaOppositeSpin(Eigen::VectorXd& sigma, const Eigen::VectorXd& guessVector) const;
  /** @} **/
  /**
   * @brief Creates an unique sorted set of onvs. Assumes that alpha and beta occuaptions are equal.
   */
  static auto getSortedOnvList(const std::vector<ElectronicDeterminant>& basis) -> std::vector<OccupationNumberVector>;

  /**
   * @brief Constructs the matrix
   * <K_\beta|E^\beta_{rs}| J_\beta>
   * This needs to be done just once in the beginning.
   * If the numer of alpha and beta electrons is the same, only one needs to be calculated.
   * NOTE ON SIGNS:
   * - < Ka | E_pq^a | Ja> is NOT just 0 or 1, it can also be -1.
   *   The sign is given by the Jordan-Wigner transformation.
   */
  auto constructKroeneckerDeltaMatrix(const std::vector<OccupationNumberVector>& onvs) const
      -> std::tuple<Eigen::MatrixXi, std::vector<std::vector<int>>>;

  auto precontractG(const Eigen::MatrixXd& twoElectronIntegrals) const -> Eigen::MatrixXd;
  static auto constructK(const boost::optional<Eigen::MatrixXd>& oneBodyIntegral,
                         const boost::optional<Eigen::MatrixXd>& twoBodyIntegral) -> boost::optional<Eigen::MatrixXd>;
  /**
   *  If the spatial parts of alpha and beta are the same this is just one matrix.
   *  Otherwise two matrices for alpha and beta (multiplied with the corresponding k)
   *  should be created.
   *  Contains the elements:
   *  K_{onv_i, onv_j} = sum_{pq} k_{pq} * <onv_i | E^{\alpha/\beta}_{pq} | onv_j>
   *
   *  The sigma vector is then calculated as
   *
   *  sigma^(1) = sigma^\alpha + sigma^\beta =
   *  \sum_{onv^\alpha_j} k_{onv^alpha_i, onv^alpha_j} * C_{onv^\alpha_j, onv^\beta_i}
   *  +
   *  \sum_{onv^\beta_j} k_{onv^\beta_i, onv^\beta_j} * C_{onv^\alpha_i, onv^\beta_j}
   */
  auto precontractK(const boost::optional<Eigen::MatrixXd>& k) const -> boost::optional<Eigen::MatrixXd>;

  auto indices(SpinComponent spin) const -> const std::vector<std::vector<int>>&;
  auto connection() const -> const std::vector<std::vector<int>>&;
  auto onvToIndexMap(const OccupationNumberVector& onv) const -> int;

  /*
   * @brief Implementation of the Jordan-Wigner transformation for obtaining the sign
   * NOTE ON SIGNS:
   * - < Ka | E_pq^a | Ja> is NOT just 0 or 1, it can also be -1.
   *   The sign is given by the Jordan-Wigner transformation.
   *
   */
  static int excitationSign(const OccupationNumberVector& onv, int toDestroy, int toCreate);

 private:
  auto constructD(const Eigen::VectorXd& guessVector) const -> std::vector<Eigen::SparseMatrix<double>>;
  /**
   * @brief Constructs the maps of the indices between ONV and allowed determinants.
   *
   * connection() returns
   * for each ONV, which ONV of opposite spin builds a determinant in the basis.
   * the int returned is the index of the ONV of opposite spin in onvs_
   *
   * indices(SpinComponent spin) returns
   * for each ONV, which ONV (as addresses in connection()) of opposite spin builds a determinant in the basis.
   * the int populating the vectors are the index of the determinant in basis_.
   */
  void constructIndicesMaps();
  const std::vector<ElectronicDeterminant>& basis_;
  // one electron matrix element: k_pq = h_pq - 0.5 * \sum_r g_prrq
  double coreEnergy_{};
  boost::optional<Eigen::MatrixXd> precontractedK_;
  Eigen::MatrixXd precontractedG_;
  boost::optional<Eigen::MatrixXd> g_;
  mutable Eigen::MatrixXd sigmaVectors_;
  Eigen::MatrixXi kroeneckerDeltaMatrix_;
  std::unordered_map<OccupationNumberVector, int, OnvHash> indexMap_;
  std::vector<OccupationNumberVector> onvs_;
  // stores for each det in basis the index of the onvs.
  // first index ist alpha onv, second index is beta onv.
  std::vector<std::pair<int, int>> detToOnvs_;
  std::vector<std::vector<int>> alphaIndices_, betaIndices_;
  std::vector<std::vector<int>> onvConnection_;
  std::vector<std::vector<int>> diagonalKDElements_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_CISIGMAVECTOREVALUATOR_H
