/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CISOLVER_H
#define UTILS_CISOLVER_H

#include "HermitianHamiltonian.h"
#include <iostream>

namespace Scine {
namespace Core {
struct Log;
}
namespace Utils {
class Settings;
class SigmaVectorEvaluator;
namespace SecondQuantization {
struct MOIndicesForReferences {
  std::vector<int> data;
};
struct MOIndicesForDoubles {
  std::vector<int> data;
};

class CISigmaVectorEvaluator;

class CISolver {
 public:
  /// @brief Constructor from an Hamiltonian.
  CISolver(const Settings& settings, const Specifier& specifier);
  /// @brief destructor.
  ~CISolver();

  /**
   * @brief Sets the indices defining the MO subspace from which reference determinants are constructed.
   */
  void setIndices(MOIndicesForReferences moIndicesReferences);
  /**
   * @brief Sets the indices defining the MO subspace from which double excitations are constructed.
   */
  void setIndices(MOIndicesForDoubles moIndicesDoubles);
  /**
   * @brief Sets whether the reference determinants be created from symmetric alpha/beta double excitations from the
   * reference.
   *
   * For example, with a reference determinant
   *
   * 222000
   *
   * and indices = {2,3}
   *
   * the reference determiants will normally be
   *
   * 222000, 220200, 22ab00, 22ba00.
   *
   * With this option set to true, the reference determinants are chosen as:
   *
   * 222000, 220200.
   *
   */
  void setReferenceFromDoubleExcitations(bool referenceFromDoubles);

  auto solve(Core::Log& log) -> EigenContainer;

  auto generateCIBasis(const std::vector<ElectronicDeterminant>& references) const -> std::vector<ElectronicDeterminant>;

  auto basis() const -> const std::vector<ElectronicDeterminant>&;
  auto basis() -> std::vector<ElectronicDeterminant>&;

  auto diagonalize(Core::Log& log, std::shared_ptr<SigmaVectorEvaluator> sve, const Eigen::VectorXd& diagonal) -> EigenContainer;

  /**
   * @brief If there are indices set in moIndicesReferences, create references accordingly,
   *        otherwise only the Hartree--Fock Determinant is returned.
   */
  auto generateReferenceDeterminants() const -> std::vector<ElectronicDeterminant>;

  /**
   * @brief Converts the sign convention of Direct CI to the one of Hamiltonian.
   *        Direct CI represents determinant as |onv_alpha>x|onv_beta>, whereas
   *        the determinant representation in Hamiltonian is
   *        |\phi_{1\alpha} \phi_{1\beta} \phi_{2\alpha}...>.
   *        The conversion of the two leads to a sign mismatch. Properties calculated
   *        on them assume the same sign convention. This function converts
   *        the representation |onv_alpha>x|onv_beta> to
   *        |\phi_{1\alpha} \phi_{1\beta} \phi_{2\alpha}...>.
   */
  void convertSigns(EigenContainer& result) const;

  /**
   * @brief returns the excitation sign for a certain determinant.
   */
  auto getExcitationSign(const ElectronicDeterminant& det, const SingleExcitation& excitation) const -> double;

 private:
  static void checkAndCorrectNumberOfRoots(int& numberOfEnergyLevels, int& initialSubspaceDimension, int nConfigurations);
  auto getDiagonal() const -> std::unordered_map<ElectronicDeterminant, double, ElectronicDeterminantHash>;
  /**
   * @brief Gets the occuparion lists to be given to Hamiltonian::generateAllConnected.
   *
   * If it is empty, only generateSingles sould be called.
   */
  auto getOccupationLists(const ElectronicDeterminant& reference) const
      -> boost::optional<std::array<std::map<SpinComponent, std::vector<int>>, 2>>;
  const Settings& settings_;
  const MoIntegrals& integrals_;
  HermitianHamiltonian hamiltonian_;
  boost::optional<MOIndicesForReferences> moIndicesReferences_;
  boost::optional<MOIndicesForDoubles> moIndicesDoubles_;
  ElectronicDeterminant hfDet_;
  bool onlyDoubleExcitations_{false};
  std::vector<ElectronicDeterminant> basis_;
  double coreEnergy_{};
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_CISOLVER_H
