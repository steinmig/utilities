/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SQ_HAMILTONIAN_H
#define UTILS_SQ_HAMILTONIAN_H

#include "HeatBathCIContainer.h"
#include "SQSpecifier.h"
#include <Eigen/Dense>
#include <array>
#include <boost/functional/hash.hpp>
#include <unordered_map>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

class ElectronicDeterminant;
class Excitation;
struct SingleExcitation;

template<class Crtp>
class Hamiltonian {
 public:
  // Types definition
  using OneBodyIdentifier = std::pair<int, int>;
  using TwoBodyIdentifier = std::array<int, 4>;
  using OneBodyContainer = Eigen::MatrixXd;
  using TwoBodyContainer = std::unordered_map<TwoBodyIdentifier, double, boost::hash<TwoBodyIdentifier>>;

  /**
   * @brief Constructor initializing the OneBodyContainer to 0.
   * TwoBodyContainer does not need initializing as it is a sparse Tensor and no Eigen::MatrixXd.
   * @param moBasisSize the size of the mo basis.
   */
  explicit Hamiltonian(int moBasisSize);

  /**
   * @brief Constructor initializing the OneBodyContainer to 0.
   * TwoBodyContainer does not need initializing as it is a sparse Tensor and no Eigen::MatrixXd.
   * @param moBasisSize the size of the mo basis.
   */
  explicit Hamiltonian(const Specifier& specifier);

  /**
   * @brief Prints the coefficients of the Hamiltonian
   */
  void print(std::ostream& out) const;
  /** @brief Add a one body term to the SQ Hamiltonian.
   * @param i Index of the first operator.
   * @param j Index of the second operator.
   * @param coeff Coefficient in the SQ representation
   */
  void addTerm(int i, int j, double coeff);
  /**
   * @brief Adds all one-body terms stored in an Eigen matrix.
   *
   * The matrix is expected to be lower or upper diagonal.
   *
   * @param InputOneBody Eigen::Matrix storing the one-body terms.
   */
  void addOneBodyContribution(const Eigen::MatrixXd& inputOneBody);
  /**
   * @brief Adds all two-body terms stored in an in-house 4-dimensional tensor.
   * @param InputTwoBody in-house Tensor with 4 indexes that stores the two-body part of the Hamiltonian.
   */
  void addTwoBodyContribution(const TwoBodyContainer& inputTwoBody, bool isPhys = false);
  /**
   * @brief Adds all two-body terms stored in an Eigen::MatrixXd in chemstry format.
   * @param InputTwoBody in-house Tensor with 4 indexes that stores the two-body part of the Hamiltonian.
   */
  void addTwoBodyContribution(const Eigen::MatrixXd& inputTwoBody);
  /**
   * @brief Setter for the core energy
   * @param coreEnergy value of the core energy
   */
  void addConstantTerm(double coreEnergy);
  /** @brief Populates the HBCI container */
  void populateHBCIContainer();
  /**
   * @brief Get the Number of one body terms.
   * @return int number of one body terms in the Hamiltonian.
   */
  int getNumberOneBody() const;
  /**
   * @brief Get the Number of two body terms.
   * @return int number of two body terms in the Hamiltonian.
   */
  int getNumberTwoBody() const;

  /**
   * @brief Calculates the matrix element of H between det1 and det2 for FH.
   *
   * The Slater-Condon rules are employed in the calculation of the matrix elements (i.e., determinants
   * that differ by more than 2 excitations are assumed to be zero). Also for FH There are some spin restrictions
   * reducing the number of integrals This routine calculates which orbitals are involved in the excitation, as well as
   * the orbitals that are occupied. If these informations are already available, use the other overloaded version of
   * CalculateMatrixElement.
   *
   * @tparam Determinant class representing the determinant. Must have two methods, specifically GetExcitations,
   *         that returns the vector of all excitations, and GetAllOccupied, that returns the indexes of the
   *         occupied orbitals for a given spin.
   * @param det1 determinant appearing in the bra.
   * @param det2 determinant appearing in the ket.
   * @return double < det1 | H | det2 >.
   */
  double calculateMatrixElement(const ElectronicDeterminant& det1, const ElectronicDeterminant& det2) const;

  /**
   * @brief Calculates the sign of a single excitation based on the Jordan-Wigner mapping.
   *
   * Note that here we employ the JW mapping to calculate the parity of the excitation.
   * Remember that the Jordan-Wigner transformation transforms a Hopping term proportional to a_i^+ a_j as (assuming
   * i>j)
   *
   * b_i^+ x f_{i-1} x f_{i+2} x ... x f_{j} x b_j (for alpha electrons)
   * b_i^+ x f_{i} x f_{i+2} x ... x f_{j+1} x b_j (for beta electrons)
   *
   * If i<j we the same equations, but they have before the creator and then the annihilator and, in addition,
   * we have a - sign coming from the anticommutation rule.
   * Note also that the sign is unique and is the same for the one- and the two-body term.
   *
   * @param lst_occupied_alpha list with the occupied alpha orbitals.
   * @param lst_occupied_beta list with the occupied beta orbitals.
   * @param exc1 identifier of the excitation.
   * @return int sign of the excitation.
   */
  int calculateExcitationSign(const std::vector<int>& lst_occupied_alpha, const std::vector<int>& lst_occupied_beta,
                              const SingleExcitation& exc1) const;

  /**
   * @brief Overload of the previous function
   *
   * @param lst_occupied_alpha
   * @param lst_occupied_beta
   * @param exc1
   * @param exc2
   * @return int
   */
  int calculateExcitationSign(const std::vector<int>& lst_occupied_alpha, const std::vector<int>& lst_occupied_beta,
                              const SingleExcitation& exc1, const SingleExcitation& exc2) const;

  /** @brief Implementation of the method to calculate < det1 | H | det2 > */
  double calculateMatrixElement(const Excitation& excitations, const std::vector<int>& i_occ_alpha,
                                const std::vector<int>& i_occ_beta) const;

  /** @brief Method to generate all excitations */
  std::vector<ElectronicDeterminant> generateAllConnected(const ElectronicDeterminant& det, bool addSame = true) const;
  /** @brief Implementation of the method to generate all excitations */
  std::vector<ElectronicDeterminant> generateAllConnectedImpl(const ElectronicDeterminant& det, bool addSame = true) const;
  /** @brief Method to generate all excitations but doubles from a subset */
  std::vector<ElectronicDeterminant>
  generateAllConnected(const ElectronicDeterminant& det, std::map<SpinComponent, std::vector<int>> listDoublesOccupied,
                       std::map<SpinComponent, std::vector<int>> listDoublesVirtual, bool addSame = true) const;
  /** @brief Implementation of the method to generate all excitations but doubles from a subset */
  std::vector<ElectronicDeterminant>
  generateAllConnectedImpl(const ElectronicDeterminant& det, std::map<SpinComponent, std::vector<int>> listDoublesOccupied,
                           std::map<SpinComponent, std::vector<int>> listDoublesVirtual, bool addSame = true) const;
  /** @brief Implementation of the method to generate single excitations */
  std::vector<ElectronicDeterminant> generateSingleExcitations(const ElectronicDeterminant& det, bool addSame = true) const;
  /** @brief Implementation of the method to generate double excitations */
  std::vector<ElectronicDeterminant> generateDoubleExcitations(const ElectronicDeterminant& det, bool addSame = true) const;
  /** @brief Implementation of the method to generate double excitations from an orbital subset */
  std::vector<ElectronicDeterminant>
  generateDoubleExcitations(const ElectronicDeterminant& det, std::map<SpinComponent, std::vector<int>> listOccupied,
                            std::map<SpinComponent, std::vector<int>> listVirtual, bool addSame = true) const;

  /** @brief Implementation of the method to count all possible excitations */
  int getNumberOfConnectedDeterminants(const ElectronicDeterminant& det, bool addSame = true) const;

  /**
   * @brief Changes the threshold epsilon for the HBCI Hamiltonian elements valuation.
   * Loops over all elements of the hbci map and updates the threshold.
   */
  void setHbciThreshold(double threshold);

 private:
  double matrixElementSameDeterminant(const std::vector<int>& lst_occupied_alpha, const std::vector<int>& lst_occupied_beta) const;

  /**
   * @brief Calculates the matrix element of the Hamiltonian for a single excitation.
   * @param exc1 ExcitationDescriptor for the target excitation.
   * @param lst_occupied_alpha list of occupied orbitals for alpha electrons.
   * @param lst_occupied_alpha list of occupied orbitals for beta electrons.
   * @return double < Psi_{ia} | H | Psi >
   */
  double matrixElementSingleExcitation(const SingleExcitation& exc1, const std::vector<int>& lst_occupied_alpha,
                                       const std::vector<int>& lst_occupied_beta) const;

  /**
   * @brief Calculates the matrix element of the Hamiltonian for a single excitation.
   * @param exc1 ExcitationDescriptor for the first excitation.
   * @param exc2 ExcitationDescriptor for the second excitation.
   * @param lst_occupied_alpha list of occupied orbitals for alpha electrons (required for the Jordan-Wigner mapping).
   * @param lst_occupied_alpha list of occupied orbitals for beta electrons (required for the Jordan-Wigner mapping).
   * @return double < Psi_{ia} | H | Psi >
   */
  double matrixElementDoubleExcitation(const SingleExcitation& exc1, const SingleExcitation& exc2,
                                       const std::vector<int>& lst_occupied_alpha,
                                       const std::vector<int>& lst_occupied_beta) const;
  /**
   * @brief Add the iOcc1 --> iVirt1 / iOcc2 --> iVirt2 excitation to the HBCI container
   * @param iOcc1 index of the first occupied orbital
   * @param iOcc2 index of the second occupied orbital
   * @param iVirt1 index of the first virtual orbital
   * @param iVirt2 index of the second virtual orbital
   * @param value coefficient of the excitation
   */
  void addTermToHBCIContainer(int iOcc1, int iOcc2, int iVir1, int iVir2, double value);

  /**
   * @brief Brings the key for the one-body term to normal form
   *
   * Here by normal form we mean that i<j. In this way, we can exploit the two-fold symmetry of the
   * one-body integral and save only one of the terms.
   *
   * @param input_ob_identifier Pair of integers labelling a one-body term.
   * @return the OneBodyIdentifier in the correct form
   */
  OneBodyIdentifier bringToNormalForm(const OneBodyIdentifier& input_ob_identifier) const;

  /**
   * @brief Brings the key for the two-body term to normal form
   *
   * Here by normal form we mean that, for a two-electron orbital < ij | kl > (in physics notation), we must
   * have i<j, k<l, i<k and j<l. In this way, we can exploit the eight-fold symmetry of the tensor
   *
   * @param input_tb_identifier Tuple of integers labelling a two-body term.
   * @return the TwoBodyIdentifier in the correct form
   */
  TwoBodyIdentifier bringToNormalForm(const TwoBodyIdentifier& input_ob_identifier) const;

  /**
   * @brief Returns the CRTP.
   *
   * Here the CRTP is returns as a static cast as const and non const function.
   * @return the CRTP
   */
  const auto& derived() const;
  auto& derived();

  // Members
  /** @brief Number of molecular orbitals in the system */
  int nMOs_;
  /** @brief Core contribution to the Hamiltonian */
  double coreEnergy_ = 0.;
  /** @brief unordered map that associates a pair, labelling a one-body term, to a scalar. */
  OneBodyContainer oneElectronMatrix_;
  /** @brief unordered map that associates a tuple of integers, labelling a two-body term, to a scalar. */
  TwoBodyContainer twoElectronMatrix_;
  static constexpr double thresholdForTerms_ = 1.0E-10;
  double hbciThreshold_ = 0.0;

 protected:
  /** @brief Data structures for HBCI calculation. */
  std::map<OrderedIndex, HeatBathCIContainer> hbciContainerSameSpin_, hbciContainerOppositeSpin_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#include "Hamiltonian.hpp"
#endif // SPARROW_SQ_HAMILTONIAN_H
