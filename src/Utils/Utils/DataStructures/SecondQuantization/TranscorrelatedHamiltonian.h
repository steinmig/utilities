/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SQ_TRANSCORRELATEDHAMILTONIAN_H
#define UTILS_SQ_TRANSCORRELATEDHAMILTONIAN_H

// #include "Hamiltonian.h"
#include "HeatBathCIContainer.h"
#include "NonHermitianHamiltonian.h"
#include <array>
#include <boost/functional/hash.hpp>
#include <unordered_map>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

class TranscorrelatedHamiltonian : public NonHermitianHamiltonian<TranscorrelatedHamiltonian> {
 public:
  // Types definition
  using Base = NonHermitianHamiltonian<TranscorrelatedHamiltonian>;
  // using hbciContainerSameSpin_ = typename Base::hbciContainerSameSpin_;
  using ThreeBodyIdentifier = std::array<int, 6>;
  using ThreeBodyContainer = std::unordered_map<ThreeBodyIdentifier, double, boost::hash<ThreeBodyIdentifier>>;

  /**
   * @brief Calls the Constructor of the baseclass.
   * @param moBasisSize the size of the mo basis.
   */
  explicit TranscorrelatedHamiltonian(int moBasisSize) : Base(moBasisSize){};

  /** @brief Implementation to generate all connect for nonhermitian Hamiltonians with a three body operator */
  std::vector<ElectronicDeterminant> generateAllConnectedImpl(const ElectronicDeterminant& det, bool addSmae = true) const;
  /** @brief Implementation of the method to generate triple excitations */
  // std::vector<ElectronicDeterminant> generateTripleExcitations(const ElectronicDeterminant& det, bool addSame = true)
  // const;

  // TODO: write this function
  /**
   * @brief Brings the key for the three-body term to normal form
   *
  //  * Here by normal form we mean that, for a three-electron orbital < ij | kl > (in physics notation), we must
  //  * have i<j, k<l, i<k and j<l. In this way, we can exploit the eight-fold symmetry of the tensor
   *
   * @param input_tb_identifier Tuple of integers labelling a three-body term.
   * @return the TwoBodyIdentifier in the correct form
   */
  ThreeBodyIdentifier bringToNormalFormImpl(const ThreeBodyIdentifier& input_ob_identifier) const;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#include "TranscorrelatedHamiltonian.hpp"

#endif