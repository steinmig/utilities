/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SQ_NONHERMITIANHAMILTONIAN_H
#define UTILS_SQ_NONHERMITIANHAMILTONIAN_H

#include "Hamiltonian.h"

namespace Scine {
namespace Utils {
namespace SecondQuantization {

template<typename T = void>
class NonHermitianHamiltonian : public Hamiltonian<NonHermitianHamiltonian<T>> {
 public:
  // Types definition
  using Base = Hamiltonian<NonHermitianHamiltonian<T>>;
  using TwoBodyIdentifier = typename Base::TwoBodyIdentifier;

  /**
   * @brief Calls the Constructor of the baseclass.
   * @param moBasisSize the size of the mo basis.
   */
  explicit NonHermitianHamiltonian(int moBasisSize) : Base(moBasisSize){};

  /**
   * @brief Brings the key for the two-body term to normal form
   *
   * Here by normal form we mean that, for a two-electron orbital < ij | kl > (in physics notation), we must
   * have i<j, k<l, i<k and j<l. In this way, we can exploit the eight-fold symmetry of the tensor
   *
   * @param input_tb_identifier Tuple of integers labelling a two-body term.
   * @return the TwoBodyIdentifier in the correct form
   */
  TwoBodyIdentifier bringToNormalFormImpl(const TwoBodyIdentifier& input_ob_identifier) const;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#include "NonHermitianHamiltonian.hpp"

#endif