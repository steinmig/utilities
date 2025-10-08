/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXCITATIONHELPER_H
#define UTILS_EXCITATIONHELPER_H

#include <cassert>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

/**
 * @brief Enum to distinguish alpha and beta electrons.
 */
enum class SpinComponent { Alpha, Beta };

enum class DoubleExcitationDescriptor { AlphaAlpha, AlphaBeta, BetaBeta };

struct SingleExcitation {
  SingleExcitation(int occ, int vir, SpinComponent orbSpin)
    : occupiedOrbital(occ), virtualOrbital(vir), orbitalSpin(orbSpin) {
  }
  int occupiedOrbital, virtualOrbital;
  SpinComponent orbitalSpin;
};

struct CrossExcitation {
  CrossExcitation(int occ, int vir, std::pair<SpinComponent, SpinComponent> orbSpins)
    : occupiedOrbital(occ), virtualOrbital(vir), orbitalSpins(std::move(orbSpins)) {
  }
  int occupiedOrbital, virtualOrbital;
  std::pair<SpinComponent, SpinComponent> orbitalSpins;
};

/**
 * @brief Simple structure to characterize an excitation.
 */
class Excitation {
 public:
  using const_iterator = std::vector<SingleExcitation>::const_iterator;
  /**
   * @brief Construct an empty new Excitation object.
   */
  Excitation() = default;

  /**
   * @brief Return the number of single excitations in the overall excitation.
   * @return an integer with the number of single excitations.
   */
  int cardinality() const {
    return excitationVector.size();
  }

  /**
   * @brief Method to add a single excitation to the generic excitation descriptor.
   * @param occupiedIndex Index of the occupied orbital from where the electron was excited.
   * @param virtualIndex Index of the virtual orbital to which the electron was excited.
   * @param orbitalSpin The spin of the two latter orbitals.
   */
  void createExcitation(int occupiedIndex, int virtualIndex, SpinComponent orbitalSpin) {
    excitationVector.emplace_back(occupiedIndex, virtualIndex, orbitalSpin);
  }

  const SingleExcitation& getSingleExcitation(int index) const {
    return excitationVector[index];
  }

  const_iterator begin() const {
    return excitationVector.begin();
  }

  const_iterator end() const {
    return excitationVector.end();
  }

 private:
  std::vector<SingleExcitation> excitationVector;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXCITATIONHELPER_H
