/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_OCCUPATIONNUMBERVECTOR_H
#define UTILS_OCCUPATIONNUMBERVECTOR_H

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

/**
 * @class OccupationNumberVector @file OccupationNumberVector.h
 * @brief Class representing an occupation number vector.
 * This header-only class has methods to manipulate efficiently the underlying
 * std::vector<bool>. boost::dynamic_bitset was also tried out but was found out
 * to be less efficient for some operations.
 */
class OccupationNumberVector {
 public:
  // Types definition
  using InputType = std::vector<int>;
  using OutputType = bool;
  /**
   * @brief Construct a new ONVInterface object.
   * @param input vector of integer that must be either 0 or 1.
   */
  explicit OccupationNumberVector(const InputType& onv) {
    occupationVector_.resize(onv.size(), false);
    int jcont = 0;
    for (int i_alpha : onv) {
      if (i_alpha == 1) {
        occupationVector_[jcont] = true;
      }
      jcont++;
    }
  }
  /**
   * @brief Move-construct a new ONVInterface object.
   * @param bool vector.
   */
  explicit OccupationNumberVector(std::vector<bool>&& onv) : occupationVector_(std::move(onv)) {
  }

  /**
   * @brief Default constructor, needed by the instantiation of std::map<OrbitalType, OccupationNumberVector>
   */
  OccupationNumberVector() = default;

  /**
   * @brief Default copy and move constructors and assignment operators.
   */
  OccupationNumberVector(const OccupationNumberVector& rhs) = default;
  OccupationNumberVector(OccupationNumberVector&& rhs) = default;
  OccupationNumberVector& operator=(const OccupationNumberVector& rhs) = default;
  OccupationNumberVector& operator=(OccupationNumberVector&& rhs) noexcept = default;
  /** @} */
  /**
   * @brief Counts how many occupied orbitals are present in the determinant.
   * @return int Number of occupied orbitals.
   */
  int countOccupied() const {
    return std::count(occupationVector_.begin(), occupationVector_.end(), 1);
  }
  /**
   * @brief Counts how many virtual orbitals are present in the determinant.
   * @return int Number of virtual orbitals.
   */
  int countVirtual() const {
    return std::count(occupationVector_.begin(), occupationVector_.end(), 0);
  }
  /**
   * @brief Checks if the orbital on position i is occupied.
   * @param i orbital to be checked.
   * @return true if the orbital is occupied.
   * @return false if the orbital is not occupied.
   */
  bool isOccupied(int i) const {
    if (i >= occupationVector_.size()) {
      throw std::runtime_error("Orbital not present in Determinant.");
    }
    return occupationVector_[i];
  }
  /**
   * @brief Changes the occupation of a orbital.
   * @param i position of the orbital to be flipped.
   */
  void flip(int i) {
    occupationVector_[i] = !occupationVector_[i];
  }
  /**
   * @brief Size of the ONV
   * @return int dimension of the data structure storing the occupation number vector.
   */
  int size() const {
    return occupationVector_.size();
  }
  /**
   * @brief Access operator
   * @param idx index of the element to be accessed.
   * @return int& value of the occupation
   */
  bool getElement(int index) const {
    return occupationVector_[index];
  }
  /**
   * @brief Access operator
   * @return the whole occupation number vector
   */
  const auto& getOccupation() const {
    return occupationVector_;
  }

  /**
   * @brief Comparison operator
   */
  bool operator==(const OccupationNumberVector& other) const {
    return (getOccupation() == other.getOccupation());
  }

 private:
  std::vector<bool> occupationVector_;
};

struct OnvHash {
  std::size_t operator()(const OccupationNumberVector& key) const {
    return hash_(key.getOccupation());
  }

 private:
  std::hash<std::vector<bool>> hash_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_OCCUPATIONNUMBERVECTOR_H
