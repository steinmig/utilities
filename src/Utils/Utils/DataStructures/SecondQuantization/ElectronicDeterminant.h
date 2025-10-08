/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ELECTRONIC_DETERMINANT_H
#define UTILS_ELECTRONIC_DETERMINANT_H

#include "ExcitationHelpers.h"
#include "OccupationNumberVector.h"
#include <algorithm>
#include <cassert>
#include <complex>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

class InvalidSpinComponentException : public std::exception {
  const char* what() const noexcept final {
    return "Invalid orbital type. Only valid types: Alpha, Beta.";
  }
};

/**
 * @brief Exception raised when the ONV is contracted with a vector with an incorrect size.
 */
class InvalidSizeException : public std::exception {
  const char* what() const noexcept final {
    return "Invalid size of vector to be contracted with the ONV.";
  }
};

/**
 * @brief Exception raised when an occupied orbital is populated or a virtual orbital
 * is depopulated.
 */
class ExcitationException : public std::exception {
  const char* what() const noexcept final {
    return "Invalid excitation requested.";
  }
};

/**
 * @brief Class representing an electronic determinant.
 */
class ElectronicDeterminant {
 public:
  // Types definition
  using InputType = OccupationNumberVector::InputType;
  using OutputType = bool;

  /**
   * @brief Default empty constructor.
   */
  ElectronicDeterminant() = default;

  /**
   * @brief Construct a new Determinant object. Must provide the ONV for the alpha and beta electrons.
   * @param input_onv_alpha Input vector of integers (0,1) for the alpha electrons.
   * @param input_onv_beta Input vector of integers (0,1) for the beta electrons.
   */
  ElectronicDeterminant(const InputType& input_onv_alpha, const InputType& input_onv_beta) {
    alphaOnv_ = OccupationNumberVector(input_onv_alpha);
    betaOnv_ = OccupationNumberVector(input_onv_beta);
  };

  /**
   * @brief Construct a new Determinant object. Must provide the ONV for the alpha and beta electrons.
   * @param input_onv_alpha ONV object for the alpha electrons.
   * @param input_onv_beta ONV object for the beta electrons.
   */
  ElectronicDeterminant(OccupationNumberVector input_onv_alpha, OccupationNumberVector input_onv_beta) {
    alphaOnv_ = std::move(input_onv_alpha);
    betaOnv_ = std::move(input_onv_beta);
  };

  /**
   * @brief Construct a new Electronic Determinant object from an ONV string
   *
   * The string must be in the format "2aab0", where 2 means that the orbital is doubly occupied,
   * a that the orbital is occupied by an alpha electorn, b that it is occupied by a beta electron
   * and 0 that the orbital is empty
   *
   * @param onv_as_string input string with the ONV formatted as explained above
   */
  ElectronicDeterminant(const std::string& onv_as_string) {
    InputType input_onv_alpha(onv_as_string.size(), 0);
    InputType input_onv_beta(onv_as_string.size(), 0);
    for (int icont = 0; icont < int(onv_as_string.size()); icont++) {
      if (onv_as_string[icont] == '2') {
        input_onv_alpha[icont] = 1;
        input_onv_beta[icont] = 1;
      }
      else if (onv_as_string[icont] == 'a' || onv_as_string[icont] == 'u') {
        input_onv_alpha[icont] = 1;
      }
      else if (onv_as_string[icont] == 'b' || onv_as_string[icont] == 'd') {
        input_onv_beta[icont] = 1;
      }
    };
    alphaOnv_ = OccupationNumberVector(input_onv_alpha);
    betaOnv_ = OccupationNumberVector(input_onv_beta);
  }

  /**
   * @brief Copy assignment operator of a determinant.
   */
  ElectronicDeterminant& operator=(const ElectronicDeterminant& rhs) = default;
  /**
   * @brief Move assignment operator of a determinant.
   */
  ElectronicDeterminant& operator=(ElectronicDeterminant&& rhs) noexcept = default;
  /**
   * @brief Copy constructor of a determinant.
   */
  ElectronicDeterminant(const ElectronicDeterminant& rhs) = default;
  /**
   * @brief Move constructor of a determinant.
   */
  ElectronicDeterminant(ElectronicDeterminant&& rhs) noexcept = default;

  /**
   * @brief Getter for the ONV.
   * @param SpinType alpha/beta.
   * @return The Occupation Number Vector.
   */
  auto getOnv(SpinComponent spin) const -> const OccupationNumberVector& {
    return (spin == SpinComponent::Alpha) ? alphaOnv_ : betaOnv_;
  }
  /**
   * @brief Size getter.
   * @param SpinType alpha/beta.
   * @return int dimension of the ONV for the given spin.
   */
  int size(SpinComponent spin) const {
    return (spin == SpinComponent::Alpha) ? alphaOnv_.size() : betaOnv_.size();
  };

  /**
   * @brief Overall size getter.
   * @return int overall (alpha + beta) dimension of the ONV.
   */
  int size() const {
    return this->size(SpinComponent::Alpha) + this->size(SpinComponent::Beta);
  }

  /**
   * @brief Get the occupation of the i-th orbital of a given spin.
   * @param i index of the orbital.
   * @param SpinType spin of the orbital to be checked.
   * @return true if the orbital is occupied.
   * @return false if the orbital is empty.
   */
  OutputType isOccupied(SpinComponent SpinType, int i) const {
    return (SpinType == SpinComponent::Alpha) ? alphaOnv_.getElement(i) : betaOnv_.getElement(i);
  }

  /**
   * @brief Get the occupation of the i-th orbital of a given spin.
   * @param i index of the orbital.
   * @param SpinType spin of the orbital to be checked.
   * @return true if the orbital is occupied.
   * @return false if the orbital is empty.
   */
  const auto& getOccupation(SpinComponent SpinType) const {
    return (SpinType == SpinComponent::Alpha) ? alphaOnv_.getOccupation() : betaOnv_.getOccupation();
  }

  /**
   * @brief Get the occupied orbitals of a given spin.
   * @param SpinType spin of the orbital to be checked.
   * @return std::vector<int>: vector with the indexes of the occupied orbitals
   */
  std::vector<int> getAllOccupied(SpinComponent SpinType) const {
    const auto& onv = SpinType == SpinComponent::Alpha ? alphaOnv_ : betaOnv_;
    std::vector<int> res;
    res.reserve(onv.countOccupied());
    for (int i = 0; i < onv.size(); i++) {
      if (onv.isOccupied(i)) {
        res.push_back(i);
      }
    }
    return res;
  }

  /**
   * @brief Get the virtual orbitals of a given spin.
   * @param SpinType spin of the orbital to be checked.
   * @return std::vector<int>: vector with the indexes of the virtual orbitals.
   */
  std::vector<int> getAllVirtual(SpinComponent SpinType) const {
    const auto& onv = SpinType == SpinComponent::Alpha ? alphaOnv_ : betaOnv_;
    std::vector<int> res;
    res.reserve(onv.countVirtual());
    for (int i = 0; i < onv.size(); i++) {
      if (!onv.isOccupied(i)) {
        res.push_back(i);
      }
    }
    return res;
  }

  /**
   * @brief Excites the i-th orbital. If i is empty, raise an ExcitationException
   * @param excitation SingleExcitation object that indentifies the excitation
   */
  void applyExcitation(const SingleExcitation& excitation) {
    auto& ref_orbital = (excitation.orbitalSpin == SpinComponent::Alpha) ? alphaOnv_ : betaOnv_;
    if (ref_orbital.isOccupied(excitation.virtualOrbital) || !ref_orbital.isOccupied(excitation.occupiedOrbital)) {
      throw ExcitationException();
    }
    ref_orbital.flip(excitation.occupiedOrbital);
    ref_orbital.flip(excitation.virtualOrbital);
  };

  /**
   * @brief Excites the i-th orbital. If i is empty, raise an ExcitationException. Also for cross excitations.
   * @param excitation SingleExcitation object that indentifies the excitation
   */
  void applyExcitation(const CrossExcitation& excitation) {
    auto& occ_orbital = (excitation.orbitalSpins.first == SpinComponent::Alpha) ? alphaOnv_ : betaOnv_;
    auto& vir_orbital = (excitation.orbitalSpins.second == SpinComponent::Alpha) ? alphaOnv_ : betaOnv_;
    if (vir_orbital.isOccupied(excitation.virtualOrbital) || !occ_orbital.isOccupied(excitation.occupiedOrbital)) {
      throw ExcitationException();
    }
    occ_orbital.flip(excitation.occupiedOrbital);
    vir_orbital.flip(excitation.virtualOrbital);
  };

  /**
   * @brief Excites the i-th orbital. If i is empty, raise an ExcitationException
   *
   * Overload of the previous function where the list provided in input is used to map the indexes contained in the
   * excitation object. This means that, if this version of the subroutine is called, the indexes contained in
   * excitation must be considered as relative, i.e., if we find 2 as occupied orbital, we must excite the second
   * orbital among the occupied ones.
   *
   * @param excitation SingleExcitation object that identifies the excitation
   * @param lstOccupiedAlpha list with the occupied alpha orbitals
   * @param lstOccupiedBeta list with the occupied beta orbitals
   * @param lstVirtualAlpha list with the virtual alpha orbitals
   * @param lstVirtualBeta list with the virtual beta orbitals
   */
  void applyExcitation(const SingleExcitation& excitation, const std::vector<int>& lstOccupiedAlpha,
                       const std::vector<int>& lstOccupiedBeta, const std::vector<int>& lstVirtualAlpha,
                       const std::vector<int>& lstVirtualBeta) {
    auto& ref_orbital = (excitation.orbitalSpin == SpinComponent::Alpha) ? alphaOnv_ : betaOnv_;
    int iVirtual = (excitation.orbitalSpin == SpinComponent::Alpha) ? lstVirtualAlpha[excitation.virtualOrbital]
                                                                    : lstVirtualBeta[excitation.virtualOrbital];
    int iOccupied = (excitation.orbitalSpin == SpinComponent::Alpha) ? lstOccupiedAlpha[excitation.occupiedOrbital]
                                                                     : lstOccupiedBeta[excitation.occupiedOrbital];
    if (ref_orbital.isOccupied(iVirtual) || !ref_orbital.isOccupied(iOccupied)) {
      throw ExcitationException();
    }
    ref_orbital.flip(iOccupied);
    ref_orbital.flip(iVirtual);
  };

  /**
   * @brief Getter for the number of electrons of the determinant.
   * @return std::pair<int, int> number of alpha and beta electrons.
   */
  std::pair<int, int> countElectrons() const {
    return {alphaOnv_.countOccupied(), betaOnv_.countOccupied()};
  }

  bool isCompatibleWith(const ElectronicDeterminant& other_det) const {
    bool is_compatible = true;
    is_compatible = is_compatible && size(SpinComponent::Alpha) == other_det.size(SpinComponent::Alpha) &&
                    size(SpinComponent::Beta) == other_det.size(SpinComponent::Beta) &&
                    countElectrons().first == other_det.countElectrons().first &&
                    countElectrons().second == other_det.countElectrons().second;
    return is_compatible;
  }

  /**
   * @brief Returns the excitation degree wrt another determinant.
   * @param other_onv: Determinant object wrt to which the excitation degree is calculated.
   * @param checkCompatibility: If true, checks if the determinants are compatible.
   * @return int: degree of excitation, stored as pair where the first element is the number of alpha
   *              excitations and the second element is the number of beta excitations.
   */
  std::pair<int, int> getExcitationDegree(const ElectronicDeterminant& other_det, bool checkCompatibility = false) const {
    if (checkCompatibility) {
      assert(isCompatibleWith(other_det));
    }
    int exc_alpha = 0;
    int exc_beta = 0;
    for (int idx = 0; idx < size(SpinComponent::Alpha); idx++) {
      if (isOccupied(SpinComponent::Alpha, idx) != other_det.isOccupied(SpinComponent::Alpha, idx)) {
        exc_alpha++;
      }
    }
    for (int idx = 0; idx < size(SpinComponent::Beta); idx++) {
      if (isOccupied(SpinComponent::Beta, idx) != other_det.isOccupied(SpinComponent::Beta, idx)) {
        exc_alpha++;
      }
    }
    // Here the division by 2 is safe if we checked the compatibility before (i.e., exc_alpha and exc_beta will
    // be for sure even numbers).
    return std::make_pair(exc_alpha / 2, exc_beta / 2);
  };

  /**
   * @brief Returns the excitation degree wrt another determinant.
   * @throw This throws a std::runtime_error if checkCompatibility is set to true and the determinants are not
   * compatible.
   * @param other_onv: Determinant object wrt to which the excitation degree is calculated.
   * @param checkCompatibility: If true, checks if the determinants are compatible.
   * @return int: degree of excitation, stored as pair where the first element is the number of alpha
   *              excitations and the second element is the number of beta excitations.
   */
  Excitation getExcitations(const ElectronicDeterminant& other_det, bool checkCompatibility = false) const {
    // Preliminary check
    if (checkCompatibility) {
      if (!isCompatibleWith(other_det)) {
        throw std::runtime_error("Determinants not compatible");
      }
    }
    std::vector<int> occupied_alpha, occupied_beta, virtual_alpha, virtual_beta;
    Excitation excitation;
    // Determines alpha excitations
    for (int orbital = 0; orbital < size(SpinComponent::Alpha); orbital++) {
      if (isOccupied(SpinComponent::Alpha, orbital) == 1 && other_det.isOccupied(SpinComponent::Alpha, orbital) == 0) {
        occupied_alpha.push_back(orbital);
      }
      else if (isOccupied(SpinComponent::Alpha, orbital) == 0 && other_det.isOccupied(SpinComponent::Alpha, orbital) == 1) {
        virtual_alpha.push_back(orbital);
      }
    }
    // Determines beta excitations
    for (int orbital = 0; orbital < size(SpinComponent::Beta); orbital++) {
      if (isOccupied(SpinComponent::Beta, orbital) == 1 && other_det.isOccupied(SpinComponent::Beta, orbital) == 0) {
        occupied_beta.push_back(orbital);
      }
      else if (isOccupied(SpinComponent::Beta, orbital) == 0 && other_det.isOccupied(SpinComponent::Beta, orbital) == 1) {
        virtual_beta.push_back(orbital);
      }
    }
    // Generates the final results
    for (int i_exc = 0; i_exc < int(occupied_alpha.size()); i_exc++) {
      excitation.createExcitation(occupied_alpha[i_exc], virtual_alpha[i_exc], SpinComponent::Alpha);
    }
    for (int i_exc = 0; i_exc < int(occupied_beta.size()); i_exc++) {
      excitation.createExcitation(occupied_beta[i_exc], virtual_beta[i_exc], SpinComponent::Beta);
    }
    return excitation;
  };

  /**
   * @brief Comparison operator
   */
  bool operator==(const ElectronicDeterminant& other) const {
    return (getOccupation(SpinComponent::Alpha) == other.getOccupation(SpinComponent::Alpha)) &&
           (getOccupation(SpinComponent::Beta) == other.getOccupation(SpinComponent::Beta));
  }
  /**
   * @brief Comparison operator >
   */
  bool operator>(const ElectronicDeterminant& other) const {
    return (getOccupation(SpinComponent::Alpha) > other.getOccupation(SpinComponent::Alpha)) ||
           ((getOccupation(SpinComponent::Alpha) == other.getOccupation(SpinComponent::Alpha)) &&
            (getOccupation(SpinComponent::Beta) > other.getOccupation(SpinComponent::Beta)));
  }
  bool operator<(const ElectronicDeterminant& other) const {
    return (getOccupation(SpinComponent::Alpha) < other.getOccupation(SpinComponent::Alpha)) ||
           ((getOccupation(SpinComponent::Alpha) == other.getOccupation(SpinComponent::Alpha)) &&
            (getOccupation(SpinComponent::Beta) < other.getOccupation(SpinComponent::Beta)));
  }

  void print(std::ostream& out, SpinComponent orbType) const {
    auto const& orbString = getOccupation(orbType);
    std::string orbName = orbType == SpinComponent::Alpha ? "Alpha" : "Beta";
    out << orbName + " ";
    for (auto orb : orbString) {
      out << orb;
    }
  }

  void print(std::ostream& out) const {
    print(out, SpinComponent::Alpha);
    out << "  ";
    print(out, SpinComponent::Beta);
  }

 private:
  /**
   * @brief Copy constructor of a determinant.
   */
  // Class members
  OccupationNumberVector alphaOnv_, betaOnv_;
};

struct ElectronicDeterminantHash {
  std::size_t operator()(const ElectronicDeterminant& key) const {
    std::size_t hash1 = hash_(key.getOccupation(SpinComponent::Alpha));
    std::size_t hash2 = hash_(key.getOccupation(SpinComponent::Beta));
    return hash1 + hash2;
  }

 private:
  std::hash<std::vector<bool>> hash_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_ELECTRONIC_DETERMINANT_H
