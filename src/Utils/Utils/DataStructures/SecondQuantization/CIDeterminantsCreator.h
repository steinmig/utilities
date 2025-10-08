/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CIDETERMINANTSCREATOR_H
#define UTILS_CIDETERMINANTSCREATOR_H

#include "ElectronicDeterminant.h"
#include <boost/math/special_functions/binomial.hpp>
namespace Scine {
namespace Utils {
namespace SecondQuantization {

struct CIDeterminantsCreator {
  CIDeterminantsCreator(ElectronicDeterminant referenceDeterminant) : ref_(std::move(referenceDeterminant)) {
  }

  std::vector<ElectronicDeterminant> generate(const std::vector<int>& indices, bool onlyDoubleExcitations = false) const {
    return onlyDoubleExcitations ? generateOnlyDoubleExcitations(indices) : generateAllReferences(indices);
  }
  /**
   * @brief Generates all reference determinants in an active space determined by indices.
   *
   * For example, with a reference determinant
   *
   * 222000
   *
   * and indices = {2,3}
   *
   * the reference determiants will be
   *
   * 222000, 220200, 22ab00, 22ba00.
   */
  std::vector<ElectronicDeterminant> generateAllReferences(const std::vector<int>& indices) const {
    checkIndices(indices);

    std::map<SpinComponent, std::vector<OccupationNumberVector>> onvs;

    for (auto spin : {SpinComponent::Alpha, SpinComponent::Beta}) {
      const auto& refOnv = ref_.getOccupation(spin);
      std::vector<bool> onvToPermute;
      onvToPermute.reserve(indices.size());
      for (auto idx : indices) {
        onvToPermute.push_back(refOnv[idx]);
      }
      std::sort(onvToPermute.begin(), onvToPermute.end());
      // RESERVE
      int nElectrons = std::count(onvToPermute.begin(), onvToPermute.end(), 1);
      onvs[spin].reserve(boost::math::binomial_coefficient<double>(indices.size(), nElectrons));
      do {
        auto newOnv = refOnv;
        int index = 0;
        for (auto idx : indices) {
          newOnv[idx] = onvToPermute[index++];
        }
        onvs[spin].emplace_back(std::move(newOnv));
      } while (std::next_permutation(onvToPermute.begin(), onvToPermute.end()));
    }

    std::vector<ElectronicDeterminant> refDets;
    refDets.reserve(onvs[SpinComponent::Alpha].size() * onvs[SpinComponent::Beta].size());
    for (const auto& onva : onvs.at(SpinComponent::Alpha)) {
      for (const auto& onvb : onvs.at(SpinComponent::Beta)) {
        refDets.emplace_back(onva, onvb);
      }
    }
    return refDets;
  }
  /**
   * @brief Generates reference determinants created from double excitations from a
   *        reference det in an active space determined by indices.
   *
   * For example, with a reference determinant
   *
   * 222000
   *
   * and indices = {2,3}
   *
   * the reference determiants will be
   *
   * 222000, 220200.
   */
  std::vector<ElectronicDeterminant> generateOnlyDoubleExcitations(const std::vector<int>& indices) const {
    checkIndices(indices);
    std::vector<bool> onvToPermute;
    onvToPermute.reserve(indices.size());
    std::vector<int> validIndices;
    for (int idx : indices) {
      if (ref_.isOccupied(SpinComponent::Alpha, idx) == ref_.isOccupied(SpinComponent::Beta, idx)) {
        validIndices.push_back(idx);
      }
    }
    const auto& refOnvAlpha = ref_.getOccupation(SpinComponent::Alpha);
    const auto& refOnvBeta = ref_.getOccupation(SpinComponent::Beta);
    for (auto idx : validIndices) {
      onvToPermute.push_back(refOnvAlpha[idx]);
    }
    std::sort(onvToPermute.begin(), onvToPermute.end());
    int nElectrons = std::count(onvToPermute.begin(), onvToPermute.end(), 1);

    std::vector<ElectronicDeterminant> refDets;
    // RESERVE
    refDets.reserve(boost::math::binomial_coefficient<double>(validIndices.size(), nElectrons));
    do {
      auto newOnvAlpha = refOnvAlpha;
      auto newOnvBeta = refOnvBeta;
      int index = 0;
      for (auto idx : validIndices) {
        newOnvAlpha[idx] = onvToPermute[index];
        newOnvBeta[idx] = onvToPermute[index++];
      }
      refDets.emplace_back(OccupationNumberVector(std::move(newOnvAlpha)), OccupationNumberVector(std::move(newOnvBeta)));
    } while (std::next_permutation(onvToPermute.begin(), onvToPermute.end()));

    return refDets;
  }

 private:
  void checkIndices(std::vector<int> indices) const {
    if (std::any_of(indices.begin(), indices.end(),
                    [&](int i) { return (i >= ref_.size(SpinComponent::Alpha) || i < 0); })) {
      throw std::runtime_error("Asked references from orbitals not present in active space.");
    }
    auto size = indices.size();
    std::sort(indices.begin(), indices.end());
    auto last = std::unique(indices.begin(), indices.end());
    indices.erase(last, indices.end());
    if (size != indices.size()) {
      throw std::runtime_error("Asked references from duplicated orbitals.");
    }
  }
  ElectronicDeterminant ref_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_UNOCI_H
