/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CASINDICESHANDLER_H
#define UTILS_CASINDICESHANDLER_H

#include "CICalculatorSettings.h"
#include "CasSpecifier.h"

namespace Scine {
namespace Utils {
namespace SecondQuantization {

struct CasIndicesHandler {
  CasIndicesHandler(const CasSpecifier<Utils::Reference::Restricted>& cas) : cas_(cas) {
  }

  /**
   * @brief Get the indices in the CAS from a number of orbitals around the Fermi level.
   * The indices refer to the CAS space, not the original space.
   * Special values: there are special values leading to special behavior.
   *
   * - DoubleExcitationSpecifier::FullSpace (aka -1) leads to all the indices
   *   to be included
   * - DoubleExcitationSpecifier::OnlySingleExcitations (aka 0) and
   *   ReferenceSpecifiers::SingleReference (aka 0) lead to an empty index vector.
   */
  auto getIndices(int orbitalsAroundFermi) const -> std::vector<int> {
    if (orbitalsAroundFermi == DoubleExcitationSpecifiers::FullSpace) {
      return getIndices(cas_.getActiveIndices().restricted);
    }
    if (orbitalsAroundFermi == ReferenceSpecifiers::SingleReference ||
        orbitalsAroundFermi == DoubleExcitationSpecifiers::OnlySingleExcitations) {
      return {};
    }
    int nElectrons = cas_.getCasElectrons().restricted;
    int homoIndex = CasGenerator::findHomo(nElectrons);
    std::vector<int> indices(orbitalsAroundFermi);

    std::iota(indices.begin(), indices.end(), homoIndex + 1 - orbitalsAroundFermi / 2);
    return indices;
  }

  /**
   * @brief Get the indices in the CAS from indices in the full space.
   * The indices refer to the CAS space, not the original space.
   */
  auto getIndices(const std::vector<int>& indices) const -> std::vector<int> {
    std::vector<int> result;
    result.reserve(indices.size());

    std::transform(indices.begin(), indices.end(), std::back_inserter<std::vector<int>>(result), [this](int idx) {
      const auto& actives = cas_.getActiveIndices().restricted;
      auto it = std::find(actives.begin(), actives.end(), idx);
      if (it == actives.end()) {
        throw std::runtime_error("Indices '" + std::to_string(idx) + "' not available in the cas");
      }
      int distance = std::distance(actives.begin(), it);
      return distance;
    });
    return result;
  }

 private:
  const CasSpecifier<Utils::Reference::Restricted>& cas_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_CASINDICESHANDLER_H
