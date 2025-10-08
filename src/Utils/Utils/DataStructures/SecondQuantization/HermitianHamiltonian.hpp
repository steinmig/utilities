/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

namespace Scine {
namespace Utils {
namespace SecondQuantization {

inline HermitianHamiltonian::TwoBodyIdentifier
HermitianHamiltonian::bringToNormalFormImpl(const TwoBodyIdentifier& input_ob_identifier) const {
  auto sorted_array = input_ob_identifier;
  // <pq|rs> : p > r -> <rq|ps>
  if (sorted_array[0] > sorted_array[2]) {
    std::swap(sorted_array[0], sorted_array[2]);
  }
  // <pq|rs> : q > s -> <ps|rq>
  if (sorted_array[1] > sorted_array[3]) {
    std::swap(sorted_array[1], sorted_array[3]);
  }
  // <pq|rs> : p > q -> <qp|sr>
  if (sorted_array[0] > sorted_array[1]) {
    std::swap(sorted_array[0], sorted_array[1]);
    std::swap(sorted_array[2], sorted_array[3]);
  }
  // <pp|rs> : p = q, r > s -> <pp|sr>
  else if (sorted_array[0] == sorted_array[1] && sorted_array[2] > sorted_array[3]) {
    std::swap(sorted_array[2], sorted_array[3]);
  }
  return sorted_array;
}

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
