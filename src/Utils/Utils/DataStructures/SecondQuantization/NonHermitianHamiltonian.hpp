/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

namespace Scine {
namespace Utils {
namespace SecondQuantization {

template<typename T>
typename NonHermitianHamiltonian<T>::TwoBodyIdentifier NonHermitianHamiltonian<T>::bringToNormalFormImpl(const TwoBodyIdentifier& input_ob_identifier) const {
  auto sorted_array = input_ob_identifier;
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

}
}
}