/**
 * @file ThermodynamicReferenceState.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ThermodynamicReferenceState.h"
#include <cmath>

namespace Scine {
namespace Utils {

ThermodynamicReferenceState::ThermodynamicReferenceState(double t, double p) : temperature(t), pressure(p) {
}
bool ThermodynamicReferenceState::operator==(const ThermodynamicReferenceState& other) const {
  return std::abs(temperature - other.temperature) < 1e-9 && std::abs(pressure - other.pressure) < 1e-9;
}
ThermodynamicReferenceState::~ThermodynamicReferenceState() = default;

} // namespace Utils
} // namespace Scine