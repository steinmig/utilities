/**
 * @file ThermodynamicReferenceState.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_THERMODYNAMICREFERENCESTATE_H
#define UTILSOS_THERMODYNAMICREFERENCESTATE_H

#include "Utils/Constants.h"

namespace Scine {
namespace Utils {
class ThermodynamicReferenceState {
 public:
  ThermodynamicReferenceState(double temperature, double pressure);
  ~ThermodynamicReferenceState();

  bool operator==(const ThermodynamicReferenceState& other) const;

  double temperature;
  double pressure;
};
namespace ReferenceStates {
constexpr double referenceTemperature = 298.15;
static const ThermodynamicReferenceState StandardStateGas(referenceTemperature, 101325.0);

constexpr double oneMolPerLiter = Constants::molarGasConstant * referenceTemperature * 1e+3;
static const ThermodynamicReferenceState StandardStateLiquid(referenceTemperature, oneMolPerLiter);

static const ThermodynamicReferenceState VacuumZeroKelvin(0.0, 0.0);
} // namespace ReferenceStates

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_THERMODYNAMICREFERENCESTATE_H
