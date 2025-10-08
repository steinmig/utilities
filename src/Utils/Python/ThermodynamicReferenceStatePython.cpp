/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h>
#include <pybind11/pybind11.h>

using namespace Scine::Utils;

ThermodynamicReferenceState getStandardStateGas() {
  return ReferenceStates::StandardStateGas;
}
ThermodynamicReferenceState getStandardStateLiquid() {
  return ReferenceStates::StandardStateLiquid;
}
ThermodynamicReferenceState getVacuumZeroKelvin() {
  return ReferenceStates::VacuumZeroKelvin;
}

void init_thermodynamic_reference_state(pybind11::module& m) {
  pybind11::class_<ThermodynamicReferenceState> thermodynamic_reference_state(m, "ThermodynamicReferenceState");

  thermodynamic_reference_state.def(pybind11::init<double, double>(), pybind11::arg("temperature"), pybind11::arg("pressure"),
                                    "The thermodynamic reference state for thermochemical calculations.");
  thermodynamic_reference_state.def_readwrite("pressure", &ThermodynamicReferenceState::pressure);
  thermodynamic_reference_state.def_readwrite("temperature", &ThermodynamicReferenceState::temperature);
  thermodynamic_reference_state.def("__eq__", &ThermodynamicReferenceState::operator==);
  thermodynamic_reference_state.def(pybind11::pickle(
      [](const ThermodynamicReferenceState& referenceState) { // __getstate__
        return pybind11::make_tuple(referenceState.temperature, referenceState.temperature);
      },
      [](pybind11::tuple tuple) { // __setstate__
        if (tuple.size() != 2) {
          throw std::runtime_error("Unable to restore reference state from pickle. Invalid state.");
        }
        return ThermodynamicReferenceState(tuple[0].cast<double>(), tuple[1].cast<float>());
      }));
}

void init_standard_states(pybind11::module& m) {
  m.def("standard_state_gas", &getStandardStateGas);
  m.def("standard_state_liquid", &getStandardStateLiquid);
  m.def("vacuum_zero_kelvin", &getVacuumZeroKelvin);
}