/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Properties/Thermochemistry/DegreesOfFreedomCollection.h>
#include <Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

std::shared_ptr<DegreesOfFreedomCollection> getDegreesOfFreedomPtr(std::shared_ptr<DegreesOfFreedomCollection> degreesOfFreedom) {
  return std::make_shared<DegreesOfFreedomCollection>(*degreesOfFreedom);
}

void init_degrees_of_freedom_collection(pybind11::module& m) {
  pybind11::class_<DegreesOfFreedomCollection, std::shared_ptr<DegreesOfFreedomCollection>> degrees_of_freedom(
      m, "DegreesOfFreedomCollection");

  degrees_of_freedom.def("get", &getDegreesOfFreedomPtr);
  degrees_of_freedom.def(pybind11::init<>(), "Initialize an empty degrees of freedom collection.");
  degrees_of_freedom.def("density_of_states", &DegreesOfFreedomCollection::densityOfState);
  degrees_of_freedom.def("sum_of_states", &DegreesOfFreedomCollection::sumOfStates);
  degrees_of_freedom.def("partition_function", &DegreesOfFreedomCollection::partitionFunction);
  degrees_of_freedom.def("enthalpy", &DegreesOfFreedomCollection::enthalpy);
  degrees_of_freedom.def("entropy", &DegreesOfFreedomCollection::entropy);
  degrees_of_freedom.def("helmholtz_free_energy", &DegreesOfFreedomCollection::helmholtzFreeEnergy);
  degrees_of_freedom.def("heat_capacity_Cp", &DegreesOfFreedomCollection::heatCapacityCp);
  degrees_of_freedom.def("heat_capacity_Cv", &DegreesOfFreedomCollection::heatCapacityCv);
}