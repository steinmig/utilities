/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

std::shared_ptr<MolecularDegreesOfFreedom>
getMolecularDegreesOfFreedomPtr(std::shared_ptr<MolecularDegreesOfFreedom> molecularDegreesOfFreedom) {
  return molecularDegreesOfFreedom->shared_from_this();
}

void init_molecular_degrees_of_freedom(pybind11::module& m) {
  pybind11::class_<MolecularDegreesOfFreedom, std::shared_ptr<MolecularDegreesOfFreedom>> molecular_degrees_of_freedom(
      m, "MolecularDegreesOfFreedom");

  molecular_degrees_of_freedom.def("get", &getMolecularDegreesOfFreedomPtr);
  molecular_degrees_of_freedom.def(
      pybind11::init<const HessianMatrix&, const AtomCollection&, int, double, unsigned int, bool>(),
      pybind11::arg("hessian"), pybind11::arg("atoms"), pybind11::arg("multiplicity"),
      pybind11::arg("electronic_energy"), pybind11::arg("symmetry_number") = 1,
      pybind11::arg("classical_vibrations") = false, "Initialize a molecular degrees of freedom object.");
  molecular_degrees_of_freedom.def("density_of_states", &MolecularDegreesOfFreedom::densityOfState);
  molecular_degrees_of_freedom.def("sum_of_states", &MolecularDegreesOfFreedom::sumOfStates);
  molecular_degrees_of_freedom.def("partition_function", &MolecularDegreesOfFreedom::partitionFunction);
  molecular_degrees_of_freedom.def("enthalpy", &MolecularDegreesOfFreedom::enthalpy);
  molecular_degrees_of_freedom.def("entropy", &MolecularDegreesOfFreedom::entropy);
  molecular_degrees_of_freedom.def("helmholtz_free_energy", &MolecularDegreesOfFreedom::helmholtzFreeEnergy);
  molecular_degrees_of_freedom.def("heat_capacity_Cp", &MolecularDegreesOfFreedom::heatCapacityCp);
  molecular_degrees_of_freedom.def("heat_capacity_Cv", &MolecularDegreesOfFreedom::heatCapacityCv);
  molecular_degrees_of_freedom.def("get_rrkm_degrees_of_freedom", &MolecularDegreesOfFreedom::getRRKMDegreesOfFreedom,
                                   pybind11::arg("include_rotation") = true);
  molecular_degrees_of_freedom.def("get_normal_modes_container", &MolecularDegreesOfFreedom::getNormalModesContainer);
}