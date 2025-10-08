/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Properties/ReactionRates/RRKMRateConstantCalculator.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

using namespace Scine::Utils;

void init_rrkm_rate_constant_calculator(pybind11::module& m) {
  pybind11::class_<RRKMRateConstantCalculator> rrkm_rate_constant_calculator(m, "RRKMRateConstantCalculator");
  rrkm_rate_constant_calculator.def(
      pybind11::init<std::shared_ptr<MolecularDegreesOfFreedom>, std::shared_ptr<MolecularDegreesOfFreedom>,
                     std::shared_ptr<MolecularDegreesOfFreedom>, bool>(),
      pybind11::arg("reactant_degrees"), pybind11::arg("ts_degrees"), pybind11::arg("product_degrees") = nullptr,
      pybind11::arg("active_rotor") = false,
      "Constructor from degrees of freedom.\n"
      "reactant_degrees: reactant molecular degrees of freedom\n"
      "ts_degrees: transition state molecular degrees of freedom\n"
      "product_degrees: product molecular degrees of freedom\n"
      "active_rotor: If true, the rotational degrees of freedom are considered active.");
  rrkm_rate_constant_calculator.def("get_rate_constants", &RRKMRateConstantCalculator::getRateConstants,
                                    pybind11::arg("energies"));
}
