/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/DataStructures/SecondQuantization/SQSpecifier.h"
#include <Utils/Pybind.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace Scine::Utils;
using namespace Scine::Utils::SecondQuantization;

void init_sq_specifier(pybind11::module& m) {
  pybind11::class_<Specifier> sq_specifier(m, "Specifier");
  sq_specifier.def(pybind11::init<>());
  sq_specifier.def("__repr__", [](const Specifier&) { return "specifier"; });
  sq_specifier.def_property_readonly("get_one_body_integrals", &Specifier::getOneBody);
  sq_specifier.def_property_readonly("get_one_body_integrals_alpha", &Specifier::getOneBodyAlpha);
  sq_specifier.def_property_readonly("get_one_body_integrals_beta", &Specifier::getOneBodyBeta);
  sq_specifier.def_property_readonly("get_two_body_integrals", &Specifier::getTwoBody);
  sq_specifier.def_property_readonly("get_two_body_integrals_alpha_alpha", &Specifier::getTwoBodyAlphaAlpha);
  sq_specifier.def_property_readonly("get_two_body_integrals_alpha_beta", &Specifier::getTwoBodyAlphaBeta);
  sq_specifier.def_property_readonly("get_two_body_integrals_beta_alpha", &Specifier::getTwoBodyBetaAlpha);
  sq_specifier.def_property_readonly("get_two_body_integrals_beta_beta", &Specifier::getTwoBodyBetaBeta);
}
