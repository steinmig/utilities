/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/DataStructures/SecondQuantization/CICalculator.h>
#include <Utils/DataStructures/SecondQuantization/ElectronicDeterminant.h>
#include <Utils/DataStructures/SecondQuantization/ExcitationHelpers.h>
#include <Utils/DataStructures/SecondQuantization/OccupationNumberVector.h>
#include <Utils/Pybind.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

void init_electronic_determinant(pybind11::module& m);
void init_ci_calculator(pybind11::module& m);

void init_sq_submodule(pybind11::module& m) {
  pybind11::module secondQuantization = m.def_submodule("sq");
  secondQuantization.doc() = R"(
    The ``sq`` submodule contains data structures and algorithms to handle problems posed in
    the second quantization formalism. Data structure examples are wavefunctions, electronic determinants,
    whereas algorithms are, for instance, a CI optimizer.
  )";

  init_electronic_determinant(secondQuantization);
  init_ci_calculator(secondQuantization);
}

void init_ci_calculator(pybind11::module& m) {
  pybind11::class_<Scine::Utils::SecondQuantization::CICalculator, Scine::Core::CalculatorWithReference,
                   std::shared_ptr<Scine::Utils::SecondQuantization::CICalculator>>
      ciCalculator(m, "CICalculator",
                   R"delim(
            A :class:`.CalculatorWithReference` decorated by methods to recover additional information.

            The CICalculator takes a :class:`.Calculator` providing a SQSpecifier in the results
            and calculates a CI expansion determined by the settings provided.
            The optimization can be done directly or, for small CI expansions, by
            constructing the CI matrix and diagonalizing it with the Davidson algorithm.

          )delim");
  ciCalculator.def(pybind11::init<>());
  ciCalculator.def_static("labels", &Scine::Utils::SecondQuantization::CICalculator::getLabels,
                          R"delim(
      Returns a list of electronic occupation as strings, one for each :class:`ElectronicDeterminant` in
      the input list.
      )delim");
  ciCalculator.def(
      "wavefunction",
      [](const Scine::Utils::SecondQuantization::CICalculator& self, int state) { return self.getWavefunction(state); },
      R"delim(
      Returns a state represented by determinant-coefficient pairs.

      The determinants are expressed as :math:`|\phi_{\alpha,1} \phi_{\beta,1} \phi_{\alpha,2}...>`
      with :math:`\phi` a molecular orbital, :math:`\alpha` and :math:`\beta`
      are spins and only occupied molecular orbitals are in the string.

  )delim");
}

void init_electronic_determinant(pybind11::module& m) {
  pybind11::enum_<Scine::Utils::SecondQuantization::SpinComponent>(m, "SpinComponent", pybind11::arithmetic(),
                                                                   "Enumeration type to list the spin, alpha or beta.")
      .value("Alpha", Scine::Utils::SecondQuantization::SpinComponent::Alpha, "Alpha spin.")
      .value("Beta", Scine::Utils::SecondQuantization::SpinComponent::Beta, "Beta spin.")
      .export_values();

  pybind11::class_<Scine::Utils::SecondQuantization::SingleExcitation>(m, "SingleExcitation",
                                                                       "Object representing a single excitation.")
      .def(pybind11::init<int, int, Scine::Utils::SecondQuantization::SpinComponent>())
      .def_readwrite("occ", &Scine::Utils::SecondQuantization::SingleExcitation::occupiedOrbital, "The index of an occupied orbital.")
      .def_readwrite("vir", &Scine::Utils::SecondQuantization::SingleExcitation::virtualOrbital, "The index of a virtual orbital.")
      .def_readwrite("spin", &Scine::Utils::SecondQuantization::SingleExcitation::orbitalSpin, "The spin of the excitation.");

  pybind11::class_<Scine::Utils::SecondQuantization::Excitation>(m, "Excitation", "A container of single excitations.")
      .def_property_readonly("cardinality", &Scine::Utils::SecondQuantization::Excitation::cardinality, "The number of excitations.")
      .def("__getitem__",
           [](const Scine::Utils::SecondQuantization::Excitation& self, int i) { return self.getSingleExcitation(i); });

  pybind11::class_<Scine::Utils::SecondQuantization::OccupationNumberVector>(m, "OccupationNumberVector")
      .def(pybind11::init<std::vector<int>>())
      .def("count_occupied", &Scine::Utils::SecondQuantization::OccupationNumberVector::countOccupied,
           "Counts the occupied orbitals in the occupation number vector.")
      .def("count_virtual", &Scine::Utils::SecondQuantization::OccupationNumberVector::countVirtual,
           "Counts the virtual orbitals in the occupation number vector.")
      .def("is_occupied", &Scine::Utils::SecondQuantization::OccupationNumberVector::isOccupied,
           "Checks whether given orbital is occupied or virtual. Throws runtime error if given orbital doesn't exist.")
      .def("flip", &Scine::Utils::SecondQuantization::OccupationNumberVector::flip,
           "Flips the i-th bit of the occupation number vector. 0 indexed.")
      .def("size", &Scine::Utils::SecondQuantization::OccupationNumberVector::size,
           "Size of the occupation number vector = number of orbitals.")
      .def(pybind11::self == pybind11::self)
      .def_property_readonly("onv", &Scine::Utils::SecondQuantization::OccupationNumberVector::getOccupation,
                             "The representation of the occupation number vector.")
      .def("__repr__",
           [](pybind11::handle h) -> std::string {
             const auto& name = Scine::Utils::qualifiedName(h);
             const Scine::Utils::SecondQuantization::OccupationNumberVector& onv =
                 h.cast<Scine::Utils::SecondQuantization::OccupationNumberVector>();
             if (onv.size() == 0) {
               return name + "()";
             }
             std::string onvStr;
             for (bool orb : onv.getOccupation()) {
               onvStr += std::to_string(int(orb));
             }

             return "<" + name + " (" + onvStr + ") instance>";
           })
      .def("__deepcopy__",
           [](const Scine::Utils::SecondQuantization::OccupationNumberVector& self, const pybind11::dict&
              /* memo */) -> Scine::Utils::SecondQuantization::OccupationNumberVector { return self; });
  ;

  pybind11::class_<Scine::Utils::SecondQuantization::ElectronicDeterminant>(m, "ElectronicDeterminant",
                                                                            R"delim(
      A class representing an electronic determinant as a series of alpha and beta orbitals.
        )delim")
      .def(pybind11::init<>())
      .def(pybind11::init<const Scine::Utils::SecondQuantization::OccupationNumberVector::InputType&,
                          const Scine::Utils::SecondQuantization::OccupationNumberVector::InputType&>())
      .def(pybind11::init<Scine::Utils::SecondQuantization::OccupationNumberVector, Scine::Utils::SecondQuantization::OccupationNumberVector>())
      .def(pybind11::init<const std::string&>())
      .def("onv", [](const Scine::Utils::SecondQuantization::ElectronicDeterminant& self,
                     Scine::Utils::SecondQuantization::SpinComponent spin) { return self.getOnv(spin); })
      .def("size",
           pybind11::overload_cast<Scine::Utils::SecondQuantization::SpinComponent>(
               &Scine::Utils::SecondQuantization::ElectronicDeterminant::size, pybind11::const_),
           "The size of the ONV indicated by the spin")
      .def("is_occupied", &Scine::Utils::SecondQuantization::ElectronicDeterminant::isOccupied,
           "Checks whether the i-th orbital of the given spin is occupied.")
      .def("occupied", &Scine::Utils::SecondQuantization::ElectronicDeterminant::getAllOccupied,
           "Returns a list of indices of the occupied orbitals of the given spin.")
      .def("virtual", &Scine::Utils::SecondQuantization::ElectronicDeterminant::getAllVirtual,
           "Returns a list of indices of the virtual orbitals of the given spin.")
      .def(
          "apply_excitation",
          [](Scine::Utils::SecondQuantization::ElectronicDeterminant& self, int occ, int vir,
             Scine::Utils::SecondQuantization::SpinComponent spin) {
            self.applyExcitation({occ, vir, spin});
          },
          "Moves an electron into a virtual orbital")
      .def(
          "apply_cross_excitation",
          [](Scine::Utils::SecondQuantization::ElectronicDeterminant& self, int occ, int vir,
             Scine::Utils::SecondQuantization::SpinComponent spin1, Scine::Utils::SecondQuantization::SpinComponent spin2) {
            self.applyExcitation({occ, vir, {spin1, spin2}});
          },
          "Moves an electron of spin ``spin1`` into a virtual orbital of spin ``spin2``.")
      .def("count_electrons", &Scine::Utils::SecondQuantization::ElectronicDeterminant::countElectrons,
           "Counts the electrons in the determinant. First number is alpha, second is beta.")
      .def("excitations", &Scine::Utils::SecondQuantization::ElectronicDeterminant::getExcitations,
           "Returns a list of excitations from another determinant.")
      .def(pybind11::self == pybind11::self)
      .def(pybind11::self < pybind11::self)
      .def(pybind11::self > pybind11::self)
      .def("__repr__",
           [](pybind11::handle h) -> std::string {
             const auto& name = Scine::Utils::qualifiedName(h);
             const Scine::Utils::SecondQuantization::ElectronicDeterminant& det =
                 h.cast<Scine::Utils::SecondQuantization::ElectronicDeterminant>();
             if (det.size() == 0) {
               return name + "()";
             }
             std::stringstream stream;
             det.print(stream);

             return "<" + name + " (" + stream.str() + ") instance>";
           })
      .def("__deepcopy__",
           [](const Scine::Utils::SecondQuantization::ElectronicDeterminant& self, const pybind11::dict&
              /* memo */) -> Scine::Utils::SecondQuantization::ElectronicDeterminant { return self; })
      .def("__hash__", [](const Scine::Utils::SecondQuantization::ElectronicDeterminant& self) {
        return Scine::Utils::SecondQuantization::ElectronicDeterminantHash{}(self);
      });
}
