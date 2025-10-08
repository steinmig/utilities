/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CICALCULATORSETTINGS_H
#define UTILS_CICALCULATORSETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

struct DoubleExcitationSpecifiers {
  static constexpr int FullSpace = -1;
  static constexpr int OnlySingleExcitations = 0;
};
struct ReferenceSpecifiers {
  static constexpr int SingleReference = 0;
};
struct CasSpecifiers {
  static constexpr int FullSpace = std::numeric_limits<int>::max();
};

static constexpr const char* convergence = "convergence";

class CICalculatorSettings : public Settings {
 public:
  CICalculatorSettings() : Settings("Settings for the CI Calculator") {
    /// CAS options

    // UNO CI
    Utils::UniversalSettings::BoolDescriptor autoCas("Determine the reference determinants with UNO-CAS.");
    autoCas.setDefaultValue(false);
    _fields.push_back(SettingsNames::casWithUnoCI, std::move(autoCas));

    Utils::UniversalSettings::DoubleDescriptor minimalOccupationUnoCi(
        "Occupation threshold for UNO-CI reference determinants detection."
        "The higher boundary is 2 - low_boundary.");
    minimalOccupationUnoCi.setMinimum(0);
    minimalOccupationUnoCi.setMaximum(1);
    minimalOccupationUnoCi.setDefaultValue(0.02);
    _fields.push_back(SettingsNames::unoCiLowThresholdOption, std::move(minimalOccupationUnoCi));

    // CAS space total
    Utils::UniversalSettings::IntDescriptor casOrbitalsAroundFermi("Number of orbitals around the Fermi level "
                                                                   "to include in the cas");
    casOrbitalsAroundFermi.setDefaultValue(CasSpecifiers::FullSpace);
    casOrbitalsAroundFermi.setMinimum(2);
    _fields.push_back(SettingsNames::casOrbitalsAroundFermiOption, std::move(casOrbitalsAroundFermi));

    Utils::UniversalSettings::StringDescriptor casOrbitalIndices(
        "Comma-separated list of orbital indices from which to build the reference functions.");
    casOrbitalIndices.setDefaultValue("");
    _fields.push_back(SettingsNames::casOrbitalIndicesOption, std::move(casOrbitalIndices));

    // MR treatment in this space
    Utils::UniversalSettings::IntDescriptor referenceCasOrbitalsAroundFermi("Number of orbitals around the Fermi level "
                                                                            "to include in the CAS");
    referenceCasOrbitalsAroundFermi.setDefaultValue(ReferenceSpecifiers::SingleReference);
    referenceCasOrbitalsAroundFermi.setMinimum(ReferenceSpecifiers::SingleReference);
    _fields.push_back(SettingsNames::referencesAroundFermiOption, std::move(referenceCasOrbitalsAroundFermi));

    Utils::UniversalSettings::StringDescriptor referenceCasOrbitalIndices(
        "Comma-separated list of orbital indices from which to build the reference functions.");
    referenceCasOrbitalIndices.setDefaultValue("");
    _fields.push_back(SettingsNames::referencesIndicesOption, std::move(referenceCasOrbitalIndices));

    Utils::UniversalSettings::BoolDescriptor referenceFromDoubleExcitations(
        "Whether to generate the reference determinants from double symmetric double excitations from the HF "
        "determinant in an active space.");
    referenceFromDoubleExcitations.setDefaultValue(false);
    _fields.push_back(SettingsNames::referenceFromDoubleExcitations, std::move(referenceFromDoubleExcitations));
    /// CAS options end
    /// Excitations options
    Utils::UniversalSettings::IntDescriptor excitationCasOrbitalsAroundFermi(
        "Number of orbitals around the Fermi level within which double excitations are constructed.");
    excitationCasOrbitalsAroundFermi.setDefaultValue(DoubleExcitationSpecifiers::FullSpace);
    excitationCasOrbitalsAroundFermi.setMinimum(DoubleExcitationSpecifiers::OnlySingleExcitations);
    _fields.push_back(SettingsNames::doublesAroundFermiOption, std::move(excitationCasOrbitalsAroundFermi));

    Utils::UniversalSettings::StringDescriptor excitationCasOrbitalIndices(
        "Comma-separated list of orbital indices to include within which double excitations are constructed.");
    excitationCasOrbitalIndices.setDefaultValue("");
    _fields.push_back(SettingsNames::doublesIndicesOption, std::move(excitationCasOrbitalIndices));
    /// Excitations options end

    Utils::UniversalSettings::BoolDescriptor calculateTransitionDipole(
        "Whether to calculate the transition dipole. Might be costly for large CI expansions.");
    calculateTransitionDipole.setDefaultValue(false);
    _fields.push_back("calculate_transition_dipole", std::move(calculateTransitionDipole));

    addDavidsonOptions();
    resetToDefaults();
  }
  void addDavidsonOptions() {
    Utils::UniversalSettings::IntDescriptor numberOfEigenstates("Sets the number of desired Eigenstates.");
    numberOfEigenstates.setMinimum(0);
    numberOfEigenstates.setDefaultValue(0);

    Utils::UniversalSettings::IntDescriptor initialSubspaceDimension(
        "Sets the initial space dimension to use in the Davidson diagonalizer.");
    initialSubspaceDimension.setMinimum(0);
    initialSubspaceDimension.setDefaultValue(0);

    Utils::UniversalSettings::IntDescriptor maxDavidsonIterations(
        "Sets the maximal iteration number for the iterative diagonalizer.");
    maxDavidsonIterations.setMinimum(0);
    maxDavidsonIterations.setDefaultValue(0);

    Utils::UniversalSettings::DoubleDescriptor convergenceCriterionDavidson(
        "The convergence threshold for the Davidson solver.");
    convergenceCriterionDavidson.setDefaultValue(1.0e-5);
    convergenceCriterionDavidson.setMinimum(std::numeric_limits<double>::min());

    Utils::UniversalSettings::OptionListDescriptor gepAlgorithm("Algorithm to compute the stable Generalize"
                                                                "Eigenvalue Problem Ax=lBx when B is almost singular.");
    gepAlgorithm.addOption("standard");
    gepAlgorithm.addOption("cholesky");
    gepAlgorithm.addOption("simultaneous_diag");
    gepAlgorithm.setDefaultOption("standard");

    Utils::UniversalSettings::BoolDescriptor directCi("Determine the reference determinants with UNO-CAS.");
    directCi.setDefaultValue(true);
    _fields.push_back(SettingsNames::directness, std::move(directCi));

    _fields.push_back(Utils::SettingsNames::numberOfEigenstates, std::move(numberOfEigenstates));
    _fields.push_back(Utils::SettingsNames::initialSubspaceDimension, std::move(initialSubspaceDimension));
    _fields.push_back(Utils::SettingsNames::maxDavidsonIterations, std::move(maxDavidsonIterations));
    _fields.push_back(convergence, std::move(convergenceCriterionDavidson));
    _fields.push_back("gep_algo", std::move(gepAlgorithm));
  }
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_UNOCI_H
