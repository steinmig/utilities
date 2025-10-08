/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CasSpecifier.h"
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
namespace Scine {
namespace Utils {
namespace SecondQuantization {
CasSpecifier<Reference::Restricted> CasGenerator::generateRestricted(int totalNumberOfElectrons, std::vector<int> indices) const {
  SpinAdaptedContainer<Reference::Restricted, int> casElectrons;
  SpinAdaptedContainer<Reference::Restricted, std::vector<int>> activeCasIndices;
  SpinAdaptedContainer<Reference::Restricted, std::vector<int>> coreIndices;
  activeCasIndices.restricted = indices;
  // Create HF determinant, first nElec/2 orbitals are occupied
  coreIndices.restricted = std::vector<int>(totalNumberOfElectrons / 2);
  std::iota(coreIndices.restricted.begin(), coreIndices.restricted.end(), 0);
  // Erease all occupied orbitals that are in the CAS
  coreIndices.restricted.erase(remove_if(coreIndices.restricted.begin(), coreIndices.restricted.end(),
                                         [&](auto x) { return find(indices.begin(), indices.end(), x) != indices.end(); }),
                               coreIndices.restricted.end());
  casElectrons.restricted = totalNumberOfElectrons - 2 * coreIndices.restricted.size();

  return {casElectrons, mos_, activeCasIndices, coreIndices};
}

CasSpecifier<Reference::Unrestricted> CasGenerator::generateUnrestricted(int totalNumberElectrons, int spinMultiplicity,
                                                                         std::vector<int> alphaIndices,
                                                                         std::vector<int> betaIndices) const {
  if (alphaIndices.size() != betaIndices.size()) {
    throw std::runtime_error("Different size for alpha and beta spaces in CAS generation.");
  }
  int alphaElectrons = (spinMultiplicity + totalNumberElectrons - 1) / 2;
  int betaElectrons = (totalNumberElectrons - spinMultiplicity + 1) / 2;
  SpinAdaptedContainer<Reference::Unrestricted, std::vector<int>> inactiveCasIndices;

  SpinAdaptedContainer<Reference::Unrestricted, std::vector<int>> activeCasIndices{std::move(alphaIndices),
                                                                                   std::move(betaIndices)};
  SpinAdaptedContainer<Reference::Unrestricted, int> casElectrons{
      static_cast<int>(std::count_if(activeCasIndices.alpha.begin(), activeCasIndices.alpha.end(),
                                     [&](int index) { return index <= alphaElectrons; })),
      static_cast<int>(std::count_if(activeCasIndices.beta.begin(), activeCasIndices.beta.end(),
                                     [&](int index) { return index <= betaElectrons; }))};

  int nCoreAlphaOrbitals = alphaElectrons - casElectrons.alpha;
  int nCoreBetaOrbitals = betaElectrons - casElectrons.beta;

  int alphaIndex = 0;
  for (int coreAlpha = 0; coreAlpha < nCoreAlphaOrbitals; ++coreAlpha) {
    while (std::find(activeCasIndices.alpha.begin(), activeCasIndices.alpha.end(), alphaIndex) != activeCasIndices.alpha.end()) {
      ++alphaIndex;
    }
    inactiveCasIndices.alpha.push_back(alphaIndex++);
  }

  int betaIndex = 0;
  for (int coreBeta = 0; coreBeta < nCoreBetaOrbitals; ++coreBeta) {
    while (std::find(activeCasIndices.beta.begin(), activeCasIndices.beta.end(), betaIndex) != activeCasIndices.beta.end()) {
      ++betaIndex;
    }
    inactiveCasIndices.beta.push_back(betaIndex++);
  }

  return {casElectrons, mos_.isRestricted() ? mos_.toUnrestricted() : mos_, activeCasIndices, inactiveCasIndices};
}

CasSpecifier<Reference::Restricted> CasGenerator::generateRestricted(int totalNumberOfElectrons,
                                                                     int orbitalsAroundFermiLevel) const {
  auto homo = findHomo(totalNumberOfElectrons);
  SpinAdaptedContainer<Reference::Restricted, std::vector<int>> inactiveCasIndices;
  SpinAdaptedContainer<Reference::Restricted, std::vector<int>> activeCasIndices;
  SpinAdaptedContainer<Reference::Restricted, int> casElectrons;

  int numberActiveOccupied = orbitalsAroundFermiLevel / 2;
  if (numberActiveOccupied == 0) {
    throw std::runtime_error("No electrons or orbitals in CAS!");
  }

  int nInactiveOrbitals = totalNumberOfElectrons / 2 - numberActiveOccupied;
  std::vector<int> coreOrbitalsIndices(nInactiveOrbitals);
  std::iota(coreOrbitalsIndices.begin(), coreOrbitalsIndices.end(), 0);
  inactiveCasIndices.restricted = std::move(coreOrbitalsIndices);

  activeCasIndices.restricted.reserve(orbitalsAroundFermiLevel);
  for (int casOrb = 0; casOrb < orbitalsAroundFermiLevel; ++casOrb) {
    int index = casOrb + homo - numberActiveOccupied + 1;
    activeCasIndices.restricted.push_back(index);
  }
  casElectrons.restricted = numberActiveOccupied * 2;

  return {casElectrons, mos_, activeCasIndices, inactiveCasIndices};
}

CasSpecifier<Reference::Unrestricted> CasGenerator::generateUnrestricted(int totalNumberElectrons, int spinMultiplicity,
                                                                         int orbitalsAroundFermiLevel) const {
  int alphaElectrons = (spinMultiplicity + totalNumberElectrons - 1) / 2;
  int betaElectrons = (totalNumberElectrons - spinMultiplicity + 1) / 2;
  auto alphaHomo = findHomo(alphaElectrons, SpinComponent::Alpha);
  auto betaHomo = findHomo(betaElectrons, SpinComponent::Beta);

  int nInactiveAlphaOrbitals = alphaElectrons - orbitalsAroundFermiLevel / 2;
  int nInactiveBetaOrbitals = betaElectrons - orbitalsAroundFermiLevel / 2;
  std::vector<int> coreAlphaOrbitalsIndices(nInactiveAlphaOrbitals);
  std::iota(coreAlphaOrbitalsIndices.begin(), coreAlphaOrbitalsIndices.end(), 0);
  std::vector<int> coreBetaOrbitalsIndices(nInactiveBetaOrbitals);
  std::iota(coreBetaOrbitalsIndices.begin(), coreBetaOrbitalsIndices.end(), 0);

  SpinAdaptedContainer<Reference::Unrestricted, int> casElectrons{alphaElectrons - nInactiveAlphaOrbitals,
                                                                  betaElectrons - nInactiveBetaOrbitals};
  SpinAdaptedContainer<Reference::Unrestricted, std::vector<int>> inactiveCasIndices{std::move(coreAlphaOrbitalsIndices),
                                                                                     std::move(coreBetaOrbitalsIndices)};
  SpinAdaptedContainer<Reference::Unrestricted, std::vector<int>> activeCasIndices;

  if (casElectrons.alpha == 0 && casElectrons.beta == 0) {
    throw std::runtime_error("No electrons in CAS!");
  }
  if (casElectrons.alpha < 0 || casElectrons.beta < 0) {
    throw std::runtime_error("Not enough electrons to form a CAS!");
  }
  activeCasIndices.alpha.reserve(orbitalsAroundFermiLevel);
  activeCasIndices.beta.reserve(orbitalsAroundFermiLevel);
  for (int casOrb = 0; casOrb < orbitalsAroundFermiLevel; ++casOrb) {
    int alphaIndex = casOrb + alphaHomo - casElectrons.alpha + 1;
    int betaIndex = casOrb + betaHomo - casElectrons.beta + 1;
    activeCasIndices.alpha.push_back(alphaIndex);
    activeCasIndices.beta.push_back(betaIndex);
  }

  return {casElectrons, mos_.isRestricted() ? mos_.toUnrestricted() : mos_, activeCasIndices, inactiveCasIndices};
}

int CasGenerator::findHomo(int nElectrons) {
  LcaoUtils::ElectronicOccupation occupation;
  occupation.fillLowestRestrictedOrbitalsWithElectrons(nElectrons);
  const auto& restrictedOrbitals = occupation.getFilledRestrictedOrbitals();
  return restrictedOrbitals[restrictedOrbitals.size() - 1];
}

int CasGenerator::findHomo(int nElectrons, SpinComponent type) {
  LcaoUtils::ElectronicOccupation occupation;
  if (type == SpinComponent::Alpha) {
    occupation.fillLowestUnrestrictedOrbitals(nElectrons, 0);
    const auto& alphaOrbitals = occupation.getFilledAlphaOrbitals();
    return alphaOrbitals[alphaOrbitals.size() - 1];
  }
  occupation.fillLowestUnrestrictedOrbitals(0, nElectrons);
  const auto& betaOrbitals = occupation.getFilledBetaOrbitals();
  return betaOrbitals[betaOrbitals.size() - 1];
}
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
