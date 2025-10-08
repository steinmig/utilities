/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CASSPECIFIER_H
#define UTILS_CASSPECIFIER_H

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SecondQuantization/ExcitationHelpers.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <numeric>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

namespace detail {
inline int total(const SpinAdaptedContainer<Reference::Restricted, int> c) {
  return c.restricted;
}
inline int total(const SpinAdaptedContainer<Reference::Unrestricted, int> c) {
  return c.alpha + c.beta;
}
} // namespace detail

/**
 * @brief Specifier of a Complete Active Space. Contains number of electrons, expansion coefficients and energies.
 * This specifier is enough to perform 4-index transformation of the integrals.
 */
template<Reference res>
class CasSpecifier {
 public:
  /**
   * @brief Constructor from the MOs in the CAS and the MO energies.
   */
  CasSpecifier(SpinAdaptedContainer<res, int> casElectrons, MolecularOrbitals coefficientMatrix,
               SpinAdaptedContainer<res, std::vector<int>> indexActiveOrbitals,
               SpinAdaptedContainer<res, std::vector<int>> indexInactiveOrbitals)
    : casElectrons_(std::move(casElectrons)),
      coefficientMatrix_(std::move(coefficientMatrix)),
      activeIndices_(std::move(indexActiveOrbitals)),
      inactiveIndices_(std::move(indexInactiveOrbitals)),
      totalNumberOfOrbitals_(coefficientMatrix_.numberOrbitals()) {
  }
  CasSpecifier() = default;
  CasSpecifier(const CasSpecifier<res>& rhs)
    : casElectrons_(rhs.casElectrons_),
      coefficientMatrix_(rhs.coefficientMatrix_),
      activeIndices_(rhs.activeIndices_),
      inactiveIndices_(rhs.inactiveIndices_),
      totalNumberOfOrbitals_(rhs.totalNumberOfOrbitals_) {
  }
  CasSpecifier(CasSpecifier<res>&& rhs) noexcept = default;
  CasSpecifier<res>& operator=(const CasSpecifier<res>& rhs) {
    casElectrons_ = rhs.casElectrons_;
    coefficientMatrix_ = rhs.coefficientMatrix_;
    activeIndices_ = rhs.activeIndices_;
    inactiveIndices_ = rhs.inactiveIndices_;
    totalNumberOfOrbitals_ = rhs.totalNumberOfOrbitals_;
    return *this;
  }
  CasSpecifier<res>& operator=(CasSpecifier<res>&& rhs) noexcept = default;

  bool empty() const {
    return detail::total(casElectrons_) == 0;
  }

  /**
   * @brief Returns the number of orbitals in the CAS.
   */
  const SpinAdaptedContainer<res, int>& getCasElectrons() const {
    return casElectrons_;
  }

  /**
   * @brief Returns the coefficient matrix of all the orbitals (not only CAS).
   */
  const MolecularOrbitals& getCasOrbitals() const {
    return coefficientMatrix_;
  }

  /**
   * @brief Returns the coefficient matrix of the active orbitals.
   */
  MolecularOrbitals getActiveOrbitalsCoefficients() const;
  /**
   * @brief Returns the coefficient matrix of the core orbitals.
   */
  MolecularOrbitals getCoreOrbitalsCoefficients() const;
  /**
   * @brief Maps the CAS indexing to the original indexing.
   *
   * Example:
   * orbitals:
   * 0 1 2 3 4 5 6 7 8
   *     |-------|
   *      active
   *
   * CAS index mappings:
   * CAS     Original
   * 0   --> 2
   * 1   --> 3
   * 2   --> 4
   * 3   --> 5
   * 4   --> 6
   *
   */
  int getOriginalIndex(int casIndex, SecondQuantization::SpinComponent spin = SecondQuantization::SpinComponent::Alpha) const;

  /**
   * @brief Returns the indices (original numbering) of the CAS orbitals.
   */
  const SpinAdaptedContainer<res, std::vector<int>>& getActiveIndices() const {
    return activeIndices_;
  };

  /**
   * @brief Returns the indices (original numbering) of the CAS orbitals.
   */
  const SpinAdaptedContainer<res, std::vector<int>>& getInactiveIndices() const {
    return inactiveIndices_;
  };

  /**
   * @brief Returns the total number (non only CAS) of orbitals in the system.
   */
  int getTotalNumberOfOrbitals() const {
    return totalNumberOfOrbitals_;
  }

 private:
  SpinAdaptedContainer<res, int> casElectrons_{};
  MolecularOrbitals coefficientMatrix_;
  SpinAdaptedContainer<res, std::vector<int>> activeIndices_;
  SpinAdaptedContainer<res, std::vector<int>> inactiveIndices_;
  int totalNumberOfOrbitals_{};
};

template<>
inline MolecularOrbitals CasSpecifier<Reference::Restricted>::getActiveOrbitalsCoefficients() const {
  Eigen::MatrixXd activeCoefficients(coefficientMatrix_.restrictedMatrix().rows(), activeIndices_.restricted.size());
  int index = 0;
  for (int activeOrbital : activeIndices_.restricted) {
    activeCoefficients.col(index++) = coefficientMatrix_.restrictedMatrix().col(activeOrbital);
  }
  return MolecularOrbitals::createFromRestrictedCoefficients(activeCoefficients);
}
template<>
inline MolecularOrbitals CasSpecifier<Reference::Unrestricted>::getActiveOrbitalsCoefficients() const {
  Eigen::MatrixXd activeAlphaCoefficients(coefficientMatrix_.alphaMatrix().rows(), activeIndices_.alpha.size());
  Eigen::MatrixXd activeBetaCoefficients(coefficientMatrix_.betaMatrix().rows(), activeIndices_.beta.size());
  int index = 0;
  for (int activeAlphaOrbital : activeIndices_.alpha) {
    activeAlphaCoefficients.col(index++) = coefficientMatrix_.alphaMatrix().col(activeAlphaOrbital);
  }
  index = 0;
  for (int activeBetaOrbital : activeIndices_.beta) {
    activeBetaCoefficients.col(index++) = coefficientMatrix_.betaMatrix().col(activeBetaOrbital);
  }
  return MolecularOrbitals::createFromUnrestrictedCoefficients(activeAlphaCoefficients, activeBetaCoefficients);
}

template<>
inline MolecularOrbitals CasSpecifier<Reference::Restricted>::getCoreOrbitalsCoefficients() const {
  Eigen::MatrixXd coreCoefficients(coefficientMatrix_.restrictedMatrix().rows(), inactiveIndices_.restricted.size());
  int index = 0;
  for (int coreOrbital : inactiveIndices_.restricted) {
    coreCoefficients.col(index++) = coefficientMatrix_.restrictedMatrix().col(coreOrbital);
  }
  return MolecularOrbitals::createFromRestrictedCoefficients(coreCoefficients);
}
template<>
inline MolecularOrbitals CasSpecifier<Reference::Unrestricted>::getCoreOrbitalsCoefficients() const {
  Eigen::MatrixXd coreAlphaCoefficients(coefficientMatrix_.alphaMatrix().rows(), inactiveIndices_.alpha.size());
  Eigen::MatrixXd coreBetaCoefficients(coefficientMatrix_.betaMatrix().rows(), inactiveIndices_.beta.size());
  int index = 0;
  for (int coreAlphaOrbital : inactiveIndices_.alpha) {
    coreAlphaCoefficients.col(index++) = coefficientMatrix_.alphaMatrix().col(coreAlphaOrbital);
  }
  index = 0;
  for (int coreBetaOrbital : inactiveIndices_.beta) {
    coreBetaCoefficients.col(index++) = coefficientMatrix_.betaMatrix().col(coreBetaOrbital);
  }
  return MolecularOrbitals::createFromUnrestrictedCoefficients(coreAlphaCoefficients, coreBetaCoefficients);
}

template<>
inline int CasSpecifier<Reference::Restricted>::getOriginalIndex(int casIndex, SecondQuantization::SpinComponent /*spin*/) const {
  if (casIndex > int(activeIndices_.restricted.size())) {
    throw std::runtime_error("casIndex out of scope!");
  }
  return activeIndices_.restricted[casIndex];
}

template<>
inline int CasSpecifier<Reference::Unrestricted>::getOriginalIndex(int casIndex, SecondQuantization::SpinComponent spin) const {
  if (spin == SecondQuantization::SpinComponent::Alpha) {
    if (casIndex > int(activeIndices_.alpha.size())) {
      throw std::runtime_error("Alpha casIndex out of scope!");
    }
    return activeIndices_.alpha[casIndex];
  }
  if (casIndex > int(activeIndices_.beta.size())) {
    throw std::runtime_error("Beta casIndex out of scope!");
  }
  return activeIndices_.beta[casIndex];
}

/**
 * @brief Generates a specifier of a Complete Active Space.
 * This specifier is enough to perform 4-index transformation of the integrals. Also contains the original
 * indices for consistency. It is templated with the restrictedness of the space (spatial or spin orbitals).
 * Right now only restricted is available.
 * @tparam r The restrictedness of the molecular orbitals.
 */
class CasGenerator {
 public:
  /**
   * @brief Constructor taking a molecular orbital as arguments.
   * As soon as the reference arguments are no more valid, this class is also no more valid.
   * @param mos A MolecularOrbitals with the expansion coefficients.
   * @param moEnergies The energies of the orbitals.
   */
  CasGenerator(const MolecularOrbitals& mos) : mos_(mos){};

  /**
   * @brief Generates a CAS specifier with desired number of electrons and MO indices.
   * Assumes an occupation according to the Aufbau principle.
   * @param totalNumberOfElectrons The total number of electrons in the molecular system.
   * @param activeIndices The indices of the (spatial) active molecular orbitals.
   */
  CasSpecifier<Reference::Restricted> generateRestricted(int totalNumberOfElectrons, std::vector<int> indices) const;

  /**
   * @brief Generates a CAS specifier with desired number of electrons and MO indices.
   * Assumes an occupation according to the Aufbau principle for both alpha and beta.
   * @param totalNumberElectrons The number of electrons in the CAS.
   * @param spinMultiplicity The spin mulitplicity M_s of the total system. 1 for singulett, 2 for dublett,...
   * @param alphaIndices The indices of the alpha (spin) active molecular orbitals.
   * @param betaIndices The indices of the beta (spin) active molecular orbitals.
   * This functions calculates the alpha and beta electrons from the total number with
   *
   * n_\alpha = (n_{tot} + M_s - 1) / 2
   * n_\beta = (n_{tot} - M_s + 1) / 2
   */
  CasSpecifier<Reference::Unrestricted> generateUnrestricted(int totalNumberElectrons, int spinMultiplicity,
                                                             std::vector<int> alphaIndices, std::vector<int> betaIndices) const;

  /**
   * @brief Generates a CAS specifier with desired number of electrons and MO indices.
   * This method assumes an occupation according to the Aufbau principle and
   * spatial orbitals.
   * @param totalNumberOfElectrons The number of electrons in the full space.
   * @param numberOfOrbitalsAroundFermiLevel number of (spatial) molecular orbitals to generate around the Fermi level.
   */
  CasSpecifier<Reference::Restricted> generateRestricted(int totalNumberOfElectrons, int orbitalsAroundFermiLevel) const;
  /**
   * @brief Generates a CAS specifier with desired number of electrons and MO indices.
   * This method assumes an occupation according to the Aufbau principle and
   * spin orbitals.
   * @param totalNumberElectrons The number of electrons in the CAS.
   * @param spinMultiplicity The spin mulitplicity M_s of the total system. 1 for singulett, 2 for dublett,...
   * @param numberOfOrbitalsAroundFermiLevel number of (spin) molecular orbitals to generate around the Fermi level.
   *        It is the same number for alpha and beta. If not enough electrons are present, this throws.
   * @throws if there are not enought alpha/beta electrons
   * This functions calculates the alpha and beta electrons from the total number with
   *
   * n_\alpha = (n_{tot} + M_s - 1) / 2
   * n_\beta = (n_{tot} - M_s + 1) / 2
   */
  CasSpecifier<Reference::Unrestricted> generateUnrestricted(int totalNumberElectrons, int spinMultiplicity,
                                                             int orbitalsAroundFermiLevel) const;

  /**
   * @brief Generates a vector of indices from a comma-separated list.
   */
  static std::vector<int> parseCasString(const std::string& casString) {
    std::vector<int> casIndices;
    std::stringstream ss(casString);
    std::string token;
    while (std::getline(ss, token, ',')) {
      casIndices.push_back(std::stoi(token));
    }
    return casIndices;
  }

  static int findHomo(int nElectrons);
  static int findHomo(int nElectrons, SecondQuantization::SpinComponent type);

 private:
  const MolecularOrbitals& mos_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_CASSPECIFIER_H
