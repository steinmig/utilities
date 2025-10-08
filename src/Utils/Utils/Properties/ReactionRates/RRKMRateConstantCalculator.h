/**
 * @file ThermochemistryContainer.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_RRKM_H
#define UTILSOS_RRKM_H

#include <Eigen/Dense>
#include <memory>

namespace Scine {
namespace Utils {
class DegreesOfFreedomCollection;
class MolecularDegreesOfFreedom;
/**
 * @class
 * @brief Calculates rate constants according to the Rice–Ramsperger–Kassel–Marcus (RRKM) theory.
 *
 * RRKM rate constants are unimolecular microcanonical rate constants, i.e., they provide a rate constant
 * for a given energy of the microstate.
 * The rough idea of the model is that the reactant and transition state are modeled as harmonic oscillators
 * on which we can partition the energy E in the system. For the reaction to occur, more energy
 * than \f$ V^\ddagger \f$ must accumulate in the mode leading to the transition state. The reaction rate constant is
 * then given as the probability that this happens assuming that every oscillation mode is equally likely to be
 * excited.
 *
 * Assumptions:
 * * Harmonic oscillators for vibrational degrees of freedom.
 * * No recrossing over the transitions state.
 * * Energy exchange between modes is much faster than the reaction.
 * * Classical dynamics for the reaction (no tunneling, calssical harmonic oscillators).
 *
 * Because of the "no-recrossing" assumption, reaction rate constants are typically overestimated if tunneling is
 * not important. Eyring's transition state theory for the canonical ensemble is recovered by integration:
 *
 * \f$
 * k(\beta) = \int_0^\infty k(E) P(E, \beta) \mathrm{d}E~,
 * \f$
 * where \f$ P(E, \beta) \f$ is the probability of finding a microstate with energy \f$ E \f$
 * (typically the Boltzman distribution) and \f$ \beta = 1/(k_B T) \f$.
 *
 * An introduction to RRKM theory is given in https://doi.org/10.1016/S0069-8040(08)70206-1
 */
class RRKMRateConstantCalculator {
 public:
  /**
   * @brief Constructor.
   * @param reactantDegreesOfFreedom  The reactant's degrees of freedom.
   * @param tsDegreesOfFreedom        The transition state's degrees of freedom..
   * @param productDegreesOfFreedom   The product's degrees of freedom.
   * @param activeRotor               If true, rotational degrees of freedom are considered
   *                                  active during RRKM evaluation, i.e., energy may be dispersed to them.
   */
  explicit RRKMRateConstantCalculator(std::shared_ptr<MolecularDegreesOfFreedom> reactantDegreesOfFreedom,
                                      std::shared_ptr<MolecularDegreesOfFreedom> tsDegreesOfFreedom,
                                      std::shared_ptr<MolecularDegreesOfFreedom> productDegreesOfFreedom = nullptr,
                                      bool activeRotor = false);
  /**
   * @brief Calculate the rate constant for the given energy.
   * @param energies               The energies for which to calculate the rate constants.
   * @return The rate constants
   */
  std::pair<Eigen::VectorXd, Eigen::VectorXd> getRateConstants(const Eigen::VectorXd& energies);
  /**
   * @brief Calculate the equilibrium constants from the density of states.
   * @param energies The energy range.
   * @return The equilibrium constants.
   */
  Eigen::VectorXd getEquilibriumConstants(const Eigen::VectorXd& energies);

 private:
  std::shared_ptr<MolecularDegreesOfFreedom> tsDegreesOfFreedom_;
  std::shared_ptr<MolecularDegreesOfFreedom> reactantDegreesOfFreedom_;
  std::shared_ptr<MolecularDegreesOfFreedom> productDegreesOfFreedom_;
  const bool activeRotor_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_RRKM_H
