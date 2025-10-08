/**
 * @file ThermochemistryContainer.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RRKMRateConstantCalculator.h"
#include <Utils/Constants.h>
#include <Utils/Properties/Thermochemistry/MolecularDegreesOfFreedom.h>
#include <Utils/Properties/Thermochemistry/RigidRotor.h>
#include <Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h>

namespace Scine::Utils {

std::pair<Eigen::VectorXd, Eigen::VectorXd> RRKMRateConstantCalculator::getRateConstants(const Eigen::VectorXd& energies) {
  /*
   * The rate expression is simply:
   *
   * k(E) = N^\ddagger(E - V^\ddagger) / (h \Omega(E))
   *
   * Note that I assume that the input energy vector are total energies. In the equation above,
   * E is the energy relative to the reactants. Therefore, we subtract the reactant energy from
   * the input vector. Furthermore, the density of states provided by the DegreeOfFreedom objects
   * is \rho(E) = \Omega(E) * \delta E, where \delta E is the energy increment on the energy grid.
   * Therefore, we have to correct for this factor.
   */
  constexpr double h = Constants::planckConstant * Constants::hartree_per_joule;
  const double eZeroTS = tsDegreesOfFreedom_->helmholtzFreeEnergy(ReferenceStates::VacuumZeroKelvin);
  const double eZeroR = reactantDegreesOfFreedom_->helmholtzFreeEnergy(ReferenceStates::VacuumZeroKelvin);
  double deltaE = (energies.size() == 1) ? 1.0 : energies[1] - energies[0];
  const Eigen::VectorXd shiftedEnergiesR = energies.array() - eZeroR;
  const Eigen::VectorXd shiftedEnergiesTS = energies.array() - eZeroTS;
  const Eigen::VectorXd sumOfStatesTS =
      tsDegreesOfFreedom_->getRRKMDegreesOfFreedom(activeRotor_)->sumOfStates(shiftedEnergiesTS);
  const Eigen::VectorXd densityOfStatesR =
      reactantDegreesOfFreedom_->getRRKMDegreesOfFreedom(activeRotor_)->densityOfState(shiftedEnergiesR);
  const Eigen::VectorXd forwardRateConstants = sumOfStatesTS.array() * deltaE / (h * densityOfStatesR.array());
  Eigen::VectorXd backwardRateConstants = Eigen::VectorXd::Zero(energies.size());
  if (productDegreesOfFreedom_) {
    const double eZeroP = productDegreesOfFreedom_->helmholtzFreeEnergy(ReferenceStates::VacuumZeroKelvin);
    const Eigen::VectorXd shiftedEnergiesP = energies.array() - eZeroP;
    const Eigen::VectorXd densityOfStatesP =
        productDegreesOfFreedom_->getRRKMDegreesOfFreedom(activeRotor_)->densityOfState(shiftedEnergiesP);
    backwardRateConstants = sumOfStatesTS.array() * deltaE / (h * densityOfStatesP.array());
  }
  return {forwardRateConstants, backwardRateConstants};
}
RRKMRateConstantCalculator::RRKMRateConstantCalculator(std::shared_ptr<MolecularDegreesOfFreedom> reactantDegreesOfFreedom,
                                                       std::shared_ptr<MolecularDegreesOfFreedom> tsDegreesOfFreedom,
                                                       std::shared_ptr<MolecularDegreesOfFreedom> productDegreesOfFreedom,
                                                       bool activeRotor)
  : tsDegreesOfFreedom_(tsDegreesOfFreedom),
    reactantDegreesOfFreedom_(reactantDegreesOfFreedom),
    productDegreesOfFreedom_(productDegreesOfFreedom),
    activeRotor_(activeRotor) {
}
Eigen::VectorXd RRKMRateConstantCalculator::getEquilibriumConstants(const Eigen::VectorXd& energies) {
  if (!productDegreesOfFreedom_) {
    throw std::runtime_error("No product degrees of freedom available to calculate the equilibrium constant");
  }
  // Note that the equilibrium constant will simply be the ratios of the density of states of reactants and products.
  const double eZeroR = reactantDegreesOfFreedom_->helmholtzFreeEnergy(ReferenceStates::VacuumZeroKelvin);
  const double eZeroP = productDegreesOfFreedom_->helmholtzFreeEnergy(ReferenceStates::VacuumZeroKelvin);
  return productDegreesOfFreedom_->getRRKMDegreesOfFreedom(activeRotor_)->densityOfState(energies.array() - eZeroP).array() /
         reactantDegreesOfFreedom_->getRRKMDegreesOfFreedom(activeRotor_)->densityOfState(energies.array() - eZeroR).array();
}

} // namespace Scine::Utils
