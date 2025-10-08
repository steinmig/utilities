/**
 * @file HarmonicOscillator.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "HarmonicOscillator.h"
#include <Utils/Constants.h>
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <boost/math/special_functions/factorials.hpp>

namespace Scine {
namespace Utils {

HarmonicOscillator::~HarmonicOscillator() = default;

HarmonicOscillator::HarmonicOscillator(std::shared_ptr<NormalModesContainer> normalModesContainer, bool classicalMechanics)
  : normalModesContainer_(normalModesContainer), classicalMechanics_(classicalMechanics) {
}
HarmonicOscillator::HarmonicOscillator(const HessianMatrix& hessian, const AtomCollection& atoms, bool classicalMechanics)
  : classicalMechanics_(classicalMechanics) {
  normalModesContainer_ = std::make_shared<NormalModesContainer>(
      Utils::NormalModeAnalysis::calculateNormalModes(hessian, atoms.getElements(), atoms.getPositions()));
}
void HarmonicOscillator::analyzeVibrations() {
  waveNumberProduct_ = std::make_unique<double>(1.0);
  zeroPointVibrationalCorrection_ = std::make_unique<double>(0.0);
  nRealModes_ = std::make_unique<unsigned int>(0);
  for (auto wn : this->normalModesContainer_->getWaveNumbers()) {
    if (wn > 0) {
      double hNu = wn * Constants::hartree_per_invCentimeter;
      *zeroPointVibrationalCorrection_ += hNu;
      *waveNumberProduct_ *= hNu;
      ++*nRealModes_;
    }
  }
  *zeroPointVibrationalCorrection_ *= 0.5;
}
Eigen::VectorXd HarmonicOscillator::densityOfState(const Eigen::VectorXd& energies) {
  if (!waveNumberProduct_) {
    this->analyzeVibrations();
  }
  Eigen::VectorXd densityOfStates = Eigen::VectorXd::Zero(energies.size());
  if (energies.size() == 0) {
    return densityOfStates;
  }
  if (*nRealModes_ < 1) {
    densityOfStates[0] = 1.0;
    return densityOfStates;
  }
  assertNonNegativeEnergies(energies);
  if (classicalMechanics_) {
    if (*nRealModes_ > boost::math::max_factorial<double>::value) {
      throw std::runtime_error("Factorial is out of bounds. Too many vibrations to enumerate the states");
    }
    long double denominator = boost::math::factorial<double>(*nRealModes_ - 1) * *waveNumberProduct_;
    double deltaE = (energies.size() == 1) ? 1.0 : energies[1] - energies[0];
    // Note that the shift by the ZPE is a correction of the classical density of states
    // (see https://doi.org/10.1016/S0069-8040(08)70206-1 Fig. 3)
    densityOfStates = (energies.array() + *zeroPointVibrationalCorrection_).pow(*nRealModes_ - 1) / denominator;
    densityOfStates *= deltaE;
  }
  else {
    /*
     * The individual oscillators are independent degrees of freedom. Use convolution to calculate the joint
     * sum of states.
     */
    densityOfStates[0] = 1.0;
    const auto waveNumbers = normalModesContainer_->getWaveNumbers();
    for (const auto& wn : waveNumbers) {
      if (wn > 0.0) {
        double energyGap = wn * Constants::hartree_per_invCentimeter;
        densityOfStates = this->bsEquidistantStates(energies, densityOfStates, energyGap, 1.0);
      }
    }
  }
  return densityOfStates;
}
Eigen::VectorXd HarmonicOscillator::sumOfStates(const Eigen::VectorXd& energies) {
  if (!waveNumberProduct_) {
    this->analyzeVibrations();
  }
  if (*nRealModes_ < 1) {
    return Eigen::VectorXd::Constant(energies.size(), 1.0);
  }
  if (*nRealModes_ > boost::math::max_factorial<double>::value) {
    throw std::runtime_error("Factorial is out of bounds. Too many vibrations to enumerate the states");
  }
  assertNonNegativeEnergies(energies);
  Eigen::VectorXd sumOverStates = Eigen::VectorXd::Zero(energies.size());
  if (energies.size() == 0) {
    return sumOverStates;
  }
  if (classicalMechanics_) {
    double denominator = boost::math::factorial<double>(*nRealModes_) * *waveNumberProduct_;
    // Note that the shift by the ZPE is a correction of the classical density of states
    // (see https://doi.org/10.1016/S0069-8040(08)70206-1 Fig. 3)
    sumOverStates.array() = (energies.array() + *zeroPointVibrationalCorrection_).pow(*nRealModes_) / denominator;
  }
  else {
    /*
     * The individual oscillators are independent degrees of freedom. Use convolution to calculate the joint
     * density of states.
     */
    sumOverStates = Eigen::VectorXd::Ones(energies.size());
    const auto waveNumbers = normalModesContainer_->getWaveNumbers();
    for (const auto& wn : waveNumbers) {
      if (wn > 0.0) {
        double energyGap = wn * Constants::hartree_per_invCentimeter;
        sumOverStates = this->bsEquidistantStates(energies, sumOverStates, energyGap, 1.0);
      }
    }
  }
  return sumOverStates;
}
double HarmonicOscillator::partitionFunction(const ThermodynamicReferenceState& state) {
  double q = 1.0;
  if (state.temperature < 1e-6) {
    return 1.0;
  }
  double beta = 1.0 / (kB * state.temperature);
  for (auto wn : this->normalModesContainer_->getWaveNumbers()) {
    if (wn > 0) {
      double hNuBeta = wn * Constants::hartree_per_invCentimeter * beta;
      // If the zero point hv/2 is included, we have this additional term.
      // However, we simply choose the zero point to be hv/2.
      //      double zeroPoint = std::exp(- hNuBeta / 2.0);
      double level = 1 - std::exp(-hNuBeta);
      //      q *= zeroPoint / level;
      q *= 1.0 / level;
    }
  }
  return q;
}
double HarmonicOscillator::entropy(const ThermodynamicReferenceState& state) {
  double s = 0.0;
  if (state.temperature < 1e-6) {
    return s;
  }
  double beta = 1.0 / (kB * state.temperature);
  for (auto wn : this->normalModesContainer_->getWaveNumbers()) {
    if (wn > 0) {
      double hNu = wn * Constants::hartree_per_invCentimeter;
      double expBetaNu = std::exp(-hNu * beta);
      s += -kB * std::log(1 - expBetaNu) + hNu / state.temperature * expBetaNu / (1 - expBetaNu);
    }
  }
  return s;
}
double HarmonicOscillator::enthalpy(const ThermodynamicReferenceState& state) {
  if (!this->nRealModes_) {
    analyzeVibrations();
  }
  double h = *zeroPointVibrationalCorrection_;
  if (state.temperature < 1e-6) {
    return h;
  }
  double beta = 1.0 / (kB * state.temperature);
  for (auto wn : this->normalModesContainer_->getWaveNumbers()) {
    if (wn > 0) {
      double hNu = wn * Constants::hartree_per_invCentimeter;
      h += hNu / (std::exp(beta * hNu) - 1);
    }
  }
  return h;
}
double HarmonicOscillator::heatCapacityCp(const ThermodynamicReferenceState& state) {
  double cp = 0.0;
  if (state.temperature < 1e-6) {
    return cp;
  }
  double beta = 1.0 / (kB * state.temperature);
  for (auto wn : this->normalModesContainer_->getWaveNumbers()) {
    if (wn > 0) {
      double hNu = wn * Constants::hartree_per_invCentimeter;
      double hNuBeta = hNu * beta;
      double exponential = std::exp(hNuBeta);
      double oneMinusExp = 1 - exponential;
      cp += hNuBeta * hNuBeta * exponential / (oneMinusExp * oneMinusExp);
    }
  }
  return cp * kB;
}
double HarmonicOscillator::heatCapacityCv(const ThermodynamicReferenceState& state) {
  return this->heatCapacityCp(state) * 3.0 / 5.0;
}
std::shared_ptr<HarmonicOscillator> HarmonicOscillator::getClassicalVariant() {
  if (classicalMechanics_) {
    return this->shared_from_this();
  }
  if (!classicalVariant_) {
    classicalVariant_ = std::make_shared<HarmonicOscillator>(normalModesContainer_, true);
  }
  return classicalVariant_;
}
const NormalModesContainer& HarmonicOscillator::getNormalModesContainer() {
  return *this->normalModesContainer_;
}

} // namespace Utils
} // namespace Scine