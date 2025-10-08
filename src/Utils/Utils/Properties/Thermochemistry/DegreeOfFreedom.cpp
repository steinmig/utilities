/**
 * @file DegreeOfFreedom.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DegreeOfFreedom.h"
#include <iostream>

namespace Scine {
namespace Utils {

DegreeOfFreedom::~DegreeOfFreedom() = default;

double DegreeOfFreedom::helmholtzFreeEnergy(const ThermodynamicReferenceState& state) {
  return this->enthalpy(state) - state.temperature * this->entropy(state);
}
unsigned int DegreeOfFreedom::degeneracy(unsigned int) {
  throw std::runtime_error("Getter for the degeneracy is not implemented by default.");
}
double DegreeOfFreedom::energyLevel(unsigned int) {
  throw std::runtime_error("Getter for the energy level is not implemented by default.");
}
Eigen::VectorXd DegreeOfFreedom::calculateQuantumDensityOfStates(const Eigen::VectorXd& energies, unsigned int minQuantumNumber,
                                                                 unsigned int maxQuantumNumber) {
  if (energies.size() == 0) {
    return Eigen::VectorXd::Zero(0);
  }
  assertNonNegativeEnergies(energies);
  Eigen::VectorXd t = Eigen::VectorXd::Zero(energies.size());
  t[0] = 1.0;
  return bsseAlgorithm(energies, t, minQuantumNumber, maxQuantumNumber);
}
Eigen::VectorXd DegreeOfFreedom::calculateQuantumSumOfStates(const Eigen::VectorXd& energies,
                                                             unsigned int minQuantumNumber, unsigned int maxQuantumNumber) {
  if (energies.size() == 0) {
    return Eigen::VectorXd::Zero(0);
  }
  assertNonNegativeEnergies(energies);
  Eigen::VectorXd t = Eigen::VectorXd::Constant(energies.size(), 1.0);
  for (unsigned int i = 0; i < energies.size(); ++i) {
    if (energies[i] < 0.0) {
      t[i] = 0.0;
    }
  }
  return bsseAlgorithm(energies, t, minQuantumNumber, maxQuantumNumber);
}
Eigen::VectorXd DegreeOfFreedom::bsseAlgorithm(const Eigen::VectorXd& energies, const Eigen::VectorXd& t,
                                               unsigned int minQuantumNumber, unsigned int maxQuantumNumber) {
  /*
   * Beyer-Swinehart-Stein-Rabinovitch enumeration
   * https://pubs.aip.org/aip/jcp/article/58/6/2438/84783/Accurate-evaluation-of-internal-energy-level-sums
   * For the density of states:
   * t(0) = 1.0
   * Search for energy gaps that match/exceed the quantum energy level. If found, add the number of states to the
   * density of states for the energy grain.
   *
   * For the sum of states:
   * t(i) = 1.0 for all i
   * Effectively sum the density of states for all energy grains.
   */
  Eigen::VectorXd result = Eigen::VectorXd::Zero(energies.size());
  unsigned int n = minQuantumNumber;
  double eN = energyLevel(n);
  double gN = degeneracy(n);
  double eMax = energies.maxCoeff();
  while (eN < eMax) {
    for (unsigned int i = 0; i < energies.size(); ++i) {
      for (unsigned int j = i; j < energies.size(); ++j) {
        if (std::abs(energies[i] - energies[j]) >= eN * 0.999999) {
          result[j] += gN * t[i];
          break;
        }
      }
    }
    ++n;
    if (n > maxQuantumNumber) {
      break;
    }
    eN = energyLevel(n);
    gN = degeneracy(n);
  }
  return result;
}
double DegreeOfFreedom::calculateQuantumPartitioningFunction(double temperature, unsigned int minQuantumNumber,
                                                             unsigned int maxQuantumNumber, double precision) {
  double q = 0.0;
  unsigned int n = minQuantumNumber;
  const double beta = 1.0 / (kB * temperature);
  while (true) {
    double eN = energyLevel(n);
    unsigned int gN = degeneracy(n);
    double delta = gN * std::exp(-beta * eN);
    q += delta;

    ++n;
    if (delta < q * precision || n > maxQuantumNumber) {
      break;
    }
  }
  return q;
}
double DegreeOfFreedom::calculateQuantumEnthalpy(double temperature, unsigned int minQuantumNumber,
                                                 unsigned int maxQuantumNumber, double precision) {
  double h = std::get<1>(calculatePartitionFunctionAndEnthalpy(temperature, minQuantumNumber, maxQuantumNumber, precision));
  return h;
}
double DegreeOfFreedom::calculateQuantumEntropy(double temperature, unsigned int minQuantumNumber,
                                                unsigned int maxQuantumNumber, double precision) {
  double q, h;
  std::tie(q, h) = calculatePartitionFunctionAndEnthalpy(temperature, minQuantumNumber, maxQuantumNumber, precision);
  if (temperature < 1e-6) {
    return kB * log(q);
  }
  double a = -kB * temperature * std::log(q);
  return (h - a) / temperature;
}
std::tuple<double, double>
DegreeOfFreedom::calculatePartitionFunctionAndEnthalpy(double temperature, unsigned int minQuantumNumber,
                                                       unsigned int maxQuantumNumber, double precision) {
  unsigned int n = minQuantumNumber;
  if (temperature < 1e-6) {
    return {degeneracy(n), 0.0};
  }
  double q = 0.0;
  double weightedEnergy = 0.0;
  const double beta = 1.0 / (kB * temperature);
  while (true) {
    double eN = energyLevel(n);
    unsigned int gN = degeneracy(n);
    double delta = gN * std::exp(-beta * eN);
    double deltaE = eN * delta;
    weightedEnergy += deltaE;
    q += delta;
    ++n;
    if ((delta < q * precision && deltaE < weightedEnergy * precision) || n > maxQuantumNumber) {
      break;
    }
  }
  double h = weightedEnergy / q;
  return {q, h};
}
Eigen::VectorXd DegreeOfFreedom::bsEquidistantStates(const Eigen::VectorXd& energies, const Eigen::VectorXd& t,
                                                     double energyGap, unsigned int degeneracy) {
  Eigen::VectorXd result = t;
  for (unsigned int i = 0; i < energies.size(); ++i) {
    for (unsigned int j = i; j < energies.size(); ++j) {
      if (std::fabs(energies[i] - energies[j]) > energyGap * 0.999999) {
        result[j] += degeneracy * result[i];
        break;
      }
    }
  }
  return result;
}
ThermochemicalContainer DegreeOfFreedom::createThermochemicalContainer(const ThermodynamicReferenceState& state,
                                                                       bool excludeZeroPoint) {
  ThermochemicalContainer container{};
  container.partitionFunction = partitionFunction(state);
  container.gibbsFreeEnergy = helmholtzFreeEnergy(state); // The contribution P * V is ignored for the time being.
  container.heatCapacityV = heatCapacityCv(state);
  container.heatCapacityP = heatCapacityCp(state);
  container.entropy = entropy(state);
  container.enthalpy = enthalpy(state);
  if (excludeZeroPoint) {
    auto zpe = this->enthalpy(ReferenceStates::VacuumZeroKelvin);
    container.enthalpy -= zpe;
    container.gibbsFreeEnergy -= zpe;
  }
  return container;
}
void DegreeOfFreedom::assertNonNegativeEnergies(const Eigen::VectorXd& energies) {
  if ((energies.array() < 0.0).count()) {
    throw std::runtime_error("Negative energy encountered for density/sum of states calculation.\n"
                             "It is impossible to distribute a negative energy on the degrees of freedom.");
  }
}

} // namespace Utils
} // namespace Scine