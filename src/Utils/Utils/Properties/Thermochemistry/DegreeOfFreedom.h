/**
 * @file DegreeOfFreedom.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_DEGREEOFFREEDOM_H
#define UTILSOS_DEGREEOFFREEDOM_H

#include "Utils/Properties/Thermochemistry/ThermochemistryCalculator.h"
#include "Utils/Properties/Thermochemistry/ThermodynamicReferenceState.h"
#include <Eigen/Dense>
#include <limits>
#include <memory>

namespace Scine {
namespace Utils {
constexpr double kB = Constants::boltzmannConstant * Constants::hartree_per_joule;
/**
 * @brief This class represent an arbitrary degree of freedom. It provides state enumeration functions
 * to calculate quantum mechanical partition functions, density of states, sum of states, and other
 * thermodynamic properties.
 */
class DegreeOfFreedom {
 public:
  /**
   * @brief Default constructor.
   */
  DegreeOfFreedom() = default;
  /**
   * @brief Default destructor.
   */
  ~DegreeOfFreedom();

  /**
   * @brief Getter for the density of states at the given energies.
   * @param energies The energies.
   * @return The density of states.
   */
  virtual Eigen::VectorXd densityOfState(const Eigen::VectorXd& energies) = 0;
  /**
   * @brief Getter for the sum of states at the given energies.
   * @param energies The energies.
   * @return The sum of states.
   */
  virtual Eigen::VectorXd sumOfStates(const Eigen::VectorXd& energies) = 0;
  /**
   * @brief Getter for the partition function.
   * @param state The thermodynamic reference state.
   * @return The partition function (in atomic units).
   */
  virtual double partitionFunction(const ThermodynamicReferenceState& state) = 0;
  /**
   * @brief Getter for the entropy.
   * @param state The thermodynamic reference state.
   * @return The entropy (in atomic units).
   */
  virtual double entropy(const ThermodynamicReferenceState& state) = 0;
  /**
   * @brief Getter for the enthalpy.
   * @param state The thermodynamic reference state.
   * @return The enthalpy (in atomic units).
   */
  virtual double enthalpy(const ThermodynamicReferenceState& state) = 0;
  /**
   * @brief Getter for the Helmholtz free energy.
   * @param state The thermodynamic reference state.
   * @return The enthalpy (in atomic units).
   */
  virtual double helmholtzFreeEnergy(const ThermodynamicReferenceState& state);
  /**
   * @brief Getter for the heat capacity at constant pressure.
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCp(const ThermodynamicReferenceState& state) = 0;
  /**
   * @brief Getter for the heat capacity at constant volume.
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCv(const ThermodynamicReferenceState& state) = 0;
  /**
   * @brief Getter for the degeneracy of the state with the given quantum number.
   * @param quantumNumber The quantum number.
   * @return The degeneracy.
   */
  virtual unsigned int degeneracy(unsigned int quantumNumber);
  /**
   * @brief Getter for the energy of the state with the given quantum number.
   * @param quantumNumber The quantum number.
   * @return The energy.
   */
  virtual double energyLevel(unsigned int quantumNumber);
  /**
   * @brief Create a thermochemistry container.
   * @param state            The thermodynamic reference state.
   * @param excludeZeroPoint If true, the zero point energy is assumed to be included in the electronic energy.
   * @return The ThermochemicalContainer.
   */
  ThermochemicalContainer createThermochemicalContainer(const ThermodynamicReferenceState& state, bool excludeZeroPoint = false);

 protected:
  /**
   * @brief Use state enumeration to calculate the density of states.
   * @param energies          The energy range.
   * @param minQuantumNumber  The minimum quantum number.
   * @param maxQuantumNumber  The maximum quantum number. If not provided, the states are enumerated until their energy
   * exceeds the maximum energy in the energy array.
   * @return The density of states at the given energies.
   */
  Eigen::VectorXd
  calculateQuantumDensityOfStates(const Eigen::VectorXd& energies, unsigned int minQuantumNumber = 0,
                                  unsigned int maxQuantumNumber = std::numeric_limits<unsigned int>::infinity() - 1);
  /**
   * @brief Use state enumeration to calculate the sum of states.
   * @param energies         The energy range.
   * @param minQuantumNumber The minimum quantum number.
   * @param maxQuantumNumber The maximum quantum number. If not provided, the states are enumerated until their energy
   * exceeds the maximum energy in the energy array.
   * @return The sum of states at the given energies.
   */
  Eigen::VectorXd calculateQuantumSumOfStates(const Eigen::VectorXd& energies, unsigned int minQuantumNumber = 0,
                                              unsigned int maxQuantumNumber = std::numeric_limits<unsigned int>::infinity() - 1);
  /**
   * @brief Calculate the partition function through state enumeration.
   * @param temperature      The temperature.
   * @param minQuantumNumber The minimum quantum number.
   * @param maxQuantumNumber The maximum quantum number. If not provided, the states are enumerated until the partition
   * function converges.
   * @param precision        The precision up to which the partition function should be converged.
   * @return The partition function.
   */
  double calculateQuantumPartitioningFunction(double temperature, unsigned int minQuantumNumber = 0,
                                              unsigned int maxQuantumNumber = std::numeric_limits<unsigned int>::infinity() - 1,
                                              double precision = 1e-9);
  /**
   * @brief Calculate the enthalpy through state enumeration.
   * @param temperature      The temperature.
   * @param minQuantumNumber The minimum quantum number.
   * @param maxQuantumNumber The maximum quantum number. If not provided, the states are enumerated until the enthalpy
   * converges.
   * @param precision        The precision up to which the enthalpy should be converged.
   * @return The enthalpy.
   */
  double calculateQuantumEnthalpy(double temperature, unsigned int minQuantumNumber = 0,
                                  unsigned int maxQuantumNumber = std::numeric_limits<unsigned int>::infinity() - 1,
                                  double precision = 1e-9);
  /**
   * @brief Calculate the entropy through state enumeration.
   * @param temperature      The temperature.
   * @param minQuantumNumber The minimum quantum number.
   * @param maxQuantumNumber The maximum quantum number. If not provided, the states are enumerated until the entropy
   * converges.
   * @param precision        The precision up to which the entropy should be converged.
   * @return The entropy.
   */
  double calculateQuantumEntropy(double temperature, unsigned int minQuantumNumber = 0,
                                 unsigned int maxQuantumNumber = std::numeric_limits<unsigned int>::infinity() - 1,
                                 double precision = 1e-9);
  /**
   * @brief Calculate the enthalpy and partition function through state enumeration.
   * @param temperature      The temperature.
   * @param minQuantumNumber The minimum quantum number.
   * @param maxQuantumNumber The maximum quantum number. If not provided, the states are enumerated until the partition
   * function converges.
   * @param precision        The precision up to which the partition function should be converged.
   * @return The partition function and enthalpy.
   */
  std::tuple<double, double>
  calculatePartitionFunctionAndEnthalpy(double temperature, unsigned int minQuantumNumber = 0,
                                        unsigned int maxQuantumNumber = std::numeric_limits<unsigned int>::infinity() - 1,
                                        double precision = 1e-9);
  /**
   * @brief Use Beyer-Swinehart-Stein-Rabinovitch enumeration.
   * @note https://pubs.aip.org/aip/jcp/article/58/6/2438/84783/Accurate-evaluation-of-internal-energy-level-sums
   *
   * @param energies         The energy range.
   * @param t                The input array with which the convolution is calculated.
   * @param minQuantumNumber The minimum quantum number.
   * @param maxQuantumNumber The maximum quantum number. If not provided, the states are enumerated until their energy
   * exceeds the maximum energy in the energy array.
   * @return The density of states or sum of states (depends on t).
   */
  Eigen::VectorXd bsseAlgorithm(const Eigen::VectorXd& energies, const Eigen::VectorXd& t, unsigned int minQuantumNumber = 0,
                                unsigned int maxQuantumNumber = std::numeric_limits<unsigned int>::infinity() - 1);
  /**
   * @brief Beyer-Swinehart enumeration for equidistant states.
   * @param energies    The energy range.
   * @param t           The input array with which the convolution is calculated.
   * @param energyGap   The energy difference between states.
   * @param degeneracy  The degeneracy of states.
   * @return The density of states or sum of states (depends on t).
   */
  Eigen::VectorXd bsEquidistantStates(const Eigen::VectorXd& energies, const Eigen::VectorXd& t, double energyGap,
                                      unsigned int degeneracy);
  /**
   * @brief Assert that all energies in the energy range are non-negative.
   * @param energies The energies.
   *
   * @throw runtime_error If an energy is negative.
   */
  void assertNonNegativeEnergies(const Eigen::VectorXd& energies);
};

} // namespace Utils
} // namespace Scine
#endif // UTILSOS_DEGREEOFFREEDOM_H
