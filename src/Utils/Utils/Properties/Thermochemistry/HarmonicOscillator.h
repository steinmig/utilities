/**
 * @file HarmonicOscillator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_HARMONICOSCILLATOR_H
#define UTILSOS_HARMONICOSCILLATOR_H

#include "Utils/Properties/Thermochemistry/DegreeOfFreedom.h"
#include "Utils/Typenames.h"
#include <memory>

namespace Scine {
namespace Utils {
class NormalModesContainer;
class AtomCollection;

/**
 * @class
 * @brief This class represents a molecule's vibrational degrees of freedom as harmonic oscillations.
 */
class HarmonicOscillator : public std::enable_shared_from_this<HarmonicOscillator>, public DegreeOfFreedom {
 public:
  /**
   * @brief Constructor from normal modes.
   * @param normalModesContainer The normal mode container.
   * @param classicalMechanics   If true, classical mechanics are assumed.
   */
  explicit HarmonicOscillator(std::shared_ptr<NormalModesContainer> normalModesContainer, bool classicalMechanics = false);
  /**
   * @brief Constructor.
   * @param hessian            The hessian.
   * @param atoms              The molecule's atom collection.
   * @param classicalMechanics If true, classical mechanics are assumed.
   */
  HarmonicOscillator(const HessianMatrix& hessian, const AtomCollection& atoms, bool classicalMechanics = false);
  /**
   * @brief Default destructor.
   */
  ~HarmonicOscillator();

  /**
   * @brief Getter for the density of states at the given energies.
   * @param energies The energies.
   * @return The density of states.
   */
  virtual Eigen::VectorXd densityOfState(const Eigen::VectorXd& energies) override;
  /**
   * @brief Getter for the sum of states at the given energies.
   * @param energies The energies.
   * @return The sum of states.
   */
  virtual Eigen::VectorXd sumOfStates(const Eigen::VectorXd& energies) override;
  /**
   * @brief Getter for the partition function.
   * @param state The thermodynamic reference state.
   * @return The partition function (in atomic units).
   */
  virtual double partitionFunction(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the entropy.
   * @param state The thermodynamic reference state.
   * @return The entropy (in atomic units).
   */
  virtual double entropy(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the enthalpy.
   * @param state The thermodynamic reference state.
   * @return The enthalpy (in atomic units).
   */
  virtual double enthalpy(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the heat capacity at constant pressure.
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCp(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Getter for the heat capacity at constant volume.
   * @param state The thermodynamic reference state.
   * @return The heat capacity (in atomic units).
   */
  virtual double heatCapacityCv(const ThermodynamicReferenceState& state) override;
  /**
   * @brief Returns a classical harmonic oscillator object. If the object is already classical, a shared pointer
   *        to *this is returned. Otherwise, a new harmonic oscillator object may be constructed.
   * @return A HarmonicOscillator with classicalMechanics_=true.
   */
  std::shared_ptr<HarmonicOscillator> getClassicalVariant();
  /**
   * @brief Getter for the normal modes container.
   * @return The normal modes container.
   */
  const NormalModesContainer& getNormalModesContainer();

 private:
  std::shared_ptr<NormalModesContainer> normalModesContainer_;
  std::shared_ptr<HarmonicOscillator> classicalVariant_;
  std::unique_ptr<double> waveNumberProduct_;
  std::unique_ptr<unsigned int> nRealModes_;
  std::unique_ptr<double> zeroPointVibrationalCorrection_;
  const bool classicalMechanics_;
  void analyzeVibrations();
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_HARMONICOSCILLATOR_H
