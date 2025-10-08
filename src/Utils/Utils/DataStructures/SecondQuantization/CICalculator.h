/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CICALCULATOR_H
#define UTILS_CICALCULATOR_H

#include <Core/Interfaces/CalculatorWithReference.h>
#include <vector>

namespace Scine {
namespace Core {
class Calculator;
}
namespace Utils {

struct EigenContainer;

namespace SecondQuantization {

class CISolver;
class ElectronicDeterminant;
struct ElectronicDeterminantHash;
class CasIndicesHandler;
struct MOIndicesForReferences;
struct MOIndicesForDoubles;

class CICalculator final : public Core::CalculatorWithReference {
 public:
  static constexpr const char* model = "ci_calculator";
  /// @brief Default constructor.
  CICalculator();
  /// @brief destructor.
  ~CICalculator() final;
  /**
   * @brief Sets the calculator to be used to perform the reference calculation.
   * @throws std::runtime_error if the referenceCalculator has no Hamiltonian possible property.
   */
  void setReferenceCalculator(std::shared_ptr<Core::Calculator> referenceCalculator) override;
  /**
   * @brief Performs a reference calculation.
   */
  void referenceCalculation() override;

  /**
   * @brief Accessor for the reference calculator.
   * @return Core::Calculator& The reference calculator.
   */
  Core::Calculator& getReferenceCalculator() override;
  /**
   * @brief Constant accessor for the reference calculator.
   * @return const Core::Calculator& The reference calculator.
   */
  const Core::Calculator& getReferenceCalculator() const override;

  /**
   * @brief The main function running the calculation with reference.
   * @returns A const-ref of stored (and newly calculated) Results.
   */
  const Results& calculate() override;

  /**
   * @brief Getter for the name of the calculator with reference.
   * @return Returns the name of the calculator with reference.
   */
  std::string name() const override;

  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Settings& settings() override;
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Settings& settings() const override;
  /**
   * @brief Method to apply the settings stored in the settings data structure.
   */
  void applySettings() override;
  /**
   * @brief Accessor for the saved instance of Utils::Results.
   * @return Utils::Results& The results of the previous calculation.
   */
  Results& results() override;
  /**
   * @brief Constant accessor for the Utils::Results.
   * @return const Utils::Results& The results of the previous calculation.
   */
  const Results& results() const override;
  /**
   * @brief Creates labels for a list of determinants.
   * The labels are in the format
   * "222ab000"
   * where
   * - 2 is a doubly occupied orbital
   * - 0 is a virtual orbital
   * - a is an occupied alpha orbital
   * - b is an occupied beta orbital
   */
  static auto getLabels(const std::vector<ElectronicDeterminant>& basis) -> std::vector<std::string>;

  /**
   * @brief Gets a state as a representation suitable for calculating operators.
   */
  auto getWavefunction(int state) const -> std::unordered_map<ElectronicDeterminant, double, ElectronicDeterminantHash>;

 private:
  auto calculateTransitionDipoles(const EigenContainer& results) const -> Eigen::Matrix3Xd;
  auto calculateSpinMultiplicity(const EigenContainer& results) const -> Eigen::VectorXd;
  void checkReferenceCalculator();
  auto getReferenceIndices(const CasIndicesHandler& handler) const -> MOIndicesForReferences;
  auto getIndices() const -> std::tuple<MOIndicesForReferences, MOIndicesForDoubles>;
  std::unique_ptr<Settings> settings_;
  std::unique_ptr<Results> results_;
  std::shared_ptr<Core::Calculator> calc_;
  std::unique_ptr<CISolver> solver_;
  std::unique_ptr<EigenContainer> eigenpairs_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_CICALCULATOR_H
