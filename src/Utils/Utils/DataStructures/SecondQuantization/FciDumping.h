/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_INTEGRALS_FCIDUMPING_H
#define SCINE_INTEGRALS_FCIDUMPING_H

#include "EriUtilities.h"
#include <Eigen/Core>
#include <string>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

struct FciDumpData;

/**
 * @class FciDumper @file FciDumper.h
 * @brief Class for input and output of FCIDump files for one- and two-bodies integrals.
 * In unrestricted calculations, the order of the integrals is
 * 1 1alpha
 * 2 1beta
 * 3 2alpha
 * 4 2beta
 * ...
 */
namespace FciDumping {

/**
 * @brief Type definitions
 * @{
 */
template<Reference restrictedness>
using OneElectronMatrixType = typename MOTypeTrait<restrictedness>::OneElectronMatrixType;
template<Reference restrictedness>
using TwoElectronMatrixType = typename MOTypeTrait<restrictedness>::Type;
template<Reference restrictedness>
using PackageType = std::tuple<OneElectronMatrixType<restrictedness>, TwoElectronMatrixType<restrictedness>, FciDumpData>;
/**
 * @}
 */
/**
 * @brief Static variable to control the threshold value for integrals to be printed.
 */
static constexpr double integralThreshold = 1e-16;
/**
 * @brief Write one- and two-bodies integrals to a file.
 *
 * @param fileName The file path.
 * @param hCore The one-body integrals matrix.
 * @param eris The two-bodies integrals supermatrix.
 * @param fciDumpData Struct containing auxiliary infos as number of orbitals, number of electrons, symmetries...
 * @throws std::runtime_error If the file could not be created.
 */
template<Reference restrictedness>
void write(const std::string& filename, const OneElectronMatrixType<restrictedness>& hCore,
           const TwoElectronMatrixType<restrictedness>& eris, const FciDumpData& fciDumpData);
/**
 * @brief Write one- and two-bodies integrals to a stream.
 *
 * @param out  The output stream.
 * @param hCore The one-body integrals matrix.
 * @param eris The two-bodies integrals supermatrix.
 * @param fciDumpData Struct containing auxiliary infos as number of orbitals, number of electrons, symmetries...
 */
template<Reference restrictedness>
void write(std::ostream& out, const OneElectronMatrixType<restrictedness>& hCore,
           const TwoElectronMatrixType<restrictedness>& eris, const FciDumpData& fciDumpData);
/**
 * @brief Read one- and two-bodies integrals from a file.
 *
 * @param fileName The file path.
 * @throws std::runtime_error If the file could not be opened.
 * @return Returns a std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> with the one-, two-bodies integrals and auxiliary
 * infos. In the unrestricted case, an unordered map with the different spin components is returned
 */
template<Reference restrictedness>
auto read(const std::string& filename) -> PackageType<restrictedness>;
/**
 * @brief Read one- and two-bodies integrals from a stream.
 *
 * @param in The input stream.
 * @return Returns a std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, FciDumpData> with one-, two-bodies integrals and
 * auxiliary infos. If unrestricted, an unordered map of Eigen::MatrixXd is returned with the different spin component
 */
template<Reference restrictedness>
auto read(std::istream& in) -> PackageType<restrictedness>;

template<Reference restrictedness>
void writeHeader(std::ostream& out, const FciDumpData& fciDumpData);
template<Reference restrictedness>
void writeOneElectronIntegrals(std::ostream& out, const OneElectronMatrixType<restrictedness>& hCore);
template<Reference restrictedness>
void writeTwoElectronIntegrals(std::ostream& out, const TwoElectronMatrixType<restrictedness>& eris);
template<Reference restrictedness>
auto contractInactiveOrbitalsInCore(const OneElectronMatrixType<restrictedness>& hCore,
                                    const TwoElectronMatrixType<restrictedness>& eris) -> double;
template<Reference restrictedness>
auto parseHeader(std::istream& in) -> FciDumpData;
template<Reference restrictedness>
auto parseIntegrals(std::istream& in, FciDumpData& fciDumpData)
    -> std::tuple<OneElectronMatrixType<restrictedness>, TwoElectronMatrixType<restrictedness>>;

} // namespace FciDumping
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#include "FciDumpingImpl.h"
#endif // SCINE_INTEGRALS_FCIDUMPING_H
