/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_FCIDUMPINGIMPL_H
#define SCINE_FCIDUMPINGIMPL_H

#include <Utils/IO/Regex.h>
#include <Utils/Technical/ScopedLocale.h>
#include <fstream>
#include <iomanip>
#include <regex>

namespace Scine {
namespace Utils {
namespace SecondQuantization {
namespace FciDumping {

template<Reference restrictedness>
void write(const std::string& filename, const OneElectronMatrixType<restrictedness>& hCore,
           const TwoElectronMatrixType<restrictedness>& eris, const FciDumpData& fciDumpData) {
  std::ofstream fout(filename);

  if (!fout.is_open()) {
    throw std::runtime_error("Problem when opening/creating file " + filename);
  }

  return write<restrictedness>(fout, hCore, eris, fciDumpData);
}

template<Reference restrictedness>
void write(std::ostream& out, const OneElectronMatrixType<restrictedness>& hCore,
           const TwoElectronMatrixType<restrictedness>& eris, const FciDumpData& fciDumpData) {
  auto scopedLocale = ScopedLocale::cLocale();
  writeHeader<restrictedness>(out, fciDumpData);
  out << std::scientific << std::left << std::setprecision(15);
  writeTwoElectronIntegrals<restrictedness>(out, eris);
  writeOneElectronIntegrals<restrictedness>(out, hCore);

  // Write core energy
  out << std::setw(25) << fciDumpData.coreEnergy << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 0
      << std::setw(10) << 0;
}

template<Reference restrictedness>
void writeHeader(std::ostream& out, const FciDumpData& fciDumpData) {
  out << "&FCI NORB=" << fciDumpData.nOrbitals << ",NELEC=" << fciDumpData.nElectrons
      << ",MS2=" << fciDumpData.spinPolarization << "," << std::endl;
  out << " ORBSYM=";
  for (auto sym : fciDumpData.orbitalSymmetry) {
    out << sym << ",";
  }
  out << std::endl;
  out << " ISYM=" << fciDumpData.wfSymIfLowestOccupied << std::endl;
  if (fciDumpData.unrestrictedReference) {
    assert(restrictedness == Reference::Unrestricted);
    out << "UHF" << std::endl;
  }
  out << "&END" << std::endl;
}

template<>
inline void writeOneElectronIntegrals<Reference::Restricted>(std::ostream& out,
                                                             const OneElectronMatrixType<Reference::Restricted>& hCore) {
  for (int row = 0; row < hCore.rows(); ++row) {
    for (int col = 0; col <= row; ++col) {
      if (std::abs(hCore(row, col)) > integralThreshold) {
        out << std::setw(25) << hCore(row, col) << std::setw(10) << row + 1 << std::setw(10) << col + 1 << std::setw(10)
            << 0 << std::setw(10) << 0 << std::endl;
      }
    }
  }
}

template<>
inline void writeOneElectronIntegrals<Reference::Unrestricted>(std::ostream& out,
                                                               const OneElectronMatrixType<Reference::Unrestricted>& hCore) {
  const Eigen::MatrixXd& alphaMatrix = hCore.at(SpinComponent::Alpha);
  const Eigen::MatrixXd& betaMatrix = hCore.at(SpinComponent::Beta);
  for (int row = 0; row < alphaMatrix.rows(); ++row) {
    int rowAlpha = 2 * row;
    int rowBeta = 2 * row + 1;
    for (int col = 0; col < row + 1; ++col) {
      int colAlpha = 2 * col;
      int colBeta = 2 * col + 1;
      out << std::setw(25) << alphaMatrix(row, col) << std::setw(10) << rowAlpha + 1 << std::setw(10) << colAlpha + 1
          << std::setw(10) << 0 << std::setw(10) << 0 << std::endl;
      out << std::setw(25) << betaMatrix(row, col) << std::setw(10) << rowBeta + 1 << std::setw(10) << colBeta + 1
          << std::setw(10) << 0 << std::setw(10) << 0 << std::endl;
    }
  }
}

template<>
inline void writeTwoElectronIntegrals<Reference::Restricted>(std::ostream& out,
                                                             const TwoElectronMatrixType<Reference::Restricted>& eris) {
  int nOrb = std::sqrt(eris.cols());
  for (int row = 0; row < eris.rows(); ++row) {
    int p = row / nOrb;
    int q = row % nOrb;
    for (int col = 0; col <= row; ++col) {
      int r = col / nOrb;
      int s = col % nOrb;
      auto index = getMappedIndex<ERISymmetry::eightfold>({{p, q, r, s}});
      if (p == index[0] && q == index[1] && r == index[2] && s == index[3] && std::abs(eris(row, col)) > integralThreshold) {
        out << std::setw(25) << eris(row, col) << std::setw(10) << p + 1 << std::setw(10) << q + 1 << std::setw(10)
            << r + 1 << std::setw(10) << s + 1 << std::endl;
      }
    }
  }
}

template<>
inline void writeTwoElectronIntegrals<Reference::Unrestricted>(std::ostream& out,
                                                               const TwoElectronMatrixType<Reference::Unrestricted>& eris) {
  const Eigen::MatrixXd& eriAA = eris.at({SpinComponent::Alpha, SpinComponent::Alpha});
  const Eigen::MatrixXd& eriAB = eris.at({SpinComponent::Alpha, SpinComponent::Beta});
  const Eigen::MatrixXd& eriBA = eris.at({SpinComponent::Beta, SpinComponent::Alpha});
  const Eigen::MatrixXd& eriBB = eris.at({SpinComponent::Beta, SpinComponent::Beta});
  int nOrb = std::sqrt(eriAA.cols());
  for (int row = 0; row < eriAA.rows(); ++row) {
    int p = row / nOrb;
    int q = row % nOrb;
    for (int col = 0; col <= row; ++col) {
      int r = col / nOrb;
      int s = col % nOrb;
      auto index = getMappedIndex<ERISymmetry::eightfold>({{p, q, r, s}});
      if (p == index[0] && q == index[1] && r == index[2] && s == index[3]) {
        if (std::abs(eriAA(row, col)) > integralThreshold) {
          out << std::setw(25) << eriAA(row, col) << std::setw(10) << p * 2 + 1 << std::setw(10) << q * 2 + 1
              << std::setw(10) << r * 2 + 1 << std::setw(10) << s * 2 + 1 << std::endl;
        }
        if (std::abs(eriAB(row, col)) > integralThreshold) {
          out << std::setw(25) << eriAB(row, col) << std::setw(10) << p * 2 + 1 << std::setw(10) << q * 2 + 1
              << std::setw(10) << r * 2 + 2 << std::setw(10) << s * 2 + 2 << std::endl;
        }
        if (std::abs(eriBA(row, col)) > integralThreshold) {
          out << std::setw(25) << eriBA(row, col) << std::setw(10) << p * 2 + 2 << std::setw(10) << q * 2 + 2
              << std::setw(10) << r * 2 + 1 << std::setw(10) << s * 2 + 1 << std::endl;
        }
        if (std::abs(eriBB(row, col)) > integralThreshold) {
          out << std::setw(25) << eriBB(row, col) << std::setw(10) << p * 2 + 2 << std::setw(10) << q * 2 + 2
              << std::setw(10) << r * 2 + 2 << std::setw(10) << s * 2 + 2 << std::endl;
        }
      }
    }
  }
}

template<Reference restrictedness>
auto read(const std::string& filename) -> PackageType<restrictedness> {
  std::ifstream fin(filename);

  if (!fin.is_open()) {
    throw std::runtime_error("Problem when opening file " + filename);
  }

  return read<restrictedness>(fin);
}

template<Reference restrictedness>
auto read(std::istream& in) -> PackageType<restrictedness> {
  OneElectronMatrixType<restrictedness> hCoreMo;
  TwoElectronMatrixType<restrictedness> eris;
  auto fciDumpData = parseHeader<restrictedness>(in);
  std::tie(hCoreMo, eris) = parseIntegrals<restrictedness>(in, fciDumpData);
  return std::make_tuple(hCoreMo, eris, fciDumpData);
}

template<Reference restrictedness>
auto parseHeader(std::istream& in) -> FciDumpData {
  auto scopedLocale = ScopedLocale::cLocale();
  std::string number = Regex::capturingIntegerNumber();
  std::regex r(number);

  FciDumpData fciDumpData{};

  std::string line;
  // skip lines until &FCI is found or eof.
  do {
    std::getline(in, line);
  } while ((line.find("&FCI") == std::string::npos) && !in.eof());

  do {
    std::string field;
    auto spaceStrippedLine = std::stringstream(line);
    if (spaceStrippedLine.str().find("ORBSYM") != std::string::npos) {
      std::smatch match;
      while (std::regex_search(line, match, r)) {
        fciDumpData.orbitalSymmetry.push_back(std::stoi(match[0]));
        line = match.suffix();
      }
    }
    else {
      while (std::getline(spaceStrippedLine, field, ',')) {
        if (field.find("NORB") != std::string::npos) {
          std::smatch match;
          std::regex_search(field, match, r);
          fciDumpData.nOrbitals = std::stoi(match[0]);
        }
        else if (field.find("NELEC") != std::string::npos) {
          std::smatch match;
          std::regex_search(field, match, r);
          fciDumpData.nElectrons = std::stoi(match[0]);
        }
        else if (field.find("MS2") != std::string::npos) {
          auto pos = field.find("MS2") + 3;
          auto substring = field.substr(pos);
          std::smatch match;
          std::regex_search(substring, match, r);
          fciDumpData.spinPolarization = std::stoi(match[0]);
        }
        else if (field.find("ISYM") != std::string::npos) {
          std::smatch match;
          std::regex_search(field, match, r);
          fciDumpData.wfSymIfLowestOccupied = std::stoi(match[0]);
        }
        else if (field.find("UHF") != std::string::npos) {
          fciDumpData.unrestrictedReference = true;
          if (restrictedness == Reference::Restricted) {
            throw std::runtime_error("UHF FciDump, but restricted data detected.");
          }
        }
      }
    }
    std::getline(in, line);
  } while (line.find("&END") == std::string::npos && !in.eof());

  return fciDumpData;
}

namespace detail {
inline void fill(OneElectronMatrixType<Reference::Restricted>& lhs, int i, int j, double value) {
  lhs(i - 1, j - 1) = value;
}
inline void fill(OneElectronMatrixType<Reference::Unrestricted>& lhs, int i, int j, double value) {
  if ((i - 1) % 2 == 0) {
    assert((j - 1) % 2 == 0);
    lhs.at(SpinComponent::Alpha)((i - 1) / 2, (j - 1) / 2) = value;
  }
  else {
    assert((i - 1) % 2 == 1);
    assert((j - 1) % 2 == 1);
    lhs.at(SpinComponent::Beta)((i - 2) / 2, (j - 2) / 2) = value;
  }
}

inline void fill(TwoElectronMatrixType<Reference::Restricted>& lhs, int nMOs, int i, int j, int k, int l, double value) {
  for (auto idx : getSymmetricIndices<ERISymmetry::eightfold>({{i - 1, j - 1, k - 1, l - 1}})) {
    lhs(idx[0] * nMOs + idx[1], idx[2] * nMOs + idx[3]) = value;
  }
}
inline void fill(TwoElectronMatrixType<Reference::Unrestricted>& lhs, int nMOs, int i, int j, int k, int l, double value) {
  SpinComponent spin1 = (i - 1) % 2 == 0 ? SpinComponent::Alpha : SpinComponent::Beta;
  SpinComponent spin2 = (k - 1) % 2 == 0 ? SpinComponent::Alpha : SpinComponent::Beta;
  int p = spin1 == SpinComponent::Alpha ? (i - 1) / 2 : (i - 2) / 2;
  int q = spin1 == SpinComponent::Alpha ? (j - 1) / 2 : (j - 2) / 2;
  int r = spin2 == SpinComponent::Alpha ? (k - 1) / 2 : (k - 2) / 2;
  int s = spin2 == SpinComponent::Alpha ? (l - 1) / 2 : (l - 2) / 2;
  for (auto idx : getSymmetricIndices<ERISymmetry::eightfold>({{p, q, r, s}})) {
    lhs.at({spin1, spin2})(idx[0] * nMOs + idx[1], idx[2] * nMOs + idx[3]) = value;
  }
}

inline void init(OneElectronMatrixType<Reference::Restricted>& oneBody,
                 TwoElectronMatrixType<Reference::Restricted>& eri, int nMOs) {
  oneBody = Eigen::MatrixXd::Zero(nMOs, nMOs);
  eri = Eigen::MatrixXd::Zero(nMOs * nMOs, nMOs * nMOs);
}
inline void init(OneElectronMatrixType<Reference::Unrestricted>& oneBody,
                 TwoElectronMatrixType<Reference::Unrestricted>& eri, int nMOs) {
  for (auto spin1 : {SpinComponent::Alpha, SpinComponent::Beta}) {
    oneBody[spin1] = Eigen::MatrixXd::Zero(nMOs, nMOs);
    for (auto spin2 : {SpinComponent::Alpha, SpinComponent::Beta}) {
      eri[{spin1, spin2}] = Eigen::MatrixXd::Zero(nMOs * nMOs, nMOs * nMOs);
    }
  }
}
} // namespace detail

template<Reference restrictedness>
auto parseIntegrals(std::istream& in, FciDumpData& fciDumpData)
    -> std::tuple<OneElectronMatrixType<restrictedness>, TwoElectronMatrixType<restrictedness>> {
  OneElectronMatrixType<restrictedness> hCoreMo;
  TwoElectronMatrixType<restrictedness> eris;
  detail::init(hCoreMo, eris, fciDumpData.nOrbitals);
  std::string line;
  while (std::getline(in, line)) {
    std::stringstream buffer(line);
    double value = 0.0;
    int i = 0, j = 0, k = 0, l = 0;
    buffer >> value >> i >> j >> k >> l;
    if (i != 0 && j != 0 && k != 0 && l != 0) {
      detail::fill(eris, fciDumpData.nOrbitals, i, j, k, l, value);
    }
    else if (i == 0 && j == 0 && k == 0 && l == 0) {
      fciDumpData.coreEnergy = value;
    }
    else if (k == 0 && l == 0) {
      detail::fill(hCoreMo, i, j, value);
      detail::fill(hCoreMo, j, i, value);
    }
  }
  return std::make_tuple(hCoreMo, eris);
}

} // namespace FciDumping
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // SCINE_FCIDUMPINGIMPL_H
