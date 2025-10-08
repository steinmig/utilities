/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

// #include "Hamiltonian.h"
#include "ElectronicDeterminant.h"
#include "EriUtilities.h"
#include "ExcitationHelpers.h"
#include "HeatBathCIContainer.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <boost/functional/hash.hpp>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <utility>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

template<class Crtp>
Hamiltonian<Crtp>::Hamiltonian(int moBasisSize) : nMOs_(moBasisSize) {
  oneElectronMatrix_ = Eigen::MatrixXd::Zero(moBasisSize, moBasisSize);
}

template<class Crtp>
Hamiltonian<Crtp>::Hamiltonian(const Specifier& specifier)
  : Hamiltonian(specifier.cas.getActiveOrbitalsCoefficients().numberOrbitals()) {
  addConstantTerm(specifier.coreEnergy);
  if (specifier.integrals.oneBodyIntegrals) {
    addOneBodyContribution(specifier.integrals.get(MoIntegrals::Type::OneBody));
  }
  if (specifier.integrals.twoBodyIntegrals) {
    addTwoBodyContribution(specifier.integrals.get(MoIntegrals::Type::TwoBody));
  }
}

template<class Crtp>
void Hamiltonian<Crtp>::print(std::ostream& out) const {
  out << "One-body term" << std::endl;
  for (int iRow = 0; iRow < oneElectronMatrix_.rows(); iRow++) {
    for (int iCol = 0; iCol < oneElectronMatrix_.cols(); iCol++) {
      out << iRow << " " << iCol << " " << oneElectronMatrix_(iRow, iCol) << std::endl;
    }
  }
  out << "Two-body term" << std::endl;
  for (const auto& iMap : twoElectronMatrix_) {
    out << iMap.first[0] << " " << iMap.first[1] << " " << iMap.first[2] << " " << iMap.first[3] << " " << iMap.second
        << std::endl;
  }
}

template<class Crtp>
void Hamiltonian<Crtp>::addTerm(int i, int j, double coeff) {
  oneElectronMatrix_(std::min(i, j), std::max(i, j)) = coeff;
}

template<class Crtp>
void Hamiltonian<Crtp>::addOneBodyContribution(const Eigen::MatrixXd& inputOneBody) {
  for (int iRow = 0; iRow < inputOneBody.rows(); iRow++) {
    for (int iCol = iRow; iCol < inputOneBody.cols(); iCol++) {
      if (iRow == iCol) {
        addTerm(iRow, iCol, inputOneBody(iRow, iCol));
      }
      else {
        if (std::fabs(inputOneBody(iRow, iCol)) > thresholdForTerms_ && std::fabs(inputOneBody(iCol, iRow)) < thresholdForTerms_) {
          addTerm(iRow, iCol, inputOneBody(iRow, iCol));
        }
        else if (std::fabs(inputOneBody(iRow, iCol)) < thresholdForTerms_ &&
                 std::fabs(inputOneBody(iCol, iRow)) > thresholdForTerms_) {
          addTerm(iRow, iCol, inputOneBody(iCol, iRow));
        }
        else if (std::fabs(inputOneBody(iRow, iCol)) > thresholdForTerms_ &&
                 std::fabs(inputOneBody(iCol, iRow)) > thresholdForTerms_) {
          if (std::fabs(inputOneBody(iCol, iRow) - inputOneBody(iRow, iCol)) < thresholdForTerms_) {
            addTerm(iCol, iRow, inputOneBody(iRow, iCol));
          }
          else {
            throw std::runtime_error("Incoherence in the input data");
          }
        }
      }
    }
  }
}

template<class Crtp>
void Hamiltonian<Crtp>::addTwoBodyContribution(const TwoBodyContainer& inputTwoBody, bool isPhys) {
  twoElectronMatrix_ = {};
  // Actual loading of the HBCI structure
  for (const auto& eri : inputTwoBody) {
    if (std::fabs(eri.second) > thresholdForTerms_) {
      auto index = eri.first;
      // Now brings everything in physics notation (this is how data MUST BE stored)
      if (!isPhys) {
        std::swap(index[1], index[2]);
      }
      // We add the exchange only if the excitation involves same spin orbitals
      twoElectronMatrix_[bringToNormalForm(index)] = eri.second;
    }
  }
  populateHBCIContainer();
}

template<class Crtp>
void Hamiltonian<Crtp>::addTwoBodyContribution(const Eigen::MatrixXd& inputTwoBody) {
  twoElectronMatrix_ = {};
  for (int i = 0; i < nMOs_; ++i) {
    for (int j = 0; j <= i; ++j) {
      int bra = i * nMOs_ + j;
      for (int k = 0; k <= i; ++k) {
        for (int l = 0; l <= k; ++l) {
          int ket = k * nMOs_ + l;
          double value = inputTwoBody(bra, ket);
          if (std::fabs(value) > thresholdForTerms_) {
            std::array<int, 4> index = {{i, j, k, l}};
            // Now brings everything in physics notation (this is how data MUST BE stored)
            std::swap(index[1], index[2]);
            // We add the exchange only if the excitation involves same spin orbitals
            twoElectronMatrix_[bringToNormalForm(index)] = value;
          }
        }
      }
    }
  }
  // Actual loading of the HBCI structure
  populateHBCIContainer();
}

template<class Crtp>
void Hamiltonian<Crtp>::addConstantTerm(double coreEnergy) {
  coreEnergy_ = coreEnergy;
}

template<class Crtp>
void Hamiltonian<Crtp>::populateHBCIContainer() {
  // Same spin
  for (auto el : twoElectronMatrix_) {
    // The HBCI container must be populated with the "true" excitations.
    // For this reason, if we have a ijkl integral, we manually add all combinations
    // i.e., ilkj, kjil, klij.
    // When we then take the data from the container to generate an excitation,
    // we consider two cases:
    // 1) alpha/beta: we manually add the symmetric excitation if i=j.
    // 2) alpha/alpha or beta/beta: we manually remove the case in which the starting
    //    or finishing orbital are the same.
    if (el.first[0] != el.first[2] && el.first[1] != el.first[3]) {
      // Prepares all permutations
      IndexType<4> idx1 = {el.first[0], el.first[1], el.first[2], el.first[3]};
      if (el.first[0] > el.first[1]) {
        std::swap(idx1[0], idx1[1]);
        std::swap(idx1[2], idx1[3]);
      }
      IndexType<4> idx2 = {el.first[0], el.first[3], el.first[2], el.first[1]};
      if (el.first[0] > el.first[3]) {
        std::swap(idx2[0], idx2[1]);
        std::swap(idx2[2], idx2[3]);
      }
      IndexType<4> idx3 = {el.first[2], el.first[1], el.first[0], el.first[3]};
      if (el.first[2] > el.first[1]) {
        std::swap(idx3[0], idx3[1]);
        std::swap(idx3[2], idx3[3]);
      }
      IndexType<4> idx4 = {el.first[2], el.first[3], el.first[0], el.first[1]};
      if (el.first[2] > el.first[3]) {
        std::swap(idx4[0], idx4[1]);
        std::swap(idx4[2], idx4[3]);
      }
      // ij --> kl
      addTermToHBCIContainer(idx1[0], idx1[1], idx1[2], idx1[3], el.second);
      // il --> kj
      if (idx2 != idx1) {
        addTermToHBCIContainer(idx2[0], idx2[1], idx2[2], idx2[3], el.second);
      }
      // kj --> il
      if (idx3 != idx2 && idx3 != idx1) {
        addTermToHBCIContainer(idx3[0], idx3[1], idx3[2], idx3[3], el.second);
      }
      // kl --> ij
      if (idx4 != idx3 && idx4 != idx2 && idx4 != idx1) {
        addTermToHBCIContainer(idx4[0], idx4[1], idx4[2], idx4[3], el.second);
      }
    }
  }
}

template<class Crtp>
int Hamiltonian<Crtp>::getNumberOneBody() const {
  return oneElectronMatrix_.size();
}

template<class Crtp>
int Hamiltonian<Crtp>::getNumberTwoBody() const {
  return twoElectronMatrix_.size();
}

template<class Crtp>
double Hamiltonian<Crtp>::calculateMatrixElement(const ElectronicDeterminant& det1, const ElectronicDeterminant& det2) const {
  // Initialization
  const auto excitation = det1.getExcitations(det2, false);
  return calculateMatrixElement(excitation, det1.getAllOccupied(SpinComponent::Alpha),
                                det1.getAllOccupied(SpinComponent::Beta));
}

template<class Crtp>
int Hamiltonian<Crtp>::calculateExcitationSign(const std::vector<int>& lst_occupied_alpha,
                                               const std::vector<int>& lst_occupied_beta, const SingleExcitation& exc1) const {
  const auto& lst_samespin = (exc1.orbitalSpin == SpinComponent::Alpha) ? lst_occupied_alpha : lst_occupied_beta;
  const auto& lst_otherspin = (exc1.orbitalSpin == SpinComponent::Alpha) ? lst_occupied_beta : lst_occupied_alpha;
  int initialPtr = exc1.orbitalSpin == SpinComponent::Alpha ? std::min(exc1.occupiedOrbital, exc1.virtualOrbital)
                                                            : std::min(exc1.occupiedOrbital, exc1.virtualOrbital) + 1;
  int finalPtr = exc1.orbitalSpin == SpinComponent::Alpha ? std::max(exc1.occupiedOrbital, exc1.virtualOrbital) - 1
                                                          : std::max(exc1.occupiedOrbital, exc1.virtualOrbital);
  int nInterSameSpin = std::count_if(lst_samespin.begin(), lst_samespin.end(),
                                     [initialPtr, finalPtr](int i) { return i >= initialPtr && i <= finalPtr; });
  int nInterOtherSpin = std::count_if(lst_otherspin.begin(), lst_otherspin.end(),
                                      [initialPtr, finalPtr](int i) { return i >= initialPtr && i <= finalPtr; });
  int coeff = (nInterSameSpin % 2 == 0) ? 1 : -1;
  coeff *= (nInterOtherSpin % 2 == 0) ? 1 : -1;
  // Corrects for the offset introduced by the action of the operators before the filling (i.e., the occupation)
  // changes by virtue of the b operators.
  if (exc1.orbitalSpin == SpinComponent::Alpha) {
    coeff *= -1;
  }
  if (exc1.occupiedOrbital > exc1.virtualOrbital) {
    coeff *= -1;
  }
  return coeff;
}

template<class Crtp>
int Hamiltonian<Crtp>::calculateExcitationSign(const std::vector<int>& lst_occupied_alpha,
                                               const std::vector<int>& lst_occupied_beta, const SingleExcitation& exc1,
                                               const SingleExcitation& exc2) const {
  int itFirstBegin{}, itFirstEnd{}, itSecondBegin{}, itSecondEnd{}, offset{};
  if (exc1.orbitalSpin == exc2.orbitalSpin) {
    // -- Alpha-alpha and Beta-beta excitations --
    itFirstBegin = (exc1.orbitalSpin == SpinComponent::Alpha) ? std::min(exc1.virtualOrbital, exc1.occupiedOrbital)
                                                              : std::min(exc1.virtualOrbital, exc1.occupiedOrbital) + 1;
    itFirstEnd = (exc1.orbitalSpin == SpinComponent::Alpha) ? std::max(exc1.virtualOrbital, exc1.occupiedOrbital) - 1
                                                            : std::max(exc1.virtualOrbital, exc1.occupiedOrbital);
    itSecondBegin = (exc2.orbitalSpin == SpinComponent::Alpha) ? std::min(exc2.virtualOrbital, exc2.occupiedOrbital)
                                                               : std::min(exc2.virtualOrbital, exc2.occupiedOrbital) + 1;
    itSecondEnd = (exc2.orbitalSpin == SpinComponent::Alpha) ? std::max(exc2.virtualOrbital, exc2.occupiedOrbital) - 1
                                                             : std::max(exc2.virtualOrbital, exc2.occupiedOrbital);
    // Commutator-related offset
    if (exc1.orbitalSpin == SpinComponent::Alpha) {
      if (exc1.occupiedOrbital < exc1.virtualOrbital)
        offset += 1;
      if (exc2.occupiedOrbital < exc2.virtualOrbital)
        offset += 1;
    }
    else {
      if (exc1.occupiedOrbital > exc1.virtualOrbital)
        offset += 1;
      if (exc2.occupiedOrbital > exc2.virtualOrbital)
        offset += 1;
    }
    // Final offset
    if (itFirstBegin <= exc2.occupiedOrbital && exc2.occupiedOrbital <= itFirstEnd)
      offset += 1;
    if (itFirstBegin <= exc2.virtualOrbital && exc2.virtualOrbital <= itFirstEnd)
      offset += 1;
  }
  else {
    // -- Mixed Alpha-Beta excitation --
    const auto& iAlphaExc = exc1.orbitalSpin == SpinComponent::Alpha ? exc1 : exc2;
    const auto& iBetaExc = exc1.orbitalSpin == SpinComponent::Beta ? exc1 : exc2;
    int iMinAlpha = std::min(iAlphaExc.occupiedOrbital, iAlphaExc.virtualOrbital);
    int iMaxAlpha = std::max(iAlphaExc.occupiedOrbital, iAlphaExc.virtualOrbital);
    int iMinBeta = std::min(iBetaExc.occupiedOrbital, iBetaExc.virtualOrbital);
    int iMaxBeta = std::max(iBetaExc.occupiedOrbital, iBetaExc.virtualOrbital);
    // Here we add the offset that are present in the definition of the JW mapping
    itFirstBegin = iMinAlpha;
    itFirstEnd = iMaxAlpha - 1;
    itSecondBegin = iMinBeta + 1;
    itSecondEnd = iMaxBeta;
    // Commutator-related offset
    // Note that in the alpha case we have an additional offset due to the fact that the JW
    // mapping adds a filling just in front of the destructor, so the filling operator acts on
    // a determinant with an electron less. This introduces the offsets.
    if (iAlphaExc.occupiedOrbital < iAlphaExc.virtualOrbital)
      offset += 1;
    if (iBetaExc.occupiedOrbital > iBetaExc.virtualOrbital)
      offset += 1;
    // Final offset
    if (itFirstBegin <= iMinBeta && iMinBeta <= itFirstEnd)
      offset += 1;
    if (itFirstBegin <= iMaxBeta && iMaxBeta <= itFirstEnd)
      offset += 1;
  }
  // Final rescaling
  auto firstCheck = [itFirstBegin, itFirstEnd](int i) -> bool { return i >= itFirstBegin && i <= itFirstEnd; };
  auto secondCheck = [itSecondBegin, itSecondEnd](int i) -> bool { return i >= itSecondBegin && i <= itSecondEnd; };
  int nToCheck = std::count_if(lst_occupied_alpha.begin(), lst_occupied_alpha.end(), firstCheck) +
                 std::count_if(lst_occupied_beta.begin(), lst_occupied_beta.end(), firstCheck) +
                 std::count_if(lst_occupied_alpha.begin(), lst_occupied_alpha.end(), secondCheck) +
                 std::count_if(lst_occupied_beta.begin(), lst_occupied_beta.end(), secondCheck);
  nToCheck += offset;
  int ret = 1;
  if (nToCheck % 2 == 1)
    ret = -1;
  return ret;
}

template<class Crtp>
double Hamiltonian<Crtp>::calculateMatrixElement(const Excitation& excitations, const std::vector<int>& i_occ_alpha,
                                                 const std::vector<int>& i_occ_beta) const {
  if (excitations.cardinality() == 0) {
    return matrixElementSameDeterminant(i_occ_alpha, i_occ_beta) + coreEnergy_;
  }
  if (excitations.cardinality() == 1) {
    return matrixElementSingleExcitation(excitations.getSingleExcitation(0), i_occ_alpha, i_occ_beta);
  }
  if (excitations.cardinality() == 2) {
    return matrixElementDoubleExcitation(excitations.getSingleExcitation(0), excitations.getSingleExcitation(1),
                                         i_occ_alpha, i_occ_beta);
  }
  return 0.0;
}

template<class Crtp>
std::vector<ElectronicDeterminant> Hamiltonian<Crtp>::generateAllConnected(const ElectronicDeterminant& det, bool addSame) const {
  return derived().generateAllConnectedImpl(det, addSame);
}

template<class Crtp>
std::vector<ElectronicDeterminant> Hamiltonian<Crtp>::generateAllConnectedImpl(const ElectronicDeterminant& det,
                                                                               bool addSame) const {
  auto singles = generateSingleExcitations(det, addSame);
  auto doubles = generateDoubleExcitations(det, false);
  std::vector<ElectronicDeterminant> ret;
  ret.reserve(singles.size() + doubles.size());
  ret.insert(ret.end(), std::make_move_iterator(singles.begin()), std::make_move_iterator(singles.end()));
  ret.insert(ret.end(), std::make_move_iterator(doubles.begin()), std::make_move_iterator(doubles.end()));
  return ret;
}

template<class Crtp>
std::vector<ElectronicDeterminant>
Hamiltonian<Crtp>::generateAllConnected(const ElectronicDeterminant& det,
                                        std::map<SpinComponent, std::vector<int>> listDoublesOccupied,
                                        std::map<SpinComponent, std::vector<int>> listDoublesVirtual, bool addSame) const {
  return derived().generateAllConnectedImpl(det, listDoublesOccupied, listDoublesVirtual, addSame);
}

template<class Crtp>
std::vector<ElectronicDeterminant>
Hamiltonian<Crtp>::generateAllConnectedImpl(const ElectronicDeterminant& det,
                                            std::map<SpinComponent, std::vector<int>> listDoublesOccupied,
                                            std::map<SpinComponent, std::vector<int>> listDoublesVirtual, bool addSame) const {
  auto singles = generateSingleExcitations(det, addSame);
  auto doubles = generateDoubleExcitations(det, listDoublesOccupied, listDoublesVirtual, false);
  std::vector<ElectronicDeterminant> ret;
  ret.reserve(singles.size() + doubles.size());
  ret.insert(ret.end(), std::make_move_iterator(singles.begin()), std::make_move_iterator(singles.end()));
  ret.insert(ret.end(), std::make_move_iterator(doubles.begin()), std::make_move_iterator(doubles.end()));
  return ret;
}

template<class Crtp>
std::vector<ElectronicDeterminant> Hamiltonian<Crtp>::generateSingleExcitations(const ElectronicDeterminant& det,
                                                                                bool addSame) const {
  std::vector<ElectronicDeterminant> ret;
  std::map<SpinComponent, std::vector<int>> listOccupied, listVirtual;
  listOccupied[SpinComponent::Alpha] = det.getAllOccupied(SpinComponent::Alpha);
  listOccupied[SpinComponent::Beta] = det.getAllOccupied(SpinComponent::Beta);
  listVirtual[SpinComponent::Alpha] = det.getAllVirtual(SpinComponent::Alpha);
  listVirtual[SpinComponent::Beta] = det.getAllVirtual(SpinComponent::Beta);
  int nAlphaOcc = listOccupied[SpinComponent::Alpha].size();
  int nAlphaVir = listVirtual[SpinComponent::Alpha].size();
  int nBetaOcc = listOccupied[SpinComponent::Beta].size();
  int nBetaVir = listVirtual[SpinComponent::Beta].size();
  int totalNumber = nAlphaOcc * nAlphaVir + nBetaOcc * nBetaVir;
  ret.reserve(totalNumber + 1);
  if (addSame) {
    ret.push_back(det);
  }
  // Singles excitations
  // number: AlphaOcc*AlphcVir + BetaOcc*BetaVir
  // NO HBCI here, otherwise CIS for example would make problems
  // (0 coupling btw HF and singly excited dets)
  for (auto spin : {SpinComponent::Alpha, SpinComponent::Beta}) {
    for (auto&& iOcc : listOccupied[spin]) {
      for (auto&& iVir : listVirtual[spin]) {
        Excitation excitation;
        excitation.createExcitation(iOcc, iVir, spin);
        auto tmp = det;
        tmp.applyExcitation({iOcc, iVir, spin});
        ret.emplace_back(std::move(tmp));
      }
    }
  }
  return ret;
}

template<class Crtp>
std::vector<ElectronicDeterminant>
Hamiltonian<Crtp>::generateDoubleExcitations(const ElectronicDeterminant& det,
                                             std::map<SpinComponent, std::vector<int>> listOccupied,
                                             std::map<SpinComponent, std::vector<int>> listVirtual, bool addSame) const {
  std::vector<ElectronicDeterminant> ret;

  int nAlphaOcc = listOccupied[SpinComponent::Alpha].size();
  int nAlphaVir = listVirtual[SpinComponent::Alpha].size();
  int nBetaOcc = listOccupied[SpinComponent::Beta].size();
  int nBetaVir = listVirtual[SpinComponent::Beta].size();
  int totalNumber = nAlphaOcc * nAlphaVir * ((nAlphaOcc - 1) * (nAlphaVir - 1) + nBetaOcc * nBetaVir) +
                    nBetaOcc * nBetaVir * ((nBetaOcc - 1) * (nBetaVir - 1));
  ret.reserve(totalNumber + 1);
  if (addSame) {
    ret.push_back(det);
  }
  // Double excitations alpha-alpa/beta-beta
  // number: nAlphaOcc * nAlphaVir * (nAlphaOcc - 1) * (nAlphaVir-1) +
  //         |--------------------|  |-----------------------------|
  //           first excitation             second excitation
  //         |-----------------|        (1 particle and hole less)
  //         |-----------------|   |----------------------------|
  //         nBetaOcc * nBetaVir * (nBetaOcc - 1) * (nBetaVir-1)
  for (auto spin : {SpinComponent::Alpha, SpinComponent::Beta}) {
    for (auto&& iOcc1 : listOccupied[spin]) {
      for (auto&& iOcc2 : listOccupied[spin]) {
        // Select "manually" only one of the two symmetric terms
        if (iOcc1 < iOcc2) {
          for (auto&& iVir1 : listVirtual[spin]) {
            for (auto&& iVir2 : listVirtual[spin]) {
              if (iVir1 < iVir2) {
                Excitation excitations;
                excitations.createExcitation(iOcc1, iVir1, spin);
                excitations.createExcitation(iOcc2, iVir2, spin);
                auto tmp = det;
                tmp.applyExcitation(excitations.getSingleExcitation(0));
                tmp.applyExcitation(excitations.getSingleExcitation(1));
                ret.emplace_back(std::move(tmp));
              }
            }
          }
        }
      }
    }
  }
  // Alpha-beta excitation
  // number: nAlphaOcc * nAlphaVir * nBetaOcc * nBetaVir
  for (auto&& iAlphaOcc : listOccupied[SpinComponent::Alpha]) {
    for (auto&& iBetaOcc : listOccupied[SpinComponent::Beta]) {
      OrderedIndex sortedIndex = {iAlphaOcc, iBetaOcc};
      SpinComponent firstOrbital = SpinComponent::Alpha;
      SpinComponent secondOrbital = SpinComponent::Beta;
      if (sortedIndex[0] > sortedIndex[1]) {
        std::swap(sortedIndex[0], sortedIndex[1]);
        std::swap(firstOrbital, secondOrbital);
      }
      // TODO: Was this just for ergodicity of selected CI?
      // if (hbciContainerOppositeSpin_.find(sortedIndex) != hbciContainerOppositeSpin_.end()) {
      // for (const auto& el : hbciContainerOppositeSpin_.at(sortedIndex)) {
      // if (!det.isOccupied(firstOrbital, el.second[0]) && !det.isOccupied(secondOrbital, el.second[1])) {
      for (auto&& iAlphaVir : listVirtual[firstOrbital]) {
        for (auto&& iBetaVir : listVirtual[secondOrbital]) {
          auto tmp = det;
          tmp.applyExcitation({sortedIndex[0], iAlphaVir, firstOrbital});
          tmp.applyExcitation({sortedIndex[1], iBetaVir, secondOrbital});
          ret.emplace_back(std::move(tmp));
          // if (iAlphaOcc == iBetaOcc && sortedIndexV[0] != sortedIndexV[1]) {
          //  auto tmp = det;
          //  tmp.applyExcitation({sortedIndex[0], sortedIndexV[0], firstOrbital});
          //  tmp.applyExcitation({sortedIndex[1], sortedIndexV[1], secondOrbital});
          //  ret.emplace_back(std::move(tmp));
          //}
          //}
        }
      }
    }
  }
  return ret;
}

template<class Crtp>
std::vector<ElectronicDeterminant> Hamiltonian<Crtp>::generateDoubleExcitations(const ElectronicDeterminant& det,
                                                                                bool addSame) const {
  std::map<SpinComponent, std::vector<int>> listOccupied{{SpinComponent::Alpha, det.getAllOccupied(SpinComponent::Alpha)},
                                                         {SpinComponent::Beta, det.getAllOccupied(SpinComponent::Beta)}};
  std::map<SpinComponent, std::vector<int>> listVirtual{{SpinComponent::Alpha, det.getAllVirtual(SpinComponent::Alpha)},
                                                        {SpinComponent::Beta, det.getAllVirtual(SpinComponent::Beta)}};
  return generateDoubleExcitations(det, listOccupied, listVirtual, addSame);
}

template<class Crtp>
int Hamiltonian<Crtp>::getNumberOfConnectedDeterminants(const ElectronicDeterminant& det, bool addSame) const {
  int nAlphaOcc = det.getAllOccupied(SpinComponent::Alpha).size();
  int nBetaOcc = det.getAllOccupied(SpinComponent::Beta).size();
  int nAlphaVir = det.getAllVirtual(SpinComponent::Alpha).size();
  int nBetaVir = det.getAllVirtual(SpinComponent::Beta).size();
  int totalNumber = nAlphaOcc * nAlphaVir * ((nAlphaOcc - 1) * (nAlphaVir - 1) + 1 + nBetaOcc * nBetaVir) +
                    nBetaOcc * nBetaVir * ((nBetaOcc - 1) * (nBetaVir - 1) + 1);
  if (addSame)
    totalNumber += 1;
  return totalNumber;
}

template<class Crtp>
void Hamiltonian<Crtp>::setHbciThreshold(double threshold) {
  hbciThreshold_ = threshold;
  for (auto& hbciElement : hbciContainerSameSpin_)
    hbciElement.second.setThreshold(hbciThreshold_);
}

template<class Crtp>
double Hamiltonian<Crtp>::matrixElementSameDeterminant(const std::vector<int>& lst_occupied_alpha,
                                                       const std::vector<int>& lst_occupied_beta) const {
  double ret = 0.;
  // One-body contribution
  for (const auto& occ_alpha : lst_occupied_alpha) {
    ret += oneElectronMatrix_(occ_alpha, occ_alpha);
  }
  for (const auto& occ_beta : lst_occupied_beta) {
    ret += oneElectronMatrix_(occ_beta, occ_beta);
  }
  // Two-bodies contribution
  // Alpha-Alpha
  for (const auto& occ_alpha_i : lst_occupied_alpha) {
    for (const auto& occ_alpha_j : lst_occupied_alpha) {
      if (occ_alpha_i != occ_alpha_j) {
        auto key1 = bringToNormalForm({occ_alpha_i, occ_alpha_j, occ_alpha_i, occ_alpha_j});
        auto key2 = bringToNormalForm({occ_alpha_i, occ_alpha_j, occ_alpha_j, occ_alpha_i});
        if (twoElectronMatrix_.find(key1) != twoElectronMatrix_.end())
          ret += twoElectronMatrix_.at(key1) / 2.;
        if (twoElectronMatrix_.find(key2) != twoElectronMatrix_.end())
          ret -= twoElectronMatrix_.at(key2) / 2.;
      }
    }
  }
  //  Beta-Beta
  for (const auto& occ_beta_i : lst_occupied_beta) {
    for (const auto& occ_beta_j : lst_occupied_beta) {
      if (occ_beta_i != occ_beta_j) {
        auto key1 = bringToNormalForm({occ_beta_i, occ_beta_j, occ_beta_i, occ_beta_j});
        auto key2 = bringToNormalForm({occ_beta_i, occ_beta_j, occ_beta_j, occ_beta_i});
        if (twoElectronMatrix_.find(key1) != twoElectronMatrix_.end())
          ret += twoElectronMatrix_.at(key1) / 2.;
        if (twoElectronMatrix_.find(key2) != twoElectronMatrix_.end())
          ret -= twoElectronMatrix_.at(key2) / 2.;
      }
    }
  }
  // Alpha-Beta
  for (const auto& occ_alpha_i : lst_occupied_alpha) {
    for (const auto& occ_beta_j : lst_occupied_beta) {
      auto&& key = bringToNormalForm({occ_alpha_i, occ_beta_j, occ_alpha_i, occ_beta_j});
      if (twoElectronMatrix_.find(key) != twoElectronMatrix_.end())
        ret += twoElectronMatrix_.at(key);
    }
  }
  return ret;
}

template<class Crtp>
double Hamiltonian<Crtp>::matrixElementSingleExcitation(const SingleExcitation& exc1, const std::vector<int>& lst_occupied_alpha,
                                                        const std::vector<int>& lst_occupied_beta) const {
  // Calculates auxiliary quantities
  double ret = 0.;
  const auto& lst_samespin = (exc1.orbitalSpin == SpinComponent::Alpha) ? lst_occupied_alpha : lst_occupied_beta;
  const auto& lst_otherspin = (exc1.orbitalSpin == SpinComponent::Alpha) ? lst_occupied_beta : lst_occupied_alpha;
  auto coeff = calculateExcitationSign(lst_occupied_alpha, lst_occupied_beta, exc1);
  // One-body term
  ret = oneElectronMatrix_(std::min(exc1.virtualOrbital, exc1.occupiedOrbital),
                           std::max(exc1.virtualOrbital, exc1.occupiedOrbital));
  // Two-body contribution
  for (const auto& i_same_spin : lst_samespin) {
    if (i_same_spin != exc1.occupiedOrbital && i_same_spin != exc1.virtualOrbital) {
      auto key1 = bringToNormalForm({exc1.virtualOrbital, i_same_spin, exc1.occupiedOrbital, i_same_spin});
      auto key2 = bringToNormalForm({exc1.virtualOrbital, i_same_spin, i_same_spin, exc1.occupiedOrbital});
      if (twoElectronMatrix_.find(key1) != twoElectronMatrix_.end())
        ret += twoElectronMatrix_.at(key1);
      if (twoElectronMatrix_.find(key2) != twoElectronMatrix_.end())
        ret -= twoElectronMatrix_.at(key2);
    }
  }
  for (const auto& i_other_spin : lst_otherspin) {
    auto key = bringToNormalForm({exc1.virtualOrbital, i_other_spin, exc1.occupiedOrbital, i_other_spin});
    if (twoElectronMatrix_.find(key) != twoElectronMatrix_.end())
      ret += twoElectronMatrix_.at(key);
  }
  return coeff * ret;
}

template<class Crtp>
double Hamiltonian<Crtp>::matrixElementDoubleExcitation(const SingleExcitation& exc1, const SingleExcitation& exc2,
                                                        const std::vector<int>& lst_occupied_alpha,
                                                        const std::vector<int>& lst_occupied_beta) const {
  // Also here we have to use the JW transformation to keep track of the correct sign for the excitation.
  // For the moment, we implement the mapping here, but in the future this should be moved maybe inside the
  // determinant class.
  double ret = 0.;
  int sign = calculateExcitationSign(lst_occupied_alpha, lst_occupied_beta, exc1, exc2);
  if (exc1.orbitalSpin == exc2.orbitalSpin) {
    auto key1 = bringToNormalForm({exc1.virtualOrbital, exc2.virtualOrbital, exc1.occupiedOrbital, exc2.occupiedOrbital});
    auto key2 = bringToNormalForm({exc1.virtualOrbital, exc2.virtualOrbital, exc2.occupiedOrbital, exc1.occupiedOrbital});
    if (twoElectronMatrix_.find(key1) != twoElectronMatrix_.end())
      ret += twoElectronMatrix_.at(key1);
    if (twoElectronMatrix_.find(key2) != twoElectronMatrix_.end())
      ret -= twoElectronMatrix_.at(key2);
  }
  else {
    auto&& key = bringToNormalForm({exc1.virtualOrbital, exc2.virtualOrbital, exc1.occupiedOrbital, exc2.occupiedOrbital});
    if (twoElectronMatrix_.find(key) != twoElectronMatrix_.end())
      ret += twoElectronMatrix_.at(key);
  }
  return ret * sign;
}

template<class Crtp>
void Hamiltonian<Crtp>::addTermToHBCIContainer(int iOcc1, int iOcc2, int iVir1, int iVir2, double value) {
  OrderedIndex ij = {iOcc1, iOcc2};
  OrderedIndex kl = {iVir1, iVir2};
  if (iOcc1 > iOcc2) {
    std::swap(ij[0], ij[1]);
    std::swap(kl[0], kl[1]);
  }
  // Sets the threshold for HBCI
  if (hbciContainerSameSpin_.find(ij) == hbciContainerSameSpin_.end())
    hbciContainerSameSpin_[ij].setThreshold(hbciThreshold_);
  if (hbciContainerOppositeSpin_.find(ij) == hbciContainerOppositeSpin_.end())
    hbciContainerOppositeSpin_[ij].setThreshold(hbciThreshold_);
  // Loads the opposite-spin container
  hbciContainerOppositeSpin_[ij].insert(value, kl);
  // Loads the same-spin container
  if (iOcc1 != iOcc2 && iVir1 != iVir2) {
    if (twoElectronMatrix_.find(bringToNormalForm({iOcc1, iOcc2, iVir2, iVir1})) != twoElectronMatrix_.end())
      value -= twoElectronMatrix_.at(bringToNormalForm({iOcc1, iOcc2, iVir2, iVir1}));
    hbciContainerSameSpin_[ij].insert(value, kl);
  }
}

template<class Crtp>
typename Hamiltonian<Crtp>::OneBodyIdentifier
Hamiltonian<Crtp>::bringToNormalForm(const OneBodyIdentifier& input_ob_identifier) const {
  if (input_ob_identifier.first > input_ob_identifier.second) {
    return std::make_pair(input_ob_identifier.second, input_ob_identifier.first);
  }
  return std::make_pair(input_ob_identifier.first, input_ob_identifier.second);
}

template<class Crtp>
typename Hamiltonian<Crtp>::TwoBodyIdentifier
Hamiltonian<Crtp>::bringToNormalForm(const TwoBodyIdentifier& input_ob_identifier) const {
  return derived().bringToNormalFormImpl(input_ob_identifier);
  // auto sorted_array = input_ob_identifier;
  // if (sorted_array[0] > sorted_array[2]) {
  //   std::swap(sorted_array[0], sorted_array[2]);
  // }
  // if (sorted_array[1] > sorted_array[3]) {
  //   std::swap(sorted_array[1], sorted_array[3]);
  // }
  // if (sorted_array[0] > sorted_array[1]) {
  //   std::swap(sorted_array[0], sorted_array[1]);
  //   std::swap(sorted_array[2], sorted_array[3]);
  // }
  // else if (sorted_array[0] == sorted_array[1] && sorted_array[2] > sorted_array[3]) {
  //   std::swap(sorted_array[2], sorted_array[3]);
  // }
  // return sorted_array;
}

template<class Crtp>
const auto& Hamiltonian<Crtp>::derived() const {
  return static_cast<const Crtp&>(*this);
}

template<class Crtp>
auto& Hamiltonian<Crtp>::derived() {
  return static_cast<Crtp&>(*this);
}
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
