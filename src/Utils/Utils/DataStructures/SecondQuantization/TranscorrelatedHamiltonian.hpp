/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

namespace Scine {
namespace Utils {
namespace SecondQuantization {

std::vector<ElectronicDeterminant> TranscorrelatedHamiltonian::generateAllConnectedImpl(const ElectronicDeterminant& det, bool addSame) const {
  auto singles = Base::generateSingleExcitations(det, addSame);
  auto doubles = Base::generateDoubleExcitations(det, false);
//   auto tripples = generateTrippleExcitations(det, false);
  std::vector<ElectronicDeterminant> ret;
  ret.reserve(singles.size() + doubles.size());
  ret.insert(ret.end(), std::make_move_iterator(singles.begin()), std::make_move_iterator(singles.end()));
  ret.insert(ret.end(), std::make_move_iterator(doubles.begin()), std::make_move_iterator(doubles.end()));
//   ret.insert(ret.end(), std::make_move_iterator(tripples.begin()), std::make_move_iterator(tripples.end()));
  return ret;
}

// std::vector<ElectronicDeterminant> TranscorrelatedHamiltonian::generateTrippleExcitations(const ElectronicDeterminant& det, bool addSame) const {
//     std::vector<ElectronicDeterminant> ret;
//     std::map<SpinComponent, std::vector<int>> listOccupied{{SpinComponent::Alpha, det.getAllOccupied(SpinComponent::Alpha)},
//                                                            {SpinComponent::Beta, det.getAllOccupied(SpinComponent::Beta)}};
//     std::map<SpinComponent, std::vector<int>> listVirtual{{SpinComponent::Alpha, det.getAllVirtual(SpinComponent::Alpha)},
//                                                           {SpinComponent::Beta, det.getAllVirtual(SpinComponent::Beta)}};
//     int nAlphaOcc = listOccupied[SpinComponent::Alpha].size();
//     int nAlphaVir = listVirtual[SpinComponent::Alpha].size();
//     int nBetaOcc = listOccupied[SpinComponent::Beta].size();
//     int nBetaVir = listVirtual[SpinComponent::Beta].size();
//     int totalNumber = nAlphaOcc * nAlphaVir * ((nAlphaOcc - 1) * (nAlphaVir - 1) + nBetaOcc * nBetaVir) +
//                       nBetaOcc * nBetaVir * ((nBetaOcc - 1) * (nBetaVir - 1));
//     ret.reserve(totalNumber + 1);
//     if (addSame) {
//         ret.push_back(det);
//     }
//     // Double excitations alpha-alpa-alpha/beta-beta-beta
//     // number: nAlphaOcc * nAlphaVir * (nAlphaOcc - 1) * (nAlphaVir-1) * (nAlphaOcc - 2) * (nAlphaVir-2) +
//     //         |--------------------|  |-----------------------------|   |-----------------------------|
//     //           first excitation             second excitation                 third excitation
//     //         |------------------|       (1 particle and hole less)      (2 particles and holes less )
//     //         |------------------|  |---------------------------|   |---------------------------|
//     //         nBetaOcc * nBetaVir * (nBetaOcc - 1) * (nBetaVir-1) * (nBetaOcc - 2) * (nBetaVir-2)
//     for (auto spin : {SpinComponent::Alpha, SpinComponent::Beta}) {
//         for (auto&& iOcc1 : listOccupied[spin]) {
//             for (auto&& iOcc2 : listOccupied[spin]) {
//                 for (auto&& iOcc3 : listOccupied[spin]) {
//                     // Select "manually" only one of the two symmetric terms
//                     if (iOcc1 < iOcc2 && iOcc1 < iOcc3 && iOcc2 < iOcc3) {
//                         if (hbciContainerSameSpin_.find({iOcc1, iOcc2, iOcc3}) != hbciContainerSameSpin_.end()) {
//                             for (const auto& el : hbciContainerSameSpin_.at({iOcc1, iOcc2, iOcc3})) {
//                                 if (!det.isOccupied(spin, el.second[0]) && !det.isOccupied(spin, el.second[1]) && !det.isOccupied(spin, el.second[1])
//                                 && el.second[0] < el.second[1] && el.second[1] < el.second[2]) {
//                                     Excitation excitations;
//                                     excitations.createExcitation(iOcc1, el.second[0], spin);
//                                     excitations.createExcitation(iOcc2, el.second[1], spin);
//                                     excitations.createExcitation(iOcc2, el.second[2], spin);
//                                     auto tmp = det;
//                                     tmp.applyExcitation(excitations.getSingleExcitation(0));
//                                     tmp.applyExcitation(excitations.getSingleExcitation(1));
//                                     tmp.applyExcitation(excitations.getSingleExcitation(2));
//                                     ret.emplace_back(std::move(tmp));
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     // Alpha-alpha-beta/Alpha-beta-beta excitation
//     // number: nAlphaOcc * nAlphaVir * nBetaOcc * nBetaVir * (nBetaOcc-1) * (nBetaVir-1)
//     //       + nAlphaOcc * nAlphaVir * nBetaOcc * nBetaVir * (nAlphaOcc-1) * (nAlphaOcc-1)
//     for (auto&& iAlphaOcc1 : listOccupied[SpinComponent::Alpha]) {
//         for (auto&& iBetaOcc1 : listOccupied[SpinComponent::Beta]) {
//             for (auto&& iBetaOcc2 : listOccupied[SpinComponent::Beta]) {

//                 OrderedIndex sortedIndex = {iAlphaOcc1, iBetaOcc1, iBetaOcc2};
//                 SpinComponent firstOrbital = SpinComponent::Alpha;
//                 SpinComponent secondOrbital = SpinComponent::Beta;
//                 SpinComponent thirdOrbital = SpinComponent::Beta;
//                 // bring everything in order 2 > 1 > 0
//                 // 0 > 1 > 2
//                 if (sortedIndex[0] > sortedIndex[1] &&
//                     sortedIndex[1] > sortedIndex[2] &&
//                     sortedIndex[0] > sortedIndex[2]) {
//                     std::swap(sortedIndex[0], sortedIndex[2]);
//                     std::swap(firstOrbital, thirdOrbital);
//                 }
//                 // 1 > 2 > 0
//                 else if (sortedIndex[1] > sortedIndex[0] &&
//                         sortedIndex[1] > sortedIndex[2] &&
//                         sortedIndex[2] > sortedIndex[0]) {
//                     std::swap(sortedIndex[1], sortedIndex[2]);
//                     std::swap(secondOrbital, thirdOrbital);
//                 }
//                 // 2 > 0 > 1
//                 else if (sortedIndex[2] > sortedIndex[0] &&
//                         sortedIndex[2] > sortedIndex[1] &&
//                         sortedIndex[0] > sortedIndex[1]) {
//                     std::swap(sortedIndex[0], sortedIndex[1]);
//                     std::swap(firstOrbital, secondOrbital);
//                 }
//                 // 1 > 0 > 2
//                 else if (sortedIndex[1] > sortedIndex[0] &&
//                         sortedIndex[1] > sortedIndex[2] &&
//                         sortedIndex[0] > sortedIndex[2]) {
//                     std::swap(sortedIndex[0], sortedIndex[2]);
//                     std::swap(firstOrbital, thirdOrbital);
//                     std::swap(sortedIndex[1], sortedIndex[2]);
//                     std::swap(secondOrbital, thirdOrbital);
//                 }
//                 // 0 > 2 > 1
//                 else if (sortedIndex[0] > sortedIndex[2] &&
//                         sortedIndex[0] > sortedIndex[1] &&
//                         sortedIndex[2] > sortedIndex[1]) {
//                     std::swap(sortedIndex[0], sortedIndex[1]);
//                     std::swap(firstOrbital, secondOrbital);
//                     std::swap(sortedIndex[1], sortedIndex[2]);
//                     std::swap(secondOrbital, thirdOrbital);
//                 }
//                 // TODO: Was this just for ergodicity of selected CI?
//                 // if (hbciContainerOppositeSpin_.find(sortedIndex) != hbciContainerOppositeSpin_.end()) {
//                 // for (const auto& el : hbciContainerOppositeSpin_.at(sortedIndex)) {
//                 // if (!det.isOccupied(firstOrbital, el.second[0]) && !det.isOccupied(secondOrbital, el.second[1])) {
//                 for (auto&& iAlphaVir1 : listVirtual[firstOrbital]) {
//                     for (auto&& iBetaVir1 : listVirtual[secondOrbital]) {
//                         for (auto&& iBetaVir2 : listVirtual[thirdOrbital]) {
//                             auto tmp = det;
//                             tmp.applyExcitation({sortedIndex[0], iAlphaVir1, firstOrbital});
//                             tmp.applyExcitation({sortedIndex[1], iBetaVir1, secondOrbital});
//                             tmp.applyExcitation({sortedIndex[2], iBetaVir2, thirdOrbital});
//                             ret.emplace_back(std::move(tmp));
//                             // if (iAlphaOcc == iBetaOcc && sortedIndexV[0] != sortedIndexV[1]) {
//                             //  auto tmp = det;
//                             //  tmp.applyExcitation({sortedIndex[0], sortedIndexV[0], firstOrbital});
//                             //  tmp.applyExcitation({sortedIndex[1], sortedIndexV[1], secondOrbital});
//                             //  ret.emplace_back(std::move(tmp));
//                             //}
//                             //}
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     return ret;
// }

}
}
}