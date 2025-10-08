/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SecondQuantization/CIDeterminantsCreator.h>
#include <Utils/DataStructures/SecondQuantization/ElectronicDeterminant.h>
#include <Utils/DataStructures/SecondQuantization/EriUtilities.h>
#include <Utils/DataStructures/SecondQuantization/HermitianHamiltonian.h>
#include <Utils/DataStructures/SecondQuantization/NonHermitianHamiltonian.h>
#include <Utils/DataStructures/SecondQuantization/TranscorrelatedHamiltonian.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <gmock/gmock.h>
#include <vector>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

using namespace testing;

/**
 * @test SQUtilsTest @file SecondQuantizationTest.cpp
 * Simple testing of excitation of a determinant.
 */

class SQUtilsTest : public Test {
 public:
  // Bare determinants
  std::vector<int> onv_2e5o, onv_2e6o, onv_4e20o, onv_1e4o, onv_10e1000o;
  // Hamiltonians
  HermitianHamiltonian hydrogenHamiltonianOneBody{4};
  HermitianHamiltonian nitrogenHamiltonianOneBody{6};
  HermitianHamiltonian hydrogenHamiltonianTwoBody{4};
  HermitianHamiltonian nitrogenHamiltonianTwoBody{6};
  NonHermitianHamiltonian<> heTranscorrelatedHamiltonian{5};
  // TranscorrelatedHamiltonian SomeTranscorrelatedHamiltonian{5}; //TODO

  // Full-CI wave functions (often then loaded with DMRG-based wave functions)
  std::unordered_map<ElectronicDeterminant, double, ElectronicDeterminantHash> H2CICoeff, N2CICoeff, FHCICoeff;

 private:
  void SetUp() final {
    // ONVs
    onv_1e4o = {1, 0, 0, 0};
    onv_2e5o = {1, 1, 0, 0, 0};
    onv_2e6o = {1, 1, 0, 0, 0, 0};
    onv_4e20o = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    onv_10e1000o = std::vector<int>(10000, 0);
    for (int i = 0; i < 100; i++)
      onv_10e1000o[i] = 1;
    // CI Wave function for H2
    H2CICoeff[ElectronicDeterminant({1, 0, 0, 0}, {1, 0, 0, 0})] = -0.9927422191;
    H2CICoeff[ElectronicDeterminant({0, 0, 1, 0}, {0, 0, 1, 0})] = 0.0767633238;
    H2CICoeff[ElectronicDeterminant({0, 1, 0, 0}, {0, 1, 0, 0})] = 0.0504776815;
    H2CICoeff[ElectronicDeterminant({0, 0, 1, 0}, {0, 0, 0, 1})] = -0.0454083555;
    H2CICoeff[ElectronicDeterminant({0, 0, 0, 1}, {0, 0, 1, 0})] = 0.0454083555;
    H2CICoeff[ElectronicDeterminant({0, 0, 0, 1}, {0, 0, 0, 1})] = 0.0428155555;
    H2CICoeff[ElectronicDeterminant({1, 0, 0, 0}, {0, 1, 0, 0})] = 0.0057128337;
    H2CICoeff[ElectronicDeterminant({0, 1, 0, 0}, {1, 0, 0, 0})] = -0.0057128337;
    // CI Wave function for FH
    // FHCICoeff[ElectronicDeterminant({1, 1, 0, 0}, {1, 0, 0, 0})] = 1.;
    FHCICoeff[ElectronicDeterminant({1, 0, 0, 0}, {1, 0, 0, 0})] = -0.9694575537;
    FHCICoeff[ElectronicDeterminant({0, 1, 0, 0}, {0, 1, 0, 0})] = 0.1648916292;
    FHCICoeff[ElectronicDeterminant({0, 0, 1, 0}, {0, 0, 1, 0})] = 0.1648916292;
    FHCICoeff[ElectronicDeterminant({0, 0, 0, 1}, {0, 0, 0, 1})] = 0.07598389752;
    FHCICoeff[ElectronicDeterminant({1, 0, 0, 0}, {1, 0, 0, 0})] = -0.9694575537;
    FHCICoeff[ElectronicDeterminant({0, 1, 0, 0}, {0, 1, 0, 0})] = 0.1648916292;
    FHCICoeff[ElectronicDeterminant({0, 0, 1, 0}, {0, 0, 1, 0})] = 0.1648916292;
    FHCICoeff[ElectronicDeterminant({0, 0, 0, 1}, {0, 0, 0, 1})] = 0.07598389752;
    // CI wave function for N2 (obtained with DMRG)
    N2CICoeff[ElectronicDeterminant("222000")] = -0.9671507429;
    N2CICoeff[ElectronicDeterminant("202020")] = 0.1214684848;
    N2CICoeff[ElectronicDeterminant("220002")] = 0.1214682754;
    N2CICoeff[ElectronicDeterminant("2ud0du")] = 0.0771249498;
    N2CICoeff[ElectronicDeterminant("2du0ud")] = 0.0771249497;
    N2CICoeff[ElectronicDeterminant("2uu0dd")] = -0.0490935732;
    N2CICoeff[ElectronicDeterminant("2dd0uu")] = -0.0490935731;
    N2CICoeff[ElectronicDeterminant("022200")] = 0.0463591830;
    N2CICoeff[ElectronicDeterminant("u2dd0u")] = 0.0428244585;
    N2CICoeff[ElectronicDeterminant("d2uu0d")] = 0.0428244585;
    N2CICoeff[ElectronicDeterminant("ud2du0")] = 0.0428244491;
    N2CICoeff[ElectronicDeterminant("du2ud0")] = 0.0428244490;
    N2CICoeff[ElectronicDeterminant("200022")] = -0.0318556912;
    N2CICoeff[ElectronicDeterminant("2ud0ud")] = -0.0280313766;
    N2CICoeff[ElectronicDeterminant("2du0du")] = -0.0280313766;
    N2CICoeff[ElectronicDeterminant("d2du0u")] = -0.0279399224;
    N2CICoeff[ElectronicDeterminant("u2ud0d")] = -0.0279399224;
    N2CICoeff[ElectronicDeterminant("dd2uu0")] = -0.0279399145;
    N2CICoeff[ElectronicDeterminant("uu2dd0")] = -0.0279399145;
    N2CICoeff[ElectronicDeterminant("220020")] = 0.0163120655;
    N2CICoeff[ElectronicDeterminant("202002")] = 0.0163120426;
    N2CICoeff[ElectronicDeterminant("u2du0d")] = -0.0148845361;
    N2CICoeff[ElectronicDeterminant("d2ud0u")] = -0.0148845361;
    N2CICoeff[ElectronicDeterminant("du2du0")] = -0.0148845345;
    N2CICoeff[ElectronicDeterminant("ud2ud0")] = -0.0148845345;
    N2CICoeff[ElectronicDeterminant("u0dd2u")] = -0.0123087568;
    N2CICoeff[ElectronicDeterminant("d0uu2d")] = -0.0123087568;
    N2CICoeff[ElectronicDeterminant("du0ud2")] = -0.0123087443;
    N2CICoeff[ElectronicDeterminant("ud0du2")] = -0.0123087443;
    N2CICoeff[ElectronicDeterminant("202200")] = 0.0106833942;
    N2CICoeff[ElectronicDeterminant("220200")] = 0.0106833938;
    N2CICoeff[ElectronicDeterminant("002220")] = -0.0099785117;
    N2CICoeff[ElectronicDeterminant("020202")] = -0.0099785021;
    N2CICoeff[ElectronicDeterminant("d0du2u")] = 0.0097324477;
    N2CICoeff[ElectronicDeterminant("u0ud2d")] = 0.0097324477;
    N2CICoeff[ElectronicDeterminant("uu0dd2")] = 0.0097324386;
    N2CICoeff[ElectronicDeterminant("dd0uu2")] = 0.0097324386;
    N2CICoeff[ElectronicDeterminant("022002")] = 0.0084695339;
    N2CICoeff[ElectronicDeterminant("022020")] = 0.0084695154;
    N2CICoeff[ElectronicDeterminant("0ud2du")] = -0.0071568181;
    N2CICoeff[ElectronicDeterminant("0du2ud")] = -0.0071568181;
    N2CICoeff[ElectronicDeterminant("0dd2uu")] = 0.0053517701;
    N2CICoeff[ElectronicDeterminant("0uu2dd")] = 0.0053517701;
    N2CICoeff[ElectronicDeterminant("000222")] = 0.0046262302;
    N2CICoeff[ElectronicDeterminant("d0ud2u")] = 0.0025763091;
    N2CICoeff[ElectronicDeterminant("u0du2d")] = 0.0025763091;
    N2CICoeff[ElectronicDeterminant("ud0ud2")] = 0.0025763056;
    N2CICoeff[ElectronicDeterminant("du0du2")] = 0.0025763056;
    N2CICoeff[ElectronicDeterminant("0ud2ud")] = 0.0018050480;
    N2CICoeff[ElectronicDeterminant("0du2du")] = 0.0018050480;
    N2CICoeff[ElectronicDeterminant("200220")] = -0.0016533117;
    N2CICoeff[ElectronicDeterminant("200202")] = -0.0016533089;
    N2CICoeff[ElectronicDeterminant("002022")] = -0.0013778651;
    N2CICoeff[ElectronicDeterminant("020022")] = -0.0013778609;
    N2CICoeff[ElectronicDeterminant("020220")] = -0.0010166412;
    N2CICoeff[ElectronicDeterminant("002202")] = -0.0010166401;
    // One-boddy H2 Hamiltonian
    hydrogenHamiltonianOneBody.addTerm(0, 0, -1.245095342160);
    hydrogenHamiltonianOneBody.addTerm(1, 0, 0.167073344018);
    hydrogenHamiltonianOneBody.addTerm(1, 1, -0.178953017984);
    hydrogenHamiltonianOneBody.addTerm(2, 2, -0.549284223834);
    hydrogenHamiltonianOneBody.addTerm(3, 2, 0.207313806957);
    hydrogenHamiltonianOneBody.addTerm(3, 3, 0.214479164774);
    // One-Body N2 Hamiltonian (calculated with cc-pVDZ at a bond distance of 2.074 Bohr)
    nitrogenHamiltonianOneBody.addTerm(0, 0, -4.32877350104);
    nitrogenHamiltonianOneBody.addTerm(1, 1, -3.44821176858);
    nitrogenHamiltonianOneBody.addTerm(2, 2, -3.44821159660);
    nitrogenHamiltonianOneBody.addTerm(3, 3, -2.73682333759);
    nitrogenHamiltonianOneBody.addTerm(4, 4, -3.10648615388);
    nitrogenHamiltonianOneBody.addTerm(5, 5, -3.10648710673);
    // Two-Body H2 Hamiltonian
    Eigen::MatrixXd twoBodyH2 = Eigen::MatrixXd::Zero(16, 16);
    std::unordered_map<std::array<int, 4>, double, boost::hash<std::array<int, 4>>> twoBodyH2Map;
    twoBodyH2Map[{0, 0, 0, 0}] = 0.649702723853;
    twoBodyH2Map[{0, 1, 0, 0}] = -0.167073344018;
    twoBodyH2Map[{1, 1, 0, 0}] = 0.109300885764;
    twoBodyH2Map[{1, 1, 1, 0}] = -0.119851267485;
    twoBodyH2Map[{0, 1, 0, 1}] = 0.531826350603;
    twoBodyH2Map[{1, 1, 1, 1}] = 0.463674286178;
    twoBodyH2Map[{2, 2, 0, 0}] = 0.801465148695E-01;
    twoBodyH2Map[{2, 2, 0, 1}] = 0.192573520766E-01;
    twoBodyH2Map[{2, 2, 1, 1}] = 0.359193058681E-01;
    twoBodyH2Map[{0, 2, 0, 2}] = 0.433764500363;
    twoBodyH2Map[{1, 2, 0, 2}] = -0.500848028627E-01;
    twoBodyH2Map[{1, 2, 1, 2}] = 0.381382374889;
    twoBodyH2Map[{2, 2, 2, 2}] = 0.385855813983;
    twoBodyH2Map[{2, 3, 0, 0}] = -0.793764546306E-01;
    twoBodyH2Map[{2, 3, 1, 0}] = 0.218346809057E-01;
    twoBodyH2Map[{3, 3, 0, 0}] = 0.137553175046;
    twoBodyH2Map[{2, 3, 0, 1}] = 0.833226832200E-01;
    twoBodyH2Map[{2, 3, 1, 1}] = -0.270771262413E-02;
    twoBodyH2Map[{3, 3, 0, 1}] = -0.123112450743;
    twoBodyH2Map[{3, 3, 1, 1}] = 0.127594100045;
    twoBodyH2Map[{0, 3, 0, 2}] = -0.143345130794;
    twoBodyH2Map[{1, 3, 0, 2}] = 0.733156940616E-01;
    twoBodyH2Map[{1, 3, 1, 2}] = -0.984145436380E-01;
    twoBodyH2Map[{2, 3, 2, 2}] = -0.548241399472E-01;
    twoBodyH2Map[{3, 3, 2, 2}] = 0.675771906037E-01;
    twoBodyH2Map[{3, 3, 3, 2}] = -0.167748159247;
    twoBodyH2Map[{0, 3, 0, 3}] = 0.662820064425;
    twoBodyH2Map[{1, 3, 0, 3}] = -0.201494810589;
    twoBodyH2Map[{1, 3, 1, 3}] = 0.552219756992;
    twoBodyH2Map[{2, 3, 2, 3}] = 0.442474250678;
    twoBodyH2Map[{3, 3, 3, 3}] = 0.740170356689;
    for (const auto& el : twoBodyH2Map) {
      for (const auto& index :
           getSymmetricIndices<ERISymmetry::eightfold>({el.first[0], el.first[2], el.first[1], el.first[3]})) {
        twoBodyH2(index[0] * 4 + index[1], index[2] * 4 + index[3]) = el.second;
      }
    }
    hydrogenHamiltonianTwoBody.addTwoBodyContribution(twoBodyH2);
    // Two-Body N2 Hamiltonian
    Eigen::MatrixXd twoBodyN2 = Eigen::MatrixXd::Zero(36, 36);
    std::unordered_map<std::array<int, 4>, double, boost::hash<std::array<int, 4>>> twoBodyN2Map;
    twoBodyN2Map[{0, 0, 0, 0}] = 0.866820952529;
    twoBodyN2Map[{1, 1, 0, 0}] = 0.504415557452E-01;
    twoBodyN2Map[{0, 1, 0, 1}] = 0.633798153839;
    twoBodyN2Map[{1, 1, 1, 1}] = 0.582266817195;
    twoBodyN2Map[{2, 2, 0, 0}] = 0.504415560309E-01;
    twoBodyN2Map[{2, 2, 1, 1}] = 0.233330020057E-01;
    twoBodyN2Map[{0, 2, 0, 2}] = 0.633798110584;
    twoBodyN2Map[{1, 2, 1, 2}] = 0.535600779651;
    twoBodyN2Map[{2, 2, 2, 2}] = 0.582266750130;
    twoBodyN2Map[{3, 3, 0, 0}] = 0.195926182697;
    twoBodyN2Map[{3, 3, 1, 1}] = 0.348061889568E-01;
    twoBodyN2Map[{3, 3, 2, 2}] = 0.348061827325E-01;
    twoBodyN2Map[{0, 3, 0, 3}] = 0.842807220119;
    twoBodyN2Map[{1, 3, 1, 3}] = 0.648313915010;
    twoBodyN2Map[{2, 3, 2, 3}] = 0.648313859378;
    twoBodyN2Map[{3, 3, 3, 3}] = 0.879004197322;
    twoBodyN2Map[{3, 4, 1, 0}] = 0.230002415882E-01;
    twoBodyN2Map[{4, 4, 0, 0}] = 0.188825936828E-01;
    twoBodyN2Map[{3, 4, 0, 1}] = 0.149722796368;
    twoBodyN2Map[{4, 4, 1, 1}] = 0.185346091589;
    twoBodyN2Map[{4, 4, 2, 2}] = 0.179848570046E-01;
    twoBodyN2Map[{1, 4, 0, 3}] = 0.449318220227E-01;
    twoBodyN2Map[{4, 4, 3, 3}] = 0.519777084859E-01;
    twoBodyN2Map[{0, 4, 0, 4}] = 0.655294777491;
    twoBodyN2Map[{1, 4, 1, 4}] = 0.606886302485;
    twoBodyN2Map[{2, 4, 2, 4}] = 0.562386019246;
    twoBodyN2Map[{3, 4, 3, 4}] = 0.697573816772;
    twoBodyN2Map[{4, 4, 4, 4}] = 0.667701348660;
    twoBodyN2Map[{3, 5, 2, 0}] = 0.230002689806E-01;
    twoBodyN2Map[{5, 5, 0, 0}] = 0.188826272363E-01;
    twoBodyN2Map[{4, 5, 2, 1}] = 0.179848657808E-01;
    twoBodyN2Map[{5, 5, 1, 1}] = 0.179848745570E-01;
    twoBodyN2Map[{3, 5, 0, 2}] = 0.149722852886;
    twoBodyN2Map[{4, 5, 1, 2}] = 0.149376395297;
    twoBodyN2Map[{5, 5, 2, 2}] = 0.185346162128;
    twoBodyN2Map[{2, 5, 0, 3}] = 0.449318669954E-01;
    twoBodyN2Map[{5, 5, 3, 3}] = 0.519778158152E-01;
    twoBodyN2Map[{2, 5, 1, 4}] = 0.222501245895E-01;
    twoBodyN2Map[{5, 5, 4, 4}] = 0.269357801382E-01;
    twoBodyN2Map[{0, 5, 0, 5}] = 0.655295288869;
    twoBodyN2Map[{1, 5, 1, 5}] = 0.562386409245;
    twoBodyN2Map[{2, 5, 2, 5}] = 0.606886624363;
    twoBodyN2Map[{3, 5, 3, 5}] = 0.697574397066;
    twoBodyN2Map[{4, 5, 4, 5}] = 0.613830245027;
    twoBodyN2Map[{5, 5, 5, 5}] = 0.667702261948;
    for (const auto& el : twoBodyN2Map) {
      for (const auto& index :
           getSymmetricIndices<ERISymmetry::eightfold>({el.first[0], el.first[2], el.first[1], el.first[3]})) {
        twoBodyN2(index[0] * 6 + index[1], index[2] * 6 + index[3]) = el.second;
      }
    }
    nitrogenHamiltonianTwoBody.addTwoBodyContribution(twoBodyN2);
    // One-Body transcorrelated He Hamiltonian
    heTranscorrelatedHamiltonian.addTerm(0, 0, -1.9410228773342584e+00);
    heTranscorrelatedHamiltonian.addTerm(0, 1, -3.1641663736652353e-01);
    heTranscorrelatedHamiltonian.addTerm(0, 2, -3.1641663736652353e-10);
    heTranscorrelatedHamiltonian.addTerm(0, 3, -3.1641663736652353e-10);
    heTranscorrelatedHamiltonian.addTerm(0, 4, -3.1641663736652353e-10);
    heTranscorrelatedHamiltonian.addTerm(2, 0, -3.1641663736652353e-10);
    heTranscorrelatedHamiltonian.addTerm(3, 0, -3.1641663736652353e-10);
    heTranscorrelatedHamiltonian.addTerm(4, 0, -3.1641663736652353e-10);
    heTranscorrelatedHamiltonian.addTerm(1, 0, -3.1641663736652381e-01);
    heTranscorrelatedHamiltonian.addTerm(1, 1, -9.0227670561568662e-02);
    heTranscorrelatedHamiltonian.addTerm(2, 2, 7.8499729043522848e-01);
    heTranscorrelatedHamiltonian.addTerm(3, 3, 7.8499729043522803e-01);
    heTranscorrelatedHamiltonian.addTerm(4, 4, 7.8499729043522803e-01);
    Eigen::MatrixXd twoBodyHe = Eigen::MatrixXd::Zero(25, 25);
    std::unordered_map<std::array<int, 4>, double, boost::hash<std::array<int, 4>>> twoBodyHeMap;
    twoBodyHeMap[{0, 0, 0, 0}] = 1.0183556583001618e+00;
    twoBodyHeMap[{0, 1, 0, 0}] = 3.0476675403355524e-01;
    twoBodyHeMap[{0, 1, 0, 1}] = 2.5174364834784252e-01;
    twoBodyHeMap[{0, 2, 0, 2}] = 1.8330678004452303e-01;
    twoBodyHeMap[{0, 3, 0, 3}] = 1.8330678004452336e-01;
    twoBodyHeMap[{0, 4, 0, 4}] = 1.8330678004452336e-01;
    twoBodyHeMap[{0, 0, 1, 0}] = 3.2114288110691158e-01;
    twoBodyHeMap[{0, 0, 1, 1}] = 8.5366263969802025e-01;
    twoBodyHeMap[{0, 1, 1, 0}] = 2.2385496280331604e-01;
    twoBodyHeMap[{0, 1, 1, 1}] = 2.4714874264182601e-01;
    twoBodyHeMap[{0, 2, 1, 2}] = 2.5134760449013090e-02;
    twoBodyHeMap[{0, 3, 1, 3}] = 2.5134760449013145e-02;
    twoBodyHeMap[{0, 4, 1, 4}] = 2.5134760449013152e-02;
    twoBodyHeMap[{0, 0, 2, 2}] = 9.5116233984038423e-01;
    twoBodyHeMap[{0, 1, 2, 2}] = 2.1082535928181326e-01;
    twoBodyHeMap[{0, 2, 2, 0}] = 1.8098179143498963e-01;
    twoBodyHeMap[{0, 2, 2, 1}] = 4.4024363522944136e-02;
    twoBodyHeMap[{0, 0, 3, 3}] = 9.5116233984038412e-01;
    twoBodyHeMap[{0, 1, 3, 3}] = 2.1082535928181315e-01;
    twoBodyHeMap[{0, 3, 3, 0}] = 1.8098179143498991e-01;
    twoBodyHeMap[{0, 3, 3, 1}] = 4.4024363522944199e-02;
    twoBodyHeMap[{0, 0, 4, 4}] = 9.5116233984038412e-01;
    twoBodyHeMap[{0, 1, 4, 4}] = 2.1082535928181315e-01;
    twoBodyHeMap[{0, 4, 4, 0}] = 1.8098179143498991e-01;
    twoBodyHeMap[{0, 4, 4, 1}] = 4.4024363522944199e-02;
    twoBodyHeMap[{1, 0, 1, 0}] = 1.9596627725879007e-01;
    twoBodyHeMap[{1, 1, 1, 0}] = 2.5546841384749081e-01;
    twoBodyHeMap[{1, 1, 1, 1}] = 7.4940475756065861e-01;
    twoBodyHeMap[{1, 2, 1, 2}] = 3.6676840894247595e-02;
    twoBodyHeMap[{1, 3, 1, 3}] = 3.6676840894247623e-02;
    twoBodyHeMap[{1, 4, 1, 4}] = 3.6676840894247630e-02;
    twoBodyHeMap[{1, 0, 2, 2}] = 2.4775473575494669e-01;
    twoBodyHeMap[{1, 1, 2, 2}] = 7.8286400008585633e-01;
    twoBodyHeMap[{1, 2, 2, 0}] = 4.7410516610265629e-03;
    twoBodyHeMap[{1, 2, 2, 1}] = 2.5466242990294319e-02;
    twoBodyHeMap[{1, 0, 3, 3}] = 2.4775473575494658e-01;
    twoBodyHeMap[{1, 1, 3, 3}] = 7.8286400008585633e-01;
    twoBodyHeMap[{1, 3, 3, 0}] = 4.7410516610265786e-03;
    twoBodyHeMap[{1, 3, 3, 1}] = 2.5466242990294343e-02;
    twoBodyHeMap[{1, 0, 4, 4}] = 2.4775473575494658e-01;
    twoBodyHeMap[{1, 1, 4, 4}] = 7.8286400008585633e-01;
    twoBodyHeMap[{1, 4, 4, 0}] = 4.7410516610265751e-03;
    twoBodyHeMap[{1, 4, 4, 1}] = 2.5466242990294340e-02;
    twoBodyHeMap[{2, 0, 2, 0}] = 1.7865680282545623e-01;
    twoBodyHeMap[{2, 1, 2, 0}] = 2.3630654734957610e-02;
    twoBodyHeMap[{2, 1, 2, 1}] = 1.4255645086341050e-02;
    twoBodyHeMap[{2, 2, 2, 2}] = 9.3158111092842799e-01;
    twoBodyHeMap[{2, 2, 3, 3}] = 9.3158111092842788e-01;
    twoBodyHeMap[{2, 2, 4, 4}] = 9.3158111092842788e-01;
    twoBodyHeMap[{3, 0, 3, 0}] = 1.7865680282545646e-01;
    twoBodyHeMap[{3, 1, 3, 0}] = 2.3630654734957651e-02;
    twoBodyHeMap[{3, 1, 3, 1}] = 1.4255645086341067e-02;
    twoBodyHeMap[{3, 3, 3, 3}] = 9.3158111092842688e-01;
    twoBodyHeMap[{3, 3, 4, 4}] = 9.3158111092842688e-01;
    twoBodyHeMap[{4, 0, 4, 0}] = 1.7865680282545646e-01;
    twoBodyHeMap[{4, 1, 4, 0}] = 2.3630654734957651e-02;
    twoBodyHeMap[{4, 1, 4, 1}] = 1.4255645086341064e-02;
    twoBodyHeMap[{4, 4, 4, 4}] = 9.3158111092842688e-01;
    for (const auto& el : twoBodyHeMap) {
      for (const auto& index :
           getSymmetricIndices<ERISymmetry::eightfold>({el.first[0], el.first[2], el.first[1], el.first[3]})) {
        twoBodyHe(index[0] * 5 + index[1], index[2] * 5 + index[3]) = el.second;
      }
    }
    heTranscorrelatedHamiltonian.addTwoBodyContribution(twoBodyHeMap);
  }
};

TEST_F(SQUtilsTest, CanInitializeABoolVecONV) {
  auto determinant = ElectronicDeterminant(onv_2e5o, onv_2e6o);
  ASSERT_EQ(determinant.size(SpinComponent::Beta), 6);
}

TEST_F(SQUtilsTest, InitializeDeterminantFromONVString) {
  auto determinant1 = ElectronicDeterminant(onv_2e6o, onv_2e6o);
  auto determinant2 = ElectronicDeterminant("220000");
  ASSERT_TRUE(determinant1.isCompatibleWith(determinant2));
}

TEST_F(SQUtilsTest, InitializeAndExciteVeryLongBoolVecONV) {
  auto determinant = ElectronicDeterminant(onv_10e1000o, onv_10e1000o);
  for (int idx = 0; idx < 10000; idx++) {
    for (int iExc = 0; iExc < 100; iExc++) {
      determinant.applyExcitation({iExc, iExc + 100, SpinComponent::Alpha});
      determinant.applyExcitation({iExc, iExc + 100, SpinComponent::Beta});
      determinant.applyExcitation({iExc + 100, iExc, SpinComponent::Alpha});
      determinant.applyExcitation({iExc + 100, iExc, SpinComponent::Beta});
    }
  }
}

TEST_F(SQUtilsTest, CorrectExcitationDegreeEvaluation) {
  auto determinant1 = ElectronicDeterminant(onv_2e5o, onv_2e6o);
  auto determinant2 = ElectronicDeterminant(onv_2e5o, onv_2e6o);
  determinant1.applyExcitation({1, 3, SpinComponent::Alpha});
  ASSERT_EQ(determinant1.getExcitationDegree(determinant2).first, 1);
}

TEST_F(SQUtilsTest, CorrectNumberOfElectrons) {
  auto determinant = ElectronicDeterminant(onv_4e20o, onv_4e20o);
  auto n_ele = determinant.countElectrons().first + determinant.countElectrons().second;
  ASSERT_EQ(n_ele, 8);
}

TEST_F(SQUtilsTest, GenerateConnectedExcitationWithHBCI) {
  auto determinantHF = ElectronicDeterminant({1, 1, 0, 1, 0, 0}, {1, 0, 1, 0, 0, 1});
  nitrogenHamiltonianTwoBody.setHbciThreshold(0);
  nitrogenHamiltonianTwoBody.populateHBCIContainer();
  auto referenceConnected = nitrogenHamiltonianTwoBody.generateAllConnected(determinantHF);

  for (const auto& detRef : referenceConnected) {
    std::stringstream stream;
    // detRef.print(stream);
    SCOPED_TRACE(stream.str() + " : " +
                 std::to_string(nitrogenHamiltonianTwoBody.calculateMatrixElement(determinantHF, detRef)));
    ASSERT_TRUE(std::find(referenceConnected.begin(), referenceConnected.end(), detRef) != referenceConnected.end() ||
                (nitrogenHamiltonianTwoBody.calculateMatrixElement(determinantHF, detRef) <= 1e-13));
  }
}

TEST_F(SQUtilsTest, CountConnectedExcitation) {
  auto determinant = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  auto numberOfConnected = hydrogenHamiltonianTwoBody.getNumberOfConnectedDeterminants(determinant);
  ASSERT_EQ(numberOfConnected, 16);
}

TEST_F(SQUtilsTest, MatrixElementOneBodyHamiltonian) {
  auto determinant = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  double matrixElement = hydrogenHamiltonianOneBody.calculateMatrixElement(determinant, determinant);
  EXPECT_THAT(matrixElement, DoubleNear(-1.245095342160 * 2., 1e-13));
}

TEST_F(SQUtilsTest, MatrixElementOneBodyHamiltonianDoubleExcitation) {
  auto determinant1 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  auto determinant2 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  determinant2.applyExcitation({0, 1, SpinComponent::Alpha});
  determinant2.applyExcitation({0, 1, SpinComponent::Beta});
  double matrixElement = hydrogenHamiltonianOneBody.calculateMatrixElement(determinant1, determinant2);
  EXPECT_THAT(matrixElement, DoubleNear(0., 1e-13));
}

TEST_F(SQUtilsTest, MatrixElementTwoBodyHamiltonian) {
  auto determinant = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  double matrixElement = hydrogenHamiltonianTwoBody.calculateMatrixElement(determinant, determinant);
  EXPECT_THAT(matrixElement, DoubleNear(0.649702723853, 1e-13));
}

TEST_F(SQUtilsTest, CheckDoubleExcitation) {
  auto determinant1 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  auto determinant2 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  determinant2.applyExcitation({0, 1, SpinComponent::Alpha});
  determinant2.applyExcitation({0, 1, SpinComponent::Beta});
  double matrixElement = hydrogenHamiltonianTwoBody.calculateMatrixElement(determinant1, determinant2);
  EXPECT_THAT(matrixElement, DoubleNear(0.109300885764, 1e-13));
}

TEST_F(SQUtilsTest, CheckDoubleExcitationWithChangeOfSign) {
  // This is an excitation with an odd parity, therefore the sign is expected to change wrt the two-electron integral.
  auto determinant1 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  auto determinant2 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  determinant1.applyExcitation({0, 3, SpinComponent::Alpha});
  determinant1.applyExcitation({0, 2, SpinComponent::Beta});
  determinant2.applyExcitation({0, 2, SpinComponent::Alpha});
  determinant2.applyExcitation({0, 3, SpinComponent::Beta});
  double matrixElement = hydrogenHamiltonianTwoBody.calculateMatrixElement(determinant1, determinant2);
  double matrixElement2 = hydrogenHamiltonianTwoBody.calculateMatrixElement(determinant2, determinant1);
  EXPECT_THAT(matrixElement, DoubleNear(-0.675771906037E-01, 1e-13));
  EXPECT_THAT(matrixElement - matrixElement2, DoubleNear(0., 1e-13));
}

TEST_F(SQUtilsTest, CheckSingleExcitationSign) {
  auto determinant1 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  auto determinant2 = ElectronicDeterminant(onv_1e4o, onv_1e4o);
  determinant2.applyExcitation({0, 1, SpinComponent::Alpha});
  double matrixElement = hydrogenHamiltonianTwoBody.calculateMatrixElement(determinant1, determinant2);
  double matrixElement2 = hydrogenHamiltonianTwoBody.calculateMatrixElement(determinant2, determinant1);
  matrixElement += hydrogenHamiltonianOneBody.calculateMatrixElement(determinant1, determinant2);
  matrixElement2 += hydrogenHamiltonianOneBody.calculateMatrixElement(determinant2, determinant1);
  EXPECT_THAT(matrixElement - matrixElement2, DoubleNear(0., 1e-13));
}

TEST_F(SQUtilsTest, CheckSingleExcitationSignAlphaAlpha) {
  auto determinant1 = ElectronicDeterminant({1, 1, 1, 0, 0, 0}, {1, 1, 1, 0, 0, 0});
  auto determinant2 = ElectronicDeterminant({0, 0, 1, 1, 1, 0}, {1, 1, 1, 0, 0, 0});
  double matrixElement = nitrogenHamiltonianTwoBody.calculateMatrixElement(determinant1, determinant2);
  EXPECT_THAT(matrixElement, DoubleNear(0.230002415882E-01 - 0.149722796368, 1.0E-13));
}

TEST_F(SQUtilsTest, CheckSingleExcitationSignBetaBeta) {
  auto determinant1 = ElectronicDeterminant({1, 1, 1, 0, 0, 0}, {1, 1, 1, 0, 0, 0});
  auto determinant2 = ElectronicDeterminant({1, 1, 1, 0, 0, 0}, {0, 0, 1, 1, 1, 0});
  double matrixElement = nitrogenHamiltonianTwoBody.calculateMatrixElement(determinant1, determinant2);
  EXPECT_THAT(matrixElement, DoubleNear(0.230002415882E-01 - 0.149722796368, 1.0E-13));
}

TEST_F(SQUtilsTest, CheckHFEnergyN2) {
  auto determinantHF = ElectronicDeterminant({1, 1, 1, 0, 0, 0}, {1, 1, 1, 0, 0, 0});
  double HFenergy = nitrogenHamiltonianOneBody.calculateMatrixElement(determinantHF, determinantHF) +
                    nitrogenHamiltonianTwoBody.calculateMatrixElement(determinantHF, determinantHF);
  HFenergy -= 95.4968869041;
  EXPECT_THAT(HFenergy, DoubleNear(-108.9515677, 1.0E-4));
}

TEST_F(SQUtilsTest, CheckSlaterCondon) {
  auto determinantHF = ElectronicDeterminant({1, 1, 1, 0, 0, 0}, {1, 1, 1, 0, 0, 0});
  auto determinantCIS = ElectronicDeterminant({1, 1, 1, 0, 0, 0}, {1, 1, 0, 1, 0, 0});
  double HFenergy = nitrogenHamiltonianOneBody.calculateMatrixElement(determinantHF, determinantCIS) +
                    nitrogenHamiltonianTwoBody.calculateMatrixElement(determinantHF, determinantCIS);
  EXPECT_THAT(HFenergy, DoubleNear(0., 1.0E-6));
}

TEST_F(SQUtilsTest, FullCIEnergyH2) {
  double energy = 0.;
  for (auto&& it1 : H2CICoeff) {
    for (auto&& it2 : H2CICoeff) {
      energy += hydrogenHamiltonianOneBody.calculateMatrixElement(it1.first, it2.first) * it1.second * it2.second;
      energy += hydrogenHamiltonianTwoBody.calculateMatrixElement(it1.first, it2.first) * it1.second * it2.second;
    }
  }
  EXPECT_THAT(energy, DoubleNear(-1.865436725511013, 1e-13));
}

TEST_F(SQUtilsTest, FullCIEnergyN2) {
  double energy = 0.;
  for (auto&& it1 : N2CICoeff) {
    for (auto&& it2 : N2CICoeff) {
      energy += nitrogenHamiltonianOneBody.calculateMatrixElement(it1.first, it2.first) * it1.second * it2.second;
      energy += nitrogenHamiltonianTwoBody.calculateMatrixElement(it1.first, it2.first) * it1.second * it2.second;
    }
  }
  // Add the core part
  energy -= 95.4968869041;
  EXPECT_THAT(energy, DoubleNear(-109.09000800, 1e-8));
}

TEST_F(SQUtilsTest, CIDeterminantsCreatorCreatesCorrectDeterminantsMrciAll) {
  CIDeterminantsCreator creator({"222000"});
  auto mrciDets23 = creator.generateAllReferences({2, 3});
  auto mrciDets13 = creator.generateAllReferences({1, 3});
  std::vector<ElectronicDeterminant> expectedDets23{{"222000"}, {"220200"}, {"22ab00"}, {"22ba00"}};
  std::vector<ElectronicDeterminant> expectedDets13{{"222000"}, {"202200"}, {"2a2b00"}, {"2b2a00"}};

  EXPECT_EQ(mrciDets23.size(), expectedDets23.size());
  EXPECT_EQ(mrciDets13.size(), expectedDets13.size());
  for (const auto& det : mrciDets23) {
    EXPECT_TRUE(std::find(expectedDets23.begin(), expectedDets23.end(), det) != expectedDets23.end());
  }
  for (const auto& det : mrciDets13) {
    EXPECT_TRUE(std::find(expectedDets13.begin(), expectedDets13.end(), det) != expectedDets13.end());
  }

  // test that code is routed correctly
  auto mrciDets = creator.generate({1, 2, 3}, false);

  for (const auto& det : mrciDets13) {
    EXPECT_TRUE(std::find(mrciDets.begin(), mrciDets.end(), det) != mrciDets.end());
  }
}

TEST_F(SQUtilsTest, CIDeterminantsCreatorCreatesCorrectDeterminantsMrciOnlyDoubles) {
  CIDeterminantsCreator creator({"222000"});
  // Test that double excitations are generated for closed shell orbitals
  auto mrciDoubleDets23 = creator.generateOnlyDoubleExcitations({2, 3});
  auto mrciDoubleDets13 = creator.generateOnlyDoubleExcitations({1, 3});
  std::vector<ElectronicDeterminant> expectedDets23{{"222000"}, {"220200"}};
  std::vector<ElectronicDeterminant> expectedDets13{{"222000"}, {"202200"}};

  EXPECT_EQ(mrciDoubleDets23.size(), expectedDets23.size());
  EXPECT_EQ(mrciDoubleDets13.size(), expectedDets13.size());
  for (const auto& det : mrciDoubleDets23) {
    EXPECT_TRUE(std::find(expectedDets23.begin(), expectedDets23.end(), det) != expectedDets23.end());
  }
  for (const auto& det : mrciDoubleDets13) {
    EXPECT_TRUE(std::find(expectedDets13.begin(), expectedDets13.end(), det) != expectedDets13.end());
  }

  // Test that open shell orbitals are not touched
  creator = CIDeterminantsCreator{{"22a000"}};
  auto mrciDoubleDetsOpenShell13 = creator.generateOnlyDoubleExcitations({1, 2, 3});
  std::vector<ElectronicDeterminant> expectedDetsOpenShell13{{"22a000"}, {"20a200"}};
  EXPECT_EQ(mrciDoubleDetsOpenShell13.size(), expectedDetsOpenShell13.size());
  for (const auto& det : mrciDoubleDetsOpenShell13) {
    EXPECT_TRUE(std::find(expectedDetsOpenShell13.begin(), expectedDetsOpenShell13.end(), det) !=
                expectedDetsOpenShell13.end());
  }

  // test that code is routed correctly
  auto mrciDets = creator.generate({1, 2, 3}, true);

  for (const auto& det : mrciDoubleDetsOpenShell13) {
    EXPECT_TRUE(std::find(mrciDets.begin(), mrciDets.end(), det) != mrciDets.end());
  }
}

TEST_F(SQUtilsTest, CIDeterminantsCreatorChecksInputIndices) {
  CIDeterminantsCreator creator({"222000"});
  // No duplicate in indices
  ASSERT_THROW(creator.generateAllReferences({1, 2, 1, 3}), std::runtime_error);
  // Checks elements
  ASSERT_THROW(creator.generateAllReferences({1, 2, -3}), std::runtime_error);
  ASSERT_THROW(creator.generateAllReferences({1, 2, 7}), std::runtime_error);
}

TEST_F(SQUtilsTest, ConstructTranscorrelatedCIMatrix) {
  auto determinantHFHe = ElectronicDeterminant({1, 0, 0, 0, 0}, {1, 0, 0, 0, 0});
  auto referenceConnected = heTranscorrelatedHamiltonian.generateAllConnected(determinantHFHe);
  Eigen::MatrixXd hamiltonianMatrix = Eigen::MatrixXd::Zero(referenceConnected.size(), referenceConnected.size());
  for (int iRow = 0; iRow < int(referenceConnected.size()); iRow++) {
    for (int iCol = 0; iCol < int(referenceConnected.size()); iCol++) {
      hamiltonianMatrix(iRow, iCol) =
          heTranscorrelatedHamiltonian.calculateMatrixElement(referenceConnected[iRow], referenceConnected[iCol]);
    }
  }

  Eigen::EigenSolver<Eigen::MatrixXd> es(hamiltonianMatrix);
  Eigen::VectorXcd eigenValues = es.eigenvalues();
  Eigen::VectorXd expectedEigenValues(25);
  expectedEigenValues << -2.89572750407, -0.960429568619, 0.587004952037, -1.401442871, 2.52225418712, -0.0264543683241,
      -0.0264543683241, -0.0264543683241, -0.391800797439, -0.391800797439, -0.391800797439, 1.50567277565,
      1.50567277565, 1.50567277565, 1.45812313591, 1.45812313591, 1.45812313591, 2.5015756918, 2.5015756918,
      2.5015756918, 2.5015756918, 2.5015756918, 2.5015756918, 2.5015756918, 2.5015756918;
  std::sort(expectedEigenValues.data(), expectedEigenValues.data() + 25);
  std::sort(eigenValues.data(), eigenValues.data() + 25, [&](auto a, auto b) { return std::real(a) < std::real(b); });
  for (int i = 0; i < 25; i++) {
    //  std::cout << eigenValues[i].real() << "  " << expectedEigenValues[i] << std::endl;
    EXPECT_NEAR(eigenValues[i].real(), expectedEigenValues[i], 1e-10);
  }
}

TEST_F(SQUtilsTest, ConstructTranscorrelatedCIMatrixWithThreeBody) {
  // auto determinantHFHe = ElectronicDeterminant({1, 0, 0, 0, 0}, {1, 0, 0, 0, 0});
  // auto referenceConnected = heTranscorrelatedHamiltonian.generateDoubleExcitations(determinantHFHe);
  // Eigen::MatrixXd hamiltonianMatrix = Eigen::MatrixXd::Zero(referenceConnected.size(), referenceConnected.size());
  // for (int iRow = 0; iRow < referenceConnected.size(); iRow++)
  //   for (int iCol = 0; iCol < referenceConnected.size(); iCol++)
  //     hamiltonianMatrix(iRow, iCol) = heTranscorrelatedHamiltonian.calculateMatrixElement(referenceConnected[iRow],
  //                                                                                         referenceConnected[iCol]);
  // Eigen::EigenSolver<Eigen::MatrixXd> es(hamiltonianMatrix);
  // Eigen::VectorXcd eigenValues = es.eigenvalues();
  // Eigen::VectorXd expectedEigenValues(25);
  // expectedEigenValues << -2.89572750407,
  //                        -0.960429568619,
  //                         0.587004952037,
  //                        -1.401442871,
  //                         2.52225418712,
  //                        -0.0264543683241,
  //                        -0.0264543683241,
  //                        -0.0264543683241,
  //                        -0.391800797439,
  //                        -0.391800797439,
  //                        -0.391800797439,
  //                         1.50567277565,
  //                         1.50567277565,
  //                         1.50567277565,
  //                         1.45812313591,
  //                         1.45812313591,
  //                         1.45812313591,
  //                         2.5015756918,
  //                         2.5015756918,
  //                         2.5015756918,
  //                         2.5015756918,
  //                         2.5015756918,
  //                         2.5015756918,
  //                         2.5015756918,
  //                         2.5015756918;
  // for(int i=0; i<25; i++){
  //   // std::cout << eigenValues[i].real() << std::endl;
  //   ASSERT_NEAR(eigenValues[i].real(), expectedEigenValues[i], 1e-10);
  // }
}

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine
