/**
 * @file constants.hpp
 * @brief Core constants for NRLMSIS 2.1 parameter handling.
 * @author Watosn
 */
#pragma once

#include <array>
#include <cstddef>

namespace msis21::detail {

inline constexpr double kPi = 3.1415926535897932384626433832795;
inline constexpr double kDeg2Rad = kPi / 180.0;
inline constexpr double kDoy2Rad = 2.0 * kPi / 365.0;
inline constexpr double kLst2Rad = kPi / 12.0;

inline constexpr std::size_t kMaxBasisFunctions = 512;
inline constexpr std::size_t kParmColumnCount = 131;
inline constexpr int kNd = 27;
inline constexpr int kSplineOrder = 4;
inline constexpr int kNl = kNd - kSplineOrder;
inline constexpr int kNdO1 = 13;
inline constexpr int kNsplO1 = kNdO1 - 5;
inline constexpr int kNdNo = 13;
inline constexpr int kNsplNo = kNdNo - 5;
inline constexpr int kIzfmx = 13;
inline constexpr int kIzfx = 14;
inline constexpr int kIzax = 17;
inline constexpr int kItex = kNl;
inline constexpr int kItgb0 = kNl - 1;
inline constexpr int kItb0 = kNl - 2;
inline constexpr double kZetaF = 70.0;
inline constexpr double kZetaA = 85.0;
inline constexpr double kZetaB = 122.5;
inline constexpr double kKb = 1.380649e-23;
inline constexpr double kNa = 6.02214076e23;
inline constexpr double kG0 = 9.80665;
inline constexpr double kG0DivKb = kG0 / kKb * 1.0e3;
inline constexpr double kMbar = 28.96546 / (1.0e3 * kNa);
inline constexpr double kMbarG0DivKb = kMbar * kG0 / kKb * 1.0e3;
// Legacy Fortran defines lnP0 without kind suffix, so it is rounded to default
// REAL(4) before promotion in DBLE mode. Keep that exact value for parity.
inline constexpr double kLnP0 = 11.515613555908203;
inline constexpr double kToa = 4000.0;
inline constexpr double kHoa =
    (kKb * kToa) / ((16.0 / (1.0e3 * kNa)) * kG0) * 1.0e-3;
inline constexpr double kDmissing = 9.999e-38;
inline constexpr std::array<double, 10> kSpecMass = {
    0.0,
    28.0134 / (1.0e3 * kNa),
    31.9988 / (1.0e3 * kNa),
    (31.9988 / 2.0) / (1.0e3 * kNa),
    4.0 / (1.0e3 * kNa),
    1.0 / (1.0e3 * kNa),
    39.948 / (1.0e3 * kNa),
    (28.0134 / 2.0) / (1.0e3 * kNa),
    (31.9988 / 2.0) / (1.0e3 * kNa),
    ((28.0134 + 31.9988) / 2.0) / (1.0e3 * kNa)};
inline constexpr std::array<double, 10> kLnVmr = {
    0.0,
    -0.24737477036295383,   // log(0.780848)
    -1.563556737177912,     // log(0.209390)
    0.0,
    -12.166851932376892,    // log(5.2e-6)
    0.0,
    -4.674305924822954,     // log(0.009332)
    0.0,
    0.0,
    0.0};
inline constexpr int kMbf = 383;
inline constexpr int kMaxLegendreDegree = 6;
inline constexpr int kAmaxN = 6;
inline constexpr int kAmaxS = 2;
inline constexpr int kTmaxL = 3;
inline constexpr int kTmaxN = 6;
inline constexpr int kTmaxS = 2;
inline constexpr int kPmaxM = 2;
inline constexpr int kPmaxN = 6;
inline constexpr int kPmaxS = 2;
inline constexpr int kNsfx = 5;
inline constexpr int kNsfxMod = 5;
inline constexpr int kNmag = 54;
inline constexpr int kNut = 12;
inline constexpr int kCtimeInd = 0;
inline constexpr int kCintAnn = kCtimeInd + (kAmaxN + 1);
inline constexpr int kCtide = kCintAnn + ((kAmaxN + 1) * 2 * kAmaxS);
inline constexpr int kCspw =
    kCtide + (4 * kTmaxS + 2) * (kTmaxL * (kTmaxN + 1) - (kTmaxL * (kTmaxL + 1)) / 2);
inline constexpr int kCsfx =
    kCspw + (4 * kPmaxS + 2) * (kPmaxM * (kPmaxN + 1) - (kPmaxM * (kPmaxM + 1)) / 2);
inline constexpr int kCextra = kCsfx + kNsfx;
inline constexpr int kCnonLin = kMbf + 1;
inline constexpr int kCsfxMod = kCnonLin;
inline constexpr int kCmag = kCsfxMod + kNsfxMod;
inline constexpr int kCut = kCmag + kNmag;
inline constexpr double kSfluxAvgReference = 150.0;
inline constexpr double kSfluxAvgQuadCutoff = 150.0;
inline constexpr double kZetagamma = 100.0;
inline constexpr double kHgamma = 1.0 / 30.0;
inline constexpr std::array<double, 30> kNodesTN = {
    -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0,
    35.0,  40.0,  45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0,
    85.0,  92.5,  102.5, 112.5, 122.5, 132.5, 142.5, 152.5, 162.5, 172.5};
inline constexpr std::array<double, 14> kNodesO1 = {
    35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 92.5, 102.5, 112.5};
inline constexpr std::array<double, 14> kNodesNO = {
    47.5, 55.0, 62.5, 70.0, 77.5, 85.0, 92.5, 100.0, 107.5, 115.0, 122.5, 130.0, 137.5, 145.0};
inline constexpr std::array<double, 3> kS4zetaA = {
    0.257142857142857, 0.653968253968254, 0.088888888888889};
inline constexpr std::array<double, 3> kWghtAxdz = {
    -0.102857142857, 0.0495238095238, 0.053333333333};
inline constexpr std::array<double, 3> kS4zetaF = {
    0.166666666666667, 0.666666666666667, 0.166666666666667};
inline constexpr std::array<double, 3> kS5zeta0 = {
    0.458333333333333, 0.458333333333333, 0.041666666666667};
inline constexpr std::array<double, 4> kS5zetaB = {
    0.041666666666667, 0.458333333333333, 0.458333333333333, 0.041666666666667};
inline constexpr std::array<double, 4> kS5zetaA = {
    0.085714285714286, 0.587590187590188, 0.313020313020313, 0.013675213675214};
inline constexpr std::array<double, 4> kS5zetaF = {
    0.041666666666667, 0.458333333333333, 0.458333333333333, 0.041666666666667};
inline constexpr std::array<double, 5> kS6zetaB = {
    0.008771929824561, 0.216228070175439, 0.55, 0.216666666666667, 0.008333333333333};
inline constexpr std::array<double, 5> kS6zetaA = {
    0.023376623376623, 0.378732378732379, 0.500743700743701, 0.095538448479625,
    0.001608848667672};
inline constexpr std::array<double, 2> kC1o1Adj = {0.257142857142857, -0.102857142686844};
inline constexpr std::array<double, 4> kC1o1 = {
    1.75, -2.916666573405061, -1.624999900076852, 21.458332647194382};
inline constexpr std::array<double, 2> kC1noAdj = {0.166666666666667, -0.066666666666667};
inline constexpr std::array<double, 4> kC1no = {1.5, -3.75, 0.0, 15.0};

}  // namespace msis21::detail
