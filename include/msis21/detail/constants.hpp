/**
 * @file constants.hpp
 * @brief Core constants for NRLMSIS 2.1 parameter handling.
 * @author Watosn
 */
#pragma once

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
inline constexpr int kIzfmx = 13;
inline constexpr int kIzfx = 14;
inline constexpr int kIzax = 17;
inline constexpr int kItex = kNl;
inline constexpr int kItgb0 = kNl - 1;
inline constexpr int kItb0 = kNl - 2;
inline constexpr double kZetaF = 70.0;
inline constexpr double kZetaA = 85.0;
inline constexpr double kZetaB = 122.5;
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

}  // namespace msis21::detail
