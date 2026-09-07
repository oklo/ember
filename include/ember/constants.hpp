#pragma once
// Physical constants in CGS.
//
// Values are CODATA 2018 and IAU 2015 Resolution B3 nominal solar quantities.
// The IAU values matter: "the solar radius" is a *defined* nominal constant
// (6.957e10 cm), not a measurement, so a model that reproduces it is
// reproducing an agreed yardstick rather than chasing a moving number.
namespace ember::constants {

inline constexpr double c        = 2.99792458e10;    // cm/s (exact)
inline constexpr double G        = 6.67430e-8;       // cm^3 g^-1 s^-2
inline constexpr double h        = 6.62607015e-27;   // erg s (exact)
inline constexpr double kB       = 1.380649e-16;     // erg/K (exact)
inline constexpr double NA       = 6.02214076e23;    // 1/mol (exact)
inline constexpr double amu      = 1.66053906660e-24;// g
inline constexpr double me       = 9.1093837015e-28; // g
inline constexpr double eV       = 1.602176634e-12;  // erg (exact)
inline constexpr double sigma_SB = 5.670374419e-5;   // erg cm^-2 s^-1 K^-4
inline constexpr double a_rad    = 4.0 * sigma_SB / c;
inline constexpr double R_gas    = NA * kB;          // erg /(g K) per mean molecular weight

// IAU 2015 nominal solar values (exact by definition)
inline constexpr double Lsun = 3.828e33;   // erg/s
inline constexpr double Rsun = 6.957e10;   // cm
inline constexpr double GMsun = 1.3271244e26; // cm^3/s^2
inline constexpr double Msun = GMsun / G;  // g  (~1.98892e33)
inline constexpr double Teff_sun = 5772.0; // K (nominal, from Lsun and Rsun)

// Solar age, Bahcall/Serenelli meteoritic value
inline constexpr double t_sun = 4.57e9 * 3.1557e7; // s

inline constexpr double yr = 3.1557e7;     // Julian year in s
} // namespace ember::constants
