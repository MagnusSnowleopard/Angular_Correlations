#ifndef RACAH_CALC_H
#define RACAH_CALC_H

#include <cmath>
#include <cstdint>
#include <array>
#include <limits>

namespace RacahCalc {

// ============================================================
// Data container for the six coefficients you want:
//   rk01, rk11, rk21, rk02, rk12, rk22
// mapped as:
//   rk01 = F2(1,1)
//   rk11 = F2(1,2)
//   rk21 = F2(2,2)
//   rk02 = F4(1,1)
//   rk12 = F4(1,2)
//   rk22 = F4(2,2)
// ============================================================
struct Coeffs {
    double rk01 = 0.0;
    double rk11 = 0.0;
    double rk21 = 0.0;
    double rk02 = 0.0;
    double rk12 = 0.0;
    double rk22 = 0.0;
};

// ============================================================
// Internal helpers
// ============================================================
namespace detail {

inline bool nearly_integer(double x, double tol = 1e-10) {
    return std::fabs(x - std::round(x)) < tol;
}

inline bool is_half_integer(double j, double tol = 1e-10) {
    return nearly_integer(2.0 * j, tol);
}

inline int to_twoj(double j) {
    return static_cast<int>(std::llround(2.0 * j));
}

inline int phase_m1_from_int(int n) {
    return (n & 1) ? -1 : +1;
}

// log(n!) for integer n >= 0
inline double log_fact_int(int n) {
    if (n < 0) return std::numeric_limits<double>::quiet_NaN();
    return std::lgamma(static_cast<double>(n) + 1.0);
}

// All arguments are in "2*j" integer representation.
// Delta(a,b,c) = sqrt( (x! y! z!) / (w!) )
// with x=(a+b-c)/2, y=(a-b+c)/2, z=(-a+b+c)/2, w=(a+b+c)/2 + 1
inline double delta_coeff_2j(int ta, int tb, int tc) {
    // Triangle + parity
    if (ta < 0 || tb < 0 || tc < 0) return 0.0;
    if ((ta + tb + tc) & 1) return 0.0;  // sum must be integer in j-space
    if (ta + tb < tc || ta + tc < tb || tb + tc < ta) return 0.0;

    const int x = (ta + tb - tc) / 2;
    const int y = (ta - tb + tc) / 2;
    const int z = (-ta + tb + tc) / 2;
    const int w = (ta + tb + tc) / 2 + 1;

    if (x < 0 || y < 0 || z < 0 || w < 0) return 0.0;

    const double logv =
        log_fact_int(x) + log_fact_int(y) + log_fact_int(z) - log_fact_int(w);

    return std::exp(0.5 * logv);
}

// ============================================================
// Wigner 3-j symbol
// Computes:
//   ( j1  j2  j3 )
//   ( m1  m2  m3 )
// using 2*j, 2*m integer notation.
//
// Returns 0 if selection rules fail.
// ============================================================
inline double wigner3j_2j(int tj1, int tj2, int tj3,
                          int tm1, int tm2, int tm3)
{
    // Selection rules
    if (tm1 + tm2 + tm3 != 0) return 0.0;
    if (tj1 < 0 || tj2 < 0 || tj3 < 0) return 0.0;

    // m bounds and parity consistency (j and m same parity)
    if (std::abs(tm1) > tj1 || std::abs(tm2) > tj2 || std::abs(tm3) > tj3) return 0.0;
    if (((tj1 - tm1) & 1) || ((tj2 - tm2) & 1) || ((tj3 - tm3) & 1)) return 0.0;

    // triangle
    if (tj1 + tj2 < tj3 || tj1 + tj3 < tj2 || tj2 + tj3 < tj1) return 0.0;
    if ((tj1 + tj2 + tj3) & 1) return 0.0;

    // Phase factor (-1)^(j1-j2-m3)
    const int phase_exp_num = (tj1 - tj2 - tm3); // equals 2*(integer exponent)
    if (phase_exp_num & 1) return 0.0; // exponent must be integer
    const int phase = phase_m1_from_int(phase_exp_num / 2);

    const double delta = delta_coeff_2j(tj1, tj2, tj3);
    if (delta == 0.0) return 0.0;

    // sqrt factorial part:
    // sqrt[(j1+m1)!(j1-m1)!...(j3-m3)!] in integer factorial arguments
    const int a1 = (tj1 + tm1) / 2;
    const int a2 = (tj1 - tm1) / 2;
    const int b1 = (tj2 + tm2) / 2;
    const int b2 = (tj2 - tm2) / 2;
    const int c1 = (tj3 + tm3) / 2;
    const int c2 = (tj3 - tm3) / 2;

    if (a1 < 0 || a2 < 0 || b1 < 0 || b2 < 0 || c1 < 0 || c2 < 0) return 0.0;

    const double log_sqrt_fac =
        0.5 * (log_fact_int(a1) + log_fact_int(a2) +
               log_fact_int(b1) + log_fact_int(b2) +
               log_fact_int(c1) + log_fact_int(c2));

    // Racah sum:
    // denominator arguments:
    // t,
    // (j1+j2-j3)-t,
    // (j1-m1)-t,
    // (j2+m2)-t,
    // (j3-j2+m1)+t,
    // (j3-j1-m2)+t
    const int A = (tj1 + tj2 - tj3) / 2;
    const int B = (tj1 - tm1) / 2;
    const int C = (tj2 + tm2) / 2;
    const int D = (tj3 - tj2 + tm1) / 2;
    const int E = (tj3 - tj1 - tm2) / 2;

    // t bounds from nonnegative denominator arguments
    const int tmin = std::max(0, std::max(-D, -E));
    const int tmax = std::min(A, std::min(B, C));

    if (tmax < tmin) return 0.0;

    double sum = 0.0;
    for (int t = tmin; t <= tmax; ++t) {
        const int d1 = t;
        const int d2 = A - t;
        const int d3 = B - t;
        const int d4 = C - t;
        const int d5 = D + t;
        const int d6 = E + t;

        if (d1 < 0 || d2 < 0 || d3 < 0 || d4 < 0 || d5 < 0 || d6 < 0) continue;

        const double log_term =
            - (log_fact_int(d1) + log_fact_int(d2) + log_fact_int(d3) +
               log_fact_int(d4) + log_fact_int(d5) + log_fact_int(d6));

        const double term = phase_m1_from_int(t) * std::exp(log_term);
        sum += term;
    }

    const double val = phase * delta * std::exp(log_sqrt_fac) * sum;
    return (std::fabs(val) < 1e-15) ? 0.0 : val;
}

// ============================================================
// Clebsch-Gordan coefficient
// <j1 m1 j2 m2 | J M>
// via relation to 3-j:
//   <j1m1 j2m2 | JM> = (-1)^(j1-j2+M) sqrt(2J+1) (j1 j2 J; m1 m2 -M)
// ============================================================
inline double clebsch_gordan_2j(int tj1, int tm1,
                                int tj2, int tm2,
                                int tJ,  int tM)
{
    if (tm1 + tm2 != tM) return 0.0;

    const int phase_exp_num = (tj1 - tj2 + tM); // 2 * exponent
    if (phase_exp_num & 1) return 0.0;
    const int phase = phase_m1_from_int(phase_exp_num / 2);

    const double threej = wigner3j_2j(tj1, tj2, tJ, tm1, tm2, -tM);
    if (threej == 0.0) return 0.0;

    const double pref = std::sqrt(static_cast<double>(tJ + 1)); // 2J+1 = tJ+1
    const double val = phase * pref * threej;
    return (std::fabs(val) < 1e-15) ? 0.0 : val;
}

// ============================================================
// Wigner 6-j symbol using Racah formula
// { a b c
//   d e f }
// all args supplied as 2*j integers.
// ============================================================
inline double wigner6j_2j(int ta, int tb, int tc,
                          int td, int te, int tf)
{
    // Quick triangle checks for the 4 required triads:
    // (a,b,c), (a,e,f), (d,b,f), (d,e,c)
    const double d1 = delta_coeff_2j(ta, tb, tc);
    const double d2 = delta_coeff_2j(ta, te, tf);
    const double d3 = delta_coeff_2j(td, tb, tf);
    const double d4 = delta_coeff_2j(td, te, tc);
    if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0 || d4 == 0.0) return 0.0;

    // z bounds (all in integer z = sums in j-units)
    const int A1 = (ta + tb + tc) / 2;
    const int A2 = (ta + te + tf) / 2;
    const int A3 = (td + tb + tf) / 2;
    const int A4 = (td + te + tc) / 2;

    const int B1 = (ta + tb + td + te) / 2;
    const int B2 = (ta + tc + td + tf) / 2;
    const int B3 = (tb + tc + te + tf) / 2;

    const int zmin = std::max(std::max(A1, A2), std::max(A3, A4));
    const int zmax = std::min(B1, std::min(B2, B3));
    if (zmax < zmin) return 0.0;

    double sum = 0.0;
    for (int z = zmin; z <= zmax; ++z) {
        const int n1 = z - A1;
        const int n2 = z - A2;
        const int n3 = z - A3;
        const int n4 = z - A4;
        const int n5 = B1 - z;
        const int n6 = B2 - z;
        const int n7 = B3 - z;

        if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 || n5 < 0 || n6 < 0 || n7 < 0) continue;

        const double log_term =
            log_fact_int(z + 1)
            - (log_fact_int(n1) + log_fact_int(n2) + log_fact_int(n3) + log_fact_int(n4)
               + log_fact_int(n5) + log_fact_int(n6) + log_fact_int(n7));

        const double term = phase_m1_from_int(z) * std::exp(log_term);
        sum += term;
    }

    const double val = d1 * d2 * d3 * d4 * sum;
    return (std::fabs(val) < 1e-15) ? 0.0 : val;
}

} // namespace detail

// ============================================================
// Public low-level wrappers using physical j,m values (double)
// Supports integer and half-integer values exactly (within tolerance).
// ============================================================

inline bool IsAllowedJ(double j) {
    return (j >= 0.0) && detail::is_half_integer(j);
}

inline double Wigner3j(double j1, double j2, double j3,
                       double m1, double m2, double m3)
{
    if (!detail::is_half_integer(j1) || !detail::is_half_integer(j2) || !detail::is_half_integer(j3)) return 0.0;
    if (!detail::is_half_integer(m1) || !detail::is_half_integer(m2) || !detail::is_half_integer(m3)) return 0.0;

    return detail::wigner3j_2j(detail::to_twoj(j1), detail::to_twoj(j2), detail::to_twoj(j3),
                               detail::to_twoj(m1), detail::to_twoj(m2), detail::to_twoj(m3));
}

inline double ClebschGordan(double j1, double m1,
                            double j2, double m2,
                            double J,  double M)
{
    if (!detail::is_half_integer(j1) || !detail::is_half_integer(j2) || !detail::is_half_integer(J)) return 0.0;
    if (!detail::is_half_integer(m1) || !detail::is_half_integer(m2) || !detail::is_half_integer(M)) return 0.0;

    return detail::clebsch_gordan_2j(detail::to_twoj(j1), detail::to_twoj(m1),
                                     detail::to_twoj(j2), detail::to_twoj(m2),
                                     detail::to_twoj(J),  detail::to_twoj(M));
}

inline double Wigner6j(double a, double b, double c,
                       double d, double e, double f)
{
    if (!IsAllowedJ(a) || !IsAllowedJ(b) || !IsAllowedJ(c) ||
        !IsAllowedJ(d) || !IsAllowedJ(e) || !IsAllowedJ(f)) return 0.0;

    return detail::wigner6j_2j(detail::to_twoj(a), detail::to_twoj(b), detail::to_twoj(c),
                               detail::to_twoj(d), detail::to_twoj(e), detail::to_twoj(f));
}

// ============================================================
// Candidate Rose/Brink-style coefficient
//
// F_k(L,L';J1,J2) = (-1)^(J1-J2-1)
//                   * sqrt((2L+1)(2L'+1)(2J1+1))
//                   * <L,1 ; L',-1 | k,0>
//                   * { J1 J1 k ; L L' J2 }
//
// This is the convention used in the Python checker I gave you.
// If your table differs by a systematic phase convention, this is
// the place to tweak.
// ============================================================
inline double Fk_candidate(int k, int L, int Lp, double J1, double J2)
{
    if (!IsAllowedJ(J1) || !IsAllowedJ(J2)) return 0.0;
    if (k < 0 || (k % 2) != 0) return 0.0;
    if (L < 0 || Lp < 0) return 0.0;

    // Angular selection for CG <L1 L'-1 | k0>
    if (k < std::abs(L - Lp) || k > (L + Lp)) return 0.0;

    const int tJ1 = detail::to_twoj(J1);
    const int tJ2 = detail::to_twoj(J2);

    // exponent = J1 - J2 - 1 ; must be integer
    const int phase_exp_num = tJ1 - tJ2 - 2; // = 2*(J1-J2-1)
    if (phase_exp_num & 1) return 0.0;
    const int phase = detail::phase_m1_from_int(phase_exp_num / 2);

    const double pref = std::sqrt((2.0*L + 1.0) * (2.0*Lp + 1.0) * (2.0*J1 + 1.0));

    // <L,1 ; L',-1 | k,0>
    const double cg = ClebschGordan((double)L, 1.0, (double)Lp, -1.0, (double)k, 0.0);

    // { J1 J1 k ; L L' J2 }
    const double sixj = Wigner6j(J1, J1, (double)k, (double)L, (double)Lp, J2);

    double val = phase * pref * cg * sixj;
    if (std::fabs(val) < 1e-14) val = 0.0;
    return val;
}

// ============================================================
// Compute the 6 coefficients matching your table column layout
//
// rk01 = F2(1,1)
// rk11 = F2(1,2)
// rk21 = F2(2,2)
// rk02 = F4(1,1)
// rk12 = F4(1,2)
// rk22 = F4(2,2)
// ============================================================
inline bool ComputeRacahCoeffs(double J1, double J2, Coeffs &out)
{
    if (!IsAllowedJ(J1) || !IsAllowedJ(J2)) return false;

    out.rk01 = Fk_candidate(2, 1, 1, J1, J2);
    out.rk11 = Fk_candidate(2, 1, 2, J1, J2);
    out.rk21 = Fk_candidate(2, 2, 2, J1, J2);

    out.rk02 = Fk_candidate(4, 1, 1, J1, J2);
    out.rk12 = Fk_candidate(4, 1, 2, J1, J2);
    out.rk22 = Fk_candidate(4, 2, 2, J1, J2);

    return true;
}

// Convenience overload
inline bool ComputeRacahCoeffs(double J1, double J2,
                               double &rk01, double &rk11, double &rk21,
                               double &rk02, double &rk12, double &rk22)
{
    Coeffs c;
    if (!ComputeRacahCoeffs(J1, J2, c)) return false;
    rk01 = c.rk01; rk11 = c.rk11; rk21 = c.rk21;
    rk02 = c.rk02; rk12 = c.rk12; rk22 = c.rk22;
    return true;
}

// ============================================================
// Utility: absolute-error comparison to a reference row
// ============================================================
inline std::array<double,6> Diff(const Coeffs &calc, const Coeffs &ref)
{
    return {
        calc.rk01 - ref.rk01,
        calc.rk11 - ref.rk11,
        calc.rk21 - ref.rk21,
        calc.rk02 - ref.rk02,
        calc.rk12 - ref.rk12,
        calc.rk22 - ref.rk22
    };
}

inline double MaxAbsDiff(const Coeffs &calc, const Coeffs &ref)
{
    const auto d = Diff(calc, ref);
    double m = 0.0;
    for (double v : d) m = std::max(m, std::fabs(v));
    return m;
}

} // namespace RacahCalc

#endif // RACAH_CALC_H
