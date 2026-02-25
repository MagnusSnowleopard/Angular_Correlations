#ifndef RACAH_H
#define RACAH_H

#include <cmath>
#include <array>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace RACAH_COEFFS {

struct RKSet {
    double rk01; // k=2, (L,L)
    double rk11; // k=2, (L,L+1)
    double rk21; // k=2, (L+1,L+1)
    double rk02; // k=4, (L,L)
    double rk12; // k=4, (L,L+1)
    double rk22; // k=4, (L+1,L+1)
    int Lbase;   // L = max(1, |J1-J2|)
};

inline bool nearly_int(double x, double eps = 1e-9) {
    return std::fabs(x - std::llround(x)) < eps;
}

inline int twoj(double x) {
    double tx = 2.0 * x;
    if (!nearly_int(tx)) {
        throw std::runtime_error("Value is not integer/half-integer.");
    }
    return static_cast<int>(std::llround(tx));
}

inline double fact_gamma(double x) {
    // x! = Gamma(x+1), valid for x = integer or half-integer (and x >= 0)
    if (x < -1e-12) return 0.0;
    return std::exp(std::lgamma(x + 1.0));
}

inline double phase_int(long long n) {
    return (n % 2LL == 0LL) ? 1.0 : -1.0;
}

inline bool triangle(double a, double b, double c, double eps = 1e-10) {
    if (a < -eps || b < -eps || c < -eps) return false;
    if ((a + b) < c - eps) return false;
    if ((a + c) < b - eps) return false;
    if ((b + c) < a - eps) return false;
    // parity condition: a+b+c integer
    return nearly_int(a + b + c);
}

inline double delta_tri(double a, double b, double c) {
    if (!triangle(a, b, c)) return 0.0;

    double n1 = a + b - c;
    double n2 = a - b + c;
    double n3 = -a + b + c;
    double d  = a + b + c + 1.0;

    if (n1 < -1e-12 || n2 < -1e-12 || n3 < -1e-12) return 0.0;

    double logv = std::lgamma(n1 + 1.0) + std::lgamma(n2 + 1.0) + std::lgamma(n3 + 1.0)
                - std::lgamma(d  + 1.0);

    return std::exp(0.5 * logv);
}

inline double wigner_3j(double j1, double j2, double j3,
                        double m1, double m2, double m3) {
    // Selection rules
    if (!nearly_int(j1) && !nearly_int(j1 - 0.5)) {} // harmless
    if (!nearly_int(m1) && !nearly_int(m1 - 0.5)) {} // harmless

    if (!nearly_int(m1 + m2 + m3)) return 0.0;
    if (std::fabs(m1 + m2 + m3) > 1e-10) return 0.0;

    if (!triangle(j1, j2, j3)) return 0.0;

    if (std::fabs(m1) > j1 + 1e-10) return 0.0;
    if (std::fabs(m2) > j2 + 1e-10) return 0.0;
    if (std::fabs(m3) > j3 + 1e-10) return 0.0;

    if (!nearly_int(j1 + m1) || !nearly_int(j1 - m1) ||
        !nearly_int(j2 + m2) || !nearly_int(j2 - m2) ||
        !nearly_int(j3 + m3) || !nearly_int(j3 - m3)) {
        return 0.0;
    }

    double pref_phase = phase_int(static_cast<long long>(std::llround(j1 - j2 - m3)));
    double pref = pref_phase * delta_tri(j1, j2, j3);

    double logfac = 0.5 * (
        std::lgamma(j1 + m1 + 1.0) + std::lgamma(j1 - m1 + 1.0) +
        std::lgamma(j2 + m2 + 1.0) + std::lgamma(j2 - m2 + 1.0) +
        std::lgamma(j3 + m3 + 1.0) + std::lgamma(j3 - m3 + 1.0)
    );

    pref *= std::exp(logfac);

    // Sum over z
    // Standard Racah formula bounds
    double zmin_d = std::max({
        0.0,
        j2 - j3 - m1,
        j1 - j3 + m2
    });

    double zmax_d = std::min({
        j1 + j2 - j3,
        j1 - m1,
        j2 + m2
    });

    int zmin = static_cast<int>(std::ceil(zmin_d - 1e-12));
    int zmax = static_cast<int>(std::floor(zmax_d + 1e-12));

    if (zmin > zmax) return 0.0;

    double sum = 0.0;
    for (int z = zmin; z <= zmax; ++z) {
        double a1 = z;
        double a2 = j1 + j2 - j3 - z;
        double a3 = j1 - m1 - z;
        double a4 = j2 + m2 - z;
        double a5 = j3 - j2 + m1 + z;
        double a6 = j3 - j1 - m2 + z;

        if (a1 < -1e-12 || a2 < -1e-12 || a3 < -1e-12 ||
            a4 < -1e-12 || a5 < -1e-12 || a6 < -1e-12) continue;

        double term_log = -(
            std::lgamma(a1 + 1.0) + std::lgamma(a2 + 1.0) + std::lgamma(a3 + 1.0) +
            std::lgamma(a4 + 1.0) + std::lgamma(a5 + 1.0) + std::lgamma(a6 + 1.0)
        );

        double term = phase_int(z) * std::exp(term_log);
        sum += term;
    }

    return pref * sum;
}

inline double clebsch_gordan(double j1, double m1,
                             double j2, double m2,
                             double J,  double M) {
    // <j1 m1 j2 m2 | J M>
    if (std::fabs((m1 + m2) - M) > 1e-10) return 0.0;

    double threej = wigner_3j(j1, j2, J, m1, m2, -M);
    if (std::fabs(threej) < 1e-300) return 0.0;

    double phase = phase_int(static_cast<long long>(std::llround(j1 - j2 + M)));
    return phase * std::sqrt(2.0 * J + 1.0) * threej;
}

inline double wigner_6j(double a, double b, double c,
                        double d, double e, double f) {
    // Selection rules (triangle conditions for each triad)
    if (!triangle(a, b, c)) return 0.0;
    if (!triangle(a, e, f)) return 0.0;
    if (!triangle(d, b, f)) return 0.0;
    if (!triangle(d, e, c)) return 0.0;

    double d1 = delta_tri(a, b, c);
    double d2 = delta_tri(a, e, f);
    double d3 = delta_tri(d, b, f);
    double d4 = delta_tri(d, e, c);

    if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0 || d4 == 0.0) return 0.0;

    double pref = d1 * d2 * d3 * d4;

    // Racah sum bounds
    double x1 = a + b + c;
    double x2 = a + e + f;
    double x3 = d + b + f;
    double x4 = d + e + c;

    double y1 = a + b + d + e;
    double y2 = a + c + d + f;
    double y3 = b + c + e + f;

    int zmin = static_cast<int>(std::ceil(std::max({x1, x2, x3, x4}) - 1e-12));
    int zmax = static_cast<int>(std::floor(std::min({y1, y2, y3}) + 1e-12));

    if (zmin > zmax) return 0.0;

    double sum = 0.0;
    for (int z = zmin; z <= zmax; ++z) {
        double t1 = z - x1;
        double t2 = z - x2;
        double t3 = z - x3;
        double t4 = z - x4;
        double t5 = y1 - z;
        double t6 = y2 - z;
        double t7 = y3 - z;

        if (t1 < -1e-12 || t2 < -1e-12 || t3 < -1e-12 || t4 < -1e-12 ||
            t5 < -1e-12 || t6 < -1e-12 || t7 < -1e-12) continue;

        double term_log =
            std::lgamma(z + 2.0) // (z+1)!
          - ( std::lgamma(t1 + 1.0) + std::lgamma(t2 + 1.0) + std::lgamma(t3 + 1.0) +
              std::lgamma(t4 + 1.0) + std::lgamma(t5 + 1.0) + std::lgamma(t6 + 1.0) +
              std::lgamma(t7 + 1.0) );

        double term = phase_int(z) * std::exp(term_log);
        sum += term;
    }

    return pref * sum;
}

inline double angular_rk_component(double J1, double J2, int k, int L, int Lp) {
    // Table convention matching your posted Rose/Brink-style rows:
    //
    // R_k = (-1)^(J1 + J2 - 1) * sqrt((2L+1)(2L'+1)(2J1+1))
    //       * <L 1 L' -1 | k 0> * { J1 J1 k ; L L' J2 }
    //
    // NOTE: This matches the sign convention in your table for tested rows.

    // quick selection rules
    if (k < 0) return 0.0;
    if (!triangle(static_cast<double>(L), static_cast<double>(Lp), static_cast<double>(k))) return 0.0;
    if (!triangle(J1, static_cast<double>(Lp), J2)) return 0.0;
    if (!triangle(static_cast<double>(L), J1, J2)) return 0.0;
    if (!triangle(J1, J1, static_cast<double>(k))) return 0.0;

    double phase = phase_int(static_cast<long long>(std::llround(J1 + J2 - 1.0)));

    double pref = std::sqrt((2.0 * L + 1.0) * (2.0 * Lp + 1.0) * (2.0 * J1 + 1.0));

    double cg = clebsch_gordan(static_cast<double>(L),  1.0,
                               static_cast<double>(Lp), -1.0,
                               static_cast<double>(k),  0.0);

    double sixj = wigner_6j(J1, J1, static_cast<double>(k),
                            static_cast<double>(L), static_cast<double>(Lp), J2);

    double val = phase * pref * cg * sixj;

    // avoid -0.0000 noise
    if (std::fabs(val) < 1e-14) val = 0.0;
    return val;
}

inline RKSet ComputeRK(double J1, double J2) {
    // J values must be integer or half-integer
    (void)twoj(J1);
    (void)twoj(J2);

    // For gamma transitions, base multipole L = max(1, |J1-J2|)
    // (L=0 monopole gamma is not allowed; 0->0 not gamma)
    double dJ = std::fabs(J1 - J2);
    if (!nearly_int(dJ)) {
        throw std::runtime_error("|J1-J2| must be an integer for electromagnetic transitions.");
    }

    int L = std::max(1, static_cast<int>(std::llround(dJ)));

    RKSet out{};
    out.Lbase = L;

    // k=2 group
    out.rk01 = angular_rk_component(J1, J2, 2, L,   L);
    out.rk11 = angular_rk_component(J1, J2, 2, L,   L+1);
    out.rk21 = angular_rk_component(J1, J2, 2, L+1, L+1);

    // k=4 group
    out.rk02 = angular_rk_component(J1, J2, 4, L,   L);
    out.rk12 = angular_rk_component(J1, J2, 4, L,   L+1);
    out.rk22 = angular_rk_component(J1, J2, 4, L+1, L+1);

    return out;
}

inline void PrintRK(double J1, double J2, std::ostream& os = std::cout) {
    RKSet r = ComputeRK(J1, J2);
    os << std::fixed << std::setprecision(4);
    os << "J1=" << J1 << " J2=" << J2 << "  (Lbase=" << r.Lbase << ")\n";
    os << "calc rk01 = " << std::setw(8) << r.rk01 << "\n";
    os << "calc rk11 = " << std::setw(8) << r.rk11 << "\n";
    os << "calc rk21 = " << std::setw(8) << r.rk21 << "\n";
    os << "calc rk02 = " << std::setw(8) << r.rk02 << "\n";
    os << "calc rk12 = " << std::setw(8) << r.rk12 << "\n";
    os << "calc rk22 = " << std::setw(8) << r.rk22 << "\n";
}

} // namespace RACAH_COEFFS

#endif // RACAH_H
