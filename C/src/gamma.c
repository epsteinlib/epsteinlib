// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file gamma.c
 * @brief Gamma functions.
 *
 * Calculates the gamma function, the incomplete upper gamma function and the
regularized lower incomplete gamma function for evaluations of Crandall's
formula.
 * @see Walter Gautschi. “A Computational Procedure for
Incomplete Gamma Func-257 tions”. In: ACM Trans. Math. Softw. 5 (1979), pp.
466–481 for basic principles
 * some improvements were made to increase accuracy
 */

#include "gamma.h"
#include <math.h>

/*!
 * @brief epsilon for cutoff around integers.
 */
#define EGF_EPS ldexp(1, -54)
/*!
 * @brief enum for the choice of algorithm for the upper incomplete gamma
 * function.
 */
enum dom { pt, qt, cf, ua, rek };

/**
 * @brief set type of algorithm to use depending on the parameters.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return enum for the type of algorithm to use.
 */
enum dom egf_domain(double a, double x) {
    double alpha;
    if (x >= 0.5) {
        alpha = x;
    } else {
        alpha = log(0.5) / log(0.5 * x);
    }
    if (a <= alpha) {
        if (x <= 1.5 && a >= -0.5) {
            return qt;
        }
        if (x <= 1.5) {
            return rek;
        }
        if (a >= 12 && a >= x / 2.35) {
            return ua;
        }
        return cf;
    }
    if (a >= 12 && x >= 0.3 * a) {
        return ua;
    }
    return pt;
}

enum dom egf_ldomain(double a, double x) {
    double alpha;
    if (x >= 0.5) {
        alpha = x;
    } else {
        alpha = log(0.5) / log(0.5 * x);
    }
    if (a <= alpha) {
        if (x <= 1.5 && (a >= -0.5 || (a >= -0.75 && x <= ldexp(2, -15)))) {
            return pt;
        }
        if (x <= 1.5) {
            return rek;
        }
        if (a >= 12 && a >= x / 2.35) {
            return ua;
        }
        return cf;
    }
    if (a >= 12 && x >= 0.3 * a) {
        return ua;
    }
    return pt;
}

/**
 * @brief calculate upper gamma function with the recursion formula.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_pt(double a, double x) {
    // optional für x >= 10000: Prüfe Restglied
    double sn = 1;
    double add = x / (a + 1);
    for (int i = 1; i < 80 && fabs(add / sn) >= EGF_EPS; i++) {
        sn += add;
        add *= (x / (a + i + 1));
    }
    return sn * exp(-x) / tgamma(a + 1);
}

/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_qt(double a, double x) {
    static double taylor[21] = {
        -0.57721566490153286061,    0.078662406618721020471,
        0.120665041652816256,       -0.045873569729475233502,
        -0.003675835173930896754,   0.0059461363539460768081,
        -0.0012728068927170227343,  -0.00010763930085795762215,
        0.00010760237325699335067,  -0.000020447909131122835485,
        -3.1305435033459682903e-7,  9.3743913180807382831e-7,
        -1.9558810017362205406e-7,  1.0045741524138656286e-8,
        3.9296464196572404677e-9,   -1.0723612248119824624e-9,
        1.0891334567503768218e-10,  4.5706745059276311356e-12,
        -3.2115889339774401184e-12, 4.8521668466476558978e-13,
        -2.4820344080682008122e-14};
    double u;
    if (fabs(a) < 0.5) {
        double u1 = taylor[0];
        double f = 1;
        for (int i = 1; i < 21; i++) {
            f *= a;
            u1 += taylor[i] * f;
        } // u1 = g(a)
        double u2 = 0;
        double y = a * log(x);
        f = 1;
        if (fabs(y) < 1) {
            for (int n = 1; n <= 30; n++) {
                f /= (double)n;
                u2 += f;
                f *= y;
            }
        } else {
            u2 = (exp(y) - 1) / y;
        }
        u = tgamma(1 + a) * (1 - a) * u1 - u2 * log(x);
    } else {
        u = tgamma(a) - pow(x, a) / a;
    }
    double v = 0;
    double f = 1;
    for (int i = 1; i <= 30; i++) {
        f *= (-1) * x / (double)i;
        v += f / (double)(a + i);
    }
    v *= -pow(x, a);
    return u + v;
}

/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_rek(double a, double x) {
    int m = (int)(0.5 - a);
    double epsilon = a + m;
    double g = egf_qt(epsilon, x) * exp(x) * pow(x, -epsilon);
    for (int n = 1; n <= m; n++) {
        g = 1. / (n - epsilon) * (1. - x * g);
    }
    return g;
}

/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_cf(double a, double x) {
    double s = 1;
    double rp = 1; // t_k-1
    double rv = 0; // rho_0
    for (int k = 1; k <= 200 && fabs(rp / s) >= EGF_EPS; k++) {
        double ak =
            k * (a - k) / (double)((x + 2 * k - 1 - a) * (x + 2 * k + 1 - a));
        rv = -ak * (1 + rv) / (1 + ak * (1 + rv));
        rp *= rv;
        s += rp;
    }
    double r = s * pow(x, a) * exp(-x) / (x + 1 - a);
    return r;
}

/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_ua_r(double a, double eta) {
    static double d[27] = {1.0,
                           -1.0 / 3.0,
                           1.0 / 12.0,
                           -2.0 / 135.0,
                           1.0 / 864.0,
                           1.0 / 2835.0,
                           -139.0 / 777600.0,
                           1.0 / 25515.0,
                           -571.0 / 261273600.0,
                           -281.0 / 151559100.0,
                           8.29671134095308601e-7,
                           -1.76659527368260793e-7,
                           6.70785354340149857e-9,
                           1.02618097842403080e-8,
                           -4.38203601845335319e-9,
                           9.14769958223679023e-10,
                           -2.55141939949462497e-11,
                           -5.83077213255042507e-11,
                           2.43619480206674162e-11,
                           -5.02766928011417559e-12,
                           1.10043920319561347e-13,
                           3.37176326240098538e-13,
                           -1.39238872241816207e-13,
                           2.85348938070474432e-14,
                           -5.13911183424257258e-16,
                           -1.97522882943494428e-15,
                           8.09952115670456133e-16};
    double beta[26];
    beta[25] = d[26];
    beta[24] = d[25];
    for (int n = 23; n >= 0; n--) {
        beta[n] = (double)(n + 2) * beta[n + 2] / a + d[n + 1];
    }
    double s = 0;
    double f = 1.;
    for (int i = 0; i <= 25; i++) {
        s += beta[i] * f;
        f *= eta;
    }
    s *= a / (a + beta[1]);
    return s * exp(-0.5 * a * eta * eta) / sqrt(2 * M_PI * a);
}

/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_ua(double a, double x) {
    double lambda = x / a;
    double eta = sqrt(2 * (lambda - 1 - log(lambda)));
    if (lambda - 1 < 0) {
        eta = -eta;
    }
    double ra = egf_ua_r(a, eta);
    return 0.5 * erfc(eta * sqrt(a / 2.)) + ra;
}

/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_ugamma(double a, double x) {
    double r = NAN;
    enum dom g = egf_domain(a, x);
    switch (g) {
    case pt:
        r = tgamma(a) * (1 - egf_pt(a, x) * pow(x, a));
        break;
    case qt:
        r = egf_qt(a, x);
        break;
    case cf:
        r = egf_cf(a, x);
        break;
    case ua:
        r = tgamma(a) * egf_ua(a, x);
        break;
    case rek:
        r = exp(-x) * pow(x, a) * egf_rek(a, x);
        break;
    }
    return r;
}

/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param[in] a: exponent of the upper incomplete gamma function.
 * @param[in] x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_gammaStar(double a, double x) {
    double r = NAN;
    if (fabs(x) < EGF_EPS) {
        if (a <= 0.1 && fabs(a - nearbyint(a)) < EGF_EPS) {
            return 0;
        }
        return 1. / tgamma(a + 1);
    }
    enum dom g = egf_ldomain(a, x);
    switch (g) {
    case pt:
    case qt:
        r = egf_pt(a, x);
        break;
    case cf:
        if (a <= 0.1 && fabs(a - nearbyint(a)) < EGF_EPS) {
            r = pow(x, -a);
        } else {
            r = (1 - egf_cf(a, x) / tgamma(a)) * pow(x, -a);
        }
        break;
    case ua:
        r = (1 - egf_ua(a, x)) * pow(x, -a);
        break;
    case rek:
        if (a <= 0.1 && fabs(a - nearbyint(a)) < EGF_EPS) {
            r = pow(x, -a);
        } else {
            r = (1. - exp(-x) * pow(x, a) * egf_rek(a, x) / tgamma(a)) * pow(x, -a);
        }
        break;
    }
    return r;
}
#undef EGF_EPS
