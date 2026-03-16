#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include "math_func.h"

const double EPS = 1e-12;
const int MAX_ITER = 1000;

int main() {
    double a, b, c, d, n, m, k;

    std::cout << "Enter coefficients a b c d n m k:\n";
    std::cin >> a >> b >> c >> d >> n >> m >> k;

    std::vector<double> coeffs = {a, b, c, d, n, m, k};
    std::vector<double> roots;

    if (std::abs(coeffs[0]) > EPS) {
        double r = combinedMethod(coeffs);
        roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    if (std::abs(coeffs[0]) > EPS) {
        double r = chordMethod(coeffs);
        roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    if (std::abs(coeffs[0]) > EPS) {
        double r = newtonMethod(coeffs);
        roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    if (coeffs.size() == 4) {

        double A = coeffs[0];
        double B = coeffs[1];
        double C = coeffs[2];
        double D = coeffs[3];

        if (std::abs(A) > EPS) {
            std::vector<double> r = cardano(A, B, C, D);
            roots.insert(roots.end(), r.begin(), r.end());
        }
        else if (std::abs(B) > EPS) {
            std::vector<double> r = viet(B, C, D);
            roots.insert(roots.end(), r.begin(), r.end());
        }
        else if (std::abs(C) > EPS) {
            roots.push_back(line_eq(C, D));
        }
    }

    std::cout << "\nRoots:\n";
    for (double r : roots) {
    if (!std::isnan(r) && std::isfinite(r))
        std::cout << r << std::endl;
    }

    return 0;
}




