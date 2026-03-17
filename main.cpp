#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include "math_func.h"

const double EPS = 1e-12;
const int MAX_ITER = 1000;

int main() {
    double a, b, c, d, n, m, k;

    while (true) {
        std::cout << "Enter coefficients a b c d n m k:\n";

        std::string line;
        std::getline(std::cin, line);

        std::istringstream iss(line);
        double extra;

        if ((iss >> a >> b >> c >> d >> n >> m >> k) && !(iss >> extra)) {
            break; 
        }

        std::cout << "Error: Enter 7 numeric values ​​correctly!\n";
    }

    std::vector<double> coeffs = {a, b, c, d, n, m, k};
    std::vector<double> roots;

    if (std::abs(coeffs[0]) > EPS) {
        double r = combined(coeffs);
        if (!std::isnan(r) && std::isfinite(r))
            roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    if (std::abs(coeffs[0]) > EPS) {
        double r = chord(coeffs);
        if (!std::isnan(r) && std::isfinite(r))
            roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    if (std::abs(coeffs[0]) > EPS) {
        double r = newton(coeffs);
        if (!std::isnan(r) && std::isfinite(r))
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
            for (double x : r) {
                if (!std::isnan(x) && std::isfinite(x))
                     roots.push_back(x);
            }
        }
        else if (std::abs(B) > EPS) {
            std::vector<double> r = viet(B, C, D);
            for (double x : r) {
                if (!std::isnan(x) && std::isfinite(x))
                     roots.push_back(x);
            }
        }
        else if (std::abs(C) > EPS) {
            double x = line_eq(C, D);
            if (!std::isnan(x) && std::isfinite(x))
                roots.push_back(x);
        }
    }

    std::cout << "\nRoots:\n";
    for (double r : roots) {
        if (!std::isnan(r) && std::isfinite(r))
            std::cout << r << std::endl;
    }


    return 0;
}




