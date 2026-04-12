#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include "math_func.h"

const double EPS = 1e-15;
const int MAX_ITER = 2000;

int main() {
    double a, b, c, d, n, m, k;

    // ввод коэффициентов (с проверкой его корректности)
    while (true) {
        std::cout << "Enter coefficients a b c d n m k:\n";
        std::string line;
        std::getline(std::cin, line);
        std::istringstream iss(line);
        double extra;
        if ((iss >> a >> b >> c >> d >> n >> m >> k) && !(iss >> extra)) {
            break;
        }
        std::cout << "Error: Enter 7 numeric values correctly!\n";
    }

    std::vector<double> coeffs = {a, b, c, d, n, m, k};
    std::vector<double> roots;

    normalize(coeffs);

    int deg = coeffs.size() - 1;
    int deg = coeffs.size() - 1;
    
    // нахождение первого корня (комбинированый метод хорд-касателных)
    if (deg == 6) {
    if (deg == 6) {
        double r = combined(coeffs);
        if (std::isfinite(r)) roots.push_back(r);
        coeffs = horner(coeffs, r);
        normalize(coeffs);
    }

    deg = coeffs.size() - 1;
    deg = coeffs.size() - 1;

    // нахождение второго корня (метод хорд)
    if (deg == 5) {
    if (deg == 5) {
        double r = chord(coeffs);
        if (std::isfinite(r)) roots.push_back(r);
        coeffs = horner(coeffs, r);
        normalize(coeffs);
    }

    deg = coeffs.size() - 1;
    deg = coeffs.size() - 1;

    // нахождение третьего корня (метод касательных)
    if (deg == 4) {
    if (deg == 4) {
        double r = newton(coeffs);
        if (std::isfinite(r)) roots.push_back(r);
        coeffs = horner(coeffs, r);
        normalize(coeffs);
    }

    deg = coeffs.size() - 1; 
    deg = coeffs.size() - 1; 

    // нахождение оставшихся корней (метод Кардано)
    if (deg == 3) {
    if (deg == 3) {
        std::vector<double> r = cardano(coeffs);
        for (double x : r)
            if (std::isfinite(x)) roots.push_back(x);
        deg = coeffs.size() - 1;

    // нахождение оставшихся корней (формула Виета)
    } else if (deg == 2) {
    } else if (deg == 2) {
        std::vector<double> r = viet(coeffs);
        for (double x : r)
            if (std::isfinite(x)) roots.push_back(x);
        deg = coeffs.size() - 1;

    // нахождение оставшихся корней (-k/m)
    } else if (deg == 1) {
        double x = line_eq(coeffs);
            if (std::isfinite(x)) roots.push_back(x);
        deg = coeffs.size() - 1;

    // проверка на наличие найденных корней 
    } else if (std::abs(coeffs[3]) > EPS) {
        if (roots.empty()) {
            std::cout << "No roots found.\n";
            return 0;
        }
    } else {
        if (roots.empty()) {
            std::cout << "Any real number is a solution.\n";
            return 0;
        }
    }

    // отсеивание повторяющихся корней (с точностью 0.01)
    roots = unique_roots(roots, 0.01);

    // проверка на наличие корней
    if (roots.empty()) {
        std::cout << "No roots found.\n";
        return 0;
    }

    // вывод корней, если они есть
    std::cout << "\nRoots:\n";
    for (double r : roots) {
        if (std::isfinite(r)) {
            r = round_root(r);
            std::cout << r << std::endl;
        }
    }
    return 0;
}