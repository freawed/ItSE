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
    
    // нахождение первого корня (комбинированый метод хорд-касателных)
    if (std::abs(coeffs[0]) > EPS) {
        double r = combined(coeffs);
        if (std::isfinite(r)) roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    // нахождение второго корня (метод хорд)
    if (std::abs(coeffs[0]) > EPS) {
        double r = chord(coeffs);
        if (std::isfinite(r)) roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    // нахождение третьего корня (метод касательных)
    if (std::abs(coeffs[0]) > EPS) {
        double r = newton(coeffs);
        if (std::isfinite(r)) roots.push_back(r);
        coeffs = horner(coeffs, r);
    }

    // нахождение оставшихся корней (метод Кардано)
    if (std::abs(coeffs[0]) > EPS) {
        std::vector<double> r = cardano(coeffs);
        for (double x : r)
            if (std::isfinite(x)) roots.push_back(x);
    // нахождение оставшихся корней (формула Виета)
    } else if (std::abs(coeffs[1]) > EPS) {
        std::vector<double> r = viet(coeffs);
        for (double x : r)
            if (std::isfinite(x)) roots.push_back(x);
    // нахождение оставшихся корней (-k/m)
    } else if (std::abs(coeffs[2]) > EPS) {
        double x = line_eq(coeffs);
            if (std::isfinite(x)) roots.push_back(x);
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
        if (std::isfinite(r))
            std::cout << r << std::endl;
    }
    return 0;
}