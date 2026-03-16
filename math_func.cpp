#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include "math_func.h"

const double EPS = 1e-12;
const int MAX_ITER = 1000;

double line_eq(double m, double k) { // cm, ck
    return -k / m;
}

std::vector<double> viet(double a, double b, double c) { // cn, cm, ck 
    std::vector<double> roots;
    double D = b * b - 4 * a * c;

    if (D > EPS) {
        double x1 = (-b + sqrt(D)) / (2 * a);
        double x2 = (-b - sqrt(D)) / (2 * a);
        roots.push_back(x1);
        roots.push_back(x2);	
    } else if (std::abs(D) <= EPS) {
        double x = -b / (2 * a); 
        roots.push_back(x);
    }
    return roots;
}


    std::vector<double> cardano(double a, double b, double c, double d) { // cd, cn, cm, ck
        std::vector<double> roots;

        double p = (3 * a * c - b * b) / (3 * a * a);
        double q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);

        double D = (q * q) / 4.0 + (p * p * p) / 27.0;
        double shift = b / (3.0 * a);

        if (D > EPS) {
            double u = cbrt(-q / 2.0 + std::sqrt(D));
            double v = cbrt(-q / 2.0 - std::sqrt(D));

            double y = u + v;
            roots.push_back(y - shift);
        } else if (std::abs(D) <= EPS) {
            double u = cbrt(-q / 2.0);

            roots.push_back(2 * u - shift);
            roots.push_back(-u - shift);
        } else {
            double r = std::sqrt(-(p*p*p) / 27.0);
            double arg = -q / (2.0 * r);
            arg = std::max(-1.0, std::min(1.0, arg));

            double phi = std::acos(arg);

            double m = 2.0 * std::sqrt(-p / 3.0);
            for (int k = 0; k < 3; ++k) {
                double y = m * std::cos((phi + 2*M_PI*k)/3);
                roots.push_back(y - shift);
            }
        }
        return roots;
    }


std::vector<double> horner(const std::vector<double>& coeffs, double r) {
    int n = coeffs.size();
    std::vector<double> new_coeffs(n - 1);
    new_coeffs[0] = coeffs[0];
    for (int i = 1; i < n - 1; i++)
        new_coeffs[i] = coeffs[i] + new_coeffs[i - 1] * r;
    return new_coeffs;
}


double evaluate(const std::vector<double>& coeffs, double x) {
    double result = 0.0;
    for (double c: coeffs)
        result = result * x + c;
    return result;
}


double firstDerivative(const std::vector<double>& coeffs, double x) {
    double result = 0.0;
    int n = coeffs.size() - 1;
    for (int i = 0; i < n; ++i)
        result = result * x + coeffs[i] * (n - i);
    return result;
}


double secondDerivative(const std::vector<double>& coeffs, double x) {
    double result = 0.0;
    int n = coeffs.size() - 1;  
    for (int i = 0; i < n - 1; ++i)
        result = result * x + coeffs[i] * (n - i) * (n - i - 1);
    return result;
}


double lagrangeBound(const std::vector<double>& coeffs) {
    double A = fabs(coeffs[0]);
    if (A < 1e-15) return 1000.0;
        
    double max_ratio = 0.0;
    for (size_t i = 1; i < coeffs.size(); i++) {
        max_ratio = std::max(max_ratio, fabs(coeffs[i]) / A);
    }
    return 1.0 + max_ratio;
}


bool findInterval(const std::vector<double>& coeffs, double& a, double& b, double& fa, double& fb) {
    double R = lagrangeBound(coeffs);
    double step = R / 5000.0;
    double x1 = -R;
    double f1 = evaluate(coeffs, x1);
        
    for (double x2 = -R + step; x2 <= R; x2 += step) {
        double f2 = evaluate(coeffs, x2);
            
        if (f1 * f2 <= 0) {
            a = x1;
            b = x2;
            fa = f1;
            fb = f2;
            return true;
        }
            
        if (std::abs(f1) < EPS) {
            a = x1 - step;
            b = x1 + step;
            fa = evaluate(coeffs, a);
            fb = evaluate(coeffs, b);
            return true;
        }
            
        x1 = x2;
        f1 = f2;
    }
    return false;
}


double chordMethod(const std::vector<double>& coeffs) {
    double a, b, fa, fb;
    if (!findInterval(coeffs, a, b, fa, fb)) {
        return NAN;
    }
        
    double x;
    double x_prev = a;
        
    for (int iter = 0; iter < MAX_ITER; iter++) {
        if (std::abs(fb - fa) < EPS) {
            x = (a + b) / 2.0;
        } else {
            x = (a * fb - b * fa) / (fb - fa);
        }
            
        double fx = evaluate(coeffs, x);
            
        if (std::abs(fx) < EPS || std::abs(x - x_prev) < EPS) {
            return x;
        }
            
        if (fa * fx <= 0) {
            b = x;
            fb = fx;
        } else {
            a = x;
            fa = fx;
        }    
        x_prev = x;
    }
    return x;
}


double newtonMethod(const std::vector<double>& coeffs) {
    double a, b, fa, fb;
    if (!findInterval(coeffs, a, b, fa, fb)) {
        return NAN;
    }
    double mid = (a + b) / 2.0;
    double f2_mid = secondDerivative(coeffs, mid);
    double xn;
        
    if (fa * f2_mid > 0) {
        xn = a;
    } else if (fb * f2_mid > 0) {
        xn = b;
    } else {
        xn = mid;
    }
        
    double x0;
    int iter = 0;        
    do {
        iter++;
        x0 = xn;    
        double fx = evaluate(coeffs, xn);
        double dfx = firstDerivative(coeffs, xn);
        double x_new;
            
        if (std::abs(dfx) < EPS) {
            x_new = (a + b) / 2.0;
        } else {
            x_new = xn - fx / dfx;
        }
            
        if (x_new < a || x_new > b) {
            x_new = (a + b) / 2.0;
        }
            
        double f_new = evaluate(coeffs, x_new);
        if (fa * f_new <= 0) {
            b = x_new;
            fb = f_new;
        } else {
            a = x_new;
            fa = f_new;
        }
        xn = x_new;    
        if (iter > MAX_ITER) break;    
        } while (std::abs(xn - x0) > EPS && std::abs(evaluate(coeffs, xn)) > EPS);
        
        return xn;
    }


double combinedMethod(const std::vector<double>& coeffs) {
    double a, b, fa, fb;
    if (!findInterval(coeffs, a, b, fa, fb)) {
        return NAN;
    }
        
    double xn = (a + b) / 2.0;
    double x_prev;
        
    for (int iter = 0; iter < MAX_ITER; iter++) {
        double x_chord = (a * fb - b * fa) / (fb - fa);
        double f_chord = evaluate(coeffs, x_chord);
            
        double df_chord = firstDerivative(coeffs, x_chord);
        double x_newton;
            
        if (std::abs(df_chord) < EPS) {
            x_newton = x_chord;
        } else {
            x_newton = x_chord - f_chord / df_chord;
        }
            
        xn = x_newton;
        double fx = evaluate(coeffs, xn);
            
        if (std::abs(fx) < EPS || (iter > 0 && std::abs(xn - x_prev) < EPS)) {
            return xn;
        }
        if (fa * fx <= 0) {
            b = xn;
            fb = fx;
        } else {
            a = xn;
            fa = fx;
        }
        x_prev = xn;
    }
        
    return xn;
}
