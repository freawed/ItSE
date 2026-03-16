#ifndef MATH_UTILS_H
#define MATH_UTILS_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>        

double line_eq(double m, double k);
std::vector<double> viet(double a, double b, double c);
std::vector<double> cardano(double a, double b, double c, double d);
std::vector<double> horner(const std::vector<double>& coeffs, double r);
double evaluate(const std::vector<double>& coeffs, double x);
double lagrangeBound(const std::vector<double>& coeffs);
double firstDerivative(const std::vector<double>& coeffs, double x);
double secondDerivative(const std::vector<double>& coeffs, double x);
double newtonMethod(const std::vector<double>& coeffs);
double chordMethod(const std::vector<double>& coeffs);
double combinedMethod(const std::vector<double>& coeffs);
bool findInterval(const std::vector<double>& coeffs, double& a, double& b, double& fa, double& fb);

#endif