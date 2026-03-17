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
double f(const std::vector<double>& coeffs, double x);
double lagrangeBound(const std::vector<double>& coeffs);
double f1(const std::vector<double>& coeffs, double x);
double f2(const std::vector<double>& coeffs, double x);
double newton(const std::vector<double>& coeffs);
double chord(const std::vector<double>& coeffs);
double combined(const std::vector<double>& coeffs);
bool findInterval(const std::vector<double>& coeffs, double& a, double& b, double& fa, double& fb);

#endif