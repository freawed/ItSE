#ifndef MATH_UTILS_H
#define MATH_UTILS_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>        

double round_root(double x);
double line_eq(const std::vector<double>& coeffs);
std::vector<double> viet(const std::vector<double>& coeffs);
std::vector<double> cardano(const std::vector<double>& coeffs);
std::vector<double> horner(const std::vector<double>& coeffs, double r);
std::vector<double> unique_roots(std::vector<double> roots, double eps);
double f(const std::vector<double>& coeffs, double x);
double lagrangeBound(const std::vector<double>& coeffs);
double df(const std::vector<double>& coeffs, double x);
double ddf(const std::vector<double>& coeffs, double x);
double newton(const std::vector<double>& coeffs);
double chord(const std::vector<double>& coeffs);
double combined(const std::vector<double>& coeffs);
bool findInterval(const std::vector<double>& coeffs, double& a, double& b, double& fa, double& fb);

#endif