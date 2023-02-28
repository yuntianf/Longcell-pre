#pragma once
#include <vector>
#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions/binomial.hpp>

using namespace std;


double qnorm(double prob, double mean, double sd);
double dnorm(double x, double mean, double sd);
double mean(vector<double> num);
double var(vector<double> num);
double update_sigma(vector<double> num, double sigma_start);
double update_prob(vector<double> m, double n);