#include "normal.h"

using namespace std;

using boost::math::normal;
using boost::math::binomial;

double qnorm(double prob, double mean, double sd) {
    normal Normal(mean, sd);
    return quantile(Normal, prob);
}

double dnorm(double x, double mean, double sd) {
    normal  Normal(mean, sd);
    return pdf(Normal, x);
}


double mean(vector<double> num) {
    double n = num.size(), sum = 0;
    for (int i = 0; i < n; i++) {
        sum += num[i];
    }
    return(sum / n);
}

double var(vector<double> num) {
    double num_m = mean(num), sum = 0, n = num.size();
    for (int i = 0; i < n; i++) {
        sum += (num[i] - num_m) * (num[i] - num_m);
    }
    return(sum / (n - 1));
}

double update_sigma(vector<double> num, double sigma_start) {
    double n = num.size();
    return(sqrt(var(num) * (n + 1) / (n + 0.5)) + sigma_start / n);
}

double update_prob(vector<double> m, double n) {
    int m_size = m.size();
    double sum = 0;
    for (int i = 0; i < m_size; i++) {
        sum += m[i];
    }
    return(sum / (n * m_size));
}