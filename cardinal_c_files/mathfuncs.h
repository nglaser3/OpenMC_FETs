#pragma once
#include <cmath>
#include <vector>
double factorial(int input)
{
    double value = 1;
    for (int i = 1; i <= input; i++)
    {
        value *= i;
    }
    return value;
}

double binomial_coeff(int n, int k)
{
    return factorial(n) / factorial(k) / factorial(n-k);
}


double calculateLegendre(int order, double _x,std::vector<double> coeffs)
{
    std::vector<double> results{1*coeffs[0], _x * coeffs[1]};
    double value;
    for (int n = 1; n <= order; n++)
    {
        results.push_back
        (
            coeffs[n+1]*(((2*n+1) * _x * results[n] - n * results[n-1])/(n+1))
        );
    }
    return results[order];

}

double calculateZernike(int n,int m, double _r,double _theta,double coeff)
{
    double angular{n%2==0 ? std::cos(m*_theta) : std::sin(m*_theta)};
    double radial{0.0};
    for (int k = 0; k < (n-m)/2; k++)
    {
        radial += std::pow(-1,k) * factorial(n-k) / (factorial(k) 
        * factorial((n+m)/2 -k) * factorial((n-m)/2 - k)) * std::pow(_r,n-2*k);
    }
    return coeff*radial*angular;
}