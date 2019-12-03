/*
This c++ code for calculation of Gauss-Legendre quadratures and weights was obtained from the following website with minor bug fixes.
https://thoughts-on-cpp.com/2019/04/25/numerical-methods-in-c-part-2-gauss-legendre-integration/
*/

#ifndef LEGENDREPOLYNOMIAL_H
#define LEGENDREPOLYNOMIAL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

class LegendrePolynomial {
public:
    LegendrePolynomial(double lowerBound, double upperBound, unsigned int numberOfIterations)
        : mLowerBound(lowerBound), mUpperBound(upperBound), mNumberOfIterations(numberOfIterations), mWeight(numberOfIterations), mRoot(numberOfIterations) {
        calculateWeightAndRoot();
    }

    const std::vector<double> & getWeight() const {
        return mWeight;
    }

    const std::vector<double> & getRoot() const {
        return mRoot;
    }

private:
    const double EPSILON = 1e-15;
    struct Result {
        double value;
        double derivative;

        Result() : value(0), derivative(0) {}
        Result(double val, double deriv) : value(val), derivative(deriv) {}
    };

    void calculateWeightAndRoot() {
        for(int step = 0; step <= mNumberOfIterations; step++) {
            double root = cos(M_PI * (step-0.25)/(mNumberOfIterations+0.5));
            Result result = calculatePolynomialValueAndDerivative(root);

            double newtonRaphsonRatio;
            do {
                newtonRaphsonRatio = result.value/result.derivative;
                root -= newtonRaphsonRatio;
                result = calculatePolynomialValueAndDerivative(root);
            } while (fabs(newtonRaphsonRatio) > EPSILON);
            if(step>=1){
              mRoot[step-1] = root;
              mWeight[step-1] = 2.0/((1-root*root)*result.derivative*result.derivative);
            }
        }
        std::reverse(mRoot.begin(), mRoot.end());
        std::reverse(mWeight.begin(), mWeight.end());
    }

    Result calculatePolynomialValueAndDerivative(double x) {
        Result result(x, 0);

        double value_minus_1 = 1;
        const double f = 1/(x*x-1);
        for(int step = 2; step <= mNumberOfIterations; step++) {
            const double value = ((2*step-1)*x*result.value-(step-1)*value_minus_1)/step;
            result.derivative = step*f*(x*value - result.value);

            value_minus_1 = result.value;
            result.value = value;
        }
        return result;
    }

    const double mLowerBound;
    const double mUpperBound;
    const int mNumberOfIterations;
    std::vector<double> mWeight;
    std::vector<double> mRoot;
};

#endif
