#ifndef COSTFUNCTION_H
#define COSTFUNCTION_H

#pragma once
#include <vector>
#include <Eigen/Dense>
#include <memory>
#include <cmath>
using namespace Eigen;

struct SimulationResults {
    double R;
    double dt = 0.001;  
};

class CostFunction{
    public:
        virtual ~CostFunction() = default;

        virtual double computeCost(const SimulationResults& results) const = 0;
};

class AeroCoefficientsCost : public CostFunction{
    public:
        double w_R_err;

        AeroCoefficientsCost(double w_R_e = 1.0) : w_R_err(w_R_e) {}

        double computeCost(const SimulationResults& res) const override{
            double cost = w_R_err*res.R;

            return cost;
        }
};

#endif