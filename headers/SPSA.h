#ifndef SPSA_H
#define SPSA_H

#pragma once
#include <functional>
#include <random>
#include <cmath>
#include <csv.h>
#include <algorithm>
#include <numeric>

class SPSA{
    public:
        struct Config{
            double a = 0.1892;
            double c = 0.05;
            double A = 10.0;
            double alpha = .8;
            double gamma = .1;

            int max_iters = 100;
            bool verbose = true;

            int ave_runs = 5;

            bool clip_bounds = false;
            std::vector<double> lower_bounds;
            std::vector<double> upper_bounds;

            bool log_to_file = false;
            std::string log_file = "C:\\Software Development\\6 DOF Sim\\Data\\SPSA_data.csv";
        };

        SPSA(size_t dim,const Config& cfg)
            : dim_(dim),cfg_(cfg),k_(0),rng_(std::random_device{}()),dist_(0.0,1.0),csv_(cfg_.log_file)
        {
            if (cfg_.log_to_file){
                std::string data = "iter,cost";
                for (size_t i=0; i < dim_; i++){
                    data += ",theta"+std::to_string(i);
                }
                csv_.write_line(data);
            }
        }

        std::vector<double> optimize(std::vector<double> theta,
                    const std::function<double(const std::vector<double>&)>& loss_fn)
                    {
                        double best_cost = std::numeric_limits<double>::infinity();
                        std::vector<double> best_theta = theta;

                        for (int iter = 0; iter < cfg_.max_iters; ++iter) {
                            theta = step(theta, loss_fn, iter);
                            double cost = averagedCost(theta, loss_fn);

                            if (cost < best_cost) {
                                best_cost = cost;
                                best_theta = theta;
                            }
                        }

                        if (cfg_.verbose) {
                            std::cout << "\nBest cost: " << best_cost << std::endl;
                        }

                        return best_theta;
                    }

    private:
        size_t dim_;
        Config cfg_;
        int k_;
        std::mt19937 rng_;
        std::uniform_real_distribution<double> dist_;
        CSV_Writer csv_;

        std::vector<double> step(const std::vector<double>& theta,
                const std::function<double(const std::vector<double>&)>& loss_fn,
                int iter)
                {
                    std::vector<double> delta(dim_);
                    for (size_t i = 0; i <dim_; ++i){
                        delta[i] = (dist_(rng_) < 0.5) ? -1.0 : 1.0;
                    }

                    double ak = cfg_.a / std::pow(cfg_.A+k_+1,cfg_.alpha);
                    double ck = cfg_.c / std::pow(k_+1,cfg_.gamma);

                    std::vector<double> theta_plus(dim_), theta_minus(dim_);
                    for (size_t i = 0; i < dim_; ++i){
                        theta_plus[i] = theta[i]+ck*delta[i];
                        theta_minus[i] = theta[i]-ck*delta[i];
                    }

                    double y_plus = averagedCost(theta_plus,loss_fn);
                    double y_minus = averagedCost(theta_minus,loss_fn);

                    std::vector<double> ghat(dim_);
                    for (size_t i = 0; i <dim_; ++i){
                        ghat[i] = (y_plus-y_minus)/(2.0*ck*delta[i]);
                    }

                    std::vector<double> theta_next(dim_);
                    for (size_t i = 0; i <dim_; ++i){
                        theta_next[i] = theta[i]-ak*ghat[i];
                        if (cfg_.clip_bounds && cfg_.lower_bounds.size() == dim_ && cfg_.upper_bounds.size() == dim_){
                            theta_next[i] = std::clamp(theta_next[i],cfg_.lower_bounds[i],cfg_.upper_bounds[i]);
                        }
                    }

                    double cost = averagedCost(theta_next,loss_fn);

                    if (cfg_.verbose){
                        std::cout << "Iter " << iter
                                  << " | Cost(avg): " << cost
                                  << " | Params: ";
                        for (auto t : theta_next){
                            std::cout << t << " ";
                        }
                        std::cout << std::endl;
                    }

                    if (cfg_.log_to_file){
                        std::string data = std::to_string(iter)+","+std::to_string(cost);
                        for (auto t : theta_next){
                            data += ","+std::to_string(t);
                        }
                        csv_.write_line(data);
                    }

                    k_++;
                    return theta_next;
                }

        double averagedCost(const std::vector<double>& theta,
                        const std::function<double(const std::vector<double>&)>& loss_fn)
                {
                    double sum = 0.0;
                    for (int i = 0; i <cfg_.ave_runs; i++){
                        sum += loss_fn(theta);
                    }

                    return sum/cfg_.ave_runs;
                }
};



#endif