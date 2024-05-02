#include "Optimizer.h"

#include <cmath>

SimpleConjugateGradient::SimpleConjugateGradient(BaseFunction &obj,
                                                 std::vector<Point2<double>> &var,
                                                 const double &alpha)
    : BaseOptimizer(obj, var),
      grad_prev_(var.size()),
      dir_prev_(var.size()),
      step_(0),
      alpha_(alpha) {}

void SimpleConjugateGradient::Initialize() {
    // Before the optimization starts, we need to initialize the optimizer.
    step_ = 0;
}

/**
 * @details Update the solution once using the conjugate gradient method.
 */
void SimpleConjugateGradient::Step() {
    const size_t &kNumModule = var_.size();

    // Compute the gradient direction
    obj_(var_);       // Forward, compute the function value and cache from the input
    obj_.Backward();  // Backward, compute the gradient according to the cache

    // Compute the Polak-Ribiere coefficient and conjugate directions
    double beta;                                  // Polak-Ribiere coefficient
    std::vector<Point2<double>> dir(kNumModule);  // conjugate directions
    if (step_ == 0) {
        // For the first step, we will set beta = 0 and d_0 = -g_0
        beta = 0.;
        for (size_t i = 0; i < kNumModule; ++i) {
            dir[i] = -obj_.grad().at(i);
        }
    } else {
        // For the remaining steps, we will calculate the Polak-Ribiere coefficient and
        // conjugate directions normally
        double t1 = 0.;  // Store the numerator of beta
        double t2 = 0.;  // Store the denominator of beta
        for (size_t i = 0; i < kNumModule; ++i) {
            Point2<double> t3 =
                obj_.grad().at(i) * (obj_.grad().at(i) - grad_prev_.at(i));
            t1 += t3.x + t3.y;
            t2 += std::abs(obj_.grad().at(i).x) + std::abs(obj_.grad().at(i).y);
        }
        beta = t1 / (t2 * t2);
        for (size_t i = 0; i < kNumModule; ++i) {
            dir[i] = -obj_.grad().at(i) + beta * dir_prev_.at(i);
        }
    }

    // Assume the step size is constant
    // TODO(Optional): Change to dynamic step-size control

    // Update the solution
    // Please be aware of the updating directions, i.e., the sign for each term.
    for (size_t i = 0; i < kNumModule; ++i) {
        var_[i] = var_[i] + alpha_ * dir[i];
    }

    // Update the cache data members
    grad_prev_ = obj_.grad();
    dir_prev_ = dir;
    step_++;
}
