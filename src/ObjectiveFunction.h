#define _GLIBCXX_USE_CXX11_ABI 0  // Align the ABI version to avoid compatibility issues with `Placment.h`
#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H


#include <vector>

#include "Placement.h"
#include "Point.h"

/**
 * @brief Base class for objective functions
 */
class BaseFunction {
   public:
    /////////////////////////////////
    // Conssutructors
    /////////////////////////////////

    BaseFunction(const size_t &input_size) : grad_(input_size) {}

    /////////////////////////////////
    // Accessors
    /////////////////////////////////

    const std::vector<Point2<double>> &grad() const { return grad_; }
    const double &value() const { return value_; }

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    // Forward pass, compute the value of the function
    virtual const double &operator()(const std::vector<Point2<double>> &input) = 0;

    // Backward pass, compute the gradient of the function
    virtual const std::vector<Point2<double>> &Backward() = 0;

   protected:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> grad_;  // Gradient of the function
    double value_;                      // Value of the function
};

/**
 * @brief Example function for optimization
 *
 * This is a simple example function for optimization. The function is defined as:
 *      f(t) = 3*t.x^2 + 2*t.x*t.y + 2*t.y^2 + 7
 */
class ExampleFunction : public BaseFunction {
   public:
    /////////////////////////////////
    // Constructors
    /////////////////////////////////

    ExampleFunction(Placement &placement);

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

   private:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> input_;  // Cache the input for backward pass
    Placement &placement_;
};

/**
 * @brief Wirelength function
 */
class Wirelength : public BaseFunction {
    // TODO: Implement the wirelength function, add necessary data members for caching
   public:
    /////////////////////////////////
    // Methods
    Wirelength(Placement &placement);
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

    private:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> input_;  // Cache the input for backward pass
    Placement &placement_;
};

/**
 * @brief Density function
 */
class Density : public BaseFunction {
    // TODO: Implement the density function, add necessary data members for caching
   public:
    /////////////////////////////////
    // Methods
    Density(Placement &placement);
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;
    const double getOverflowRatio() const { return overflowRatio_; }
    /////////////////////////////////
    // Data members
    /////////////////////////////////
    

    std::vector<Point2<double>> input_;  // Cache the input for backward pass
    Placement &placement_;
    int binSize_;
    int widthBinNum_;
    int heightBinNum_;
    double Mb ;
    double overflowRatio_;
    vector<vector<double>> binDensity;
    
};

/**
 * @brief Objective function for global placement
 */
class ObjectiveFunction : public BaseFunction {
    // TODO: Implement the objective function for global placement, add necessary data
    // members for caching
    //
    // Hint: The objetive function of global placement is as follows:
    //       f(t) = wirelength(t) + lambda * density(t),
    // where t is the positions of the modules, and lambda is the penalty weight.
    // You may need an interface to update the penalty weight (lambda) dynamically.
   public:
    /////////////////////////////////
    // Methods
    ObjectiveFunction(Placement &placement)
        : BaseFunction(placement.numModules()), placement_(placement),  wirelength_(placement), density_(placement),iterNum_(0),lambda(0),spreadEnough_(0),currentvalue(0),lastvalue(0)
    {

    }

    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;
    const double getOverflowRatio() const { return density_.getOverflowRatio(); }

   private:
    Placement &placement_;
    Wirelength wirelength_;
    Density density_;
    int iterNum_;
    double lambda;
    bool spreadEnough_;
    double currentvalue;
    double lastvalue ;
    
    std::vector<Point2<double>> input_;
};

#endif  // OBJECTIVEFUNCTION_H
