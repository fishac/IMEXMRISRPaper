#ifndef ADAPTIVESTEPMULTIRATEMETHOD_DEFINED__
#define ADAPTIVESTEPMULTIRATEMETHOD_DEFINED__

#include "MRIGARKAdaptiveStep.hpp"
#include "Controller.hpp"

using namespace arma;

struct AdaptiveMultiRateMethodReturnValue {
	mat Y;
	std::vector<double> ts;
	std::vector<double> Hs;
	std::vector<int> Ms;
	std::vector<double> tolfs;
	std::vector<int> full_function_evals;
	std::vector<int> fast_function_evals;
	std::vector<int> slow_function_evals;
	std::vector<int> implicit_function_evals;
	std::vector<int> explicit_function_evals;
	std::vector<int> full_jacobian_evals;
	std::vector<int> fast_jacobian_evals;
	std::vector<int> slow_jacobian_evals;
	std::vector<int> implicit_jacobian_evals;
	int total_timesteps = 0;
	int total_successful_timesteps = 0;
	int total_microtimesteps = 0;
	int total_successful_microtimesteps = 0;
	int status = 0;
};

class AdaptiveStepMultiRateMethod {
public:
	const char* name;
	virtual void solve(double t_0, double H_0, int M_0, vec* y_0, vec* output_tspan, MRIGARKAdaptiveStep* mrigark_step, Controller* controller, AdaptiveMultiRateMethodReturnValue* ret) {}
};

#endif