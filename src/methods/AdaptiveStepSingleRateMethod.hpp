#ifndef ADAPTIVESTEPSINGLERATEMETHOD_DEFINED__
#define ADAPTIVESTEPSINGLERATEMETHOD_DEFINED__

using namespace arma;

struct AdaptiveSingleRateMethodReturnValue {
	mat Y;
	std::vector<double> ts;
	std::vector<double> hs;
	int total_timesteps = 0;
	int total_successful_timesteps = 0;
	std::vector<int> full_function_evals;
	std::vector<int> fast_function_evals;
	std::vector<int> slow_function_evals;
	std::vector<int> implicit_function_evals;
	std::vector<int> explicit_function_evals;
	std::vector<int> full_jacobian_evals;
	std::vector<int> fast_jacobian_evals;
	std::vector<int> slow_jacobian_evals;
	std::vector<int> implicit_jacobian_evals;
	int status = 0;
};

class AdaptiveStepSingleRateMethod {
public:
	const char* name;
	virtual void solve(double t_0, double t_f, double h_, vec* y_0, vec* output_tspan, AdaptiveSingleRateMethodReturnValue* ret) {}
};

#endif