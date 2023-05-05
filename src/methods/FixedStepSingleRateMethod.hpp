#ifndef FIXEDSTEPSINGLERATEMETHOD_DEFINED__
#define FIXEDSTEPSINGLERATEMETHOD_DEFINED__

using namespace arma;

class FixedStepSingleRateMethod {
public:
	const char* name;
	virtual vec solve(double t_0, double t_f, double h_, vec* y_0) {
		vec x(1);
		return x;
	}
};

#endif