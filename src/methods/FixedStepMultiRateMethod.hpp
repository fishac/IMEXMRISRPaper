#ifndef FIXEDSTEPMULTIRATEMETHOD_DEFINED__
#define FIXEDSTEPMULTIRATEMETHOD_DEFINED__

#include "FixedStepMultiRateStep.hpp"

using namespace arma;

class FixedStepMultiRateMethod {
public:
	const char* name;
	virtual mat solve(double t_0, double H_, int M_, vec* y_0, vec* output_tspan, FixedStepMultiRateStep* step) {
		mat x(1,1);
		return x;
	}
};

#endif