#ifndef FIXEDSTEPMULTIRATESTEP_DEFINED__
#define FIXEDSTEPMULTIRATESTEP_DEFINED__

using namespace arma;

class FixedStepMultiRateStep {
public:
	const char* name;
	virtual void step_solution(double t, double H, int M, vec* y_prev, vec* solution_vec) {};
	virtual void set_coeffs() {};
	virtual void refresh_coupling_coeffs(int M) {};
};

#endif