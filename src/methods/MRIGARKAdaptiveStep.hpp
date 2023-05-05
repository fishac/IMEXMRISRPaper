#ifndef MRIGARKADAPTIVESTEP_DEFINED__
#define MRIGARKADAPTIVESTEP_DEFINED__

using namespace arma;

struct MRIGARKAdaptiveStepReturnValue {
	vec y;
	double ess;
	double esf;
	int total_microtimesteps = 0;
	int total_successful_microtimesteps = 0;
	int status = 0;
};

class MRIGARKAdaptiveStep {
public:
	const char* name;
	virtual void step_solution(double t, double H, int M, vec* y_0, const char* measurement_type, MRIGARKAdaptiveStepReturnValue* ret) {}
	virtual void step_solution(double t, double H, double tolf, vec* y_0, const char* measurement_type, MRIGARKAdaptiveStepReturnValue* ret) {}
};

#endif