#ifndef SINGLERATEMETHODCOEFFICIENTS_DEFINED__
#define SINGLERATEMETHODCOEFFICIENTS_DEFINED__

using namespace arma;

class SingleRateMethodCoefficients {
public:
	const char* name;
	mat A;
	vec b;
	vec d;
	vec c;
	int num_stages;
	double primary_order;
	double secondary_order;
};

#endif