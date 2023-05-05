#ifndef SDIRK21COEFFICIENTS_DEFINED__
#define SDIRK21COEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class SDIRK21Coefficients: public SingleRateMethodCoefficients {
public:
	SDIRK21Coefficients() {
		name = "SDIRK-2-1";
		
		num_stages = 2;
		primary_order = 2.0;
		secondary_order = 1.0;
		
		A = mat(num_stages,num_stages,fill::zeros);
		A(0,0) = 1.0;
		A(1,0) = -1.0;
		A(1,1) = 1.0;

		b = vec(num_stages,fill::zeros);
		b(0) = 0.5;
		b(1) = 0.5;

		d = vec(num_stages,fill::zeros);
		d(0) = 1.0;
		d(1) = 0.0;

		c = A*vec(num_stages,fill::ones);
	}
};

#endif