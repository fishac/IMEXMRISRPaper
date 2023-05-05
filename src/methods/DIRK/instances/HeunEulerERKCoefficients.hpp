#ifndef HEUNEULERERKCOEFFICIENTS_DEFINED__
#define HEUNEULERERKCOEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class HeunEulerERKCoefficients: public SingleRateMethodCoefficients {
public:
	HeunEulerERKCoefficients() {
		name = "HeunEulerERK";
		
		num_stages = 2;
		primary_order = 2.0;
		secondary_order = 1.0;

		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = 1.0;

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