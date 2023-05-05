#ifndef ZONNEVELDERKCOEFFICIENTS_DEFINED__
#define ZONNEVELDERKCOEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class ZonneveldERKCoefficients: public SingleRateMethodCoefficients {
public:
	ZonneveldERKCoefficients() {
		name = "ZonneveldERK";
		
		num_stages = 5;
		primary_order = 4.0;
		secondary_order = 3.0;
		
		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = 1.0/2.0;

		A(2,0) = 0.0;
		A(2,1) = 1.0/2.0;

		A(3,0) = 0.0;
		A(3,1) = 0.0;
		A(3,2) = 1.0;

		A(4,0) = 5.0/32.0;
		A(4,1) = 7.0/32.0;
		A(4,2) = 13.0/32.0;
		A(4,3) = -1.0/32.0;

		b = vec(num_stages,fill::zeros);
		b(0) = 1.0/6.0;
		b(1) = 1.0/3.0;
		b(2) = 1.0/3.0;
		b(3) = 1.0/6.0;
		b(4) = 0.0;

		d = vec(num_stages,fill::zeros);
		d(0) = -1.0/2.0;
		d(1) = 7.0/3.0;
		d(2) = 7.0/3.0;
		d(3) = 13.0/6.0;
		d(4) = -16.0/3.0;

		c = A*vec(num_stages,fill::ones);
	}
};

#endif