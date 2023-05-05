#ifndef BOGACKISHAMPINEERKCOEFFICIENTS_DEFINED__
#define BOGACKISHAMPINEERKCOEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class BogackiShampineERKCoefficients: public SingleRateMethodCoefficients {
public:
	BogackiShampineERKCoefficients() {
		name = "BogackiShampineERK";
		
		num_stages = 4;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = 0.5;

		A(2,1) = 0.75;

		A(3,0) = 2.0/9.0;
		A(3,1) = 1.0/3.0;
		A(3,2) = 4.0/9.0;

		b = vec(num_stages,fill::zeros);
		b(0) = 2.0/9.0;
		b(1) = 1.0/3.0;
		b(2) = 4.0/9.0;
		b(3) = 0.0;

		d = vec(num_stages,fill::zeros);
		d(0) = 7.0/24.0;
		d(1) = 1.0/4.0;
		d(2) = 1.0/3.0;
		d(3) = 1.0/8.0;

		c = A*vec(num_stages,fill::ones);
	}
};

#endif