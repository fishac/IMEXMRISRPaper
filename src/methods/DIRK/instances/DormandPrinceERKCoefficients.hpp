#ifndef DormandPrinceERKCOEFFICIENTS_DEFINED__
#define DormandPrinceERKCOEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class DormandPrinceERKCoefficients: public SingleRateMethodCoefficients {
public:
	DormandPrinceERKCoefficients() {
		name = "DormandPrinceERK";
		
		num_stages = 7;
		primary_order = 5.0;
		secondary_order = 4.0;
		
		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = 1.0/5.0;

		A(2,0) = 3.0/40.0;
		A(2,1) = 9.0/40.0;

		A(3,0) = 44.0/45.0;
		A(3,1) = -56.0/15.0;
		A(3,2) = 32.0/9.0;

		A(4,0) = 19372.0/6561.0;
		A(4,1) = -25360.0/2187.0;
		A(4,2) = 64448.0/6561.0;
		A(4,3) = -212.0/729.0;

		A(5,0) = 9017.0/3168.0;
		A(5,1) = -355.0/33.0;
		A(5,2) = 46732.0/5247.0;
		A(5,3) = 49.0/176.0;
		A(5,4) = -5103.0/18656.0;

		A(6,0) = 35.0/384.0;
		A(6,1) = 0.0;
		A(6,2) = 500.0/1113.0;
		A(6,3) = 125.0/192.0;
		A(6,4) = -2187.0/6784.0;
		A(6,5) = 11.0/84.0;

		b = vec(num_stages,fill::zeros);
		b(0) = 35.0/384.0;
		b(1) = 0.0;
		b(2) = 500.0/1113.0;
		b(3) = 125.0/192.0;
		b(4) = -2187.0/6784.0;
		b(5) = 11.0/84.0;
		b(6) = 0.0;

		d = vec(num_stages,fill::zeros);
		d(0) = 5179.0/57600.0;
		d(1) = 0.0;
		d(2) = 7571.0/16695.0;
		d(3) = 393.0/640.0;
		d(4) = -92097.0/339200.0;
		d(5) = 187.0/2100.0;
		d(6) = 1.0/40.0;

		c = A*vec(num_stages,fill::ones);
	}
};

#endif