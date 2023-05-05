#ifndef VERNERERKCOEFFICIENTS_DEFINED__
#define VERNERERKCOEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class VernerERKCoefficients: public SingleRateMethodCoefficients {
public:
	VernerERKCoefficients() {
		name = "VernerERK";
		
		num_stages = 8;
		primary_order = 6.0;
		secondary_order = 5.0;

		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = 1.0/6.0;
		A(2,0) = 4.0/75.0;
		A(2,1) = 16.0/75.0;
		A(3,0) = 5.0/6.0;
		A(3,1) = -8.0/3.0;
		A(3,2) = 5.0/2.0;
		A(4,0) = -165.0/64.0;
		A(4,1) = 55.0/6.0;
		A(4,2) = -425.0/64.0;
		A(4,3) = 85.0/96.0;
		A(5,0) = 12.0/5.0;
		A(5,1) = -8.0;
		A(5,2) = 4015.0/612.0;
		A(5,3) = -11.0/36.0;
		A(5,4) = 88.0/255.0;
		A(6,0) = -8263.0/15000.0;
		A(6,1) = 124.0/75.0;
		A(6,2) = -643.0/680.0;
		A(6,3) = -81.0/250.0;
		A(6,4) = 2484.0/10625.0;
		A(7,0) = 3501.0/1720.0;
		A(7,1) = -300.0/43.0;
		A(7,2) = 297275.0/52632.0;
		A(7,3) = -319.0/2322.0;
		A(7,4) = 24068.0/84065.0;
		A(7,6) = 3850.0/26703.0;

		b = vec(num_stages,fill::zeros);
		b(0) = 3.0/40.0;
		b(2) = 875.0/2244.0;
		b(3) = 23.0/72.0;
		b(4) = 264.0/1955.0;
		b(6) = 125.0/11592.0;
		b(7) = 43.0/616.0;

		d = vec(num_stages,fill::zeros);
		d(0) = 13.0/160.0;
		d(2) = 2375.0/5984.0;
		d(3) = 5.0/16.0;
		d(4) = 12.0/85.0;
		d(5) = 3.0/44.0;

		c = A*vec(num_stages,fill::ones);
	}
};

#endif