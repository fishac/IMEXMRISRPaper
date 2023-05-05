#ifndef ESDIRK325COEFFICIENTS_DEFINED__
#define ESDIRK325COEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"
#include <cmath>

using namespace arma;

class ESDIRK325Coefficients: public SingleRateMethodCoefficients {
public:
	ESDIRK325Coefficients() {
		name = "ESDIRK3(2)5L[2]SA";
		
		num_stages = 5;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		double sqrt2 = std::sqrt(2.0);

		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = 9.0/40.0;
		A(1,1) = 9.0/40.0;

		A(2,0) = 9.0*(1.0+sqrt2)/80.0;
		A(2,1) = 9.0*(1.0+sqrt2)/80.0;
		A(2,2) = 9.0/40.0;

		A(3,0) = (22.0+15.0*sqrt2)/(80.0*(1.0+sqrt2));
		A(3,1) = (22.0+15.0*sqrt2)/(80.0*(1.0+sqrt2));
		A(3,2) = -7.0/(40.0*(1.0+sqrt2));
		A(3,3) = 9.0/40.0;
		
		A(4,0) = (2398.0+1205.0*sqrt2)/(2835.0*(4.0+3.0*sqrt2));
		A(4,1) = (2398.0+1205.0*sqrt2)/(2835.0*(4.0+3.0*sqrt2));
		A(4,2) = -2374.0*(1.0+2.0*sqrt2)/(2835.0*(5.0+3.0*sqrt2));
		A(4,3) = 5827.0/7560.0;
		A(4,4) = 9.0/40.0;

		b = vec(num_stages,fill::zeros);
		b(0) = (2398.0+1205.0*sqrt2)/(2835.0*(4.0+3.0*sqrt2));
		b(1) = (2398.0+1205.0*sqrt2)/(2835.0*(4.0+3.0*sqrt2));
		b(2) = -2374.0*(1.0+2.0*sqrt2)/(2835.0*(5.0+3.0*sqrt2));
		b(3) = 5827.0/7560.0;
		b(4) = 9.0/40.0;

		d = vec(num_stages,fill::zeros);
      	d(0) = 4555948517383.0/24713416420891.0;
		d(1) = 4555948517383.0/24713416420891.0;
		d(2) = -7107561914881.0/25547637784726.0;
		d(3) = 30698249.0/44052120.0;
		d(4) = 49563.0/233080.0;
		
		c = A*vec(num_stages,fill::ones);
	}
};

#endif