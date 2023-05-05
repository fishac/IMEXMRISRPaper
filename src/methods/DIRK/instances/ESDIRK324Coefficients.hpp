#ifndef ESDIRK324COEFFICIENTS_DEFINED__
#define ESDIRK324COEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class ESDIRK324Coefficients: public SingleRateMethodCoefficients {
public:
	ESDIRK324Coefficients() {
		name = "ESDIRK3(2)4L[2]SA";
		
		num_stages = 4;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		double gamma = 0.43586652150845899941601945;
		double gamma2 = gamma*gamma;
		double gamma3 = gamma*gamma*gamma;
		double gamma4 = gamma*gamma*gamma*gamma;
		double gamma5 = gamma*gamma*gamma*gamma*gamma;
		
		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 2.0*gamma;
		c(2) = 3.0/5.0;
		c(3) = 1.0;
		
		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = gamma;
		A(1,1) = gamma;

		A(2,2) = gamma;
		A(2,1) = c(2)*(c(2)-2.0*gamma)/(4.0*gamma);
		A(2,0) = c(2)-A(2,1)-gamma;

		A(3,3) = gamma;
		A(3,2) = (1.0-6.0*gamma+6.0*gamma2)/(3.0*c(2)*(c(2)-2.0*gamma));
		A(3,1) = (-2.0+3.0*c(2)+6.0*gamma*(1.0-c(2)))/(12.0*gamma*(c(2)-2.0*gamma));
		A(3,0) = 1.0-A(3,1)-A(3,2)-gamma;

		b = vec(num_stages,fill::zeros);
		b(3) = gamma;
		b(2) = (1.0-6.0*gamma+6.0*gamma2)/(3.0*c(2)*(c(2)-2.0*gamma));
		b(1) = (-2.0+3.0*c(2)+6.0*gamma*(1.0-c(2)))/(12.0*gamma*(c(2)-2.0*gamma));
		b(0) = 1.0-A(3,1)-A(3,2)-gamma;

		d = vec(num_stages,fill::zeros);
      	d(3) = -3.0*gamma2*(-1.0+4.0*gamma-2.0*gamma2+gamma3)/(1.0-6.0*gamma+6.0*gamma2);
      	d(2) = -gamma*(-2.0+21.0*gamma-68.0*gamma2+79.0*gamma3-33.0*gamma4+12.0*gamma5)/(c(2)*(c(2)-2.0*gamma)*(1.0-6.0*gamma+6.0*gamma2));
		d(1) = c(2)*(-1.0 + 6.0*gamma - 24.0*gamma3 + 12.0*gamma4 - 6.0*gamma5)/(4.0*gamma*(2.0*gamma-c(2))*(1.0-6.0*gamma+6.0*gamma2)) + (3.0-27.0*gamma+68.0*gamma2-55.0*gamma3+21.0*gamma4-6.0*gamma5)/(2.0*(2.0*gamma-c(2))*(1.0-6.0*gamma+6.0*gamma2));                                                       
		d(0) = 1.0-b(1)-b(2)-b(3);
	}
};

#endif