#ifndef ESDIRK436COEFFICIENTS_DEFINED__
#define ESDIRK436COEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"
#include <cmath>

using namespace arma;

class ESDIRK436Coefficients: public SingleRateMethodCoefficients {
public:
	ESDIRK436Coefficients() {
		name = "ESDIRK4(3)6L[2]SA";
		
		num_stages = 6;
		primary_order = 4.0;
		secondary_order = 3.0;
		
		A = mat(num_stages,num_stages,fill::zeros);
		A(1,0) = 31.0/125.0;
		A(1,1) = 31.0/125.0;

		A(2,0) = -360286518617.0/7014585480527.0;
		A(2,1) = -360286518617.0/7014585480527.0;
		A(2,2) = 31.0/125.0;

		A(3,0) = -506388693497.0/5937754990171.0;
		A(3,1) = -506388693497.0/5937754990171.0;
		A(3,2) = 7149918333491.0/13390931526268.0;
		A(3,3) = 31.0/125.0;
		
		A(4,0) = -7628305438933.0/11061539393788.0;
		A(4,1) = -7628305438933.0/11061539393788.0;
		A(4,2) = 21592626537567.0/14352247503901.0;
		A(4,3) = 11630056083252.0/17263101053231.0;
		A(4,4) = 31.0/125.0;
		
		A(5,0) = -12917657251.0/5222094901039.0;
		A(5,1) = -12917657251.0/5222094901039.0;
		A(5,2) = 5602338284630.0/15643096342197.0;
		A(5,3) = 9002339615474.0/18125249312447.0;
		A(5,4) = -2420307481369.0/24731958684496.0;
		A(5,5) = 31.0/125.0;

		b = vec(num_stages,fill::zeros);
		b(0) = -12917657251.0/5222094901039.0;
		b(1) = -12917657251.0/5222094901039.0;
		b(2) = 5602338284630.0/15643096342197.0;
		b(3) = 9002339615474.0/18125249312447.0;
		b(4) = -2420307481369.0/24731958684496.0;
		b(5) = 31.0/125.0;

		d = vec(num_stages,fill::zeros);
      	d(0) = -1007911106287.0/12117826057527.0;
		d(1) = -1007911106287.0/12117826057527.0;
		d(2) = 17694008993113.0/35931961998873.0;
		d(3) = 5816803040497.0/11256217655929.0;
		d(4) = -538664890905.0/7490061179786.0;
		d(5) = 2032560730450.0/8872919773257.0;
		
		c = A*vec(num_stages,fill::ones);
	}
};

#endif