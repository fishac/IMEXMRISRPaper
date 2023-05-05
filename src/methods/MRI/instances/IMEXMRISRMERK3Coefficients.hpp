#ifndef IMEXMRISRMERK3COEFFICIENTS_DEFINED__
#define IMEXMRISRMERK3COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISRMERK3Coefficients : public MRICoefficients {
public:
	IMEXMRISRMERK3Coefficients() {
		name = "IMEXMRISRMERK3";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 2;
		num_stages = 4;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 0.5;
		O1(2,0) = 2.0/3.0;
		O1(3,0) = 1.0;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(2,0) = -4.0/(9.0*0.5);
		O2(2,1) = 4.0/(9.0*0.5);
		O2(3,0) = -3.0/2.0;
		O2(3,1) = 0.0;
		O2(3,2) = 3.0/2.0;
		omegas.push_back(O2);
		

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.5;
		c(2) = 2.0/3.0;
		c(3) = 1.0;
		
		num_groups = 3;
		stage_groups = {{1}, {2}, {3}};
	}
};

#endif