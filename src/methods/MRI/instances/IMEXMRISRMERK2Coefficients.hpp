#ifndef IMEXMRISRMERK2COEFFICIENTS_DEFINED__
#define IMEXMRISRMERK2COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISRMERK2Coefficients : public MRICoefficients {
public:
	IMEXMRISRMERK2Coefficients() {
		name = "IMEXMRISRMERK2";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 2;
		num_stages = 3;
		primary_order = 2.0;
		secondary_order = 1.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 0.5;
		O1(2,0) = 1.0;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(2,0) = -2.0;
		O2(2,1) = 2.0;
		omegas.push_back(O2);
		

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.5;
		c(2) = 1.0;
		
		num_groups = 2;
		stage_groups = {{1}, {2}};
	}
};

#endif