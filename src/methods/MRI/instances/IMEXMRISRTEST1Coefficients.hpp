#ifndef IMEXMRISRTEST1COEFFICIENTS_DEFINED__
#define IMEXMRISRTEST1COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISRTEST1Coefficients : public MRICoefficients {
public:
	IMEXMRISRTEST1Coefficients() {
		name = "IMEXMRISRTEST1";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 1;
		num_stages = 3;
		primary_order = 2.0;
		secondary_order = 1.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = 1.0/4.0;
		G1(1,1) = 1.0/4.0;
		G1(2,0) = 1.0/4.0;
		G1(2,1) = 1.0/2.0;
		G1(2,2) = 1.0/4.0;
		G1(3,0) = 0.0;
		G1(3,1) = 0.0;
		G1(3,2) = 0.0;

		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 1.0/2.0;
		O1(2,1) = 1.0;
		
		omegas.push_back(O1);
		

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 1.0/2.0;
		c(2) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1}, {2}};
	}
};

#endif