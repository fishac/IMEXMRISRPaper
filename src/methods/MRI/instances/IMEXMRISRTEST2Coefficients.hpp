#ifndef IMEXMRISRTEST2COEFFICIENTS_DEFINED__
#define IMEXMRISRTEST2COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISRTEST2Coefficients : public MRICoefficients {
public:
	IMEXMRISRTEST2Coefficients() {
		name = "IMEXMRISRTEST2";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 2;
		num_stages = 4;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = -1.0/6.0;
		G1(1,1) = 1.0/2.0;
		G1(2,0) = 1.0/12.0;
		G1(2,1) = 1.0/2.0;
		G1(2,2) = 1.0/12.0;
		G1(3,0) = 1.0/4.0;
		G1(3,1) = 0.0;
		G1(3,2) = 3.0/4.0;
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 1.0/3.0;
		O1(2,0) = 1.0/6.0;
		O1(2,1) = 1.0/2.0;
		O1(3,0) = 3.0/2.0;
		O1(3,1) = -1.0;
		O1(3,2) = 1.0/2.0;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(2,0) = -1.0/3.0;
		O2(2,1) = 1.0/3.0;
		O2(3,0) = -5.0/2.0;
		O2(3,1) = 2.0;
		O2(3,2) = 1.0/2.0;
		omegas.push_back(O2);
		

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 1.0/3.0;
		c(2) = 2.0/3.0;
		c(3) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1}, {2}, {3}};
	}
};

#endif