#ifndef MRIGARKERK22ACOEFFICIENTS_DEFINED__
#define MRIGARKERK22ACOEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MRIGARKERK22aCoefficients: public MRICoefficients {
public:
	MRIGARKERK22aCoefficients() {
		name = "MRIGARKERK22a";
		explicit_mrigark = true;
		method_type = 0;
		
		num_gammas = 1;
		num_omegas = 0;
		num_stages = 3;
		primary_order = 2.0;
		secondary_order = 1.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = 0.5;
		G1(2,0) = -0.5;
		G1(2,1) = 1.0;
		G1(3,0) = 0.5;
		gammas.push_back(G1);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.5;
		c(2) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2}};
	}
};

#endif