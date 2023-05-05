#ifndef MRIGARKIRK21ACOEFFICIENTS_DEFINED__
#define MRIGARKIRK21ACOEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MRIGARKIRK21aCoefficients: public MRICoefficients {
public:
	MRIGARKIRK21aCoefficients() {
		name = "MRIGARKIRK21a";
		explicit_mrigark = false;
		method_type = 0;
		
		num_gammas = 1;
		num_omegas = 0;
		num_stages = 4;
		primary_order = 2.0;
		secondary_order = 1.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = 1.0;
		G1(2,0) = -0.5;
		G1(2,2) = 0.5;
		G1(4,0) = -0.5;
		G1(4,2) = 0.5;
		gammas.push_back(G1);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 1.0;
		c(2) = 1.0;
		c(3) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2, 3}};
	}
};

#endif