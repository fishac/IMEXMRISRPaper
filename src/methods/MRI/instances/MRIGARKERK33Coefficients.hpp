#ifndef MRIGARKERK33COEFFICIENTS_DEFINED__
#define MRIGARKERK33COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MRIGARKERK33Coefficients: public MRICoefficients {
public:
	MRIGARKERK33Coefficients() {
		name = "MRIGARKERK33";
		explicit_mrigark = true;
		method_type = 0;
		
		num_gammas = 2;
		num_omegas = 0;
		num_stages = 4;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = 1.0/3.0;
		G1(2,0) = -1.0/3.0;
		G1(2,1) = 2.0/3.0;
		G1(3,1) = -2.0/3.0;
		G1(3,2) = 1.0;
		G1(4,0) = 1.0/12.0;
		G1(4,1) = -1.0/3.0;
		G1(4,2) = 7.0/12.0;

		mat G2(num_stages+1,num_stages,fill::zeros);
		G2(3,0) = 0.5;
		G2(3,2) = -0.5;

		gammas.push_back(G1);
		gammas.push_back(G2);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 1.0/3.0;
		c(2) = 2.0/3.0;
		c(3) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2, 3}};
	}
};

#endif