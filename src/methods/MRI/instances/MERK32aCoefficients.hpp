#ifndef MERK32aCOEFFICIENTS_DEFINED__
#define MERK32aCOEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MERK32aCoefficients: public MRICoefficients {
public:
	MERK32aCoefficients() {
		name = "MERK32a";
		explicit_mrigark = true;
		method_type = 1;
		
		num_gammas = 2;
		num_omegas = 0;
		num_stages = 5;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		c = vec(num_stages+1,fill::zeros);
		c(0) = 0.0;
		c(1) = 1.0/2.0;
		c(2) = 1.0/2.0;
		c(3) = 1.0/3.0;
		c(4) = 1.0;
		c(5) = 1.0;
		
		mat G1(num_stages+1,num_stages-2,fill::zeros);
		G1(2,0) = 1.0/c(1);
		G1(3,0) = 1.0/c(1);
		G1(4,1) = -c(3)/(c(2)*(c(2)-c(3)));
		G1(4,2) = c(2)/(c(3)*(c(2)-c(3)));
		//
		G1(5,1) = -c(3)/(c(2)*(c(2)-c(3)));
		G1(5,2) = c(2)/(c(3)*(c(2)-c(3)));

		mat G2(num_stages+1,num_stages-2,fill::zeros);
		G2(4,1) = 1.0/(c(2)*(c(2)-c(3)));
		G2(4,2) = -1.0/(c(3)*(c(2)-c(3)));
		//G2(5,1) = 1.0/(c(2)*(c(2)-c(3)));
		//G2(5,2) = -1.0/(c(3)*(c(2)-c(3)));

		gammas.push_back(G1);
		gammas.push_back(G2);
		
		num_groups = 3;
		stage_groups = {{1}, {3,2}, {4}};
	}
};

#endif