#ifndef IMEXMRIGARK21COEFFICIENTS_DEFINED__
#define IMEXMRIGARK21COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRIGARK21Coefficients: public MRICoefficients {
public:
	IMEXMRIGARK21Coefficients() {
		name = "IMEXMRIGARK21";
		explicit_mrigark = false;
		method_type = 0;
		
		num_gammas = 2;
		num_omegas = 2;
		num_stages = 4;
		primary_order = 2.0;
		secondary_order = 1.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = 1.0;
		G1(2,0) = -567.0/290.0;
		G1(2,2) = 567.0/290.0;
		G1(3,0) = -119623.0/80910.0;
		G1(3,2) = -133.0/279.0;
		G1(3,3) = 567.0/290.0;
		G1(4,0) = 1678.0/735.0;
		gammas.push_back(G1);
		
		mat G2(num_stages+1,num_stages,fill::zeros);
		G2(2,0) = -567.0/290.0;
		G2(2,2) = 567.0/290.0;
		G2(3,0) = 126583.0/16182.0;
		G2(3,2) = -395554.0/40455.0;
		G2(3,3) = 567.0/290.0;
		G2(4,0) = 121.0/106.0;
		G2(4,2) = -444671.0/77910.0;
		gammas.push_back(G2);

		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 1.0;
		O1(3,0) = 80.0/311.0;
		O1(3,2) = -80.0/311.0;
		O1(4,0) = -17.0/269.0;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(3,0) = -471.0/311.0;
		O2(3,2) = 471.0/311.0;
		O2(4,0) = -11.0/348.0;
		O2(4,2) = 14791.0/93612.0;
		omegas.push_back(O2);

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