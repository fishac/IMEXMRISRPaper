#ifndef IMEXMRISRMERK4COEFFICIENTS_DEFINED__
#define IMEXMRISRMERK4COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISRMERK4Coefficients : public MRICoefficients {
public:
	IMEXMRISRMERK4Coefficients() {
		name = "IMEXMRISRMERK4";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 3;
		num_stages = 7;
		primary_order = 4.0;
		secondary_order = 3.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 0.5;
		O1(2,0) = 0.5;
		O1(3,0) = 1.0/3.0;
		O1(4,0) = 5.0/6.0;
		O1(5,0) = 1.0/3.0;
		O1(6,0) = 1.0;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(2,0) = -0.5;
		O2(2,1) = 0.5;
		O2(3,0) = -2.0/9.0;
		O2(3,1) = 2.0/9.0;
		O2(4,0) = -125.0/36.0;
		O2(4,2) = -25.0/9.0;
		O2(4,3) = 25.0/4.0;
		O2(5,0) = -5.0/9.0;
		O2(5,2) = -4.0/9.0;
		O2(5,3) = 1.0;
		O2(6,0) = -21.0/5.0;
		O2(6,4) = -4.0/5.0;
		O2(6,5) = 5.0;
		omegas.push_back(O2);
		
		mat O3(num_stages+1,num_stages,fill::zeros);
		O3(4,0) = 125.0/36.0;
		O3(4,2) = 125.0/18.0;
		O3(4,3) = -125.0/12.0;
		O3(5,0) = 2.0/9.0;
		O3(5,2) = 4.0/9.0;
		O3(5,3) = -2.0/3.0;
		O3(6,0) = 18.0/5.0;
		O3(6,4) = 12.0/5.0;
		O3(6,5) = -6.0;
		omegas.push_back(O3);
		

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.5;
		c(2) = 0.5;
		c(3) = 1.0/3.0;
		c(4) = 5.0/6.0;
		c(5) = 1.0/3.0;
		c(6) = 1.0;
		
		num_groups = 4;
		stage_groups = {{1}, {3,2}, {5,4}, {6}};
		//num_groups = 6;
		//stage_groups = {{1}, {2}, {3}, {4}, {5}, {6}};
	}
};

#endif