#ifndef IMEXMRISRTEST4COEFFICIENTS_DEFINED__
#define IMEXMRISRTEST4COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISRTEST4Coefficients : public MRICoefficients {
public:
	IMEXMRISRTEST4Coefficients() {
		name = "IMEXMRISRTEST4";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 2;
		num_stages = 7;
		primary_order = 4.0;
		secondary_order = 3.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = -0.25;
		G1(1,1) = 0.25;
		G1(2,0) = 0.25;
		G1(2,1) = -0.5;
		G1(2,2) = 0.25;
		G1(3,0) = 13.0/100.0;
		G1(3,1) = -7.0/30.0;
		G1(3,2) = -11.0/75.0;
		G1(3,3) = 0.25;
		G1(4,0) = 6.0/85.0;
		G1(4,1) = -301.0/1360.0;
		G1(4,2) = -99.0/544.0;
		G1(4,3) = 45.0/544.0;
		G1(4,4) = 0.25;
		G1(5,0) = 0.0;
		G1(5,1) = -9.0/4.0;
		G1(5,2) = -19.0/48.0;
		G1(5,3) = -75.0/16.0;
		G1(5,4) = 85.0/12.0;
		G1(5,5) = 0.25;
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 0.25;
		O1(2,0) = 0.75;
		O1(3,0) = 11.0/20.0;
		O1(4,0) = 1.0/2.0;
		O1(5,0) = 1.0;
		O1(6,2) = -269.0/24.0;
		O1(6,3) = 1175.0/72.0;
		O1(6,4) = -85.0/12.0;
		O1(6,5) = 107.0/36.0;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(2,0) = -2.0;
		O2(2,1) = 2.0;
		O2(3,0) = -34.0/25.0;
		O2(3,1) = 86.0/75.0;
		O2(3,2) = 16.0/75.0;
		O2(4,0) = -97.0/85.0;
		O2(4,1) = 84.0/85.0;
		O2(4,2) = 179.0/680.0;
		O2(4,3) = -15.0/136.0;
		O2(5,0) = -2.0;
		O2(5,1) = 79.0/12.0;
		O2(5,2) = -5.0/4.0;
		O2(5,3) = 25.0;
		O2(5,4) = -85.0/3.0;
		O2(6,1) = 25.0/12.0;
		O2(6,2) = 163.0/8.0;
		O2(6,3) = -1225.0/72.0;
		O2(6,4) = 0.0;
		O2(6,5) = -49.0/9.0;
		omegas.push_back(O2);
		

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.25;
		c(2) = 0.75;
		c(3) = 11.0/20.0;
		c(4) = 1.0/2.0;
		c(5) = 1.0;
		c(6) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1}, {2}, {3}, {4}, {5}, {6}};
	}
};

#endif