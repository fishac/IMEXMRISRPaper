#ifndef MERK5COEFFICIENTS_DEFINED__
#define MERK5COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MERK5Coefficients: public MRICoefficients {
public:
	MERK5Coefficients() {
		name = "MERK5";
		explicit_mrigark = true;
		method_type = 1;
		
		num_gammas = 3;
		num_omegas = 0;
		num_stages = 11;
		primary_order = 5.0;
		secondary_order = 4.0;
		
		mat G1(num_stages+1,num_stages-2,fill::zeros);
		G1(2,0) = 2.0;
		G1(3,0) = 2.0;
		G1(4,1) = -4.0;
		G1(4,2) = 9.0;
		G1(5,1) = -4.0;
		G1(5,2) = 9.0;
		G1(6,1) = -4.0;
		G1(6,2) = 9.0;
		G1(7,3) = 4.0;
		G1(7,4) = -27.0;
		G1(7,5) = 32.0;
		G1(8,3) = 4.0;
		G1(8,4) = -27.0;
		G1(8,5) = 32.0;
		G1(9,3) = 4.0;
		G1(9,4) = -27.0;
		G1(9,5) = 32.0;
		G1(10,6) = 500.0/7.0;
		G1(10,7) = 28.0;
		G1(10,8) = -189.0/2.0;

		mat G2(num_stages+1,num_stages-2,fill::zeros);
		G2(4,1) = 12.0;
		G2(4,2) = -18.0;
		G2(5,1) = 12.0;
		G2(5,2) = -18.0;
		G2(6,1) = 12.0;
		G2(6,2) = -18.0;
		G2(7,3) = -28.0;
		G2(7,4) = 162.0;
		G2(7,5) = -160.0;
		G2(8,3) = -28.0;
		G2(8,4) = 162.0;
		G2(8,5) = -160.0;
		G2(9,3) = -28.0;
		G2(9,4) = 162.0;
		G2(9,5) = -160.0;
		G2(10,6) = -250.0;
		G2(10,7) = -82.0;
		G2(10,8) = 324.0;

		mat G3(num_stages+1,num_stages-2,fill::zeros);
		G3(7,3) = 48.0;
		G3(7,4) = -216.0;
		G3(7,5) = 192.0;
		G3(8,3) = 48.0;
		G3(8,4) = -216.0;
		G3(8,5) = 192.0;
		G3(9,3) = 48.0;
		G3(9,4) = -216.0;
		G3(9,5) = 192.0;
		G3(10,6) = 1500.0/7.0;
		G3(10,7) = 60.0;
		G3(10,8) = -270.0;

		gammas.push_back(G1);
		gammas.push_back(G2);
		gammas.push_back(G3);
		
		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.5;
		c(2) = 0.5;
		c(3) = 1.0/3.0;
		c(4) = 0.5;
		c(5) = 1.0/3.0;
		c(6) = 1.0/4.0;
		c(7) = 7.0/10.0;
		c(8) = 0.5;
		c(9) = 2.0/3.0;
		c(10) = 1.0;
		
		num_groups = 5;
		stage_groups = {{1}, {3,2}, {6,5,4}, {8,9,7}, {10}};
		//num_groups = 10;
		//stage_groups = {{1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}};
	}
};

#endif