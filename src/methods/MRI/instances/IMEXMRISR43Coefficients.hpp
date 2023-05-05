#ifndef IMEXMRISR43COEFFICIENTS_DEFINED__
#define IMEXMRISR43COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISR43Coefficients : public MRICoefficients {
public:
	IMEXMRISR43Coefficients() {
		name = "IMEXMRISR43";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 2;
		num_stages = 7;
		primary_order = 4.0;
		secondary_order = 3.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = -1.0/4.0;
		G1(1,1) = 1.0/4.0;
		G1(2,0) = 1.0/4.0;
		G1(2,1) = -1.0/2.0;
		G1(2,2) = 1.0/4.0;
		G1(3,0) = 13.0/100.0;
		G1(3,1) = -7.0/30.0;
		G1(3,2) = -11.0/75.0;
		G1(3,3) = 1.0/4.0;
		G1(4,0) = 6.0/85.0;
		G1(4,1) = -301.0/1360.0;
		G1(4,2) = -99.0/544.0;
		G1(4,3) = 45.0/544.0;
		G1(4,4) = 1.0/4.0;
		G1(5,0) = 0.0;
		G1(5,1) = -9.0/4.0;
		G1(5,2) = -19.0/48.0;
		G1(5,3) = -75.0/16.0;
		G1(5,4) = 85.0/12.0;
		G1(5,5) = 1.0/4.0;
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 1.0/4.0;
		O1(2,0) = 9.0/8.0;
		O1(2,1) = -3.0/8.0;
		O1(3,0) = 187.0/2340.0;
		O1(3,1) = 7.0/9.0;
		O1(3,2) = -4.0/13.0;
		O1(4,0) = 64.0/165.0;
		O1(4,1) = 1.0/6.0;
		O1(4,2) = -3.0/5.0;
		O1(4,3) = 6.0/11.0;
		O1(5,0) = 1816283.0/549120.0;
		O1(5,1) = -2.0/9.0;
		O1(5,2) = -4.0/11.0;
		O1(5,3) = -1.0/6.0;
		O1(5,4) = -2561809.0/1647360.0;
		O1(6,0) = 0.0;
		O1(6,1) = 7.0/11.0;
		O1(6,2) = -2203.0/264.0;
		O1(6,3) = 10825.0/792.0;
		O1(6,4) = -85.0/12.0;
		O1(6,5) = 841.0/396.0;
		O1(7,0) = 1.0/400.0;
		O1(7,1) = 49.0/12.0;
		O1(7,2) = 43.0/6.0;
		O1(7,3) = -7.0/10.0;
		O1(7,4) = -85.0/12.0;
		O1(7,5) = -2963.0/1200.0;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(2,0) = -11.0/4.0;
		O2(2,1) = 11.0/4.0;
		O2(3,0) = -1228.0/2925.0;
		O2(3,1) = -92.0/225.0;
		O2(3,2) = 808.0/975.0;
		O2(4,0) = -2572.0/2805.0;
		O2(4,1) = 167.0/255.0;
		O2(4,2) = 199.0/136.0;
		O2(4,3) = -1797.0/1496.0;
		O2(5,0) = -1816283.0/274560.0;
		O2(5,1) = 253.0/36.0;
		O2(5,2) = -23.0/44.0;
		O2(5,3) = 76.0/3.0;
		O2(5,4) = -20775791.0/823680.0;
		O2(6,0) = 0.0;
		O2(6,1) = 107.0/132.0;
		O2(6,2) = 1289.0/88.0;
		O2(6,3) = -9275.0/792.0;
		O2(6,4) = 0.0;
		O2(6,5) = -371.0/99.0;
		O2(7,0) = -1.0/200.0;
		O2(7,1) = -137.0/24.0;
		O2(7,2) = -235.0/16.0;
		O2(7,3) = 1237.0/80.0;
		O2(7,4) = 0.0;
		O2(7,5) = 2963.0/600.0;
		omegas.push_back(O2);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 1.0/4.0;
		c(2) = 3.0/4.0;
		c(3) = 11.0/20.0;
		c(4) = 1.0/2.0;
		c(5) = 1.0;
		c(6) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1}, {2}, {3}, {4}, {5}, {6}};
	}
};

#endif