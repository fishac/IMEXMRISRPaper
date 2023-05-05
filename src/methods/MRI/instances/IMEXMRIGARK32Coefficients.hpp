#ifndef IMEXMRIGARK32COEFFICIENTS_DEFINED__
#define IMEXMRIGARK32COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRIGARK32Coefficients: public MRICoefficients {
public:
	IMEXMRIGARK32Coefficients() {
		name = "IMEXMRIGARK32";
		explicit_mrigark = false;
		method_type = 0;
		
		num_gammas = 1;
		num_omegas = 1;
		num_stages = 8;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		mat G(num_stages+1,num_stages,fill::zeros);
		G(1,0) = 3.0/7.0;
		G(2,0) = -1.0;
		G(2,2) = 1.0;
		G(3,0) = -4.0/105.0;
		G(3,2) = 1.0/7.0;
		G(4,0) = 388.0/315.0;
		G(4,2) = -703.0/315.0;
		G(4,4) = 1.0;
		G(5,0) = -33997.0/92610.0;
		G(5,2) = 6178.0/9261.0;
		G(5,4) = 1.0/6.0;
		G(6,0) = 43461623.0/23245110.0;
		G(6,2) = -38315719.0/11622555.0;
		G(6,4) = 643.0/1506.0;
		G(6,6) = 1.0;
		G(7,0) = -14243.0/7530.0;
		G(7,2) = 99701.0/30120.0;
		G(7,4) = -1835.0/3514.0;
		G(7,6) = -531.0/280.0;
		G(7,7) = 1.0;
		G(8,0) = -2569.0/4518.0;
		G(8,2) = 2569.0/4518.0;
		gammas.push_back(G);

		mat O(num_stages+1,num_stages,fill::zeros);
		O(1,0) = 3.0/7.0;
		O(3,0) = -73.0/105.0;
		O(3,2) = 4.0/5.0;
		O(4,0) = 3119.0/7425.0;
		O(4,2) = -3119.0/7425.0;
		O(5,0) = -225086.0/1091475.0;
		O(5,2) = 491891.0/1091475.0;
		O(5,4) = 2.0/9.0;
		O(6,0) = -9450719.0/91320075.0;
		O(6,2) = -32058406.0/91320075.0;
		O(6,4) = 5.0/11.0;
		O(7,0) = 680411548.0/2132416935.0;
		O(7,2) = -401996371.0/1550848680.0;
		O(7,4) = -295928.0/1321551.0;
		O(7,6) = 87599.0/533960.0;
		O(8,0) = 2569.0/16566.0;
		O(8,2) = -2569.0/16566.0;
		omegas.push_back(O);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 3.0/7.0;
		c(2) = 3.0/7.0;
		c(3) = 8.0/15.0;
		c(4) = 8.0/15.0;
		c(5) = 1.0;
		c(6) = 1.0;
		c(7) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2, 3, 4, 5, 6, 7}};
	}
};

#endif