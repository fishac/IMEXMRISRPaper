#ifndef IMEXMRISR21COEFFICIENTS_DEFINED__
#define IMEXMRISR21COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISR21Coefficients : public MRICoefficients {
public:
	IMEXMRISR21Coefficients() {
		name = "IMEXMRISR21";
		explicit_mrigark = false;
		method_type = 2;
		
		num_gammas = 1;
		num_omegas = 1;
		num_stages = 4;
		primary_order = 2.0;
		secondary_order = 1.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = -11.0/23.0;
		G1(1,1) = 11.0/23.0;
		G1(2,0) = -6692.0/52371.0;
		G1(2,1) = -18355.0/52371.0;
		G1(2,2) = 11.0/23.0;
		G1(3,0) = 11621.0/90666.0;
		G1(3,1) = -215249.0/226665.0;
		G1(3,2) = 17287.0/50370.0;
		G1(3,3) = 11.0/23.0;
		G1(4,0) = -31.0/12.0;
		G1(4,1) = -1.0/6.0;
		G1(4,2) = 11.0/4.0;
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 3.0/5.0;
		O1(2,0) = 14.0/165.0;
		O1(2,1) = 2.0/11.0;
		O1(3,0) = -13.0/54.0;
		O1(3,1) = 137.0/270.0;
		O1(3,2) = 11.0/15.0;
		O1(4,0) = -1.0/4.0;
		O1(4,1) = 1.0/2.0;
		O1(4,2) = 3.0/4.0;
		omegas.push_back(O1);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 3.0/5.0;
		c(2) = 4.0/15.0;
		c(3) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1}, {2}, {3}};
	}
};

#endif