#ifndef IMEXMRISRTEST3COEFFICIENTS_DEFINED__
#define IMEXMRISRTEST3COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRISRTEST3Coefficients : public MRICoefficients {
public:
	IMEXMRISRTEST3Coefficients() {
		name = "IMEXMRISRTEST3";
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
		G1(2,0) = -0.084;
		G1(2,1) = -0.166;
		G1(2,2) = 0.25;
		G1(3,0) = 0.19348346118010076;
		G1(3,1) = -0.046211255286943741;
		G1(3,2) = -0.39727220589315702;
		G1(3,3) = 0.25;
		G1(4,0) = 0.25367564170848026;
		G1(4,1) = -0.23483923299747126;
		G1(4,2) = -0.24860482604014310;
		G1(4,3) = -0.020231582670865905;
		G1(4,4) = 0.25;
		G1(5,0) = -0.043508055511004974;
		G1(5,1) = -0.0087420578429041841;
		G1(5,2) = 0.026818983452319619;
		G1(5,3) = 0.27673623478725708;
		G1(5,4) = -0.50130510488566755;
		G1(5,5) = 0.25;
		gammas.push_back(G1);
		
		mat O1(num_stages+1,num_stages,fill::zeros);
		O1(1,0) = 0.5;
		O1(2,0) = 0.332;
		O1(3,0) = 0.62;
		O1(4,0) = 0.85;
		O1(5,0) = 1.0;
		O1(6,0) = 0.80174392120469429;
		O1(6,3) = 0.021178493241912901;
		O1(6,4) = 1.2680550090891921;
		O1(6,5) = -1.0909774235357992;
		omegas.push_back(O1);
		
		mat O2(num_stages+1,num_stages,fill::zeros);
		O2(2,0) = -0.22044800000000000;
		O2(2,1) = 0.22044800000000000;
		O2(3,0) = -1.3376931903062372;
		O2(3,1) = -0.35544130465280200;
		O2(3,2) = 1.6931344949590392;
		O2(4,0) = -2.0108337168498310;
		O2(4,1) = -0.71341001964439826;
		O2(4,2) = 2.1174517597368854;
		O2(4,3) = 0.60679197675734383;
		O2(5,0) = -1.5971512986546473;
		O2(5,1) = 0.017484115685808368;
		O2(5,2) = 0.31987991414336229;
		O2(5,3) = 0.80765812104415499;
		O2(5,4) = 0.45212914778132169;
		O2(6,0) = -1.2876552520860459;
		O2(6,2) = 0.37351788104800153;
		O2(6,3) = 1.3187736041348434;
		O2(6,4) = -3.0865910801683975;
		O2(6,5) = 2.6819548470715985;
		omegas.push_back(O2);
		

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.5;
		c(2) = 0.332;
		c(3) = 0.62;
		c(4) = 0.85;
		c(5) = 1.0;
		c(6) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1}, {2}, {3}, {4}, {5}, {6}};
	}
};

#endif