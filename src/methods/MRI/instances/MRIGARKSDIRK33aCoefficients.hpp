#ifndef MRIGARKSDIRK33aCOEFFICIENTS_DEFINED__
#define MRIGARKSDIRK33aCOEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MRIGARKSDIRK33aCoefficients: public MRICoefficients {
public:
	MRIGARKSDIRK33aCoefficients() {
		name = "MRIGARKSDIRK33a";
		explicit_mrigark = false;
		method_type = 0;
		
		num_gammas = 1;
		num_omegas = 0;
		num_stages = 8;
		primary_order = 3.0;
		secondary_order = 2.0;

		double lambda = 0.435866521508458999416019;
		double lambda4 = std::pow(lambda,4.0);
		double lambda3 = std::pow(lambda,3.0);
		double lambda2 = std::pow(lambda,2.0);
		
		mat G1(num_stages+1,num_stages,fill::zeros);
	    G1(1,0) = 1.0/3.0;
	    G1(2,0) = -lambda;
	    G1(2,2) = lambda;
	    G1(3,0) = (3.0-10.0*lambda)/(24.0*lambda-6.0);
	    G1(3,2) = (5.0-18.0*lambda)/(6.0-24.0*lambda);
	    G1(4,0) = (-24.0*lambda2+6.0*lambda+1.0)/(6.0-24.0*lambda);
	    G1(4,2) = (-48.0*lambda2+12.0*lambda+1.0)/(24.0*lambda-6.0);
	    G1(4,4) = lambda;
	    G1(5,0) = (3.0-16.0*lambda)/(12.0-48.0*lambda);
	    G1(5,2) = (48.0*lambda2-21.0*lambda+2.0)/(12.0*lambda-3.0);
	    G1(5,4) = (3.0-16.0*lambda)/4.0;
	    G1(6,0) = -lambda;
	    G1(6,6) = lambda;
	    G1(8,0) = 0.0;
	    G1(8,2) = (-6.0*lambda3+14.0*lambda2-7.0*lambda+1.0)/(12.0*lambda4-48.0*lambda3+60.0*lambda2-28.0*lambda+4.0);
	    G1(8,4) = 3.0*std::pow(2.0*lambda2-4.0*lambda+1.0,2.0)/(4.0*(3.0*lambda3-9.0*lambda2+6.0*lambda-1.0));
	    G1(8,6) = (2.0*lambda2-4.0*lambda+1.0)/(2.0-2.0*lambda);
		gammas.push_back(G1);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = lambda;
	    c(2) = lambda;
	    c(3) = (6.0*lambda2-9.0*lambda+2.0)/(6.0*lambda2-12.0*lambda+3.0);
	    c(4) = (6.0*lambda2-9.0*lambda+2.0)/(6.0*lambda2-12.0*lambda+3.0);
	    c(5) = 1.0;
	    c(6) = 1.0;
		c(7) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2, 3, 4, 5, 6, 7}};
	}
};

#endif