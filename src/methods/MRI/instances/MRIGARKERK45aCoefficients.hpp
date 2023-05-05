#ifndef MRIGARKERK45aCOEFFICIENTS_DEFINED__
#define MRIGARKERK45aCOEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MRIGARKERK45aCoefficients: public MRICoefficients {
public:
	MRIGARKERK45aCoefficients() {
		name = "MRIGARKERK45a";
		explicit_mrigark = true;
		method_type = 0;
		
		num_gammas = 2;
		num_omegas = 0;
		num_stages = 6;
		primary_order = 4.0;
		secondary_order = 3.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = 0.2;
	    G1(2,0) = -53.0/16.0;
	    G1(2,1) = 281.0/80.0;
	    G1(3,0) = -36562993.0/71394880.0;
	    G1(3,1) = 34903117.0/17848720.0;
	    G1(3,2) = -88770499.0/71394880.0;
	    G1(4,0) = -7631593.0/71394880.0;
	    G1(4,1) = -166232021.0/35697440.0;
	    G1(4,2) = 6068517.0/1519040.0;
	    G1(4,3) = 8644289.0/8924360.0;
	    G1(5,0) = 277061.0/303808.0;
	    G1(5,1) = -209323.0/1139280.0;
	    G1(5,2) = -1360217.0/1139280.0;
	    G1(5,3) = -148789.0/56964.0;
	    G1(5,4) = 147889.0/45120.0;
	    G1(6,0) = -1482837.0/759520.0;
	    G1(6,1) = 175781.0/71205.0;
	    G1(6,2) = -790577.0/1139280.0;
	    G1(6,3) = -6379.0/56964.0;
	    G1(6,4) = 47.0/96.0;
		gammas.push_back(G1);

		mat G2(num_stages+1,num_stages,fill::zeros);
	    G2(2,0) = 503.0/80.0;
	    G2(2,1) = -503.0/80.0;
	    G2(3,0) = -1365537.0/35697440.0;
	    G2(3,1) = 4963773.0/7139488.0;
	    G2(3,2) = -1465833.0/2231090.0;
	    G2(4,0) = 66974357.0/35697440.0;
	    G2(4,1) = 21445367.0/7139488.0;
	    G2(4,2) = -3.0;
	    G2(4,3) = -8388609.0/4462180.0;
	    G2(5,0) = -18227.0/7520.0;
	    G2(5,1) = 2.0;
	    G2(5,2) = 1.0;
	    G2(5,3) = 5.0;
	    G2(5,4) = -41933.0/7520.0;
	    G2(6,0) = 6213.0/1880.0;
	    G2(6,1) = -6213.0/1880.0;
		gammas.push_back(G2);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.2;
		c(2) = 0.4;
		c(3) = 0.6;
		c(4) = 0.8;
		c(5) = 1.0;
		
		num_groups = num_stages;
		stage_groups = {{1, 2, 3, 4, 5}};
	}
};

#endif