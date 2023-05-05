#ifndef MRIGARKESDIRK46aCOEFFICIENTS_DEFINED__
#define MRIGARKESDIRK46aCOEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MRIGARKESDIRK46aCoefficients: public MRICoefficients {
public:
	MRIGARKESDIRK46aCoefficients() {
		name = "MRIGARKESDIRK46a";
		explicit_mrigark = false;
		method_type = 0;
		
		num_gammas = 2;
		num_omegas = 0;
		num_stages = 12;
		primary_order = 3.0;
		secondary_order = 2.0;

		mat G1(num_stages+1,num_stages,fill::zeros);
	    G1(1,0) = 1.0/5.0;
	    G1(2,0) = -1.0/4.0;
	    G1(2,2) = 1.0/4.0;
	    G1(3,0) = 1771023115159.0/1929363690800.0;
	    G1(3,2) = -1385150376999.0/1929363690800.0;
	    G1(4,0) = 914009.0/345800.0;
	    G1(4,2) = -1000459.0/345800.0;
	    G1(4,4) = 1.0/4.0;
	    G1(5,0) = 18386293581909.0/36657910125200.0;
	    G1(5,2) = 5506531089.0/80566835440.0;
	    G1(5,4) = -178423463189.0/482340922700.0;
	    G1(6,0) = 36036097.0/8299200.0;
		G1(6,2) = 4621.0/118560.0;
		G1(6,4) = -38434367.0/8299200.0;
	    G1(6,6) = 1.0/4.0;
	    G1(7,0) = -247809665162987.0/146631640500800.0;
	    G1(7,2) = 10604946373579.0/14663164050080.0;
	    G1(7,4) = 10838126175385.0/5865265620032.0;
	    G1(7,6) = -24966656214317.0/36657910125200.0;
		G1(8,0) = 38519701.0/11618880.0;
		G1(8,2) = 10517363.0/9682400.0;
		G1(8,4) = -23284701.0/19364800.0;
		G1(8,6) = -10018609.0/2904720.0;
		G1(8,8) = 1.0/4.0;
		G1(9,0) = -52907807977903.0/33838070884800.0;
		G1(9,2) = 74846944529257.0/73315820250400.0;
		G1(9,4) = 365022522318171.0/146631640500800.0;
		G1(9,6) = -20513210406809.0/109973730375600.0;
		G1(9,8) = -2918009798.0/1870301537.0;
		G1(10,0) = 19.0/100.0;
		G1(10,2) = -73.0/300.0;
		G1(10,4) = 127.0/300.0;
		G1(10,6) = 127.0/300.0;
		G1(10,8) = -313.0/300.0;
		G1(10,10) = 1.0/4.0;
		G1(12,0) = -1.0/4.0;
		G1(12,2) = 5595.0/8804.0;
		G1(12,4) = -2445.0/8804.0;
		G1(12,6) = -4225.0/8804.0;
		G1(12,8) = 2205.0/4402.0;
		G1(12,10) = -567.0/4402.0;
		gammas.push_back(G1);
		
		mat G2(num_stages+1,num_stages,fill::zeros);
	    G2(3,0) = -1674554930619.0/964681845400.0;
	    G2(3,2) = 1674554930619.0/964681845400.0;
	    G2(4,0) = -1007739.0/172900.0;
	    G2(4,2) = 1007739.0/172900.0;
	    G2(5,0) = -8450070574289.0/18328955062600.0;
	    G2(5,2) = -39429409169.0/40283417720.0;
	    G2(5,4) = 173621393067.0/120585230675.0;
	    G2(6,0) = -122894383.0/16598400.0;
		G2(6,2) = 14501.0/237120.0;
		G2(6,4) = 121879313.0/16598400.0;
	    G2(7,0) = 32410002731287.0/15434909526400.0;
	    G2(7,2) = -46499276605921.0/29326328100160.0;
	    G2(7,4) = -34914135774643.0/11730531240064.0;
	    G2(7,6) = 45128506783177.0/18328955062600.0;
		G2(8,0) = -128357303.0/23237760.0;
		G2(8,2) = -35433927.0/19364800.0;
		G2(8,4) = 71038479.0/38729600.0;
		G2(8,6) = 8015933.0/1452360.0;
		G2(9,0) = 136721604296777.0/67676141769600.0;
		G2(9,2) = -349632444539303.0/146631640500800.0;
		G2(9,4) = -1292744859249609.0/293263281001600.0;
		G2(9,6) = 8356250416309.0/54986865187800.0;
		G2(9,8) = 17282943803.0/3740603074.0;
		G2(10,0) = 3.0/25.0;
		G2(10,2) = -29.0/300.0;
		G2(10,4) = 71.0/300.0;
		G2(10,6) = 71.0/300.0;
		G2(10,8) = -149.0/300.0;
		gammas.push_back(G2);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 1.0/5.0;
	    c(2) = 1.0/5.0;
	    c(3) = 2.0/5.0;
	    c(4) = 2.0/5.0;
		c(5) = 3.0/5.0;
		c(6) = 3.0/5.0;
		c(7) = 4.0/5.0;
		c(8) = 4.0/5.0;
		c(9) = 1.0;
		c(10) = 1.0;
		c(11) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};
	}
};

#endif