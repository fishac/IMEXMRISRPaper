#ifndef MRIGARKTEST2COEFFICIENTS_DEFINED__
#define MRIGARKTEST2COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class MRIGARKTEST2Coefficients: public MRICoefficients {
public:
	MRIGARKTEST2Coefficients() {
		name = "MRIGARKTEST2";
		explicit_mrigark = true;
		method_type = 0;
		
		num_gammas = 2;
		num_omegas = 0;
		num_stages = 4;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		mat G1(num_stages+1,num_stages,fill::zeros);
		G1(1,0) = 0.5807851970537204;
		G1(2,0) = -0.08803151336500341;
		G1(2,1) = 0.22187482965628463;
		G1(3,0) = -0.20577456214553907;
		G1(3,1) = -1.9394255390998725;
		G1(3,2) = 2.4305715879004337;
		G1(4,0) = 1.0212332618433528;
		G1(4,1) = -2.0934260666270346;
		G1(4,2) = 1.35756429143868;

		mat G2(num_stages+1,num_stages,fill::zeros);
		G2(2,0) = -0.8346621210667525;
		G2(2,1) = 0.8346621210667525;
		G2(3,0) = 0.7458174416586598;
		G2(3,1) = 3.2174371738696617;
		G2(3,2) = -3.9632546155283426;
		G2(4,0) = -1.7081982063191243;
		G2(4,1) = 3.5254382289239787;
		G2(4,2) = -1.817240022604861;

		gammas.push_back(G1);
		gammas.push_back(G2);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.5807851970537204;
		c(2) = 0.7146285133450017;
		c(3) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2, 3}};
	}
};

#endif