#ifndef IMEXMRITest1COEFFICIENTS_DEFINED__
#define IMEXMRITest1COEFFICIENTS_DEFINED__

#include "MRICoefficients.hpp"

using namespace arma;

class IMEXMRITest1Coefficients: public MRICoefficients {
public:
	IMEXMRITest1Coefficients() {
		name = "IMEXMRITest1";
		explicit_mrigark = false;
		method_type = 0;
		
		num_gammas = 1;
		num_omegas = 1;
		num_stages = 8;
		primary_order = 3.0;
		secondary_order = 2.0;
		
		mat G(num_stages+1,num_stages,fill::zeros);
		G(1,0) = 0.1458858459507417;
		G(2,0) = -0.8424794539487297;
		G(2,2) = 0.8424794539487297;
		G(3,0) = 0.09556159730913771;
		G(3,2) = 0.399008980903116;
		G(4,0) = 1.6241679444646575;
		G(4,2) = -2.466647398413387;
		G(4,4) = 0.8424794539487297;
		G(5,0) = 0.0015854893558506947;
		G(5,2) = 0.8937617664977989;
		G(5,4) = -0.5358036800166469;
		G(6,0) = 3.8687581642396314;
		G(6,2) = -3.982704458261963;
		G(6,4) = -0.7285331599263944;
		G(6,6) = 0.8424794539487297;
		G(7,0) = -9.074947007248928;
		G(7,2) = 11.017408837536848;
		G(7,4) = -2.481607433141217;
		G(7,6) = -0.3033338510954216;
		G(7,7) = 0.8424794539487297;
		G(8,0) = 3.0854848312877583;
		G(8,1) = -15.532711966649634;
		G(8,2) = 10.823835619553055;
		G(8,3) = 3.7430150146786687;
		G(8,4) = -2.687885835815083;
		G(8,5) = 0.391346863822345;
		G(8,6) = 0.17691547312288855;
		/*G(7,0) = 3.0854848312877583;
		G(7,1) = -15.532711966649634;
		G(7,2) = 10.823835619553055;
		G(7,3) = 3.7430150146786687;
		G(7,4) = -2.687885835815083;
		G(7,5) = 0.391346863822345;
		G(7,6) = 0.17691547312288855;*/
		gammas.push_back(G);

		mat O(num_stages+1,num_stages,fill::zeros);
		O(1,0) = 0.1458858459507417;
		O(3,0) = 1.0686231936286892;
		O(3,2) = -0.5740526154164356;
		O(4,0) = -1.8786984593376057;
		O(4,2) = 1.8786984593376057;
		O(5,0) = -4.066756330132985;
		O(5,2) = 4.3941402293314535;
		O(5,4) = 0.032159676638538054;
		O(6,0) = 5.149299364853728;
		O(6,2) = -6.063312981549586;
		O(6,4) = 0.9140136166958586;
		O(7,0) = -0.9059963833122844;
		O(7,2) = 1.4674003180402257;
		O(7,4) = -0.8191706654689832;
		O(7,6) = 0.25776673074103884;
		O(8,0) = -0.9059964141348158;
		O(8,2) = 1.4674003505877729;
		O(8,4) = -0.8191706570605083;
		O(8,6) = 0.25776672060755074;
		/*O(7,0) = -0.9059964141348158;
		O(7,2) = 1.4674003505877729;
		O(7,4) = -0.8191706570605083;
		O(7,6) = 0.25776672060755074;*/
		omegas.push_back(O);

		c = vec(num_stages,fill::zeros);
		c(0) = 0.0;
		c(1) = 0.1458858459507417;
		c(2) = 0.1458858459507417;
		c(3) = 0.6404564241629954;
		c(4) = 0.6404564241629954;
		c(5) = 1.0;
		c(6) = 1.0;
		c(7) = 1.0;
		
		num_groups = num_stages-1;
		stage_groups = {{1, 2, 3, 4, 5, 6, 7}};
	}
};

#endif