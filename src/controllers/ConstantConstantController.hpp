#ifndef CONSTANTCONSTANTCONTROLLER_DEFINED__
#define CONSTANTCONSTANTCONTROLLER_DEFINED__

#include "Controller.hpp"

using namespace arma;

class ConstantConstantController : public Controller {
public:
	double Hnm1_pow;
	double Hn_pow;
	double H_essn_pow;
	double Mnm2_pow;
	double Mnm1_pow;
	double Mn_pow;
	double M_essnm1_pow;
	double M_essn_pow;
	double M_esfnm1_pow;
	double M_esfn_pow;
	double ess_split = 0.5;
	double esf_split = 0.5;

	ConstantConstantController(double P_, double p_, double tol_, double safety_factor_, double* k1_, double* k2_) {
		name = "ConstantConstant";
		P = P_;
		p = p_;
		tol = tol_;
		safety_factor = safety_factor_;
		k1 = k1_;
		k2 = k2_;
		update_exponent_terms();
	}

	void update_exponent_terms() {
		Hn_pow = 1.0; 

		Mn_pow = 1.0;

		H_essn_pow = k1[0]/P;

		M_essn_pow = k1[0]*(1.0+p)/(p*P);
		M_esfn_pow = -k2[0]/p;
	}

	double get_new_H() override {
		double Hn = H_array[0];
		//printf("Hn: %.16f, essn: %.16f, H_essn_pow: %.16f\n",Hn,err1_array[0],H_essn_pow);
		double H_new = safety_factor*std::pow(Hn,Hn_pow)*std::pow(ess_split*tol/err1_array[0],H_essn_pow);
		return H_new;
	}

	int get_new_M() {
		int Mn = M_array[0];
		double M_terms = std::pow(Mn,Mn_pow);
		double err1_terms = std::pow(ess_split*tol/err1_array[0],M_essn_pow);
		double err2_terms = std::pow(esf_split*tol/err2_array[0],M_esfn_pow);
		int M_new = std::ceil(safety_factor*M_terms*err1_terms*err2_terms);
		return M_new;
	}

	void initialize(double H, int M) {
		H_array = { H };
		M_array = { M };
		err1_array = { ess_split*tol };
		err2_array = { esf_split*tol };
	}

	void reset() {
		initialize(1.0, 1);
		iteration = 0;
	}

	void replace_last_errors(double err1, double err2) {
		err1_array[0] = err1;
		err2_array[0] = err2;
	}

	void replace_last_H(double H) {
		H_array[0] = H;
	}

	void replace_last_M(int M) {
		M_array[0] = M;
	}
};

#endif