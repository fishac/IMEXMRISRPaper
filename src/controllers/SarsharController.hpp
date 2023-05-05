#ifndef SARSHARCONTROLLER_DEFINED__
#define SARSHARCONTROLLER_DEFINED__

#include "Controller.hpp"

using namespace arma;

class SarsharController : public Controller {
public:
	SarsharController(double P_, double p_, double tol_, double safety_factor_, double k1_[1], double k2_[1]) {
		name = "SarsharController";
		P = P_;
		p = p_;
		tol = tol_;
		safety_factor = safety_factor_;
	}

	double get_new_H() {
		double Hn = H_array[0];
		double Hnp1 = safety_factor*Hn*std::pow(err1_array[0],-1.0/P);
		return Hnp1;
	}

	int get_new_M() {
		int Mn = M_array[0];
		int Mnp1 = (int)ceil(Mn*pow(err2_array[0]/err1_array[0],1.0/p));
		return Mnp1;
	}

	void initialize(double H, int M) {
		H_array = { H };
		M_array = { M };
		err1_array = { 1.0 };
		err2_array = { 1.0 };
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