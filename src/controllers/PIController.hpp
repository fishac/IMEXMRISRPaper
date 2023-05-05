#ifndef PICONTROLLER_DEFINED__
#define PICONTROLLER_DEFINED__

#include "Controller.hpp"

using namespace arma;

class PIController : public Controller {
public:
	PIController(double P_, double p_, double tol_, double safety_factor_, double k1_[2], double k2_[2])
	{
		name = "PI";
		P = P_;
		p = p_;
		tol = tol_;
		safety_factor = safety_factor_;
		k1 = k1_;
		k2 = k2_;
		num_H = 1;
		num_M = 1;
		num_err1 = 2;
		num_err2 = 2;
	}

	void update_errors(double err1, double err2) {
		if (err1_array.size() > 0) {
			err1_array.pop_front();
		}
		if (err2_array.size() > 0) {
			err2_array.pop_front();
		}
		err1_array.push_back(err1);
		err2_array.push_back(err2);
	}

	double get_new_H() {
		double errn_pow = k1[0]/P;
		double errnm1_pow = -k2[0]/P;
		double Hn = H_array[0];
		double Hnp1 = safety_factor*Hn*std::pow(tol/err1_array[0],errnm1_pow)*std::pow(tol/err1_array[1],errn_pow);
		return Hnp1;
	}

	int get_new_M() {
		return M_array[0];
	}

	void initialize(double H, int M) {
		H_array = { H };
		M_array = { M };
		err1_array = { 1.0, 1.0 };
		err2_array = { 1.0, 1.0 };
	}
	
	void replace_last_errors(double err1, double err2) {
		err1_array[1] = err1;
		err2_array[1] = err2;
	}

	void replace_last_H(double H) {
		H_array[0] = H;
	}

	void replace_last_M(int M) {
		M_array[0] = M;
	}
};

#endif