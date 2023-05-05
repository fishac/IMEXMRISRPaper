#ifndef GUSTAFSSONCONTROLLER_DEFINED__
#define GUSTAFSSONCONTROLLER_DEFINED__

#include "Controller.hpp"

using namespace arma;

class GustafssonController : public Controller {
public:
	double Hnm1_pow;
	double Hn_pow;
	double Hn_pow_initial;
	double H_essnm1_pow;
	double H_essn_pow;
	double Mnm2_pow;
	double Mnm1_pow;
	double Mn_pow;
	double Mn_pow_initial;
	double M_essnm1_pow;
	double M_essn_pow;
	double M_esfnm1_pow;
	double M_esfn_pow;
	double ess_split = 0.5;
	double esf_split = 0.5;
	Controller* initial_controller;

	GustafssonController(double P_, double p_, double tol_, double safety_factor_, double k1_[2], double k2_[1], Controller* initial_controller_) {
		name = "Gustafsson";
		P = P_;
		p = p_;
		tol = tol_;
		safety_factor = safety_factor_;
		k1 = k1_;
		k2 = k2_;
		initial_controller = initial_controller_;
		update_exponent_terms();
	}

	void update_H(double H_new) {
		H_array.push_back(H_new);
		H_array.pop_front();
		if (iteration <= 1) {
			initial_controller->update_H(H_new);
		}	
	}

	void update_M(int M_new) {
		M_array.push_back(M_new);
		M_array.pop_front();
		if (iteration <= 1) {
			initial_controller->update_M(M_new);
		}	
	}

	void update_exponent_terms() {
		Hnm1_pow = -1.0;
		Hn_pow = 2.0; 

		H_essnm1_pow = -k1[0]/P;
		H_essn_pow = (k1[0]+k1[1])/(P);
	}

	double get_new_H() {
		double Hnm1 = H_array[0];
		double Hn = H_array[1];
		double H_new;
		if (iteration == 1) {
			H_new = initial_controller->get_new_H();
		} else {
			H_new = safety_factor*std::pow(Hnm1,Hnm1_pow)*std::pow(Hn,Hn_pow)*std::pow(tol/err1_array[0],H_essnm1_pow)*std::pow(tol/err1_array[1],H_essn_pow);
		} 
		return H_new;
	}

	int get_new_M() {
		return M_array[0];
	}

	void initialize(double H, int M) {
		initial_controller->initialize(H, M);
		H_array = { 1.0, H };
		M_array = { M };
		err1_array = { tol, tol };
		err2_array = { tol };
		iteration = 0;
	}

	void set_orders(double P_, double p_) {
		P = P_;
		p = p_;
		update_exponent_terms();
		initial_controller->set_orders(P,p);
	}

	void update_errors(double err1, double err2) {
		err1_array.push_back(err1);
		err2_array.push_back(err2);
		err1_array.pop_front();
		err2_array.pop_front();
		initial_controller->update_errors(err1,err2);
		iteration++;
	}

	void reset() {
		H_array = {};
		M_array = {};
		err1_array = {};
		err2_array = {};
		iteration = 0;
		initial_controller->reset();
	}
};

#endif