#ifndef PIDMRCONTROLLER_DEFINED__
#define PIDMRCONTROLLER_DEFINED__

#include "Controller.hpp"

using namespace arma;

class PIDMRController : public Controller {
public:
	double Hnm2_pow;
	double Hnm1_pow;
	double Hn_pow;
	double H_essnm2_pow;
	double H_essnm1_pow;
	double H_essn_pow;
	double Mnm2_pow;
	double Mnm1_pow;
	double Mn_pow;
	double M_essnm2_pow;
	double M_essnm1_pow;
	double M_essn_pow;
	double M_esfnm2_pow;
	double M_esfnm1_pow;
	double M_esfn_pow;
	double ess_split = 0.5;
	double esf_split = 0.5;
	Controller* initial_controller;

	PIDMRController(double P_, double p_, double tol_, double safety_factor_, double k1_[3], double k2_[3], Controller* initial_controller_) {
		name = "PIDMR";
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
		if (iteration <= 3) {
			initial_controller->update_H(H_new);
		}	
	}

	void update_M(int M_new) {
		M_array.push_back(M_new);
		M_array.pop_front();
		if (iteration <= 3) {
			initial_controller->update_M(M_new);
		}	
	}

	void update_exponent_terms() {
		
		H_essnm2_pow = k1[0]/(3.0*P);
		H_essnm1_pow = -(k1[0]+k1[1])/(3.0*P);
		H_essn_pow = (k1[0]+k1[1]+k1[2])/(3.0*P);

		//printf("k1[0]: %.3f, k1[1]: %.3f, k1[2]: %.3f, P: %.0f\n",k1[0],k1[1],k1[2],P);
		//printf("H_essnm2_pow: %.4f, H_essnm1_pow: %.4f, H_essn_pow: %.4f\n",H_essnm2_pow,H_essnm1_pow,H_essn_pow);
		
		M_essnm2_pow = k1[0]*(1.0+p)/(3.0*P*p);
		M_essnm1_pow = -(k1[0]+k1[1])*(1.0+p)/(3.0*P*p);
		M_essn_pow = (k1[0]+k1[1]+k1[2])*(1.0+p)/(3.0*P*p);
		
		M_esfnm2_pow = -k2[0]/(3.0*p);
		M_esfnm1_pow = (k2[0]+k2[1])/(3.0*p);
		M_esfn_pow = -(k2[0]+k2[1]+k2[2])/(3.0*p);
		
		/*
		H_essnm2_pow = k1[0]/(3.0*P);
		H_essnm1_pow = (k1[0]+k1[1]-8.0)/(3.0*P);
		H_essn_pow = -(k1[0]+k1[1]+k1[2]-17.0/2.0)/(3.0*P);

		//printf("k1[0]: %.3f, k1[1]: %.3f, k1[2]: %.3f, P: %.0f\n",k1[0],k1[1],k1[2],P);
		//printf("H_essnm2_pow: %.4f, H_essnm1_pow: %.4f, H_essn_pow: %.4f\n",H_essnm2_pow,H_essnm1_pow,H_essn_pow);
		
		M_essnm2_pow = -k1[0]*(1.0+p)/(3.0*P*p);
		M_essnm1_pow = -(k1[0]+k1[1])*(1.0+p)/(3.0*P*p);
		M_essn_pow = (k1[0]+k1[1]+k1[2])*(1.0+p)/(3.0*P*p);
		
		M_esfnm2_pow = k2[0]/(3.0*p);
		M_esfnm1_pow = (k2[0]+k2[1])/(3.0*p);
		M_esfn_pow = -(k2[0]+k2[1]+k2[2])/(3.0*p);
		*/
	}

	double get_new_H() {
		double Hn = H_array[0];
		double H_new;
		if (iteration < 3) {
			H_new = initial_controller->get_new_H();
		} else {
			//printf("\tHn: %.16f, Hnm1: %.16f, Hnm2: %.16f, essn: %.16f, essnm1: %.16f, essnm2: %.16f\n",Hn,Hnm1,Hnm2,err1_array[2],err1_array[1],err1_array[0]);
			double H_terms = Hn;
			double err1_terms = std::pow(ess_split*tol/err1_array[0],H_essnm2_pow)*std::pow(ess_split*tol/err1_array[1],H_essnm1_pow)*std::pow(ess_split*tol/err1_array[2],H_essn_pow);
			H_new = safety_factor*H_terms*err1_terms;
			//printf("%.16f = %.2f(%.16f)^%.0f(%.16f)^%.0f(%.16f)^%.0f(%.1f*%.1f/%.16f)^%.3f(%.1f*%.1f/%.16f)^%.3f(%.1f*%.1f/%.16f)^%.3f\n",
			//	H_new,safety_factor,Hn,Hn_pow,Hnm1,Hnm1_pow,Hnm2,Hnm2_pow,ess_split,tol,err1_array[0],H_essnm2_pow,ess_split,tol,err1_array[1],H_essnm1_pow,ess_split,tol,err1_array[2],H_essn_pow);
		} 
		return H_new;
	}

	int get_new_M() {
		int Mn = M_array[0];
		int M_new;
		if (iteration < 3) {
			M_new = initial_controller->get_new_M();
		} else {
			double M_terms = (double) Mn;
			double err1_terms = std::pow(ess_split*tol/err1_array[0],M_essnm2_pow)*std::pow(ess_split*tol/err1_array[1],M_essnm1_pow)*std::pow(ess_split*tol/err1_array[2],M_essn_pow);
			double err2_terms = std::pow(esf_split*tol/err2_array[0],M_esfnm2_pow)*std::pow(esf_split*tol/err2_array[1],M_esfnm1_pow)*std::pow(esf_split*tol/err2_array[2],M_esfn_pow);
			M_new = std::ceil(safety_factor*M_terms*err1_terms*err2_terms);
		}
		return M_new;
	}

	void initialize(double H, int M) {
		H_array = { H };
		M_array = { M };
		err1_array = { ess_split*tol, ess_split*tol, ess_split*tol };
		err2_array = { esf_split*tol, esf_split*tol, esf_split*tol };
		iteration = 0;
		initial_controller->initialize(H, M);
	}

	void set_orders(double P_, double p_) {
		P = P_;
		p = p_;
		update_exponent_terms();
		initial_controller->set_orders(P,p);
	}

	void set_tol(double tol_) {
		tol = tol_;
		initial_controller->set_tol(tol_);
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
		initialize(1.0, 1);
		iteration = 0;
		initial_controller->reset();
	}

	void replace_last_errors(double err1, double err2) {
		err1_array[2] = err1;
		err2_array[2] = err2;
		if (iteration < 3) {
			initial_controller->replace_last_errors(err1, err2);
		}
	}

	void replace_last_H(double H) {
		H_array[0] = H;
		if (iteration < 3) {
			initial_controller->replace_last_H(H);
		}
	}

	void replace_last_M(int M) {
		M_array[0] = M;
		if (iteration < 3) {
			initial_controller->replace_last_M(M);
		}
	}
};

#endif