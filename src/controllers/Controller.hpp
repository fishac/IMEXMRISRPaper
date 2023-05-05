#ifndef CONTROLLER_DEFINED__
#define CONTROLLER_DEFINED__

#include <deque>

using namespace arma;

class Controller {
public:
	const char* name;
	double P; // The outer method/primary order.
	double p; // The inner method/embedding order.
	double tol; // The tolerance goal.
	double safety_factor; // Safety parameter.
	double* k1; // Array of parameter values for phi1.
	double* k2; // Array of parameter values for phi2.
	int num_H; // Number of H values to store.
	int num_M; // Number of M values to store.
	int num_tolf;
	int num_err1; // Number of error 1 values to store. 
	int num_err2; // Number of error 2 values to store.
	int iteration; // Iteration number. Useful for building up extended-history controllers.
	std::deque<double> H_array = {}; // H array. Deque was chosen for ease of cycling out old values, and element access by index.
	std::deque<int> M_array = {}; // M array. Deque was chosen for ease of cycling out old values, and element access by index.
	std::deque<double> tolf_array = {}; // M array. Deque was chosen for ease of cycling out old values, and element access by index.
	std::deque<double> err1_array = {}; // Error 1 array. Deque was chosen for ease of cycling out old values, and element access by index.
	std::deque<double> err2_array = {}; // Error 2 array. Deque was chosen for ease of cycling out old values, and element access by index.
	bool hm_controller = true;
	bool htol_controller = false;

	virtual void initialize(double H, int M) {};
	virtual void update_exponent_terms() {};
	virtual void update_values(double err1, double err2) {};
	virtual double get_new_H() { return 0.0; };
	virtual int get_new_M() { return 0; };
	virtual double get_new_tolf() { return 0.0; };
	virtual void replace_last_errors(double err1, double err2) {};
	virtual void replace_last_H(double H) {};
	virtual void replace_last_M(int M) {};
	virtual void replace_last_tolf(double tolf) {};

	virtual void update_H(double H_new) {
		H_array.push_back(H_new);
		H_array.pop_front();
	}

	virtual void update_M(int M_new) {
		M_array.push_back(M_new);
		M_array.pop_front();
	}
	
	virtual void update_tolf(double tolf_new) {
		tolf_array.push_back(tolf_new);
		tolf_array.pop_front();
	}

	virtual void update_errors(double err1, double err2) {
		if (err1_array.size() > 0) {
			err1_array.pop_front();
		}
		if (err2_array.size() > 0) {
			err2_array.pop_front();
		}
		err1_array.push_back(err1);
		err2_array.push_back(err2);
		iteration++;
	}

	virtual void reset() {
		if (hm_controller) {
			initialize(1.0, 1);
		} else if (htol_controller) {
			initialize(1.0, 1.0);
		}
		iteration = 0;
	}

	virtual void set_orders(double P_, double p_) {
		P = P_;
		p = p_;
		update_exponent_terms();
	}

	virtual void set_parameters(double* k1_, double* k2_) {
		k1 = k1_;
		k2 = k2_;
		update_exponent_terms();
	}

	virtual void set_tol(double tol_) {
		tol = tol_;
	}
	
	virtual void print_status() {
		printf("H:\n");
		for(double H : H_array) {
			printf("\t%.16f\n",H);
		}
		if (hm_controller) {
			printf("M:\n");
			for(int M : M_array) {
				printf("\t%d\n",M);
			}
		} else if (htol_controller) {
			printf("tolf:\n");
			for(double tolf : tolf_array) {
				printf("\t%.16f\n",tolf);
			}
		}
		printf("err1:\n");
		for(double err1 : err1_array) {
			printf("\t%.16f\n",err1);
		}
		printf("err2:\n");
		for(double err2 : err2_array) {
			printf("\t%.16f\n",err2);
		}
	}
};

#endif