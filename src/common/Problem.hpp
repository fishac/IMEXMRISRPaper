#ifndef PROBLEM_DEFINED__
#define PROBLEM_DEFINED__

using namespace std;
using namespace arma;

class Problem {
public:
	const char* name;
	int problem_dimension = 0;
	double default_H;
	double t_0;
	double t_f;
	bool has_true_solution;
	bool explicit_only;
	vec y_0;
	vec yhat;
	double that;
	int full_function_evals = 0;
	int fast_function_evals = 0;
	int slow_function_evals = 0;
	int implicit_function_evals = 0;
	int explicit_function_evals = 0;
	int full_jacobian_evals = 0;
	int fast_jacobian_evals = 0;
	int slow_jacobian_evals = 0;
	int implicit_jacobian_evals = 0;
	int slow_nonlinear_solves = 0;
	double time_step_multiplier = 1.0;
	
	virtual void full_rhs_custom(double t, vec* y, vec* f) {}
	
	virtual void fast_rhs_custom(double t, vec* y, vec* f) {}
	
	virtual void slow_rhs_custom(double t, vec* y, vec* f) {}
	
	virtual void implicit_rhs_custom(double t, vec* y, vec* f) {}
	
	virtual void explicit_rhs_custom(double t, vec* y, vec* f) {}
	
	virtual void linear_rhs_custom(double t, vec* y, vec* f) {}
	
	virtual void nonlinear_rhs_custom(double t, vec* y, vec* f) {}
	
	virtual void full_rhsjacobian_custom(double t, vec* y, mat* j) {}
	
	virtual void fast_rhsjacobian_custom(double t, vec* y, mat* j) {}
	
	virtual void slow_rhsjacobian_custom(double t, vec* y, mat* j) {}
	
	virtual void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {}
	
	virtual void linear_rhsjacobian_custom(double t, vec* y, mat* j) {}
	
	virtual void true_solution(double t, vec* y) {}

	void full_rhs(double t, vec* y, vec* f) {
		full_function_evals++;
		full_rhs_custom(t,y,f);
	}
	
	void fast_rhs(double t, vec* y, vec* f) {
		fast_function_evals++;
		fast_rhs_custom(t,y,f);
	}
	
	void slow_rhs(double t, vec* y, vec* f) {
		slow_function_evals++;
		slow_rhs_custom(t,y,f);
	}
	
	void implicit_rhs(double t, vec* y, vec* f) {
		implicit_function_evals++;
		implicit_rhs_custom(t,y,f);
	}
	
	void explicit_rhs(double t, vec* y, vec* f) {
		explicit_function_evals++;
		explicit_rhs_custom(t,y,f);
	}
	
	void linear_rhs(double t, vec* y, vec* f) {
		fast_function_evals++;
		linear_rhs_custom(t,y,f);
	}
	
	void nonlinear_rhs(double t, vec* y, vec* f) {
		slow_function_evals++;
		nonlinear_rhs_custom(t,y,f);
	}
	
	void full_rhsjacobian(double t, vec* y, mat* j) {
		full_jacobian_evals++;
		full_rhsjacobian_custom(t,y,j);
	}
	
	void fast_rhsjacobian(double t, vec* y, mat* j) {
		fast_jacobian_evals++;
		fast_rhsjacobian_custom(t,y,j);
	}
	
	void slow_rhsjacobian(double t, vec* y, mat* j) {
		slow_jacobian_evals++;
		slow_rhsjacobian_custom(t,y,j);
	}
	
	void implicit_rhsjacobian(double t, vec* y, mat* j) {
		implicit_jacobian_evals++;
		implicit_rhsjacobian_custom(t,y,j);
	}
	
	void linear_rhsjacobian(double t, vec* y, mat* j) {
		fast_jacobian_evals++;
		linear_rhsjacobian_custom(t,y,j);
	}

	void reset_eval_counts() {
		full_function_evals = 0;
		fast_function_evals = 0;
		slow_function_evals = 0;
		implicit_function_evals = 0;
		explicit_function_evals = 0;
		full_jacobian_evals = 0;
		fast_jacobian_evals = 0;
		slow_jacobian_evals = 0;
		implicit_jacobian_evals = 0;
		slow_nonlinear_solves = 0;
	}
	
	void increment_slow_nonlinear_solves() {
		slow_nonlinear_solves++;
	}

	void set_yhat_that(vec yhat_, double that_) {
		yhat = yhat_;
		that = that_;
	}
};

#endif