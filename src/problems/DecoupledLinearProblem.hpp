#ifndef DECOUPLEDLINEARPROBLEM_DEFINED__
#define DECOUPLEDLINEARPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class DecoupledLinearProblem : public Problem {
public:
	double lambda_f = -1.0;
	double lambda_s = -10.0;

	DecoupledLinearProblem() {
		name = "DecoupledLinear";
		problem_dimension = 2;
		default_H = std::pow(2.0,-6.0);
		t_0 = 0.0;
		t_f = 1.0;
		has_true_solution = true;
		explicit_only = false;
		y_0 = { 1.0, 1.0 };
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);

		(*f)(0) = lambda_f*y0;
		(*f)(1) = lambda_s*y1;
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);

		(*f)(0) = lambda_f*y0;
		(*f)(1) = 0.0;
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = lambda_s*y1;
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = (lambda_s/2.0)*y1;
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = (lambda_s/2.0)*y1;
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		
		(*j)(0,0) = lambda_f;
		(*j)(0,1) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = lambda_s;
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		
		(*j)(0,0) = lambda_f;
		(*j)(0,1) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
	}

	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		
		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = lambda_s;
		}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		
		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = lambda_s/2.0;
	}

	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		
		(*j)(0,0) = 0;
		(*j)(0,1) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
	}
	
	void true_solution(double t, vec* y) {
		(*y)(0) = exp(lambda_f*t);
		(*y)(1) = exp(lambda_s*t);
	}
};
#endif