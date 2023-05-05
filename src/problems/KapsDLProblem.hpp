#ifndef KAPSDLPROBLEM_DEFINED__
#define KAPSDLPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class KapsDLProblem : public Problem {
public:
	KapsDLProblem() {
		name = "KapsDL";
		problem_dimension = 2;
		default_H = std::pow(2.0,-6.0);
		t_0 = 0.0;
		t_f = 2.0;
		has_true_solution = true;
		explicit_only = false;
		y_0 = { 1.0, 1.0 };
	}

	void full_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = -102.0*u+100.0*v*v;
		(*f)(1) = -v*v + u - v;
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = u*(-102.0)
		+ v*(100.0*2.0*vhat);

		(*f)(1) = u*1.0
		+ v*(-2.0*vhat - 1.0);
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = 100.0*v*(v-2.0*vhat);
		(*f)(1) = v*(2.0*vhat-v);
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = 100.0*v*(v-2.0*vhat);
		(*f)(1) = 0.0;
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = 0.0;
		(*f)(1) = v*(2.0*vhat-v);
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = u*(-102.0)
		+ v*(100.0*2.0*vhat);

		(*f)(1) = u*1.0
		+ v*(-2.0*vhat - 1.0);
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = 100.0*v*(v-2.0*vhat);
		(*f)(1) = v*(2.0*vhat-v);
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		
		(*j)(0,0) = -102.0;
		(*j)(0,1) = 100.0*2.0*v;

		(*j)(1,0) = 1.0;
		(*j)(1,1) = -2.0*v - 1.0;
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = -102.0;
		(*j)(0,1) = 100.0*2.0*vhat;

		(*j)(1,0) = 1.0;
		(*j)(1,1) = -2.0*vhat - 1.0;
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 100.0*2.0*(v-vhat);
		
		(*j)(1,0) = 0.0;
		(*j)(1,1) = 2.0*(vhat-v);
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 100.0*2.0*(v-vhat);

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = -102.0;
		(*j)(0,1) = 100.0*2.0*vhat;

		(*j)(1,0) = 1.0;
		(*j)(1,1) = -2.0*vhat - 1.0;
	}
	
	void true_solution(double t, vec* y) {
		(*y)(0) = exp(-2.0*t);
		(*y)(1) = exp(-t);
	}
};

#endif