#ifndef LORENZ63PROBLEM_DEFINED__
#define LORENZ63PROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class Lorenz63Problem : public Problem {
public:
	double sigma = 10.0;
	double beta = 8.0/3.0;
	double rho = 28.0;
	
	Lorenz63Problem() {
		name = "Lorenz63";
		problem_dimension = 3;
		default_H = std::pow(2.0,-8.0);
		t_0 = 0.0;
		t_f = 2.0;
		has_true_solution = false;
		explicit_only = false;
		y_0 = { 1.0, 1.0, 1.0 };
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = sigma*(v-u);
		(*f)(1) = u*(rho-w)-v;
		(*f)(2) = u*v-beta*w;
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 0.0;
		(*f)(1) = u*rho;
		(*f)(2) = 0.0;
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = sigma*(v-u);
		(*f)(1) = -u*w-v;
		(*f)(2) = u*v-beta*w;
	}

	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 0.0;
		(*f)(1) = -u*w-v;
		(*f)(2) = u*v;
	}

	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = sigma*(v-u);
		(*f)(1) = -v;
		(*f)(2) = -beta*w;
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
		(*f)(2) = 0.0;
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
		(*f)(2) = 0.0;
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = -sigma;
		(*j)(0,1) = sigma;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = rho-w;
		(*j)(1,1) = -1.0;
		(*j)(1,2) = -u;
		
		(*j)(2,0) = v;
		(*j)(2,1) = u;
		(*j)(2,2) = -beta;
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = rho;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;
		
		(*j)(2,0) = 0.0;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = 0.0;
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = -sigma;
		(*j)(0,1) = sigma;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = -w;
		(*j)(1,1) = -1.0;
		(*j)(1,2) = -u;
		
		(*j)(2,0) = v;
		(*j)(2,1) = u;
		(*j)(2,2) = -beta;
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = -w;
		(*j)(1,1) = -1.0;
		(*j)(1,2) = -u;
		
		(*j)(2,0) = v;
		(*j)(2,1) = u;
		(*j)(2,2) = 0.0;
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;
		
		(*j)(2,0) = 0.0;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = 0.0;
	}
	
	void true_solution(double t, vec* y) {}
};

#endif