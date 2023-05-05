#ifndef OREGONATORPROBLEM_DEFINED__
#define OREGONATORPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class OregonatorProblem : public Problem {
public:
	OregonatorProblem() {
		name = "Oregonator";
		problem_dimension = 3;
		default_H = std::pow(2.0,-12.0);
		t_0 = 0.0;
		t_f = 360.0;
		has_true_solution = false;
		explicit_only = false;
		y_0 = { 1.0, 2.0, 3.0 };
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = -77.27*(8.375e-6*u*u + u*v) + 77.27*(u + v);
		(*f)(1) = -u*v/77.27 + (w-v)/77.27;
		(*f)(2) = 0.161*(u-w);
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = -77.27*(8.375e-6*u*u + u*v);
		(*f)(1) = -u*v/77.27;
		(*f)(2) = 0.0;
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 77.27*(u + v);
		(*f)(1) = (w-v)/77.27;
		(*f)(2) = 0.161*(u-w);
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 77.27*(u + v);
		(*f)(1) = 0.0;
		(*f)(2) = 0.0;
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 0.0;
		(*f)(1) = (w-v)/77.27;
		(*f)(2) = 0.161*(u-w);
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = 77.27*(u + v);
		(*f)(1) = (w-v)/77.27;
		(*f)(2) = 0.161*(u-w);
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = -77.27*(8.375e-6*u*u + u*v);
		(*f)(1) = -u*v/77.27;
		(*f)(2) = 0.0;
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = -77.27*(1.675e-5*u + v) + 77.27;
		(*j)(0,1) = -77.27*u + 77.27;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = -v/77.27;
		(*j)(1,1) = -u/77.27 - 1.0/77.27;
		(*j)(1,2) = 0.0 + 1.0/77.27;

		(*j)(2,0) = 0.161;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -0.161;
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = -77.27*(1.675e-5*u + v);
		(*j)(0,1) = -77.27*u;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = -v/77.27;
		(*j)(1,1) = -u/77.27;
		(*j)(1,2) = 0.0;

		(*j)(2,0) = 0.0;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = 0.0;
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = 77.27;
		(*j)(0,1) = 77.27;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = -1.0/77.27;
		(*j)(1,2) = 1.0/77.27;

		(*j)(2,0) = 0.161;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -0.161;
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = 77.27;
		(*j)(0,1) = 77.27;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;

		(*j)(2,0) = 0.0;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = 0.0;
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = 77.27;
		(*j)(0,1) = 77.27;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = -1.0/77.27;
		(*j)(1,2) = 1.0/77.27;

		(*j)(2,0) = 0.161;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -0.161;
	}
};

#endif