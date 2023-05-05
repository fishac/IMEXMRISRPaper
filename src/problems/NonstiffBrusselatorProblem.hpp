#ifndef NONSTIFFBRUSSELATORPROBLEM_DEFINED__
#define NONSTIFFBRUSSELATORPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class NonstiffBrusselatorProblem : public Problem {
public:
	double a = 1.0;
	double b = 3.0;
	
	NonstiffBrusselatorProblem() {
		name = "NonstiffBrusselator";
		problem_dimension = 2;
		default_H = std::pow(2.0,-6.0);
		t_0 = 0.0;
		t_f = 2.0;
		has_true_solution = false;
		explicit_only = false;
		y_0 = { 1.0, 1.0 };
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = a - u*(b+1.0) + u*u*v;
		(*f)(1) = b*u - u*u*v;
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = u*u*v;
		(*f)(1) = -u*u*v;
	}
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = a - u*(b+1.0);
		(*f)(1) = b*u;
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = -u*b;
		(*f)(1) = b*u;
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = a-u;
		(*f)(1) = 0.0;
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = -(b+1.0)+2.0*u*v;
		(*j)(0,1) = u*u;

		(*j)(1,0) = b-2.0*u*v;
		(*j)(1,1) = -u*u;
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = 2.0*u*v;
		(*j)(0,1) = u*u;
		
		(*j)(1,0) = -u*u;
		(*j)(1,1) = -u*u;
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = -b;
		(*j)(0,1) = 0.0;
		
		(*j)(1,0) = b;
		(*j)(1,1) = 0.0;
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = -b;
		(*j)(0,1) = 0.0;

		(*j)(1,0) = b;
		(*j)(1,1) = 0.0;
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
	}
};

#endif