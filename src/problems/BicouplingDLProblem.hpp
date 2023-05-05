#ifndef BICOUPLINGDLPROBLEM_DEFINED__
#define BICOUPLINGDLPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class BicouplingDLProblem : public Problem {
public:
	double a = 1.0;
	double b = 20.0;
	double w = 100.0;
	double l = 5.0;
	double p = 0.01;

	BicouplingDLProblem() {
		name = "BicouplingDL";
		problem_dimension = 3;
		default_H = std::pow(2.0,-12.0);
		t_0 = 0.0;
		t_f = 1.0;
		has_true_solution = true;
		explicit_only = false;
		y_0 = { 1.0 + a, b, a*l+b*w };
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);

		(*f)(0) = w*y1 - y2 - p*t;
		(*f)(1) = -w*y0;
		(*f)(2) = -l*y2 - l*p*t - p*std::pow(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w), 2.0) - p*std::pow(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w), 2.0);
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);

		(*f)(0) = y1*w - y2;
		(*f)(1) = -w*y0;
		(*f)(2) = y0*(-2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)))
			+ y1*(-2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
			+ y2*(-l + (2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)));
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);

		(*f)(0) = -p*t;
		(*f)(1) = 0.0;
		(*f)(2) = -l*p*t - p*std::pow(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w), 2.0) - p*std::pow(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w), 2.0)
		- (
			y0*(-2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)))
			+ y1*(-2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
			+ y2*((2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
		);
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
		(*f)(2) =  - p*std::pow(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w), 2.0) - p*std::pow(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w), 2.0)
		- (
			y0*(-2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)))
			+ y1*(-2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
			+ y2*((2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
		);
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);

		(*f)(0) = -p*t;
		(*f)(1) = 0.0;
		(*f)(2) = -l*p*t;
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);

		(*f)(0) = y1*w - y2;
		(*f)(1) = -w*y0;
		(*f)(2) = y0*(-2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)))
			+ y1*(-2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
			+ y2*(-l + (2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)));
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);

		(*f)(0) = -p*t;
		(*f)(1) = 0.0;
		(*f)(2) = -l*p*t - p*std::pow(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w), 2.0) - p*std::pow(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w), 2.0)
		- (
			y0*(-2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)))
			+ y1*(-2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
			+ y2*((2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)))
		);
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		
		(*j)(0,0) = 0.0;
		(*j)(0,1) = w;
		(*j)(0,2) = -1.0;

		(*j)(1,0) = -w;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;

		(*j)(2,0) = -2.0*p*(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w));
		(*j)(2,1) = -2.0*p*(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w));
		(*j)(2,2) = -l + (2.0*a*p)/(a*l + b*w)*(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w));
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);

		(*j)(0,0) = 0;
		(*j)(0,1) = w;
		(*j)(0,2) = -1.0;

		(*j)(1,0) = -w;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;

		(*j)(0,0) = -2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w));
		(*j)(0,1) = -2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w));
		(*j)(0,2) = -l + (2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w));
	}

	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);
		
		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;

		(*j)(2,0) = -2.0*p*(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w)) + 2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w));
		(*j)(2,1) = -2.0*p*(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w)) + 2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w));
		(*j)(2,2) = (2.0*a*p)/(a*l + b*w)*(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w))
		- ((2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)));
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);
		
		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;

		(*j)(2,0) = -2.0*p*(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w)) + 2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w));
		(*j)(2,1) = -2.0*p*(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w)) + 2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w));
		(*j)(2,2) = (2.0*a*p)/(a*l + b*w)*(y0 - a*y2/(a*l + b*w) - a*p*t/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(y1 - b*y2/(a*l + b*w) - b*p*t/(a*l + b*w))
		- ((2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w)));
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double y0 = (*y)(0);
		double y1 = (*y)(1);
		double y2 = (*y)(2);
		double yhat0 = yhat(0);
		double yhat1 = yhat(1);
		double yhat2 = yhat(2);

		(*j)(0,0) = 0;
		(*j)(0,1) = w;
		(*j)(0,2) = -1.0;

		(*j)(1,0) = -w;
		(*j)(1,1) = 0.0;
		(*j)(1,2) = 0.0;

		(*j)(0,0) = -2.0*p*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w));
		(*j)(0,1) = -2.0*p*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w));
		(*j)(0,2) = -l + (2.0*a*p)/(a*l + b*w)*(yhat0 - a*yhat2/(a*l + b*w) - a*p*that/(a*l + b*w)) + (2.0*b*p)/(a*l + b*w)*(yhat1 - b*yhat2/(a*l + b*w) - b*p*that/(a*l + b*w));
	}

	void true_solution(double t, vec* y) {
		(*y)(0) = cos(w*t) + a*exp(-l*t);
		(*y)(1) = -sin(w*t) + b*exp(-l*t);
		(*y)(2) = (a*l+b*w)*exp(-l*t) - p*t;
	}
};

#endif