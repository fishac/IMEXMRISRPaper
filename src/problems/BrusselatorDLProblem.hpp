#ifndef BRUSSELATORDLPROBLEM_DEFINED__
#define BRUSSELATORDLPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class BrusselatorDLProblem : public Problem {
public:
	double a = 1.0;
	double b = 3.5;
	double epsilon = 0.01;
	
	BrusselatorDLProblem() {
		name = "BrusselatorDL";
		problem_dimension = 3;
		default_H = std::pow(2.0,-6.0);
		t_0 = 0.0;
		t_f = 2.0;
		has_true_solution = false;
		explicit_only = false;
		y_0 = { 1.2, 3.1, 3.0 };
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*f)(0) = a - (w + 1.0)*u + u*u*v;
		(*f)(1) = w*u - u*u*v;
		(*f)(2) = (b-w)/epsilon - u*w;
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = u*(-(what + 1.0) + 2.0*uhat*vhat)
		+ v*(uhat*uhat)
		+ w*(-uhat);

		(*f)(1) = u*(what - 2.0*uhat*vhat)
		+ v*(-uhat*uhat)
		+ w*(uhat);

		(*f)(2) = u*(-what)
		+ w*(-1.0/epsilon - uhat);
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = a + (what - w)*u + v*(u*u-uhat*uhat) + uhat*(w - u*2.0*vhat);
		(*f)(1) = (w-what)*u + u*2.0*uhat*vhat + v*(uhat*uhat-u*u) - w*(uhat);
		(*f)(2) = b/epsilon + u*(what-w) + w*uhat;
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = a + (what - w)*u + v*(u*u-uhat*uhat) + uhat*(w - u*2.0*vhat);
		(*f)(1) = (w-what)*u + u*2.0*uhat*vhat + v*(uhat*uhat-u*u) - w*(uhat);
		(*f)(2) = 0.0;
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
		(*f)(2) = b/epsilon + u*(what-w) + w*uhat;
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = u*(-(what + 1.0) + 2.0*uhat*vhat)
		+ v*(uhat*uhat)
		+ w*(-uhat);

		(*f)(1) = u*(what - 2.0*uhat*vhat)
		+ v*(-uhat*uhat)
		+ w*(uhat);

		(*f)(2) = u*(-what)
		+ w*(-1.0/epsilon - uhat);
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = a + (what - w)*u + v*(u*u-uhat*uhat) + uhat*(w - u*2.0*vhat);
		(*f)(1) = (w-what)*u + u*2.0*uhat*vhat + v*(uhat*uhat-u*u) - w*(uhat);
		(*f)(2) = b/epsilon + u*(what-w) + w*uhat;
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);

		(*j)(0,0) = -(w + 1.0) + 2.0*u*v;
		(*j)(0,1) = u*u;
		(*j)(0,2) = -u;

		(*j)(1,0) = w - 2.0*u*v;
		(*j)(1,1) = -u*u;
		(*j)(1,2) = u;

		(*j)(2,0) = -w;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -1.0/epsilon -u;
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = -(what + 1.0) + 2.0*uhat*vhat;
		(*j)(0,1) = uhat*uhat;
		(*j)(0,2) = -uhat;

		(*j)(1,0) = what - 2.0*uhat*vhat;
		(*j)(1,1) = -uhat*uhat;
		(*j)(1,2) = uhat;

		(*j)(2,0) = -what;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -1.0/epsilon - uhat;
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = (what - w) + 2.0*u*v - 2.0*uhat*vhat;
		(*j)(0,1) = u*u-uhat*uhat;
		(*j)(0,2) = uhat-u;
		
		(*j)(1,0) = w-what + 2.0*(uhat*vhat-u*v);
		(*j)(1,1) = uhat*uhat-u*u;
		(*j)(1,2) = uhat-u;
		
		(*j)(2,0) = what-w;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = uhat-u;
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = (what - w) + 2.0*u*v - 2.0*uhat*vhat;
		(*j)(0,1) = u*u-uhat*uhat;
		(*j)(0,2) = uhat-u;
		
		(*j)(1,0) = w-what + 2.0*(uhat*vhat-u*v);
		(*j)(1,1) = uhat*uhat-u*u;
		(*j)(1,2) = uhat-u;

		(*j)(2,0) = 0.0;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = 0.0;
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = -(what + 1.0) + 2.0*uhat*vhat;
		(*j)(0,1) = uhat*uhat;
		(*j)(0,2) = -uhat;

		(*j)(1,0) = what - 2.0*uhat*vhat;
		(*j)(1,1) = -uhat*uhat;
		(*j)(1,2) = uhat;

		(*j)(2,0) = -what;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -1.0/epsilon - uhat;
	}
};

#endif