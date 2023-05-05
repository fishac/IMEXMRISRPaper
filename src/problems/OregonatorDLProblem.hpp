#ifndef OREGONATORDLPROBLEM_DEFINED__
#define OREGONATORDLPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class OregonatorDLProblem : public Problem {
public:
	OregonatorDLProblem() {
		name = "OregonatorDL";
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
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = u*(-77.27*(1.675e-5*uhat + vhat) + 77.27)
		+ v*(-77.27*uhat + 77.27);

		(*f)(1) = u*(-vhat/77.27)
		+ v*(-uhat/77.27 - 1.0/77.27)
		+ w*(1.0/77.27);

		(*f)(2) = u*(0.161)
		+ w*(-0.161);
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = 77.27*v*(uhat-u) + u*(-77.27*8.375e-6*u + 77.27*1.675e-5 + 77.27*vhat);
		(*f)(1) = v*1.0/77.27*(uhat-u) + 1.0/77.27*u*vhat;
		(*f)(2) = 0.0;
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = 77.27*v*(uhat-u) + u*(-77.27*8.375e-6*u + 77.27*1.675e-5 + 77.27*vhat);
		(*f)(1) = 0.0;
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
		(*f)(1) = v*1.0/77.27*(uhat-u) + 1.0/77.27*u*vhat;
		(*f)(2) = 0.0;
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = u*(-77.27*(1.675e-5*uhat + vhat) + 77.27)
		+ v*(-77.27*uhat + 77.27);

		(*f)(1) = u*(-vhat/77.27)
		+ v*(-uhat/77.27 - 1.0/77.27)
		+ w*(1.0/77.27);

		(*f)(2) = u*(0.161)
		+ w*(-0.161);
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*f)(0) = 77.27*v*(uhat-u) + u*(-77.27*8.375e-6*u + 77.27*1.675e-5 + 77.27*vhat);
		(*f)(1) = v*1.0/77.27*(uhat-u) + 1.0/77.27*u*vhat;
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
		(*j)(1,1) = -u/77.27 - 1.0/77.26;
		(*j)(1,2) = 0.0 + 1.0/77.27;

		(*j)(2,0) = 0.161;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -0.161;
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = -77.27*(1.675e-5*uhat + vhat) + 77.27;
		(*j)(0,1) = -77.27*uhat + 77.27;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = -vhat/77.27;
		(*j)(1,1) = -uhat/77.27 - 1.0/77.26;
		(*j)(1,2) = 0.0 + 1.0/77.27;

		(*j)(2,0) = 0.161;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -0.161;
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = -77.27*(1.675e-5*(uhat-u) + (vhat-v));
		(*j)(0,1) = 77.27*(uhat-u);
		(*j)(0,2) = 0.0;

		(*j)(1,0) = (vhat-v)/77.27;
		(*j)(1,1) = (uhat-u)/77.27;
		(*j)(1,2) = 0.0;

		(*j)(2,0) = 0.0;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = 0.0;
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double w = (*y)(2);
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = -77.27*(1.675e-5*(uhat-u) + (vhat-v));
		(*j)(0,1) = 77.27*(uhat-u);
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
		double uhat = yhat(0);
		double vhat = yhat(1);
		double what = yhat(2);

		(*j)(0,0) = -77.27*(1.675e-5*uhat + vhat) + 77.27;
		(*j)(0,1) = -77.27*uhat + 77.27;
		(*j)(0,2) = 0.0;

		(*j)(1,0) = -vhat/77.27;
		(*j)(1,1) = -uhat/77.27 - 1.0/77.26;
		(*j)(1,2) = 0.0 + 1.0/77.27;

		(*j)(2,0) = 0.161;
		(*j)(2,1) = 0.0;
		(*j)(2,2) = -0.161;
	}
};

#endif