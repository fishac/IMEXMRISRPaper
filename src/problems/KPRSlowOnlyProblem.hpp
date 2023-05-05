#ifndef KPRSLOWONLYPROBLEM_DEFINED__
#define KPRSLOWONLYPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class KPRSlowOnlyProblem : public Problem {
public:
	double lambda_f = -10.0;
	double lambda_s = -1.0;
	double alpha = 1.0;
	double beta = 20.0;
	double epsilon = 0.1;
	
	KPRSlowOnlyProblem() {
		name = "KPRSlowOnly";
		problem_dimension = 2;
		default_H = M_PI * std::pow(2.0,-6.0);
		t_0 = 0.0;
		t_f = 5.0*M_PI/2.0;
		has_true_solution = true;
		explicit_only = false;
		y_0 = { 2.0, sqrt(3.0) };
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = (1.0 - epsilon)*(lambda_f - lambda_s)*(-2.0 + v*v - cos(t))/(2.0*alpha*v) + lambda_f*(-3.0 + u*u - cos(beta*t))/(2.0*u) - beta*sin(beta*t)/(2.0*u);
		(*f)(1) = lambda_s*(-2.0 + v*v - cos(t))/(2.0*v) - alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u) - sin(t)/(2.0*v);
	}
		
	void fast_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = (1.0 - epsilon)*(lambda_f - lambda_s)*(-2.0 + v*v - cos(t))/(2.0*alpha*v) + lambda_f*(-3.0 + u*u - cos(beta*t))/(2.0*u) - beta*sin(beta*t)/(2.0*u);
		(*f)(1) = lambda_s*(-2.0 + v*v - cos(t))/(2.0*v) - alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u) - sin(t)/(2.0*v);
	}


	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
	}

	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double v = (*y)(1);

		(*f)(0) = 0.0;
		(*f)(1) = 0.0;
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = (1.0 - epsilon)*(lambda_f - lambda_s)/(2.0*alpha)*v + lambda_f/2.0*u;
		(*f)(1) = lambda_s/2.0*v - alpha*epsilon*(lambda_f - lambda_s)/2.0*u;
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*f)(0) = (1.0 - epsilon)*(lambda_f - lambda_s)*(-2.0 - cos(t))/(2.0*alpha*v) + lambda_f*(-3.0 - cos(beta*t))/(2.0*u) - beta*sin(beta*t)/(2.0*u);
		(*f)(1) = lambda_s*(-2.0 - cos(t))/(2.0*v) - alpha*epsilon*(lambda_f - lambda_s)*(-3.0 - cos(beta*t))/(2.0*u) - sin(t)/(2.0*v);
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = lambda_f - lambda_f*(u*u-cos(beta*t)-3.0)/(2.0*u*u) + beta*sin(beta*t)/(2.0*u*u);
		(*j)(0,1) = (1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(v*v-cos(t)-2.0)/(2.0*alpha*v*v);

		(*j)(1,0) = -alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u*u);
		(*j)(1,1) = lambda_s - lambda_s*(-2.0 + v*v - cos(t))/(2.0*v*v) + sin(t)/(2.0*v*v);
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = lambda_f - lambda_f*(u*u-cos(beta*t)-3.0)/(2.0*u*u) + beta*sin(beta*t)/(2.0*u*u);
		(*j)(0,1) = (1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(v*v-cos(t)-2.0)/(2.0*alpha*v*v);
		
		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;
		
		(*j)(1,0) = -alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u*u);
		(*j)(1,1) = lambda_s - lambda_s*(-2.0 + v*v - cos(t))/(2.0*v*v) + sin(t)/(2.0*v*v);
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = 0.0;
		(*j)(0,1) = 0.0;
		
		(*j)(1,0) = -alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u*u);
		(*j)(1,1) = lambda_s - lambda_s*(-2.0 + v*v - cos(t))/(2.0*v*v);
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);

		(*j)(0,0) = lambda_f/2.0;
		(*j)(0,1) = (1.0 - epsilon)*(lambda_f - lambda_s)/(2.0*alpha);

		(*j)(1,0) = -alpha*epsilon*(lambda_f - lambda_s)/2.0;
		(*j)(1,1) = lambda_s/2.0;
	}
	
	void true_solution(double t, vec* y) {
		(*y)(0) = sqrt(3.0 + cos(beta*t));
		(*y)(1) = sqrt(2.0 + cos(t));
	}
};

#endif