#ifndef KPRDLDLPROBLEM_DEFINED__
#define KPRDLDLPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class KPRDLProblem : public Problem {
public:
	double lambda_f = -10.0;
	double lambda_s = -1.0;
	double alpha = 1.0;
	double beta = 20.0;
	double epsilon = 0.1;
	
	KPRDLProblem() {
		name = "KPRDL";
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
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = u*(lambda_f - lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*that)/(2.0*uhat*uhat))
		+ v*((1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat));

		(*f)(1) = u*(-alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat))
		+ v*(lambda_s - lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat) + sin(that)/(2.0*vhat*vhat));
	}

	void slow_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = ( (1.0 - epsilon)*(lambda_f - lambda_s)*(-2.0 + v*v - cos(t))/(2.0*alpha*v) + lambda_f*(-3.0 + u*u - cos(beta*t))/(2.0*u) - beta*sin(beta*t)/(2.0*u)) 
		- (u*(lambda_f - lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*that)/(2.0*uhat*uhat))
		+ v*((1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat)));
		(*f)(1) = (lambda_s*(-2.0 + v*v - cos(t))/(2.0*v) - alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u) - sin(t)/(2.0*v))
		- (u*(-alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat))
		+ v*(lambda_s - lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat) + sin(that)/(2.0*vhat*vhat)));
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = ((1.0 - epsilon)*(lambda_f - lambda_s)*(-2.0 + v*v - cos(t))/(2.0*alpha*v) + lambda_f*(-3.0 + u*u - cos(beta*t))/(2.0*u) - beta*sin(beta*t)/(2.0*u)) 
		- (u*(lambda_f - lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*that)/(2.0*uhat*uhat))
		+ v*((1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat)));
		(*f)(1) = 0.0;
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = 0.0;
		(*f)(1) = (lambda_s*(-2.0 + v*v - cos(t))/(2.0*v) - alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u) - sin(t)/(2.0*v))
		- (u*(-alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat))
		+ v*(lambda_s - lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat) + sin(that)/(2.0*vhat*vhat)));
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = u*(lambda_f - lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*that)/(2.0*uhat*uhat))
		+ v*((1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat));

		(*f)(1) = u*(-alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat))
		+ v*(lambda_s - lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat) + sin(that)/(2.0*vhat*vhat));
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*f)(0) = ( (1.0 - epsilon)*(lambda_f - lambda_s)*(-2.0 + v*v - cos(t))/(2.0*alpha*v) + lambda_f*(-3.0 + u*u - cos(beta*t))/(2.0*u) - beta*sin(beta*t)/(2.0*u)) 
		- (u*(lambda_f - lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*that)/(2.0*uhat*uhat))
		+ v*((1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat)));
		(*f)(1) = (lambda_s*(-2.0 + v*v - cos(t))/(2.0*v) - alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u) - sin(t)/(2.0*v))
		- (u*(-alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat))
		+ v*(lambda_s - lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat) + sin(that)/(2.0*vhat*vhat)));
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
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = lambda_f - lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*that)/(2.0*uhat*uhat);
		(*j)(0,1) = (1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat);

		(*j)(1,0) = -alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat);
		(*j)(1,1) = lambda_s - lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat) + sin(that)/(2.0*vhat*vhat);
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = (- lambda_f*(u*u-cos(beta*t)-3.0)/(2.0*u*u)) + (lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*t)/(2.0*u*u) - beta*sin(beta*that)/(2.0*uhat*uhat));
		(*j)(0,1) = (- (1.0-epsilon)*(lambda_f - lambda_s)*(v*v-cos(t)-2.0)/(2.0*alpha*v*v)) + ((1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat));

		(*j)(1,0) = (alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + u*u - cos(beta*t))/(2.0*u*u)) - (alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat));
		(*j)(1,1) = (lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat)) - (lambda_s*(-2.0 + v*v - cos(t))/(2.0*v*v)) + sin(t)/(2.0*v*v) - sin(that)/(2.0*vhat*vhat);
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = (- lambda_f*(u*u-cos(beta*t)-3.0)/(2.0*u*u)) + (lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*t)/(2.0*u*u) - beta*sin(beta*that)/(2.0*uhat*uhat));
		(*j)(0,1) = (- (1.0-epsilon)*(lambda_f - lambda_s)*(v*v-cos(t)-2.0)/(2.0*alpha*v*v)) + ((1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat));

		(*j)(1,0) = 0.0;
		(*j)(1,1) = 0.0;
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		double u = (*y)(0);
		double v = (*y)(1);
		double uhat = yhat(0);
		double vhat = yhat(1);

		(*j)(0,0) = lambda_f - lambda_f*(uhat*uhat-cos(beta*that)-3.0)/(2.0*uhat*uhat) + beta*sin(beta*that)/(2.0*uhat*uhat);
		(*j)(0,1) = (1.0-epsilon)*(lambda_f - lambda_s)/alpha - (1.0-epsilon)*(lambda_f - lambda_s)*(vhat*vhat-cos(that)-2.0)/(2.0*alpha*vhat*vhat);

		(*j)(1,0) = -alpha*epsilon*(lambda_f - lambda_s) + alpha*epsilon*(lambda_f - lambda_s)*(-3.0 + uhat*uhat - cos(beta*that))/(2.0*uhat*uhat);
		(*j)(1,1) = lambda_s - lambda_s*(-2.0 + vhat*vhat - cos(that))/(2.0*vhat*vhat) + sin(that)/(2.0*vhat*vhat);
	}
	
	void true_solution(double t, vec* y) {
		(*y)(0) = sqrt(3.0 + cos(beta*t));
		(*y)(1) = sqrt(2.0 + cos(t));
	}
};

#endif