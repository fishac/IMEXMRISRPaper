#ifndef BRUSSELATORPDEPROBLEM_DEFINED__
#define BRUSSELATORPDEPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class BrusselatorPDEProblem : public Problem {
public:
	int n = 0;
	double dx = 0.0;
	int u_start = 0;
	int v_start = 0;
	int w_start = 0;
	int u_end = 0;
	int v_end = 0;
	int w_end = 0;
	vec x;
	double a = 0.6;
	double b = 2.0;
	double epsilon = 0.01;
	double alpha_u = 0.01;
	double alpha_v = 0.01;
	double alpha_w = 0.01;
	double rho_u = 0.001;
	double rho_v = 0.001;
	double rho_w = 0.001;
	
	BrusselatorPDEProblem(int n_) {
		name = "BrusselatorPDE";
		n = n_;
		time_step_multiplier = 0.1;
		problem_dimension = 3*n;
		u_start = 0;
		v_start = n;
		w_start = 2*n;
		u_end = n-1;
		v_end = 2*n-1;
		w_end = 3*n-1;
		default_H = std::pow(2.0,-9.0);
		t_0 = 0.0;
		t_f = 10.0;
		has_true_solution = false;
		explicit_only = false;
		x = linspace(0.0, 1.0, n_);
		dx = x(1)-x(0);
		y_0 = get_y_0();
	}
	
	vec get_y_0() {
		vec y_0_(problem_dimension,fill::zeros);
		for(int i=0; i<n; i++) {
			y_0_(u_start + i) = a + 0.1*sin(M_PI*x(i));
			y_0_(v_start + i) = b/a + 0.1*sin(M_PI*x(i));
			y_0_(w_start + i) = b + 0.1*sin(M_PI*x(i));
		}
		return y_0_;
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
		
		for(int i=1; i<n-1; i++) {
			double un = (*y)(u_start+i);
			double unm1 = (*y)(u_start+i-1);
			double unp1 = (*y)(u_start+i+1);

			double vn = (*y)(v_start+i);
			double vnm1 = (*y)(v_start+i-1);
			double vnp1 = (*y)(v_start+i+1);
			
			double wn = (*y)(w_start+i);
			double wnm1 = (*y)(w_start+i-1);
			double wnp1 = (*y)(w_start+i+1);
			
			(*f)(u_start+i) = alpha_u*(1.0*unm1 -2.0*un + 1.0*unp1)/(dx*dx) - rho_u*(-1.0*unm1 + 1.0*unp1)/(2.0*dx) + (a - (wn+1.0)*un + un*un*vn);
			(*f)(v_start+i) = alpha_v*(1.0*vnm1 -2.0*vn + 1.0*vnp1)/(dx*dx) - rho_v*(-1.0*vnm1 + 1.0*vnp1)/(2.0*dx) + (un*wn - un*un*vn);
			(*f)(w_start+i) = alpha_w*(1.0*wnm1 -2.0*wn + 1.0*wnp1)/(dx*dx) - rho_w*(-1.0*wnm1 + 1.0*wnp1)/(2.0*dx) + ((b - wn)/epsilon - un*wn);
		}
	}
	
	void slow_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
		
		for(int i=1; i<n-1; i++) {
			double un = (*y)(u_start+i);
			double unm1 = (*y)(u_start+i-1);
			double unp1 = (*y)(u_start+i+1);

			double vn = (*y)(v_start+i);
			double vnm1 = (*y)(v_start+i-1);
			double vnp1 = (*y)(v_start+i+1);
			
			double wn = (*y)(w_start+i);
			double wnm1 = (*y)(w_start+i-1);
			double wnp1 = (*y)(w_start+i+1);
			
			(*f)(u_start+i) = alpha_u*(1.0*unm1 -2.0*un + 1.0*unp1)/(dx*dx) - rho_u*(-1.0*unm1 + 1.0*unp1)/(2.0*dx);
			(*f)(v_start+i) = alpha_v*(1.0*vnm1 -2.0*vn + 1.0*vnp1)/(dx*dx) - rho_v*(-1.0*vnm1 + 1.0*vnp1)/(2.0*dx);
			(*f)(w_start+i) = alpha_w*(1.0*wnm1 -2.0*wn + 1.0*wnp1)/(dx*dx) - rho_w*(-1.0*wnm1 + 1.0*wnp1)/(2.0*dx);
		}
	}

	void fast_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
		
		for(int i=1; i<n-1; i++) {
			double un = (*y)(u_start+i);
			double unm1 = (*y)(u_start+i-1);
			double unp1 = (*y)(u_start+i+1);

			double vn = (*y)(v_start+i);
			double vnm1 = (*y)(v_start+i-1);
			double vnp1 = (*y)(v_start+i+1);
			
			double wn = (*y)(w_start+i);
			double wnm1 = (*y)(w_start+i-1);
			double wnp1 = (*y)(w_start+i+1);
			
			(*f)(u_start+i) = a - (wn+1.0)*un + un*un*vn;
			(*f)(v_start+i) = un*wn - un*un*vn;
			(*f)(w_start+i) = (b - wn)/epsilon - un*wn;
		}
	}
	
	void implicit_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
		
		for(int i=1; i<n-1; i++) {
			double un = (*y)(u_start+i);
			double unm1 = (*y)(u_start+i-1);
			double unp1 = (*y)(u_start+i+1);

			double vn = (*y)(v_start+i);
			double vnm1 = (*y)(v_start+i-1);
			double vnp1 = (*y)(v_start+i+1);
			
			double wn = (*y)(w_start+i);
			double wnm1 = (*y)(w_start+i-1);
			double wnp1 = (*y)(w_start+i+1);
			
			(*f)(u_start+i) = alpha_u*(1.0*unm1 -2.0*un + 1.0*unp1)/(dx*dx);
			(*f)(v_start+i) = alpha_v*(1.0*vnm1 -2.0*vn + 1.0*vnp1)/(dx*dx);
			(*f)(w_start+i) = alpha_w*(1.0*wnm1 -2.0*wn + 1.0*wnp1)/(dx*dx);
		}
	}
	
	void explicit_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
		
		for(int i=1; i<n-1; i++) {
			double un = (*y)(u_start+i);
			double unm1 = (*y)(u_start+i-1);
			double unp1 = (*y)(u_start+i+1);

			double vn = (*y)(v_start+i);
			double vnm1 = (*y)(v_start+i-1);
			double vnp1 = (*y)(v_start+i+1);
			
			double wn = (*y)(w_start+i);
			double wnm1 = (*y)(w_start+i-1);
			double wnp1 = (*y)(w_start+i+1);
			
			(*f)(u_start+i) = -rho_u*(-1.0*unm1 + 1.0*unp1)/(2.0*dx);
			(*f)(v_start+i) = -rho_v*(-1.0*vnm1 + 1.0*vnp1)/(2.0*dx);
			(*f)(w_start+i) = -rho_w*(-1.0*wnm1 + 1.0*wnp1)/(2.0*dx);
		}
	}
	
	void linear_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
	}
	
	void full_rhsjacobian_custom(double t, vec* y, mat* j) {
		j->zeros();
		
		for(int i=1; i<n-1; i++) {
			(*j)(u_start+i,u_start+i-1) = alpha_u/(dx*dx)+rho_u/(2.0*dx);
			(*j)(u_start+i,u_start+i) = alpha_u*(-2.0)/(dx*dx) + ((*y)(w_start+i)+1) + 2.0*(*y)(u_start+i)*(*y)(v_start+i);
			(*j)(u_start+i,u_start+i+1) = alpha_u/(dx*dx)-rho_w/(2.0*dx);
			(*j)(u_start+i,v_start+i) = (*y)(u_start+i)*(*y)(u_start+i);
			(*j)(u_start+i,w_start+i) = (*y)(u_start+i);
			
			(*j)(v_start+i,v_start+i-1) = alpha_v/(dx*dx)+rho_v/(2.0*dx);
			(*j)(v_start+i,v_start+i) = alpha_v*(-2.0)/(dx*dx) + -(*y)(u_start+i)*(*y)(u_start+i);
			(*j)(v_start+i,v_start+i+1) = alpha_v/(dx*dx)-rho_v/(2.0*dx);
			(*j)(v_start+i,u_start+i) = ((*y)(w_start+i) - 2.0*(*y)(u_start+i)*(*y)(v_start+i));
			(*j)(v_start+i,w_start+i) = (*y)(u_start+i);
			
			(*j)(w_start+i,w_start+i-1) = alpha_w/(dx*dx)+rho_w/(2.0*dx);
			(*j)(w_start+i,w_start+i) = alpha_w*(-2.0)/(dx*dx) + -1.0/epsilon - (*y)(u_start+i);
			(*j)(w_start+i,w_start+i+1) = alpha_w/(dx*dx)-rho_w/(2.0*dx);
			(*j)(w_start+i,u_start+i) = -(*y)(w_start+i);
			(*j)(w_start+i,v_start+i) = 0.0;
		}
	}
	
	void fast_rhsjacobian_custom(double t, vec* y, mat* j) {
		j->zeros();
		
		for(int i=1; i<n-1; i++) {
			(*j)(u_start+i,u_start+i) = ((*y)(w_start+i)+1) + 2.0*(*y)(u_start+i)*(*y)(v_start+i);
			(*j)(u_start+i,v_start+i) = (*y)(u_start+i)*(*y)(u_start+i);
			(*j)(u_start+i,w_start+i) = (*y)(u_start+i);
			
			(*j)(v_start+i,v_start+i) = -(*y)(u_start+i)*(*y)(u_start+i);
			(*j)(v_start+i,u_start+i) = (*y)(w_start+i) - 2.0*(*y)(u_start+i)*(*y)(v_start+i);
			(*j)(v_start+i,w_start+i) = (*y)(u_start+i);
			
			(*j)(w_start+i,w_start+i) = -1.0/epsilon - (*y)(u_start+i);
			(*j)(w_start+i,u_start+i) = -(*y)(w_start+i);
			(*j)(w_start+i,v_start+i) = 0.0;
		}
	}
	
	void slow_rhsjacobian_custom(double t, vec* y, mat* j) {
		j->zeros();
		
		for(int i=1; i<n-1; i++) {
			(*j)(u_start+i,u_start+i-1) = alpha_u/(dx*dx)+rho_u/(2.0*dx);
			(*j)(u_start+i,u_start+i) = alpha_u*(-2.0)/(dx*dx);
			(*j)(u_start+i,u_start+i+1) = alpha_u/(dx*dx)-rho_u/(2.0*dx);
			
			(*j)(v_start+i,v_start+i-1) = alpha_v/(dx*dx)+rho_v/(2.0*dx);
			(*j)(v_start+i,v_start+i) = alpha_v*(-2.0)/(dx*dx);
			(*j)(v_start+i,v_start+i+1) = alpha_v/(dx*dx)-rho_v/(2.0*dx);
			
			(*j)(w_start+i,w_start+i-1) = alpha_w/(dx*dx)+rho_w/(2.0*dx);
			(*j)(w_start+i,w_start+i) = alpha_w*(-2.0)/(dx*dx);
			(*j)(w_start+i,w_start+i+1) = alpha_w/(dx*dx)-rho_w/(2.0*dx);
		}
	}
	
	void implicit_rhsjacobian_custom(double t, vec* y, mat* j) {
		j->zeros();
		
		for(int i=1; i<n-1; i++) {
			(*j)(u_start+i,u_start+i-1) = alpha_u/(dx*dx);
			(*j)(u_start+i,u_start+i) = alpha_u*(-2.0)/(dx*dx);
			(*j)(u_start+i,u_start+i+1) = alpha_u/(dx*dx);
			
			(*j)(v_start+i,v_start+i-1) = alpha_v/(dx*dx);
			(*j)(v_start+i,v_start+i) = alpha_v*(-2.0)/(dx*dx);
			(*j)(v_start+i,v_start+i+1) = alpha_v/(dx*dx);
			
			(*j)(w_start+i,w_start+i-1) = alpha_w/(dx*dx);
			(*j)(w_start+i,w_start+i) = alpha_w*(-2.0)/(dx*dx);
			(*j)(w_start+i,w_start+i+1) = alpha_w/(dx*dx);
		}
	}
	
	void linear_rhsjacobian_custom(double t, vec* y, mat* j) {
		j->zeros();
	}
};

#endif