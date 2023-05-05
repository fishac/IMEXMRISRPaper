#ifndef STRANGMARCHUKIMPLICITSOLVERESIDUAL_DEFINED__
#define STRANGMARCHUKIMPLICITSOLVERESIDUAL_DEFINED__

#include "Residual.hpp"
#include "Problem.hpp"

using namespace arma;

class StrangMarchukImplicitSolveResidual: public Residual {
public:
	Problem* problem;
	double H;
	vec f_temp;
	mat jac_temp;
	mat I;
	vec explicit_data_local;
	
	StrangMarchukImplicitSolveResidual(Problem* problem_) {
		problem = problem_;
		f_temp = vec(problem->problem_dimension,fill::zeros);
		jac_temp = mat(problem->problem_dimension,problem->problem_dimension,fill::zeros);
		I = eye(problem->problem_dimension, problem->problem_dimension);
	}

	void residual(double t, vec* explicit_data, vec* y_prev, vec* y, vec* f) {
		problem->implicit_rhs(t, y, &f_temp);
		*f = *y - *y_prev - H/4.0*(*explicit_data + f_temp);
	}
	
	void residual_jacobian(double t, vec* y, mat* jac) {
		problem->implicit_rhsjacobian(t, y, &jac_temp);
		*jac = I - H/4.0*jac_temp;
	}

	void evaluate_explicit_data(vec* explicit_data) {
		*explicit_data = explicit_data_local;
	}
	
	void set_explicit_data(vec* explicit_data_) {
		explicit_data_local = *explicit_data_;
	}
	
	void set_step_dependent_data(double H_) {
		H = H_;
	}
};

#endif