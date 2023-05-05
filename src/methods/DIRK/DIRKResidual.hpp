#ifndef DIRKRESIDUAL_DEFINED__
#define DIRKRESIDUAL_DEFINED__

#include "SingleRateMethodCoefficients.hpp"
#include "Residual.hpp"
#include "Problem.hpp"

using namespace arma;

class DIRKResidual: public Residual {
public:
	Problem* problem;
	SingleRateMethodCoefficients* coeffs;
	vec y_temp;
	mat jac_temp;
	mat I;
	int problem_dimension;
	vec* explicit_data;
	double h;
	int stage_index;

	DIRKResidual(SingleRateMethodCoefficients* coeffs_, Problem* problem_, int problem_dimension_) {
		problem = problem_;
		coeffs = coeffs_;
		problem_dimension = problem_dimension_;
		y_temp = vec(problem_dimension, fill::zeros);
		jac_temp = mat(problem_dimension, problem_dimension, fill::zeros);
		I = eye(problem_dimension,problem_dimension);
	}

	void residual(double t, vec* explicit_data, vec* y_0, vec* y, vec* f) {
		problem->full_rhs(t+(coeffs->c(stage_index))*h, y, &y_temp);
		*f = *y - *y_0 - h*(*explicit_data + (coeffs->A(stage_index,stage_index))*y_temp);
	}
	
	void residual_jacobian(double t, vec* y, mat* jac) {
		problem->full_rhsjacobian(t+(coeffs->c(stage_index))*h, y, &jac_temp);
		*jac = I - h*(coeffs->A(stage_index,stage_index))*jac_temp;
	}

	void evaluate_explicit_data(vec* explicit_data_) {
		*explicit_data_ = *explicit_data;
	}

	void set_problem_dependent_data(double h_) {
		h = h_;
	}

	void set_explicit_data_pointer(vec* explicit_data_) {
		explicit_data = explicit_data_;
	}

	void set_function_dependent_data(int stage_index_) {
		stage_index = stage_index_;
	}

};

#endif