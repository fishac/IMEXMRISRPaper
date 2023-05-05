#ifndef NEWTONSOLVER_DEFINED__
#define NEWTONSOLVER_DEFINED__

#include "Residual.hpp"
#include "WeightedErrorNorm.hpp"

using namespace arma;

struct NewtonSolverReturnValue {
	vec y;
	int status = 0;
};

class NewtonSolver {
public:
	Residual* residual;
	WeightedErrorNorm* err_norm;
	vec y;
	vec* y_prev;
	mat jac;
	vec f;
	vec subtraction_term;
	vec explicit_data;
	double err;
	double H;
	int max_iter;
	double tol;
	int problem_dimension;
	int status = 0;

	NewtonSolver(Residual* residual_, int max_iter_, double tol_, int problem_dimension_, WeightedErrorNorm* err_norm_) {
		residual = residual_;
		max_iter = max_iter_;
		tol = tol_;
		problem_dimension = problem_dimension_;
		err_norm = err_norm_;

		jac = mat(problem_dimension, problem_dimension, fill::zeros);
		f = vec(problem_dimension, fill::zeros);
		subtraction_term = vec(problem_dimension, fill::zeros);
		explicit_data = vec(problem_dimension, fill::zeros);
	}

	void solve(double t, vec* y_prev_, NewtonSolverReturnValue* ret) {
		status = 0;
		y_prev = y_prev_;
		y = *y_prev_;

		err = tol + 1.0;
		int iter = 0;
		residual->evaluate_explicit_data(&explicit_data);
		while(err > tol && iter < max_iter) {
			if (status == 0) {
				calculate_subtraction_term(t);
				y -= subtraction_term;
				err = err_norm->compute_norm(subtraction_term);
				iter++;
			}
		}

		if(iter == max_iter) {
			//printf("!! Newton did not converge !!\n\tError: %.16f\n", err);
			status = 1;
		} else {
			//printf("Convergence at iter: %d, error: %.16f\n", iter, err);
		}

		ret->y = y;
		ret->status = status;
	}

	void calculate_subtraction_term(double t) {
		residual->residual(t, &explicit_data, y_prev, &y, &f);
		residual->residual_jacobian(t, &y, &jac);

		if (arma::solve(subtraction_term, jac, f) == false) {
			printf("Linear solver failure.\n");
			status = 2;
		} 
	}

	void set_problem_dependent_data(double H_) {
		H = H_;
	}
};

#endif