#ifndef LIETROTTERMETHOD_DEFINED__
#define LIETROTTERMETHOD_DEFINED__

#include "LieTrotterExplicitSolveProblem.hpp"
#include "LieTrotterImplicitSolveResidual.hpp"
#include "NewtonSolver.hpp"
#include "FixedDIRKMethod.hpp"
#include "WeightedErrorNorm.hpp"

using namespace arma;

class LieTrotterMethod {
public:
	vec y;
	vec v_0;
	vec v_H;
	vec f_temp;
	mat Y;
	double effective_H;
	int total_output_points;
	int problem_dimension;
	Problem* problem;
	LieTrotterExplicitSolveProblem explicit_solve_problem;
	LieTrotterImplicitSolveResidual implicit_solve_residual;
	NewtonSolver newton_solver;
	struct NewtonSolverReturnValue newton_ret;
	FixedDIRKMethod dirk;
	WeightedErrorNorm* err_norm;

	LieTrotterMethod(SingleRateMethodCoefficients* inner_coeffs_, Problem* problem_, int problem_dimension_, WeightedErrorNorm* err_norm_) : 
	implicit_solve_residual(problem_),
	explicit_solve_problem(problem_),
	dirk(inner_coeffs_, &(LieTrotterMethod::explicit_solve_problem), problem_dimension_, err_norm_, false),
	newton_solver(&(LieTrotterMethod::implicit_solve_residual), 20, 1.0, problem_dimension_, err_norm_) 
	{
		problem = problem_;
		problem_dimension = problem_dimension_;
		err_norm = err_norm_;
	}

	mat solve(double t_0, double H, int M, vec* y_0, vec* output_tspan) {
		prepare_solve(H, output_tspan);
		y = *y_0;
		problem->set_yhat_that(*y_0,t_0);
		err_norm->set_weights(y_0);

		int output_index = 0;
		if((*output_tspan)(0) == t_0) {
			Y.col(0) = *y_0;
			output_index++;
		}

		double t = t_0;
		while(output_index < total_output_points) {
			if(t + H - (*output_tspan)(output_index) >= 1e-14) {
				effective_H = (*output_tspan)(output_index) - t;
			} else {
				effective_H = H;
			}
			step_solution(t, effective_H, M, &y, &y);
			if (abs(t + effective_H - (*output_tspan)(output_index)) < 1e-14) {
				Y.col(output_index) = y;
				output_index++;
			}
			t += effective_H;
			problem->set_yhat_that(y,t);
			err_norm->set_weights(&y);
		}
		return Y;
	}
	
	void step_solution(double t, double H, int M, vec* y_n, vec* y_np1) {
		implicit_solve_residual.set_step_dependent_data(H);
		explicit_solve_problem.set_step_dependent_data(t);
		
		v_0 = *y_n;
		problem->explicit_rhs(t,&v_0,&f_temp);
		v_H = v_0 + H*f_temp;
		
		v_0 = v_H;
		newton_solver.solve(t+H, &v_0, &newton_ret);
		problem->increment_slow_nonlinear_solves();
		v_H = newton_ret.y;
		
		v_0 = v_H;
		vec output_tspan = {H};
		mat Y = dirk.solve(0.0, H/M, &v_0, &output_tspan);
		v_H = Y.col(0);
		
		*y_np1 = v_H;
	}

	void prepare_solve(double H_, vec* output_tspan) {
		effective_H = H_;
		total_output_points = output_tspan->n_elem;
		y = vec(problem_dimension, fill::zeros);
		v_0 = vec(problem_dimension, fill::zeros);
		v_H = vec(problem_dimension, fill::zeros);
		f_temp = vec(problem_dimension, fill::zeros);
		Y = mat(problem_dimension, total_output_points, fill::zeros);
	}
};

#endif