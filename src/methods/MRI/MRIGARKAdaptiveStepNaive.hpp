#ifndef MRIGARKADAPTIVESTEPNAIVE_DEFINED__
#define MRIGARKADAPTIVESTEPNAIVE_DEFINED__

#include "MRIGARKAdaptiveStep.hpp"
#include "MRIGARKCoefficients.hpp"
#include "MRIGARKInnerRHSFunctions.hpp"
#include "MRIGARKExplicitSolveProblem.hpp"
#include "MRIGARKImplicitSolveResidual.hpp"
#include "AdaptiveDIRKMethod.hpp"
#include "SingleRateMethodCoefficients.hpp"
#include "AdaptiveDIRKMethod.hpp"
#include "NewtonSolver.hpp"
#include "Controller.hpp"
#include "WeightedErrorNorm.hpp"
#include <set>

using namespace arma;

class MRIGARKAdaptiveStepNaive : public MRIGARKAdaptiveStep {
public:
	MRIGARKCoefficients* coeffs;
	MRIGARKInnerRHSFunctions inner_rhs_funcs;
	MRIGARKExplicitSolveProblem explicit_solve_problem;
	NewtonSolver newton_solver;
	struct NewtonSolverReturnValue newton_ret;
	AdaptiveDIRKMethod dirk;
	struct AdaptiveSingleRateMethodReturnValue dirk_ret;
	Controller* controller;
	WeightedErrorNorm* err_norm;
	vec output_tspan;
	mat y_stages;
	vec y;
	vec y_hat;
	vec v_0;
	vec v_H;
	double H;
	int M;
	int total_microtimesteps;
	int total_successful_microtimesteps;
	int problem_dimension;
	int num_stages;
	int status;
	int inner_status;

	MRIGARKAdaptiveStepNaive(Problem* problem, MRIGARKCoefficients* coeffs_, SingleRateMethodCoefficients* inner_coeffs_, Controller* controller_, RHS* fast_func_, RHS* slow_func_, RHSJacobian* fast_func_jac_, RHSJacobian* slow_func_jac_, int problem_dimension_, WeightedErrorNorm* err_norm_) :
	inner_rhs_funcs(coeffs_, fast_func_, slow_func_, fast_func_jac_, slow_func_jac_, problem_dimension_),
	implicit_solve_residual(&(MRIGARKAdaptiveStepNaive::inner_rhs_funcs)),
	explicit_solve_problem(&(MRIGARKAdaptiveStepNaive::inner_rhs_funcs)),
	dirk(inner_coeffs_, &(MRIGARKAdaptiveStepNaive::explicit_solve_problem), problem_dimension_, err_norm_, false),
	newton_solver(&(MRIGARKAdaptiveStepNaive::implicit_solve_residual), 20, 0.1, problem_dimension_, err_norm_)
	{
		coeffs = coeffs_;
		num_stages = coeffs->num_stages;
		problem_dimension = problem_dimension_;
		controller = controller_;
		err_norm = err_norm_;
		declare_vectors();
	}

	// Adaptive step solution
	void step_solution(double t, double H, int M, vec* y_0, int esf_measurement_type, MRIGARKAdaptiveStepReturnValue* ret) {
		prepare_solve(H);
		total_microtimesteps = 0;
		total_successful_microtimesteps = 0;
		status = 0;

		y_stages.col(0) = *y_0;
		for(int stage_index=1; stage_index<num_stages; stage_index++) {
			if (status == 0) {
				inner_rhs_funcs.set_function_dependent_data(H, t, stage_index, false);
				v_0 = y_stages.col(stage_index-1);

				if (coeffs->method_type == 0) {
					double gbar = inner_rhs_funcs.gamma_bar(stage_index,stage_index);
					double delta_c = coeffs->c(stage_index) - coeffs->c(stage_index-1);
					if (delta_c == 0) {
						if(gbar != 0.0) {
							inner_rhs_funcs.implicit_set_previous_terms();
							implicit_step(t);
						} else  {
							erk_step();
						}
					} else {
						// No support for methods where delta_C != 0 and gbar != 0 simultaneously.
						// Therefore assume all cases where delta_C != 0, gbar = 0, corresponding to ODE stage solve.
						explicit_solve(H,M,interval_start,interval_end);
					}
				} else if (coeffs->method_type == 1) {
					explicit_solve(H,M,interval_start,interval_end);
				}
				y_stages.col(stage_index) = v_H;
			}
		}
		y = y_stages.col(num_stages-1);

		int stage_index = num_stages-1;
		if (status == 0) {
			inner_rhs_funcs.set_function_dependent_data(H, t, stage_index, true);
			double gbar = inner_rhs_funcs.gamma_bar(stage_index+1,stage_index);
			if(gbar != 0.0) {
				implicit_step(t);
			} else {
				explicit_step(H, M);
			}
			y_hat = v_H;
		}

		ret->y = y;
		ret->ess = err_norm->compute_norm(y-y_hat);
		ret->esf = 0.0;
		ret->total_microtimesteps = total_microtimesteps;
		ret->total_successful_microtimesteps = total_successful_microtimesteps;
		ret->status = status;
	}

	void implicit_step(double t) {
		newton_solver.solve(t, &v_0, &newton_ret);
		problem->increment_slow_nonlinear_solves();
		v_H = newton_ret.y;
		inner_status = newton_ret.status;
		total_microtimesteps += 1;
		if (inner_status == 0) {
			total_successful_microtimesteps += 1;
		} else if (inner_status == 1) {
			// Newton nonconvergence.
			status = 3;
		} else if (inner_status == 2) {
			// Newton linear solver failure.
			status = 4;
		}
	}

	void explicit_step(double H, int M) {
		inner_rhs_funcs.explicit_set_previous_terms();
		dirk.solve(0, (double) H/M, &v_0, &output_tspan, controller, &dirk_ret);
		v_H = dirk_ret.Y.tail_cols(1);
		inner_status = dirk_ret.status;
		total_microtimesteps += dirk_ret.ts.size();
		if (inner_status == 0) {
			total_successful_microtimesteps += std::set<double>(dirk_ret.ts.begin(), dirk_ret.ts.end()).size();
		} else if (inner_status == 1) {
			// Inner ODE h_new too small.
			status = 1;
		} else if (inner_status == 2) {
			// Inner ODE h_new nonfinite.
			status = 2;
		}
	}

	void declare_vectors() {
		output_tspan = vec(1,fill::zeros);
		y = vec(problem_dimension, fill::zeros);
		y_hat = vec(problem_dimension, fill::zeros);
		v_0 = vec(problem_dimension, fill::zeros);
		v_H = vec(problem_dimension, fill::zeros);
		y_stages = mat(problem_dimension, num_stages, fill::zeros);
		inner_rhs_funcs.set_problem_dependent_data(&y_stages);
	}

	void prepare_solve(double H) {
		output_tspan(0) = H;
		v_0.zeros();
		v_H.zeros();
		y_stages.zeros();
	}
};

#endif