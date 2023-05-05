#ifndef MRIGARKADAPTIVESTEPFASTMEASUREMENT_DEFINED__
#define MRIGARKADAPTIVESTEPFASTMEASUREMENT_DEFINED__

#include "MRIGARKAdaptiveStep.hpp"
#include "MRICoefficients.hpp"
#include "MRIGARKInnerRHSFunctions.hpp"
#include "MRIGARKExplicitSolveProblem.hpp"
#include "MRIGARKImplicitSolveResidual.hpp"
#include "AdaptiveDIRKMethod.hpp"
#include "SingleRateMethodCoefficients.hpp"
#include "FixedDIRKMethodFastMeasurement.hpp"
#include "NewtonSolver.hpp"
#include "WeightedErrorNorm.hpp"
#include "FastErrorMeasurementTypes.hpp"
#include <set>

using namespace arma;

class MRIGARKAdaptiveStepFastMeasurement : public MRIGARKAdaptiveStep {
public:
	Problem* problem;
	MRICoefficients* coeffs;
	MRIGARKInnerRHSFunctions inner_rhs_funcs;
	MRIGARKExplicitSolveProblem explicit_solve_problem;
	MRIGARKImplicitSolveResidual implicit_solve_residual;
	NewtonSolver newton_solver;
	struct NewtonSolverReturnValue newton_ret;
	FixedDIRKMethodFastMeasurement dirk;
	WeightedErrorNorm* err_norm;
	mat y_stages;
	vec y;
	vec y_hat;
	vec v_0;
	vec v_H;
	vec errs;
	vec err_vec;
	double H;
	int M;
	int total_microtimesteps;
	int total_successful_microtimesteps;
	int problem_dimension;
	int num_stages;
	int status;
	int inner_status;

	MRIGARKAdaptiveStepFastMeasurement(MRICoefficients* coeffs_, SingleRateMethodCoefficients* inner_coeffs_, Problem* problem_, int problem_dimension_, WeightedErrorNorm* err_norm_) :
	inner_rhs_funcs(coeffs_, problem_, problem_dimension_),
	implicit_solve_residual(&(MRIGARKAdaptiveStepFastMeasurement::inner_rhs_funcs)),
	explicit_solve_problem(&(MRIGARKAdaptiveStepFastMeasurement::inner_rhs_funcs)),
	dirk(inner_coeffs_, &(MRIGARKAdaptiveStepFastMeasurement::explicit_solve_problem), problem_dimension_, err_norm_),
	newton_solver(&(MRIGARKAdaptiveStepFastMeasurement::implicit_solve_residual), 20, 0.1, problem_dimension_, err_norm_)
	{
		problem = problem_;
		coeffs = coeffs_;
		num_stages = coeffs->num_stages;
		problem_dimension = problem_dimension_;
		err_norm = err_norm_;
		declare_vectors();
	}

	// Adaptive step solution
	void step_solution(double t, double H, int M, vec* y_0, const char* measurement_type, MRIGARKAdaptiveStepReturnValue* ret) {
		prepare_solve(H);
		inner_rhs_funcs.reset_stage_func_eval_storage();
		inner_rhs_funcs.set_step_dependent_data(H, t);
		y_stages.col(0) = *y_0;
		inner_rhs_funcs.store_stage_func_eval(y_0,0);
		//printf("step_solution. t: %.16f, H: %.16f, M: %d\n",t,H,M);
		total_microtimesteps = 0;
		total_successful_microtimesteps = 0;
		status = 0;
		
		double interval_start;
		double interval_end;
		int stage_initial_condition_index;
		int stage_index;
		
		double ess;
		
		for(std::vector<int> stage_group : coeffs->stage_groups) {
			interval_start = t;
			interval_end = 0.0;
			stage_initial_condition_index = 0;
			for(int stage_group_index=0; stage_group_index<stage_group.size(); stage_group_index++) {
				if (status == 0) {
					stage_index = stage_group[stage_group_index];
					inner_rhs_funcs.set_stage_dependent_data(stage_index, false);
						
					if (coeffs->method_type == 0) {
						interval_start = t+coeffs->c(stage_index-1)*H;
						interval_end = t+coeffs->c(stage_index)*H;
						stage_initial_condition_index = stage_index-1;
					} else if (coeffs->method_type == 1) {
						if (stage_group_index == 0) {
							interval_start = 0.0;
						}
						interval_end = coeffs->c(stage_index)*H;
					} else if (coeffs->method_type == 2) {
						//interval_start = t;
						interval_end = t+coeffs->c(stage_index)*H;
					} 
					v_0 = y_stages.col(stage_initial_condition_index);
					
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
							explicit_solve(H,M,stage_index,interval_start,interval_end);
						}
					} else if (coeffs->method_type == 1) {
						explicit_solve(H,M,stage_index,interval_start,interval_end);
					} else if (coeffs->method_type == 2) { // IMEXMRISR
						double c = coeffs->c(stage_index);
						explicit_solve(H,M,stage_index,interval_start,interval_end);
						v_0 = v_H;
						inner_rhs_funcs.implicit_set_previous_terms();
						double gdiag = inner_rhs_funcs.gamma_bar(stage_index,stage_index);
						if (gdiag > 1e-8) {
							implicit_step(interval_end);
						} else {
							erk_step();
						}
					}
					y_stages.col(stage_index) = v_H;
					
					interval_start = interval_end;
					stage_initial_condition_index = stage_index;
					
					// If not last stage, store function evals of the stage
					if (stage_index < num_stages-1) {
						inner_rhs_funcs.store_stage_func_eval(&v_H,stage_index);
					}
				}
			}
		}
		y = y_stages.col(num_stages-1);

		// Compute embedding
		if (status == 0) {
			int stage_index = num_stages-1;
			inner_rhs_funcs.set_stage_dependent_data(stage_index, true);
			if (coeffs->method_type == 0) {
				interval_start = t+coeffs->c(stage_index-1)*H;
				interval_end = t+coeffs->c(stage_index)*H;
				stage_initial_condition_index = stage_index-1;
				v_0 = y_stages.col(stage_initial_condition_index);
			} else if (coeffs->method_type == 1) {
				interval_start = 0.0;
				interval_end = H;
				stage_initial_condition_index = 0;
				v_0 = y_stages.col(stage_initial_condition_index);
			} else if (coeffs->method_type == 2) {
				interval_start = t;
				interval_end = t+H;
				stage_initial_condition_index = 0;
				v_0 = y_stages.col(stage_initial_condition_index);
			}
			
			if (coeffs->method_type == 0) {
				double gbar = inner_rhs_funcs.gamma_bar(stage_index+1,stage_index);
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
					explicit_solve(H,M,-1,interval_start,interval_end);
				}
			} else if (coeffs->method_type == 1) {
				explicit_solve(H,M,-1,interval_start,interval_end);
			} else if (coeffs->method_type == 2) { // IMEXMRISR
				double c = coeffs->c(stage_index);
				explicit_solve(H,M,-1,interval_start,interval_end);
				v_0 = v_H;
				inner_rhs_funcs.implicit_set_previous_terms();
				double gdiag = inner_rhs_funcs.gamma_bar(stage_index+1,stage_index);
				if (gdiag > 1e-8) {
					implicit_step(interval_end);
				} else {
					erk_step();
				}
			}
			y_hat = v_H;
			
			//if (coeffs->method_type == 0) {
				ess = err_norm->compute_norm(y-y_hat);
			//} else if (coeffs->method_type == 1) {
			//	ess = err_norm->compute_norm(y_hat);
			//}
		}

		ret->y = y;
		ret->ess = ess;
		ret->esf = compute_esf(measurement_type);
		ret->total_microtimesteps = total_microtimesteps;
		ret->total_successful_microtimesteps = total_successful_microtimesteps;
		ret->status = status;

		//printf("step slow measurement ess: %.16f, esf: %.16f, t: %.16f, H: %.16f, M: %d\n",ret->ess,ret->esf,t,H,M);
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
	
	void erk_step() {
		inner_rhs_funcs.erk_step(&v_0, &v_H);
	}

	void explicit_solve(double H, int M, int stage_index, double interval_start, double interval_end) {
		struct FixedDIRKMethodFastMeasurementReturnValue ret;
		vec output_tspan = { interval_end };
		dirk.solve(interval_start, H/M, &v_0, &output_tspan, &ret);
		v_H = ret.y;
		if (stage_index > 0) {
			errs(stage_index) = ret.err_estimate;
		}
		inner_status = ret.status;
		total_microtimesteps += M;
		total_successful_microtimesteps += M;
	}

	double compute_esf(const char* measurement_type) {
		if (FastError::is_LASAsum(measurement_type)) {
			// Aggregator: sum
			return accu(errs);
		} else if (FastError::is_LASAmean(measurement_type)) {
			// Aggregator: average
			return 1.0/num_stages*accu(errs);
		} else if (FastError::is_LASAmax(measurement_type)) {
			// Aggregator: max
			return max(errs);
		} else {
			return 0.0;
		}
	}

	void declare_vectors() {
		y = vec(problem_dimension, fill::zeros);
		y_hat = vec(problem_dimension, fill::zeros);
		v_0 = vec(problem_dimension, fill::zeros);
		v_H = vec(problem_dimension, fill::zeros);
		errs = vec(num_stages, fill::zeros);
		y_stages = mat(problem_dimension, num_stages, fill::zeros);
	}

	void prepare_solve(double H) {
		v_0.zeros();
		v_H.zeros();
		y_stages.zeros();
	}
};

#endif