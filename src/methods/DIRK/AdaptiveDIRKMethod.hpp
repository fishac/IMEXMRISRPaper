#ifndef ADAPTIVEDIRKMETHOD_DEFINED__
#define ADAPTIVEDIRKMETHOD_DEFINED__

#include "AdaptiveStepSingleRateMethod.hpp"
#include "Residual.hpp"
#include "NewtonSolver.hpp"
#include "DIRKResidual.hpp"
#include "Controller.hpp"
#include "Problem.hpp"
#include <set>

using namespace arma;
using namespace std;

class AdaptiveDIRKMethod : AdaptiveStepSingleRateMethod {
public:
	SingleRateMethodCoefficients* coeffs;
	DIRKResidual dirk_residual;
	NewtonSolver newton_solver;
	struct NewtonSolverReturnValue newton_ret;
	WeightedErrorNorm* err_norm;
	bool is_top_level;
	int problem_dimension;
	vec explicit_data;
	mat Y;
	vec y;
	vec y_temp;
	vec y_stage;
	vec y_hat;
	vec y_accepted;
	mat y_stages;
	vec t_span;
	double h;
	double h_new;
	double err;
	double tol = 1.0;
	double safety_fac = 0.85;
	int output_index;
	int status = 0;
	int newton_status = 0;
	Problem* problem;

	AdaptiveDIRKMethod(SingleRateMethodCoefficients* coeffs_, Problem* problem_, int problem_dimension_, WeightedErrorNorm* err_norm_, bool is_top_level_) :
	dirk_residual(coeffs_, problem_, problem_dimension_),
	newton_solver(&(AdaptiveDIRKMethod::dirk_residual), 20, 0.1, problem_dimension_, err_norm_)
	{
		problem = problem_;
		coeffs = coeffs_;
		problem_dimension = problem_dimension_;
		err_norm = err_norm_;
		is_top_level = is_top_level_;
		
		declare_vectors();
	}

	void solve(double t_0, double h_, vec* y_0, vec* output_tspan, Controller* controller, AdaptiveSingleRateMethodReturnValue* ret) {
		prepare_solve(h_, output_tspan);
		output_index = 0;
		if (t_0 == (*output_tspan)(0)) {
			Y.col(0) = *y_0;
			output_index++;
		}
		y_accepted = *y_0;
		if (is_top_level) {
			err_norm->set_weights(y_0);
		}

		std::vector<double> ts;
		std::vector<double> hs;
		std::vector<int> full_function_evals;
		std::vector<int> fast_function_evals;
		std::vector<int> slow_function_evals;
		std::vector<int> implicit_function_evals;
		std::vector<int> explicit_function_evals;
		std::vector<int> full_jacobian_evals;
		std::vector<int> fast_jacobian_evals;
		std::vector<int> slow_jacobian_evals;
		std::vector<int> implicit_jacobian_evals;
		ts.push_back(t_0);
		hs.push_back(h);

		full_function_evals.push_back(0);
		fast_function_evals.push_back(0);
		slow_function_evals.push_back(0);
		implicit_function_evals.push_back(0);
		explicit_function_evals.push_back(0);
		full_jacobian_evals.push_back(0);
		fast_jacobian_evals.push_back(0);
		slow_jacobian_evals.push_back(0);
		implicit_jacobian_evals.push_back(0);

		controller->initialize(h,0);
		status = 0;
		newton_status = 0;

		double t = t_0;
		bool continue_computation = true;
		while(continue_computation) {
			y = y_accepted;
			compute_stages(t,h);
			if (newton_status == 0) {
				compute_solutions(h);
				err = err_norm->compute_norm(y-y_hat);
				
				controller->update_errors(err,0.0);
				h_new = controller->get_new_H();

				if (abs(t+h - (*output_tspan)(output_index)) > 1e-15 && t+h+h_new > (*output_tspan)(output_index)) {
					h_new = (*output_tspan)(output_index)-t-h;
				} else if (abs(t+h - (*output_tspan)(output_index)) < 1e-15 && output_index < output_tspan->n_elem-1) {
					if (t+h+h_new > (*output_tspan)(output_index+1)) {
						h_new = (*output_tspan)(output_index+1)-t-h;
					}
				}

				//Accept step.
				if (err < tol) {
					if(abs(t+h - (*output_tspan)(output_index)) < 1e-15 ) {
						Y.col(output_index) = y;
						output_index++;

						if (output_index == output_tspan->n_elem) {
							continue_computation = false;
						}
					} 
					
					if (is_top_level) {
						err_norm->set_weights(&y);
					}
					y_accepted = y;
					t += h;
				}	
				if (continue_computation && (abs(h_new) < 1e-15 || t+h_new == t)) {
					if (abs(h_new) < 1e-15) {
						//printf("DIRK H_new too small at t: %.16f\n", t);
					} else if (t+h_new == t) {
						//printf("DIRK t+H_new == t at t: %.16f\n",t);
					}
					continue_computation = false;
					status = 1;
				}
				if (continue_computation && !isfinite(h_new)) {
					//printf("DIRK Nonfinite H_new. Stopping solve at t: %.16f.\n",t);
					continue_computation = false;
					status = 2;
				}
			} else if (newton_status == 1) {
				h_new = h/2.0;
				newton_status = 0;
			} else if (newton_status == 2) {
				status = 3;
				continue_computation = false;
			}

			ts.push_back(t);
			hs.push_back(h);
			full_function_evals.push_back(problem->full_function_evals - full_function_evals.back());
			fast_function_evals.push_back(problem->fast_function_evals - fast_function_evals.back());
			slow_function_evals.push_back(problem->slow_function_evals - slow_function_evals.back());
			implicit_function_evals.push_back(problem->implicit_function_evals - implicit_function_evals.back());
			explicit_function_evals.push_back(problem->explicit_function_evals - explicit_function_evals.back());
			full_jacobian_evals.push_back(problem->full_jacobian_evals - full_jacobian_evals.back());
			fast_jacobian_evals.push_back(problem->fast_jacobian_evals - fast_jacobian_evals.back());
			slow_jacobian_evals.push_back(problem->slow_jacobian_evals - slow_jacobian_evals.back());
			implicit_jacobian_evals.push_back(problem->implicit_jacobian_evals - implicit_jacobian_evals.back());
				
			controller->update_H(h_new);
			h = h_new;
			set_problem_dependent_data(h);
		}
		
		ret->Y = Y;
		ret->ts = ts;
		ret->hs = hs;
		ret->total_timesteps = ts.size();
		ret->total_successful_timesteps = std::set<double>(ts.begin(), ts.end()).size();
		ret->full_function_evals = full_function_evals;
		ret->fast_function_evals = fast_function_evals;
		ret->slow_function_evals = slow_function_evals;
		ret->implicit_function_evals = implicit_function_evals;
		ret->explicit_function_evals = explicit_function_evals;
		ret->full_jacobian_evals = full_jacobian_evals;
		ret->fast_jacobian_evals = fast_jacobian_evals;
		ret->slow_jacobian_evals = slow_jacobian_evals;
		ret->implicit_jacobian_evals = implicit_jacobian_evals;
		ret->status = status;
	}

	void compute_stages(double t, double h) {
		for(int stage_idx=0; stage_idx<coeffs->num_stages; stage_idx++) {
			if (status == 0) {
				dirk_residual.set_function_dependent_data(stage_idx);

				y_stage.zeros();
				y_temp.zeros();
				compute_explicit_data(stage_idx);

				if (coeffs->A(stage_idx,stage_idx) != 0.0) {
					newton_solver.solve(t, &y, &newton_ret);
					y_temp = newton_ret.y;
					newton_status = newton_ret.status;
				} else {
					y_temp = y + h*explicit_data;
				}

				problem->full_rhs(t+(coeffs->c(stage_idx))*h,&y_temp,&y_stage);
				y_stages.col(stage_idx) = y_stage;
			}
		}
	}

	void compute_explicit_data(int stage_idx) {
		explicit_data.zeros();
		for(int inner_stage_idx = 0; inner_stage_idx<stage_idx; inner_stage_idx++) {
			explicit_data += (coeffs->A(stage_idx,inner_stage_idx))*y_stages.col(inner_stage_idx);
		}
	}

	void compute_solutions(double h) {
		y_temp.zeros();
		for(int stage_idx=0; stage_idx<(coeffs->num_stages); stage_idx++) {
			y_temp += (coeffs->b(stage_idx))*y_stages.col(stage_idx);
		}
		y = y_accepted + h*y_temp;

		y_temp.zeros();
		for(int stage_idx=0; stage_idx<(coeffs->num_stages); stage_idx++) {
			y_temp += (coeffs->d(stage_idx))*y_stages.col(stage_idx);
		}
		y_hat = y_accepted + h*y_temp;
	}

	void declare_vectors() {
		y_stages = mat(problem_dimension, coeffs->num_stages, fill::zeros);
		y = vec(problem_dimension, fill::zeros);
		y_temp = vec(problem_dimension, fill::zeros);
		y_stage = vec(problem_dimension, fill::zeros);
		y_hat = vec(problem_dimension, fill::zeros);
		y_accepted = vec(problem_dimension, fill::zeros);
		explicit_data = vec(problem_dimension, fill::zeros);
		dirk_residual.set_explicit_data_pointer(&explicit_data);
	}

	void prepare_solve(double h_, vec* output_tspan) {
		h = h_;
		set_problem_dependent_data(h_);

		Y = mat(problem_dimension, output_tspan->n_elem, fill::zeros);
		y.zeros();
		y_stages.zeros();
		y_temp.zeros();
		y_stage.zeros();
		y_hat.zeros();
		y_accepted.zeros();
		explicit_data.zeros();
	}

	void set_problem_dependent_data(double h_) {
		dirk_residual.set_problem_dependent_data(h_);
		newton_solver.set_problem_dependent_data(h_);
	}
};

#endif