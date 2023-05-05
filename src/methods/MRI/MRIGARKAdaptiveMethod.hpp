#ifndef MRIGARKADAPTIVEMETHOD_DEFINED__
#define MRIGARKADAPTIVEMETHOD_DEFINED__

#include "MRIGARKAdaptiveStep.hpp"
#include "AdaptiveStepMultiRateMethod.hpp"
#include "MRIGARKAdaptiveStep.hpp"
#include "MRIGARKFixedStep.hpp"
#include "Problem.hpp"
#include "WeightedErrorNorm.hpp"
#include "Controller.hpp"
#include "FastErrorMeasurementTypes.hpp"
#include <set>

using namespace arma;

class MRIGARKAdaptiveMethod : AdaptiveStepMultiRateMethod {
public:
	mat Y;
	vec y_accepted;
	vec y;
	vec y_hat;
	double H;
	double H_new;
	int M_new;
	double ess;
	double esf;
	double tol = 1.0;
	double H_recovery_factor = 0.3;
	int M;
	int output_index;
	int problem_dimension;
	WeightedErrorNorm* err_norm;
	Problem* problem;
	int status;

	MRIGARKAdaptiveMethod(Problem* problem_, int problem_dimension_, WeightedErrorNorm* err_norm_) {
		problem = problem_;
		problem_dimension = problem_dimension_;
		err_norm = err_norm_;
	}

	void solve(double t_0, double H_0, int M_0, vec* y_0, vec* output_tspan, MRIGARKAdaptiveStep* mrigark_step, Controller* controller, const char* measurement_type, AdaptiveMultiRateMethodReturnValue* ret) {
		solve(t_0, H_0, M_0, y_0, output_tspan, mrigark_step, controller, measurement_type, true, ret);
	}

	void solve(double t_0, double H_0, int M_0, vec* y_0, vec* output_tspan, MRIGARKAdaptiveStep* mrigark_step, Controller* controller, const char* measurement_type, bool record_step_data, AdaptiveMultiRateMethodReturnValue* ret) {
		prepare_solve(t_0, H_0, M_0, output_tspan);
		output_index = 0;
		if (t_0 == (*output_tspan)(0)) {
			Y.col(0) = *y_0;
			output_index++;
		}
		y_accepted = *y_0;
		err_norm->set_weights(y_0);
		problem->set_yhat_that(*y_0,t_0);

		std::deque<int> recent_statuses;
		std::deque<int> recent_Ms;

		int iterations_at_t = 0;
		int iterations_at_min_esf = 0;
		int iterations_at_max_M = 0;

		bool simple_step_failure = false;
		double previous_successful_H = -1.0;
		int previous_successful_M = -1;
		status = 0;
		int total_timesteps = 1;
		int total_microtimesteps = 0;
		int total_successful_microtimesteps = 0;
		std::vector<double> ts;
		std::vector<double> Hs;
		std::vector<int> effective_Ms;
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
		Hs.push_back(H);
		effective_Ms.push_back(M);

		full_function_evals.push_back(0);
		fast_function_evals.push_back(0);
		slow_function_evals.push_back(0);
		implicit_function_evals.push_back(0);
		explicit_function_evals.push_back(0);
		full_jacobian_evals.push_back(0);
		fast_jacobian_evals.push_back(0);
		slow_jacobian_evals.push_back(0);
		implicit_jacobian_evals.push_back(0);

		controller->initialize(H, M);

		double t = t_0;
		bool continue_computation = true;
		struct MRIGARKAdaptiveStepReturnValue step_ret;
		while(continue_computation) {
			simple_step_failure = false;
			y = y_accepted;
			mrigark_step->step_solution(t, H, M, &y, measurement_type, &step_ret);

			if (t > (*output_tspan)(output_tspan->n_elem-1)) {
				break;
			}
			if (step_ret.status == 0 && status == 0) {
				y = step_ret.y;
				ess = std::max(step_ret.ess,1e-6);
				esf = std::max(step_ret.esf,1e-6);
				total_microtimesteps += step_ret.total_microtimesteps;
				total_successful_microtimesteps += step_ret.total_successful_microtimesteps;

				//printf("t: %.16f, H: %.16f, t+H: %.16f, M: %d, ess: %.16f, esf: %.16f\n",t,H,t+H,M,ess,esf);

				if (iterations_at_t == 0) {
					controller->update_H(H);
					controller->update_M(M);
					controller->update_errors(ess,esf);
					H_new = controller->get_new_H();
					M_new = controller->get_new_M();
				} else {
					controller->replace_last_H(H);
					controller->replace_last_M(M);
					controller->replace_last_errors(ess,esf);
					H_new = controller->get_new_H();
					M_new = controller->get_new_M();
				}

				M_new = std::max(1,std::min(M_new,10000));
				//Accept step.
				if (ess+esf < tol) {
					if (H_new < 0) {
						break;
					}
					if(abs(t+H - (*output_tspan)(output_index)) < 1e-15 ) {
						//printf("Storing solution at t: %.16f\n",t+H);
						Y.col(output_index) = y;
						output_index++;

						if (output_index == output_tspan->n_elem) {
							continue_computation = false;
						}
					} 

					iterations_at_t = 0;
					y_accepted = y;
					t += H;
					err_norm->set_weights(&y);
					problem->set_yhat_that(y,t);
				} else {
					simple_step_failure = true;
				}

				if (continue_computation && (abs(H_new) < 1e-15 || t+H_new == t)) {
					continue_computation = false;
					status = 1;
				} else if (!isfinite(H_new)) {
					H_new = H/2.0;
					//continue_computation = false;
					//status = 2;
					//controller->print_status();
					//y.print("y");
					//y_accepted.print("y_accepted");
				} else if (iterations_at_t > 10) {
					continue_computation = false;
					status = 3;
				}

			} else if (step_ret.status == 1) {
				simple_step_failure = true;
			} else if (step_ret.status == 2) {
				continue_computation = false;
				status = 4;
			} else if (step_ret.status == 3) {
				simple_step_failure = true;
			} else if (step_ret.status == 4) {
				continue_computation = false;
				status = 5;
			}

			//printf("H_new: %.16f, M_new: %d\n",H_new,M_new);
			if (continue_computation) {
				H_new = std::min(H_new, (*output_tspan)(output_index)-t);
				if (simple_step_failure) {
					recent_statuses.push_back(1);
					iterations_at_t++;
					//printf("simple_step_failure\n");
					if (previous_successful_H > 0.0) {
						H_new = std::min(H_new, previous_successful_H);
					}
					if (previous_successful_M > 0) {
						M_new = std::max(M_new, previous_successful_M);
					}

					if (iterations_at_t > 0) {
						//H_new = std::min(H_new, H_recovery_factor*H);
						//int M_temp = std::ceil(M_new/H_recovery_factor);
						//M_new = std::max(M_new, M_temp);
					}
				} else {
					recent_statuses.push_back(0);
					previous_successful_H = H;
					previous_successful_M = M;
				}
				M_new = std::max(1,std::min(M_new,10000));
				recent_Ms.push_back(M);

				if (recent_Ms.size() > 20) {
					recent_Ms.pop_front();
				}
				if (recent_statuses.size() > 100) {
					recent_statuses.pop_front();
				}

				if (std::count(recent_Ms.begin(),recent_Ms.end(),10000) >= 10) {
					//continue_computation = false;
					//status = 6;
				} else if (ts.size() > 1e5) {
					// Only necessary in testing usually.
					//continue_computation = false;
					//status = 7;
				} else if (std::count(recent_statuses.begin(),recent_statuses.end(),1) >= 30) {
					//continue_computation = false;
					//status = 8;
				} else if (H_new/M_new < 1e-14) {
					continue_computation = false;
					status = 9;
				}
				int recent_failures = std::count(recent_statuses.begin(),recent_statuses.end(),1);
				//printf("t: %.16f, H_new: %.16f, M_new: %d, H_new/M_new: %.16f, recent failures: %d\n",t,H_new,M_new,H_new/M_new,recent_failures);

				total_timesteps++;
				if (record_step_data) {
					ts.push_back(t);
					Hs.push_back(H);
					full_function_evals.push_back(problem->full_function_evals - full_function_evals.back());
					fast_function_evals.push_back(problem->fast_function_evals - fast_function_evals.back());
					slow_function_evals.push_back(problem->slow_function_evals - slow_function_evals.back());
					implicit_function_evals.push_back(problem->implicit_function_evals - implicit_function_evals.back());
					explicit_function_evals.push_back(problem->explicit_function_evals - explicit_function_evals.back());
					full_jacobian_evals.push_back(problem->full_jacobian_evals - full_jacobian_evals.back());
					fast_jacobian_evals.push_back(problem->fast_jacobian_evals - fast_jacobian_evals.back());
					slow_jacobian_evals.push_back(problem->slow_jacobian_evals - slow_jacobian_evals.back());
					implicit_jacobian_evals.push_back(problem->implicit_jacobian_evals - implicit_jacobian_evals.back());
					effective_Ms.push_back(M);
				}

				H = H_new;
				M = M_new;
			}
		}
		ret->Y = Y;
		ret->ts = ts;
		ret->Hs = Hs;
		ret->Ms = effective_Ms;
		ret->total_timesteps = total_timesteps;
		ret->total_successful_timesteps = std::set<double>(ts.begin(), ts.end()).size();
		ret->total_microtimesteps = total_microtimesteps;
		ret->total_successful_microtimesteps = total_successful_microtimesteps;
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

	void declare_vectors() {
		y_accepted = vec(problem_dimension, fill::zeros);
		y = vec(problem_dimension, fill::zeros);
		y_hat = vec(problem_dimension, fill::zeros);
	}

	void prepare_solve(double t_0, double H_0, int M_0, vec* output_tspan) {
		H = H_0;
		M = M_0;
		y_accepted.zeros();
		y.zeros();
		y_hat.zeros();
		Y = mat(problem_dimension, output_tspan->n_elem, fill::zeros);
	}
};

#endif