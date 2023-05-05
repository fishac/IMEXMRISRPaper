#ifndef ADAPTIVESTEPDRIVER_DEFINED__
#define ADAPTIVESTEPDRIVER_DEFINED__

#include "Problem.hpp"
#include "AdaptiveStepMultiRateMethod.hpp"
#include "MRIGARKAdaptiveMethod.hpp"
#include "MRIGARKAdaptiveStepSlowMeasurement.hpp"
#include "MRIGARKAdaptiveStepFastMeasurement.hpp"
#include "MRICoefficients.hpp"
#include "MRIGARKERK33Coefficients.hpp"
#include "MRIGARKIRK21aCoefficients.hpp"
#include "MRIGARKERK45aCoefficients.hpp"
#include "MRIGARKESDIRK34aCoefficients.hpp"
#include "MERK32aCoefficients.hpp"
#include "SingleRateMethodCoefficients.hpp"
#include "HeunEulerERKCoefficients.hpp"
#include "BogackiShampineERKCoefficients.hpp"
#include "ZonneveldERKCoefficients.hpp"
#include "WeightedErrorNorm.hpp"
#include "ConstantConstantController.hpp"
#include "LinearLinearController.hpp"
#include "PIMRController.hpp"
#include "QuadraticQuadraticController.hpp"
#include "PIDMRController.hpp"
#include "IController.hpp"
#include "PIController.hpp"
#include "PIDController.hpp"
#include "GustafssonController.hpp"
#include "Controller.hpp"
#include "FastErrorMeasurementTypes.hpp"

using namespace std;
using namespace arma;
using namespace std::chrono;

struct stats {
	int total_timesteps; int total_successful_timesteps; int total_microtimesteps; 
	int total_successful_microtimesteps; double rel_err; double abs_err;
	int full_function_evals; int fast_function_evals; int slow_function_evals; 
	int implicit_function_evals;  int explicit_function_evals; 
	int full_jacobian_evals; int fast_jacobian_evals; int slow_jacobian_evals; 
	int implicit_jacobian_evals; int slow_nonlinear_solves; double runtime;
	int status;
};

struct stats_over_time {
	std::vector<double> ts; std::vector<double> Hs; std::vector<int> Ms;
	std::vector<int> full_function_evals; std::vector<int> fast_function_evals; std::vector<int> slow_function_evals;
	std::vector<int> implicit_function_evals; std::vector<int> explicit_function_evals;
	std::vector<int> full_jacobian_evals; std::vector<int> fast_jacobian_evals; std::vector<int> slow_jacobian_evals;
	std::vector<int> implicit_jacobian_evals;
};

class AdaptiveStepDriver {
public:
	void save_stats(const char* problem_name, const char* method_name, const char* tol_string, const char* controller_name, const char* measurement_type, stats* solve_stats) {
		char filename[200];
		sprintf(filename, "./output/%s/%s_AdaptiveStep_%s_%s_%s_%s_stats.csv", problem_name, problem_name, tol_string, controller_name, measurement_type, method_name);
		vec output = {
			(double) solve_stats->total_timesteps, 
			(double) solve_stats->total_successful_timesteps, 
			(double) solve_stats->total_microtimesteps, 
			(double) solve_stats->total_successful_microtimesteps, 
			solve_stats->rel_err, 
			solve_stats->abs_err,
			(double) solve_stats->full_function_evals, 
			(double) solve_stats->fast_function_evals, 
			(double) solve_stats->slow_function_evals,
			(double) solve_stats->implicit_function_evals, 
			(double) solve_stats->explicit_function_evals, 
			(double) solve_stats->full_jacobian_evals, 
			(double) solve_stats->fast_jacobian_evals, 
			(double) solve_stats->slow_jacobian_evals, 
			(double) solve_stats->implicit_jacobian_evals, 
			(double) solve_stats->slow_nonlinear_solves, 
			(double) solve_stats->runtime,
			(double) solve_stats->status
		};
		output.save(filename, csv_ascii);
	}
	
	void save_stats_over_time(const char* problem_name, const char* method_name, const char* tol_string, const char* controller_name, const char* measurement_type, stats_over_time* solve_stats) {
		char filename[200];
		sprintf(filename, "./output/%s/%s_AdaptiveStep_%s_%s_%s_%s_SOT.csv", problem_name, problem_name, tol_string, controller_name, measurement_type, method_name);
		int n_elem = (solve_stats->ts).size();
		mat output(n_elem,12,fill::zeros);
		for(int i=0; i<n_elem; i++) {
			output(i,0) = (solve_stats->ts)[i];
			output(i,1) = (solve_stats->Hs)[i];
			output(i,2) = (solve_stats->Ms)[i];
			output(i,3) = (solve_stats->full_function_evals)[i];
			output(i,4) = (solve_stats->fast_function_evals)[i];
			output(i,5) = (solve_stats->slow_function_evals)[i];
			output(i,6) = (solve_stats->implicit_function_evals)[i];
			output(i,7) = (solve_stats->explicit_function_evals)[i];
			output(i,8) = (solve_stats->full_jacobian_evals)[i];
			output(i,9) = (solve_stats->fast_jacobian_evals)[i];
			output(i,10) = (solve_stats->slow_jacobian_evals)[i];
			output(i,11) = (solve_stats->implicit_jacobian_evals)[i];
		}
		output.save(filename, csv_ascii);
	}

	void run(Problem* problem, MRICoefficients* coeffs, SingleRateMethodCoefficients* inner_coeffs, const char* input_controller_name, const char* tol_string, WeightedErrorNorm* err_norm, double H_0, int M_0, vec* output_tspan, mat* Y_true) {
		
		double k1_CC[1] = { 0.42 }; 
		double k2_CC[1] = { 0.44 }; 
		ConstantConstantController CCcontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_CC,
			k2_CC
		);
		
		if (strcmp(input_controller_name,"ConstantConstant")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &CCcontroller, H_0, M_0, Y_true, output_tspan);
		}
		double k1_LL[2] = { 0.82, 0.54 }; 
		double k2_LL[2] = { 0.94, 0.9 }; 
		LinearLinearController LLcontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_LL,
			k2_LL,
			&CCcontroller
		);

		if (strcmp(input_controller_name,"LinearLinear")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &LLcontroller, H_0, M_0, Y_true, output_tspan);
		}
		
		double k1_PIMR[2] = { 0.18, 0.86 };
		double k2_PIMR[2] = { 0.34, 0.8 }; 
		PIMRController pimrcontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_PIMR,
			k2_PIMR,
			&CCcontroller
		);

		if (strcmp(input_controller_name,"PIMR")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &pimrcontroller, H_0, M_0, Y_true, output_tspan);
		}
		
		double k1_PIDMR[3] = { 0.34, 0.1, 0.78 }; 
		double k2_PIDMR[3] = { 0.46, 0.42, 0.74 }; 
		PIDMRController pidmrcontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_PIDMR,
			k2_PIDMR,
			&CCcontroller
		);

		if (strcmp(input_controller_name,"PIDMR")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &pidmrcontroller, H_0, M_0, Y_true, output_tspan);
		}
		
		double k1_i[1] = { 1.0 };
		double k2_i[1] = { 0.0 }; 
		IController icontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_i,
			k2_i
		);
		
		if (strcmp(input_controller_name,"I")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &icontroller, H_0, M_0, Y_true, output_tspan);
		}
		
		double k1_pi[2] = { 0.6, 0.2 };
		double k2_pi[2] = { 0.0, 0.0 }; 
		PIController picontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_pi,
			k2_pi
		);
		
		
		if (strcmp(input_controller_name,"PI")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &picontroller, H_0, M_0, Y_true, output_tspan);
		}
		
		double k1_pid[3] = { 0.49, 0.34, 0.1 };
		double k2_pid[3] = { 0.0, 0.0, 0.0 }; 
		PIDController pidcontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_pid,
			k2_pid
		);
		
		if (strcmp(input_controller_name,"PID")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &pidcontroller, H_0, M_0, Y_true, output_tspan);
		}
		
		double k1_g[2] = { 0.6, 0.2 };
		double k2_g[2] = { 0.0, 0.0 }; 
		GustafssonController gcontroller(
			1.0,
			1.0,
			1.0,
			0.85,
			k1_g,
			k2_g,
			&icontroller
		);
		
		if (strcmp(input_controller_name,"Gustafsson")==0) {
			run_single_mrigark_method(problem, coeffs, inner_coeffs, tol_string, err_norm, &gcontroller, H_0, M_0, Y_true, output_tspan);
		}
	}

	void run_single_mrigark_method(Problem* problem, MRICoefficients* coeffs, SingleRateMethodCoefficients* inner_coeffs, const char* tol_string, WeightedErrorNorm* err_norm,  Controller* controller, double H_0, int M_0, mat* Y_true, vec* output_tspan) {
		controller->reset();
		double P = std::min(coeffs->primary_order,coeffs->secondary_order);
		double p = std::min(inner_coeffs->primary_order,inner_coeffs->secondary_order);
		controller->set_orders(P,p);
		const char* measurement_type = "LASA-mean";
		//const char* measurement_type = "FS";
		//const char* measurement_type = "LASA-max";
		AdaptiveMultiRateMethodReturnValue ret;
		MRIGARKAdaptiveMethod mrigark_method(
			problem,
			problem->problem_dimension,
			err_norm
		);
		MRIGARKAdaptiveStepFastMeasurement mrigark_step_fm(
			coeffs, 
			inner_coeffs, 
			problem,
			problem->problem_dimension, 
			err_norm
		);
		*output_tspan = { problem->t_f };
		auto start = high_resolution_clock::now();
		mrigark_method.solve(problem->t_0, H_0, M_0, &(problem->y_0), output_tspan, &mrigark_step_fm, controller, measurement_type, &ret);
		auto stop = high_resolution_clock::now();
		double runtime = duration_cast<milliseconds>(stop-start).count();
		process_multirate_method(&ret, problem, coeffs->name, tol_string, controller->name, measurement_type, Y_true, runtime);

	}

	void process_multirate_method(AdaptiveMultiRateMethodReturnValue* ret, Problem* problem, const char* instance_name, const char* tol_string, const char* controller_name, const char* measurement_type, mat* Y_true, double runtime) {
		std::vector<double> ts = ret->ts;
		std::vector<double> Hs = ret->Hs;
		std::vector<int> Ms = ret->Ms;
		mat Y = ret->Y;
		int status = ret->status;
		
		double abs_err = abs(Y_true->col(Y_true->n_cols-1)-Y).max();
		double rel_err = norm(Y_true->col(Y_true->n_cols-1)-Y,2)/norm(Y_true->col(Y_true->n_cols-1),2);

		if (status == 0) {
			printf("%s, %s. Measurement Type: %s. Total timesteps: %d, total microtimesteps: %d, total successful timesteps: %d, rel err: %.16f, abs err: %.16f\n", instance_name, controller_name, measurement_type, ret->total_timesteps, ret->total_microtimesteps, ret->total_successful_timesteps, rel_err, abs_err);
		} else if (status == 1) {
			printf("%s, %s. Measurement Type: %s. Solver failure: H_new too small.\n", instance_name, controller_name, measurement_type);
		} else if (status == 2) {
			printf("%s, %s. Measurement Type: %s. Solver failure: H_new nonfinite.\n", instance_name, controller_name, measurement_type);
		} else if (status == 3) {
			printf("%s, %s. Measurement Type: %s. Solver failure: Excessive iterations without progressing in time.\n", instance_name, controller_name, measurement_type);
		} else if (status == 4) {
			printf("%s, %s. Measurement Type: %s. Solver failure: Inner solver h_new nonfinite.\n", instance_name, controller_name, measurement_type);
		} else if (status == 5) {
			printf("%s, %s. Measurement Type: %s. Solver failure: NewtonSolver linear solver failure.\n", instance_name, controller_name, measurement_type);
		} else if (status == 6) {
			printf("%s, %s. Measurement Type: %s. Solver failure: Excessive iterations with M=M_max, esf=esf_min.\n", instance_name, controller_name, measurement_type);
		} else if (status == 7) {
			printf("%s, %s. Measurement Type: %s. Solver failure: Excessive iterations.\n", instance_name, controller_name, measurement_type);
		} else if (status == 8) {
			printf("%s, %s. Measurement Type: %s. Solver failure: Excessive failures in recent steps.\n", instance_name, controller_name, measurement_type);
		}
		
		struct stats solve_stats = {
			ret->total_timesteps, ret->total_successful_timesteps, ret->total_microtimesteps, 
			ret->total_successful_microtimesteps, rel_err, abs_err,
			problem->full_function_evals, problem->fast_function_evals, problem->slow_function_evals,
			problem->implicit_function_evals, problem->explicit_function_evals, 
			problem->full_jacobian_evals, problem->fast_jacobian_evals, problem->slow_jacobian_evals, 
			problem->implicit_jacobian_evals, problem->slow_nonlinear_solves, runtime, status
		};
		save_stats(problem->name, instance_name, tol_string, controller_name, measurement_type, &solve_stats);
			
		struct stats_over_time solve_stats_over_time = {
			ret->ts, ret->Hs, ret->Ms,
			ret->full_function_evals, ret->fast_function_evals, ret->slow_function_evals,
			ret->implicit_function_evals, ret->explicit_function_evals,
			ret->full_jacobian_evals, ret->fast_jacobian_evals, ret->slow_jacobian_evals,
			ret->implicit_jacobian_evals
		};
		save_stats_over_time(problem->name, instance_name, tol_string, controller_name, measurement_type, &solve_stats_over_time);
		problem->reset_eval_counts();
	}
};

#endif