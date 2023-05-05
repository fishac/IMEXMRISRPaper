#ifndef FIXEDSTEPDRIVER_DEFINED__
#define FIXEDSTEPDRIVER_DEFINED__

#include <armadillo>
#include <math.h>
#include <chrono>

#include "MRIGARKERK33Coefficients.hpp"
#include "MRIGARKIRK21aCoefficients.hpp"
#include "MRIGARKERK45aCoefficients.hpp"
#include "MRIGARKESDIRK34aCoefficients.hpp"
#include "MRIGARKESDIRK46aCoefficients.hpp"
#include "MRIGARKSDIRK33aCoefficients.hpp"
#include "MRIGARKERK22aCoefficients.hpp"
#include "MRIGARKTEST1Coefficients.hpp"
#include "MRIGARKTEST2Coefficients.hpp"
#include "IMEXMRI3aCoefficients.hpp"
#include "IMEXMRI3bCoefficients.hpp"
#include "IMEXMRI4Coefficients.hpp"
#include "IMEXMRITest1Coefficients.hpp"
#include "IMEXMRISRTEST1Coefficients.hpp"
#include "IMEXMRISRTEST2Coefficients.hpp"
#include "IMEXMRISRTEST3Coefficients.hpp"
#include "IMEXMRISRTEST4Coefficients.hpp"
#include "IMEXMRISRMERK2Coefficients.hpp"
#include "IMEXMRISRMERK3Coefficients.hpp"
#include "IMEXMRISRMERK4Coefficients.hpp"
#include "IMEXMRISRMERK5Coefficients.hpp"
#include "IMEXMRISR21Coefficients.hpp"
#include "IMEXMRISR32Coefficients.hpp"
#include "IMEXMRISR43Coefficients.hpp"
#include "IMEXMRIGARK32Coefficients.hpp"
#include "IMEXMRIGARK21Coefficients.hpp"
#include "MERK32aCoefficients.hpp"
#include "MERK5Coefficients.hpp"
#include "MRIGARKFixedMethod.hpp"
#include "MRIGARKFixedStep.hpp"
#include "Problem.hpp"
#include "FixedStepMultiRateMethod.hpp"
#include "ARK54ESDIRKCoefficients.hpp"
#include "ARK43ESDIRKCoefficients.hpp"
#include "ARK32ESDIRKCoefficients.hpp"
#include "SDIRK21Coefficients.hpp"
#include "ESDIRK324Coefficients.hpp"
#include "ESDIRK325Coefficients.hpp"
#include "ESDIRK436Coefficients.hpp"
#include "ESDIRK437Coefficients.hpp"
#include "DormandPrinceERKCoefficients.hpp"
#include "ZonneveldERKCoefficients.hpp"
#include "BogackiShampineERKCoefficients.hpp"
#include "HeunEulerERKCoefficients.hpp"
#include "WeightedErrorNorm.hpp"
#include "LieTrotterMethod.hpp"
#include "StrangMarchukMethod.hpp"

using namespace std;
using namespace arma;
using namespace std::chrono;

void save_stats(const char* problem_name, const char* method_name, int M, vec* errs, vec* H_vec, vec* fast_function_evals, vec* slow_function_evals, vec* implicit_function_evals, vec* explicit_function_evals, vec* fast_jacobian_evals, vec* slow_jacobian_evals, vec* implicit_jacobian_evals, vec* slow_nonlinear_solves, vec* runtimes) {
	char filename[120];
	sprintf(filename, "./output/%s/%s_PaperRuns_%s_M%d_stats.csv",problem_name, problem_name,method_name,M);
	mat output = join_rows(*H_vec, *errs);
	output = join_rows(output, *fast_function_evals);
	output = join_rows(output, *slow_function_evals);
	output = join_rows(output, *implicit_function_evals);
	output = join_rows(output, *explicit_function_evals);
	output = join_rows(output, *fast_jacobian_evals);
	output = join_rows(output, *slow_jacobian_evals);
	output = join_rows(output, *implicit_jacobian_evals);
	output = join_rows(output, *slow_nonlinear_solves);
	output = join_rows(output, *runtimes);
	output.save(filename, csv_ascii);
}

class PaperRunsDriver {
public:
	void run(Problem* problem, MRICoefficients* coeffs, SingleRateMethodCoefficients* inner_coeffs, vec* H_vec, int M, vec* output_tspan, mat* Y_true) {
		vec atol(problem->problem_dimension,fill::ones);
		atol = 1e-12*atol;
		double rtol = 1e-12;
		WeightedErrorNorm err_norm(&atol, rtol);
	
		MRIGARKFixedMethod mrigark_method(
			problem,
			problem->problem_dimension
		);
		MRIGARKFixedStep mrigark_step(
			coeffs, 
			inner_coeffs,
			problem,
			problem->problem_dimension, 
			&err_norm
		);
		run_solves(problem, &mrigark_method, &mrigark_step, coeffs->name, H_vec, M, output_tspan, Y_true);
	}

	void run_solves(Problem* problem, FixedStepMultiRateMethod* method, FixedStepMultiRateStep* step, const char* instance_name, vec* H_vec, int M, vec* output_tspan, mat* Y_true) {
		int total_Hs = H_vec->n_elem;

		vec errs(total_Hs, fill::zeros);
		vec fast_function_evals(total_Hs, fill::zeros);
		vec slow_function_evals(total_Hs, fill::zeros);
		vec implicit_function_evals(total_Hs, fill::zeros);
		vec explicit_function_evals(total_Hs, fill::zeros);
		vec fast_jacobian_evals(total_Hs, fill::zeros);
		vec slow_jacobian_evals(total_Hs, fill::zeros);
		vec implicit_jacobian_evals(total_Hs, fill::zeros);
		vec slow_nonlinear_solves(total_Hs, fill::zeros);
		vec runtimes(total_Hs, fill::zeros);

		vec p;

		printf("\n%s.\n",instance_name);
		for(int ih=0; ih<total_Hs; ih++) {
			double H = (*H_vec)(ih);

			problem->reset_eval_counts();
			auto start = high_resolution_clock::now();
			mat Y = method->solve(problem->t_0, H, M, &(problem->y_0), output_tspan, step);
			auto stop = high_resolution_clock::now();
			errs(ih) = abs(Y-*Y_true).max();
			fast_function_evals(ih) = problem->fast_function_evals;
			slow_function_evals(ih) = problem->slow_function_evals;
			implicit_function_evals(ih) = problem->implicit_function_evals;
			explicit_function_evals(ih) = problem->explicit_function_evals;
			fast_jacobian_evals(ih) = problem->fast_jacobian_evals;
			slow_jacobian_evals(ih) = problem->slow_jacobian_evals;
			implicit_jacobian_evals(ih) = problem->implicit_jacobian_evals;
			slow_nonlinear_solves(ih) = problem->slow_nonlinear_solves;
			runtimes(ih) = duration_cast<milliseconds>(stop-start).count();

			printf("H = %.16f, M = %d, Err = %.16f\n", H, M, errs(ih));
		}
		
		p = polyfit(arma::log10(*H_vec),arma::log10(errs),1);
		printf("\tOrder estimate: %f\n",p(0));

		save_stats(problem->name, instance_name, M, &errs, H_vec, &fast_function_evals, &slow_function_evals, &implicit_function_evals, &explicit_function_evals, &fast_jacobian_evals, &slow_jacobian_evals, &implicit_jacobian_evals, &slow_nonlinear_solves, &runtimes);
		
		errs.zeros();
		fast_function_evals.zeros();
		slow_function_evals.zeros();
		implicit_function_evals.zeros();
		explicit_function_evals.zeros();
		fast_jacobian_evals.zeros();
		slow_jacobian_evals.zeros();
		implicit_jacobian_evals.zeros();
		slow_nonlinear_solves.zeros();
		runtimes.zeros();
	}
	
	void run_lietrotter_method(Problem* problem, SingleRateMethodCoefficients* inner_coeffs, vec* H_vec, int M, vec* output_tspan, mat* Y_true) {
		vec atol(problem->problem_dimension,fill::ones);
		atol = 1e-14*atol;
		double rtol = 1e-14;
		WeightedErrorNorm err_norm(&atol, rtol);
		LieTrotterMethod method(inner_coeffs, problem, problem->problem_dimension, &err_norm);
		const char* instance_name = "LieTrotter";
		
		int total_Hs = H_vec->n_elem;

		vec errs(total_Hs, fill::zeros);
		vec fast_function_evals(total_Hs, fill::zeros);
		vec slow_function_evals(total_Hs, fill::zeros);
		vec implicit_function_evals(total_Hs, fill::zeros);
		vec explicit_function_evals(total_Hs, fill::zeros);
		vec fast_jacobian_evals(total_Hs, fill::zeros);
		vec slow_jacobian_evals(total_Hs, fill::zeros);
		vec implicit_jacobian_evals(total_Hs, fill::zeros);
		vec slow_nonlinear_solves(total_Hs, fill::zeros);
		vec runtimes(total_Hs, fill::zeros);

		vec p;

		printf("\n%s.\n",instance_name);
		for(int ih=0; ih<total_Hs; ih++) {
			double H = (*H_vec)(ih);
			double h = H/M;

			problem->reset_eval_counts();
			auto start = high_resolution_clock::now();
			mat Y = method.solve(problem->t_0, H, M, &(problem->y_0), output_tspan);
			auto stop = high_resolution_clock::now();
			errs(ih) = abs(Y-*Y_true).max();
			fast_function_evals(ih) = problem->fast_function_evals;
			slow_function_evals(ih) = problem->slow_function_evals;
			implicit_function_evals(ih) = problem->implicit_function_evals;
			explicit_function_evals(ih) = problem->explicit_function_evals;
			fast_jacobian_evals(ih) = problem->fast_jacobian_evals;
			slow_jacobian_evals(ih) = problem->slow_jacobian_evals;
			implicit_jacobian_evals(ih) = problem->implicit_jacobian_evals;
			slow_nonlinear_solves(ih) = problem->slow_nonlinear_solves;
			runtimes(ih) = duration_cast<milliseconds>(stop-start).count();

			printf("H = %.16f, M = %d, Err = %.16f\n", H, M, errs(ih));
		}
		
		p = polyfit(arma::log10(*H_vec),arma::log10(errs),1);
		printf("\tOrder estimate: %f\n",p(0));

		save_stats(problem->name, instance_name, M, &errs, H_vec, &fast_function_evals, &slow_function_evals, &implicit_function_evals, &explicit_function_evals, &fast_jacobian_evals, &slow_jacobian_evals, &implicit_jacobian_evals, &slow_nonlinear_solves, &runtimes);
		
		errs.zeros();
		fast_function_evals.zeros();
		slow_function_evals.zeros();
		implicit_function_evals.zeros();
		explicit_function_evals.zeros();
		fast_jacobian_evals.zeros();
		slow_jacobian_evals.zeros();
		implicit_jacobian_evals.zeros();
		slow_nonlinear_solves.zeros();
		runtimes.zeros();
	}
	
	void run_strangmarchuk_method(Problem* problem, SingleRateMethodCoefficients* inner_coeffs, vec* H_vec, int M, vec* output_tspan, mat* Y_true) {
		vec atol(problem->problem_dimension,fill::ones);
		atol = 1e-14*atol;
		double rtol = 1e-14;
		WeightedErrorNorm err_norm(&atol, rtol);
		StrangMarchukMethod method(inner_coeffs, problem, problem->problem_dimension, &err_norm);
		const char* instance_name = "StrangMarchuk";
		
		int total_Hs = H_vec->n_elem;

		vec errs(total_Hs, fill::zeros);
		vec fast_function_evals(total_Hs, fill::zeros);
		vec slow_function_evals(total_Hs, fill::zeros);
		vec implicit_function_evals(total_Hs, fill::zeros);
		vec explicit_function_evals(total_Hs, fill::zeros);
		vec fast_jacobian_evals(total_Hs, fill::zeros);
		vec slow_jacobian_evals(total_Hs, fill::zeros);
		vec implicit_jacobian_evals(total_Hs, fill::zeros);
		vec slow_nonlinear_solves(total_Hs, fill::zeros);
		vec runtimes(total_Hs, fill::zeros);

		vec p;

		printf("\n%s.\n",instance_name);
		for(int ih=0; ih<total_Hs; ih++) {
			double H = (*H_vec)(ih);
			double h = H/M;

			problem->reset_eval_counts();
			auto start = high_resolution_clock::now();
			mat Y = method.solve(problem->t_0, H, M, &(problem->y_0), output_tspan);
			auto stop = high_resolution_clock::now();
			errs(ih) = abs(Y-*Y_true).max();
			fast_function_evals(ih) = problem->fast_function_evals;
			slow_function_evals(ih) = problem->slow_function_evals;
			implicit_function_evals(ih) = problem->implicit_function_evals;
			explicit_function_evals(ih) = problem->explicit_function_evals;
			fast_jacobian_evals(ih) = problem->fast_jacobian_evals;
			slow_jacobian_evals(ih) = problem->slow_jacobian_evals;
			implicit_jacobian_evals(ih) = problem->implicit_jacobian_evals;
			slow_nonlinear_solves(ih) = problem->slow_nonlinear_solves;
			runtimes(ih) = duration_cast<milliseconds>(stop-start).count();

			printf("H = %.16f, M = %d, Err = %.16f\n", H, M, errs(ih));
		}
		
		p = polyfit(arma::log10(*H_vec),arma::log10(errs),1);
		printf("\tOrder estimate: %f\n",p(0));

		save_stats(problem->name, instance_name, M, &errs, H_vec, &fast_function_evals, &slow_function_evals, &implicit_function_evals, &explicit_function_evals, &fast_jacobian_evals, &slow_jacobian_evals, &implicit_jacobian_evals, &slow_nonlinear_solves, &runtimes);
		
		errs.zeros();
		fast_function_evals.zeros();
		slow_function_evals.zeros();
		implicit_function_evals.zeros();
		explicit_function_evals.zeros();
		fast_jacobian_evals.zeros();
		slow_jacobian_evals.zeros();
		implicit_jacobian_evals.zeros();
		slow_nonlinear_solves.zeros();
		runtimes.zeros();
	}
};

#endif