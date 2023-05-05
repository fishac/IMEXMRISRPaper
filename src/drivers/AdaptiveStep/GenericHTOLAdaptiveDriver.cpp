#include <armadillo>
#include <math.h>
#include <stdio.h>

#include "BicouplingProblem.hpp"
#include "BicouplingDLProblem.hpp"
#include "BicouplingLNProblem.hpp"
#include "BrusselatorProblem.hpp"
#include "BrusselatorDLProblem.hpp"
#include "FourBody3dProblem.hpp"
#include "KapsProblem.hpp"
#include "KapsDLProblem.hpp"
#include "KapsLNProblem.hpp"
#include "KPRProblem.hpp"
#include "KPRDLProblem.hpp"
#include "LienardProblem.hpp"
#include "LienardDLProblem.hpp"
#include "LienardLNProblem.hpp"
#include "OregonatorProblem.hpp"
#include "OregonatorDLProblem.hpp"
#include "PleiadesProblem.hpp"
#include "Brusselator1DProblem.hpp"
#include "Bicoupling1DProblem.hpp"
#include "Problem.hpp"
#include "AdaptiveHTOLStepDriver.hpp"
#include "WeightedErrorNorm.hpp"
#include "MRIGARKERK33Coefficients.hpp"
#include "MRIGARKIRK21aCoefficients.hpp"
#include "MRIGARKERK45aCoefficients.hpp"
#include "MRIGARKESDIRK34aCoefficients.hpp"
#include "MRIGARKERK22aCoefficients.hpp"

using namespace std;
using namespace arma;

mat load_true_sol(const char* problem_name) {
	char filename[100];
	sprintf(filename,"./resources/%s/%s_truesol.csv",problem_name,problem_name);
	mat Y_true;
	bool success = Y_true.load(filename, csv_ascii);
	if (!success) {
		printf("Failed to load true solution.\n");
	}
	return Y_true;
}

mat get_true_sol(vec* output_tspan, Problem* problem) {
	if (problem->has_true_solution) {
		mat Y_true(problem->problem_dimension, output_tspan->n_elem, fill::zeros);
		vec y_true(problem->problem_dimension, fill::zeros);
		double t = 0.0;
		for(int it=0; it<output_tspan->n_elem; it++) {
			t = (*output_tspan)(it);
			problem->true_solution(t, &y_true);
			Y_true.col(it) = y_true;
		}
		return Y_true;
	} else {
		return load_true_sol(problem->name);
	}
}

void run(Problem* problem, MRICoefficients* coeffs, const char* tol_string, WeightedErrorNorm* err_norm, double tol) {
	AdaptiveStepDriver driver;
	int problem_dimension = problem->problem_dimension;
	double H_0 = std::pow(10.0,-3.0);
	double tolf_0 = 0.5*tol;
	vec output_tspan = linspace(problem->t_0, problem->t_f, 11);
	mat Y_true = get_true_sol(&output_tspan, problem);

	printf("\n%s Problem.\n", problem->name);
	driver.run(problem, coeffs, tol_string, err_norm, H_0, tolf_0, &output_tspan, &Y_true);
}

void setup_and_run(Problem* problem, const char* input_method_name, const char* tol_string, double tol) {
	vec atol = tol*vec(problem->problem_dimension,fill::ones);
	double rtol = tol;
	WeightedErrorNorm err_norm(&atol, rtol);

	if(strcmp("MRIGARKIRK21a",input_method_name) == 0) {
		MRIGARKIRK21aCoefficients mrigarkirk21a;
		run(problem, &mrigarkirk21a, tol_string, &err_norm, tol);
	} else if(strcmp("MRIGARKERK33",input_method_name) == 0) {
		MRIGARKERK33Coefficients mrigarkerk33;
		run(problem, &mrigarkerk33, tol_string, &err_norm, tol);
	} else if(strcmp("MRIGARKERK45a",input_method_name) == 0) {
		MRIGARKERK45aCoefficients mrigarkerk45a;
		run(problem, &mrigarkerk45a, tol_string, &err_norm, tol);
	} else if(strcmp("MRIGARKESDIRK34a",input_method_name) == 0) {
		MRIGARKESDIRK34aCoefficients mrigarkesdirk34a;
		run(problem, &mrigarkesdirk34a, tol_string, &err_norm, tol);
	} else if(strcmp("MRIGARKERK22a",input_method_name) == 0) {
		MRIGARKERK22aCoefficients mrigarkerk22a;
		run(problem, &mrigarkerk22a, tol_string, &err_norm, tol);
	} else {
		printf("Error: Did not recognize method name: %s\n", input_method_name);
	}
}

int main(int argc, char* argv[]) {
	if(argc != 4) {
		printf("Error: Requires 3 command-line arguments. Ex: GenericAdaptiveDriver.exe <ProblemName> <MethodName> <tol>\n");
		return 1;
	} else {
		double tol = 0.0;
		const char* input_problem_name = argv[1];
		const char* input_method_name = argv[2];
		const char* tol_string = argv[3];
		sscanf(argv[3], "%lf", &tol);

		if(strcmp("Bicoupling",input_problem_name) == 0) {
			BicouplingProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("BicouplingDL",input_problem_name) == 0) {
			BicouplingDLProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("BicouplingLN",input_problem_name) == 0) {
			BicouplingLNProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("Brusselator",input_problem_name) == 0) {
			BrusselatorProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("BrusselatorDL",input_problem_name) == 0) {
			BrusselatorDLProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("FourBody3d",input_problem_name) == 0) {
			FourBody3dProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("Kaps",input_problem_name) == 0) {
			KapsProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("KapsDL",input_problem_name) == 0) {
			KapsDLProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("KapsLN",input_problem_name) == 0) {
			KapsLNProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("KPR",input_problem_name) == 0) {
			KPRProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("KPRDL",input_problem_name) == 0) {
			KPRProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("Lienard",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("LienardDL",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("LienardLN",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("Oregonator",input_problem_name) == 0) {
			OregonatorProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("OregonatorDL",input_problem_name) == 0) {
			OregonatorProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("Pleiades",input_problem_name) == 0) {
			PleiadesProblem problem;
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("Brusselator1D",input_problem_name) == 0) {
			Brusselator1DProblem problem(100);
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else if(strcmp("Bicoupling1D",input_problem_name) == 0) {
			Bicoupling1DProblem problem(100);
			setup_and_run(&problem, input_method_name, tol_string, tol);
		} else {
			printf("Error: Did not recognize problem name: %s\n", input_problem_name);
			return 1;
		}
		return 0;
	}
}