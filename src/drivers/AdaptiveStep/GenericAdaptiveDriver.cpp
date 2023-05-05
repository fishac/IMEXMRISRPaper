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
#include "Brusselator1DIMEXProblem.hpp"
#include "Bicoupling1DProblem.hpp"
#include "Bicoupling1DIMEXProblem.hpp"
#include "Problem.hpp"
#include "AdaptiveStepDriver.hpp"
#include "WeightedErrorNorm.hpp"
#include "MRIGARKERK33Coefficients.hpp"
#include "MRIGARKIRK21aCoefficients.hpp"
#include "MRIGARKERK45aCoefficients.hpp"
#include "MRIGARKESDIRK34aCoefficients.hpp"
#include "MRIGARKERK22aCoefficients.hpp"
#include "IMEXMRITest1Coefficients.hpp"
#include "IMEXMRISR21Coefficients.hpp"
#include "IMEXMRISR32Coefficients.hpp"
#include "IMEXMRISR43Coefficients.hpp"
#include "IMEXMRIGARK21Coefficients.hpp"
#include "IMEXMRIGARK32Coefficients.hpp"
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

void run(Problem* problem, MRICoefficients* coeffs, SingleRateMethodCoefficients* inner_coeffs, const char* input_controller_name, const char* tol_string, WeightedErrorNorm* err_norm) {
	AdaptiveStepDriver driver;
	int problem_dimension = problem->problem_dimension;
	double H_0 = std::pow(10.0,-3.0);
	int M_0 = 7;
	vec output_tspan = linspace(problem->t_0, problem->t_f, 11);
	mat Y_true = get_true_sol(&output_tspan, problem);

	printf("\n%s Problem.\n", problem->name);
	driver.run(problem, coeffs, inner_coeffs, input_controller_name, tol_string, err_norm, H_0, M_0, &output_tspan, &Y_true);
}

void setup_and_run(Problem* problem, const char* input_method_name, SingleRateMethodCoefficients* inner_coeffs, const char* input_controller_name, const char* tol_string, double tol) {
	vec atol = tol*vec(problem->problem_dimension,fill::ones);
	double rtol = tol;
	WeightedErrorNorm err_norm(&atol, rtol);

	MRIGARKAdaptiveMethod mrigark_method(
		problem,
		problem->problem_dimension,
		&err_norm
	);
	if(strcmp("MRIGARKIRK21a",input_method_name) == 0) {
		MRIGARKIRK21aCoefficients mrigarkirk21a;
		run(problem, &mrigarkirk21a, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("MRIGARKERK33",input_method_name) == 0) {
		MRIGARKERK33Coefficients mrigarkerk33;
		run(problem, &mrigarkerk33, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("MRIGARKERK45a",input_method_name) == 0) {
		MRIGARKERK45aCoefficients mrigarkerk45a;
		run(problem, &mrigarkerk45a, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("MRIGARKESDIRK34a",input_method_name) == 0) {
		MRIGARKESDIRK34aCoefficients mrigarkesdirk34a;
		run(problem, &mrigarkesdirk34a, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("MRIGARKERK22a",input_method_name) == 0) {
		MRIGARKERK22aCoefficients mrigarkerk22a;
		run(problem, &mrigarkerk22a, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("IMEXMRITest1",input_method_name) == 0) {
		IMEXMRITest1Coefficients imexmritest1;
		run(problem, &imexmritest1, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("IMEXMRISR21",input_method_name) == 0) {
		IMEXMRISR21Coefficients imexmrisr21;
		run(problem, &imexmrisr21, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("IMEXMRISR32",input_method_name) == 0) {
		IMEXMRISR32Coefficients imexmrisr32;
		run(problem, &imexmrisr32, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("IMEXMRISR43",input_method_name) == 0) {
		IMEXMRISR43Coefficients imexmrisr43;
		run(problem, &imexmrisr43, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("IMEXMRIGARK21",input_method_name) == 0) {
		IMEXMRIGARK21Coefficients imexmrigark21;
		run(problem, &imexmrigark21, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else if(strcmp("IMEXMRIGARK32",input_method_name) == 0) {
		IMEXMRIGARK32Coefficients imexmrigark32;
		run(problem, &imexmrigark32, inner_coeffs, input_controller_name, tol_string, &err_norm);
	} else {
		printf("Error: Did not recognize method name: %s\n", input_method_name);
	}
}

void setup_and_run_with_problem(Problem* problem, const char* input_method_name, const char* input_inner_method_name, const char*  input_controller_name, const char* tol_string, double tol) {
	if(strcmp("SDIRK21",input_inner_method_name) == 0) {
		SDIRK21Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("ESDIRK324",input_inner_method_name) == 0) {
		ESDIRK324Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("ESDIRK325",input_inner_method_name) == 0) {
		ESDIRK325Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("ESDIRK436",input_inner_method_name) == 0) {
		ESDIRK436Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("ESDIRK437",input_inner_method_name) == 0) {
		ESDIRK437Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("ARK32ESDIRK",input_inner_method_name) == 0) {
		ARK32ESDIRKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("ARK43ESDIRK",input_inner_method_name) == 0) {
		ARK32ESDIRKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("HeunEulerERK",input_inner_method_name) == 0) {
		HeunEulerERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("BogackiShampineERK",input_inner_method_name) == 0) {
		BogackiShampineERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("ZonneveldERK",input_inner_method_name) == 0) {
		ZonneveldERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else if(strcmp("DormandPrinceERK",input_inner_method_name) == 0) {
		DormandPrinceERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, input_controller_name, tol_string, tol);
	} else {
		printf("Error: Did not recognize method name: %s\n", input_inner_method_name);
	}
}

int main(int argc, char* argv[]) {
	if(argc != 6) {
		printf("Error: Requires 5 command-line arguments. Ex: GenericAdaptiveDriver.exe <ProblemName> <MethodName> <InnerMethodName> <Controller> <tol>\n");
		return 1;
	} else {
		double tol = 0.0;
		const char* input_problem_name = argv[1];
		const char* input_method_name = argv[2];
		const char* input_inner_method_name = argv[3];
		const char* input_controller_name = argv[4];
		const char* tol_string = argv[5];
		sscanf(argv[5], "%lf", &tol);

		if(strcmp("Bicoupling",input_problem_name) == 0) {
			BicouplingProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("BicouplingDL",input_problem_name) == 0) {
			BicouplingDLProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("BicouplingLN",input_problem_name) == 0) {
			BicouplingLNProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Brusselator",input_problem_name) == 0) {
			BrusselatorProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("BrusselatorDL",input_problem_name) == 0) {
			BrusselatorDLProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("FourBody3d",input_problem_name) == 0) {
			FourBody3dProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Kaps",input_problem_name) == 0) {
			KapsProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("KapsDL",input_problem_name) == 0) {
			KapsDLProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("KapsLN",input_problem_name) == 0) {
			KapsLNProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("KPR",input_problem_name) == 0) {
			KPRProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("KPRDL",input_problem_name) == 0) {
			KPRProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Lienard",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("LienardDL",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("LienardLN",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Oregonator",input_problem_name) == 0) {
			OregonatorProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("OregonatorDL",input_problem_name) == 0) {
			OregonatorProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Pleiades",input_problem_name) == 0) {
			PleiadesProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Brusselator1D",input_problem_name) == 0) {
			Brusselator1DProblem problem(101);
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Bicoupling1D",input_problem_name) == 0) {
			Bicoupling1DProblem problem(101);
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Brusselator1DIMEX",input_problem_name) == 0) {
			Brusselator1DIMEXProblem problem(101);
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else if(strcmp("Bicoupling1DIMEX",input_problem_name) == 0) {
			Bicoupling1DIMEXProblem problem(101);
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, input_controller_name, tol_string, tol);
		} else {
			printf("Error: Did not recognize problem name: %s\n", input_problem_name);
			return 1;
		}
		return 0;
	}
}