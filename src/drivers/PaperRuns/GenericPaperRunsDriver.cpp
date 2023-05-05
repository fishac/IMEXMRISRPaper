#include <armadillo>
#include <math.h>
#include <stdio.h>

#include "BicouplingProblem.hpp"
#include "BicouplingDLProblem.hpp"
#include "BicouplingLNProblem.hpp"
#include "BrusselatorProblem.hpp"
#include "BrusselatorDLProblem.hpp"
#include "BrusselatorPDEProblem.hpp"
#include "BrusselatorPDE201Problem.hpp"
#include "BrusselatorPDE801Problem.hpp"
#include "FourBody3dProblem.hpp"
#include "KapsProblem.hpp"
#include "KapsDLProblem.hpp"
#include "KapsLNProblem.hpp"
#include "KPRProblem.hpp"
#include "KPRSlowOnlyProblem.hpp"
#include "KPRFastOnlyProblem.hpp"
#include "KPRDLProblem.hpp"
#include "LienardProblem.hpp"
#include "LienardDLProblem.hpp"
#include "LienardLNProblem.hpp"
#include "OregonatorProblem.hpp"
#include "OregonatorDLProblem.hpp"
#include "PleiadesProblem.hpp"
#include "DecoupledLinearProblem.hpp"
#include "Lorenz63Problem.hpp"
#include "NonstiffBrusselatorProblem.hpp"
#include "Problem.hpp"
#include "PaperRunsDriver.hpp"

using namespace std;
using namespace arma;

mat load_true_sol(const char* problem_name) {
	char filename[100];
	sprintf(filename,"./resources/%s/%s_fixed_truesol_11.csv",problem_name,problem_name);
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

vec generate_H_vec(Problem* problem, int two_exp_steps_min, int two_exp_steps_max) {
	int total_Hs = two_exp_steps_max - two_exp_steps_min + 1;
	vec H_vec(total_Hs, fill::zeros);
	for(int ih=0; ih<total_Hs; ih++) {
		H_vec(ih) = (problem->time_step_multiplier)*pow(2.0, -(ih+two_exp_steps_min));
	}
	return H_vec;
}

void setup_and_run(Problem* problem, const char* input_method_name, SingleRateMethodCoefficients* inner_coeffs, int two_exp_steps_min, int two_exp_steps_max, int M) {
	vec H_vec = generate_H_vec(problem, two_exp_steps_min, two_exp_steps_max);
	
	vec output_tspan = linspace(problem->t_0, problem->t_f, 11);
	mat Y_true = get_true_sol(&output_tspan, problem);

	PaperRunsDriver driver;
	
	if(strcmp("MRIGARKIRK21a",input_method_name) == 0) {
		MRIGARKIRK21aCoefficients mrigarkirk21a;
		driver.run(problem, &mrigarkirk21a, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("MRIGARKERK33",input_method_name) == 0) {
		MRIGARKERK33Coefficients mrigarkerk33;
		driver.run(problem, &mrigarkerk33, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("MRIGARKERK45a",input_method_name) == 0) {
		MRIGARKERK45aCoefficients mrigarkerk45a;
		driver.run(problem, &mrigarkerk45a, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("MRIGARKESDIRK34a",input_method_name) == 0) {
		MRIGARKESDIRK34aCoefficients mrigarkesdirk34a;
		driver.run(problem, &mrigarkesdirk34a, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("MRIGARKESDIRK46a",input_method_name) == 0) {
		MRIGARKESDIRK46aCoefficients mrigarkesdirk46a;
		driver.run(problem, &mrigarkesdirk46a, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("MRIGARKERK22a",input_method_name) == 0) {
		MRIGARKERK22aCoefficients mrigarkerk22a;
		driver.run(problem, &mrigarkerk22a, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRITest1",input_method_name) == 0) {
		IMEXMRITest1Coefficients imexmritest1;
		driver.run(problem, &imexmritest1, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISR21",input_method_name) == 0) {
		IMEXMRISR21Coefficients imexmrisr21;
		driver.run(problem, &imexmrisr21, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISR32",input_method_name) == 0) {
		IMEXMRISR32Coefficients imexmrisr32;
		driver.run(problem, &imexmrisr32, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISR43",input_method_name) == 0) {
		IMEXMRISR43Coefficients imexmrisr43;
		driver.run(problem, &imexmrisr43, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISRTEST4",input_method_name) == 0) {
		IMEXMRISRTEST4Coefficients imexmrisrtest4;
		driver.run(problem, &imexmrisrtest4, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRI4",input_method_name) == 0) {
		IMEXMRI4Coefficients imexmri4;
		driver.run(problem, &imexmri4, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRI3a",input_method_name) == 0) {
		IMEXMRI3aCoefficients imexmri3a;
		driver.run(problem, &imexmri3a, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRI3b",input_method_name) == 0) {
		IMEXMRI3bCoefficients imexmri3b;
		driver.run(problem, &imexmri3b, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISRMERK2",input_method_name) == 0) {
		IMEXMRISRMERK2Coefficients imexmrisrmerk2;
		driver.run(problem, &imexmrisrmerk2, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISRMERK3",input_method_name) == 0) {
		IMEXMRISRMERK3Coefficients imexmrisrmerk3;
		driver.run(problem, &imexmrisrmerk3, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISRMERK4",input_method_name) == 0) {
		IMEXMRISRMERK4Coefficients imexmrisrmerk4;
		driver.run(problem, &imexmrisrmerk4, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRISRMERK5",input_method_name) == 0) {
		IMEXMRISRMERK5Coefficients imexmrisrmerk5;
		driver.run(problem, &imexmrisrmerk5, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if (strcmp("LieTrotter", input_method_name) == 0) {
		driver.run_lietrotter_method(problem, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if (strcmp("StrangMarchuk", input_method_name) == 0) {
		driver.run_strangmarchuk_method(problem, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRIGARK21",input_method_name) == 0) {
		IMEXMRIGARK21Coefficients imexmrigark21;
		driver.run(problem, &imexmrigark21, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	} else if(strcmp("IMEXMRIGARK32",input_method_name) == 0) {
		IMEXMRIGARK32Coefficients imexmrigark32;
		driver.run(problem, &imexmrigark32, inner_coeffs, &H_vec, M, &output_tspan, &Y_true);
	}  else {
		printf("Error: Did not recognize method name: %s\n", input_method_name);
	}
}

void setup_and_run_with_problem(Problem* problem, const char* input_method_name, const char* input_inner_method_name, int two_exp_steps_min, int two_exp_steps_max, int M) {
	if(strcmp("SDIRK21",input_inner_method_name) == 0) {
		SDIRK21Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("ESDIRK324",input_inner_method_name) == 0) {
		ESDIRK324Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("ESDIRK325",input_inner_method_name) == 0) {
		ESDIRK325Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("ESDIRK436",input_inner_method_name) == 0) {
		ESDIRK436Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("ESDIRK437",input_inner_method_name) == 0) {
		ESDIRK437Coefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("ARK32ESDIRK",input_inner_method_name) == 0) {
		ARK32ESDIRKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("ARK43ESDIRK",input_inner_method_name) == 0) {
		ARK32ESDIRKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("HeunEulerERK",input_inner_method_name) == 0) {
		HeunEulerERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("BogackiShampineERK",input_inner_method_name) == 0) {
		BogackiShampineERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("ZonneveldERK",input_inner_method_name) == 0) {
		ZonneveldERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else if(strcmp("DormandPrinceERK",input_inner_method_name) == 0) {
		DormandPrinceERKCoefficients inner_coeffs;
		setup_and_run(problem, input_method_name, &inner_coeffs, two_exp_steps_min, two_exp_steps_max, M);
	} else {
		printf("Error: Did not recognize method name: %s\n", input_inner_method_name);
	}
}


int main(int argc, char* argv[]) {
	if(argc != 7) {
		printf("Error: Requires 6 command-line arguments. Ex: GenericPaperRunsDriver.exe <ProblemName:string> <MethodName:string> <InnerMethodName:string> <TwoExpStepsMin:int> <TwoExpStepsMax:int> <M:int>\n");
		printf("Ex: GenericPaperRunsDriver.exe KPR IMEXMRISR21 SDIRK21 5 12 10\n");
		return 1;
	} else {
		const char* input_problem_name = argv[1];
		const char* input_method_name = argv[2];
		const char* input_inner_method_name = argv[3];
		const char* input_two_exp_steps_min = argv[4];
		const char* input_two_exp_steps_max = argv[5];
		const char* input_m = argv[6];
		
		int two_exp_steps_min = 0;
		int two_exp_steps_max = 0;
		int M = 0;
		
		sscanf(input_two_exp_steps_min, "%d", &two_exp_steps_min);
		sscanf(input_two_exp_steps_max, "%d", &two_exp_steps_max);
		sscanf(input_m, "%d", &M);

		if(strcmp("Bicoupling",input_problem_name) == 0) {
			BicouplingProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("BicouplingDL",input_problem_name) == 0) {
			BicouplingDLProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("BicouplingLN",input_problem_name) == 0) {
			BicouplingLNProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("Brusselator",input_problem_name) == 0) {
			BrusselatorProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("BrusselatorDL",input_problem_name) == 0) {
			BrusselatorDLProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("FourBody3d",input_problem_name) == 0) {
			FourBody3dProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("Kaps",input_problem_name) == 0) {
			KapsProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("KapsDL",input_problem_name) == 0) {
			KapsDLProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("KapsLN",input_problem_name) == 0) {
			KapsLNProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("KPR",input_problem_name) == 0) {
			KPRProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("KPRFastOnly",input_problem_name) == 0) {
			KPRFastOnlyProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("KPRSlowOnly",input_problem_name) == 0) {
			KPRSlowOnlyProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("KPRDL",input_problem_name) == 0) {
			KPRProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("Lienard",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("LienardDL",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("LienardLN",input_problem_name) == 0) {
			LienardProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("Oregonator",input_problem_name) == 0) {
			OregonatorProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("OregonatorDL",input_problem_name) == 0) {
			OregonatorProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("Pleiades",input_problem_name) == 0) {
			PleiadesProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("DecoupledLinear",input_problem_name) == 0) {
			DecoupledLinearProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("Lorenz63",input_problem_name) == 0) {
			Lorenz63Problem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("NonstiffBrusselator",input_problem_name) == 0) {
			NonstiffBrusselatorProblem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("BrusselatorPDE",input_problem_name) == 0) {
			BrusselatorPDEProblem problem(201);
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("BrusselatorPDE201",input_problem_name) == 0) {
			BrusselatorPDE201Problem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else if(strcmp("BrusselatorPDE801",input_problem_name) == 0) {
			BrusselatorPDE801Problem problem;
			setup_and_run_with_problem(&problem, input_method_name, input_inner_method_name, two_exp_steps_min, two_exp_steps_max, M);
		} else {
			printf("Error: Did not recognize problem name: %s\n", input_problem_name);
			return 1;
		}
		return 0;
	}
}