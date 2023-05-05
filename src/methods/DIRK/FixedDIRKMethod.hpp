#ifndef FIXEDDIRKMETHOD_DEFINED__
#define FIXEDDIRKMETHOD_DEFINED__

#include "FixedStepSingleRateMethod.hpp"
#include "Residual.hpp"
#include "NewtonSolver.hpp"
#include "WeightedErrorNorm.hpp"
#include "DIRKResidual.hpp"

using namespace arma;
using namespace std;

class FixedDIRKMethod : FixedStepSingleRateMethod {
public:
	SingleRateMethodCoefficients* coeffs;
	Problem* problem;
	DIRKResidual dirk_residual;
	NewtonSolver newton_solver;
	struct NewtonSolverReturnValue newton_ret;
	WeightedErrorNorm* err_norm;
	double h;
	double effective_h;
	int total_output_points;
	int problem_dimension;
	vec explicit_data;
	vec y;
	vec y_temp;
	vec y_stage;
	mat y_stages;
	mat Y;
	bool is_top_level;
	
	FixedDIRKMethod(SingleRateMethodCoefficients* coeffs_, Problem* problem_, int problem_dimension_, WeightedErrorNorm* err_norm_, bool is_top_level_) :
	dirk_residual(coeffs_, problem_, problem_dimension_),
	newton_solver(&(FixedDIRKMethod::dirk_residual), 20, 1.0, problem_dimension_, err_norm_)
	{
		coeffs = coeffs_;
		problem = problem_;
		problem_dimension = problem_dimension_;
		err_norm = err_norm_;
		is_top_level = is_top_level_;
		
		declare_vectors();
	}

	mat solve(double t_0, double h_, vec* y_0, vec* output_tspan) {
		return solve(t_0, h_, y_0, output_tspan, true);
	}

	mat solve(double t_0, double h_, vec* y_0, vec* output_tspan, bool use_primary_coeffs) {
		prepare_solve(h_, output_tspan);
		y = *y_0;
		if (is_top_level) {
			err_norm->set_weights(y_0);
		}

		int output_index = 0;
		if((*output_tspan)(0) == t_0) {
			Y.col(0) = *y_0;
			output_index++;
		}

		double t = t_0;
		while(output_index < total_output_points) {
			if(t + h_ - (*output_tspan)(output_index) >= 1e-14) {
				effective_h = (*output_tspan)(output_index) - t;
			} else {
				effective_h = h_;
			}
			set_problem_dependent_data(effective_h);
			compute_stages(t, effective_h);
			compute_solution(effective_h, use_primary_coeffs);
			if (abs(t + effective_h - (*output_tspan)(output_index)) < 1e-14) {
				Y.col(output_index) = y;
				output_index++;
			}
			t += effective_h;
			if (is_top_level) {
				err_norm->set_weights(&y);
			}
		}
		return Y;
	}

	void compute_stages(double t, double h) {
		y_stages.zeros();
		for(int stage_idx=0; stage_idx<coeffs->num_stages; stage_idx++) {

			y_stage.zeros();
			y_temp.zeros();
			compute_explicit_data(stage_idx);
			//printf("Calculating stage\n");
			if ((coeffs->A(stage_idx,stage_idx)) != 0.0) {
				dirk_residual.set_function_dependent_data(stage_idx);
				newton_solver.solve(t, &y, &newton_ret);
				y_temp = newton_ret.y;
			} else {
				y_temp = y + h*explicit_data;
			}
			//printf("Finished calculating stage\n");
			problem->full_rhs(t+(coeffs->c(stage_idx))*h,&y_temp,&y_stage);
			y_stages.col(stage_idx) = y_stage;
		}
	}

	void compute_explicit_data(int stage_idx) {
		explicit_data.zeros();
		for(int inner_stage_idx=0; inner_stage_idx<stage_idx; inner_stage_idx++) {
			explicit_data += (coeffs->A(stage_idx,inner_stage_idx))*y_stages.col(inner_stage_idx);
		}
	}

	void compute_solution(double h, bool use_primary_coeffs) {
		y_temp.zeros();
		for(int stage_idx=0; stage_idx<coeffs->num_stages; stage_idx++) {
			if (use_primary_coeffs) {
				y_temp += (coeffs->b(stage_idx))*y_stages.col(stage_idx);
			} else {
				y_temp += (coeffs->d(stage_idx))*y_stages.col(stage_idx);
			}
		}
		y += h*y_temp;
	}

	void declare_vectors() {
		y_stages = mat(problem_dimension, coeffs->num_stages, fill::zeros);
		y = vec(problem_dimension, fill::zeros);
		y_temp = vec(problem_dimension, fill::zeros);
		y_stage = vec(problem_dimension, fill::zeros);
		explicit_data = vec(problem_dimension, fill::zeros);
		dirk_residual.set_explicit_data_pointer(&explicit_data);
		Y = mat(problem_dimension, 1, fill::zeros);
	}

	void prepare_solve(double h, vec* output_tspan) {
		effective_h = h;
		total_output_points = output_tspan->n_elem;

		y.zeros();
		y_stages.zeros();
		y_temp.zeros();
		y_stage.zeros();
		if (is_top_level) {
			Y = mat(problem_dimension, total_output_points, fill::zeros);
		}
		explicit_data.zeros();
	}

	void set_problem_dependent_data(double h_) {
		dirk_residual.set_problem_dependent_data(h_);
		newton_solver.set_problem_dependent_data(h_);
	}
};

#endif