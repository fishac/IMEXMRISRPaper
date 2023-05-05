#ifndef FIXEDDIRKMETHODFASTMEASUREMENT_DEFINED__
#define FIXEDDIRKMETHODFASTMEASUREMENT_DEFINED__

#include "FixedStepSingleRateMethod.hpp"
#include "Residual.hpp"
#include "NewtonSolver.hpp"
#include "WeightedErrorNorm.hpp"
#include "DIRKResidual.hpp"

using namespace arma;
using namespace std;

struct FixedDIRKMethodFastMeasurementReturnValue {
	vec y;
	double err_estimate;
	int status;
};

class FixedDIRKMethodFastMeasurement {
public:
	SingleRateMethodCoefficients* coeffs;
	Problem* problem;
	DIRKResidual dirk_residual;
	NewtonSolver newton_solver;
	WeightedErrorNorm* err_norm;
	struct NewtonSolverReturnValue newton_ret;
	double h;
	double effective_h;
	int total_output_points;
	int problem_dimension;
	vec explicit_data;
	vec y;
	vec y2;
	vec y_temp;
	vec y_temp2;
	vec y_stage;
	mat y_stages;
	mat Y;
	double err_estimate;

	FixedDIRKMethodFastMeasurement(SingleRateMethodCoefficients* coeffs_, Problem* problem_, int problem_dimension_, WeightedErrorNorm* err_norm_) :
	dirk_residual(coeffs_, problem_, problem_dimension_),
	newton_solver(&(FixedDIRKMethodFastMeasurement::dirk_residual), 20, 1.0, problem_dimension_, err_norm_)
	{
		coeffs = coeffs_;
		problem = problem_;
		problem_dimension = problem_dimension_;
		err_norm = err_norm_;
		
		declare_vectors();
	}

	void solve(double t_0, double h_, vec* y_0, vec* output_tspan, FixedDIRKMethodFastMeasurementReturnValue* ret) {
		solve(t_0, h_, y_0, output_tspan, true, ret);
	}

	void solve(double t_0, double h_, vec* y_0, vec* output_tspan, bool useRelNorm, FixedDIRKMethodFastMeasurementReturnValue* ret) {
		prepare_solve(h_, output_tspan);
		y = *y_0;

		int output_index = 0;
		if((*output_tspan)(0) == t_0) {
			Y.col(0) = *y_0;
			output_index++;
		}

		double t = t_0;
		while(output_index < total_output_points) {
			if(t + h_ - (*output_tspan)(output_index) > 0.0) {
				effective_h = (*output_tspan)(output_index) - t;
			} else {
				effective_h = h_;
			}
			set_problem_dependent_data(effective_h);
			compute_stages(t, effective_h);
			compute_solution(effective_h, useRelNorm);
			if (t + effective_h == (*output_tspan)(output_index)) {
				Y.col(output_index) = y;
				output_index++;
			}
			t += effective_h;
		}
		ret->y = y;
		ret->err_estimate = err_estimate;
		ret->status = 0;
	}

	void compute_stages(double t, double h) {
		for(int stage_idx=0; stage_idx<coeffs->num_stages; stage_idx++) {
			dirk_residual.set_function_dependent_data(stage_idx);

			y_stage.zeros();
			y_temp.zeros();
			compute_explicit_data(stage_idx);
			//printf("Calculating stage\n");
			if ((coeffs->A(stage_idx,stage_idx)) != 0.0) {
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

	void compute_solution(double h, bool useRelNorm) {
		y_temp.zeros();
		y_temp2.zeros();
		for(int stage_idx=0; stage_idx<coeffs->num_stages; stage_idx++) {
			y_temp += (coeffs->b(stage_idx))*y_stages.col(stage_idx);
			y_temp2 += (coeffs->d(stage_idx))*y_stages.col(stage_idx);
		}

		y2 = y + h*y_temp2;
		y += h*y_temp;
		if (useRelNorm) {
			err_estimate += err_norm->compute_norm(y-y2);
		} else {
			err_estimate += norm(y-y2,2);
		}
	}

	void declare_vectors() {
		y_stages = mat(problem_dimension, coeffs->num_stages, fill::zeros);
		y = vec(problem_dimension, fill::zeros);
		y2 = vec(problem_dimension, fill::zeros);
		y_temp = vec(problem_dimension, fill::zeros);
		y_temp2 = vec(problem_dimension, fill::zeros);
		y_stage = vec(problem_dimension, fill::zeros);
		explicit_data = vec(problem_dimension, fill::zeros);
		dirk_residual.set_explicit_data_pointer(&explicit_data);
	}

	void prepare_solve(double h, vec* output_tspan) {
		err_estimate = 0.0;
		effective_h = h;
		total_output_points = output_tspan->n_elem;

		y.zeros();
		y2.zeros();
		y_stages.zeros();
		y_temp.zeros();
		y_temp2.zeros();
		y_stage.zeros();
		Y = mat(problem_dimension, total_output_points, fill::zeros);
		explicit_data.zeros();
	}

	void set_problem_dependent_data(double h_) {
		dirk_residual.set_problem_dependent_data(h_);
		newton_solver.set_problem_dependent_data(h_);
	}
};

#endif