#ifndef MRIGARKFIXEDMETHOD_DEFINED__
#define MRIGARKFIXEDMETHOD_DEFINED__

#include "FixedStepMultiRateMethod.hpp"
#include "SingleRateMethodCoefficients.hpp"
#include "MRIGARKFixedStep.hpp"
#include "WeightedErrorNorm.hpp"
#include "FixedStepMultiRateStep.hpp"

using namespace arma;

class MRIGARKFixedMethod : public FixedStepMultiRateMethod {
public:
	vec y;
	mat Y;
	double effective_H;
	int total_output_points;
	int problem_dimension;
	Problem* problem;

	MRIGARKFixedMethod(Problem* problem_, int problem_dimension_) {
		problem = problem_;
		problem_dimension = problem_dimension_;
	}

	mat solve(double t_0, double H, int M, vec* y_0, vec* output_tspan, FixedStepMultiRateStep* mrigark_step) {
		prepare_solve(H, output_tspan);
		y = *y_0;
		problem->set_yhat_that(*y_0,t_0);

		int output_index = 0;
		if((*output_tspan)(0) == t_0) {
			Y.col(0) = *y_0;
			output_index++;
		}

		double t = t_0;
		while(output_index < total_output_points) {
			if(t + H - (*output_tspan)(output_index) >= 1e-14) {
				effective_H = (*output_tspan)(output_index) - t;
			} else {
				effective_H = H;
			}
			mrigark_step->step_solution(t, effective_H, M, &y, &y);
			if (abs(t + effective_H - (*output_tspan)(output_index)) < 1e-14) {
				Y.col(output_index) = y;
				output_index++;
			}
			t += effective_H;
			problem->set_yhat_that(y,t);
		}
		return Y;
	}

	void prepare_solve(double H_, vec* output_tspan) {
		effective_H = H_;
		total_output_points = output_tspan->n_elem;
		y = vec(problem_dimension, fill::zeros);
		Y = mat(problem_dimension, total_output_points, fill::zeros);
	}
};

#endif