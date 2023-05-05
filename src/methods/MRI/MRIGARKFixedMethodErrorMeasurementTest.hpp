#ifndef MRIGARKFIXEDMETHODERRORMEASUREMENTTEST_DEFINED__
#define MRIGARKFIXEDMETHODERRORMEASUREMENTTEST_DEFINED__

#include "FixedStepMultiRateMethod.hpp"
#include "MRIGARKCoefficients.hpp"
#include "SingleRateMethodCoefficients.hpp"
#include "MRIGARKFixedStep.hpp"
#include "WeightedErrorNorm.hpp"
#include "RHS.hpp"
#include "RHSJacobian.hpp"
#include "MRIGARKFixedStepErrorMeasurementTest.hpp"

using namespace arma;

struct MRIGARKFixedMethodErrorMeasurementTestReturnValue {
	mat Y;
	std::vector<double> err_true;
	std::vector<double> err0;
	std::vector<double> err1;
	std::vector<double> err2;
	std::vector<double> err3;
	std::vector<double> err4;
	std::vector<double> err5;
	std::vector<double> err6;
};

class MRIGARKFixedMethodErrorMeasurementTest : public FixedStepMultiRateMethod {
public:
	vec y;
	mat Y;
	double effective_H;
	int total_output_points;
	int problem_dimension;
	MRIGARKFixedStepErrorMeasurementTestReturnValue step_ret;

	MRIGARKFixedMethodErrorMeasurementTest(int problem_dimension_) {
		problem_dimension = problem_dimension_;
	}

	void solve(double t_0, double H, int M, vec* y_0, vec* output_tspan, MRIGARKFixedStepErrorMeasurementTest* mrigark_step, MRIGARKFixedMethodErrorMeasurementTestReturnValue* ret) {
		prepare_solve(H, output_tspan);
		y = *y_0;

		std::vector<double> err_true;
		std::vector<double> err0;
		std::vector<double> err1;
		std::vector<double> err2;
		std::vector<double> err3;
		std::vector<double> err4;
		std::vector<double> err5;
		std::vector<double> err6;

		int output_index = 0;
		if((*output_tspan)(0) == t_0) {
			Y.col(0) = *y_0;
			output_index++;
		}

		double t = t_0;
		while(output_index < total_output_points) {
			if(t + H - (*output_tspan)(output_index) > 0.0) {
				effective_H = (*output_tspan)(output_index) - t;
			} else {
				effective_H = H;
			}
			mrigark_step->step_solution(t, effective_H, M, &y, &step_ret);
			y = step_ret.y;
			err_true.push_back(step_ret.err_true);
			err0.push_back(step_ret.err0);
			err1.push_back(step_ret.err1);
			err2.push_back(step_ret.err2);
			err3.push_back(step_ret.err3);
			err4.push_back(step_ret.err4);
			err5.push_back(step_ret.err5);
			err6.push_back(step_ret.err6);

			if (t + effective_H == (*output_tspan)(output_index)) {
				Y.col(output_index) = y;
				output_index++;
			}
			t += effective_H;
		}
		ret->Y = Y;
		ret->err_true = err_true;
		ret->err0 = err0;
		ret->err1 = err1;
		ret->err2 = err2;
		ret->err3 = err3;
		ret->err4 = err4;
		ret->err5 = err5;
		ret->err6 = err6;
	}

	void prepare_solve(double H_, vec* output_tspan) {
		effective_H = H_;
		total_output_points = output_tspan->n_elem;
		y = vec(problem_dimension, fill::zeros);
		Y = mat(problem_dimension, total_output_points, fill::zeros);
	}
};

#endif