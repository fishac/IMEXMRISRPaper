#ifndef MRIGARKFIXEDSTEPERRORMEASUREMENTTEST_DEFINED__
#define MRIGARKFIXEDSTEPERRORMEASUREMENTTEST_DEFINED__

#include "MRIGARKCoefficients.hpp"
#include "MRIGARKInnerRHSFunctions.hpp"
#include "MRIGARKExplicitRHS.hpp"
#include "MRIGARKExplicitRHSJacobian.hpp"
#include "MRIGARKImplicitResidual.hpp"
#include "MRIGARKImplicitResidualJacobian.hpp"
#include "SingleRateMethodCoefficients.hpp"
#include "FixedDIRKMethod.hpp"
#include "NewtonSolver.hpp"
#include "WeightedErrorNorm.hpp"
#include "FixedStepMultiRateStep.hpp"
#include "VernerERKCoefficients.hpp"
#include "FixedDIRKMethodFastMeasurement.hpp"

using namespace arma;

struct MRIGARKFixedStepErrorMeasurementTestReturnValue {
	vec y;
	double err_true;
	double err0;
	double err1;
	double err2;
	double err3;
	double err4;
	double err5;
	double err6;
};

class MRIGARKFixedStepErrorMeasurementTest : public FixedStepMultiRateStep {
public:
	MRIGARKCoefficients* coeffs;
	MRIGARKInnerRHSFunctions inner_rhs_funcs;
	MRIGARKExplicitRHS explicit_rhs;
	MRIGARKExplicitRHSJacobian explicit_rhs_jacobian;
	MRIGARKImplicitResidual implicit_residual;
	MRIGARKImplicitResidualJacobian implicit_residual_jacobian;
	NewtonSolver newton_solver;
	struct NewtonSolverReturnValue newton_ret;
	FixedDIRKMethod dirk;
	VernerERKCoefficients inner_coeffs_ref;
	FixedDIRKMethod dirk_reference;
	FixedDIRKMethodFastMeasurement dirk_fast_measurement;
	mat y_stages;
	mat y_stages_star;
	mat y_stages_fast;
	mat y_stages_reference;
	vec y;
	vec y_ref;
	vec y_star;
	vec v_0;
	vec v_H;
	double H;
	int M;
	int problem_dimension;
	int num_stages;

	MRIGARKFixedStepErrorMeasurementTest(MRIGARKCoefficients* coeffs_, SingleRateMethodCoefficients* inner_coeffs_, RHS* fast_func_, RHS* slow_func_, RHSJacobian* fast_func_jac_, RHSJacobian* slow_func_jac_, int problem_dimension_, WeightedErrorNorm* err_norm) :
	inner_rhs_funcs(coeffs_, fast_func_, slow_func_, fast_func_jac_, slow_func_jac_, problem_dimension_),
	implicit_residual(&(MRIGARKFixedStepErrorMeasurementTest::inner_rhs_funcs)),
	implicit_residual_jacobian(&(MRIGARKFixedStepErrorMeasurementTest::inner_rhs_funcs)),
	explicit_rhs(&(MRIGARKFixedStepErrorMeasurementTest::inner_rhs_funcs)),
	explicit_rhs_jacobian(&(MRIGARKFixedStepErrorMeasurementTest::inner_rhs_funcs)),
	dirk(inner_coeffs_, &(MRIGARKFixedStepErrorMeasurementTest::explicit_rhs), &(MRIGARKFixedStepErrorMeasurementTest::explicit_rhs_jacobian), problem_dimension_, err_norm),
	dirk_fast_measurement(inner_coeffs_, &(MRIGARKFixedStepErrorMeasurementTest::explicit_rhs), &(MRIGARKFixedStepErrorMeasurementTest::explicit_rhs_jacobian), problem_dimension_, err_norm),
	newton_solver(&(MRIGARKFixedStepErrorMeasurementTest::implicit_residual), &(MRIGARKFixedStepErrorMeasurementTest::implicit_residual_jacobian), 20, 1.0, problem_dimension_, err_norm),
	dirk_reference(&inner_coeffs_ref, &(MRIGARKFixedStepErrorMeasurementTest::explicit_rhs), &(MRIGARKFixedStepErrorMeasurementTest::explicit_rhs_jacobian), problem_dimension_, err_norm)
	{
		coeffs = coeffs_;
		problem_dimension = problem_dimension_;
		num_stages = coeffs_->num_stages;
		
		declare_vectors();
	}

	void step_solution(double t, double H, int M, vec* y_prev, MRIGARKFixedStepErrorMeasurementTestReturnValue* ret) {
		y_stages.col(0) = *y_prev;
		y_stages_star.col(0) = *y_prev;
		y_stages_fast.col(0) = *y_prev;
		y_stages_reference.col(0) = *y_prev;
		vec fast_errs(num_stages,fill::zeros);
		for(int stage_index=1; stage_index<num_stages; stage_index++) {
			inner_rhs_funcs.set_function_dependent_data(H, t, stage_index, false);
			

			double gbar = inner_rhs_funcs.gamma_bar(stage_index,stage_index);
			if(gbar != 0.0) {
				v_0 = y_stages.col(stage_index-1);
				implicit_step(t);
				y_stages.col(stage_index) = v_H;

				v_0 = y_stages_star.col(stage_index-1);
				implicit_step(t);
				y_stages_star.col(stage_index) = v_H;

				v_0 = y_stages_fast.col(stage_index-1);
				implicit_step(t);
				y_stages_fast.col(stage_index) = v_H;

				v_0 = y_stages_reference.col(stage_index-1);
				implicit_step(t);
				y_stages_reference.col(stage_index) = v_H;
			} else {
				//cout.precision(16);
				v_0 = y_stages.col(stage_index-1);
				explicit_step(H,M);
				y_stages.col(stage_index) = v_H;
				//v_H.raw_print(cout,"computed solution");

				v_0 = y_stages_star.col(stage_index-1);
				explicit_step_star(H,M);
				y_stages_star.col(stage_index) = v_H;

				v_0 = y_stages_fast.col(stage_index-1);
				double fast_err = explicit_step_fast(H,M);
				fast_errs(stage_index) = fast_err;
				y_stages_fast.col(stage_index) = v_H;

				v_0 = y_stages_reference.col(stage_index-1);
				explicit_step_reference(H,M);
				y_stages_reference.col(stage_index) = v_H;
				//v_H.raw_print(cout,"reference solution");

				//vec stage_err = y_stages.col(stage_index) - y_stages_reference.col(stage_index);
				//double stage_err_norm = norm(stage_err,2);
				//printf("stage_err_norm: %.16f, H: %.16f, t: %.16f\n",stage_err_norm,H,t);
			}
		}

		//cout.precision(16);
		//mat err_mat = y_stages - y_stages_reference;
		//err_mat.raw_print("err_mat");
		y = y_stages.col(num_stages-1);
		y_star = y_stages_star.col(num_stages-1);
		y_ref = y_stages_reference.col(num_stages-1);

		double err_true = norm(y-y_ref,2);
		double err0 = norm(y-y_star,2);
		double err1 = 0.0;
		double err2 = 0.0;
		double err3 = 0.0;
		double err4 = 0.0;
		double err5 = 0.0;
		double err6 = 0.0;

		mat temp1 = y_stages-y_stages_star;
		vec temp2(problem_dimension,fill::zeros);
		vec temp3(problem_dimension,fill::zeros);
		for (int i=0; i<num_stages; i++) {
			temp2 = temp1.col(i);
			temp3 += temp1.col(i);
			err3 = std::max(err3, norm(temp1,2));

			err6 = std::max(err6, fast_errs(i));
		}
		err1 = norm(temp3,2);
		err2 = norm(temp3,2)/num_stages;
		err4 = accu(fast_errs);
		err5 = accu(fast_errs)/num_stages;

		ret->y = y;
		ret->err_true = err_true;
		ret->err0 = err0;
		ret->err1 = err1;
		ret->err2 = err2;
		ret->err3 = err3;
		ret->err4 = err4;
		ret->err5 = err5;
		ret->err6 = err6;
	}

	void implicit_step(double t) {
		newton_solver.solve(t, &v_0, &newton_ret);
		problem->increment_slow_nonlinear_solves();
		v_H = newton_ret.y;
	}

	void explicit_step(double H, int M) {
		inner_rhs_funcs.explicit_set_previous_terms();
		vec output_tspan = {H};
		mat Y = dirk.solve(0.0, H/M, &v_0, &output_tspan);
		v_H = Y.col(0);
	}

	void explicit_step_star(double H, int M) {
		inner_rhs_funcs.explicit_set_previous_terms();
		vec output_tspan = {H};
		mat Y = dirk.solve(0.0, H/M, &v_0, &output_tspan, false);
		v_H = Y.col(0);
	}

	double explicit_step_fast(double H, int M) {
		inner_rhs_funcs.explicit_set_previous_terms();
		vec output_tspan = {H};
		FixedDIRKMethodFastMeasurementReturnValue ret;
		dirk_fast_measurement.solve(0.0, H/M, &v_0, &output_tspan, false, &ret);
		v_H = ret.y;
		return ret.err_estimate;
	}

	void explicit_step_reference(double H, int M) {
		inner_rhs_funcs.explicit_set_previous_terms();
		vec output_tspan = {H};
		mat Y = dirk_reference.solve(0.0, H/M, &v_0, &output_tspan);
		v_H = Y.col(0);
	}

	void declare_vectors() {
		v_0 = vec(problem_dimension, fill::zeros);
		v_H = vec(problem_dimension, fill::zeros);
		y_ref = vec(problem_dimension, fill::zeros);
		y_stages = mat(problem_dimension, num_stages, fill::zeros);
		y_stages_star = mat(problem_dimension, num_stages, fill::zeros);
		y_stages_reference = mat(problem_dimension, num_stages, fill::zeros);
		y_stages_fast = mat(problem_dimension, num_stages, fill::zeros);
		inner_rhs_funcs.set_problem_dependent_data(&y_stages);
	}
};

#endif