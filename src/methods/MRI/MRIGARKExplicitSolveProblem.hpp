#ifndef MRIGARKEXPLICITSOLVEPROBLEM_DEFINED__
#define MRIGARKEXPLICITSOLVEPROBLEM_DEFINED__

#include "MRIGARKInnerRHSFunctions.hpp"

using namespace arma;

class MRIGARKExplicitSolveProblem: public Problem {
public:
	MRIGARKInnerRHSFunctions* inner_rhs_funcs;
	MRIGARKExplicitSolveProblem(MRIGARKInnerRHSFunctions* inner_rhs_funcs_) {
		name = "MRIGARKExplicitSolve";
		inner_rhs_funcs = inner_rhs_funcs_;
	}
	
	void full_rhs_custom(double theta, vec* v, vec* f) {
		inner_rhs_funcs->explicit_solve_rhs(theta, v, f);
	}
	
	void full_rhsjacobian_custom(double theta, vec* v, mat* jac) {
		inner_rhs_funcs->explicit_solve_rhsjacobian(theta, v, jac);
	}
};

#endif