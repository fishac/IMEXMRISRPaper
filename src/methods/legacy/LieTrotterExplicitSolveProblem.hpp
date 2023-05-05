#ifndef LIETROTTEREXPLICITSOLVEPROBLEM_DEFINED__
#define LIETROTTEREXPLICITSOLVEPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace arma;

class LieTrotterExplicitSolveProblem: public Problem {
public:
	Problem* problem;
	double t;
	
	LieTrotterExplicitSolveProblem(Problem* problem_) {
		name = "LieTrotterExplicitSolve";
		problem = problem_;
	}
	
	void full_rhs_custom(double theta, vec* v, vec* f) {
		problem->fast_rhs(t+theta, v, f);
	}
	
	void full_rhsjacobian_custom(double theta, vec* v, mat* jac) {
		problem->fast_rhsjacobian(t+theta, v, jac);
	}
	
	void set_step_dependent_data(double t_) {
		t = t_;
	}
};

#endif