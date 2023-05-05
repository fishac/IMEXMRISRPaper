#ifndef MRIGARKINNERRHSFUNCTIONS_DEFINED__
#define MRIGARKINNERRHSFUNCTIONS_DEFINED__

#include "MRICoefficients.hpp"
#include "Problem.hpp"

using namespace arma;

class MRIGARKInnerRHSFunctions {
public:
	Problem* problem;
	
	vec* c;
	std::vector<mat>* gammas;
	std::vector<mat>* omegas;

	int method_type;
	int problem_dimension;
	int num_stages;
	int num_gammas;
	int num_omegas;
	double H;
	double t;
	int stage_index;
	bool embedding;
	int embedding_shift;

	vec f_temp;
	vec f_temp2;
	vec implicit_previous_terms;
	mat mri_eval_stages1;
	mat mri_eval_stages2;
	vec y_temp;
	mat jac_temp;
	mat jac_temp2;
	mat I;

	MRIGARKInnerRHSFunctions(MRICoefficients* coeffs, Problem* problem_, int problem_dimension_) {
		c = &(coeffs->c);
		gammas = &(coeffs->gammas);
		omegas = &(coeffs->omegas);
		problem = problem_;
		problem_dimension = problem_dimension_;
		
		num_stages = coeffs->num_stages;
		num_gammas = coeffs->num_gammas;
		num_omegas = coeffs->num_omegas;
		method_type = coeffs->method_type;

		f_temp = vec(problem_dimension, fill::zeros);
		f_temp2 = vec(problem_dimension, fill::zeros);
		mri_eval_stages1 = mat(problem_dimension, coeffs->num_stages, fill::zeros);
		mri_eval_stages2 = mat(problem_dimension, coeffs->num_stages, fill::zeros);
		implicit_previous_terms = vec(problem_dimension, fill::zeros);
		y_temp = vec(problem_dimension, fill::zeros);
		jac_temp = mat(problem_dimension, problem_dimension, fill::zeros);
		jac_temp2 = mat(problem_dimension, problem_dimension, fill::zeros);
		I = eye(problem_dimension, problem_dimension);
	}

	void explicit_solve_rhs(double theta, vec* v, vec* f) {		
		y_temp.zeros();
		f_temp.zeros();
		f->zeros();

		if (method_type == 0) {
			if (num_omegas == 0) {
				explicit_solve_rhs_mrigark(theta,v,f);
			} else {
				explicit_solve_rhs_imexmrigark(theta,v,f);
			}
		} else if (method_type == 1) {
			explicit_solve_rhs_mrexp(theta,v,f);
		} else if (method_type == 2) {
			explicit_solve_rhs_imexmrisr2(theta,v,f);
		}
	}
	
	void explicit_solve_rhs_mrigark(double theta, vec* v, vec* f) {
		double delta_c = (*c)(stage_index) - (*c)(stage_index-1); 
		double T_prev = t + (*c)(stage_index-1)*H;
		double tau = (theta-T_prev)/(delta_c*H);
		
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			(*f) += gamma(stage_index + embedding_shift, sub_stage_index, tau)*f_temp;
		}
		(*f) /= delta_c;
		problem->fast_rhs(theta, v, &f_temp);
		(*f) += f_temp;
	}
	
	void explicit_solve_rhs_imexmrigark(double theta, vec* v, vec* f) {
		double delta_c = (*c)(stage_index) - (*c)(stage_index-1); 
		double T_prev = t + (*c)(stage_index-1)*H;
		double tau = (theta-T_prev)/(delta_c*H);
		
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			(*f) += gamma(stage_index + embedding_shift, sub_stage_index, tau)*f_temp;
			
			f_temp = mri_eval_stages2.col(sub_stage_index);
			(*f) += omega(stage_index + embedding_shift, sub_stage_index, tau)*f_temp;
		}
		(*f) /= delta_c;
		problem->fast_rhs(theta, v, &f_temp);
		(*f) += f_temp;
	}
	
	void explicit_solve_rhs_mrexp(double theta, vec* v, vec* f) {
		problem->linear_rhs(t+theta,v,f);
		
		(*f) += mri_eval_stages1.col(0);
		
		if (stage_index > 1) {
			f_temp.zeros();
			for(int sub_stage_index=0; sub_stage_index<stage_index-1; sub_stage_index++) {
				f_temp = mri_eval_stages1.col(sub_stage_index+1);
				(*f) += gamma_mrexp(stage_index + embedding_shift,sub_stage_index,theta/H)*f_temp;
			}
		}
	}
	
	void explicit_solve_rhs_mrigarksr(double theta, vec* v, vec* f) {
		double tau = (theta-t)/((*c)(stage_index)*H);
		
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			(*f) += gamma(stage_index + embedding_shift, sub_stage_index, tau)*f_temp;
		}
		
		(*f) /= (*c)(stage_index);
		problem->fast_rhs(theta, v, &f_temp);
		(*f) += f_temp;
	}
	
	void explicit_solve_rhs_imexmrisr(double theta, vec* v, vec* f) {
		double tau = (theta-t)/((*c)(stage_index)*H);
		
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages2.col(sub_stage_index);
			(*f) += omega(stage_index + embedding_shift, sub_stage_index, tau)*f_temp;
		}
		
		(*f) /= (*c)(stage_index);
		problem->fast_rhs(theta, v, &f_temp);
		(*f) += f_temp;
	}
	
	void explicit_solve_rhs_imexmrisr2(double theta, vec* v, vec* f) {
		double tau = (theta-t)/((*c)(stage_index)*H);
		
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			(*f) += omega(stage_index + embedding_shift, sub_stage_index, tau)*f_temp;
			
			f_temp = mri_eval_stages2.col(sub_stage_index);
			(*f) += omega(stage_index + embedding_shift, sub_stage_index, tau)*f_temp;
		}
		
		(*f) /= (*c)(stage_index);
		problem->fast_rhs(theta, v, &f_temp);
		(*f) += f_temp;
	}

	void explicit_solve_rhsjacobian(double theta, vec* v, mat* jac) {
		if (method_type == 0) {
			problem->fast_rhsjacobian(theta, v, jac);
		} else if (method_type == 1) {
			problem->linear_rhsjacobian(theta, v, jac);
		} else if (method_type == 2) {
			problem->fast_rhsjacobian(theta, v, jac);
		}
	}

	void store_stage_func_eval(vec* y_stage, int stage_index) {
		f_temp.zeros();

		if (method_type == 0) {
			if (num_omegas == 0) {
				store_stage_func_eval_mrigark(y_stage, stage_index);
			} else {
				store_stage_func_eval_imexmrigark(y_stage, stage_index);
			}
		} else if (method_type == 1) {
			store_stage_func_eval_mrexp(y_stage, stage_index);
		} else if (method_type == 2) {
			store_stage_func_eval_imexmrigark(y_stage, stage_index);
		}
	}
	
	void store_stage_func_eval_mrigark(vec* y_stage, int stage_index) {
		double tj = t + (*c)(stage_index)*H;
		
		// Set f_temp = f^S(t+c_j*H,Y_j)
		problem->slow_rhs(tj, y_stage, &f_temp);
		mri_eval_stages1.col(stage_index) = f_temp;
	}
	
	void store_stage_func_eval_imexmrigark(vec* y_stage, int stage_index) {
		double tj = t + (*c)(stage_index)*H;
		
		// Set f_temp = f^I(t+c_j*H,Y_j)
		problem->implicit_rhs(tj, y_stage, &f_temp);
		mri_eval_stages1.col(stage_index) = f_temp;
		
		// Set f_temp = f^E(t+c_j*H,Y_j)
		problem->explicit_rhs(tj, y_stage, &f_temp);
		mri_eval_stages2.col(stage_index) = f_temp;
	}
	
	void store_stage_func_eval_mrexp(vec* y_stage, int stage_index) {
		if (stage_index == 0) {
			// Set col = f^N(t,Y_0)
			problem->nonlinear_rhs(t, y_stage, &f_temp);
			mri_eval_stages1.col(stage_index) = f_temp;
		} else {
			// Set col = f^N(t+c_j*H,Y_j) - f^N(t,Y_0)
			double tj = t + (*c)(stage_index)*H;
			problem->nonlinear_rhs(tj, y_stage, &f_temp2);
			mri_eval_stages1.col(stage_index) = f_temp2 - mri_eval_stages1.col(0);
		}
	}
	
	void reset_stage_func_eval_storage() {
		mri_eval_stages1.zeros();
		mri_eval_stages2.zeros();
	}

	void implicit_set_previous_terms() {
		implicit_previous_terms.zeros();
		
		if (method_type == 0) {
			if (num_omegas == 0) {
				implicit_set_previous_terms_mrigark();
			} else {
				implicit_set_previous_terms_imexmrigark();
			}
		} else if (method_type == 2) {
			implicit_set_previous_terms_imexmrisr2();
		}
	}
	
	void implicit_set_previous_terms_mrigark() {
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			f_temp *= gamma_bar(stage_index + embedding_shift, sub_stage_index);
			implicit_previous_terms += f_temp;
		}
	}
	
	void implicit_set_previous_terms_imexmrigark() {
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			implicit_previous_terms += gamma_bar(stage_index + embedding_shift, sub_stage_index)*f_temp;
			
			f_temp = mri_eval_stages2.col(sub_stage_index);
			implicit_previous_terms += omega_bar(stage_index + embedding_shift, sub_stage_index)*f_temp;
		}
	}
	
	void implicit_set_previous_terms_imexmrisr() {
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			implicit_previous_terms += (gammas->at(0))(stage_index + embedding_shift, sub_stage_index)*f_temp;
		}
	}
	
	void implicit_set_previous_terms_imexmrisr2() {
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp = mri_eval_stages1.col(sub_stage_index);
			implicit_previous_terms += (gammas->at(0))(stage_index + embedding_shift, sub_stage_index)*f_temp;
		}
	}

	vec implicit_get_previous_terms() {
		return implicit_previous_terms;
	}
	
	void implicit_solve_residual(vec* explicit_data, vec* y_prev, vec* y, vec* f) {
		implicit_solve_current_term(y, &f_temp);

		// Compute result vector
		*f = *y - *y_prev - H*(*explicit_data + f_temp);
	}

	void implicit_solve_current_term(vec* y, vec* f) {
		double ti = t + (*c)(stage_index)*H;

		if (num_omegas == 0) {
			problem->slow_rhs(ti, y, &f_temp2);
		} else {
			problem->implicit_rhs(ti, y, &f_temp2);
		}
		if (method_type == 0) {
			f_temp2 *= gamma_bar(stage_index + embedding_shift, stage_index);
		} else if (method_type == 1) {
			// Do nothing, MERK cannot be implicit
		} else if (method_type == 2) {
			f_temp2 *= (gammas->at(0))(stage_index + embedding_shift, stage_index);
		}

		(*f) = f_temp2;
	}
	
	void implicit_solve_residual_jacobian(vec* y, mat* jac) {
		if (method_type == 0) {
			implicit_solve_residual_jacobian_mrigark(y, jac);
		} else if (method_type == 1) {
			// Do nothing. MrExp methods are not implicit.
		} else if (method_type == 2) {
			implicit_solve_residual_jacobian_imexmrisr(y, jac);
		}
	}

	void implicit_solve_residual_jacobian_mrigark(vec* y, mat* jac) {
		double gbar = gamma_bar(stage_index + embedding_shift, stage_index);
		if (gbar <= 2.0*1e-8) {
			(*jac) = I;
		} else {
			implicit_solve_jacobian(y, &jac_temp);
			(*jac) = I - H*gbar*jac_temp;
		}
	}
	
	void implicit_solve_residual_jacobian_imexmrisr(vec* y, mat* jac) {
		double g = (gammas->at(0))(stage_index + embedding_shift, stage_index);
		if (g <= 2.0*1e-8) {
			(*jac) = I;
		} else {
			implicit_solve_jacobian(y, &jac_temp);
			(*jac) = I - H*g*jac_temp;
		}
	}

	void implicit_solve_jacobian(vec* y, mat* jac) {
		// Set ti = t + c_i*H
		double ti = t + (*c)(stage_index)*H;

		if (num_omegas == 0) {
			problem->slow_rhsjacobian(ti, y, &jac_temp2);
		} else {
			problem->implicit_rhsjacobian(ti, y, &jac_temp2);
		}
		(*jac) = jac_temp2;
	}
	
	void erk_step(vec* y_prev, vec* f) {	
		f_temp.zeros();
		if (method_type == 0) {
			if (num_omegas == 0) {
				erk_step_mrigark(y_prev,f);
			} else {
				erk_step_imexmrigark(y_prev,f);
			}
		} else if (method_type == 2) {
			erk_step_imexmrisr(y_prev,f);
		}
	}
	
	void erk_step_mrigark(vec* y_prev, vec* f) {
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp2 = mri_eval_stages1.col(sub_stage_index);
			f_temp += gamma_bar(stage_index + embedding_shift, sub_stage_index)*f_temp2;
		}
		(*f) = *y_prev + H*f_temp;
	}
	
	void erk_step_imexmrigark(vec* y_prev, vec* f) {
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp2 = mri_eval_stages1.col(sub_stage_index);
			f_temp += gamma_bar(stage_index + embedding_shift, sub_stage_index)*f_temp2;
			
			f_temp2 = mri_eval_stages2.col(sub_stage_index);
			f_temp += omega_bar(stage_index + embedding_shift, sub_stage_index)*f_temp2;
		}
		(*f) = *y_prev + H*f_temp;
	}
	
	void erk_step_imexmrisr(vec* y_prev, vec* f) {
		for(int sub_stage_index=0; sub_stage_index<stage_index; sub_stage_index++) {
			f_temp2 = mri_eval_stages1.col(sub_stage_index);
			f_temp += (gammas->at(0))(stage_index + embedding_shift, sub_stage_index)*f_temp2;
		}
		(*f) = *y_prev + H*f_temp;
	}
	
	double gamma(int i, int j, double tau) {
		double gamma_val = 0.0;
		for(int k=0; k<num_gammas; k++) {
			gamma_val += (gammas->at(k))(i,j)*pow(tau,k);
		}
		return gamma_val;
	}
	
	double gamma_mrexp(int i, int j, double tau) {
		double gamma_val = 0.0;
		for(int k=0; k<num_gammas; k++) {
			gamma_val += (gammas->at(k))(i,j)*pow(tau,k+1);
		}
		return gamma_val;
	}
	
	double omega(int i, int j, double tau) {
		double omega_val = 0.0;
		for(int k=0; k<num_omegas; k++) {
			omega_val += (omegas->at(k))(i,j)*pow(tau,k);
		}
		return omega_val;
	}

	double gamma_bar(int i, int j) {
		double gamma_bar_val = 0.0;
		for(int k=0; k<num_gammas; k++) {
			gamma_bar_val += (gammas->at(k))(i,j)/(k+1);
		}
		return gamma_bar_val;
	}
	
	double omega_bar(int i, int j) {
		double omega_bar_val = 0.0;
		for(int k=0; k<num_omegas; k++) {
			omega_bar_val += (omegas->at(k))(i,j)/(k+1);
		}
		return omega_bar_val;
	}

	void set_stage_dependent_data(int stage_index_, bool embedding_) {
		stage_index = stage_index_;
		embedding = embedding_;
		if (embedding_) {
			embedding_shift = 1;
		} else {
			embedding_shift = 0;
		}
	}
	
	void set_step_dependent_data(double H_, double t_) {
		H = H_;
		t = t_;
	}

};

#endif