#ifndef PLEIADESPROBLEM_DEFINED__
#define PLEIADESPROBLEM_DEFINED__

#include "Problem.hpp"

using namespace std;
using namespace arma;

class PleiadesProblem : public Problem {
public:
	double g = 1.0;
	double softening_length_squared = 0.0;
	double masses[7] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
	double n = 7;
	double d = 2;
	vec delta_pos;
	vec delta_accel;
	int i_start;
	int j_start;
	double square_sum;

	PleiadesProblem() {
		name = "Pleiades";
		problem_dimension = 28;
		default_H = std::pow(2.0,-9.0);
		t_0 = 0.0;
		t_f = 3.0;
		has_true_solution = false;
		explicit_only = true;
		y_0 = { 3.0, 3.0, 3.0, -3.0, -1.0, 2.0, -3.0, 0.0, 2.0, 0.0, -2.0, -4.0, 2.0, 4.0,
    		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.25, 0.0, 1.0, 1.75, 0.0, -1.5, 0.0};
		
		delta_pos = vec(d,fill::zeros);
		delta_accel = vec(d,fill::zeros);
	}
	
	void full_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();

		for(int i=0; i<n*d; i++) {
			(*f)(i) = (*y)(n*d + i);
		}

		for(int i=0; i<n; i++) {
			i_start = i*d;

			for(int j=i+1; j<n; j++) {
				j_start = j*d;

				for(int k=0; k<d; k++) {
					delta_pos(k) = (*y)(j_start + k) - (*y)(i_start + k);
				}

				square_sum = sum(square(delta_pos));
				delta_accel = g*delta_pos/(std::pow(square_sum+softening_length_squared, 1.5));
				for(int k=0; k<d; k++) {
					(*f)(n*d + i_start + k) += masses[j]*delta_accel(k);
				}
				for(int k=0; k<d; k++) {
					(*f)(n*d + j_start + k) -= masses[i]*delta_accel(k);
				}
			}
		}
	}
	
	void fast_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();

		for(int i=0; i<n; i++) {
			i_start = i*d;

			for(int j=i+1; j<n; j++) {
				j_start = j*d;

				for(int k=0; k<d; k++) {
					delta_pos(k) = (*y)(j_start + k) - (*y)(i_start + k);
				}

				square_sum = sum(square(delta_pos));
				delta_accel = g*delta_pos/(std::pow(square_sum+softening_length_squared, 1.5));
				for(int k=0; k<d; k++) {
					(*f)(n*d + i_start + k) += masses[j]*delta_accel(k);
				}
				for(int k=0; k<d; k++) {
					(*f)(n*d + j_start + k) -= masses[i]*delta_accel(k);
				}
			}
		}
	}
	void slow_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
		for(int i=0; i<n*d; i++) {
			(*f)(i) = (*y)(n*d + i);
		}
	}

	void linear_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();
		for(int i=0; i<n*d; i++) {
			(*f)(i) = (*y)(n*d + i);
		}
	}
	
	void nonlinear_rhs_custom(double t, vec* y, vec* f) {
		f->zeros();

		for(int i=0; i<n; i++) {
			i_start = i*d;

			for(int j=i+1; j<n; j++) {
				j_start = j*d;

				for(int k=0; k<d; k++) {
					delta_pos(k) = (*y)(j_start + k) - (*y)(i_start + k);
				}

				square_sum = sum(square(delta_pos));
				delta_accel = g*delta_pos/(std::pow(square_sum+softening_length_squared, 1.5));
				for(int k=0; k<d; k++) {
					(*f)(n*d + i_start + k) += masses[j]*delta_accel(k);
				}
				for(int k=0; k<d; k++) {
					(*f)(n*d + j_start + k) -= masses[i]*delta_accel(k);
				}
			}
		}
	}
};

#endif