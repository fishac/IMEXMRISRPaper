#ifndef WEIGHTEDERRORNORM_DEFINED__
#define WEIGHTEDERRORNORM_DEFINED__

using namespace arma;

class WeightedErrorNorm {
public:
	vec* atol;
	double rtol;
	int dim;
	vec weight_vec;
	
	WeightedErrorNorm(int dim_) {
		dim = dim_;
		weight_vec = vec(dim_, fill::zeros);
	}

	WeightedErrorNorm(vec* atol_, double rtol_) {
		atol = atol_;
		rtol = rtol_;
		dim = atol->n_elem;
		weight_vec = vec(dim, fill::zeros);
	}

	void set_weights(vec* base_vec) {
		for(int i=0; i<dim; i++) {
			weight_vec(i) = 1.0/((*atol)(i) + rtol*std::abs((*base_vec)(i)));
			//printf("weight_vec(%d): %.16f\n",i,weight_vec(i));
		}
	}

	double compute_norm(vec norm_vec) {
		return std::max(norm(norm_vec % weight_vec,2.0),1e-6);
	}

	double compute_norm_nosafe(vec norm_vec) {
		return norm(norm_vec % weight_vec,2.0);
	}

	void reset_weights() {
		weight_vec.zeros();
	}
	
	void set_atol_rtol(vec* atol_, double rtol_) {
		atol = atol_;
		rtol = rtol_;
	}
};

#endif