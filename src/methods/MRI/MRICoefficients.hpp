#ifndef MRICOEFFICIENTS_DEFINED__
#define MRICOEFFICIENTS_DEFINED__

using namespace arma;

class MRICoefficients {
public:
	const char* name;
	bool explicit_mrigark;
	int method_type; // 0: MRIGARK, 1: MERK/MERB
	std::vector<mat> gammas;
	std::vector<mat> omegas;
	vec c;
	int num_stages;
	int num_gammas;
	int num_omegas;
	int num_groups;
	std::vector<std::vector<int>> stage_groups;
	double primary_order;
	double secondary_order;
};

#endif