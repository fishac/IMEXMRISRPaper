#ifndef FASTERRORMEASUREMENTTYPES_DEFINED__
#define FASTERRORMEASUREMENTTYPES_DEFINED__

namespace FastError {
	static const char* types[] = { "FS", "SA-mean", "SA-max", "LASA-mean", "LASA-max" };
	std::vector<std::string> types_strings = { "FS", "SA-mean", "SA-max", "LASA-mean", "LASA-max" };
	// "Slow" types because the error estimation happens primarily at the slow scale (the multirate method stages)
	static const char* slow_types[] = { "FS", "SA-mean", "SA-max" };
	// "Fast" types because the error estimation happens primarily at the fast scale (the inner method steps)
	static const char* fast_types[] = { "LASA-mean", "LASA-max" };

	static const char* all_types[] = { "FS", "SA-sum", "SA-mean", "SAmax", "LASA-sum", "LASA-mean", "LASA-max" };
	static const char* all_slow_types[] = { "FS", "SA-sum", "SA-mean", "SAmax" };
	static const char* all_fast_types[] = { "LASA-sum", "LASA-mean", "LASA-max" };

	int is_FS(const char* measurement_type) {
		return strcmp(measurement_type,"FS") == 0;
	}

	int is_SAsum(const char* measurement_type) {
		return strcmp(measurement_type,"SA-sum") == 0;
	}

	int is_SAmean(const char* measurement_type) {
		return strcmp(measurement_type,"SA-mean") == 0;
	}

	int is_SAmax(const char* measurement_type) {
		return strcmp(measurement_type,"SA-max") == 0;
	}

	int is_LASAsum(const char* measurement_type) {
		return strcmp(measurement_type,"LASA-sum") == 0;
	}

	int is_LASAmean(const char* measurement_type) {
		return strcmp(measurement_type,"LASA-mean") == 0;
	}

	int is_LASAmax(const char* measurement_type) {
		return strcmp(measurement_type,"LASA-max") == 0;
	}

	int is_slow_type(const char* measurement_type) {
		return ((strcmp(measurement_type,"FS") == 0)
		|| (strcmp(measurement_type,"SA-sum") == 0)
		|| (strcmp(measurement_type,"SA-mean") == 0)
		|| (strcmp(measurement_type,"SA-max") == 0));
	}

	int is_fast_type(const char* measurement_type) {
		return ((strcmp(measurement_type,"LASA-sum") == 0)
		|| (strcmp(measurement_type,"LASA-sum") == 0)
		|| (strcmp(measurement_type,"LASA-max") == 0));
	}
}

#endif