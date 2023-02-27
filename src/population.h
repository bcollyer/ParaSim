#pragma once
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <process.h>
#include "person.h"

struct Population {
	std::vector<Person> hosts;
	double mean_burden;
	double prevalence;
	double reservoir;
	double location_x;
	double location_y;
	gsl_rng * rando;
	const gsl_rng_type * T;
	std::vector< std::vector<int> > temp_hosts; // temp hosts store migrants from other locations, elements are [pop_index,host_index]
	std::vector< std::vector<double> > migration_probs; // 2D vector. Each row {source_location_ID, migration_prob}
};

Population create_population(int nhosts, double m, double k, double x, double y, int index);
//void evolve_population(Population &population, int timesteps);
void set_prevalence(Population &population);
void set_meanburden(Population &population);