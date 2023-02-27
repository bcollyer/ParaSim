#pragma once
#include <vector>
#include "parameters.h"
#include "population.h"

struct Simulation
{
	std::vector<Population> populations;
	std::vector<Parameters> parameters;
	std::vector<double> mda_group;
	int total_time;
};

Simulation create_simulation(std::string filename, int n_villages);
Simulation load_from_map(std::string map_filename, std::string param_filename);
void evolve_population(Population &population, Parameters &parameters, int timesteps);
//void evolve_population(Simulation &sim, int timesteps);
void evolve_populations(Simulation &sim, int timesteps);
void reservoir_coupling(Simulation &sim, int &n);
void output_prevalence(Simulation &sim);
void output_meanburden(Simulation &sim);
void sac_mda_pop(Population &population, Parameters &parameters, double sac_frac, double sac_cov);
void sac_mda_sim(Simulation &sim);
void adult_mda_pop(Population &population, Parameters &parameters, double adult_frac, double adult_cov);
void adult_mda_sim(Simulation &sim,  std::vector<double> group);
double age_dep_contact(Parameters &parameters, double age);
