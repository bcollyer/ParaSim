#include <vector>
#include "Parameters.h"
#include "Population.h"

struct Simulation
{
  std::vector<Population> populations;
  std::vector<Parameters> parameters;
};

Simulation create_simulation(std::string filename,int n_villages);
Simulation load_from_map(std::string map_filename, std::string param_filename);
void evolve_population(Population &population,Parameters &parameters, int timesteps);
void evolve_populations(Simulation &sim,int timesteps);
void output_prevalence(Simulation &sim);
void output_meanburden(Simulation &sim);
