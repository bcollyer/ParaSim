#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"


struct Population{
  std::vector<Person> hosts;
  double mean_burden;
  double prevalence;
  double reservoir;
  double location_x;
  double location_y;
};

Population create_population(int nhosts, double m, double k, double x, double y);
//void evolve_population(Population &population, int timesteps);
void set_prevalence(Population &population);
void set_meanburden(Population &population);
