#include "Population.h"




void set_prevalence(Population &population)
{
  int n_hosts;
  double prevalence;

  n_hosts = population.hosts.size();
  prevalence = 0.0;

  for (int i = 0; i<n_hosts;i++)
  {
      prevalence += (int) (population.hosts[i].burden.eggs_test > 0);
  }

  population.prevalence = 1.0*prevalence/n_hosts;
};

void set_meanburden(Population &population){
  int n_hosts;
  double meanburden;

  n_hosts = population.hosts.size();
  meanburden = 0.0;

  for (int i = 0; i<n_hosts;i++)
  {
      meanburden += population.hosts[i].burden.eggs;
  }

  population.mean_burden = 1.0*meanburden/n_hosts;
};

Population create_population(int nhosts, double m, double k, double x, double y)
{
  // Set up random number generation
  gsl_rng * rando;
  const gsl_rng_type * T;
  Population population;
  std::vector<Person> hosts;
  double rate;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rando = gsl_rng_alloc(T);
  gsl_rng_set(rando,time(0));

  for (int i = 0; i < nhosts; i++){
    Person person;
    Burden burden;

    person.age = 10.0;
    person.risk = gsl_ran_gamma(rando,k,1.0/k);

    rate = 0.5 * person.risk * m;
    burden.male_worms = gsl_ran_poisson(rando,rate);
    burden.female_worms = gsl_ran_poisson(rando,rate);
    burden.eggs = 0;
    burden.eggs_test = 0;


    person.burden = burden;
    hosts.push_back(person);
  };
  population.hosts = hosts;
  population.reservoir = 1000.0;
  population.location_x = x;
  population.location_y = y;


  return population;
};
