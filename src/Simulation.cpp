//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <ctime>
#include <cmath>
#include <sstream>
#include "Simulation.h"

Simulation create_simulation(std::string filename,int n_villages)
{
  Parameters parameters;
  Population population;
  std::vector<Population> villages;
  std::vector<Parameters> params;
  Simulation sim;


  for (int i =0; i < n_villages; i++){
    parameters = read_from_file(filename);
    population = create_population(parameters.nhosts,2.0,parameters.k,0.0,0.0);
    set_prevalence(population);
    set_meanburden(population);


    villages.push_back(population);
    params.push_back(parameters);
  }

  sim.populations = villages;
  sim.parameters = params;

  return sim;
}

Simulation load_from_map(std::string map_filename, std::string param_filename)
{
  int pop;
  double x;
  double y;
  double contact_rate;
  Parameters parameters;
  Population population;
  std::vector<Population> villages;
  std::vector<Parameters> params;
  Simulation sim;
  std::ifstream myfile (map_filename);
  gsl_rng * rando;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rando = gsl_rng_alloc(T);
  gsl_rng_set(rando,123);

  if (myfile.is_open())
  {
    std::string line;
    while (getline(myfile, line))
    {
      std::istringstream iss(line);
      iss >> x >> y >> pop;

      //std::cout << pop << "\n";
      parameters = read_from_file(param_filename);
      contact_rate = exp(-1.0*pop/2000.0)*gsl_ran_flat(rando,0.0,0.00035);
      parameters.contact_rate = contact_rate;
      //std::cout << parameters.nhosts << "\n";
      //population = create_population(parameters.nhosts,2.0,parameters.k,0.0,0.0);
      population = create_population(pop,2.0,parameters.k,x,y);
      set_prevalence(population);
      set_meanburden(population);

      villages.push_back(population);
      params.push_back(parameters);

      //std::cout << population.mean_burden << "\n";
    }
    sim.populations = villages;
    sim.parameters = params;
  }



  return sim;
}

void evolve_population(Population &population,Parameters &parameters, int timesteps)
{
  // Set up random number generation
  gsl_rng * rando;
  const gsl_rng_type * T;
  int index;
  double rate, m, p, k;
  int deaths, births, total_worms;
  int nhosts;
  double dt = 1.0/12.0;
  double t = 0.0;

  //index = population.location_x + 10*population.location_y;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rando = gsl_rng_alloc(T);
  gsl_rng_set(rando,index+time(0));


  nhosts = population.hosts.size();
  for ( int i=0; i < timesteps; i++)
  {

    // depletion of reservoir
    population.reservoir = (population.reservoir*(1.0-1.0*dt));

    // loop over hosts
    for (int j=0; j < nhosts; j++)
    {

      population.hosts[j].age += dt;

      // human death and replacement birth
      if (gsl_ran_flat(rando,0,1) < (1.0 - exp(-1.0 * parameters.human_death_rate * dt))) {
        population.hosts[j].age = 0.0;
        Burden burden;
        burden.male_worms = 0;
        burden.female_worms = 0;
        burden.eggs = 0;
        burden.eggs_test = 0;
        population.hosts[j].burden = burden;
        population.hosts[j].risk = gsl_ran_gamma(rando,parameters.k,1.0/parameters.k);
      }

      // parasite death
      rate = population.hosts[j].burden.male_worms * parameters.parasite_death_rate * dt;
      if (rate < 10E3) {
        deaths = gsl_ran_poisson(rando,rate);
      } else {
        deaths = (int) (rate + gsl_ran_gaussian(rando,sqrt(rate)));
      }
      population.hosts[j].burden.male_worms -= deaths;
      population.hosts[j].burden.male_worms = std::max(0,population.hosts[j].burden.male_worms);

      rate = population.hosts[j].burden.female_worms * parameters.parasite_death_rate * dt;
      if (rate < 10E3) {
        deaths = gsl_ran_poisson(rando,rate);
      } else {
        deaths = (int) (rate + gsl_ran_gaussian(rando,sqrt(rate)));
      }
      population.hosts[j].burden.female_worms -= deaths;
      population.hosts[j].burden.female_worms = std::max(0,population.hosts[j].burden.female_worms);



      //parasite aquisition with equal sex probability
      rate = 0.5*population.reservoir * parameters.contact_rate * population.hosts[j].risk * dt;
      //rate = 10.0;
      if (rate < 10E3) {
        births = gsl_ran_poisson(rando,rate);
      } else {
        births = (int) (rate + gsl_ran_gaussian(rando,sqrt(rate)));
      }
      population.hosts[j].burden.female_worms += births;

      if (rate < 10E3) {
        births = gsl_ran_poisson(rando,rate);
      } else {
        births = (int) (rate + gsl_ran_gaussian(rando,sqrt(rate)));
      }
      population.hosts[j].burden.male_worms += births;



      // egg/larvae production
      total_worms = population.hosts[j].burden.male_worms + population.hosts[j].burden.female_worms;
      population.hosts[j].burden.eggs = (int) parameters.fecundity * std::min(population.hosts[j].burden.male_worms,population.hosts[j].burden.female_worms) * exp(-1.0 * parameters.gamma * total_worms);

      m = parameters.fecundity * std::min(population.hosts[j].burden.male_worms,population.hosts[j].burden.female_worms) * exp(-1.0 * parameters.gamma * total_worms);
      p = parameters.k_egg/ (m + parameters.k_egg);

      /*
      if (p <10E-6){
        population.hosts[j].burden.eggs_test = 0;
      } else {
        //population.hosts[j].burden.eggs = (int) m;//gsl_ran_negative_binomial (rando, p, parameters.k_egg);

        population.hosts[j].burden.eggs_test  = (int) ceil((gsl_ran_negative_binomial (rando, p, parameters.k_egg)
                                              + gsl_ran_negative_binomial (rando, p, parameters.k_egg)
                                              + gsl_ran_negative_binomial (rando, p, parameters.k_egg))/3.0);
      }
      */

      population.hosts[j].burden.eggs_test  =  0.5*(gsl_ran_poisson(rando,m)+gsl_ran_poisson(rando,m));
      // contribution to reservoir
      population.reservoir += population.hosts[j].burden.eggs*dt;
    }



    t+=dt;
    set_meanburden(population);
    set_prevalence(population);
    if ((i%12 == 0) && (i>1))
    {
      std::cout << " " << t << " ";
      std::cout << population.prevalence;
    }

    //std::cout << std::fixed;
    //std::cout << std::setprecision(2);
    //std::cout << t <<"  " <<population.reservoir;
    //std::cout << "  " << population.mean_burden << " "<<population.prevalence<<"\n";

  }
}

void evolve_populations(Simulation &sim,int timesteps)
{
  int n;
  double temp, dist;

  n = sim.populations.size();
  for (int i = 0;i < timesteps; i++)
  {
    for (int j = 0; j < n; j++)
    {
      evolve_population(sim.populations[j],sim.parameters[j],1);
    }

    for (int j = 0; j < n; j++)
    {
      temp = 0.0;
      for (int k = 0; k < n; k++)
      {
        if (k != j){
          dist = pow(pow(sim.populations[j].location_x - sim.populations[k].location_x,2.0) + pow(sim.populations[j].location_y -sim.populations[k].location_y,2.0),0.5);
          temp += sim.populations[k].reservoir * pow(dist,-2.0);
        }
      }
      sim.populations[j].reservoir += sim.parameters[j].spat_kernel_magnitude * temp;
    }

    if ((i%12 == 0))
    {
      std::cout << i << " ";
      output_prevalence(sim);
    }


  }
}

void output_prevalence(Simulation &sim)
{
  int n;
  n = sim.populations.size();
  for (int j = 0; j < n; j++)
  {
    std::cout << sim.populations[j].prevalence << " ";
  }
  std::cout << "\n";
}

void output_meanburden(Simulation &sim)
{
  int n;
  n = sim.populations.size();
  for (int j = 0; j < n; j++)
  {
    std::cout << sim.populations[j].mean_burden<< " ";
  }
  std::cout << "\n";
}
