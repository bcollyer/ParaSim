//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <sstream>
#include <numeric>
//#include <numeric>
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
    population = create_population(parameters.nhosts,2.0,parameters.k,0.0,0.0,i);
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
  int index = 0;
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

  parameters = read_from_file(param_filename);


  if (myfile.is_open())
  {
    std::string line;
    while (getline(myfile, line))
    {
      std::istringstream iss(line);
      iss >> x >> y >> pop;

      //std::cout << pop << "\n";
      //parameters = read_from_file(param_filename);
      contact_rate = exp(-1.0*pop/4000.0)*gsl_ran_flat(rando,0.0,0.00015);
      //contact_rate = exp(-1.0*pop/2000.0)*0.0002;
      //contact_rate = 0.00002;
      parameters.contact_rate = contact_rate;
      //std::cout << parameters.nhosts << "\n";
      //population = create_population(parameters.nhosts,2.0,parameters.k,0.0,0.0);
      population = create_population(pop,2.0,parameters.k,x,y,index);
      //std::cout << "created pop" << "\n";

      set_prevalence(population);
      set_meanburden(population);

      villages.push_back(population);
      params.push_back(parameters);
      index++;
      //std::cout << population.mean_burden << "\n";
    }
    sim.populations = villages;
    sim.parameters = params;
  }


  gsl_rng_free(rando);
  return sim;
}

void evolve_population(Population &population,Parameters &parameters, int timesteps)
{
  // Set up random number generation
  //gsl_rng * rando;
  //const gsl_rng_type * T;
  double rate = 0.0;
  double m, p, k;
  int deaths, births, total_worms;
  int nhosts;
  double dt = 1.0/12.0;
  double t = 0.0;
  double death_check = 0.0;

  //index = population.location_x + 10*population.location_y;

  //gsl_rng_env_setup();
  //T = gsl_rng_default;
  //rando = gsl_rng_alloc(T);
  //gsl_rng_set(rando,index+time(0));


  nhosts = population.hosts.size();
  for ( int i=0; i < timesteps; i++)
  {

    // depletion of reservoir
    population.reservoir = (population.reservoir*(1.0-1.0*dt));

    // loop over hosts
    for (int j=0; j < nhosts; j++)
    {

      population.hosts[j].age += dt;
      //population.hosts[j].b += dt;


      // human death and replacement birth
      death_check = gsl_ran_flat(population.rando,0,1);
      if (death_check < (1.0 - exp(-1.0 * parameters.human_death_rate * dt))) {
        population.hosts[j].age = 0.0;
        Burden burden;
        burden.male_worms = 0;
        burden.female_worms = 0;
        burden.eggs = 0;
        burden.eggs_test = 0;
        population.hosts[j].burden = burden;
        population.hosts[j].risk = gsl_ran_gamma(population.rando,parameters.k,1.0/parameters.k);
      }

      // parasite death
      rate = population.hosts[j].burden.male_worms * parameters.parasite_death_rate * dt;
      if (rate < 10E3) {
        deaths = gsl_ran_poisson(population.rando,rate);
      } else {
        deaths = (int) (rate + gsl_ran_gaussian(population.rando,sqrt(rate)));
      }
      population.hosts[j].burden.male_worms -= deaths;
      population.hosts[j].burden.male_worms = std::max(0,population.hosts[j].burden.male_worms);

      rate = population.hosts[j].burden.female_worms * parameters.parasite_death_rate * dt;
      if (rate < 10E3) {
        deaths = gsl_ran_poisson(population.rando,rate);
      } else {
        deaths = (int) (rate + gsl_ran_gaussian(population.rando,sqrt(rate)));
      }
      population.hosts[j].burden.female_worms -= deaths;
      population.hosts[j].burden.female_worms = std::max(0,population.hosts[j].burden.female_worms);



      //parasite aquisition with equal sex probability
      rate = 0.5*population.reservoir * parameters.contact_rate * population.hosts[j].risk * dt;
      rate = rate * population.hosts[j].home;
      //rate = 10.0;
      if (rate < 10E3) {
        births = gsl_ran_poisson(population.rando,rate);
      } else {
        births = (int) (rate + gsl_ran_gaussian(population.rando,sqrt(rate)));
      }
      population.hosts[j].burden.female_worms += births;

      if (rate < 10E3) {
        births = gsl_ran_poisson(population.rando,rate);
      } else {
        births = (int) (rate + gsl_ran_gaussian(population.rando,sqrt(rate)));
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

      population.hosts[j].burden.eggs_test  =  0.5*(gsl_ran_poisson(population.rando,m)+gsl_ran_poisson(population.rando,m));
      // contribution to reservoir
      population.reservoir +=  population.hosts[j].home * population.hosts[j].burden.eggs*dt;
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
  //gsl_rng_free(rando);
}

void evolve_populations(Simulation &sim,int timesteps)
{
  int n, n_temp, n_hosts, n_mig;
  int pop_idx, host_idx;
  int births;
  double rate;
  std::vector< std::vector<double> > migration_probs;
  std::vector<double>  temp_probs;
  std::vector<int> temp_v;
  std::vector<double> temp_v2;
  double temp_p, dist;
  double dt = 1.0/12.0;
  double total = 0;
  int out_freq = 0;
  gsl_rng * rando;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rando = gsl_rng_alloc(T);
  gsl_rng_set(rando,1);


  n = sim.populations.size();

  // calculate migration probabilities for each subpopulation using gravity model
  for (int i = 0;i < n; i++)
  {
    for (int j = 0;j < n; j++)
    {
      dist = pow(pow(sim.populations[j].location_x - sim.populations[i].location_x,2.0) + pow(sim.populations[j].location_y -sim.populations[i].location_y,2.0),0.5);
      dist = pow((1+dist/10.0),-1.5);
      temp_p = (1-(i==j)) * 0.0000001 * sim.populations[j].hosts.size() * dist; //migration from i to j
      temp_p = 1.0 - exp (-1.0*temp_p/12.0);
      temp_probs.push_back(temp_p);
    }
    migration_probs.push_back(temp_probs);
    temp_probs.clear();
  }
  //migration_probs[source_ID][dest_ID] = prob

  // store migration rates into each sub-pop and threshold
  for (int i = 0; i < n; i++)
  {
    for (int j = 0;j < n; j++)
    {
      if ( migration_probs[i][j] * sim.populations[i].hosts.size() * 12 * 10 > 1.0 )
      {
        temp_v2.push_back(1.0*i); //source index
        temp_v2.push_back(migration_probs[i][j]); //probability
        //std::cout << j << " "<< temp_v2[0] << " " << temp_v2[1] << "\n";
        sim.populations[j].migration_probs.push_back(temp_v2);
        temp_v2.clear();
      }
    }
  }



  // loop through time
  for (int i = 0;i < timesteps; i++)
  {

    // evolve each population 1 time-period
    for (int j = 0; j < n; j++)
    {
      evolve_population(sim.populations[j],sim.parameters[j],1);

      for (int k = 0; k < n; k++)
      {
        // number of hosts going from pop j to pop k
        n_hosts = sim.populations[j].hosts.size();
        n_mig = gsl_ran_binomial(rando, migration_probs[j][k], n_hosts); //number of hosts going from j -> k
        //std::cout << n_mig << "\n";
        for (int itr = 0; itr < n_mig; itr++)
        {
          temp_v.push_back(j);
          temp_v.push_back(itr);
          sim.populations[k].temp_hosts.push_back(temp_v);
          temp_v.clear();
        }

      }
    }

    reservoir_coupling(sim, n);

    //gravity model
    for (int j = 0; j < n; j++)
    {
      n_temp = sim.populations[j].temp_hosts.size();

      for (int k = 0; k < n_temp; k++)
      {
          pop_idx = sim.populations[j].temp_hosts[k][0];
          host_idx = sim.populations[j].temp_hosts[k][1];

          //infect hosts from other locations
          rate = 0.5*sim.populations[j].reservoir * sim.parameters[j].contact_rate * dt;
          if (rate < 10E3) {
            births = gsl_ran_poisson(sim.populations[j].rando,rate);
          } else {
            births = (int) (rate + gsl_ran_gaussian(sim.populations[j].rando,sqrt(rate)));
          }
          sim.populations[ pop_idx ].hosts[ host_idx ].burden.female_worms += births;

          if (rate < 10E3) {
            births = gsl_ran_poisson(sim.populations[j].rando,rate);
          } else {
            births = (int) (rate + gsl_ran_gaussian(sim.populations[j].rando,sqrt(rate)));
          }
          sim.populations[ pop_idx ].hosts[ host_idx ].burden.male_worms += births;
          // contributions to reservoir
      }
      sim.populations[j].temp_hosts.clear();
    }


    // output prevalence every 12 months
    out_freq = i%12;
    if (out_freq == 0)
    {
      std::cout << i << " ";
      output_prevalence(sim);
    }
  }
}

void reservoir_coupling(Simulation &sim, int &n)
{
  double temp, dist;
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
