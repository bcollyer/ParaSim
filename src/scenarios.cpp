#include <iostream>
#include <filesystem>
#include "Simulation.h"
#include "scenarios.h"

//#include "Simulation.h"

void scenario_square()
{
  int nhosts;
  int itr;
  int nx, ny;
  double rate;
  double s;
  Simulation sim;
  std::vector<Simulation> simulations;


  nx = 16;
  ny = 16;
  s = 0.8;
  sim = create_simulation("./data/params.txt",nx*ny);


  for (int i = 0; i < nx; i++){
      for ( int j=0; j <ny; j++){
        itr = i + j*nx;
        rate = sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double) i-11,(double) j-11, s,s, 0.0)/gsl_ran_bivariate_gaussian_pdf(0.0,0.0, s,s, 0.0)
              + sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double) i-4,(double) j-11, s,s, 0.0)/gsl_ran_bivariate_gaussian_pdf(0.0,0.0, s,s, 0.0)
              + sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double) i-4,(double) j-4, s,s, 0.0)/gsl_ran_bivariate_gaussian_pdf(0.0,0.0, s,s, 0.0)
              + sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double) i-11,(double) j-4, s,s, 0.0)/gsl_ran_bivariate_gaussian_pdf(0.0,0.0, s,s, 0.0);

        sim.populations[itr].location_y = 1.0*j;
        sim.populations[itr].location_x = 1.0*i;
        sim.parameters[itr].contact_rate = rate;
        //cout << "0";
        //cout << rate;
        //cout << " ";
      }
      //cout << "\n";
  }


  evolve_populations(sim,12*50);
  for (int i = 0; i < nx; i++){
      for ( int j=0; j <ny; j++){
        itr = i + j*nx;
        if ((i<6)&&(j<6))
        {
          for (int k=0;k<sim.populations[itr].hosts.size();k++)
          {
            sim.populations[itr].hosts[k].burden.male_worms = 0;
            sim.populations[itr].hosts[k].burden.female_worms = 0;
            sim.populations[itr].hosts[k].burden.eggs = 0;
          }
        }
        sim.populations[itr].reservoir = 0;
      }
  }
  evolve_populations(sim,12*50);



  //sim = create_simulation("../data/params.txt");
  //evolve_population(sim.populations[0],sim.parameters,100);
}

void scenario_swz()
{
  Simulation sim;
  std::vector<Simulation> simulations;
  int n;

  sim = load_from_map("./data/swz.txt","./data/params.txt");
  std::cout<<"loaded map\n";
  evolve_populations(sim,12*50);

  // free rngs
  n = sim.populations.size();
  for(int i=0; i < n ; i++){
    gsl_rng_free(sim.populations[i].rando);
  }
}
