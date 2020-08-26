#include <iostream>
//#include <vector>
//#include "Person.h"
//#include "Population.h"
//#include "Parameters.h"
#include "Simulation.h"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    int nhosts;
    int itr;
    int nx, ny;
    double rate;
    double s;
    Simulation sim;
    std::vector<Simulation> simulations;


    //cout << "Hello, World!\n";
    //cout << "Enter number of hosts: " << flush;
    //cin >> nhosts;

    nx = 16;
    ny = 16;
    s = 0.8;
    sim = create_simulation("params.txt",nx*ny);


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
    //cout << "\n";

    for (int i = 0; i < nx; i++){
        for ( int j=0; j <ny; j++){
          itr = i + j*10;

          //cout << "0";
          //cout << sim.parameters[itr].contact_rate ;
          //cout << " ";
        }
        //cout << "\n";
    }

    //for  (int i=0; i< nx*ny; i++){
    //    evolve_population(sim.populations[i],sim.parameters[i],1000);
    //}

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


    for (int i = 0; i < nx; i++){
        for ( int j=0; j <ny; j++){
          itr = i + j*10;
          //std::cout << std::fixed;
          //std::cout << std::setprecision(2);
          //cout << sim.populations[itr].prevalence;
          //cout << " ";
        }
        //cout << "\n";
    }
    //cout << "\n";

    for (int i = 0; i < nx; i++){
        for ( int j=0; j <ny; j++){
          itr = i + j*10;
          //std::cout << std::fixed;
          //std::cout << std::setprecision(2);
          //cout << sim.populations[itr].mean_burden;
          //cout << " ";
        }
        //cout << "\n";
    }


    //sim = create_simulation("../data/params.txt");
    //evolve_population(sim.populations[0],sim.parameters,100);

    return 0;
}
