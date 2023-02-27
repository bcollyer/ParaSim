#include "stdafx.h"
#include <iostream>
#include <filesystem>
#include <fstream>
#include "simulation.h"
#include "scenarios.h"

void scenario_square(double contact_rate)
{

	//
	// Toy example, simulating square shaped meta-pop 
	//

	int itr;
	int nx, ny;
	double rate;
	double s;
	Simulation sim;
	std::vector<Simulation> simulations;


	nx = 16;
	ny = 16;
	s = 0.8;
	//sim = create_simulation("./data/params.txt", nx*ny);
	sim = create_simulation("\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\params.txt", nx*ny);


	for (int i = 0; i < nx; i++) {
		for (int j = 0; j <ny; j++) {
			itr = i + j*nx;
			rate = sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double)i - 11, (double)j - 11, s, s, 0.0) / gsl_ran_bivariate_gaussian_pdf(0.0, 0.0, s, s, 0.0)
				+ sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double)i - 4, (double)j - 11, s, s, 0.0) / gsl_ran_bivariate_gaussian_pdf(0.0, 0.0, s, s, 0.0)
				+ sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double)i - 4, (double)j - 4, s, s, 0.0) / gsl_ran_bivariate_gaussian_pdf(0.0, 0.0, s, s, 0.0)
				+ sim.parameters[itr].contact_rate *gsl_ran_bivariate_gaussian_pdf((double)i - 11, (double)j - 4, s, s, 0.0) / gsl_ran_bivariate_gaussian_pdf(0.0, 0.0, s, s, 0.0);

			sim.populations[itr].location_y = 1.0*j;
			sim.populations[itr].location_x = 1.0*i;
			sim.parameters[itr].contact_rate = contact_rate*rate;
			//cout << "0";
			//cout << rate;
			//cout << " ";
		}
		//cout << "\n";
	}


	evolve_populations(sim, 12 * 10);
	int nhosts = sim.populations.size();

	// free rngs
	for (int i = 0; i < nhosts; i++) {
		gsl_rng_free(sim.populations[i].rando);
	}

	/*
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j <ny; j++) {
			itr = i + j*nx;
			if ((i<6) && (j<6))
			{
				for (int k = 0; k<sim.populations[itr].hosts.size(); k++)
				{
					sim.populations[itr].hosts[k].burden.male_worms = 0;
					sim.populations[itr].hosts[k].burden.female_worms = 0;
					sim.populations[itr].hosts[k].burden.eggs = 0;
				}
			}
			sim.populations[itr].reservoir = 0;
		}
	}
	*/
	//evolve_populations(sim, 12 * 50);



	//sim = create_simulation("../data/params.txt");
	//evolve_population(sim.populations[0],sim.parameters,100);
}

void scenario_swz(double cr1, double cr2, double cr3, double cr4, double cr5)
{
	Simulation sim;
	std::vector<Simulation> simulations;
	

	sim = load_from_map("\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\swz.txt", "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\params.txt");
	int n = sim.populations.size();
	
	if (cr1 > 0)
	{
		
		for (int i = 0; i < n; i++) {
			
			int nhosts = sim.populations[i].hosts.size();
			//sim.parameters[i].contact_rate = exp(-1.0*nhosts / 4000.0)*0.5*contact_rate;
			if (sim.populations[i].location_x < 11) {
				sim.parameters[i].contact_rate = cr1;
			}
			else if (sim.populations[i].location_x < 22) {
				sim.parameters[i].contact_rate = cr2;
			}
			else if (sim.populations[i].location_x < 33) {
				sim.parameters[i].contact_rate = cr3;
			}
			else if (sim.populations[i].location_x < 44) {
				sim.parameters[i].contact_rate = cr4;
			}
			else {
				sim.parameters[i].contact_rate = cr5;
			}
			//sim.parameters[i].contact_rate = 0.5*contact_rate;
		}
	}


	std::cout << "loaded map\n";
	evolve_populations(sim, 12 * 50);

	// free rngs
	for (int i = 0; i < n; i++) {
		gsl_rng_free(sim.populations[i].rando);
	}
}

void scenario_tumikia(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> ks, double theta)
{
	Simulation sim;
	std::vector<Simulation> simulations;


	sim = load_from_map("\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\TUMIKIA_POP.txt", "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\params.txt");

	int n = sim.populations.size();


	//std::cout << n << "\n";
	//std::cout << contact_rates.size() << "\n";
	for (int i = 0; i < n; i++)
	{
		sim.parameters[i].contact_rate = contact_rates[i];
		//sim.parameters[i].k = std::max(-1.0 * log(1.0 - prevalences[0]) / log(1.0 + 1 / (k1/10.0)), 0.1);
		//sim.parameters[i].k = -1.0 * log(1.0 - prevalences[0]) / log(1.0 + 1 / (k1 / 10.0));
		sim.parameters[i].k = ks[i];
		sim.parameters[i].theta = theta;
		sim.parameters[i].presac_contact = 1.0;
		sim.parameters[i].sac_contact = 1.0;
		sim.parameters[i].adult_contact = 1.0;
	}

	std::cout << "loaded map\n";
	evolve_populations(sim, 12 * 75);

	// free rngs
	for (int i = 0; i < n; i++) {
		gsl_rng_free(sim.populations[i].rando);
	}
}

void scenario_tumikia_mda(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> ks, double theta)
{
	Simulation sim;
	std::vector<Simulation> simulations;
	std::vector<int> mda_group;
	std::vector<double> annual = {1.0,3.0}; // 1.0 2.0
	std::vector<double> biannual = {3.0}; //3

	sim = load_from_map("\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\TUMIKIA_POP.txt", "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\params.txt");

	int n = sim.populations.size();
	for (int i = 0; i < n; i++)
	{
		sim.parameters[i].contact_rate = contact_rates[i];
		sim.parameters[i].k = ks[i];
		sim.parameters[i].theta = theta;
	}

	//sim.mda_group = contact_rates;

	if (sim.mda_group.size() != contact_rates.size())
	{
		return;
	}

	std::cout << "loaded map\n";

	// burn in time
	evolve_populations(sim, 12 * 75);

	// sac mda all
	sac_mda_sim(sim);

	// evolve two months
	//evolve_populations(sim, 2);
	
	// community mop up mda
	adult_mda_sim(sim,annual);

	// evolve 4 months
	evolve_populations(sim, 6);

	// community mda in binannual group
	adult_mda_sim(sim,biannual);

	// evolve population 6 months
	evolve_populations(sim, 6);

	// mda sac all and community mop up
	sac_mda_sim(sim);
	adult_mda_sim(sim,annual);

	// evolve 5 months
	evolve_populations(sim, 5);

	// community mda in biannual group
	adult_mda_sim(sim, biannual);
	//sac_mda_sim(sim); // delete later

	// evolve 2 months
	evolve_populations(sim, 7 + 23*12);


	// free rngs
	for (int i = 0; i < n; i++) {
		gsl_rng_free(sim.populations[i].rando);
	}
}


void scenario_tumikia_mda_long(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> ks, double theta)
{
	Simulation sim;
	std::vector<Simulation> simulations;
	std::vector<int> mda_group;
	std::vector<double> annual = { 1.0,3.0 }; // 1.0 2.0
	std::vector<double> biannual = {1.0, 2.0, 3.0}; //3

	sim = load_from_map("\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\TUMIKIA_POP.txt", "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\params.txt");

	int n = sim.populations.size();
	for (int i = 0; i < n; i++)
	{
		sim.parameters[i].contact_rate = contact_rates[i];
		sim.parameters[i].k = ks[i];
		sim.parameters[i].theta = theta;
	}

	//sim.mda_group = contact_rates;
	
	if (sim.mda_group.size() != contact_rates.size())
	{
		return;
	}
	
	std::cout << "loaded map\n";

	// burn in time
	evolve_populations(sim, 12 * 76);

	
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);



	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);

	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);

	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);

	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);

	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);

	// evolve two months
	evolve_populations(sim, 6);
	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);
	// evolve two months
	evolve_populations(sim, 6);

	adult_mda_sim(sim, biannual);
	sac_mda_sim(sim);

	
	
	// evolve 2 months
	evolve_populations(sim, 12 * 25);


	// free rngs
	for (int i = 0; i < n; i++) {
		gsl_rng_free(sim.populations[i].rando);
	}
}

void scenario_gesh(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> ks, double theta, double fileid)
{
	Simulation sim;
	std::vector<Simulation> simulations;
	int nhosts;
	std::ofstream myfile;

	sim = load_from_map("\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\SINGLE_POP.txt", "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\ParaSim\\data\\params_geshas.txt");

	int n = sim.populations.size();


	//std::cout << n << "\n";
	//std::cout << contact_rates.size() << "\n";
	for (int i = 0; i < n; i++)
	{
		sim.parameters[i].contact_rate = contact_rates[i];
		//sim.parameters[i].k = std::max(-1.0 * log(1.0 - prevalences[0]) / log(1.0 + 1 / (k1/10.0)), 0.1);
		//sim.parameters[i].k = -1.0 * log(1.0 - prevalences[0]) / log(1.0 + 1 / (k1 / 10.0));
		sim.parameters[i].k = ks[i];
		sim.parameters[i].theta = 0;
		sim.parameters[i].presac_contact = 0.2;
		sim.parameters[i].sac_contact = 1.5;
		sim.parameters[i].adult_contact = 0.5;
	}

	std::cout << "loaded map\n";
	evolve_populations(sim, 12 * 75);
	
	//output
	nhosts = sim.populations[0].hosts.size();
	std::string text = "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\RParasim\\finaloutput";
	text += std::to_string(fileid);	
	text += ".txt";
	//myfile.open("\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\RParasim\\finaloutput1.txt");
	myfile.open(text);
	if (myfile.is_open())
	{
		for (int count = 0; count < nhosts; count++) {
			myfile << sim.populations[0].hosts[count].age << " " << sim.populations[0].hosts[count].burden.eggs_test << "\n";
		}
		myfile.close();
	}
	
	//output
	for (int i = 0; i < n; i++) {
		gsl_rng_free(sim.populations[i].rando);
	}
}
