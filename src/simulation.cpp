#include "stdafx.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <sstream>
#include <numeric>
//#include <numeric>
#include "simulation.h"


Simulation create_simulation(std::string filename, int n_villages)
{
	// Creates simulation object for n_villages with given parameters

	Parameters parameters;
	Population population;
	std::vector<Population> villages;
	std::vector<Parameters> params;
	std::vector<double> mda_group;
	Simulation sim;


	// create villages and store in vector
	for (int i = 0; i < n_villages; i++) {
		parameters = read_from_file(filename);
		population = create_population(parameters.nhosts, 2.0, parameters.k, 0.0, 0.0, i);
		set_prevalence(population);
		set_meanburden(population);

		villages.push_back(population);
		params.push_back(parameters);
	}

	//read in mda group list
	read_k_from_file(mda_group, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\group_list.txt");

	// set sim
	sim.populations = villages;
	sim.parameters = params;
	sim.mda_group = mda_group;

	return sim;
}

Simulation load_from_map(std::string map_filename, std::string param_filename)
{
	// loads simulation object from map of prevalences

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
	std::vector<double> mda_group;
	std::ifstream myfile(map_filename);

	// initialise RNG
	gsl_rng * rando;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rando = gsl_rng_alloc(T);
	gsl_rng_set(rando, 123);

	parameters = read_from_file(param_filename);
	// hard coded age contact rates TODO change
	parameters.presac_contact = 1.0;
	parameters.sac_contact = 1.0;
	parameters.adult_contact = 1.0;

	if (myfile.is_open())
	{
		std::string line;
		while (getline(myfile, line))
		{
			std::istringstream iss(line);
			iss >> x >> y >> pop;

			//std::cout << pop << "\n";
			//parameters = read_from_file(param_filename);
			contact_rate = gsl_ran_flat(rando, 0.0, 0.15); //0.00015
			//contact_rate = exp(-1.0*pop / 4000.0)*0.5*parameters.contact_rate;

			parameters.contact_rate = contact_rate;
			//std::cout << parameters.nhosts << "\n";
			//population = create_population(parameters.nhosts,2.0,parameters.k,0.0,0.0);
			//population = create_population(pop, 0.75, parameters.k, x, y, index); // THIS USED FOR TUMIKIA
			//std::cout << "created pop" << "\n";
			population = create_population(pop, 5.75, parameters.k, x, y, index); 

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

	//read in mda group list
	read_k_from_file(mda_group, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\Rparasim\\group.txt");
	sim.mda_group = mda_group;

	sim.total_time = 0;
	gsl_rng_free(rando);
	return sim;
}

double age_dep_contact(Parameters &parameters, double age)
{
	// sets age dependent contact rates in parameters
	if (age < 5.0){
		return parameters.presac_contact;
	}
	else if (age < 15.0) {
		return parameters.sac_contact;
	}
	else{
		return parameters.adult_contact;
	}
}

void evolve_population(Population &population, Parameters &parameters, int timesteps)
{
	// evolves population for a number of timesteps

	double rate = 0.0;
	double m, p;
	int deaths, births, total_worms;
	int nhosts;
	double dt = 1.0 / 12.0;  // unit = years
	double t = 0.0;
	double death_check = 0.0;
	double temp_age_contact;

	nhosts = population.hosts.size();
	for (int i = 0; i < timesteps; i++)
	{

		// depletion of reservoir
		//population.reservoir = (population.reservoir*(1.0 - 1.0*dt));
		//population.reservoir = (population.reservoir*0.115); // correct


		// loop over hosts
		for (int j = 0; j < nhosts; j++)
		{

			// increment host agee
			population.hosts[j].age += dt;
			//population.hosts[j].b += dt;

			// set temp dependent contact rate
			if (population.hosts[j].age < 5) {
				temp_age_contact = parameters.presac_contact;
			}
			else if (population.hosts[j].age < 15) {
				temp_age_contact = parameters.sac_contact;
			}
			else {
				temp_age_contact = parameters.adult_contact;
			}

			// human death and replacement birth
			death_check = gsl_ran_flat(population.rando, 0, 1);
			if (death_check < (1.0 - exp(-1.0 * parameters.human_death_rate * dt))) {
				population.hosts[j].age = 0.0;
				Burden burden;
				burden.male_worms = 0;
				burden.female_worms = 0;
				burden.eggs = 0;
				burden.eggs_test = 0;
				population.hosts[j].burden = burden;
				population.hosts[j].risk = gsl_ran_gamma(population.rando, parameters.k, 1.0 / parameters.k);
			}

			if (population.hosts[j].age > 70.1) {
				population.hosts[j].age = 0.0;
				Burden burden;
				burden.male_worms = 0;
				burden.female_worms = 0;
				burden.eggs = 0;
				burden.eggs_test = 0;
				population.hosts[j].burden = burden;
				population.hosts[j].risk = gsl_ran_gamma(population.rando, parameters.k, 1.0 / parameters.k);

			}

			// parasite death
			rate = population.hosts[j].burden.male_worms * parameters.parasite_death_rate * dt;
			if (rate < 10E3) {
				deaths = gsl_ran_poisson(population.rando, rate);
			}
			else {
				deaths = (int)(rate + gsl_ran_gaussian(population.rando, sqrt(rate)));
			}
			population.hosts[j].burden.male_worms -= deaths;
			population.hosts[j].burden.male_worms = std::max(0, population.hosts[j].burden.male_worms);

			rate = population.hosts[j].burden.female_worms * parameters.parasite_death_rate * dt;
			if (rate < 10E3) {
				deaths = gsl_ran_poisson(population.rando, rate);
			}
			else {
				deaths = (int)(rate + gsl_ran_gaussian(population.rando, sqrt(rate)));
			}
			population.hosts[j].burden.female_worms -= deaths;
			population.hosts[j].burden.female_worms = std::max(0, population.hosts[j].burden.female_worms);



			//parasite aquisition with equal sex probability
			rate = 0.5*population.reservoir * parameters.contact_rate * temp_age_contact * population.hosts[j].risk * dt;

			// if the host is not home rate is zero
			rate = rate * population.hosts[j].home;

			// generate poisson if mean low enough, else more efficient to use gaussian approximation
			if (rate < 10E3) {
				births = gsl_ran_poisson(population.rando, rate);
			}
			else {
				births = (int)(rate + gsl_ran_gaussian(population.rando, sqrt(rate)));
			}
			population.hosts[j].burden.female_worms += births;

			if (rate < 10E3) {
				births = gsl_ran_poisson(population.rando, rate);
			}
			else {
				births = (int)(rate + gsl_ran_gaussian(population.rando, sqrt(rate)));
			}
			population.hosts[j].burden.male_worms += births;



			// egg/larvae production
			total_worms = population.hosts[j].burden.male_worms + population.hosts[j].burden.female_worms;
			//population.hosts[j].burden.eggs = (int)parameters.fecundity * std::min(population.hosts[j].burden.male_worms, population.hosts[j].burden.female_worms) * exp(-1.0 * parameters.gamma * population.hosts[j].burden.female_worms);
			population.hosts[j].burden.eggs = (int)parameters.fecundity * (population.hosts[j].burden.female_worms > 0) * population.hosts[j].burden.male_worms * exp(-1.0 * parameters.gamma * population.hosts[j].burden.female_worms);


			//m = parameters.fecundity * std::min(population.hosts[j].burden.male_worms, population.hosts[j].burden.female_worms) * exp(-1.0 * parameters.gamma * population.hosts[j].burden.female_worms);
			m = population.hosts[j].burden.eggs = parameters.fecundity * (population.hosts[j].burden.female_worms > 0) * population.hosts[j].burden.male_worms * exp(-1.0 * parameters.gamma * population.hosts[j].burden.female_worms);

			// p parameter of negbin distribution
			p = parameters.k_egg / (m + parameters.k_egg);


			if (p > 10E6) {
				population.hosts[j].burden.eggs_test = m;
			}
			else {
				population.hosts[j].burden.eggs_test = gsl_ran_negative_binomial(population.rando, p, parameters.k_egg);
				//population.hosts[j].burden.eggs_test = gsl_ran_poisson(population.rando, m);
			}
		}

		population.reservoir = (population.reservoir*0.115); //decay of two weeks hookworm
		//population.reservoir = (population.reservoir*0.606); //decay of two months ascaris

		for (int j = 1;j < nhosts; j++){
			// contribution to reservoir
			population.reservoir += population.hosts[j].home * population.hosts[j].burden.eggs*dt / nhosts;
		}

		// update time
		t += dt;

		// update mean burden and prevalence for output
		set_meanburden(population);
		set_prevalence(population);

		// this output doesnt get called
		if ((i % 12 == 0) && (i>1))
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

void evolve_populations(Simulation &sim, int timesteps)
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
	double dt = 1.0 / 12.0; //unit = years
	double trip_time;
	double total = 0;
	int out_freq = 0;
	gsl_rng * rando;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rando = gsl_rng_alloc(T);
	gsl_rng_set(rando, 1);


	n = sim.populations.size();

	// calculate migration probabilities for each subpopulation using gravity model
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			dist = pow(pow(sim.populations[j].location_x - sim.populations[i].location_x, 2.0) + pow(sim.populations[j].location_y - sim.populations[i].location_y, 2.0), 0.5);
			dist = pow((1 + dist*3.0 / 4.3), -1.22); //each pixel has size 3km, 
			temp_p = (1 - (i == j)) * sim.parameters[j].theta * pow(sim.populations[j].hosts.size(),1.5) * dist; //migration from i to j in 1 year
			temp_p = 1.0 - exp(-1.0*temp_p / 12.0); // migration prob per month
			temp_probs.push_back(temp_p);
		}
		migration_probs.push_back(temp_probs);
		temp_probs.clear();
	}
	//migration_probs[source_ID][dest_ID] = prob

	// store migration rates into each sub-pop and threshold
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (migration_probs[i][j] * sim.populations[i].hosts.size() * 12 > 1.0) // if expected number of people going from i to j per year greater than 0.1 
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
	for (int i = 0; i < timesteps; i++)
	{

		// evolve each population 1 time-period
		for (int j = 0; j < n; j++)
		{
			evolve_population(sim.populations[j], sim.parameters[j], 1);

			// select hosts migrating
			for (int k = 0; k < n; k++)
			{
				// number of hosts going from pop j to pop k
				n_hosts = sim.populations[j].hosts.size();
				n_mig = gsl_ran_binomial(rando, migration_probs[j][k], n_hosts); //number of hosts going from j -> k
																				 //std::cout << n_mig << "\n";

				// for each person going from j to k, add record in population k
				for (int itr = 0; itr < n_mig; itr++)
				{
					temp_v.push_back(j);
					temp_v.push_back(itr); // temp_v = [source, index]
					sim.populations[k].temp_hosts.push_back(temp_v); // temp_hosts = [[source1 index1] , [source2 index2] , ... , [source_N indexN]]
					temp_v.clear();
				}

			}
		}

		//reservoir_coupling(sim, n);

		//gravity infection model
		for (int j = 0; j < n; j++)
		{
			// number of migrants in cell j
			n_temp = sim.populations[j].temp_hosts.size();


			for (int k = 0; k < n_temp; k++)
			{
			
				pop_idx = sim.populations[j].temp_hosts[k][0];
				host_idx = sim.populations[j].temp_hosts[k][1];

				//time spent in cell j
				trip_time = gsl_ran_exponential(sim.populations[j].rando, 3.0)/30.0;

				//infect each migrant in cell j
				rate = 0.5*sim.populations[j].reservoir * sim.parameters[j].contact_rate * dt;
				rate = rate * trip_time;
				if (rate < 10E3) {
					births = gsl_ran_poisson(sim.populations[j].rando, rate);
				}
				else {
					births = (int)(rate + gsl_ran_gaussian(sim.populations[j].rando, sqrt(rate))); //gaussian quicker to generate than Poisson for large mean
				}
				sim.populations[pop_idx].hosts[host_idx].burden.female_worms += births;

				if (rate < 10E3) {
					births = gsl_ran_poisson(sim.populations[j].rando, rate);
				}
				else {
					births = (int)(rate + gsl_ran_gaussian(sim.populations[j].rando, sqrt(rate)));
				}
				sim.populations[pop_idx].hosts[host_idx].burden.male_worms += births;


				// contributions to reservoir
				sim.populations[j].reservoir += trip_time * sim.populations[pop_idx].hosts[host_idx].burden.eggs / sim.populations[j].hosts.size();

			}

			//migrants leave after 1 month
			sim.populations[j].temp_hosts.clear();
		}

		// output prevalence every 12 months
		//out_freq = sim.total_time % 12;
		out_freq = 0;// sim.total_time % 12;
		if (out_freq == 0)
		{
			std::cout << sim.total_time << " ";
			output_prevalence(sim);
		}
		sim.total_time += 1;

	}
}

void reservoir_coupling(Simulation &sim, int &n)
{
	// NOT USED YET couple reservoirs, maybe useful for schistosomiasis

	double temp, dist;
	for (int j = 0; j < n; j++)
	{
		temp = 0.0;
		for (int k = 0; k < n; k++)
		{
			if (k != j) {
				dist = pow(pow(sim.populations[j].location_x - sim.populations[k].location_x, 2.0) + pow(sim.populations[j].location_y - sim.populations[k].location_y, 2.0), 0.5);
				temp += sim.populations[k].reservoir * pow(dist, -2.0);
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
		std::cout << sim.populations[j].mean_burden << " ";
	}
	std::cout << "\n";
}

void sac_mda_pop(Population &population, Parameters &parameters, double sac_frac,double sac_cov)
{
	int nhosts = population.hosts.size();
	for (int i = 0; i < nhosts; i++)
	{
		if (population.hosts[i].age < 15.0) {
			if (gsl_ran_bernoulli(population.rando, sac_cov) == 1)
			{
				//population.hosts[i].burden.male_worms = 0;//
				population.hosts[i].burden.male_worms = gsl_ran_binomial(population.rando, 0.01, population.hosts[i].burden.male_worms);
				//population.hosts[i].burden.female_worms = 0;//
				population.hosts[i].burden.female_worms = gsl_ran_binomial(population.rando, 0.01, population.hosts[i].burden.female_worms);
				population.hosts[i].burden.eggs = 0;
				population.hosts[i].burden.eggs_test = 0;
			}
		}
	}
	population.reservoir = population.reservoir*0.1;
}


void sac_mda_sim(Simulation &sim)
{
	double sac_frac = 0.5;
	double sac_cov = 0.859;//0.87;
	int n = sim.populations.size();


	for (int j = 0; j < n; j++)
	{
		sac_mda_pop(sim.populations[j], sim.parameters[j],sac_frac,sac_cov);
	}
}

void adult_mda_pop(Population &population, Parameters &parameters, double adult_frac, double adult_cov)
{
	int nhosts = population.hosts.size();
	for (int i = 0; i < nhosts; i++)
	{
		if (population.hosts[i].age > 15.0) {
			if (gsl_ran_bernoulli(population.rando, adult_cov) == 1)
			{
				//population.hosts[i].burden.male_worms = 0;//
				population.hosts[i].burden.male_worms = gsl_ran_binomial(population.rando, 0.01, population.hosts[i].burden.male_worms);
				//population.hosts[i].burden.female_worms = 0;//
				population.hosts[i].burden.female_worms = gsl_ran_binomial(population.rando, 0.01, population.hosts[i].burden.female_worms);
				population.hosts[i].burden.eggs = 0;
				population.hosts[i].burden.eggs_test = 0;

			}

		}
		//population.hosts[i].burden.male_worms = 0;// gsl_ran_binomial(population.rando, 0.01, population.hosts[i].burden.male_worms);
		//population.hosts[i].burden.female_worms = 0;// gsl_ran_binomial(population.rando, 0.01, population.hosts[i].burden.female_worms);
		//population.hosts[i].burden.eggs = 0;
		//population.hosts[i].burden.eggs_test = 0;
	}
	population.reservoir = population.reservoir*(1 - adult_cov);
}


void adult_mda_sim(Simulation &sim,std::vector<double> group)
{
	double sac_frac = 0.5;
	double adult_cov = 0.859;
	int n = sim.populations.size();

	// apply mda to each sub pop
	for (int j = 0; j < n; j++)
	{
		// if population group is in groups-to-treat, treat
		if (std::find(group.begin(), group.end(), sim.mda_group[j]) != group.end())
		{
			adult_mda_pop(sim.populations[j], sim.parameters[j], 1 - sac_frac, adult_cov); 
		}
		
	}
}