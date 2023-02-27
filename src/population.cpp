#include "stdafx.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <process.h>
#include "person.h"
#include "population.h"

void set_prevalence(Population &population)
{
	// records prevalence

	int n_hosts;
	double prevalence;

	n_hosts = population.hosts.size();
	prevalence = 0.0;

	for (int i = 0; i<n_hosts; i++)
	{
		prevalence += (int)(population.hosts[i].burden.eggs_test > 0);
	}

	population.prevalence = 1.0*prevalence / n_hosts;
};

void set_meanburden(Population &population) {

	// record mean egg burden

	int n_hosts;
	double meanburden;

	n_hosts = population.hosts.size();
	meanburden = 0.0;

	for (int i = 0; i<n_hosts; i++)
	{
		meanburden += population.hosts[i].burden.eggs;
	}

	population.mean_burden = 1.0*meanburden / n_hosts;
};

Population create_population(int nhosts, double m, double k, double x, double y, int index)
{
	Population pop;
	std::vector<Person> hosts;
	std::vector< std::vector<int> > temp_hosts;
	double rate;
	double temp_age;

	gsl_rng_env_setup();
	pop.T = gsl_rng_default;
	pop.rando = gsl_rng_alloc(pop.T);
	gsl_rng_set(pop.rando, _getpid() + index); //TODO: fix

	for (int i = 0; i < nhosts; i++) {
		Person person;
		Burden burden;
		//std::cout << i << "\n";
		temp_age = gsl_ran_exponential(pop.rando, 40.0);
		while (temp_age > 80) {
			temp_age = gsl_ran_exponential(pop.rando, 40.001);
		}
		person.age = temp_age;
		person.home = 1; // person resides in home cell

		person.risk = gsl_ran_gamma(pop.rando, k, 1.0 / k);
		rate = 0.5 * person.risk * m;
		burden.male_worms = gsl_ran_poisson(pop.rando, rate);
		burden.female_worms = gsl_ran_poisson(pop.rando, rate);
		burden.eggs = gsl_ran_poisson(pop.rando, rate);
		burden.eggs_test = gsl_ran_poisson(pop.rando, rate);
		person.burden = burden;
		hosts.push_back(person);
	};
	pop.hosts = hosts;
	pop.reservoir = 0.0;
	pop.location_x = x;
	pop.location_y = y;
	pop.temp_hosts = temp_hosts;

	//gsl_rng_free(rando);
	return pop;
};