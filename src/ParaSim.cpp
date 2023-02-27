// ParaSim.cpp : Defines the entry point for the console application.


#include "stdafx.h"
#include <iostream>
#include "burden.h"
#include "person.h"
#include "parameters.h"
#include "population.h"
#include "simulation.h"
#include "scenarios.h"

using namespace std;

int main(int argc, const char * argv[])
{

	if (argc == 1) {
		cout << "No scenario selected";
		return 0;
	}
	if (atoi(argv[1]) == 0) {
		double contact_rate = atof(argv[2]);
		scenario_square(contact_rate);
	}

	
	if (atoi(argv[1]) == 1) {
		double contact_rate = -1;
		scenario_swz(-1, 0.01, 0.01, 0.01, 0.01);
	}

	/*
	if (atoi(argv[1]) == 2) {
		double contact_rate = atof(argv[2]);
		scenario_swz(contact_rate);
	}
	*/
	if (atoi(argv[1]) == 2) {

		double cr1 = atof(argv[2]);
		double cr2 = atof(argv[3]);
		double cr3 = atof(argv[4]);
		double cr4 = atof(argv[5]);
		double cr5 = atof(argv[6]);
		scenario_swz(cr1,cr2,cr3,cr4,cr5);
	}

	if (atoi(argv[1]) == 4) {


		std::vector<double> contact_rates;
		std::vector<double> prevalences;
		std::vector<double> ks;
		double k1;
		double theta;

		//read_cr_from_file(contact_rates, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\temp_crs.txt");
		read_cr_from_file(contact_rates, argv[3]);
		read_k_from_file(ks, argv[4]);
		read_cr_from_file(prevalences, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\prev_list.txt");
		k1 = atof(argv[2]);
		theta = atof(argv[2]);

		scenario_tumikia(contact_rates,prevalences,ks,theta);
	}

	// burn in and intervention (mda)
	if (atoi(argv[1]) == 5) {


		std::vector<double> contact_rates;
		std::vector<double> prevalences;
		std::vector<double> ks; 
		double k1;
		double theta;

		//read_cr_from_file(contact_rates, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\temp_crs.txt");
		read_cr_from_file(contact_rates, argv[3]);
		read_k_from_file(ks, argv[4]);
		read_cr_from_file(prevalences, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\prev_list.txt");
		k1 = atof(argv[2]);
		theta = atof(argv[2]);

		scenario_tumikia_mda(contact_rates, prevalences, ks, theta);
	}

	
	// burn in and intervention (mda)
	if (atoi(argv[1]) == 6) {


		std::vector<double> contact_rates;
		std::vector<double> prevalences;
		std::vector<double> ks;
		double k1;
		double theta;

		//read_cr_from_file(contact_rates, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\temp_crs.txt");
		read_cr_from_file(contact_rates, argv[3]);
		read_k_from_file(ks, argv[4]);
		read_cr_from_file(prevalences, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\prev_list.txt");
		k1 = atof(argv[2]);
		theta = atof(argv[2]);

		scenario_tumikia_mda_long(contact_rates, prevalences, ks, theta);
	}

	if (atoi(argv[1]) == 7) {


		std::vector<double> contact_rates;
		std::vector<double> prevalences;
		std::vector<double> ks;
		double k1;
		double theta;
		std::string filename = argv[3];
		//read_cr_from_file(contact_rates, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\temp_crs.txt");
		read_cr_from_file(contact_rates, argv[3]);
		read_k_from_file(ks, argv[4]);
		read_cr_from_file(prevalences, "\\\\qdrive.dide.ic.ac.uk\\homes\\bcollyer\\parasim\\data\\prev_list.txt");
		k1 = atof(argv[2]);
		theta = atof(argv[2]);

		// For atoi, the input string has to start with a digit, so lets search for the first digit
		size_t i = 0;
		for (; i < filename.length(); i++) { if (isdigit(filename[i])) break; }

		// remove the first chars, which aren't digits
		filename = filename.substr(i, filename.length() - i);

		// convert the remaining text to an integer
		int fileid = atoi(filename.c_str());

		scenario_gesh(contact_rates, prevalences, ks, theta, fileid);
	}
	return 0;

	//scenario_swz();
	//return 0;
}

