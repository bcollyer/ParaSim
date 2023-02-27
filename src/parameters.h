#pragma once
#include <vector>

struct Parameters
{
	int nhosts;
	double human_death_rate;
	double parasite_death_rate;
	double contact_rate;
	double gamma;
	double k;
	double k_egg;
	double fecundity;
	double spat_kernel_magnitude;
	double theta;
	double presac_contact;
	double sac_contact;
	double adult_contact;
};

double read_line(std::ifstream& myfile);
Parameters read_from_file(std::string filename);
void read_cr_from_file(std::vector<double> &x, const std::string &file_name);
void read_k_from_file(std::vector<double> &x, const std::string &file_name);
