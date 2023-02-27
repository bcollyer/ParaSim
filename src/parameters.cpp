#include "stdafx.h"
#include <fstream>
#include <string>
#include <iostream>
#include <cassert>
#include <iterator>
#include "parameters.h"

double read_line(std::ifstream& myfile)
{
	std::string line;
	std::string token;
	std::size_t found;
	std::string delim;
	double out = 0.0;

	getline(myfile, line);
	delim = ":";
	found = line.find(delim);
	token = line.substr(0, found);
	out = stod(token);
	//cout << "token :" << token << '\n';
	return out;
};

Parameters read_from_file(std::string filename)
{
	Parameters parameters;
	std::ifstream myfile(filename);
	if (myfile.is_open())
	{
		parameters.nhosts = (int)read_line(myfile);
		parameters.human_death_rate = read_line(myfile);
		parameters.parasite_death_rate = read_line(myfile);
		parameters.contact_rate = read_line(myfile);
		parameters.gamma = read_line(myfile);
		parameters.k = read_line(myfile);
		parameters.k_egg = read_line(myfile);
		parameters.fecundity = read_line(myfile);
		parameters.spat_kernel_magnitude = read_line(myfile);
	};
	return parameters;
};

void read_cr_from_file(std::vector<double> &x, const std::string &file_name)
{
	std::ifstream read_file(file_name);
	assert(read_file.is_open());

	std::copy(std::istream_iterator<double>(read_file), std::istream_iterator<double>(),
		std::back_inserter(x));

	read_file.close();
}


void read_k_from_file(std::vector<double> &x, const std::string &file_name)
{
	std::ifstream read_file(file_name);
	assert(read_file.is_open());

	std::copy(std::istream_iterator<double>(read_file), std::istream_iterator<double>(),
		std::back_inserter(x));

	read_file.close();
}


