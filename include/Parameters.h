#ifndef PARAMETERS_H
#define PARAMETERS_H


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
};

double read_line(std::ifstream& myfile);
Parameters read_from_file(std::string filename);

#endif
