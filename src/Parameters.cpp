#include "Parameters.h"

double read_line(std::ifstream& myfile)
{
    std::string line;
    getline (myfile,line);
    //cout << "reading line :" << line << '\n';
    std::string token = line.substr(0, line.find(':'));
    //cout << "token :" << token << '\n';
    return stod(token);
};

Parameters read_from_file(std::string filename)
{
  Parameters parameters;
  std::ifstream myfile (filename);
  if (myfile.is_open())
  {
    parameters.nhosts = (int) read_line(myfile);
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
