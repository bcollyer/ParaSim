#include "Parameters.h"

double read_line(std::ifstream& myfile)
{
    std::string line;
    std::string token;
    std::size_t found;
    std::string delim;
    double out = 0.0;

    getline (myfile,line);
    //cout << "reading line :" << line << '\n';
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
