#pragma once
void scenario_square(double contact_rate);
void scenario_swz(double cr1, double cr2, double cr3, double cr4, double cr5);
void scenario_tumikia(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> prevalence, double theta);
void scenario_tumikia_mda(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> prevalence, double theta);
void scenario_tumikia_mda_long(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> prevalence, double theta);
void scenario_gesh(std::vector<double> contact_rates, std::vector<double> prevalences, std::vector<double> prevalence, double theta, double fileid);