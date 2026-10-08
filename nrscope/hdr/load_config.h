#ifndef LOAD_CONFIG_H
#define LOAD_CONFIG_H

#include <iostream>
#include <yaml-cpp/yaml.h>

#include "nrscope/hdr/radio_nr.h"

// using namespace std;

int get_nof_usrp(std::string file_name);
int load_config(std::vector<Radio>&, std::string file_name);
/* The sensing: block alone, into nrscope_sensing_args; load_config() calls it, and
   sensing_offline uses it to take the same settings without the radios. */
void load_sensing_config(const YAML::Node& config_yaml);

#endif