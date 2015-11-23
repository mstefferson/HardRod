// ParseParams.h
#ifndef _HARDROD_PARSEPARAMS_H_
#define _HARDROD_PARSEPARAMS_H_

#include <yaml-cpp/yaml.h>
#include <string.h>
#include "Parameters.h"
#include <iostream>
#include <cstring>
#include <string>

void parse_params(std::string param_file, system_params &params); 

#endif
