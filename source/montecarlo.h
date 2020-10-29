#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <iostream>
#include <cstring>
#include "o2.h"

void metropolis(struct O2* Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
#endif