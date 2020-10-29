//
// Created by ilaria on 2019-11-13.
//

#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "main.h"
#include "robust_filesystem.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#define C_TWO_PI (6.2831853071795864769252867665590058L)
#define C_PI (3.1415926535897932384626433832795029L)
#define Annealing (0) //To be implemented

//static constexpr double C_TWO_PI = 6.2831853071795864769252867665590058;



struct H_parameters{
    /*These values are externally given by an input file*/
    double dx;
    double dy;
    double b_low; //lowest beta for the parallel tempering procedure
    double b_high; // highest beta for the parallel tempering procedure
    int init;
};

struct PT_parameters{
    /*Parallel Tempering parameters*/
    int np;
    int rank;
    int root=0;
};

struct PTroot_parameters{
    /*Arrays root needs to handle the swaps*/
    std::vector <double> beta;
    std::vector <double> All_Energies;
    std::vector <int> ind_to_rank;
    std::vector <int> rank_to_ind;
};

struct MC_parameters{
    /*These values are externally given by an input file*/
    int tau; // estimate of the auto-correlation time
    int nmisu; //total number of independent measures
    int n_autosave; //frequency at which intermediate configuration are saved
    double acc; //It is a parameter which set the acceptance rate of the Metropolis sweep
    double a_T; //This parameter is temperature dependent. The larger is aT the bigger is the change in orientation proposed to the single spin.

};

void initialize_lattice(struct O2* Spin, const fs::path & directory_read, int RESTART, struct H_parameters &Hp);
void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters);
void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters);
void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp);

#endif //INITIALIZATION_H
