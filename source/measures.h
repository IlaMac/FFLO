//
// Created by ilaria on 2019-11-16.
//

#ifndef MEASURES_H
#define MEASURES_H

#include "main.h"
#include <fstream>

struct Measures{

/*******THE ORDER SHOULD BE THE SAME OF THE H5T INSERT REGISTRATION**************/

    double E=0; //Energy
    double M[2]={0};
    //superfluid stiffness Js= <jd>- \beta*[<ip^2> - <ip>^2]
    double jd[2]={0}; //diamagnetic term: jdx=jd[0]; jdy=jd[1]
    double ip[2]={0}; //paramagnetic current: ipx=ip[0]; ipy=ip[1]
    //vortex density
    double rho_vplus=0.;
    double rho_vminus=0.;

    int my_rank = 0;
    void reset(){
        *this = Measures();
    }
};


void save_lattice(struct O2* Spin, const fs::path & directory_write, std::string configuration);
void all_measures(struct Measures &mis, struct H_parameters &Hp, double my_beta, struct O2* Spin);

#endif //MEASURES_H
