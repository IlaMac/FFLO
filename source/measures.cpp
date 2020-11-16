//
// Created by ilaria on 2019-11-16.
//

#include "measures.h"

void all_measures(struct Measures &mis, struct H_parameters &Hp, double my_beta, struct O2* Spin) {

    unsigned ix, iy, i;
    unsigned ipx, ipy, imy, nn_ipx, nn_ipy, nn_imy;
    double Kx;
    double phi_1=0., phi_2=0., phi_3=0., phi_4=0.;
    double n=0.;
    double invV=1./N;

    Kx=2*(power(Hp.dy,4)/ power(Hp.dx,2));

    /*Ridefinisco tutti gli angoli compresi tra -pi e pi*/
    for(ix=0; ix<Lx; ix++){
        for (iy=0; iy<Ly; iy++) {
            i=ix+iy*Lx;
            cartesian_to_polar(Spin[i]);
            if(Spin[i].t >C_PI){ Spin[i].t-= C_TWO_PI;}
        }
    }


    for (iy=0; iy<Ly; iy++) {
        for (ix=0; ix<Lx; ix++) {

            i=ix+iy*Lx;
            ipx=mod(ix+1,Lx);
            nn_ipx = ipx + Lx*iy;
            ipy=mod(iy+1,Ly);
            nn_ipy = ix + Lx*ipy;
            imy = mod(iy-1,Ly);
            nn_imy = ix+ Lx*imy;
            /*Energy*/
            mis.E+=Hp.dx*Hp.dy*( Kx*(1. -O2prod(Spin[i], Spin[nn_ipx])) + (1. -O2prod(Spin[i], Spin[nn_ipy]))
                    - 0.5*sin(-Spin[i].t + Spin[nn_ipy].t)*sin(-Spin[nn_imy].t + Spin[i].t) );

            /*Magetization*/
            mis.M[0]+=Spin[i].x;
            mis.M[1]+=Spin[i].y;

            /*Helicity (Stiffness)*/

            //Ix=Kx*sin(\phi_{i+x} - \phi_i)
            mis.ip[0]+= Kx*sin(-Spin[i].t +Spin[nn_ipx].t);
            //Jdx=Kx*cos(\phi_{i+x} - \phi_i)
            mis.jd[0]+= Kx*O2prod(Spin[i], Spin[nn_ipx]);

            //Ipy= Ky*(sin(\phi_{i+y} - \phi_i) - 0.5*sin(\phi_{i+y} - \phi_{i-y})
            mis.ip[1]+= (sin(-Spin[i].t + Spin[nn_ipy].t) - 0.5*sin(-Spin[nn_imy].t + Spin[nn_ipy].t));
            //Jdy= Ky*(cos(\phi_{i+y} - \phi_i) - cos(\phi_{i+y} - \phi_{i-y})
            mis.jd[1]+=(O2prod(Spin[i], Spin[nn_ipy]) - O2prod(Spin[nn_imy], Spin[nn_ipy]));

            /*Numero vortici*/
            /*Circuitazione*/
            n=0.;
            phi_1= Spin[i].t -Spin[nn_ipx].t;
            if(phi_1< -C_PI){ phi_1+= C_TWO_PI; }
            if(phi_1> C_PI){ phi_1-= C_TWO_PI;}
            phi_2= Spin[nn_ipx].t -Spin[ipx+imy*Lx].t;
            if(phi_2< -C_PI){ phi_2+= C_TWO_PI; }
            if(phi_2> C_PI){ phi_2-= C_TWO_PI;}
            phi_3= Spin[ipx+imy*Lx].t -Spin[nn_imy].t;
            if(phi_3< -C_PI){ phi_3+= C_TWO_PI; }
            if(phi_3> C_PI){ phi_3-= C_TWO_PI;}
            phi_4= Spin[nn_imy].t -Spin[i].t;
            if(phi_4< -C_PI){ phi_4+= C_TWO_PI; }
            if(phi_4> C_PI){ phi_4-= C_TWO_PI;}

            n=(phi_1+phi_2+phi_3+phi_4)/C_TWO_PI;

            if(n>0.1){
                mis.rho_vplus+=n;
            }
            else if (n<-0.1){
                mis.rho_vminus+=n;
            }
        }
    }


    mis.E*= invV;
    mis.M[0]*=invV;
    mis.M[1]*=invV;
    mis.jd[0]*=invV;
    mis.jd[1]*=invV;
    mis.ip[0]*=invV;
    mis.ip[1]*=invV;
    mis.rho_vplus*= invV;
    mis.rho_vminus*= invV;

    return;

}


void save_lattice(struct O2* Spin, const fs::path & directory_write, std::string configuration){

    std::string sPsi;
    sPsi= std::string("Psi_")+ configuration + std::string(".bin");
    fs::path psi_init_file = directory_write / sPsi;

    FILE *fPsi= nullptr;
    unsigned int i=0;

    if((fPsi=fopen(psi_init_file.c_str(), "w"))) {
        fwrite(Spin, sizeof(struct O2), N, fPsi);
        fclose(fPsi);
    }


}

