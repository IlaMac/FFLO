#include "montecarlo.h"
#include "main.h"
#include "rng.h"
#include "class_tic_toc.h"

void metropolis( struct O2* Spin, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){

    int ix, iy, i, preso=0;
    unsigned int nn_ipx, ipx, nn_ipy, ipy, nn_ippy, ippy;
    unsigned int nn_imx, imx, nn_imy, imy, nn_immy, immy;
    double rand_epsilon=0, rand=0;
    double Epsilon=0, minusdeltaE2=0, rnorm=0, norm=0;
    double acc_rate=0.5, acc_T=0.;
    double E_old=0., E_new=0.;
    double Kx=0.;
    struct O2 Sp, newS;

    double invV=1./N;

    Sp.x=0.;
    Sp.y=0.;

    newS.x=0.;
    newS.y=0.;
    newS.t=0.;

    Kx=2*((Hp.dy*Hp.dy*Hp.dy*Hp.dy)/ (Hp.dx*Hp.dx));
    /*Procedo col metropolis in ordine, sito per sito, perdo la reversibilità del processo, ma guadagno in efficienza e velocità */
    for(iy=0; iy<Ly; iy++){
        for(ix=0; ix<Lx; ix++){
            i=ix+iy*Lx;
            ipx = (ix == Lx - 1 ? 0 : ix + 1);
            nn_ipx = ipx + Lx*iy;
            imx = (ix == 0? Lx-1 : ix - 1);
            nn_imx = imx + Lx*iy;
            ipy = (iy == Ly - 1 ? 0 : iy + 1);
            nn_ipy = ix + Lx*ipy;
            ippy = (ipy == Ly - 1 ? 0 : ipy + 1);
            nn_ippy = ix + Lx*ippy;
            imy = (iy == 0? Ly-1 : iy - 1);
            nn_imy = ix+ Lx*imy;
            immy = (imy == 0? Ly-1 : imy - 1);
            nn_immy = ix+ Lx*immy;

            /*Calcolo il nuovo spin*/
            rand_epsilon=rn::uniform_real_box(-1, 1);
            Epsilon=((double)MCp.a_T/MCp.acc)*rand_epsilon;
            O2perp(Spin[i],Sp);
            norm= 1./O2norm(Sp);
            Sp.x*=norm;
            Sp.y*=norm;
            Sp.x*=Epsilon;
            Sp.y*=Epsilon;

            O2sum(newS, Spin[i], Sp);
            rnorm=1./O2norm(newS);
            newS.x*=rnorm;
            newS.y*=rnorm;

            E_old=-(Kx*O2prod(Spin[i], Spin[nn_ipx]) + Kx*O2prod(Spin[i], Spin[nn_imx]) +
                    O2prod(Spin[i], Spin[nn_ipy]) + O2prod(Spin[i], Spin[nn_imy])
                    +0.5*sin(Spin[nn_ipy].t - Spin[i].t)*sin(Spin[i].t - Spin[nn_imy].t)
                    + 0.5*sin(Spin[i].t - Spin[nn_imy].t)*sin(Spin[nn_imy].t - Spin[nn_immy].t)
                    +0.5*sin(Spin[nn_ippy].t - Spin[nn_ipy].t)*sin(Spin[nn_ipy].t + Spin[i].t));

            E_new=-(Kx*O2prod(newS, Spin[nn_ipx]) + Kx*O2prod(newS, Spin[nn_imx]) +
                    O2prod(newS, Spin[nn_ipy]) + O2prod(newS, Spin[nn_imy])
                    +0.5*sin(Spin[nn_ipy].t - newS.t)*sin(newS.t - Spin[nn_imy].t)
                    + 0.5*sin(newS.t - Spin[nn_imy].t)*sin(Spin[nn_imy].t - Spin[nn_immy].t)
                    +0.5*sin(Spin[nn_ippy].t - Spin[nn_ipy].t)*sin(Spin[nn_ipy].t + newS.t));

            minusdeltaE2= Hp.dx*Hp.dy*(E_old-E_new);

            if(minusdeltaE2>0){
                Spin[i].x=newS.x;
                Spin[i].y=newS.y;
                preso++;
            }else{
                rand=rn::uniform_real_box(0, 1);
                if ( rand < exp((my_beta)*minusdeltaE2) ) {
                    Spin[i].x=newS.x;
                    Spin[i].y=newS.y;
                    preso++;
                }
            }
        }
    }
    MCp.a_T= (double)preso*invV;

}



