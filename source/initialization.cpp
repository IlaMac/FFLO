//
// Created by ilaria on 2019-11-13.
//
#include "initialization.h"
#include "robust_filesystem.h"

void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters){

    printf("initialize_Hparameters\n");
    int aux_dx, aux_dy;
    fs::path hp_init_file = directory_parameters / "HP_init.txt";
    if(fs::exists(hp_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(hp_init_file.c_str(), "r"))) {
            fscanf(fin, "%d" , &aux_dx);
            Hp.dx=(double) C_TWO_PI/aux_dx;
            fscanf(fin, "%d" , &aux_dy );
            Hp.dy=(double)C_TWO_PI/aux_dy;
            fscanf(fin, "%lf" , &Hp.b_low);
            fscanf(fin, "%lf" , &Hp.b_high);
            fscanf(fin, "%d" , &Hp.init);
            fclose(fin);
            //With this modification Hp.beta in not anymore part of the Hamiltonian parameters list
        }
    }else{
        Hp.dx=C_TWO_PI/4;
        Hp.dy=C_TWO_PI/4;
        Hp.b_low=0.1;
        Hp.b_high=0.6;
        Hp.init=1;
    }
    printf("dx %lf dy %lf blow %lf bhigh %lf init %d\n", Hp.dx, Hp.dy, Hp.b_low, Hp.b_high, Hp.init);

}

void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters){
    printf("initialize_MCparameters\n");

    fs::path mc_init_file = directory_parameters / "MC_init.txt";
    if(fs::exists(mc_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(mc_init_file.c_str(), "r"))) {
            fscanf(fin, "%d", &MCp.nmisu);
            fscanf(fin, "%d", &MCp.tau);
            fscanf(fin, "%d", &MCp.n_autosave);
            fscanf(fin, "%lf", &MCp.acc);
            fscanf(fin, "%lf", &MCp.a_T);
            fclose(fin);
        }
    }else{
        MCp.nmisu=100;
        MCp.tau=10;
        MCp.n_autosave=200; //not used now
        MCp.acc=0.5;
        MCp.a_T=0.5;
    }

}

void initialize_lattice(struct O2* Spin, const fs::path & directory_read, int RESTART, struct H_parameters &Hp){

    unsigned int i;
    fs::path psi_init_file = directory_read / "Psi_restart.bin";

    if(RESTART==1){
        fs::path psi_init_file = directory_read / "Psi_restart.bin";
    }
    else if (RESTART==2){
        fs::path psi_init_file = directory_read / "Psi_final.bin";
    }

    if((fs::exists(psi_init_file)) and (RESTART!=0)){

        FILE *fPsi= nullptr;
        if((fPsi=fopen(psi_init_file.c_str(), "r"))) {
            fread(Spin, sizeof(struct O2), N, fPsi);
            fclose(fPsi);
        }

    }else{
        if(Hp.init==0) {
            for (i = 0; i < N; i++) {
                Spin[i].r =1;
                Spin[i].t = 0.;
                polar_to_cartesian(Spin[i]);
            }
        }
        else if(Hp.init!=0) {
            for (i = 0; i < N; i++) {
                    Spin[i].r = 1;
                    Spin[i].t = rn::uniform_real_box(0, C_TWO_PI);
                    polar_to_cartesian(Spin[i]);
            }
        }

    }

}



void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp){
    int p;
    double  beta_low, beta_high, delta_beta;

    if(Hp.b_high>Hp.b_low){ //Paranoic check
        beta_low=Hp.b_low;
        beta_high=Hp.b_high;
    }else{
        beta_low=Hp.b_high;
        beta_high=Hp.b_low;
    }
    PTroot.beta.resize(PTp.np, 0.0);
    PTroot.All_Energies.resize(PTp.np, 0.0);
    PTroot.ind_to_rank.resize(PTp.np, 0);
    PTroot.rank_to_ind.resize(PTp.np, 0);
    delta_beta=(beta_high-beta_low)/(PTp.np-1);
    for(p=0; p<PTp.np; p++){
        PTroot.rank_to_ind[p]=p;
        PTroot.ind_to_rank[p]=p;
        PTroot.beta[p]=beta_low + p*delta_beta;
    }
}
