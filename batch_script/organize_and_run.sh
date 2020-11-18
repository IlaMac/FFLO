#!/bin/bash

BASEDIR=${HOME}/FFLO
SCRIPT_DIR=${BASEDIR}/batch_script

cd /tmp/

if [ ! -d ./SOutput_x_ilaria ]; then
   mkdir -p Output_x_ilaria
fi

#RESTART=0-> Restart from scratch
#RESTART=1-> Restart from interrupted run
#RESTART=2-> Restart from the previois final scenario

RESTART=0
pi_greek="3.14159265359"
LLIST="8 16 32 64 128"
############ Parameters of the Hamiltonian ---> HP_init.txt in a directory whose name contains the main parameters values##################
H_dx=$(echo '3.14159265359/8' | bc -l)
H_dy=$(echo '3.14159265359/8' | bc -l)
H_blow=0.1
H_bhigh=4.0
H_init=1 #If H_init=0: phases initialized to zero; H_init=1: phases initialized randomly

############ Parameters for the Monte Carlo simulations --> MC_init.txt#####################

Nmisu=200000
ntau=32
n_autosave=200
acc=0.5
a_T=0.5

for L in $LLIST; do

############Creation of the output folder and of the two files of initialization####################

cd ${BASEDIR}/Output_FFLO

if [ ! -d ./Sdx_2pi_8_dy_2pi_8 ]; then
   mkdir -p dx_2pi_8_dy_2pi_8
fi

cd dx_2pi_8_dy_2pi_8


if [ ! -d ./SL${L}_dx_2pi_8_dy_2pi_8_bmin${H_blow}_bmax${H_bhigh}_init${H_init} ]; then
   mkdir -p L${L}_dx_2pi_8_dy_2pi_8_bmin${H_blow}_bmax${H_bhigh}_init${H_init}
fi

OUTPUT=${BASEDIR}/Output_FFLO/dx_2pi_8_dy_2pi_8/L${L}_dx_2pi_8_dy_2pi_8_bmin${H_blow}_bmax${H_bhigh}_init${H_init}

cd /tmp/Output_x_ilaria

if [ ! -d ./Sdx_2pi_8_dy_2pi_8 ]; then
   mkdir -p dx_2pi_8_dy_2pi_8
fi

cd dx_2pi_8_dy_2pi_8

if [ ! -d ./SL${L}_dx_2pi_8_dy_2pi_8_bmin${H_blow}_bmax${H_bhigh}_init${H_init} ]; then
   mkdir -p L${L}_dx_2pi_8_dy_2pi_8_bmin${H_blow}_bmax${H_bhigh}_init${H_init}
fi


OUTPUT_TEMP=/tmp/Output_x_ilaria/dx_2pi_8_dy_2pi_8/L${L}_dx_2pi_8_dy_2pi_8_bmin${H_blow}_bmax${H_bhigh}_init${H_init}

cd ${OUTPUT}

#THE ORDER OF WRITING DOES MATTER
echo $H_dx > HP_init.txt
echo $H_dy >> HP_init.txt
echo $H_blow >> HP_init.txt
echo $H_bhigh >> HP_init.txt
echo $H_init >> HP_init.txt

#THE ORDER OF WRITING DOES MATTER
echo $Nmisu > MC_init.txt
echo $ntau >> MC_init.txt
echo $n_autosave>> MC_init.txt
echo $acc >> MC_init.txt
echo $a_T >> MC_init.txt

#################Creation of the submit_runs script#########################

jobname="L${L}_dx_2pi_8_dy_2pi_8_bmin${H_blow}_bmax${H_bhigh}_init${H_init}"
nnodes=2
ntasks=64 #parallel tempering over ntasks temperatures

#I create ntasks folder: one for each rank.

cd ${OUTPUT}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${OUTPUT_TEMP}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${SCRIPT_DIR}
DIR_PAR="${OUTPUT}"
DIR_PAR_TEMP="${OUTPUT_TEMP}"

#SEED= If I want to repeat exactly a simulation I could initialize the random number generator exactly at the same way

EXECUTE_DIR="../build/Release"

#SBATCH --nodes=${nnodes}               # Number of nodes

echo "#!/bin/bash
#SBATCH --job-name=${jobname}          # Name of the job
#SBATCH --time=7-00:00:00               # Allocation time
#SBATCH --mem-per-cpu=2000              # Memory per allocated cpu
#SBATCH --nodes=${nnodes}               # Number of nodes
#SBATCH --ntasks=${ntasks}
#SBATCH --output=${DIR_PAR}/logs/log_${jobname}.o
#SBATCH --error=${DIR_PAR}/logs/log_${jobname}.e

srun ${EXECUTE_DIR}/FFLO ${L} ${DIR_PAR} ${DIR_PAR} ${RESTART} &> ${DIR_PAR}/logs/log_${jobname}.o

" >  submit_run

#Submission of the work --> sbatch submit_runs

mkdir -p ${DIR_PAR}/logs

sbatch submit_run

done
