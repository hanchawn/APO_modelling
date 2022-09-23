#!/bin/sh
#PBS -l select=1:ncpus=4:mem=80gb
#PBS -l walltime=72:00:00

echo Current Date and Time is: `date +"%Y-%m-%d %T"`
echo Running on host "$(hostname)"
echo Directory is "$(pwd)"

source ~/.bashrc

conda info --envs

conda activate acrg

echo Varying the fossil fuel oxidative ratio

python /user/home/vf20487/code/APO_modelling/SensitivityStudy/FossilFuelRatio.py
