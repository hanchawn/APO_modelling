#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=60gb
#SBATCH --account=CHEM007981

source ~/.bashrc

conda info --envs

conda activate acrg

echo Current Date and Time is: `date +"%Y-%m-%d %T"`

declare -a sites=("WAO" "HFD" "RGL")
year=2015
month='all'
climatology='false'

sector='bc'
echo $year
echo $month
echo $sector


for site in "${sites[@]}"
do
    echo 'Site: '$site
    if [ $month == 'all' ]; then
        mth=1
        while [ $mth -le 12 ]
        do
            echo 'Creating timeseries for month: '$mth
            python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_bc.py $year $mth $site

            ((mth++))
        done

        python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_join_months.py $year $site $sector $climatology
    else
        if [ $month == 'none' ]; then
            echo 'Creating timeseries for whole year'
        else
            echo 'Creating timeseries for month: '$month
        fi
        python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_bc.py $year $month $site

    fi
done