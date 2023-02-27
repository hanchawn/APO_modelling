#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=50gb
#SBATCH --account=CHEM007981

#100gb

source ~/.bashrc

conda info --envs

conda activate acrg

echo Current Date and Time is: `date +"%Y-%m-%d %T"`

year=2015
site='WAO'

month='all' #'none'

echo 'Year: '$year
echo 'Month: '$month
echo 'Site: '$Site

declare -a sectors=("domcom" "energyprod" "indcom" "indproc" "offshore" "othertrans" "roadtrans") #"solvents" "waste")
#declare -a sectors=("waste")

for sector in "${sectors[@]}"
do
    if [ $month == 'all' ]; then
        mth=1
        while [ $mth -le 12 ]
        do
            echo 'Creating timeseries for '$month
            python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_sectors.py $year $mth $site $sector

            ((mth++))
        done

        climatology='false'
        sector='ff-ukghg-'$sector

        echo 'Joining monthly time series'

        python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_join_months.py $year $site $sector $climatology
    
    else
        if [ $month == 'none' ]; then
            echo 'Creating timeseries for whole year'
        else
            echo 'Creating timeseries for month:' $month
        fi
        python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_sectors.py $year $month $site $sector

    fi
done