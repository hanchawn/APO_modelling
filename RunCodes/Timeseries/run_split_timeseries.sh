#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00
#SBATCH --mem=30gb
#SBATCH --account=CHEM007981

#22gb

source ~/.bashrc
conda info --envs
conda activate acrg

echo Current Date and Time is: `date +"%Y-%m-%d %T"`

declare -a years=("2014" "2015" "2016" "2017" "2018" "2019" "2020" "2021") # ("2015") #
declare -a sites=("HFD" "RGL" "WAO") # "TAC" "MHD"

sector='ocean' #'ffcomb' #'bc' #'bio' #

# whether to use averaged fluxes across multiple years
climatology='false' #'true' #

# oxidativeratio should be 'None' if using the sector-derived O2 fluxes
oxidativeratio='None' #'gridfed' #'gridfed-ukghg' #
# if ff_model='edgar_ukghg' the inventory names are not added to the output filename
ff_model='edgar-ukghg' #'edgar' #
bio_model='orchidee'

# if months='all' all months throughout the year are run and then joined
# if months='year' the whole year is run in one go - not recommended due to high memory use
month='all' #9 #'year' #

echo 'Years: '$years
echo 'Sites: '$sites
echo 'Sector: '$sector
if [ $sector == 'bio' ]; then
    echo 'Bio model: '$bio_model
elif [ $sector == 'ff' ]; then
    echo 'FF model: '$ff_model
    echo 'oxidative ratio: '$oxidativeratio
fi
echo 'Climatology: '$climatology

if [ $month == 'all' ]; then
    mths=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12")
elif [ $month == 'none' ]; then
    mths='none'
else
    mths=$month
fi

for year in "${years[@]}"
do
    echo 'Year: '$year
    echo 'Month: '$month
    for site in "${sites[@]}"
    do
        echo 'Site: '$site

        if [ $month == 'year' ]; then
            echo 'Creating timeseries for whole year'

            if [ $sector == 'bc' ]; then
                python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_bc.py $year $month $site
            else
                python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_split.py $year $month $site $sector $bio_model $ff_model $oxidativeratio $climatology
            fi
        
        else
            for mth in "${mths[@]}"
            do
                echo 'Creating timeseries for month: '$mth

                if [ $sector == 'bc' ]; then
                    python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_bc.py $year $mth $site
                else
                    python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_split.py $year $mth $site $sector $bio_model $ff_model $oxidativeratio $climatology
                fi

            done

            if [ $month == 'all' ]; then
                # join all months for the year
                if [ $sector == 'bio' ]; then
                    sec='bio_'$bio_model
                elif [ $sector == 'ff' ]; then
                    sec='ff_gridfed-ukghg'
                    if [ $ff_model == 'edgar-ukghg' ]; then
                        inv_str=''
                    else
                        inv_str='_'$ff_model
                    fi
                
                    if [ $oxidativeratio == 'None' ]; then
                        ox_ratio_str=''
                    else
                        if [ inv_str == '' ]; then
                            ox_ratio_str='_'$oxidativeratio
                        else
                            ox_ratio_str='-'$oxidativeratio
                        fi
                    fi
                    sec='ff'$inv_str$ox_ratio_str
                else
                    sec=$sector
                fi
                echo $sec
                python /user/home/vf20487/code/APO_modelling/Timeseries/Timeseries_join_months.py $year $site $sec $climatology
            fi
        fi
    done
done

        # else
        #     if [ $month == 'none' ]; then
        #         echo 'Creating timeseries for whole year'
        #     else
        #         echo 'Creating timeseries for month: '$month
        #     fi

        #     if [ $sector == 'bc' ]; then
        #         python /user/home/vf20487/code/Hannah/Timeseries/Timeseries_bc.py $year $month $site
        #     else
        #         python /user/home/vf20487/code/Hannah/Timeseries/Timeseries_split.py $year $month $site $sector $bio_model $ff_model $oxidativeratio $climatology
        #     fi
        # fi