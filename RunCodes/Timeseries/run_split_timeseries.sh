#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=20gb

#100gb

source ~/.bashrc
conda info --envs
conda activate acrg

declare -a years=("2014") # "2015" "2017" "2018" "2019" "2020"
declare -a sites=("WAO") #"RGL" "HFD" "MHD")

sector='ff' #'bc' #'ocean' #'bio' #

# whether to use averaged fluxes across multiple years
climatology='false' #'true' #

# oxidativeratio should be 'None' if using the sector-derived O2 fluxes
oxidativeratio='gridfed' #'None' #'gridfed-ukghg'
# if ff_model='edgar_ukghg' the inventory names are not added to the output filename
ff_model='edgar-ukghg'
bio_model='orchidee'

# if months='all' all months throughout the year are run and then joined
# if months='year' the whole year is run in one go - not recommended due to high memory use
month='all' #("5") #1 #'year' #

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
                if [ $sector == 'ff' ]; then
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
                    sec='ff'$inv_str$ox_ratio_str
                elif [ $sector == 'bio' ]; then
                    sec='bio_'$bio_model
                else
                    sec=$sector
                fi
                python /user/home/vf20487/code/Hannah/Timeseries/Timeseries_join_months.py $year $site $sec $climatology
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