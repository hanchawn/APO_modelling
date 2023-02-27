# APO_modelling

Scripts for modelling atmospheric potential oxygen (APO) at UK sites, using inventory emissions estimates and the NAME atmospheric model.

We investigate the model sensitivities to its components: the biospheric and fossil fuel oxidative ratios, the fossil fuel fluxes, and the ocean fluxes.

## Making the figures

This is where each figure in the paper is made:

- Figure 1: calculating ff fluxes\
    powerpoint
- Figure 2: NAEI within EDGAR & NAEI minus EDGAR\
    SensitivityStudy/EDGAR_UKGHG_compare.ipynb (top half)
- Figure 3: ocean flux fields\
    SensitivityStudy/Ocean_Compare.ipynb
- Figure 4: component APO contributions\
    Timeseries/Plot_all.ipynb
- Figure 5: total APO model\
    Timeseries/Plot_all.ipynb
- Figure 6: APO model-obs correlations\
    Timeseries/Plot_all.ipynb
- Figure 7: $\alpha$ uncertainties\
    SensitivityStudy/Plot_all.ipynb
- Figure 8: ff uncertainties\
    SensitivityStudy/EDGAR_UKGHG_compare.ipynb (bottom half)
- Figure 9: $\Delta$ APO\
    Timeseries/bkg_subtracted.ipynb
- Figure 10: $\Delta$ APO correlations\
    Timeseries/bkg_subtracted.ipynb
- Figure 11: model ffCO2\
    ffCO2/ffCO2.ipynb
- Figure 12: obs ffCO2\
    ffCO2/ffCO2_obs.ipynb
- Figure 13: ffCO2 model-obs correlations\
    ffCO2/ffCO2_obs.ipynb

And in the supplement:

- Figure 2: total 2015 APO with bcs\
    ffCO2/ffCO2_obs.ipynb
- Figure 3: total modelled APO\
    Timeseries/Plot_all.ipynb
- Figure 4: obs vs model scatter\
    Timeseries/Plot_all.ipynb
- Figure 5: CO2 & O2 correlations\
    Timeseries/O2_CO2_plot.ipynb
- Figure 6: total modelled APO, 2021\
    Timeseries/Plot_all.ipynb
- Figure 7: ocean O2 flux correlations\
    Timeseries/Plot_all.ipynb
- Figure 8: fluxes - year vs climatology\
    SensitivityStudy/Ocean_yrVclim.ipynb
- Figure 9: REBS modelled APO baselines\
    ffCO2/ffCO2_obs.ipynb
- Figure 10: model APO deirved ffCO2\
    ffCO2/ffCO2_obs.ipynb

## Estimating the APO forward model

The forward model is run using several scripts:

- Timeseries/Timeseries_split.py\
    Runs the forward model for each APO component.
    Run using RunCodes/Timeseries/run_split_timeseries.sh.

    In the bash script we specify the:
    - year(s);
    - site(s);
    - sector: either bio, ff, or ocean;
    - climatology: true or false;
    - oxidative ratio source;
    - ff model: edgar or edgar-ukghg;
    - month: either all (runs each month separately before joining to reduce memory use), a single month, or year (runs whole year in one go).

    This outputs a netcdf file containing the forward model for each species required within the sector specified:
    - bio: o2 & co2     (not needed for APO but run to model CO2 and O2 separately)
    - ff: o2 & co2
    - ocean: co2, o2, n2    (co2 & o2 are calculated for all ocean models)

- Timeseries/Timeseries_join_months.py\
    Joins the forward models for each month produced above into a model for the whole year.
    This is automatically run when running RunCodes/Timeseries/run_split_timeseries.sh.

- Timeseries/Timeseries_join_months.ipynb\
    As above but to be run as needed.

- Timeseries/Timeseries_sectors.py\
    Runs the forward model for each ff sector separately.

- Timeseries/Timeseries_sectors_join.ipynb\
    Joins the above into one netcdf file with variables for each sector.

- Timeseries/Timeseries_bc.py\
    Runs the forward model for the boundary conditions.
    Run using RunCodes/Timeseries/run_bc_timeseries.sh.

    In the bash script we specify the:
    - year;
    - site(s);
    - climatology: true or false;
    - month: either all (runs each month separately before joining to reduce memory use), a single month, or year (runs whole year in one go).
    
    This outputs a netcdf file containing the baseline for:
    - delta APO
    - CO2
    - O2
    - delta O2/N2'

- Timeseries/JC_baseline_adjust.ipynb\
    Calculates an adjustment to the Jena Carboscope (JC) baseline:
    - compares the no-ocean APO model with the osb for each month;
    - subtracts the difference from the JC baseline;
    - saves the adjusted baseline to the bc timeseries netcdf file.

## Plotting the APO forward model

- Timeseries/APO_plot.ipynb
    - Plots the APO components and APO models
    - Calculates & plots correlations between the APO models and observations

- Timeseries/O2_CO2.ipynb\
    Plots the CO2 and O2 forward models

## Running sensitivity tests

- SensitivityStudy/BiosphericRatio.ipynb
    - Creates a normally distributed random sample of oxidative ratio values,
    - Calculates the total APO model for each values,
    - Finds the standard deviation etc,
    - Saves the 1 sigma uncertainty timeseries to a netcdf file.

- SensitivityStudy/EDGAR_UKGHG_compare.ipynb\
    Compares the EDGAR and UKGHG ffCO2 fields and APO mdoels

    Plots both fields and the difference between them,
    Plots the APO model from each field & obs, as well as the ff uncertainty from the Monte Carlo (see below)

- SensitivityStudy/FossilFuelMonteCarlo.ipynb
    - Randomly varies the ffCO2 within a range,
    - Calculates the total APO model for each ffCO2 timeseries,
    - Finds the standard deviation etc,
    - Saves the 1 sigma uncertainy to a netcdf file.

- SensitivityStudy/FossilFuelOR_sectors.ipynb
    - Creates a normally distributed random sample of oxidative ratio values,
    - Multiplies the ffCO2 by each oxidative ratio to get ffO2 for each sector,
    - Gets a total ffO2,
    - Finds the standard deviation etc,
    - Converts to uncertainty in APO,
    - Saves the 1 sigma uncertainy to a netcdf file.

- SensitivityStudy/Ocean_Compare.ipynb\
    Compare the ocean fields for each flux estimate

    Creates a figure with a field for each ocean o2 estimate and a footprint overlaid

- SensitivityStudy/Ocean_yrVclim.ipynb\
    Compare the ocean models using:
    - yearly data vs climatology,
    - daily vs monthly resolution.

- SensitivityStudy/Plot_all.ipynb\
    Creates a figure showing the total APO model and the oxidative ratio uncertainties

## Estimating the APO-derived ffCO2

- ffCO2/ffCO2_obs.ipynb\
    Estimate ffCO2 from the observations

- ffCO2/ffCO2.ipynb\
    Derive fCO2 from the APO model

- ffCO2/weighted_oxidative_ratio_timeseries.ipynb\
    Estimate a timeseries of footprint weighted oxidative ratios and APO:ffCO2 ratios


    
