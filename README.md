## The impact of sea surface temperatures on a severe convective wind event in the East Australian Current system

This repository contains scripts and notebooks for analysis and modelling of a severe wind event in the Tasman Sea. This analysis is intended for inclusion in a research thesis (in preparation). The marine severe wind event was observed during the RV Investigator voyage [Understanding eddy interactions and their impacts in the EAC](https://www.csiro.au/en/about/facilities-collections/mnf/voyages-schedules/voyages/2023/october/in2023_v06).

## Log
1. Version 4.3 of the Weather Research and Forecasting (WRF) model is [setup on the Gadi HPC system](https://github.com/coecms/WRF) and run for a control case using [this script](scripts/wrf_setup/tasman_d04_ramssa_shortlead_20231015_20231016/run_tasman_d04) (noting that one restart is also required using [this script](scripts/wrf_setup/tasman_d04_ramssa_shortlead_20231015_20231016/run_tasman_d04_restart)). Pre-processing and runtime WRF namelists for WRF simulations are available in the [namelists](scripts/wrf_setup/namelists) directory.
2. Three experiments are performed using WRF, including with [no latent heating from microphysics](scripts/wrf_setup/tasman_d04_ramssa_shortlead_no_mp_heating_20231015_20231016), a [warm SST perturbation](scripts/wrf_setup/tasman_d04_ramssa_shortlead_3Kplus_20231015_20231016), and a [cold SST perturbation](scripts/wrf_setup/tasman_d04_ramssa_shortlead_3Kminus_20231015_20231016). [Several other WRF experiments](scripts/wrf_setup) were performed but are not presented in the research thesis.
3. Vertical cross sections are taken along several transects from WRF model output using [this script](scripts/transects/compute_transects.py). Options for transect bounds and orientation are determined by manual inspection of WRF output. These cross sections are [averaged and visualised](scripts/pub_notebooks/control_transects.ipynb), and used to [compare density current statistics](scripts/pub_notebooks/compare_transects.ipynb) between WRF simulations.
4. WRF output is compared with [wind speed](scripts/pub_notebooks/time_series.ipynb) and [radar observations](scripts/pub_notebooks/plan_view.ipynb) from oboard the RV Investigator and [cloud-top temperatures from the Himawari satellite](scripts/pub_notebooks/convection_d03.ipynb).
5. Vertical profiles of temperature, wind, and moisture from the WRF model are compared with [radiosonde observations](scripts/pub_notebooks/plot_soundings.ipynb)
