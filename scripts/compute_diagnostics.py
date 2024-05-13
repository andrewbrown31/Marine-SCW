from read_model_data.era5_read import read_era5_rt52, read_era5_cds
import datetime as dt
from diagnostic_driver import run_diagnostics

if __name__ == "__main__":

    start_lat = -44.525; end_lat = -9.975; start_lon = 111.975; end_lon = 160
    t1 = "2023101500"; t2="2023101623"
    domain = [start_lat,end_lat,start_lon,end_lon]
    time = [dt.datetime.strptime(t1,"%Y%m%d%H"),dt.datetime.strptime(t2,"%Y%m%d%H")]

    ta,dp,hur,hgt,terrain,p,ps,ua,va,uas,vas,tas,ta2d,cp,tp,wg10,cape,lon,lat,date_list = read_era5_cds(
	"/g/data/w40/ab4502/IN2023_V06/data/era5/era5_pl.nc", 
	"/g/data/w40/ab4502/IN2023_V06/data/era5/era5_sfc.nc", 
	domain,time,delta_t=1)

    #Run the diagnostic suite and save the output
    run_diagnostics(ta,hur,hgt,terrain,p,ps,ua,va,uas,vas,tas,ta2d,wg10,lon,lat,date_list,
					params="full",mdl_lvl=False,is_dcape=True,
					issave=True,out_name="era5",out_path="/g/data/w40/ab4502/IN2023_V06/data/era5/")
