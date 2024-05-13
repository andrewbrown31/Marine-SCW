import xarray as xr

if __name__ == "__main__":

	base = "https://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L4/RAMSSA/%4d/%4d10%02d120000-ABOM-L4_GHRSST-SSTfnd-RAMSSA_09km-AUS-v02.0-fv01.0.nc?time[0:1:0],lat[300:1:540],lon[960:1:1320],analysed_sst[0:1:0][300:1:540][960:1:1320]"
	url_list = []
	for d in range(2006,2023):
	    for dd in range(11,22):
		url_list.append(base % (d,d,dd))
	ramssa_climo = xr.open_mfdataset(url_list)

	ramssa_climo_mean = ramssa_climo.analysed_sst.mean("time").persist()
	ramssa_climo_mean.to_netcdf("/g/data/eg3/ab4502/WRF/WPS/ramssa_data/ramssa_climo.nc")
