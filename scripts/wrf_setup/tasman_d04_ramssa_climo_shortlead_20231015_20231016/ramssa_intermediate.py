import numpy as np
import pywinter.winter as pyw
import xarray as xr

def create_winter(fname,outname):

    f = xr.open_dataset(fname)
    lat = f.lat.values; lon=f.lon.values
    sst = f.analysed_sst.values
    winter_geo = pyw.Geo0(lat[0],lon[0],lat[1]-lat[0],lon[1]-lon[0])
    winter_sst = pyw.V2d("SST",sst)
    pyw.cinter("RAMSSA",outname,winter_geo,[winter_sst],"/g/data/eg3/ab4502/WRF/WPS/ramssa_data/")

if __name__ == "__main__":

    create_winter("/g/data/eg3/ab4502/WRF/WPS/ramssa_data/ramssa_climo.nc","2023-10-15_12")
    create_winter("/g/data/eg3/ab4502/WRF/WPS/ramssa_data/ramssa_climo.nc","2023-10-16_12")
