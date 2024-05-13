import tqdm
import datetime as dt
import wrf
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np
import xarray as xr
import metpy.calc as mpcalc

def plot_transects(path,f,xlim,ylim,no_of_transects=30,d_km=50,grad_theta_thresh=2):

    #Get variables with wrf-python
    print("Plotting transects...")
    wrflist = nc.Dataset(path+f)
    theta = wrf.getvar(wrflist, "theta").isel(bottom_top=0)
    grad_theta = calc_grad_theta(theta)
    wspd_wdir = wrf.getvar(wrflist,"wspd_wdir10")
    wspd = wspd_wdir.isel(wspd_wdir=0)
    wdir = wspd_wdir.isel(wspd_wdir=1)

    line_inds = (grad_theta.values >= grad_theta_thresh) \
                & (grad_theta.XLONG >= xlim[0]) & (grad_theta.XLONG <= xlim[1]) \
                & (grad_theta.XLAT >= ylim[0]) & (grad_theta.XLAT <= ylim[1])
    #step = int(np.round((line_inds.sum() / no_of_transects).values)) #Plot around 20 transects along the line of grad_theta >= 2
    x = theta.west_east*wrflist.DX
    y = theta.south_north*wrflist.DY
    x,y=np.meshgrid(x,y)
    steps = np.linspace(0,x[line_inds].shape[0]-1,no_of_transects).round().astype(int)
    line_wdir = 270 - wdir.rolling({"south_north":11,"west_east":11}).mean().values[line_inds][steps]
    line_x_cart = x[line_inds][steps]
    line_y_cart = y[line_inds][steps]    
    print(len(line_x_cart)," transects")
    start_points = []
    end_points = []
    angles = []
    D = (d_km/2) * 1000 
    for i in range(len(line_x_cart)):

        #Angle of transect
        angle = line_wdir[i]
        angles.append(angle)

        #Start and end points
        start_point = [line_y_cart[i] + D * np.sin(np.deg2rad(angle-180)), line_x_cart[i] + D * np.cos(np.deg2rad(angle-180))]
        end_point = [line_y_cart[i] + D * np.sin(np.deg2rad(angle)), line_x_cart[i] + D * np.cos(np.deg2rad(angle))]

        #Transform back into lat-lon coords
        start_ind = [np.argmin(np.abs(start_point[0] - theta.west_east.values*wrflist.DX)), 
                            np.argmin(np.abs(start_point[1] - theta.south_north.values*wrflist.DY))]
        end_ind = [np.argmin(np.abs(end_point[0] - theta.west_east.values*wrflist.DX)), 
                            np.argmin(np.abs(end_point[1] - theta.south_north.values*wrflist.DY))]    
        start_point_lon = grad_theta.XLONG.values[start_ind[0],start_ind[1]]
        start_point_lat = grad_theta.XLAT.values[start_ind[0],start_ind[1]]
        end_point_lon = grad_theta.XLONG.values[end_ind[0],end_ind[1]]
        end_point_lat = grad_theta.XLAT.values[end_ind[0],end_ind[1]]

        #Keep track of lat lon coords
        start_points.append([start_point_lon,start_point_lat])
        end_points.append([end_point_lon,end_point_lat])    

    #Plot
    plt.figure(figsize=[10,6])
    ax = plt.axes(projection=ccrs.PlateCarree())
    for i in range(len(line_x_cart)):
	    start_point_lon,start_point_lat = start_points[i]
	    end_point_lon,end_point_lat = end_points[i]
	    ax.plot([start_point_lon,end_point_lon],[start_point_lat,end_point_lat],lw=2,color="grey",ls="--")
	    ax.plot(start_point_lon,start_point_lat,marker="o",color="grey",mec="k")
	    ax.text(end_point_lon,end_point_lat,str(i),fontdict={"color":"k","fontweight":"bold","fontsize":"medium"})
    line_x = grad_theta.XLONG.values[line_inds][steps]
    line_y = grad_theta.XLAT.values[line_inds][steps]
    uvmet10 = wrf.getvar(wrflist, "uvmet10")
    u10 = uvmet10.isel(u_v=0).drop_vars("u_v")
    v10 = uvmet10.isel(u_v=1).drop_vars("u_v")
    line_u = u10.values[line_inds][steps]
    line_v = v10.values[line_inds][steps]
    extent = ax.get_extent()

    wspd.plot(ax=ax,x="XLONG",y="XLAT",cmap="Reds",levels=np.linspace(2,22,11))
    uv_ds = xr.Dataset({"u":u10,"v":v10}).coarsen(dim={"south_north":20,"west_east":20},boundary="trim").mean()
    uv_ds.plot.quiver("XLONG","XLAT","u","v",ax=ax,scale=200,width=0.002)  
    xr.plot.contour(grad_theta,levels=[grad_theta_thresh],colors=["tab:blue"],x="XLONG",y="XLAT")
    #ax.set_extent(extent)
    ax.gridlines(draw_labels=["left","bottom"],ls=":")
    plt.savefig("/g/data/w40/ab4502/WRF_simulations/transects/"+\
		    path.split("/")[-2] + "_" + f.split("/")[-1].split("_")[-1] + ".jpeg",dpi=300)

def convert_xyloc(ds):

    xy_loc = ds.xy_loc.values
    xy_lat = np.zeros(xy_loc.shape)
    xy_lon = np.zeros(xy_loc.shape)
    for i in range(xy_loc.shape[0]):
        for j in range(xy_loc.shape[1]):    
            xy_lat[i,j] = xy_loc[i,j].lat
            xy_lon[i,j] = xy_loc[i,j].lon        

    ds = ds.drop_vars("xy_loc")
    ds = ds.assign(xy_lat=(["transect","cross_line_idx"],xy_lat), xy_lon=(["transect","cross_line_idx"],xy_lon))
    return ds

def drop_attrs(ds):

    for d in ds.data_vars:
        try:
            del ds[d].attrs["projection"]
        except:
            pass

    return ds.drop_vars(["XTIME","wspd_wdir"])

def drop_projection(ds):
    del ds.attrs["projection"]
    return ds

def calc_grad_theta(theta):
    theta_np = theta.values
    dx, dy = mpcalc.lat_lon_grid_deltas(theta.XLONG.values, theta.XLAT.values)
    dtheta_dy, dtheta_dx = mpcalc.gradient(theta_np, deltas=[dy.to("km"), dx.to("km")])
    grad_theta = dtheta_dy + dtheta_dx
    grad_theta = xr.DataArray(grad_theta.data,dims=["south_north","west_east"])
    grad_theta["XLONG"] = theta.XLONG
    return grad_theta

def latlon_dist(lat, lon, lats, lons):

        #Calculate great circle distance (Harversine) between a lat lon point (lat, lon) and a list of lat lon
        # points (lats, lons)

        R = 6373.0

        lat1 = np.deg2rad(lat)
        lon1 = np.deg2rad(lon)
        lat2 = np.deg2rad(lats)
        lon2 = np.deg2rad(lons)

        dlon = lon2 - lon1
        dlat = lat2 - lat1

        a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

        return (R * c)
    
def cross_line_km(cross,lats,lons):
    
    xy_lats = [xy.lat for xy in cross.xy_loc.values]
    xy_lons = [xy.lon for xy in cross.xy_loc.values]
    xy = np.round(latlon_dist(lats[0],lons[0],xy_lats,xy_lons),1)#.astype(str)

    return cross.assign_coords({"cross_line_idx":xy})

def get_base_state(path):
    #Get base states for perturbation calculations
    #Define base state as the horzontal average of all hourly times
    
    initial_time = dt.datetime(2023,10,15,15)
    times = []
    while initial_time <= dt.datetime(2023,10,16,12):
        times.append(initial_time.strftime("%Y-%m-%d_%H:%M:%S"))
        initial_time = initial_time+dt.timedelta(hours=1)
        
    path = path + "highres_d04_"
    
    base = [nc.Dataset(path+times[t]) for t in range(len(times))]
    T_base = wrf.getvar(base,"T",timeidx=wrf.ALL_TIMES).mean(("west_east","south_north","Time"))
    theta_base = wrf.getvar(base,"theta",timeidx=wrf.ALL_TIMES).mean(("west_east","south_north","Time"))
    p_base = wrf.getvar(base,"p",timeidx=wrf.ALL_TIMES).mean(("west_east","south_north","Time"))
    qv_base = wrf.getvar(base,"QVAPOR",timeidx=wrf.ALL_TIMES).mean(("west_east","south_north","Time"))
    
    return T_base, theta_base, p_base, qv_base

def get_transects(path,f,xlim,ylim,no_of_transects=30,d_km=50,top=12000,mean_motion_top=6000,grad_theta_thresh=2):

    #Get variables with wrf-python
    print("Doing cross-sections...")
    wrflist = nc.Dataset(path+f)
    #2d for defining transects
    theta = wrf.getvar(wrflist, "theta").isel(bottom_top=0)
    grad_theta = calc_grad_theta(theta)
    wspd_wdir = wrf.getvar(wrflist,"wspd_wdir10")
    wspd = wspd_wdir.isel(wspd_wdir=0)
    wdir = wspd_wdir.isel(wspd_wdir=1)
    #3d for taking cross-sections
    T_base, theta_base, p_base, qv_base = get_base_state(path)
    heights = wrf.getvar(wrflist,"z")
    wspd_wdir_3d = wrf.getvar(wrflist, "wspd_wdir")
    wspd_3d = wspd_wdir_3d.isel(wspd_wdir=0)    
    theta_3d = wrf.getvar(wrflist, "theta")
    theta_pert_3d = theta_3d - theta_base
    T_3d = wrf.getvar(wrflist,"T")
    Tv_3d = wrf.getvar(wrflist,"tv")
    thetae_3d = wrf.getvar(wrflist,"theta_e")
    P_3d = wrf.getvar(wrflist,"p")
    p_pert_3d = P_3d - p_base
    U_3d = wrf.destagger(wrf.getvar(wrflist,"U"),2,meta=True)
    V_3d = wrf.destagger(wrf.getvar(wrflist,"V"),1,meta=True)
    QR_3d = wrf.getvar(wrflist,"QRAIN")
    QC_3d = wrf.getvar(wrflist,"QCLOUD")
    QV_3d = wrf.getvar(wrflist,"QVAPOR")
    QI_3d = wrf.getvar(wrflist,"QICE")
    W_3d = wrf.destagger(wrf.getvar(wrflist,"W"),0,meta=True)
    dbz_3d = wrf.getvar(wrflist,"dbz")
    #Save dbz for later
    drop_projection(dbz_3d).to_netcdf(path+"dbz_"+f)
    B = 9.8 * ( (theta_pert_3d / theta_base) + 0.61*(QV_3d - qv_base) - (QC_3d + QR_3d + QI_3d) )
    rho = mpcalc.density(P_3d,Tv_3d,0)
    #Add coordinate info to destaggered variables
    U_3d["XLONG"] = T_3d["XLONG"]
    V_3d["XLONG"] = T_3d["XLONG"]
    W_3d["XLONG"] = T_3d["XLONG"]
    
    #Set up transect points along squall line, and angle of transects.
    #Defined by number of transects (no_of_transects) along a line of grad_theta greater than 2 K/km at the lowest model level.
    #The angle of the transects are defined by the local 10 m wind direction at each point.
    #Note that a rolling 5 km mean is applied along the lat and lon direction before calculating the local wind
    line_inds = (grad_theta.values >= grad_theta_thresh) \
                & (grad_theta.XLONG >= xlim[0]) & (grad_theta.XLONG <= xlim[1]) \
                & (grad_theta.XLAT >= ylim[0]) & (grad_theta.XLAT <= ylim[1])
    #step = int(np.round((line_inds.sum() / no_of_transects).values)) #Plot around 20 transects along the line of grad_theta >= 2
    
    #Calculate the start and end points for each transect
    #Set up x and y grid in km
    x = theta.west_east*wrflist.DX
    y = theta.south_north*wrflist.DY
    x,y=np.meshgrid(x,y)
    steps = np.linspace(0,x[line_inds].shape[0]-1,no_of_transects).round().astype(int)
    line_x_cart = x[line_inds][steps]
    line_y_cart = y[line_inds][steps]    
    line_wdir = 270 - wdir.rolling({"south_north":11,"west_east":11}).mean().values[line_inds][steps]
    print(len(line_x_cart)," transects")
    start_points = []
    end_points = []
    angles = []
    D = (d_km/2) * 1000 
    for i in range(len(line_x_cart)):

        #Angle of transect
        angle = line_wdir[i]
        angles.append(angle)

        #Start and end points
        start_point = [line_y_cart[i] + D * np.sin(np.deg2rad(angle-180)), line_x_cart[i] + D * np.cos(np.deg2rad(angle-180))]
        end_point = [line_y_cart[i] + D * np.sin(np.deg2rad(angle)), line_x_cart[i] + D * np.cos(np.deg2rad(angle))]

        #Transform back into lat-lon coords
        start_ind = [np.argmin(np.abs(start_point[0] - theta.west_east.values*wrflist.DX)), 
                            np.argmin(np.abs(start_point[1] - theta.south_north.values*wrflist.DY))]
        end_ind = [np.argmin(np.abs(end_point[0] - theta.west_east.values*wrflist.DX)), 
                            np.argmin(np.abs(end_point[1] - theta.south_north.values*wrflist.DY))]    
        start_point_lon = grad_theta.XLONG.values[start_ind[0],start_ind[1]]
        start_point_lat = grad_theta.XLAT.values[start_ind[0],start_ind[1]]
        end_point_lon = grad_theta.XLONG.values[end_ind[0],end_ind[1]]
        end_point_lat = grad_theta.XLAT.values[end_ind[0],end_ind[1]]

        #Keep track of lat lon coords
        start_points.append([start_point_lon,start_point_lat])
        end_points.append([end_point_lon,end_point_lat])    
        
    #Now take cross-sections through the transects, and average
    theta_crosses = []
    wspd_crosses = []
    theta_pert_crosses = []
    p_pert_crosses = []
    u_crosses = []
    sru_crosses = []
    w_crosses = []
    dbz_crosses = []
    b_crosses = []
    p_pert_crosses = []
    rho_crosses = []
    tv_crosses = []
    for i in tqdm.tqdm(range(len(line_x_cart))):

        start_point = wrf.CoordPair(lat=start_points[i][1],lon=start_points[i][0])
        end_point = wrf.CoordPair(lat=end_points[i][1],lon=end_points[i][0])    

        theta_cross = wrf.vertcross(theta_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        theta_crosses.append( cross_line_km(theta_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))  

        wspd_cross = wrf.vertcross(wspd_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        wspd_crosses.append( cross_line_km(wspd_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))    

        theta_pert_cross = wrf.vertcross(theta_pert_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        theta_pert_crosses.append( cross_line_km(theta_pert_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))  

        u_cross = wrf.vertcross(U_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        v_cross = wrf.vertcross(V_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
	#https://en.wikipedia.org/wiki/Rotation_of_axes_in_two_dimensions
        along_transect_wind = np.cos(np.deg2rad(angles[i])) * u_cross + np.sin(np.deg2rad(angles[i])) * v_cross
        umean = along_transect_wind.sel({"vertical":slice(0,mean_motion_top)}).mean()
        u_crosses.append(cross_line_km(along_transect_wind,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))
        sru_crosses.append(cross_line_km(along_transect_wind - umean,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))

        w_cross = wrf.vertcross(W_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        w_crosses.append( cross_line_km(w_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))  

        dbz_cross = wrf.vertcross(dbz_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        dbz_crosses.append( cross_line_km(dbz_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))  

        b_cross = wrf.vertcross(B,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        b_crosses.append( cross_line_km(b_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))  

        p_cross = wrf.vertcross(P_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        p_pert_cross = ((p_cross - p_cross.mean("cross_line_idx"))/100)
        p_pert_crosses.append(cross_line_km(p_pert_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))

        rho_cross = wrf.vertcross(rho,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        rho_crosses.append( cross_line_km(rho_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))
	
        tv_cross = wrf.vertcross(Tv_3d,heights,wrfin=wrflist,levels=np.linspace(50,top,100),start_point=start_point,end_point=end_point,latlon=True)
        tv_crosses.append( cross_line_km(tv_cross,[start_point.lat,end_point.lat],[start_point.lon,end_point.lon]))
        
    #The transects are on slightly different grids due to cross-section geometry.
    #So, interpolate all the grids and combine, using the first cross-section as a reference
    theta_cross = xr.concat([theta_cross.interp({"cross_line_idx":theta_crosses[0].cross_line_idx},method="nearest") for theta_cross in theta_crosses],"transect")#.mean("transect")
    wspd_cross = xr.concat([wspd_cross.interp({"cross_line_idx":wspd_crosses[0].cross_line_idx},method="nearest") for wspd_cross in wspd_crosses],"transect")#.mean("transect")
    theta_pert_cross = xr.concat([theta_pert_cross.interp({"cross_line_idx":theta_pert_crosses[0].cross_line_idx},method="nearest") for theta_pert_cross in theta_pert_crosses],"transect")#.mean("transect")
    u_cross = xr.concat([u_cross.interp({"cross_line_idx":u_crosses[0].cross_line_idx},method="nearest") for u_cross in u_crosses],"transect")#.mean("transect")
    sru_cross = xr.concat([sru_cross.interp({"cross_line_idx":sru_crosses[0].cross_line_idx},method="nearest") for sru_cross in sru_crosses],"transect")#.mean("transect")
    w_cross = xr.concat([w_cross.interp({"cross_line_idx":w_crosses[0].cross_line_idx},method="nearest") for w_cross in w_crosses],"transect")#.mean("transect")
    dbz_cross = xr.concat([dbz_cross.interp({"cross_line_idx":dbz_crosses[0].cross_line_idx},method="nearest") for dbz_cross in dbz_crosses],"transect")#.mean("transect")
    b_cross = xr.concat([b_cross.interp({"cross_line_idx":b_crosses[0].cross_line_idx},method="nearest") for b_cross in b_crosses],"transect")#.mean("transect")
    p_pert_cross = xr.concat([p_pert_cross.interp({"cross_line_idx":p_pert_crosses[0].cross_line_idx},method="nearest") for p_pert_cross in p_pert_crosses],"transect")#.mean("transect")
    rho_cross = xr.concat([rho_cross.interp({"cross_line_idx":rho_crosses[0].cross_line_idx},method="nearest") for rho_cross in rho_crosses],"transect")#.mean("transect")
    tv_cross = xr.concat([tv_cross.interp({"cross_line_idx":tv_crosses[0].cross_line_idx},method="nearest") for tv_cross in tv_crosses],"transect")#.mean("transect")
    
    #Comine into a dataset and save
    ds = xr.Dataset({"theta":theta_cross,\
                "ground_rel_ws":wspd_cross,\
                "theta_pert":theta_pert_cross,\
                "along_transect_wind":u_cross,\
                "system_relative_along_transect_wind":sru_cross,\
                "w":w_cross,\
                "reflectivity":dbz_cross,\
                "p_pert":p_pert_cross,\
                "rho":rho_cross,\
                "tv":tv_cross,\
                "buoyancy":b_cross},
                   attrs={"transect_starts_lon":[s[0] for s in start_points],
			  "transect_starts_lat":[s[1] for s in start_points],
                          "transect_ends_lon":[e[0] for e in end_points],
                          "transect_ends_lat":[e[1] for e in end_points],
                          "transect_angles":angles,
                          "transect_length_km":d_km,
                          "system_motion_top_m":mean_motion_top,
                          "output_path":path,
                          "fname":f,
                          "line_xlim":xlim,
                          "line_ylim":ylim})

    convert_xyloc(drop_attrs(ds)).to_netcdf("/g/data/w40/ab4502/WRF_simulations/transects/"+\
			    path.split("/")[-2] + "_" + f.split("/")[-1].split("_")[-1] + ".nc")

if __name__ == "__main__":
    
    ################
    #Control run
    ################

    #Load WRF output
    path = "/g/data/w40/ab4502/WRF_simulations/tasman_d04_20231015_20231016_ramssa_shortlead/"
    fs = ["highres_d04_2023-10-16_02:30:00",\
	"highres_d04_2023-10-16_03:00:00",\
	"highres_d04_2023-10-16_03:30:00",\
	"highres_d04_2023-10-16_04:00:00"]

    #Define the region where the squall line is. Needs to be done manually at this stage
    xlims = [[152.0,156],[152.0,156],[153.2,156],[153.4,156]]
    ylims = [[-36.8,-35],[-36.8,-35],[-37,-35],[-37,-35]]

    #for f, xlim, ylim in zip(fs, xlims, ylims):
    #    print(path+f)
    #    plot_transects(path,f,xlim,ylim)
    #    ds = get_transects(path,f,xlim,ylim)

    ################
    #Nov run
    ################

    #Load WRF output
    path = "/g/data/w40/ab4502/WRF_simulations/tasman_d04_20231015_20231016_ramssa_nov23_shortlead/"
    f = ["highres_d04_2023-10-16_02:30:00",\
	"highres_d04_2023-10-16_03:00:00",\
	"highres_d04_2023-10-16_03:30:00",\
	"highres_d04_2023-10-16_04:00:00"]

    #Define the region where the squall line is. Needs to be done manually at this stage
    xlims = [[153.5,154.6],[153.5,155.5],[153.9,155.5],[154,156]]
    ylims = [[-37,-36],[-37,-36],[-37,-36],[-37,-36]]
    
    #for f, xlim, ylim in zip(fs, xlims, ylims):
        #print(path+f)
        #plot_transects(path,f,xlim,ylim)
        #ds = get_transects(path,f,xlim,ylim)

    ################
    #Climo run
    ################

    #Load WRF output
    path = "/g/data/w40/ab4502/WRF_simulations/tasman_d04_20231015_20231016_ramssa_climo_shortlead/"
    f = ["highres_d04_2023-10-16_02:30:00",\
	"highres_d04_2023-10-16_03:00:00",\
	"highres_d04_2023-10-16_03:30:00",\
	"highres_d04_2023-10-16_04:00:00"]

    #Define the region where the squall line is. Needs to be done manually at this stage
    xlims = [[153.5,154.6],[153.5,154.6],[153.5,154.6],[153.5,154.6]]
    ylims = [[-36.8,-36],[-36.8,-36],[-37,-36],[-37,-36]]
    
    #for f, xlim, ylim in zip(fs, xlims, ylims):
    #    print(path+f)
        #plot_transects(path,f,xlim,ylim)
        #ds = get_transects(path,f,xlim,ylim)

    ################
    #+3K run
    ################

    #Load WRF output
    path = "/g/data/w40/ab4502/WRF_simulations/tasman_d04_20231015_20231016_ramssa_shortlead_3Kplus/"
    fs = ["highres_d04_2023-10-16_02:30:00",\
	"highres_d04_2023-10-16_03:00:00",\
	"highres_d04_2023-10-16_03:30:00",\
	"highres_d04_2023-10-16_04:00:00"]

    #Define the region where the squall line is. Needs to be done manually at this stage
    xlims = [[153.0,154.5],[153.5,155],[153.5,156],[154,156]]
    ylims = [[-36.8,-35],[-36.8,-35],[-36.8,-35],[-36.7,-35]]

    #for f, xlim, ylim in zip(fs, xlims, ylims):
        #print(path+f)
        #plot_transects(path,f,xlim,ylim)
        #ds = get_transects(path,f,xlim,ylim)

    ################
    #-3K run
    ################

    #Load WRF output
    path = "/g/data/w40/ab4502/WRF_simulations/tasman_d04_20231015_20231016_ramssa_shortlead_3Kminus/"
    fs = ["highres_d04_2023-10-16_02:30:00",\
	"highres_d04_2023-10-16_03:00:00",\
	"highres_d04_2023-10-16_03:30:00",\
	"highres_d04_2023-10-16_04:00:00"]

    #Define the region where the squall line is. Needs to be done manually at this stage
    xlims = [[152,154],[152.0,154],[153.0,154],[153.4,155]]
    ylims = [[-36.85,-35],[-36.8,-35],[-36.8,-35],[-36.5,-35]]

    for f, xlim, ylim in zip(fs, xlims, ylims):
        print(path+f)
        plot_transects(path,f,xlim,ylim,grad_theta_thresh=1)
        ds = get_transects(path,f,xlim,ylim,grad_theta_thresh=1)
