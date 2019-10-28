"""
This script reduces (levels and variables) and
outputs reduced data in the final format
of the NEWA atlas (as of 04-2018).

authors: Martin Doerenkaemper (MD)

V0: 28.03.2017 - MD - basic scripting
V1: 25.07.2017 - MD - small bugfixes
V2: 30.01.2018 - BW - modification for final NEWA Light output
V3: 28.03.2018 - MD - Modification for NEWA production runs
V4: 26.04.2018 - MD - Bug Fixes according to Slack Discussions
V5: 28.05.2018 - MD - WRF-Python Interpolation from neda,DTU

USAGE: python3 newaoutfinal.py FILENAME
              (Optional: PREFIXOFNEWFILE)

(c) Fraunhofer IWES
(c) Carl von Ossietzky Universität Oldenburg
"""
import re
import os
import sys
import glob
import copy
import datetime
import netCDF4 as nc
import newa2nc
import numpy as np
import wrf

# Definition of Constants
RLAIR=287.0    # WRF Dry air gas constant
Rw = 461.4     # WRF Water gas constant
CP=1004.5      # WRF Specific heat capacity air
HUMCONST=0.622 # Ideal Gas law for Water Vapour and Dry Air
G = 9.81       # Gravitational Constant

def interp_4d(ipvar, wrf_hgts, out_hgts):
    """ Interpolates vertically to heights specified by
    out_hgts using the interplevel function from wrf-python
    -----------------------------------------------------
    In: 4D field of variable to interpolate (ipvar)
        Heights of the WRF simulation (wrf_hgts)
        Heights to interpolate to (out_hgts)
    -----------------------------------------------------
    Out: Interpolated 4D field (out)
    --------------------------------------
    """
    shape = list(ipvar.shape)
    shape[1] = len(out_hgts)
    out = np.empty(shape, ipvar.dtype)
    for k, out_hgt in enumerate(out_hgts):
        out[:,k,:,:] = wrf.interplevel(ipvar, wrf_hgts, out_hgt)
    return(out)

def createdatv(times):
     """ Creates a continuous time vector in pythons datetime format
     ---------------------------------------------------
     In: RAW Time Vector from wrfout*.nc file (Times)
     --------------------------------------------------
     Out: Time Vector in Datetime format (datev)
     """
     t0=datetime.datetime.strptime(times[0].tostring().decode('utf-8'), '%Y-%m-%d_%H:%M:%S')
     t1=datetime.datetime.strptime(times[1].tostring().decode('utf-8'), '%Y-%m-%d_%H:%M:%S')
     ts=(t1-t0).total_seconds()
     datev=[]
     for i in range(0,len(times)):
        datev.append(t0+datetime.timedelta(seconds=i*ts))
     return(datev)

def interphgt(ipvar,iphgt,hgts):
    """ Interpolates any 3D field [z,y,x] from WRF linearly to desired height
    ---------------------------------------------------
    In: Desired Height (iphgt [m])
        3D Field [z,y,x] (ipvar)
        Geometric Corresponding Height of the 3D field [z,y,x] (hgts)
    --------------------------------------------------
    Out: Height interpolated 2D field at iphgt
    """
    dims=hgts.shape
    ipmat=np.empty([dims[1],dims[2]])
    for j in range(0,dims[1]):
        for i in range(0,dims[2]):
            ipmat[j,i]=np.interp(iphgt,hgts[:,j,i],ipvar[:,j,i])
    return(ipmat)

def twrftot(twrf,p,pb):
    """ Converts WRF Temperature to normal Temperature
    ---------------------------------------------------
    In: WRF Temperature (twrf [K])
        Pressure (P [Pa])
        Base Pressure (PB [Pa])
    --------------------------------------------------
    Out: Regular Temperature [K]
    """

    PREF=100000.0       # Reference Pressure in Pa (1000.0 hPa)
    KAPPA=RLAIR/CP      # Meteorological Constant KAPPA
    TOFF=300.0          # Constant WRF Temperature Offset
    tink=(twrf+TOFF)*(((p+pb)/PREF)**KAPPA)
    return(tink)

def uv2wswd(u,v):
    """ Converts wind speed in vector format to horizontal
    wind speed and direction
    ----------------------------------------------------
    In:  u - eastward component of the wind speed
         v - northward component of the wind speed
    ----------------------------------------------------
    Out: ws - horizontal absolute wind speed
         wd - horizontal wind direction [degrees]
    """
    ws=np.sqrt(u*u+v*v)
    wdrad = np.rad2deg(np.arctan2(v/ws,u/ws))
    wd = np.mod((270.0 - wdrad),360.0)
    return([ws,wd])

def erelwind(uw,vw,cosalph,sinalph):
    """ Converts WRF variables U and V to earth relative winds
        and calculates horizontal wind speed and direction
    ---------------------------------------------------
    In: WRF relative wind speed components (U [m/s],V[m/s])
        sin and cos of WRF grid (COSALP, SINALPH)
    --------------------------------------------------
    Out: Horizontal Wind (WS[m/s], WD[degree])
    """
    ue = uw*cosalph - vw*sinalph    # Rotate Coordinates to Earth grid
    ve = vw*cosalph + uw*sinalph    # uw,vw are coordinates in/on WRF Grid
    [ws,wd]=uv2wswd(ue,ve)
    return([ws,wd])

def airdens(p,tk,q=None,dryfl=False):
    """ Calculates the density WRF variables U and V to earth relative winds
        and calculates horizontal wind speed and direction
    ---------------------------------------------------
    In: p - Pressure [Pa]
        t - Temperature [K]
        q - Specific Humidity [kg/kg]
        dryfl - Flag to derive dry or humid[True/False]
    --------------------------------------------------
    Out: Density of  Horizontal Wind (WS[m/s], WD[degree])
    """
    if dryfl or (q is None):
        dens=(p/(RLAIR*tk))
    else:
        vp = (q * p) / (HUMCONST + q)
        densdry=(p - vp) / (RLAIR * tk)
        densvp=vp/(Rw * tk)
        dens=densdry + densvp
    return(dens)

def interpnewa(ncin,ncfout,vdict,iphgt): # Interpolate 3D Fields
     """ Main Function to Interpolate WRF Variables
     ---------------------------------------------------
     In: In File Name (ncin - wrfout Netcdf-File)
         Variable Dictioary (vdict - contains: 'units, 'name', 'dims')
         Desired Interpolation Heights (iphgt)
     --------------------------------------------------
     Out: Status
     --------------------------------------------------
     Writes: Data to existing NETCDF File specified by ncfout
     """
     wsswitch=True
     nz1=0
     nz2=20 # Take first 25 Levels only, needs to be changed if max(iphgt) >> 500m
     # Infile
     nc_in=nc.Dataset(ncin,'r')
     # Get attributes
     attributes = {a:nc_in.getncattr(a) for a in nc_in.ncattrs()}
     attributes.update({
                "CONST_R_DRY_AIR": RLAIR,
                "CONST_R_WATER": Rw,
                "CONST_RATIO_Rd_Rw": HUMCONST,
                "CONST_Cp_DRY_AIR": CP,
                "CONST_GRAVITY": G
     })
     # Outvars
     oldtime = nc_in.variables['Times'][:]
     lats = nc_in.variables['XLAT'][0,:,:]
     lons = nc_in.variables['XLONG'][0,:,:]
     dims=[len(iphgt),len(lats[:,0]),len(lats[0,:])]
     cstat=newa2nc.createnc(ncfout,lats,lons,oldtime,iphgt,attbts=attributes,ftype="timeseries",dims=dims) # Create Basic NetCFD4 File
     # Calculate Height Variables:
     [yd,xd]=lats.shape
     td=len(oldtime)

     # Prepare Height variable for Interpolation
     topo = nc_in.variables['HGT']  # Height of the ground
     znt = nc_in.variables['ZNT'][:]
     phb = nc_in.variables['PHB'][:,nz1:(nz2+1),:,:] # extract/copy the data
     ph = nc_in.variables['PH'][:,nz1:(nz2+1),:,:]   # extract/copy the data
     geopt_stag = phb + ph
     geopt = 0.5 * (geopt_stag[:,:-1, ...] + geopt_stag[:,1:, ...])
     z = geopt / G
     topo_4d = np.expand_dims(topo, 1)
     hgts = z - topo_4d

     # Delete Fields to save Memory
     del phb,ph,geopt_stag, geopt,z,topo_4d

     # Calculate the Alpha Cone
     cosa = nc_in.variables['COSALPHA'][:,:,:]
     sina = nc_in.variables['SINALPHA'][:,:,:]

     # Copy or Interpolate Fields to desired heights and write to NetCDF4 File
     for key in sorted(vdict):
        # First, take care of the 3D Fields that can be copied
        if vdict[key]["wrfname"] is not None and vdict[key]["dims"]==3:
           print("     Copying 3D Field: "+key)
           try:
               ifield = nc_in.variables[vdict[key]["wrfname"]][:,:,:]
           except KeyError:
               print("     +Error:"+vdict[key]["wrfname"]+" not in File "+ncin+"!")

        # Then: Interpolate the 4D Fields that do not need to be derived:
        elif vdict[key]["wrfname"] is not None and vdict[key]["dims"]==4:
           print("     Interpolating 4D Field: "+key)
           try:
               ifield = np.empty([td,len(iphgt),yd,xd])
               twrf = nc_in.variables[vdict[key]["wrfname"]][:,nz1:nz2,:,:]
               ifield = interp_4d(twrf, hgts, iphgt)
           except KeyError:
               print("     Error:"+vdict[key]["wrfname"]+" not in File "+ncin+" setting to -999.0!")
               ifield = -999. + np.zeros((td, len(iphgt), yd, xd))
        # Add 2D fields that don't change over time
        elif vdict[key]["wrfname"] is not None and vdict[key]["dims"]==2:
            print("     Subsetting 2D Field: " + key)
            try:
                ifield = nc_in.variables[vdict[key]["wrfname"]][1,:,:]
            except KeyError:
                print("     +Error:"+vdict[key]["wrfname"]+" not in File "+ncin+"!")

        # Last Step: Special Treetment of Variables that need to be derived
        elif key == "ALPHA":
           # Alpha Cone
           print("     Calculating Custom Field: "+key)
           ifield=np.degrees(np.arctan2(sina[0,:,:],cosa[0,:,:]))
        elif key == "HGT":
           print("     Calculating 2D Custom Field: "+key)
           ifield=topo[0,:,:]
        elif key == "T":
           print("     Calculating 4D Custom Field: T from WRF Potential temperature")
           ifield = np.empty([td,len(iphgt),yd,xd])
           potwrf = nc_in.variables["T"][:,nz1:nz2,:,:]  # Pot. Temperature - 300.0
           pwrf = nc_in.variables["P"][:,nz1:nz2,:,:]    # Pressure WRF
           pbwrf = nc_in.variables["PB"][:,nz1:nz2,:,:]  # Base Pressure WRF
           twrf = twrftot(potwrf,pwrf,pbwrf)
           ifield = interp_4d(twrf, hgts, iphgt)
        elif key == "PRECIP":
            print("     Calculating Custom Field: "+key)
            try:
                # Precipication (Sum of Convective and Non-Convective)
                prec1=nc_in.variables["PREC_ACC_NC"][:,:,:]
                prec2=nc_in.variables["PREC_ACC_C"][:,:,:]
                ifield=prec1[:,:,:]+prec2[:,:,:]
            except KeyError:
                print("     Error:"+key+" not in File "+ncin+"!")
        elif key == "RHO":
           print("     Calculating Custom Field: "+key)
           # Density of Dry Air
           t2=nc_in.variables["T2"][:,:,:]
           q2=nc_in.variables["Q2"][:,:,:]
           psfc=nc_in.variables["PSFC"][:,:,:]
           ifield=airdens(psfc,t2,q2,dryfl=False)
        elif key == "WS10":
           print("     Calculating Custom Field: "+key)
	   # Horizontal Wind Speed at 10m
           ufield = nc_in.variables["U10"][:,:,:]
           vfield = nc_in.variables["V10"][:,:,:]
           [ifield, wd] = erelwind(ufield,vfield,cosa,sina)
        elif key == "WD10":
           print("     Calculating Custom Field: "+key)
           # Horizontal Wind Direction at 10m
           ufield = nc_in.variables["U10"][:,:,:]
           vfield = nc_in.variables["V10"][:,:,:]
           [ws, ifield] = erelwind(ufield,vfield,cosa,sina)
        elif key == "TKE":
            if nc_in.BL_PBL_PHYSICS in (2, 4): # MYJ, QNSE
                print("     Calculating 4D Custom Field: TKE from TKE_PBL")
                tke = nc_in.variables['TKE_PBL'][:,nz1:nz2,...]
                ifield = interp_4d(tke, hgts, iphgt)
                del tke
            elif nc_in.BL_PBL_PHYSICS in (5, 6): # MYNN2/3
                print("     Calculating 4D Custom Field: TKE as 0.5 * QKE")
                qke = nc_in.variables['QKE'][:,nz1:nz2,...]
                qke = interp_4d(qke, hgts, iphgt)
                ifield = qke * 0.5
                del qke
            elif nc_in.BL_PBL_PHYSICS in (1, 7): # YSU, ACM2
                print("     There is no TKE for this scheme, filling with -999.")
                ifield = -999. + np.zeros((td, len(iphgt), yd, xd))
            else:
                raise ValueError("     Don't know how to process TKE for PBL Scheme %s." % nc_in.BL_PBL_PHYSICS)

        elif (key == "WS" or key == "WD" or key == "PD") and wsswitch:
           print("     Calculating Custom Fields: WS, WD, & PD")
           wsswitch=False

           # Read wind fields
           uwrf = nc_in.variables['U'][:,nz1:nz2,:,:]     # U-Component of the Flow
           vwrf = nc_in.variables['V'][:,nz1:nz2,:,:]     # V-Component of the Flow

           # Unstager
           uen = 0.5 * (uwrf[...,0:-1] + uwrf[...,1:])
           ven = 0.5 * (vwrf[:,:,0:-1,:] + vwrf[:,:,1:,:])

           # Wind speed
           wsf = uv2wswd(uen,ven)[0]
           wsfield = interp_4d(wsf, np.log(hgts), np.log(iphgt))

           # Wind Direction
           # First convert to earth relative
           uenorig=copy.copy(uen)
           for k in range(0,(nz2-nz1)):
               uen[:,k,:,:] = uen[:,k,:,:]*cosa[:,:,:] - ven[:,k,:,:]*sina[:,:,:]         # Rotate Coordinates to Earth grid
               ven[:,k,:,:] = ven[:,k,:,:]*cosa[:,:,:] + uenorig[:,k,:,:]*sina[:,:,:]

           ueip = interp_4d(uen, hgts, iphgt)
           veip = interp_4d(ven, hgts, iphgt)
           wdfield = uv2wswd(ueip,veip)[1]

           # RHO
           t2 = np.expand_dims(nc_in.variables["T2"][:], 1)
           q2 = np.expand_dims(nc_in.variables["Q2"][:], 1)
           psfc = np.expand_dims(nc_in.variables["PSFC"][:], 1)
           pdfield = 0.5 * airdens(psfc,t2,q2,dryfl=False) * wsfield**3.0

           astat = newa2nc.addtonc(ncfout,"WS",vdict["WS"],wsfield,ftype="timeseries")
           astat = newa2nc.addtonc(ncfout,"WD",vdict["WD"],wdfield,ftype="timeseries")
           astat = newa2nc.addtonc(ncfout,"PD",vdict["PD"],pdfield,ftype="timeseries")
           # Cleanup
           del wsfield,wdfield,pdfield
        # Write Field to File
        if key not in ["WS","WD","PD"]:
           astat = newa2nc.addtonc(ncfout, key, vdict[key], ifield, ftype="timeseries")
           del ifield
        if astat==False:
            print("Something went wrong in writing variables to the files!")
     # Close Infile
     nc_in.close()
     return(True)

def main():

    # Specify Heights to Interpolate to
    fhgt=[50.0,75.0,100.0,150.0,200.0,250.0,500.0]

    # Specify Parameters to
    vardict={
     "ABLAT_CYL":{"name":"Ice Ablation on Standard Cylinder","std_name":None,"units":"kg","dims":4,"wrfname":"ABLAT_CYL"},
     "ACCRE_CYL":{"name":"Ice Accretion on Standard Cylinder","std_name":None,"units":"kg","dims":4,"wrfname":"ACCRE_CYL"},
     "ALPHA":{"name":"Alpha Cone","std_name":None,"units":"1","dims":2,"wrfname":None},
     "HFX":{"name":"Sensible Heat Flux","std_name":"surface_upward_sensible_heat_flux","units":"W m-2","dims":3,"wrfname":"HFX"},
     "HGT":{"name":"Elevation","std_name":"ground_level_altitude","units":"m","dims":2,"wrfname":"HGT"},
     "LANDMASK":{"name":"Landmask (1 FOR LAND, 0 FOR WATER)","std_name":"land_binary_mask","units":"1","dims":2, "wrfname": "LANDMASK"},
     "LH":{"name":"Latent Heat Flux","std_name":"surface_upward_latent_heat_flux","units":"W m-2","dims":3,"wrfname":"LH"},
     "LU_INDEX":{"name":"Dominant Land Use Category (USGS)","std_name":None,"units":"","dims":2, "wrfname": "LU_INDEX"},
     "PBLH":{"name":"PBL Height","std_name":"atmosphere_boundary_layer_thickness","units":"m","dims":3,"wrfname":"PBLH"},
     "PD":{"name":"Wind Power Density","std_name":None,"units":"W m-2","dims":4,"wrfname":None},
     "PSFC":{"name":"Surface Pressure","std_name":"air_pressure","units":"Pa","dims":3,"wrfname":"PSFC"},
     "PRECIP":{"name":"Total Precipitation for last 30 minutes","std_name":"precipitation_amount","units":"kg m-2","dims":3,"wrfname":None},
     "Q2":{"name":"Specific Humidity at 2m","std_name":"specific_humidity","units":"1","dims":3,"wrfname":"Q2"},
     "QVAPOR":{"name":"Humidity Mixing Ratio","std_name":"humidity_mixing_ratio","units":"1","dims":4,"wrfname":"QVAPOR"},
     "RHO":{"name":"Air Density","std_name":"air_density","units":"kg m-3","dims":3,"wrfname":None},
     "RMOL":{"name":"Inverse Obukhov Length","std_name":None,"units":"m-1","dims":3,"wrfname":"RMOL"},
     "SEAICE":{"name":"Sea Ice Fraction","std_name":"sea_ice_area_fraction", "units":"1", "dims": 3, "wrfname": "SEAICE"},
     #Identical_to_TSK_over_water "SST":{"name":"Sea Surface Temperature","std_name":"sea_surface_temperature","units":"K","dims":3,"wrfname":"SST"},
     "SWDDIR":{"name":"Shortwave Diffuse Incident Radiation","std_name":None,"units":"W m-2","dims":3,"wrfname":"SWDDIR"},
     "SWDDNI":{"name":"Shortwave Direct Normal Radiation","std_name":None,"units":"W m-2","dims":3,"wrfname":"SWDDNI"},
     "T":{"name":"Air Temperature","std_name":"air_temperature","units":"K","dims":4,"wrfname":None},
     "T2":{"name":"Air Temperature at 2m","std_name":"air_temperature","units":"K","dims":3,"wrfname":"T2"},
     "TKE":{"name":"Turbulent Kinetic Energy","std_name":None,"units":"m2 s-2","dims":4,"wrfname":None},
     "TSK":{"name":"Skin Temperature","std_name":"surface_temperature","units":"K","dims":3,"wrfname":"TSK"},
     "UST":{"name":"Friction Velocity","std_name":None,"units":"m s-1","dims":3,"wrfname":"UST"},
     "WD":{"name":"Wind Direction","std_name":"wind_from_direction","units":"degree","dims":4,"wrfname":None},
     "WD10":{"name":"Wind Direction at 10m","std_name":"wind_from_direction","units":"degree","dims":3,"wrfname":None},
     "WS":{"name":"Wind Speed","units":"m s-1","std_name":"wind_speed","dims":4,"wrfname":None},
     "WS10":{"name":"Wind Speed at 10m","std_name":"wind_speed","units":"m s-1","dims":3,"wrfname":None},
     "ZNT":{"name":"Aerodynamic Roughness Length","std_name":"surface_roughness_length","units":"m","dims":3,"wrfname":"ZNT"}
     }


    # Measure Execution Time
    dt1=datetime.datetime.now()

    # Start Loop over simulation period
    fna=sys.argv[1]
    fls=sorted(glob.glob(fna+"/wrfout_d03*"))

    prefix=os.getenv('NEWARUNNAME')

    if prefix is None:
    	prefix="NEWA-"

    spinup=24


    print("-------------------------------------------------------------------------")
    print("--- Starting TimeSeries interpolation in: "+fna+"  ---")
    print("-------------------------------------------------------------------------")

    for i,fl in enumerate(fls[:-1]): # Discard last Timestep
        print(" +++ Interpolating TimeSeries for: "+fl)

        # Generate Filename from old Filename + Spinup-Marker
        dat=re.split('_|\.',fl)
        for elem in dat:
            if re.match('\d{4}-\d{2}-\d{2}',elem):
                dstr=elem
                stime=datetime.datetime.strptime(dstr,"%Y-%m-%d")
                if i==0:
                    s0t=stime
                if stime<(s0t+datetime.timedelta(hours=spinup)):
                    fend="-SP.nc"
                else:
                    fend=".nc"

        flout=prefix+"-"+dstr+fend

        #  RUN the Post-Processing
        istat=interpnewa(fl,flout,vardict,fhgt)
        dt2=datetime.datetime.now()
        print(" +++ Interpolating TimeSeries for: "+fl)
        print(" +++ Execution Time: "+str((dt2-dt1).seconds)+" seconds.")
        print("-------------------------------------------------------------------------")

    # Measure Execution Time of the Script
    dt3=datetime.datetime.now()
    print("-------------------------------------------------------------------------")
    print("Total Execution Time of Time Series Interpolation: "+str((dt3-dt1).seconds)+" seconds.")
    print("---------------------------------Done.-----------------------------------")


if __name__ == '__main__':
  main()
