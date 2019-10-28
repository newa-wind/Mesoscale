"""
PYTHON3 Basic Function to Create NEWA Output NETCDF-Files

V0 - MD - 28.03.2018 - Moved functions to this module
V1 - MD - 16.04.2018 - Final Check of the Files
V2 - MD - 26.04.2018 - Bug Fixes According to Slack Discussion

(c) Fraunhofer IWES
"""
import datetime
import numpy as np
import netCDF4 as nc
from pyproj import Proj

INSTITUTION="The NEWA Consortium"
EXPERIMENT="Production Run"

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

def createnc(ncfout,xlat,xlon,times=None,zvals=None,wsvals=None,\
             wdvals=None,olvals=None,attbts=None,ftype="timeseries",dims=[7,180,180]):
    """ Creates a NETCDF Output file and assignes dimensions to it
    ---------------------------------------------------
    In: File Name (ncfout)
        Time Vector (times)
        Height Vector (zvals)
        xlat, xlon (3D Fields)
        Optional: dims - Dimensions in z,y,x directions
    --------------------------------------------------
    Out: Status
    --------------------------------------------------
    Writes: New Output NETCDF File With Dimenions
    """
    nc_out=nc.Dataset(ncfout,'w',clobber=True)

    # Set Attributes to the File
    if attbts is not None:
        final_attbts={}
        # Define projection
        proj_lcc  = pj_lcc = Proj("+proj=lcc +lat_1={TRUELAT1} +lat_2={TRUELAT2} +lat_0={MOAD_CEN_LAT} +lon_0={STAND_LON} +x_0=0 +y_0=0 +a=6370000 +b=6370000".format(**attbts))

        # Get x&y of domain center
        xcen, ycen = pj_lcc(attbts['CEN_LON'], attbts['CEN_LAT'])

        for key in attbts:
            if str(key).find("STAG") <= 0 : # Remove Staggered Grid Information
                final_attbts.update({key:attbts[key]})
        nc_out.setncatts(final_attbts)
        # Create a CRS Variable for the Projection (GIS Readability)
        crsv=nc_out.createVariable('crs','c')
        crsv.semi_major_axis = 6370000.0
        crsv.inverse_flattening = 0.0
        crsv.grid_mapping_name = "lambert_conformal_conic"
        crsv.longitude_of_central_meridian = attbts["STAND_LON"]
        crsv.false_easting = 0.0
        crsv.false_northing = 0.0
        crsv.latitude_of_projection_origin = attbts["MOAD_CEN_LAT"]
        crsv.standard_parallel = [attbts["TRUELAT1"],attbts["TRUELAT2"]]
        crsv.longitude_of_prime_meridian = 0.0
        crsv.proj = proj_lcc.srs



    # Override Institution and Experiment
    nc_out.INSTITUTION=INSTITUTION
    nc_out.EXPERIMENT=EXPERIMENT
    nc_out.Conventions="CF-1.6"

    # Create Dimensions First
    if ftype=="timeseries":
        nc_out.TITLE='Timeseries of the New European Wind Atlas from WRF V3.8.1'
        nc_out.createDimension('time',None)
        nc_out.createDimension('DateStrLen',19)
        nc_out.createDimension('height',dims[0])
        nc_out.createDimension('south_north',dims[1])
        nc_out.createDimension('west_east',dims[2])
        # Create Time Vector as Integer
        timesn = nc_out.createVariable('time','i8',('time',))
        timesn.units = "minutes since 1900-01-01 00:00:00.0"
        timesn.calendar = "gregorian"
        timesn.long_name = "Time"
        timesn.standard_name = "time"
        timesn[:] = nc.date2num(createdatv(times),units=timesn.units,calendar=timesn.calendar)
        # Create additional Time Vector as Character
        timesc = nc_out.createVariable('Times', 'c', ('time','DateStrLen'))
        timesc.format = "YYYY-MM-DD_HH:MM:SS"
        timesc.long_name = "Time"
        timesc[:] = times[:]
        # Height
        hgts = nc_out.createVariable('height','f4',('height',))
        hgts.units="m"
        hgts.long_name="Height above Ground"
        hgts.standard_name="height"
        hgts[:] = zvals
        # y
        south_north = nc_out.createVariable('south_north','f4',('south_north',))
        south_north.long_name = "y-coordinate in Cartesian system"
        south_north.units = "m"

        dy = attbts["DY"]
        ny = attbts["SOUTH-NORTH_PATCH_END_UNSTAG"]
        ymin = ycen - dy * (ny - 1) / 2
        s_n = np.linspace(0, ny-1, ny) * dy + ymin
        south_north[:] = s_n

        # x
        west_east = nc_out.createVariable('west_east','f4',('west_east',))
        west_east.long_name = "x-coordinate in Cartesian system"
        west_east.units = "m"

        dx = attbts["DX"]
        nx = attbts["WEST-EAST_PATCH_END_UNSTAG"]
        xmin = xcen - dx * (nx - 1) / 2
        e_w = np.linspace(0, nx-1, nx) * dx + xmin
        west_east[:] = e_w

    elif ftype=="roughness":
        nc_out.title='NEWA Roughness'
        nc_out.createDimension('south_north',dims[0])
        nc_out.createDimension('west_east',dims[1])

    elif ftype=="tabfile":
        nc_out.title='NEWA WasP Tab File'
        nc_out.createDimension('south_north',dims[0])
        nc_out.createDimension('west_east',dims[1])
        nc_out.createDimension('sector',dims[2])
        nc_out.createDimension('wind',dims[3])
        nc_out.createDimension('stab',dims[4])

        # Wind Speed Class
        wscl = nc_out.createVariable('wspdCl','f4',('wind',))
        wscl.units="ms-1"
        wscl.long_name="Velocity of bin centre"
        wscl[:] = wsvals

        # Wind Speed Class
        wdcl = nc_out.createVariable('wdirCl','f4',('sector',))
        wdcl.units="ms-1"
        wdcl.long_name="Velocity of bin centre"
        wdcl[:] = wdvals

        # Stability
        lcl = nc_out.createVariable('Ltypical','f4',('stab',))
        lcl.units="m"
        lcl.long_name="L typical"
        lcl[:] = olvals

    # Lat and Lon
    lats = nc_out.createVariable("XLAT", 'f4', ('south_north','west_east'), zlib=True,complevel=9)
    lats[:] = xlat[:]
    lats.units="degree_north"
    lats.long_name="Center Latitude of Grid Cell"
    lats.standard_name="latitude"
    lons = nc_out.createVariable("XLON", 'f4', ('south_north','west_east'), zlib=True,complevel=9)
    lons[:] = xlon[:]
    lons.units="degree_east"
    lons.long_name="Center Longitude of Grid Cell"
    lons.standard_name="longitude"
    nc_out.close()
    return(None)

def addtonc(ncfout,key,vd,ofield,ftype="timeseries"):
    """ Adds Variables to an existing NETCDF File
    ---------------------------------------------------
    In: File Name (ncfout)
        Variable Name (key)
        Variable Dictioary (vd - contains: 'units, 'name', 'dims')
        Field to write to file (ofield)
    --------------------------------------------------
    Out: Status
    --------------------------------------------------
    Writes: Data to existing NETCDF File specified by ncfout
    """
    nc_out=nc.Dataset(ncfout,'r+')
    if ftype=="timeseries":
        diml=['time','height','south_north','west_east'] # Tuple of Dimensions
        if vd['dims']==4:
            dimtup=tuple(diml)
        elif vd['dims']==3:
            dimtup = tuple([c for c in diml if c != "height"])
        elif vd['dims']==2:
            dimtup = tuple([c for c in diml if c not in ["height","time"]])
    elif ftype=="roughness":
        diml=['south_north','west_east']
        dimtup=tuple(diml)
    elif ftype=="tabfile":
        diml=['south_north','west_east','sector','wind','stab']
        if vd['dims']==3:
            dimtup=tuple(diml.remove('wind').remove('stab'))
        if vd['dims']==2:
            dimtup=tuple(diml.remove('wind').remove('stab').remove('sector'))
    if key in ("TKE", "ABLAT_CYL", "ACCRE_CYL"):
        outv=nc_out.createVariable(key, 'f4', dimtup, zlib=True,
                                   complevel=9, fill_value=-999.)
    else:
        outv=nc_out.createVariable(key,'f4',dimtup,zlib=True,complevel=9)
    outv.units=vd['units']
    outv.long_name=vd['name']
    if vd['std_name'] is not None:
        outv.standard_name=vd['std_name']
    if key=="PRECIP":
        outv.cell_methods="time: sum"
    outv.grid_mapping="crs"
    outv.coordinates="XLAT XLON"
    outv[:]=ofield[:]
    nc_out.close()
    return(None)
