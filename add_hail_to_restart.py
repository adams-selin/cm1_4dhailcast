import numpy as np
import netCDF4 as nc4
import pandas as pd

#open the restart file with hail data in it
nc = nc4.Dataset('/glade/scratch/radams/hail/20120529_quickrestart/cm1out_rst_000001.nc','r')
nhailloc = nc.dimensions['nhailloc'].size
nhailtrajs = nc.dimensions['nhailtrajs'].size
numhailtrajs = np.array(nc.variables['numhailtrajs'][0])  #int
hailtrajtim = np.array(nc.variables['hailtrajtim'][0])  #double or float64
hailloc = np.array(nc.variables['hailloc'])  #float32, dimensions 
nc.close()

#hailtrajtim need to be evenly divisible by hailfrq, set in the namelist. We want that to 
#be 1 second. So pull the time from the 'clean' restart file and round up.
nohail_nc = nc4.Dataset('/glade/scratch/radams/hail/20120529/cm1out_rst_000001.nc','r')
restart_time = np.ceil(np.array(nohail_nc.variables['time']))
nc.close()


#pull in haildata from the actual haildata file to see what it should look like at time 0
nc = nc4.Dataset('/glade/scratch/radams/hail/20120529_quickrestart/cm1out_haildata.nc','r')
hailx_from_t0 = np.array(nc.variables['x'][0,:])
haily_from_t0 = np.array(nc.variables['y'][0,:])
hailz_from_t0 = np.array(nc.variables['z'][0,:])
hailu_from_t0 = np.array(nc.variables['u'][0,:])
hailv_from_t0 = np.array(nc.variables['v'][0,:])
hailw_from_t0 = np.array(nc.variables['w'][0,:])
haild_from_t0 = np.array(nc.variables['d'][0,:])
haildense_from_t0 = np.array(nc.variables['dense'][0,:])
hailvt_from_t0 = np.array(nc.variables['tv'][0,:])
hailts_from_t0 = np.array(nc.variables['ts'][0,:])
hailfw_from_t0 = np.array(nc.variables['fw'][0,:])
hailitype_from_t0 = np.array(nc.variables['itype'][0,:])
nc.close()


#re-assign hailloc values to those from time 0
hailloc[0,0,:] = hailx_from_t0
hailloc[0,1,:] = haily_from_t0
hailloc[0,2,:] = hailz_from_t0
hailloc[0,3,:] = hailu_from_t0
hailloc[0,4,:] = hailv_from_t0
hailloc[0,5,:] = hailw_from_t0
hailloc[0,6,:] = haild_from_t0
hailloc[0,7,:] = haildense_from_t0
hailloc[0,8,:] = hailvt_from_t0
hailloc[0,9,:] = hailts_from_t0
hailloc[0,10,:] = hailfw_from_t0
hailloc[0,11,:] = hailitype_from_t0


#append the t0 hail data to a restart file without hail data in it.
nohail_nc = nc4.Dataset('/glade/scratch/radams/hail/20120529/cm1out_rst_000001.nc','a',format='NETCDF4_CLASSIC')
nhailloc_dim = nohail_nc.createDimension('nhailloc',nhailloc)
nhailtrajs_dim = nohail_nc.createDimension('nhailtrajs',nhailtrajs)

#numhailtrajs_var = nohail_nc.createVariable('numhailtrajs','i4',('time',))
#numhailtrajs_var[:] = numhailtrajs
nohail_nc['numhailtrajs'][:] = numhailtrajs

#hailtrajtim_var = nohail_nc.createVariable('hailtrajtim',np.float64,('time',))
#hailtrajtim_var[:] = restart_time
nohail_nc['hailtrajtim'][:] = restart_time

hailloc_var = nohail_nc.createVariable('hailloc',np.float32,('time','nhailloc','nhailtrajs'))
hailloc_var[:,:,:] = hailloc

nohail_nc.close
















