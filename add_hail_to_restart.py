import numpy as np
import netCDF4 as nc4
import pandas as pd

#As described in the readme doc, this script is designed to add hail embryos to a CM1 restart file
# in order to start modeling their hail trajectories.

#This script assumes you have done the following things:
# 1) Run a "clean" version of your storm with no hail trajectories (ihailtraj = 0).
#    Start your simulation at time 0 until the storm has reached a reasonable 
#    level of intensity at which hail could start.
#    Write out a restart file at that point. I call this the "clean_run".
# 2) Run a "quick restart" version of your storm. Start your simulation at time 0;
#    stop it one model timestep later and write out a restart file. Here you *include*
#    ihailtraj=1 and set your number of hail trajectories (nhailtraj = XX in the namelist,
#    also see init3d.F) to what you would like. The restart file will have starting 
#    conditions for all your hail embryos nicely written out in the correct format.

# This script basically copies the initial hail embryo from the quick_restart run to 
#  your restart file from the "clean" run. Now you can restart your storm, and hail
#  embryos will be initialized in your mature storm!


#open the restart file with hail data in it
nc = nc4.Dataset('/glade/derecho/scratch/radams/20120529/quick_restart/cm1out_rst_000001.nc','r')
nhailloc = nc.dimensions['nhailloc'].size
nhailtrajs = nc.dimensions['nhailtrajs'].size
numhailtrajs = np.array(nc.variables['numhailtrajs'][0])  #int
hailtrajtim = np.array(nc.variables['hailtrajtim'][0])  #double or float64
hailloc = np.array(nc.variables['hailloc'])  #float32, dimensions 
nc.close()

#hailtrajtim need to be evenly divisible by hailfrq, set in the namelist. We want that to 
#be 1 second. So pull the time from the 'clean' restart file and round up.
nohail_nc = nc4.Dataset('/glade/derecho/scratch/radams/20120529/clean_run/cm1out_rst_000004.nc','r')
restart_time = np.ceil(np.array(nohail_nc.variables['time']))
nc.close()


#pull in haildata from the quick_restart haildata file to see what it should look like at time 0
nc = nc4.Dataset('/glade/derecho/scratch/radams/20120529/quick_restart/cm1out_haildata.nc','r')
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
nohail_nc = nc4.Dataset('/glade/derecho/scratch/radams/20120529/rst000025/cm1out_rst_000004.nc','a')
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
















