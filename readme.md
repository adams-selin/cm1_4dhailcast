# 4D HAILCAST hail trajectory model, incorporated into CM1

This code is a four-dimensional, fully integrated version of the Adams-Selin and Ziegler 2016 hail trajectory model in version 19.8 of CM1. The physics of the model is  more fully described in Pounds et al. 2024.

**Please cite the following references**: *Adams-Selin, R.D., 2023: A three-dimensional trajectory clustering technique. Mon. Wea. Rev., 151, 2361–2375, https://doi.org/10.1175/MWR-D-22-0345.1.* and *Pounds, L. E., C. L. Ziegler, R. D. Adams-Selin, and M. I. Biggerstaff, 2024: Analysis of hail production via simulated hailstone trajectories in the 29 May 2012 Kingfisher, Oklahoma, supercell. Mon. Wea. Rev., 152, 245–276, https://doi.org/10.1175/MWR-D-23-0073.1.*

Additionally, this code has been uploaded to Zenodo and given the following DOI: [![DOI](https://zenodo.org/badge/577845382.svg)](https://zenodo.org/doi/10.5281/zenodo.13243355)


 
## Implementation details

A few differences with the ASZ16 hail trajectory model: this subroutine allows for melting at every timestep.  If inside cloud, the melting is determined via the RH87 heat balance equation. If outside or below cloud, melting is determined using Eq. 3 of Goyer et al. (1969), assuming a **spherical** hailstone melting in dry air. This replaces the mean below- cloud melting calculations performed in ASZ16.
Code for oblate spheroidal hailstones is available in the _zshape_ branch.

The CM1 parcel code was repurposed here for advecting the hailstones. Hailstones are first advected using a two-step Runge-Kutta process. Hailstone terminal velocity is calculated and incorporated at both steps.  Following the advection, physical variables are interpolated at the new hailstone location and hailstone growth/melting is determined via the hailstone_trajectory subroutine.  Once the hailstone reaches one of the edges of the domain, hailstone_trajectory is no longer called and the hailstone remains there for the rest of the simulation.

Due to the nature of hailstones, some of the largest trajectory speeds occur near the surface.  I had to add a check for below ground heights within the 1st rk loop, not just after the 2nd as in the parcel code.

One possibility I haven't accounted for: I'm still allowing the hailstone to fall and be advected even if it has reached the edge of the domain. That could have an impact if it shot off the top, as it might eventually fall back down into the simulation and start growing again.  Thus far, stratospheric hailstones don't seem to be much of a problem.

Finally, **the implementation of this code for sigma vertical coordinates is not complete and should not be used.**

For additional documentation, see George Bryan’s README for his parcel code (README.parcel), or the main parts of the code itself in hail.F.


## Namelist options

### &param1:  
* hailfrq – Frequency to output hail data (s). The hail calculations are done every model timestep, but this value tells the model how frequently to output the information.  Set to a negative number to output every timestep.
 
### &param2:

* ihailtraj – turn on HAILCAST hail trajectories (1 = on)
* nhailtrajs – total number of hailstones.
   Note:  you’ll have to define all initial embryo locations and sizes in init3d.F.  Make sure the number you are defining there matches nhailtrajs.
 
### &param9:
* output_format – this is the generic output_format namelist selector. The hail trajectory data follows this selection.
 
### &param17:  (new section)

#### Hail-related parameters to include in the output:
* hail_dice: output ice-only diameter size
* hail_ts: output hailstone temperature
* hail_fw: output hailstone water fraction
* hail_vt: output hailstone terminal velocity
* hail_qice: output environmental cloud ice. (includes snow but no graupel)
* hail_qliq: output environment cloud liquid water (does not include rain)
* hail_tc: output environmental temperature
* hail_u: output environmental u wind
* hail_v: output environmental v wind
* hail_w: output environmental w wind 

Hail diameter, density, x, y, z location, and itype (if in wet or dry growth),  are always output.


## Run Instructions

I’ve included a sample namelist in the run directory as well as the header of a cm1out_haildata.nc file so you can get an idea of how the output data is structured.  I recommend keeping the model timestep (dtl) low, around a second, to avoid instabilities.

Depending on what you are simulating, it is likely you won't want to initialize the hail embryos right at the simulation start, like parcels are. I've played around with adding a delayed start option to the namelist, but recommend instead the following:
1. Compile the CM1 code with the number/locations of hail embryos you want declared in init3d.F.
2. Start the model run with _no_ hail embryos (ihailtraj=0), and run the simulation until you a realistic looking storm has developed. Write out a restart file.
3. Add the appropriate number of hail embryos/locations to the restart file, using the included python script add_hail_to_restart.py. Make sure these match what you coded in init3d.F.
4. Start the model again with ihailtraj=1 and all other hail-related namelist parameters set. Run the model until all hail has fallen out (no more than 45-60 min).
5. As desired, you can repeat steps 2-4 with incrementally later restart files. E.g., Run a no-hail simulation, write out a restart file at 60 min, add hail embryos, run for an hour. Write out a restart file from the no-hail simulation at 65 min, add hail embryos, run for an hour. Etc.

