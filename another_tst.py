import misc_tools
from mtools import *
from numpy import arange,asarray
from antelope.stock import pfread, pfin
params=pfin('eqloc3d.pf')
loc_params = params['location_parameters']
#These should go in parameter file.
#nr = int(loc_params['nr'])
#nlat = int(loc_params['nlat'])
#nlon = int(loc_params['nlon'])
#nx, ny, nz = nlon, nlat, nr
earth_rad=6371

#Load events
fnam='pha.phase_210to9_july2011.dat'
pha=Phalist(fnam,'scedc')

for ev in pha:
    if ev.lon<-117 or ev.lon>-116 or ev.lat<33.0 or ev.lat>34.0:
        continue
    misc_tools.locate_eq(ev)
