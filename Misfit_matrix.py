import misc_tools
from mtools import *
from numpy import arange,asarray,delete
from antelope.stock import pfread, pfin
from matplotlib import pyplot as plt
params=pfin('eqloc3d.pf')
loc_params = params['location_parameters']
prop_params= params['propagation_grid']
#These should go in parameter file.
#nr = int(loc_params['nr'])
#nlat = int(loc_params['nlat'])
#nlon = int(loc_params['nlon'])
#nx, ny, nz = nlon, nlat, nr
earth_rad=6371

#Load events
fnam='pha.phase_210to9_july2011.dat'
pha=Phalist(fnam,'scedc')

#Get Propagation grid paramters
nlat=int(prop_params['nlat'])
nlon=int(prop_params['nlon'])
nz=int(prop_params['nr'])
li = misc_tools.Linear_index(nlon, nlat, nz)
olon=float(prop_params['minlon'])
olat=float(prop_params['minlat'])
oz=float(prop_params['minz'])
dlon=float(prop_params['dlon'])
dlat=float(prop_params['dlat'])
dz=  float(prop_params['dr'])

#Build vectors of geographic coordinates
qlon = arange(olon, dlon * nlon + olon, dlon)
qlat = arange(olat, dlat * nlat + olat, dlat)
qdep = arange(earth_rad-dz*nz,earth_rad,dz)+oz+dz
qlon = delete(qlon,-1,0)#Remove last element
qlat = delete(qlat,-1,0)
qdep = delete(qdep,-1,0)

ev=pha[0]
#Get station names for arrivals and a traveltime vector
absvec=[]
arrvec=[]
arrsta=[]        #a list of station names
for arrival in ev.arrivals:
    if arrival.phase is 'P':
        arrvec.append(arrival.ttime)
        absvec.append(arrival.epoch)
        arrsta.append(arrival.staname)
    if not os.path.isfile(arrival.staname+'traveltime'):
        continue
absvec=asarray(absvec)
arrvec=asarray(arrvec)

#Return the entire misfit matrix
dstep = int(loc_params['dstep2'])
dx, dy, dz = nlon / dstep, nlat / dstep, nz / dstep
dx,dy,dz=1,1,1
qx, qy, qz = range(0, nlon, dx), range(0, nlat, dy), range(0, nz, dz);
#minx, miny, minz, orgmin ,orimean,oristd= misc_tools.exp_grid_search(arrsta, qx, qy,qz, absvec, li,ev.epoch)
#orimean,oristd= misc_tools._grid_search_traveltimes_origin(arrsta, qx, qy,qz, absvec, li)
minx,miny,minz,orimean,oristd=misc_tools.exp_grid_search(arrsta,qx,qy,qz,absvec,li)
#Plot
idepth=misc_tools.find_nearest_index(earth_rad-qdep,ev.depth)
faultlon,faultlat=misc_tools.load_faults()
figid=plt.figure(1)
figid.clf();
plt.plot(faultlon,faultlat,'k')
plt.axis([min(qlon),max(qlon),min(qlat),max(qlat)])
pc_id=plt.pcolor(qlon,qlat,oristd[:,:,idepth] )
figid.colorbar(pc_id)
figid.show()
