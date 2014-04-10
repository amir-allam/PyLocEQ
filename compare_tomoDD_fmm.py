#Compare the traveltime calculations of FMM and TomoDD
from numpy import asarray
from misc_tools import *
from mtools import *
import scipy.io
import numpy as np
import math
import utm

#Load the model we'll use; this could be a class
mat = scipy.io.loadmat('combined_vel_models.mat')
C=mat['C'][0][0]
qx=C[0]
qy=C[1]
qz=C[8]
nx=uint(C[3][0][0])
ny=uint(C[4][0][0])
nz=uint(qz.size)
Vp=C[13]
Vs=C[14]
#Compute the origin and spacing of the grids in fucking RADIANS
earth_rad=6371.0
deg_to_rad=math.acos(-1)/180.0 #math.acos(-1) is pi
rad_to_deg=180.0/math.acos(-1)
olat,olon=utm.to_latlon(qx[0][0],qy[0][0],11,'S')
olat_rad=olat*deg_to_rad
olon_rad=olon*deg_to_rad
orad=earth_rad-qz[-1][0]
drad=(qz[1]-qz[0])[0] #Radius spacing
lat2,lon2=utm.to_latlon(qx[1][0],qy[0][0],11,'S')
lat3,lon3=utm.to_latlon(qx[0][0],qy[1][0],11,'S')
dlon=(lon2-olon)*deg_to_rad #lon spacing in RADIANS
dlat=(lat3-olat)*deg_to_rad #lat spacing in RADIANS
print 'olon,olat = %8.4f , %8.4f\n' %(olon,olat)
#Create qx,qy,qz vectors for the MOD file
qx_MOD=asarray(range(nx))*1.0
qy_MOD=asarray(range(ny))*1.0
qz_MOD=asarray([ii[0] for ii in qz])

#Parse the model into useful numpy vectors
outVp_FMM=np.empty((nx,ny,nz))
outVp_TomoDD=np.empty((nx,ny,nz))
for ix in range(nx):#Just be explicit and use brute force
    for iy in range(ny):
        for iz in range(nz):
            outVp_TomoDD[ix,iy,iz]=Vp[iy,ix,iz]
            outVp_FMM[ix,iy,iz]=Vp[iy,ix,nz-iz-1]

#Create and write out TomoDD MOD
MOD=TomoDD_MOD()
MOD.read() #Just creates a dummy MOD class
MOD.nx,MOD.ny,MOD.nz=nx,ny,nz
MOD.qx,MOD.qy,MOD.qz=qx_MOD,qy_MOD,qz_MOD
MOD.vel=outVp_TomoDD
MOD.velrat=outVp_TomoDD*0+1.73
MOD.write('poopMOD')
#Create and write out FMM vgrids.in
vg=Fmm_vgrids()
vg.read()
vg.nrad,vg.nlat,vg.nlon=nz,ny,nx
vg.drad,vg.dlat,vg.dlon=drad,dlat,dlon
vg.orad,vg.olat,vg.olon=orad,olat_rad,olon_rad
vg.vel=outVp_FMM
vg.write('poop_vgrids')
#Create and write interfaces.in
endlon=vg.olon+vg.dlon*(vg.nlon-1)
endlat=vg.olat+vg.dlat*(vg.nlat-1)
qlon=np.linspace(vg.olon,endlon,vg.nlon)*rad_to_deg
qlat=np.linspace(vg.olat,endlat,vg.nlat)*rad_to_deg

#topography up tup
topofile='/Users/aallam/gmt/100m_poop.grd'
topo=Interface()
topo.read_topo_netcdf(topofile)
topo.interp(qlon,qlat)
topo.z=topo.z/1000 #m to Km
topo.z=topo.z+earth_rad #elevation to radius
topo.write_fmm()

#Flat bottom
flat_z=vg.vel[:,:,0]*0+vg.orad+2 #Make interface the same size as vel
flat_bottom=Interface()
flat_bottom.set_data(qlon,qlat,flat_z)
flat_bottom.write_fmm(ifappend=True)

#Load station list and write to receivers.in
stafnam='sta.phase_220to9_july13.dat'
stalist=StationList(stafnam,False)
write_receivers_fmm(stalist)

#Run the two codes

#Read abs_tt_calc.dat

#Read FMM output

#Compare station by station


