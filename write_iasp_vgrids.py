#Read a csv file defining iasp91, write it out to vgrids format
import misc_tools
import numpy as np

#Load the parameters of the output vgrids file
outvel=misc_tools.Fmm_vgrids()
outvel.read()

#Read isap csv file
fnam='iasp.csv'
dep,rad,vp,vs=np.genfromtxt(fnam,delimiter=',').transpose()

#Build matrix from iasp
vvec=[]
radvec=np.arange(outvel.orad,outvel.orad+outvel.nrad)
for i in radvec:
    ind=misc_tools.find_nearest_index(rad,i)
    if i<=rad[ind]:
        vvec=np.append(vvec,vp[ind])
    else:
        vvec=np.append(vvec,vp[ind+1])
#Replace 3D velocity with iasp
for ilat in np.arange(outvel.nlat):
    for ilon in np.arange(outvel.nlon):
        outvel.vel[ilon][ilat]=vvec
outvel.write()
