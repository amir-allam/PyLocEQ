import time
import struct
from numpy import * #Should probably be more specific
from matplotlib import pyplot as plt
from array import array #This is for fast io when writing binary

#Read ASCII arrival time file output by FMM code

def num(s):
#Convert a string to a number, choosing int or float automatically
# SLOW, don't use for large lists
    try:
        return int(s)
    except ValueError:
        return float(s)

def load_faults():
#Load the california fault map
    fnam='cal_faults.dat'
    fid=open(fnam,'r')
    a=fid.readlines()
    faultlon=[];faultlat=[]
    for line in a:
        tmp=line.strip().split()
        if isnan(num(tmp[0])):
            faultlon.append(nan)
            faultlat.append(nan)
        else:
            faultlon.append(num(tmp[0]))
            faultlat.append(num(tmp[1]))
    print 'Faults loaded...'
    return faultlon,faultlat

class Fmm_vgrids():
# vgrids.in is the FMM velocity file format. Looks like:
#    ngrids ntypes
#    nradius nlat nlon     !number of points
#    dradius dlat dlon     !grid spacings; uniform in each direction
#    oradius olat olon     !origin of the grid
#    V(1,1,1)
#    V(2,1,1)
#    V(1,2,1)
#    V(1,1,2)
#    etc...
#
#       NOTE: Right-handed; oradius is somewhere deep in the earth
#
# ngrids is the number of volumes. We generally use 1
# ntypes will be 1 for just P, 2 for P and S
    def __init__(self):
        pass
    def read(self,fnam='vgrids.in'): #Create a matrix of velocities
        self.fnam=fnam
        fid = open(fnam,'r')
        tmp=fid.readline().strip().split()
        self.ngrids,self.ntypes=num(tmp[0]),num(tmp[1])
        tmp=fid.readline().strip().split()
        self.nrad,self.nlat,self.nlon=num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp=fid.readline().strip().split()
        self.drad,self.dlat,self.dlon=num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp=fid.readline().strip().split()
        self.orad,self.olat,self.olon=num(tmp[0]),num(tmp[1]),num(tmp[2])
        #Loop to create a numpy matrix
        self.vel=empty((self.nlon,self.nlat,self.nrad))
        for irad in range(self.nrad):
            for ilat in range(self.nlat):
                for ilon in range(self.nlon):
                    tmp=fid.readline().strip().split()
                    self.vel[ilon,ilat,irad]=float(tmp[0])
        #THERE SHOULD BE MORE STATEMENTS HERE IN CASE NTYPES,NGRIDS!= 1
        fid.close()
    def write(self,outfnam='out.vgrids'): #Write to vgrids format
        self.outfnam=outfnam
        fid=open(outfnam,'w')
        outs=str(self.ngrids)+' '+str(self.ntypes)+'\n'
        fid.write(outs)
        outs=str(self.nrad)+' '+str(self.nlat)+' '+str(self.nlon)+'\n'
        fid.write(outs)
        outs=str(self.drad)+' '+str(self.dlat)+' '+str(self.dlon)+'\n'
        fid.write(outs)
        outs=str(self.orad)+' '+str(self.olat)+' '+str(self.olon)+'\n'
        fid.write(outs)
        for ilon in range(self.nlon):
            for ilat in range(self.nlat):
                for irad in range(self.nrad):
                    fid.write('{1:.{0}f}'.format(3,self.vel[ilon,ilat,irad])+'\n')
        fid.close()

class Interface():
#A 2D matrix containing interface information, e.g., topography or moho depth
    def __init__(self):
        pass
    def set_data(self,x,y,z):
    #This builds the class ad-hoc. x,y,z are all numpy arrays. x,y are vectors, z is 2D
        self.x,self.y,self.z=x,y,z
        deg_to_rad=math.acos(-1)/180.0  #math.acos(-1) is pi
        self.nx,self.ny=len(x),len(y) #total number of x and y coordinates
        self.ox,self.oy=min(x),min(y)   #origins in lat,lon
        self.dx=self.x[1]-self.x[0]     #Spacings in degrees
        self.dy=self.y[1]-self.y[0]
        self.ox_rad=self.ox*deg_to_rad  #Origins in radians
        self.oy_rad=self.oy*deg_to_rad
        self.dx_rad=self.dx*deg_to_rad  #Spacings in radians
        self.dy_rad=self.dy*deg_to_rad
    def read_topo_netcdf(self,fnam,subsamp=1):
    #Read a netcdf formatted topography file
        from scipy.io import netcdf_file as netcdf
        topo=netcdf(fnam,'r')
        x=topo.variables['x'][::subsamp] #vector
        y=topo.variables['y'][::subsamp] #vector
        z=topo.variables['z'][::subsamp,::subsamp] #matrix
        self.subsamp=subsamp
        self.set_data(x,y,z)
    def interp(self,xi,yi):
    #Find the values at points xi,yi by cubic spline
        import scipy.interpolate
        gnn=scipy.interpolate.RectBivariateSpline(self.x,self.y,self.z.transpose())
        self.x_unint=self.x #Remember non-interpolated values
        self.y_unint=self.y
        self.z_unint=self.z
        zi=gnn.__call__(xi,yi) #the actual interpolation
        self.set_data(xi,yi,zi)
    def write_fmm(self,fnam='out_interfaces',ifappend=False):
    #Write to the fmm ascii format. ifappend should be set to True for any interface after the first; this will only write the z values and not the header
        if not ifappend:
            fid=open(fnam,'w')
            fid.write('2\n') # THIS ASSUMES THAT THERE ARE ONLY 2 INTERFACES
            fid.write('%u %u\n'%(self.ny,self.nx) )
            fid.write('%f %f\n'%(self.dy_rad,self.dx_rad) )
            fid.write('%f %f\n'%(self.oy_rad,self.ox_rad) )
        else:
            fid=open(fnam,'a')
        for iy in range(self.ny):#Write the interface as a vector
            outs=''.join('{1:.{0}f}\n'.format(3,ii)  for ii in self.z[:,iy])
            fid.write(outs)
        fid.close()

def find_containing_cube(px,py,pz,xvec,yvec,zvec):
#Find the 8 endpoints for the cell which contains point px,py
#  We take advantage of the regular grid
#  Assumes the point is inside the volume defined by xvec,yvec,zvec
#  Returns an array of size 8,3 where the rows contain x,y,z coordinates of the cubes endpoints
#  Also returns indexes of endpoints
    #Find the nearest node point and indexes <--"indices" sounds stupid to me
    xind,xnode=_find_nearest(px,xvec)
    yind,ynode=_find_nearest(py,yvec)
    zind,znode=_find_nearest(pz,zvec)
    #Now check if the 3 coordinates of p are greater or less than the node it is nearest
    if px>=xnode: #px is east of the nearest node
        xi=xind+1
        xn=xvec[xi]
    else:        #px is west of the nearest node
        xi=xind-1
        xn=xvec[xi]
    if py>=ynode: #px is north of the nearest node
        yi=yind+1
        yn=yvec[yi]
    else:        #px is south of the nearest node
        yi=yind-1
        yn=yvec[yi]
    if pz>=znode: #px is above the nearest node
        zi=zind+1
        zn=zvec[zi]
    else:        #px is below the nearest node
        zi=zind-1
        zn=zvec[zi]
    #Add new endpoints to define the cube
    endpoints=[]
    endpoints.append( [xnode,ynode,znode] )
    endpoints.append( [xn,ynode,znode] )
    endpoints.append( [xn,yn,znode]    )
    endpoints.append( [xnode,yn,znode] )
    endpoints.append( [xnode,ynode,zn] )
    endpoints.append( [xn,ynode,zn] )
    endpoints.append( [xn,yn,zn]    )
    endpoints.append( [xnode,yn,zn] )
    #Add indices
    indexes=[]
    indexes.append( [xind,yind,zind] )
    indexes.append( [xi,yind,zind] )
    indexes.append( [xi,yi,zind]    )
    indexes.append( [xind,yi,zind] )
    indexes.append( [xind,yind,zi] )
    indexes.append( [xi,yind,zi] )
    indexes.append( [xi,yi,zi]    )
    indexes.append( [xind,yi,zi] )

    return endpoints,indexes

def find_nearest(nparray,value):
    #Returns the nearest item in nparray to value
    idx= (abs(nparray-value)).argmin()
    return nparray.flat[idx]

def find_nearest_index(nparray,value):
    idx= (abs(nparray-value)).argmin()
    return idx

def _find_nearest(px,xvec):
#Find the nearest x in xvec 
#  returns index
    best_ind=0
    shortest=100000000.0;
    for ii in range(len(xvec)):
        if abs(xvec[ii]-px)<shortest:
            shortest=abs(xvec[ii]-px)
            best_ind=ii
    return best_ind,xvec[best_ind]

def read_binary_float(fid,n=0,precision='double'):
#read the nth float value from a binary file at the with the given precision
# following python conventions, 0 is the index of the first value
    if precision is 'single':
        numbytes=4;packstr='f'
    else:
        numbytes=8;packstr='d'
    offset=n*numbytes #Go to the right spot in the file
    fid.seek(offset)
    raw=fid.read(numbytes)
    tmp=struct.unpack(packstr,raw)
    val=tmp[0]
    return val

def grid_search_traveltimes_rms(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of some criterion by performing a grid search
# We aren't necessarily searching the whole grid; we may be skipping values
#   sta         list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    rms=array([])
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in qx:  #Loop over the three vectors, searching every point
        for iy in qy:
            for iz in qz:
                calctt=array([]); #initialize the calculated tt vector
                ind=li.get_1D(ix,iy,iz) #Find the vector index
                for sta in arrsta: #Build vector of calculated ttimes
                    #print 'bin.'+sta+'.traveltime'
                    fid = open('bin.'+sta+'.traveltime')
                    calctt=append(calctt, read_binary_float(fid,ind) )
                    #print read_binary_float(fid,ind)
                    fid.close()
                rms=append(rms, sqrt(mean( (arrvec-calctt)**2 )) ) #Root-mean-square, yo
    #Now find the 3D index of the best point so far
    min_ind=rms.argmin()
    (minx,miny,minz)=search_inds.get_3D(min_ind)
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return minx,miny,minz

def grid_search_traveltimes_origin(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
#  sta          list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of absolute arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    origin_std=array([])
    origin_mean=array([])
    #origin_std=empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
    #origin_mean=empty([ len(qy),len(qx),len(qz)] )
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in qx: #range(len(qx)):  #Loop over the three vectors, searching every point
        print 'On ix: ', ix,'/',len(qx),'\n'
        for iy in qy: #range(len(qy)):
            for iz in qz: #range(len(qz)):
                calctt=array([]); #initialize the calculated tt vector
                #ind=li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                ind=li.get_1D(ix,iy,iz)
                calctt=read_tt_vector(arrsta,ind) #Make a traveltime vector from calculated times
                orivec=arrvec-calctt #Take the difference
                origin_mean=append( origin_mean,orivec.mean()) #find the mean origin time
                if min(calctt)<0: #If the traveltime <0, this gridpoint is null
                    origin_std=append(origin_std,10000) #A large dummy value
                    continue
                origin_std=append(origin_std,orivec.std)
                #origin_std[iy,ix,iz]=orivec.std()
                #origin_mean[iy,ix,iz]=orivec.mean()
    #Now find the 3D index of the best point so far
    #new_origin=origin_mean.flatten()[oristd.argmin()] #The origin time with lowest std
    min_ind=origin_std.argmin()
    (minx,miny,minz)=search_inds.get_3D(min_ind)
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return minx,miny,minz,origin_mean[min_ind]

def _grid_search_traveltimes_origin(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
#  sta          list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of absolute arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    rms=array([])
    #origin_std=array([])
    origin_std=empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
    origin_mean=empty([ len(qy),len(qx),len(qz)] )
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in range(len(qx)):  #Loop over the three vectors, searching every point
        for iy in range(len(qy)):
            for iz in range(len(qz)):
                calctt=array([]); #initialize the calculated tt vector
                ind=li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                for sta in arrsta: #Build vector of calculated ttimes
                    fid = open('bin.'+sta+'.traveltime')
                    tmp=read_binary_float(fid,ind)
                    #if tmp<0: tmp=NaN #Kill negative travel times
                    calctt=append(calctt, tmp )
                    fid.close()
                orivec=arrvec-calctt
                origin_std[iy,ix,iz]=orivec.std()
                origin_mean[iy,ix,iz]=orivec.mean()
    #Now find the 3D index of the best point so far
    min_ind=origin_std.argmin()
    (minx,miny,minz)=search_inds.get_3D(min_ind)
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return origin_mean,origin_std

def exp_grid_search(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
#  sta          list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of absolute arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array,indices #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    #origin_std=array([])
    #origin_mean=array([])
    origin_std=empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
    origin_mean=empty([ len(qy),len(qx),len(qz)] )
    search_inds=Linear_index(len(qx),len(qy),len(qz))
    for ix in range(len(qx)):  #Loop over the three vectors, searching every point
        print 'On ix: ', ix,'/',len(qx),'\n'
        for iy in range(len(qy)):
            for iz in range(len(qz)):
                calctt=array([]); #initialize the calculated tt vector
                ind=li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                calctt=read_tt_vector(arrsta,ind) #Make a traveltime vector from calculated times
                orivec=arrvec-calctt #Take the difference
                if min(calctt)<0: #If the traveltime <0, this gridpoint is null
                    continue
                origin_std[iy,ix,iz]=orivec.std()
                origin_mean[iy,ix,iz]=orivec.mean()
    #Now find the 3D index of the best point so far
    #new_origin=origin_mean.flatten()[oristd.argmin()] #The origin time with lowest std
    besti=origin_std.argmin() #1D index of best point
    helper_i=indices(origin_std.shape) #3D indices
    miny= helper_i[0,:,:,:].flatten()[besti]
    minx= helper_i[1,:,:,:].flatten()[besti]
    minz= helper_i[2,:,:,:].flatten()[besti]
    minx=qx[minx]; miny=qy[miny]; minz=qz[minz];
    return minx,miny,minz,origin_mean,origin_std

def locate_eq(ev):
    #Locate an earthquake based on the arrivals in ev, traveltime files which are already saved
    from antelope.stock import pfread, pfin
    from mtools import Phalist,Stalist
    params=pfin('eqloc3d.pf')
    loc_params = params['location_parameters']
    prop_params= params['propagation_grid']
    earth_rad=6371.0

    #Get Propagation grid paramters
    nlat=int(prop_params['nlat'])
    nlon=int(prop_params['nlon'])
    nz=int(prop_params['nr'])
    li = Linear_index(nlon, nlat, nz)
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

    #Grid search for best location
    start_time=time.time()
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
    if len(arrvec)<6: #About this many phases are needed to get a decent result
        return None
    absvec=asarray(absvec)
    arrvec=asarray(arrvec)
    #print 'Number of phases used: ',len(arrvec)

    #Search coarsely
    #dstep should go in parameter file.
    dstep = int(loc_params['dstep2'])
    dx, dy, dz = nlon / dstep, nlat / dstep, nz / dstep
    dx, dy, dz = 1,1,1 #Remove this later
    qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nz, dz);
    minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
                                                                qz, absvec, li)
    #minx, miny, minz, orgmin = exp_grid_search(arrsta,qx,qy,qz, absvec, li,ev.epoch)
    #Finer search
#    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
#                                                            arrvec,li)
    buff = int(loc_params['buff2'])
    qx = range(minx - buff, minx + buff)
    qy = range(miny - buff, miny + buff)
    qz = range(minz - buff, minz + buff);
    qx = fix_boundary_search(qx, li.nx)
    qy = fix_boundary_search(qy, li.ny)
    qz = fix_boundary_search(qz, li.nz)
    #minx, miny, minz, orgmin = grid_search_traveltimes_origin(arrsta, qx, qy,
    #                                                            qz, absvec, li)
    #minx, miny, minz, orgmin = exp_grid_search(arrsta,qx,qy,qz, absvec, li,ev.epoch)
#    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
#                                                            arrvec,li)
    orgmin=0

    #Find the best subgrid location
    c,resid=get_subgrid_loc(minx,miny,minz,arrvec,arrsta,li)
    delta_x=qlon[1]-qlon[0]
    delta_y=qlat[1]-qlat[0]
    delta_z=qdep[1]-qdep[0]
    loc_change=c*[delta_x,delta_y,delta_z] #Subgrid location change in lon/lat/depth

    #Find the best-fit source location in geographic coordinates
    newloc=[newlon,newlat,newz]=[ qlon[minx],qlat[miny],qdep[minz] ] #+loc_change
    elapsed_time=time.time()-start_time
    print ev.id,len(arrvec),newlon,newlat,earth_rad-newz,ev.lon,ev.lat,ev.depth,ev.epoch-orgmin,elapsed_time,resid
    #This function will return location, origin time, and error estimates for both

def get_subgrid_loc(ix,iy,iz,arrvec,arrsta,li):
    #Test least squares on real data
    import numpy as np
    from scipy import linalg
    import matplotlib.pyplot as plt

    #Get traveltime vectors for the closest point and its neighbors
    ind=li.get_1D(ix,iy,iz)
    tt000= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix+1,iy,iz)
    tt100= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy+1,iz)
    tt010= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy,iz+1)
    tt001= read_tt_vector(arrsta,ind)
    #backwards
    ind=li.get_1D(ix-1,iy,iz)
    btt100= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy-1,iz)
    btt010= read_tt_vector(arrsta,ind)
    ind=li.get_1D(ix,iy,iz-1)
    btt001= read_tt_vector(arrsta,ind)

    #Calculate forward derivatives
    dt_dx=tt100-tt000
    dt_dy=tt010-tt000
    dt_dz=tt001-tt000
    #backwards
    bdt_dx=tt000-btt100
    bdt_dy=tt000-btt010
    bdt_dz=tt000-btt001
    #Build A matrix and r in r=Ax  (x is the spatial vector here [x,y,z]
    A=c_[dt_dx,dt_dy,dt_dz]
    r=arrvec-tt000
    c,resid,rank,sigma=linalg.lstsq(A,r)
    #backwards
    bA=c_[bdt_dx,bdt_dy,bdt_dz]
    bc,resid,rank,sigma=linalg.lstsq(bA,r)
    return c,resid



def read_tt_vector(stanames,ind):
    #Read from binary files to create a vector of traveltimes for the stations in stanames
    #   at the index location ind, which is the 1D index
    #   stanames is a list of station names only
    from numpy import array
    ttvec=array([])
    for sta in stanames:
        fid=open('bin.'+sta+'.traveltime')
        ttvec=append(ttvec, read_binary_float(fid,ind) )
        fid.close()
    return ttvec

def fix_boundary_search(qx,nx):
#When performing a grid search on a subgrid, make sure you don't go off the edges
#  qx         search vectors, these will be modified then returned
#  nx         max index [li.nx]
    for ix in range(len(qx)):
        if qx[ix]<0:
            qx[ix]=0
        if qx[ix]>=nx:
            qx[ix]=nx-1
    newqx=uniq(qx)
    return newqx

def uniq(input):
#Remove duplicate items from a list. Preserves order.
  output = []
  for x in input:
    if x not in output:
      output.append(x)
  return output

class Linear_index():
    #Holds a 1D list of 3D indices and a 3D list of 1D indices
    # where iz varies fastest, then iy, then ix
    #  The speed of this can certainly be improved
    def __init__(self,nx,ny,nz):
        from numpy import empty
        self.nx=nx; self.ny=ny; self.nz=nz;
        self.i1D=[]
        self.i3D=empty((nx,ny,nz)) #python has weird index conventions
        ic=0
        for ix in range(nx): #Fuck it, just be explicit
            for iy in range(ny):
                for iz in range(nz):
                    self.i1D.append( (ix,iy,iz) )
                    self.i3D[ix,iy,iz]=ic
                    ic=ic+1
        self.i3D=self.i3D.astype(int)
    def get_1D(self,ix,iy,iz):
        return self.i3D[ix,iy,iz]
    def get_3D(self,iv):
        return self.i1D[iv]
