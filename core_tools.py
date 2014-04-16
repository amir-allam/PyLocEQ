import os
import time
from numpy import * #Should probably be more specific
#from matplotlib import pyplot as plt
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
    fnam = 'cal_faults.dat'
    fid = open(fnam,'r')
    a = fid.readlines()
    faultlon = [];faultlat = []
    for line in a:
        tmp = line.strip().split()
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
        self.fnam = fnam
        fid = open(fnam,'r')
        tmp = fid.readline().strip().split()
        self.ngrids,self.ntypes = num(tmp[0]),num(tmp[1])
        tmp = fid.readline().strip().split()
        self.nrad,self.nlat,self.nlon = num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp = fid.readline().strip().split()
        self.drad,self.dlat,self.dlon = num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp = fid.readline().strip().split()
        self.orad,self.olat,self.olon = num(tmp[0]),num(tmp[1]),num(tmp[2])
        #Loop to create a numpy matrix
        self.vel = empty((self.nlon,self.nlat,self.nrad))
        for irad in range(self.nrad):
            for ilat in range(self.nlat):
                for ilon in range(self.nlon):
                    tmp = fid.readline().strip().split()
                    self.vel[ilon,ilat,irad] = float(tmp[0])
        #THERE SHOULD BE MORE STATEMENTS HERE IN CASE NTYPES,NGRIDS!= 1
        fid.close()
    def write(self,outfnam='out.vgrids'): #Write to vgrids format
        self.outfnam = outfnam
        fid = open(outfnam,'w')
        outs = str(self.ngrids)+' '+str(self.ntypes)+'\n'
        fid.write(outs)
        outs = str(self.nrad)+' '+str(self.nlat)+' '+str(self.nlon)+'\n'
        fid.write(outs)
        outs = str(self.drad)+' '+str(self.dlat)+' '+str(self.dlon)+'\n'
        fid.write(outs)
        outs = str(self.orad)+' '+str(self.olat)+' '+str(self.olon)+'\n'
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
        self.x,self.y,self.z = x,y,z
        deg_to_rad = math.acos(-1)/180.0  #math.acos(-1) is pi
        self.nx,self.ny = len(x),len(y) #total number of x and y coordinates
        self.ox,self.oy = min(x),min(y)   #origins in lat,lon
        self.dx = self.x[1]-self.x[0]     #Spacings in degrees
        self.dy = self.y[1]-self.y[0]
        self.ox_rad = self.ox*deg_to_rad  #Origins in radians
        self.oy_rad = self.oy*deg_to_rad
        self.dx_rad = self.dx*deg_to_rad  #Spacings in radians
        self.dy_rad = self.dy*deg_to_rad
    def read_topo_netcdf(self,fnam,subsamp=1):
    #Read a netcdf formatted topography file
        from scipy.io import netcdf_file as netcdf
        topo = netcdf(fnam,'r')
        x = topo.variables['x'][::subsamp] #vector
        y = topo.variables['y'][::subsamp] #vector
        z = topo.variables['z'][::subsamp,::subsamp] #matrix
        self.subsamp = subsamp
        self.set_data(x,y,z)
    def interp(self,xi,yi):
    #Find the values at points xi,yi by cubic spline
        import scipy.interpolate
        gnn = scipy.interpolate.RectBivariateSpline(self.x,self.y,self.z.transpose())
        self.x_unint = self.x #Remember non-interpolated values
        self.y_unint = self.y
        self.z_unint = self.z
        zi = gnn.__call__(xi,yi) #the actual interpolation
        self.set_data(xi,yi,zi)
    def write_fmm(self,fnam='out_interfaces',ifappend=False):
    #Write to the fmm ascii format. ifappend should be set to True for any interface after the first; this will only write the z values and not the header
        if not ifappend:
            fid = open(fnam,'w')
            fid.write('2\n') # THIS ASSUMES THAT THERE ARE ONLY 2 INTERFACES
            fid.write('%u %u\n'%(self.ny,self.nx) )
            fid.write('%f %f\n'%(self.dy_rad,self.dx_rad) )
            fid.write('%f %f\n'%(self.oy_rad,self.ox_rad) )
        else:
            fid = open(fnam,'a')
        for iy in range(self.ny):#Write the interface as a vector
            outs = ''.join('{1:.{0}f}\n'.format(3,ii)  for ii in self.z[:,iy])
            fid.write(outs)
        fid.close()

def find_containing_cube(px,py,pz,xvec,yvec,zvec):
#Find the 8 endpoints for the cell which contains point px,py
#  We take advantage of the regular grid
#  Assumes the point is inside the volume defined by xvec,yvec,zvec
#  Returns an array of size 8,3 where the rows contain x,y,z coordinates of the cubes endpoints
#  Also returns indexes of endpoints
    #Find the nearest node point and indexes <--"indices" sounds stupid to me
    xind,xnode = _find_nearest(px,xvec)
    yind,ynode = _find_nearest(py,yvec)
    zind,znode = _find_nearest(pz,zvec)
    #Now check if the 3 coordinates of p are greater or less than the node it is nearest
    if px>= xnode: #px is east of the nearest node
        xi = xind+1
        xn = xvec[xi]
    else:        #px is west of the nearest node
        xi = xind-1
        xn = xvec[xi]
    if py >= ynode: #px is north of the nearest node
        yi = yind+1
        yn = yvec[yi]
    else:        #px is south of the nearest node
        yi = yind-1
        yn = yvec[yi]
    if pz >= znode: #px is above the nearest node
        zi = zind+1
        zn = zvec[zi]
    else:        #px is below the nearest node
        zi = zind-1
        zn = zvec[zi]
    #Add new endpoints to define the cube
    endpoints = []
    endpoints.append( [xnode,ynode,znode] )
    endpoints.append( [xn,ynode,znode] )
    endpoints.append( [xn,yn,znode]    )
    endpoints.append( [xnode,yn,znode] )
    endpoints.append( [xnode,ynode,zn] )
    endpoints.append( [xn,ynode,zn] )
    endpoints.append( [xn,yn,zn]    )
    endpoints.append( [xnode,yn,zn] )
    #Add indices
    indexes = []
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
    idx = (abs(nparray-value)).argmin()
    return nparray.flat[idx]

def find_nearest_index(nparray,value):
    idx = (abs(nparray-value)).argmin()
    return idx

def _find_nearest(px,xvec):
#Find the nearest x in xvec 
#  returns index
    best_ind = 0
    shortest = 100000000.0;
    for ii in range(len(xvec)):
        if abs(xvec[ii]-px)<shortest:
            shortest = abs(xvec[ii]-px)
            best_ind = ii
    return best_ind,xvec[best_ind]

def read_binary_float(fid,n=0,precision='double'):
    import struct
#read the nth float value from a binary file at the with the given precision
# following python conventions, 0 is the index of the first value
    if precision is 'single':
        numbytes = 4;packstr = 'f'
    else:
        numbytes = 8;packstr = 'd'
    offset = n*numbytes #Go to the right spot in the file
    fid.seek(offset)
    raw = fid.read(numbytes)
    tmp = struct.unpack(packstr,raw)
    val = tmp[0]
    return val



def _grid_search_traveltimes_origin(arrsta,qx,qy,qz,arrvec,li):
#Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
#  sta          list of station names; strings
#   qx,qy,qz    vectors of indices to search through
#   arrvec      vector of absolute arrivals in the same order as sta
#   li          Linear_index class for the entire traveltime grid
    from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
    rms = array([])
    #origin_std = array([])
    origin_std = empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
    origin_mean = empty([ len(qy),len(qx),len(qz)] )
    search_inds = Linear_index(len(qx),len(qy),len(qz))
    for ix in range(len(qx)):  #Loop over the three vectors, searching every point
        for iy in range(len(qy)):
            for iz in range(len(qz)):
                calctt = array([]); #initialize the calculated tt vector
                ind = li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                for sta in arrsta: #Build vector of calculated ttimes
                    fid = open('bin.'+sta+'.traveltime')
                    tmp = read_binary_float(fid,ind)
                    #if tmp<0: tmp = NaN #Kill negative travel times
                    calctt = append(calctt, tmp )
                    fid.close()
                orivec = arrvec-calctt
                origin_std[iy,ix,iz] = orivec.std()
                origin_mean[iy,ix,iz] = orivec.mean()
    #Now find the 3D index of the best point so far
    min_ind = origin_std.argmin()
    (minx,miny,minz) = search_inds.get_3D(min_ind)
    minx = qx[minx]; miny = qy[miny]; minz = qz[minz];
    return origin_mean,origin_std


class Locator:
    def __init__(self, config_file):
        import ConfigParser
        config = ConfigParser.RawConfigParser()
        config.read(config_file)
        for section in config.sections():
            tmp_dict = {}
            for option in config.options(section):
                tmp_dict[option] = config.get(section, option)
            setattr(self, section, tmp_dict)

    def locate_eq(self, ev):
        #Locate an earthquake based on the arrivals in ev, traveltime files which are already saved
        from temp_core_tools import Origin #This line will be deleted once things are reorganized.
        loc_params = self.location_parameters
        prop_params = self.propagation_grid
        earth_rad = 6371.0

        #Get Propagation grid paramters
        nlat = int(prop_params['nlat'])
        nlon = int(prop_params['nlon'])
        nz = int(prop_params['nr'])
        li  =  Linear_index(nlon, nlat, nz)
        olon = float(prop_params['minlon'])
        olat = float(prop_params['minlat'])
        oz = float(prop_params['minz'])
        dlon = float(prop_params['dlon'])
        dlat = float(prop_params['dlat'])
        dz = float(prop_params['dr'])

        #Build vectors of geographic coordinates
        qlon = arange(olon, dlon * nlon + olon, dlon)
        qlat = arange(olat, dlat * nlat + olat, dlat)
        qdep = arange(earth_rad - dz*nz,earth_rad,dz)+oz+dz

        #Grid search for best location
        start_time = time.time()
        absvec = []
        arrvec = []
        arrsta = []        #a list of station names
        arrpha = []        #List of phases
        for arrival in ev.arrivals:
            if arrival.phase is 'P':
                arrvec.append(arrival.time - ev.time)
                absvec.append(arrival.time)
                arrsta.append(arrival.sta)
                arrpha.append(arrival.phase)
            if not os.path.isfile(arrival.sta + 'traveltime'):
                continue
        if len(arrvec)<6: #About this many phases are needed to get a decent result
            return None
        absvec = asarray(absvec)
        arrvec = asarray(arrvec)
        #print 'Number of phases used: ',len(arrvec)

        #Search coarsely
        #dstep should go in parameter file.
        dstep = int(loc_params['dstep2'])
        dx, dy, dz = nlon / dstep, nlat / dstep, nz / dstep
        #dx, dy, dz = 1,1,1 #Remove this later
        qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nz, dz);
        #minx, miny, minz, orgmin = self.grid_search_traveltimes_origin(arrsta, qx, qy,qz, absvec, li)
        minx, miny, minz, origin_mean, origin_std= self.exp_grid_search(arrsta,qx,qy,qz, absvec, li)
        #Finer search
    #    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,
    #                                                            arrvec,li)
        buff = int(loc_params['buff2'])
        qx = range(minx - buff, minx + buff)
        qy = range(miny - buff, miny + buff)
        qz = range(minz - buff, minz + buff);
        qx = self.fix_boundary_search(qx, li.nx)
        qy = self.fix_boundary_search(qy, li.ny)
        qz = self.fix_boundary_search(qz, li.nz)
        #minx, miny, minz, orgmin = self.grid_search_traveltimes_origin(arrsta, qx, qy,qz, absvec, li)
        minx, miny, minz, origin_mean, origin_std= self.exp_grid_search(arrsta,qx,qy,qz, absvec, li)
    #    minx, miny, minz = grid_search_traveltimes_rms(arrsta, qx, qy, qz,arrvec,li)
        #orgmin = 0

        #Find the best subgrid location
        c,resid,tt_updated = self.get_subgrid_loc(minx,miny,minz,arrvec,arrsta,li)
        delta_x = qlon[1] - qlon[0]
        delta_y = qlat[1] - qlat[0]
        delta_z = qdep[1] - qdep[0]
        loc_change = c*[delta_x,delta_y,delta_z] #Subgrid location change in lon/lat/depth

        #Add calculated travel times to Origin
        ic=0 #Just a counter
        for arrival in ev.arrivals: #Loop over all input phases
            if (arrsta[ic] is arrival.sta) and (arrpha[ic] is arrival.phase):
                arrival.tt_calc=tt_updated[ic] #Only update the ones used in the relocation
            ic+=1

        #Find the best-fit source location in geographic coordinates
        newloc = [newlon,newlat,newz] = [ qlon[minx],qlat[miny],qdep[minz] ] #+loc_change
        best_otime=origin_mean.flatten()[origin_std.argmin()]
        elapsed_time = time.time() - start_time
        print ev.orid,len(arrvec),newlon,newlat,earth_rad - newz,ev.lon,ev.lat,ev.depth,ev.time -best_otime,elapsed_time,resid
        return Origin(newlat,
                      newlon,
                      earth_rad - newz,
                      best_otime,
                      'PyLocEQ',
                      arrivals=ev.arrivals,
                      evid=ev.evid)
        #This function will return location, origin time, and error estimates for both

    def grid_search_traveltimes_origin(self, arrsta, qx, qy, qz, arrvec, li):
    #Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
    #  sta          list of station names; strings
    #   qx,qy,qz    vectors of indices to search through
    #   arrvec      vector of absolute arrivals in the same order as sta
    #   li          Linear_index class for the entire traveltime grid
        from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
        origin_std = array([])
        origin_mean = array([])
        #origin_std = empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
        #origin_mean = empty([ len(qy),len(qx),len(qz)] )
        search_inds = Linear_index(len(qx),len(qy),len(qz))
        for ix in qx: #range(len(qx)):  #Loop over the three vectors, searching every point
            #print 'On ix: ', ix,'/',len(qx),'\n'
            for iy in qy: #range(len(qy)):
                for iz in qz: #range(len(qz)):
                    calctt = array([]); #initialize the calculated tt vector
                    #ind = li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                    ind = li.get_1D(ix,iy,iz)
                    calctt = self.read_tt_vector(arrsta,ind) #Make a traveltime vector from calculated times
                    orivec = arrvec - calctt #Take the difference
                    origin_mean = append( origin_mean,orivec.mean()) #find the mean origin time
                    if min(calctt)<0: #If the traveltime <0, this gridpoint is null
                        origin_std = append(origin_std,10000) #A large dummy value
                        continue
                    origin_std = append(origin_std,orivec.std)
                    #origin_std[iy,ix,iz] = orivec.std()
                    #origin_mean[iy,ix,iz] = orivec.mean()
        #Now find the 3D index of the best point so far
        #new_origin = origin_mean.flatten()[oristd.argmin()] #The origin time with lowest std
        min_ind = origin_std.argmin()
        (minx,miny,minz) = search_inds.get_3D(min_ind)
        minx = qx[minx]; miny = qy[miny]; minz = qz[minz];
        return minx,miny,minz,origin_mean[min_ind]

    def grid_search_traveltimes_rms(self, arrsta, qx, qy, qz, arrvec, li):
    #Find the minimum value of some criterion by performing a grid search
    # We aren't necessarily searching the whole grid; we may be skipping values
    #   sta         list of station names; strings
    #   qx,qy,qz    vectors of indices to search through
    #   arrvec      vector of arrivals in the same order as sta
    #   li          Linear_index class for the entire traveltime grid
        from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
        rms = array([])
        search_inds = Linear_index(len(qx),len(qy),len(qz))
        for ix in qx:  #Loop over the three vectors, searching every point
            for iy in qy:
                for iz in qz:
                    calctt = array([]); #initialize the calculated tt vector
                    ind = li.get_1D(ix,iy,iz) #Find the vector index
                    for sta in arrsta: #Build vector of calculated ttimes
                        #print 'bin.'+sta+'.traveltime'
                        fid = open('bin.'+sta+'.traveltime')
                        calctt = append(calctt, read_binary_float(fid,ind) )
                        #print read_binary_float(fid,ind)
                        fid.close()
                    rms = append(rms, sqrt(mean( (arrvec - calctt)**2 )) ) #Root-mean-square, yo
        #Now find the 3D index of the best point so far
        min_ind = rms.argmin()
        (minx,miny,minz) = search_inds.get_3D(min_ind)
        minx = qx[minx]; miny = qy[miny]; minz = qz[minz];
        return minx,miny,minz

    def read_tt_vector(self, stanames, ind):
        #Read from binary files to create a vector of traveltimes for the stations in stanames
        #   at the index location ind, which is the 1D index
        #   stanames is a list of station names only
        from numpy import array
        ttvec = array([])
        ttdir=self.misc['tt_map_dir']
        for sta in stanames:
            fid = open(ttdir+'bin.'+sta+'.traveltime')
            ttvec = append(ttvec, read_binary_float(fid,ind) )
            fid.close()
        return ttvec

    def get_subgrid_loc(self, ix, iy, iz, arrvec, arrsta, li):
        #Test least squares on real data
        import numpy as np
        from scipy import linalg
        import matplotlib.pyplot as plt

        #Cut off derivative calculations at model boundaries
        endx=ix+1; endy=iy+1; endz=iz+1
        if endx==li.nx: endx=ix 
        if endy==li.ny: endy=iy 
        if endz==li.nz: endz=iz 
        #Get traveltime vectors for the closest point and its neighbors
        ind = li.get_1D(ix,iy,iz)
        tt000 = self.read_tt_vector(arrsta, ind)
        ind = li.get_1D(endx, iy, iz)
        tt100 = self.read_tt_vector(arrsta, ind)
        ind = li.get_1D(ix, endy, iz)
        tt010 = self.read_tt_vector(arrsta, ind)
        ind = li.get_1D(ix, iy, endz)
        tt001 = self.read_tt_vector(arrsta, ind)

        #Calculate forward derivatives
        dt_dx = tt100 - tt000
        dt_dy = tt010 - tt000
        dt_dz = tt001 - tt000
        #Build A matrix and r in r=Ax  (x is the spatial vector here [x,y,z]
        A = c_[dt_dx, dt_dy, dt_dz]
        r = arrvec - tt000
        c, resid, rank, sigma = linalg.lstsq(A, r)

        #Compute updated travel times
        tt_updated=tt000+(A*c).sum(axis=1) #c has independent changes for x,y,z, so sum them
        return c, resid,tt_updated

    def fix_boundary_search(self, qx, nx):
    #When performing a grid search on a subgrid, make sure you don't go off the edges
    #  qx         search vectors, these will be modified then returned
    #  nx         max index [li.nx]
        for ix in range(len(qx)):
            if qx[ix] < 0:
                qx[ix] = 0
            if qx[ix] >= nx:
                qx[ix] = nx - 1
        newqx = uniq(qx)
        return newqx

    def exp_grid_search(self, arrsta, qx, qy, qz, arrvec, li):
    #Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
    #  sta          list of station names; strings
    #   qx,qy,qz    vectors of indices to search through
    #   arrvec      vector of absolute arrivals in the same order as sta
    #   li          Linear_index class for the entire traveltime grid
        from numpy import array, indices #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
        #origin_std = array([])
        #origin_mean = array([])
        origin_std = empty([len(qy), len(qx), len(qz)]) + 1000 #Give large starting values
        origin_mean = empty([len(qy), len(qx), len(qz)])
        search_inds = Linear_index(len(qx), len(qy), len(qz))
        for ix in range(len(qx)):  #Loop over the three vectors, searching every point
            print 'On ix: ', ix,'/',len(qx),'\n'
            for iy in range(len(qy)):
                for iz in range(len(qz)):
                    calctt = array([]) #initialize the calculated tt vector
                    ind  = li.get_1D(qx[ix], qy[iy], qz[iz]) #Find the vector index
                    calctt = self.read_tt_vector(arrsta,ind) #Make a traveltime vector from calculated times
                    orivec = arrvec - calctt #Take the difference
                    if min(calctt) < 0: #If the traveltime <0, this gridpoint is null
                        continue
                    origin_std[iy, ix, iz] = orivec.std()
                    origin_mean[iy, ix, iz] = orivec.mean()
        #Now find the 3D index of the best point so far
        #new_origin = origin_mean.flatten()[oristd.argmin()] #The origin time with lowest std
        besti = origin_std.argmin() #1D index of best point
        helper_i = indices(origin_std.shape) #3D indices
        miny = helper_i[0, :, :, :].flatten()[besti]
        minx = helper_i[1, :, :, :].flatten()[besti]
        minz = helper_i[2, :, :, :].flatten()[besti]
        minx = qx[minx]
        miny = qy[miny]
        minz = qz[minz];
        return minx, miny, minz, origin_mean, origin_std

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
        self.nx = nx; self.ny = ny; self.nz = nz;
        self.i1D = []
        self.i3D = empty((nx,ny,nz)) #python has weird index conventions
        ic = 0
        for ix in range(nx): #Fuck it, just be explicit
            for iy in range(ny):
                for iz in range(nz):
                    self.i1D.append( (ix,iy,iz) )
                    self.i3D[ix,iy,iz] = ic
                    ic = ic+1
        self.i3D = self.i3D.astype(int)
    def get_1D(self,ix,iy,iz):
        return self.i3D[ix,iy,iz]
    def get_3D(self,iv):
        return self.i1D[iv]

################################################
class Station:
    """
    A containter class for station location data.
    """
    def __init__(self, sta, lat, lon, elev):
        self.sta = sta
        self.lat = lat
        self.lon = lon
        self.elev = elev

    def __str__(self):
        ret = 'Station Object\n--------------\n'
        ret += 'sta:\t\t%s' % name
        ret += 'lat:\t\t%s' % lat
        ret += 'lon:\t\t%s' % lon
        ret += 'elev:\t\t%s' % elev
        return ret

class Event():
    """
    A container class for earthquake event metadata.
    """
    #def __init__(self, time, lat, lon, depth, mag, magtype=None, evid=None):
    def __init__(self,
                 prefor,
                 evid=None,
                 evname=None,
                 auth=None,
                 commid=None,
                 lddate=None,
                 origins=None):
        """
        Initialize Event object using one of two possible inputs.
        """
        import time as pytime
        self.evid = evid
        self.evname = evname
        self.prefor = prefor
        self.auth = auth
        self.commid = commid
        self.lddate = lddate
        self.preferred_origin = None
        if origins == None: self.origins = []
        else: self.origins = origins

    def __str__(self):
        ret = 'Event Object\n------------\n'
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'evname:\t\t%s\n' % self.evname
        ret += 'prefor:\t\t%s\n' % self.prefor
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'origins:\n'
        if len(self.origins) == 0:
            ret += '\t\tNone\n'
        else:
            for i in range(len(self.origins)):
                for line in  ('%s' % self.origins[i]).split('\n'):
                    ret += '\t\t%s\n' % line
        return ret

    def set_preferred_origin(self, prefor):
        """
        Set self.preferred_origin to equal origin with orid == prefor.
        """
        for i in range(len(self.origins)):
            if self.origins[i].orid == prefor:
                self.preferred_origin = self.origins[i]
                return 0

    def add_origin(self,
                   lat,
                   lon,
                   depth,
                   time,
                   auth,
                   arrivals=[],
                   orid=None,
                   evid=None,
                   jdate=None,
                   nass=None,
                   ndef=None,
                   ndp=None,
                   grn=None,
                   srn=None,
                   etype=None,
                   review=None,
                   depdp=None,
                   dtype=None,
                   mb=None,
                   mbid=None,
                   ms=None,
                   msid=None,
                   ml=None,
                   mlid=None,
                   algorithm=None,
                   commid=None,
                   lddate=None):
        """
        Add an Origin object to the list of origins associated with this event.
        """
        self.origins += [Origin(lat,
                                lon,
                                depth,
                                time,
                                auth,
                                orid=orid,
                                evid=evid,
                                arrivals=arrivals,
                                jdate=jdate,
                                nass=nass,
                                ndef=ndef,
                                ndp=ndp,
                                grn=grn,
                                srn=srn,
                                etype=etype,
                                review=review,
                                depdp=depdp,
                                dtype=dtype,
                                mb=mb,
                                mbid=mbid,
                                ms=ms,
                                msid=msid,
                                ml=ml,
                                mlid=mlid,
                                algorithm=algorithm,
                                commid=commid,
                                lddate=lddate)]
    def write(self, out, fmt):
        """
        Write out newly authored data pertaining to event.

        Arguments:
        out - A datascope db pointer to an open CSS3.0 database for output.
        Alternatively the path to an output SCEDC format flat file.

        fmt - The format of the 'out' argument; 'CSS3.0' or 'SCEDC'.

        Return Values:
        0 - Success
        -1 - Failure
        """
        if fmt == 'CSS3.0':
            return self._write_CSS()
        elif fmt == 'SCEDC':
            #Do it the SCEDC way
            pass
        else:
            raise Exception('Output format %s not recognized' % fmt)

class Origin():
    """
    A container class for origin data.
    """
    def __init__(self,
                 lat,
                 lon,
                 depth,
                 time,
                 auth,
                 arrivals=[],
                 orid=None,
                 evid=None,
                 jdate=None,
                 nass=None,
                 ndef=None,
                 ndp=None,
                 grn=None,
                 srn=None,
                 etype=None,
                 review=None,
                 depdp=None,
                 dtype=None,
                 mb=None,
                 mbid=None,
                 ms=None,
                 msid=None,
                 ml=None,
                 mlid=None,
                 algorithm=None,
                 commid=None,
                 lddate=None):
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.time = time
        self.orid = orid
        self.evid = evid
        self.auth = auth
        self.arrivals = arrivals
        self.jdate = jdate
        self.nass = nass
        self.ndef = ndef
        self.ndp = ndp
        self.grn = grn
        self.srn = srn
        self.etype = etype
        self.review = review
        self.depdp = depdp
        self.dtype = dtype
        self.mb = mb
        self.mbid = mbid
        self.ms = ms
        self.msid = msid
        self.ml = ml
        self.mlid = mlid
        self.algorithm = algorithm
        self.commid = commid
        self.lddate = lddate

    def __str__(self):
        """
        Return string representation of Origin object.
        """
        ret = 'Origin Object\n-------------\n'
        ret += 'lat:\t\t%s\n' % self.lat
        ret += 'lon:\t\t%s\n' % self.lon
        ret += 'depth:\t\t%s\n' % self.depth
        ret += 'time:\t\t%s\n' % self.time
        ret += 'orid:\t\t%s\n' % self.orid
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'jdate:\t\t%s\n' % self.jdate
        ret += 'nass:\t\t%s\n' % self.nass
        ret += 'ndef:\t\t%s\n' % self.ndef
        ret += 'ndp:\t\t%s\n' % self.ndp
        ret += 'grn:\t\t%s\n' % self.grn
        ret += 'srn:\t\t%s\n' % self.srn
        ret += 'etype:\t\t%s\n' % self.etype
        ret += 'review:\t\t%s\n' % self.review
        ret += 'depdp:\t\t%s\n' % self.depdp
        ret += 'dtype:\t\t%s\n' % self.dtype
        ret += 'mb:\t\t%s\n' % self.mb
        ret += 'mbid:\t\t%s\n' % self.mbid
        ret += 'ms:\t\t%s\n' % self.ms
        ret += 'msid:\t\t%s\n' % self.msid
        ret += 'ml:\t\t%s\n' % self.ml
        ret += 'mlid:\t\t%s\n' % self.mlid
        ret += 'algorithm:\t\t%s\n' % self.algorithm
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'arrivals:\n'
        for i in range(len(self.arrivals)):
            ret += '%s' % self.arrivals[i]
        return ret

class Phase():
    """
    A container class for phase data.
    """
    def __init__(self, sta, time, phase, chan=None, qual=None, arid=None):
        self.sta = sta
        self.time = time
        self.phase = phase
        self.chan = chan
        self.qual = qual
        self.arid = arid
        self.tt_calc=None #Calculated travel time

    def __str__(self):
        ret = 'Arrival Object\n--------------\n'
        ret += 'sta:\t\t%s\n' % self.sta
        ret += 'time:\t\t%s\n' % self.time
        ret += 'phase:\t\t%s\n' % self.phase
        ret += 'arid:\t\t%s\n' % self.arid
        ret += 'qual:\t\t%s\n'  % self.qual
        return ret

