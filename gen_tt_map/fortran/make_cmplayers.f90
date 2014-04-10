program make_interfaces

! this program produces a file containing the position of the interfaces that 
! can be used as input for the 3-D fast marching code. The positions of the
! interfaces are given as values on a grid of nodes that is used as a basis 
! for cubic spline interpolation in the FM code. The resolution of the interface
! grids is typically lower than on the grid used for the fast marching.
! Interfaces are assumed to be written the file in order from top to bottom


implicit none 
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

integer        :: nv0 = 6             ! # of interfaces
integer        :: size_ratio = 4     ! ratio of velocity grid interval size 
                                     ! to propagation grid interval size
integer        :: nr,nlat,nlong      ! # of grid points of velocity grid
integer        :: n,i,j,k

real(kind=dp)  :: latmin,longmin,deg_to_rad,dr,dlat,dlong,r,lat,long
real(kind=dp)  :: prmin,platmin,plongmin,pdr,pdlat,pdlong,twopi,h,wlat,wlong
integer        :: pnr,pnlat,pnlong,nv
real(kind=dp)  :: stretch = 1.01_dp
REAL(KIND=dp)      :: earth_radius = 6371.0_dp
! some useful constants

deg_to_rad=acos(-1.0_dp)/180.0_dp
twopi=2.0_dp*acos(-1.0_dp)



! read properties of the FM propagation grid

open(1,file='propgrid.in')

read(1,*) pnr,pnlat,pnlong
read(1,*) pdr,pdlat,pdlong
read(1,*) prmin,platmin,plongmin

close(1)

! convert from degrees to radians and depth to radius

prmin =  earth_radius + prmin - dble(pnr-1)*pdr
pdlat=pdlat*deg_to_rad
pdlong=pdlong*deg_to_rad
platmin=platmin*deg_to_rad
plongmin=plongmin*deg_to_rad

! wavelength of interface perturbations (a specific example)
wlat=(pnlat-1)*pdlat*4.0
wlong=(pnlong-1)*pdlong*4.0


! # of points in interface grid

i = mod(pnlat-1,size_ratio)
if (i>0) pnlat = pnlat + (size_ratio-i)
i = mod(pnlong-1,size_ratio)
if (i>0) pnlong = pnlong + (size_ratio-i)


nlat  = (pnlat-1)/size_ratio + 3 
nlong = (pnlong-1)/size_ratio + 3


! interval size on interface grid. Grid is slightly stretched to make sure the propagation
! grid falls completely inside it

dlat   =stretch*pdlat*size_ratio
dlong  =stretch*pdlong*size_ratio


! origin of interface grid

latmin   =  platmin - dlat - (nlat-1)*dlat*(stretch-1.0)/2
longmin  =  plongmin -dlong - (nlong-1)*dlong*(stretch-1.0)/2


! write the grids to a file taken as input by the FM code

nv=nv0

open(11,file='interfaces.in')

write (11,*) nv

! if making changes here, remember that interfaces have to be written in order
! from top to bottom

   write (11,*) nlat,nlong
   write (11,*) dlat,dlong
   write (11,*) latmin,longmin


do n=1,nv

   do j=1,nlat

      lat=latmin+(j-1)*dlat

      do k=1,nlong

         long=longmin+(k-1)*dlong

         if (n == nv) then
            write (11,*) rbot(lat,long)
         else            

            select case(n)
            case(1)
               write (11,*) rsurf(lat,long)
            case(2)
               write (11,*) ints(lat,long)
            case(3)
               write (11,*) int1(lat,long)
            case(4)
               write (11,*) int2(lat,long)
            case(5)
               write (11,*) int3(lat,long)
            case(6)
               write (11,*) rbot(lat,long)
            end select

         endif

      end do

   end do
 
end do 

close(11)


! below are functions returning the !! RADIUS !! of an interface as a function of latitude
! and longitude

contains

  function rsurf(latev,longev)
    real(kind=dp)  :: rsurf,latev,longev
    rsurf=earth_radius
    return
  end function rsurf

  function rbot(latev,longev)
    real(kind=dp)  :: rbot,latev,longev
    rbot=prmin
    return
  end function rbot

  function ints(latev,longev)
    real(kind=dp)  :: ints,latev,longev,x,y,z,a,b,h1,h2
    x=(latev-latmin)/((nlat-1)*dlat)
    y=(longev-longmin)/((nlong-1)*dlong)
    z=(pnr-1)*pdr
    a=(pnr-1)*pdr/15.
    b=z/10.
    h1=prmin + z*0.95_dp
    h2=prmin + z*0.85_dp+a*(1.0-cos(twopi*latev/(wlat))*cos(twopi*longev/(wlong)))+(b*(x+y))
    ints=max(h1,h2)
    return
  end function ints

  function int1(latev,longev)
    real(kind=dp)  :: int1,latev,longev,x,y,z,a,b,h1,h2
    x=(latev-latmin)/((nlat-1)*dlat)
    y=(longev-longmin)/((nlong-1)*dlong)
    z=(pnr-1)*pdr
    a=(pnr-1)*pdr/10.
    b=z/10.
    h1=prmin + z*0.80_dp+a*(1.0 - cos(twopi*latev/(wlat))*cos(twopi*longev/(wlong)))+(b*(x+y))
    h2=prmin + z*0.90_dp
    int1=max(h1,h2)
    return
  end function int1

  function int2(latev,longev)
    real(kind=dp)  :: int2,latev,longev,x,y,z,a,b
    x=(latev-latmin)/((nlat-1)*dlat)
    y=(longev-longmin)/((nlong-1)*dlong)
    z=(pnr-1)*pdr
    a=z/3.
    b=z/5.
    int2=prmin + z*0.92_dp-a*(1.0-cos(twopi*latev/wlat)*cos(twopi*longev/wlong))-b*(x**2+y**2)
    return
  end function int2

  function int3(latev,longev)
    real(kind=dp)  :: int3,latev,longev,x,y,z,a,b
    x=(latev-latmin)/((nlat-1)*dlat)
    y=(longev-longmin)/((nlong-1)*dlong)
    z=(pnr-1)*pdr
    a=z/3.
    b=z/5.
    int3=prmin + z*0.72_dp-a*(1.0-cos(twopi*latev/wlat)*cos(twopi*longev/wlong))-b*(x**2+y**2)
    return
  end function int3

end program make_interfaces



