program make_interfaces

! this program produces a file containing the position of the interfaces that 
! can be used as input for the 3-D fast marching code. The positions of the
! interfaces are given as values on a grid of nodes that is used as a basis 
! for cubic spline interpolation in the FM code. The resolution of the interface
! grids is typically lower than on the grid used for the fast marching.
! Interfaces are assumed to be written the file in order from top to bottom


implicit none 
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

integer        :: nv0 = 7             ! # of interfaces
integer        :: size_ratio = 4     ! ratio of velocity grid interval size 
                                     ! to propagation grid interval size
integer        :: nr,nlat,nlong      ! # of grid points of velocity grid
integer        :: n,i,j,k

real(kind=dp)  :: latmin,longmin,deg_to_rad,dr,dlat,dlong,r,lat,long
real(kind=dp)  :: prmin,platmin,plongmin,pdr,pdlat,pdlong,twopi,h
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


! convert from degrees to radians

pdlat=pdlat*deg_to_rad
pdlong=pdlong*deg_to_rad
platmin=platmin*deg_to_rad
plongmin=plongmin*deg_to_rad

prmin =  earth_radius + prmin - dble(pnr-1)*pdr

! # of points in interface grid

i = mod(pnlat-1,size_ratio)
if (i>0) pnlat = pnlat + (size_ratio-i)
i = mod(pnlong-1,size_ratio)
if (i>0) pnlong = pnlong + (size_ratio-i)


nlat  = (pnlat-1)/size_ratio + 3   !size_ratio-1
nlong = (pnlong-1)/size_ratio + 3  !size_ratio-1


! interval size on interface grid 

dlat   =stretch*pdlat*size_ratio
dlong  =stretch*pdlong*size_ratio


! origin of interface grid

latmin   =  platmin - dlat - (nlat-1)*dlat*(stretch-1.0)/2
longmin  =  plongmin -dlong - (nlong-1)*dlong*(stretch-1.0)/2


if (prmin < rl6(0._dp,0._dp)) nv=7
if (prmin >= rl6(0._dp,0._dp)) nv=6
if (prmin >= rl5(0._dp,0._dp)) nv=5
if (prmin >= rl4(0._dp,0._dp)) nv=4
if (prmin >= rl3(0._dp,0._dp)) nv=3
if (prmin >= rl2(0._dp,0._dp)) nv=2


! write the grids to a file taken as input by the FM code

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
               write (11,*) rl2(lat,long)
            case(3)
               write (11,*) rl3(lat,long)
            case(4)
               write (11,*) rl4(lat,long)
            case(5)
               write (11,*) rl5(lat,long)
            case(6)
               write (11,*) rl6(lat,long)
            case(7)
               write (11,*) rbot(lat,long)
            end select

         endif

      end do

   end do
 
end do 

close(11)


! below are functions returning the radius of an interface as a function of latitude
! and longitude

contains

  function rsurf(latev,longev)
    real(kind=dp)  :: rsurf,latev,longev
    rsurf=6371.0_dp
    return
  end function rsurf

  function rbot(latev,longev)
    real(kind=dp)  :: rbot,latev,longev
    rbot=prmin
    return
  end function rbot

  function rl2(latev,longev)
    real(kind=dp)  :: rl2,latev,longev
    rl2=earth_radius - 20.0_dp
    return
  end function rl2

  function rl3(latev,longev)
    real(kind=dp)  :: rl3,latev,longev
    rl3=earth_radius - 35.0_dp
    return
  end function rl3

  function rl4(latev,longev)
    real(kind=dp)  :: rl4,latev,longev
    rl4=earth_radius - 410.0_dp
    return
  end function rl4

  function rl5(latev,longev)
    real(kind=dp)  :: rl5,latev,longev
    rl5=earth_radius - 660.0_dp
    return
  end function rl5


  function rl6(latev,longev)
    real(kind=dp)  :: rl6,latev,longev
    rl6=earth_radius - 2891.5_dp
    return
  end function rl6



  function rmoho(latev,longev)
    real(kind=dp)  :: rmoho,latev,longev
    rmoho=6300.0_dp
    return
  end function rmoho

end program make_interfaces



