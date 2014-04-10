program make_interfaces

! this program produces a file containing the position of the interfaces that 
! can be used as input for the 3-D fast marching code. The positions of the
! interfaces are given as values on a grid of nodes that is used as a basis 
! for cubic spline interpolation in the FM code. The resolution of the interface
! grids is typically lower than on the grid used for the fast marching.
! Interfaces are assumed to be written the file in order from top to bottom


implicit none 

integer        :: nif             ! # of interfaces
integer        :: size_ratio = 2     ! ratio of velocity grid interval size 
                                     ! to propagation grid interval size
integer        :: nr,nlat,nlong      ! # of grid points of velocity grid
integer        :: n,i,j,k

real(kind=8)  :: latmin,longmin,deg_to_rad,dr,dlat,dlong,r,lat,long,latdeg,longdeg
real(kind=8)  :: prmin,platmin,plongmin,pdr,pdlat,pdlong,h
integer        :: pnr,pnlat,pnlong,nv
real(kind=8)  :: stretch = 1.01
REAL(kind=8)      :: earth_radius = 6371.0
real(kind=8)  :: iface1,iface2,iface3,iface4,iface5,iface6
real(kind=8)  :: iface7,iface8,iface9,iface10,iface11

! some useful constants

deg_to_rad=acos(-1.0)/180.0

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


! get the number of interfaces

call get_n_interfaces(nif) 


! write the grids to a file taken as input by the FM code

open(11,file='interfaces.in')

write (11,*) nif

! if making changes here, remember that interfaces have to be written in order
! from top to bottom

   write (11,*) nlat,nlong
   write (11,*) dlat,dlong
   write (11,*) latmin,longmin

do n=1,nif

   do j=1,nlat
      lat=latmin+(j-1)*dlat
      latdeg=lat/deg_to_rad
      do k=1,nlong
         long=longmin+(k-1)*dlong
         longdeg=long/deg_to_rad

            select case(n)

            case(1)
               write (11,*) earth_radius - iface1(latdeg,longdeg)
            case(2)
               write (11,*) earth_radius - iface2(latdeg,longdeg)
            case(3)
               write (11,*) earth_radius - iface3(latdeg,longdeg)
            case(4)
               write (11,*) earth_radius - iface4(latdeg,longdeg)
            case(5)
               write (11,*) earth_radius - iface5(latdeg,longdeg)
            case(6)
               write (11,*) earth_radius - iface6(latdeg,longdeg)
            case(7)
               write (11,*) earth_radius - iface7(latdeg,longdeg)
            case(8)
               write (11,*) earth_radius - iface8(latdeg,longdeg)
            case(9)
               write (11,*) earth_radius - iface9(latdeg,longdeg)
            case(10)
               write (11,*) earth_radius - iface10(latdeg,longdeg)
            case(11)
               write (11,*) earth_radius - iface11(latdeg,longdeg)
            case default
               stop 'edit program if more than 11 interfaces required'

            end select

      end do

   end do
 
end do 

close(11)

end program make_interfaces



