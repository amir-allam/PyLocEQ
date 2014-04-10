program make_vgrids

! this program produces a file containing the velocity as a function of position that 
! can be used as input for the 3-D fast marching code. The velocities
! are given as values on a grid of nodes that is used as a basis 
! for cubic spline interpolation in the FM code. The resolution of the velocity
! grids is typically lower than on the grid used for the fast marching.
! If more velocity regions are present, they are assumed to be written to the file in 
! order from top to bottom, i.e. the first grid refers to the region just below the
! surface, the next one to the one below that etc. 


implicit none 
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

real(kind=dp)  :: rmin, latmin,longmin

real(kind=dp)  :: stretch = 1.01_dp
REAL(KIND=dp)  :: earth_radius = 6371.0_dp
integer        :: nv0 = 5
integer        :: n_vtypes = 2
integer        :: size_ratio = 4
integer        :: nr,nlat,nlong 
integer        :: n,m,i,j,k,nv

real(kind=dp)  :: dr,dlat,dlong,r,lat,long,deg_to_rad,vel,latc,longc,rc,vfac
real(kind=dp)  :: prmin,platmin,plongmin,pdr,pdlat,pdlong
real(kind=dp)  :: rak(136),vpak(136),vsak(136),dak(136)

integer        :: pnr,pnlat,pnlong,n_interfaces


open(1,file='propgrid.in')

read(1,*) pnr,pnlat,pnlong
read(1,*) pdr,pdlat,pdlong
read(1,*) prmin,platmin,plongmin

close(1)

open(1,file='interfaces.in')
read(1,*) n_interfaces
close(1)

nv = n_interfaces - 1

deg_to_rad=acos(-1.0_dp)/180.0_dp

pdlat=pdlat*deg_to_rad
pdlong=pdlong*deg_to_rad
platmin=platmin*deg_to_rad
plongmin=plongmin*deg_to_rad

prmin =  earth_radius + prmin - dble(pnr-1)*pdr


! # of points in velocity grid

! first extend grid to an integer number times the size ration

i = mod(pnr-1,size_ratio)
if (i>0) pnr = pnr + (size_ratio-i)
i = mod(pnlat-1,size_ratio)
if (i>0) pnlat = pnlat + (size_ratio-i)
i = mod(pnlong-1,size_ratio)
if (i>0) pnlong = pnlong + (size_ratio-i)


nr    = (pnr-1)/size_ratio + 3  !size_ratio-1
nlat  = (pnlat-1)/size_ratio + 3 !size_ratio-1
nlong = (pnlong-1)/size_ratio + 3  !size_ratio-1


! interval size on velocity grid 

dr     =stretch*pdr*size_ratio
dlat   =stretch*pdlat*size_ratio
dlong  =stretch*pdlong*size_ratio

! origin of velocity grid 

rmin     =  prmin - dr - (nr-1)*dr*(stretch-1.0)/2
latmin   =  platmin - dlat - (nlat-1)*dlat*(stretch-1.0)/2
longmin  =  plongmin -dlong - (nlong-1)*dlong*(stretch-1.0)/2


call random_seed

open(11,file='vgrids.in')

write (11,*) nv,n_vtypes
do m=1,n_vtypes
if (m==1) vfac= 1.0_dp
if (m==2) vfac= 0.6_dp 
do n=1,nv

   write (11,*) nr,nlat,nlong
   write (11,*) dr,dlat,dlong
   write (11,*) rmin,latmin,longmin

   do i=1,nr
      r=rmin+(i-1)*dr
      do j=1,nlat
         lat=latmin+(j-1)*dlat
         do k=1,nlong
            long = longmin+(k-1)*dlong
            select case(n)
               case(1)
                  vel=vel1(r,lat,long)*vfac
                  write (11,*) vel
               case(2)
                  vel=vel2(r,lat,long)*vfac
                  write (11,*) vel
               case(3)
                  vel=vel3(r,lat,long)*vfac
                  write (11,*) vel
               case(4)
                  vel=vel4(r,lat,long)*vfac
                  write (11,*) vel
               case(5)
                  vel=vel5(r,lat,long)*vfac
                  write (11,*) vel
            end select
         end do
      end do
   end do

end do 

end do  ! vtypes

close(11)

contains

  function vel1(rev,latev,longev)
    real(kind=dp)  :: vel1,rev,latev,longev
    vel1=2.0_dp
    return
  end function vel1

  function vel2(rev,latev,longev)
    real(kind=dp)  :: vel2,rev,latev,longev
    vel2=3.0_dp
    return
  end function vel2

  function vel3(rev,latev,longev)
    real(kind=dp)  :: vel3,rev,z,latev,longev,rad
    z=earth_radius-rev
    call random_number(rad)
    vel3=4.0_dp + (z-5.0_dp)*0.05_dp + (rad-0.5_dp)*0.0_dp
    return
  end function vel3

  function vel4(rev,latev,longev)
    real(kind=dp)  :: vel4,rev,z,latev,longev
    z=earth_radius-rev
    vel4=5.0_dp + (z-5.0_dp)*0.05_dp
    return
  end function vel4

  function vel5(rev,latev,longev)
    real(kind=dp)  :: vel5,rev,z,latev,longev
    vel5 = 8.0_dp
    return
  end function vel5

end program make_vgrids



