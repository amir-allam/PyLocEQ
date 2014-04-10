program make_vgrids

! this program produces a file containing the velocity as a function of position that 
! can be used as input for the 3-D fast marching code. The velocities
! are given as values on a grid of nodes that is used as a basis 
! for cubic spline interpolation in the FM code. The resolution of the velocity
! grids is typically lower than on the grid used for the fast marching.
! If more velocity regions are present, they are assumed to be written to the file in 
! order from top to bottom, i.e. the first grid refers to the region just below the
! surface, the next one to the one below that etc. 
!
! This program generates P velocities specific to the ak135 reference model. 
! Discontinuities inside the core are 
! not treated as explicit discontinuities. Only use in conjunction with the
! program make_aklayers that defines the position of the velocity discontinuities 
! treated explicitly


implicit none 
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

real(kind=dp)  :: rmin, latmin,longmin

real(kind=dp)  :: stretch = 1.01_dp
REAL(KIND=dp)      :: earth_radius = 6371.0_dp
integer        :: nv0 = 6
integer        :: n_vtypes = 1
integer        :: size_ratio = 2  ! the ratio of the grid spacing on the velocity grid
                                  ! and the propagation grid
integer        :: nr,nlat,nlong 
integer        :: n,m,i,j,k,nv

real(kind=dp)  :: dr,dlat,dlong,r,lat,long,deg_to_rad,vel
real(kind=dp)  :: prmin,platmin,plongmin,pdr,pdlat,pdlong
real(kind=dp),target  :: rak(136),vpak(136),vsak(136),dak(136)
real(kind=dp),dimension(:),pointer :: vak
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

open(1,file='ak135.dat')
do i=1,136
   read(1,*) rak(i),vpak(i),vsak(i),dak(i)
end do
rak=earth_radius - rak
close(1)

deg_to_rad=acos(-1.0_dp)/180.0_dp

pdlat=pdlat*deg_to_rad
pdlong=pdlong*deg_to_rad
platmin=platmin*deg_to_rad
plongmin=plongmin*deg_to_rad

prmin =  earth_radius + prmin - dble(pnr-1)*pdr


! # of points in velocity grid
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

open(11,file='vgrids.in')

write (11,*) nv, n_vtypes
do m=1,n_vtypes

   if (m==1) vak => vpak
   if (m==2) vak => vsak

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
                  vel=vel1(r)
                  write (11,*) vel
               case(2)
                  vel=vel2(r)
                  write (11,*) vel
               case(3)
                  vel=vel3(r)
                  write (11,*) vel
               case(4)
                  vel=vel4(r)
                  write (11,*) vel
               case(5)
                  vel=vel5(r)
                  write (11,*) vel
               case(6)
                  vel=vel6(r)
                  write (11,*) vel
               end select
            end do
         end do
      end do

   end do

end do

close(11)

contains

  function vel1(rev)
    real(kind=dp)  :: vel1,rev
    vel1=vak(1)
    return
  end function vel1

  function vel2(rev)
    real(kind=dp)  :: vel2,rev
    vel2=vak(3)
    return
  end function vel2

  function vel3(rev)
    real(kind=dp)  :: vel3,rev
    integer :: ii
    if (rev <= rak(14)) then
       vel3=vak(14)
       return
    endif
    if (rev >= rak(5)) then
       vel3=vak(5)
       return
    endif
    do ii=5,14
       if (rev >= rak(ii+1)) then
          vel3 = vak(ii)+(vak(ii+1)-vak(ii))*(rev-rak(ii))/(rak(ii+1)-rak(ii))
          return
       endif
    end do
    stop 'vel3 error'
  end function vel3

  function vel4(rev)
    real(kind=dp)  :: vel4,rev
    integer :: ii
    if (rev <= rak(20)) then
       vel4=vak(20)
       return
    endif
    if (rev >= rak(15)) then
       vel4=vak(15)
       return
    endif
    do ii=15,20
       if (rev >= rak(ii+1)) then
          vel4 = vak(ii)+(vak(ii+1)-vak(ii))*(rev-rak(ii))/(rak(ii+1)-rak(ii))
          return
       endif
    end do
    stop 'vel4 error'
  end function vel4

  function vel5(rev)
    real(kind=dp)  :: vel5,rev
    integer :: ii
    if (rev <= rak(67)) then
       vel5=vak(67)
       return
    endif
    if (rev >= rak(21)) then
       vel5=vak(21)
       return
    endif
    do ii=21,67
       if (rev >= rak(ii+1)) then
          vel5= vak(ii)+(vak(ii+1)-vak(ii))*(rev-rak(ii))/(rak(ii+1)-rak(ii))
          return
       endif
    end do
    stop 'vel5 error'
  end function vel5

  function vel6(rev)
    real(kind=dp)  :: vel6,rev
    integer :: ii
    if (rev >= rak(68)) then
       vel6=vak(68)
       return
    endif

    do ii=68,135
       if (rev >= rak(ii+1)) then
          vel6 = vak(ii)+(vak(ii+1)-vak(ii))*(rev-rak(ii))/(rak(ii+1)-rak(ii))
          return
       endif
    end do
    stop 'vel6 error'
  end function vel6


end program make_vgrids



