program make_vgrids

! this program produces a file containing the velocity as a function of position that 
! can be used as input for the 3-D fast marching code. The velocities
! are given as values on a grid of nodes that is used as a basis 
! for cubic spline interpolation in the FM code. The resolution of the velocity
! grids is typically lower than on the grid used for the fast marching.
! If more velocity regions are present, they are assumed to be written to the file in 
! order from top to bottom, i.e. the first grid refers to the region just below the
! surface, the next one to the one below that etc. 

! In general the user must provide a set of functions vel1(dep,lat,long,vtype),
! vel2(dep,lat,long,vtype), vel3()... that return the velocity as a function of position
! in each region (vel1 returns velocities in region 1, vel2 in region 2 etc.)
! vtype refers to the type of velocity, typically vtype=1 corresponds to P velocity
! and vtype=2 to S-velocity. Note that the units of the position arguments of the 
! vel1 etc functions are km, degrees and degrees.
! It is somwhat inelegant but 10 velocity functions vel1-vel10 have to provided even
! if only a smaller number is used. The unused ones will be simply ignored but
! have to be present. 

! Provide these functions in a separate file, then compile them together with this
! program. Examples are provided.
! e.g. "f90 -o mvg make_vgrid.f90 vdefs_ak135.f90" creates the executable mvg that when
! executed creates a file vgrids.in corresponding to the ak136 model that can be used 
! as input for the fast marching code




implicit none 

real(kind=8)  :: rmin, latmin,longmin

real(kind=8)  :: stretch = 1.01
REAL(kind=8)  :: earth_radius = 6371.0
integer        :: n_vtypes = 2
integer        :: size_ratio = 2
integer        :: nr,nlat,nlong 
integer        :: n,m,i,j,k,nv,vtype

real(kind=8)  :: dr,dlat,dlong,r,lat,long,deg_to_rad,vel,latc,latdeg,longdeg
real(kind=8)  :: prmin,platmin,plongmin,pdr,pdlat,pdlong,depth
real(kind=8)  :: rak(136),vpak(136),vsak(136),dak(136)

integer        :: pnr,pnlat,pnlong,n_interfaces

real(kind=8)  :: vel1,vel2,vel3,vel4,vel5,vel6,vel7,vel8,vel9,vel10

! read in the parameters of the propagation grid
open(1,file='propgrid.in')
read(1,*) pnr,pnlat,pnlong
read(1,*) pdr,pdlat,pdlong
read(1,*) prmin,platmin,plongmin
close(1)

! convert lat/long from degrees to radians
deg_to_rad=acos(-1.0d0)/180.0d0
pdlat=pdlat*deg_to_rad
pdlong=pdlong*deg_to_rad
platmin=platmin*deg_to_rad
plongmin=plongmin*deg_to_rad

! convert origin of the propagation grid from depth to radius 
prmin =  earth_radius + prmin - dble(pnr-1)*pdr

! read in the number of interfaces to get the number
! of velocity grids required
open(1,file='interfaces.in')
read(1,*) n_interfaces
close(1)
nv = n_interfaces - 1

! # of nodes required in velocity grid in each direction
i = mod(pnr-1,size_ratio)
if (i>0) pnr = pnr + (size_ratio-i)
i = mod(pnlat-1,size_ratio)
if (i>0) pnlat = pnlat + (size_ratio-i)
i = mod(pnlong-1,size_ratio)
if (i>0) pnlong = pnlong + (size_ratio-i)


nr    = (pnr-1)/size_ratio + 3 
nlat  = (pnlat-1)/size_ratio + 3
nlong = (pnlong-1)/size_ratio + 3


! interval size on velocity grid 
! slightly stretched to ensure the propagation grid
! is entirely within the second layer of velocity nodes
dr     =stretch*pdr*size_ratio
dlat   =stretch*pdlat*size_ratio
dlong  =stretch*pdlong*size_ratio

! origin of velocity grid 
rmin     =  prmin - dr - (nr-1)*dr*(stretch-1.0)/2
latmin   =  platmin - dlat - (nlat-1)*dlat*(stretch-1.0)/2
longmin  =  plongmin -dlong - (nlong-1)*dlong*(stretch-1.0)/2

! initialize random number generator in case it is to be used
! by the velocity functions
call random_seed

! open the file that will be the input file for the fast marching
! code
open(11,file='vgrids.in')

! nr. of velocity grids, nr, of velocity types
write (11,*) nv, n_vtypes

do vtype=1,n_vtypes

   do n=1,nv
   ! the parameters of this velocity grid
      write (11,*) nr,nlat,nlong
      write (11,*) dr,dlat,dlong
      write (11,*) rmin,latmin,longmin

   ! the actual values at each of the grid nodes
      do i=1,nr
         r=rmin+(i-1)*dr
         depth=earth_radius-r
         do j=1,nlat
            lat=latmin+(j-1)*dlat
            latdeg=lat/deg_to_rad
            do k=1,nlong
               long = longmin+(k-1)*dlong
               longdeg=long/deg_to_rad
               select case(n)
                  case(1)
                     vel=vel1(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(2)
                     vel=vel2(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(3)
                     vel=vel3(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(4)
                     vel=vel4(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(5)
                     vel=vel5(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(6)
                     vel=vel6(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(7)
                     vel=vel7(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(8)
                     vel=vel8(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(9)
                     vel=vel9(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case(10)
                     vel=vel10(depth,latdeg,longdeg,vtype)
                     write (11,*) vel
                  case default
                     stop 'error:modify program if more than 10 regions required'
              end select

           end do  !long loop
        end do  ! lat loop
     end do  ! r loop

  end do  ! regions loop

end do ! vtypes

close(11)

end program make_vgrids



