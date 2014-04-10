3dfm_main.f90
!************************************************************************************************
! This is the main program unit for the 3d-multistage FMM code
! The purpose of this unit is to 
! - read the input files,
! - perform a range of initializations (including the wave front propagation
!   on the refined source grid) by calling routines located in in the file 3dfmlib.f90
! - go through the multistage FMM sequences on the main propagation grid as 
!   defined in the input files by initializing the appropriate intersections 
!   and calling the fast marching routines located in the file propagate.f90
! - call the ray tracing routines for all paths by calling routines from the file rays.f90
! - finish by calling the Frechet derivatives routines in frechet.f90
! - along the way, create the OpenDx visualisation files by calling routines from visual.f90

program fm3d
  use mod_3dfm
  implicit none

  integer                            :: n,i,j,k,m,i1,i2,ns,path_id,nsec,top_id,bot_id,vtype
  integer                            :: ctype,step,prev_tf,n_tf,recpath_id,direction
  integer                            :: t1,t2,t3,t4,t5,count_rate,count_max,io
  integer                            :: is_teleseismic,n_teleseismic,nfile,nsave,save_grid_counter
!  integer                            :: i3,n_tot_ref
  type(Tintersection),pointer        :: istart,inext
  type(Tregion),pointer              :: reg
  type(Tsource),pointer              :: s,ss
  type(Treceiver),pointer            :: rec
  type(Tray),pointer                 :: ray
!  type(Ttime_field),pointer          :: tf1,tf2,tf3
  real(kind=dp)                      :: xr(3),xs(3),xcp(3)
  real(kind=dp)                      :: deg_to_rad,r
!  real(kind=dp)                      :: dist,hb,h,p_arrival,mean_error
  real(kind=dp)                      :: interpolate_interface
  real(kind=dp)                      :: t_arrival
  integer,dimension(:),allocatable   :: npp_receiver,ipp_receiver,iarr,jarr
  logical                            :: do_frechet_derivatives
  real                               :: evlat,evlon,ndlat,ndlon,deltas,cazim,bazim,azima
  integer                            :: myid,nproc


!  logical,dimension(:),allocatable   :: logarr

! read the modes under which the program is to run from file

  open(1,file='mode_set.in')

  read (1,*) file_mode
  read (1,*) no_pp_mode
  read (1,*) parallel_mode
  read (1,*) display_mode
  read (1,*) save_rays_mode
  read (1,*) save_timefields_mode

  close(1)


  if (no_pp_mode .and. file_mode) &
       stop 'if no_pp_mode enabled file_mode must be disabled'

  if (.not.no_pp_mode .and. parallel_mode) &
       stop 'currently parrallel_mode is only possible in no_pp_mode'

  if (parallel_mode) then
     open(1,file='myid')
     read(1,*) nproc,myid
     close(1)
  endif

  call system_clock(t1,count_rate,count_max)

  deg_to_rad=acos(-1.0_dp)/180._dp

!-------------------
! initialize objects that are defined in input files

  write(unit=*,fmt='(a30)',advance='no') 'initializing propagation grid'
  call initialize_propagation_grid
  print *,'......finished'

  write(unit=*,fmt='(a30)',advance='no')' initializing velocity grids '
  call initialize_velocity_grids
  print *,'......finished'

  write(unit=*,fmt='(a30)',advance='no')' initializing interfaces'
  call initialize_interfaces
  print *,'......finished'


!---------------------------------------------------------------------------------------------------
! initialize intersections (interface grid + navigation pointers) and regions (collection of grid
! points between two interfaces, including the bounding intersections ). The fast marching
! solution is evaluated region by region.  

!---------
!  calculate the intersection points, i.e. where the surfaces cut through the regular propagation grid

  n_intersections=n_interfaces
  allocate(intersection(n_intersections))

  do n=1,n_intersections ; call intersection_defaults(intersection(n)) ; intersection(n)%id=n ; end do

  do n=1, n_intersections
     call find_intersection(intersection(n),intrface(n),pgrid)
  end do

  print *,'main grid intersections initialized'


!----------
! set up the 3-D regions that lie between the interfaces 

  n_regions=n_interfaces-1
  allocate(region(n_regions))

  do n=1,n_regions

     call region_defaults(region(n))
     region(n)%id=n

! pointers to the intersections that are the boundaries of the region
     region(n)%itop => intersection(n)
     region(n)%ibot => intersection(n+1)

     if (region(n)%itop%nnode == 0) then
        print *,'region',n,'does not exist'
        region(n)%nnode = 0
        cycle
     endif

     call define_region(region(n),region(n)%itop,region(n)%ibot,pgrid)

  end do

  print *,'main grid regions initialized'


!---------

! set a flag at each node of the main grid that has only regular connected cells
! This is used later to speed up calculations. Intersections have to be initialized to do this

  call tag_regular_nodes(pgrid)


!---------
! calculate the velocities on the regular grid and intersection grid points
! and transfer them to the corresponding regional nodes

  call velocities_to_grid(pgrid)

  do n=1,n_intersections
     call velocities_to_intersection(intersection(n))
  end do

  do n=1,n_regions
     call velocities_to_region(region(n),pgrid)   
  end do

  print *,'velocities on main grid, its intersections and regions evaluated'


!------------------
! read source location from file, and determine some of its properties

  open(1,file='sources.in')

  read (1,*) n_sources             ! the number of sources

!  print *,'nsources is', n_sources

 ! we already need to know the number of receivers since these can potentially become sources
 ! if a phase containing late reflections (pp type) is requested

  open(2,file='receivers.in')

  read(2,*) n_receivers            ! the number of receivers

!  print *,'nreceivers is', n_receivers

  allocate(receiver(n_receivers))
  do n=1,n_receivers ; call receiver_defaults(receiver(n)) ; receiver(n)%id=n ; end do

  allocate(source(n_sources+n_receivers))

  do n=1,n_sources+n_receivers
     call source_defaults(source(n)) 
     source(n)%id=n 
  end do

  n_sources_ppinc=n_sources   ! will contain the number of sources including receivers of pp phases

  n_teleseismic = 0           ! will contain the number of teleseismic sources

  do ns=1,n_sources

     print *,'reading source',ns

     s => source(ns) 

     read (1,*) is_teleseismic

     if (is_teleseismic == 1) then
        
        s%is_teleseismic = .true.
        n_teleseismic = n_teleseismic+1

        read (1,*) s%teleseismic_phase
        read (1,*) s%r,s%lat,s%long

        s%r = earth_radius - s%r
        s%lat=deg_to_rad*s%lat
        s%long=deg_to_rad*s%long
        s%coslat=cos(s%lat)

        s%teleseismic_id = n_teleseismic

     else

        s%is_local = .true.

        read (1,*) s%r,s%lat,s%long

        s%r = earth_radius - s%r
        s%lat=deg_to_rad*s%lat
        s%long=deg_to_rad*s%long
        s%coslat=cos(s%lat)

        call initialize_source(s,pgrid)! determines where the source is located wrt grid/intersections


        if (.not. s%on_interface) print *,'source',ns,' is located in region ',s%region_id
        if (s%on_interface)  print *,'source',ns,' is located on interfaces ',s%topint_id,s%botint_id

     endif


!---------------------------------------------------------
! read in the required paths (reflections/transmissions)

     read(1,*) s%n_paths    ! the number of paths to this source

!     print *,s%n_paths,'is npath'

 ! allocate the array containing path information for this source

     allocate(s%path(s%n_paths))
     do n=1,s%n_paths ; call path_defaults(s%path(n)) ; end do

  ! read the paths

     do n=1,s%n_paths

!        print *,'reading path',n

 ! register the path id

        s%path(n)%id  = n 


 ! read the number of timefields on this path (= number of steps in the sequence)
        read(1,*) n_tf

        s%path(n)%n_tf = n_tf

!        print *,'ntf',s%path(n)%n_tf

        allocate(s%path(n)%sequence(2*n_tf))  ! contains the sequence definition
        allocate(s%path(n)%tf_sequence(n_tf)) ! indices of the time fields corresponding to each step
        allocate(s%path(n)%vtype_sequence(n_tf)) ! velocity type of each step

        read (1,*) s%path(n)%sequence(1:2*n_tf)
        read (1,*) s%path(n)%vtype_sequence(1:n_tf)

        if (count(s%path(n)%vtype_sequence(1:n_tf) == 2) > 0 .and. n_vtypes == 1) then
           print *
           print *,'***** ERROR  ***********************************'
           print *, 'two velocity types specified in path but only one defined'
           stop
        endif

!        print *,s%path(n)%sequence(1:2*s%path(n)%n_tf)

     ! check whether the path contains a reflection fit step. Store the number of the step if so

        do step=1,n_tf

           if (s%path(n)%sequence(2*step-1) == s%path(n)%sequence(2*step)) then
              s%path(n)%refstep=step
              s%path(n)%fitting_interface=s%path(n)%sequence(2*step)
              print *,'path',n,'step',step,' is a reflection step'
              if (no_pp_mode) then
                 print *,'when no_pp_mode is enabled late reflections are not allowed'
                 stop 'illegal request for late reflection fit while no_pp_mode enabled'
              endif
           endif

        end do

     end do

!     print *,s%n_paths,' paths read from file'


! first test paths for consistency with source position 

     if (s%is_local) then

        do n=1,s%n_paths
           if (s%path(n)%sequence(1) /= 0) then
              print *
              print *,'***** ERROR  ***********************************'
              print *,'first item in path definition must be 0 for local sources'
              print *,'source',s%id, ' path',n
              stop
           endif
        end do

        if (s%on_interface) then

           do n=1,s%n_paths

              if (.not.((s%path(n)%sequence(2) == s%botint_id+1  &
                   .or. s%path(n)%sequence(2) == s%topint_id-1))) then
                 print *
                 print *,'***** ERROR  ***********************************'
                 print *,'path',s%path(n)%id,' of source',s%id,' inconsistent with source position'
                 stop 
              endif

           end do

        else

           do n=1,s%n_paths

              if (.not.(s%path(n)%sequence(2) == region(s%region_id)%itop%id .or. &
                   s%path(n)%sequence(2) == region(s%region_id)%ibot%id)) then
                 print *
                 print *,'***** ERROR  ***********************************'
                 print *,n,s%path(n)%sequence(2),region(s%region_id)%itop%id, &
                      region(s%region_id)%ibot%id
                 print *,'path',s%path(n)%id,' of source',s%id,' inconsistent with source position'
                 stop 
              endif
           
           end do

        endif

     endif

     if (s%is_teleseismic) then

        do n=1,s%n_paths

           if (.not.(s%path(n)%sequence(2) == n_interfaces-1  &
                .and. s%path(n)%sequence(1) == n_interfaces)) then
              print *
              print *,'***** ERROR  ***********************************'
              print *,'path',s%path(n)%id,' path from teleseismic source does not start at bottom'
              stop 
           endif

        end do

     endif


! then test paths for allowable sequence


     do n=1,s%n_paths

        do m=3,2*s%path(n)%n_tf,2

           if (abs(s%path(n)%sequence(m+1)-s%path(n)%sequence(m))/= 1.and.s%path(n)%sequence(m) /= 0 &
                .and. m /= 2*s%path(n)%refstep-1) then
              print *
              print *,'***** ERROR  ***********************************'
              print *, 'illegal path sequence for source',s%id,' path',n
              stop
           endif

           if ((s%path(n)%sequence(m+1) > n_interfaces .or. s%path(n)%sequence(m+1) < 1) .and. &
                m /= 2*s%path(n)%refstep-1) then
              print *
              print *,'***** ERROR  ***********************************'
              print *, 'non-existing interface in path sequence',m,s%path(n)%refstep
              stop
           endif

           if ((s%path(n)%sequence(m) /= s%path(n)%sequence(m-1)) .and. &
                (s%path(n)%sequence(m) /= s%path(n)%sequence(m-2))) then
              print *
              print *,'***** ERROR  ***********************************'
              print *, 'illegal path sequence for source',s%id,' path',n
              stop
           endif              

        end do
     end do


  end do  ! loop over input sources

close(1)

print *,'finished reading sources'


!---------------------------------------------------------
! read in the receiver properties


do n=1,n_receivers

   rec => receiver(n)

   read(2,*) rec%r,rec%lat,rec%long    ! the position of the receiver

   rec%r = earth_radius - rec%r
   rec%lat = rec%lat*deg_to_rad
   rec%long = rec%long*deg_to_rad

   if ( (rec%lat < pgrid%lat(1) .or. rec%lat > pgrid%lat(pgrid%nlat)) .or. &
        (rec%long < pgrid%long(1) .or. rec%long > pgrid%long(pgrid%nlong))) then
      print *
      print *,'error: receiver',n,' lies outside propagation grid in lat or long'
      stop

   else

      if ((rec%r < interpolate_interface(rec%lat,rec%long,intrface(n_interfaces))-0.1_dp*pgrid%tolerance) .or. &
        (rec%r >  interpolate_interface(rec%lat,rec%long,intrface(1))+0.1_dp*pgrid%tolerance) ) then 
         print *
         print *,'error: receiver',n,' lies above or below the propagation grid'
         stop
      endif

   endif

   read(2,*) rec%n_rays                               ! the number of paths to this receiver

   allocate(receiver(n)%ray(receiver(n)%n_rays))
   do i=1,receiver(n)%n_rays ; call ray_defaults(receiver(n)%ray(i)) ; end do

   read(2,*) rec%ray(1:rec%n_rays)%source_id  ! the index of the source of the rays

! verify that sources are valid

   if (count(rec%ray(1:rec%n_rays)%source_id > n_sources) > 0) then
      print *
      print *,'****** ERROR: INCONSISTENT INPUT ********'
      print *,'receiver',n,' is requesting paths from a non-existent source'
      stop
   endif

   do i=1,rec%n_rays ; rec%ray(i)%source => source(rec%ray(i)%source_id) ; end do

   read(2,*) rec%ray(1:rec%n_rays)%raypath_id  ! the index in the source path list of the rays

! verify that path references are valid

   do i=1,rec%n_rays 
      if (rec%ray(i)%raypath_id > rec%ray(i)%source%n_paths) then
         print *
         print *,'***** ERROR: INCONSISTENT INPUT *******************'
         print *,'the list of rays to receiver ',rec%id,' contains a path that is not defined'
         print *,'the number of the offending ray is ',i
         print *,'it refers to path',rec%ray(i)%raypath_id,' from source ',rec%ray(i)%source%id
         stop
      endif
   end do

end do

close(2)


! test for consistency between receiver positions and requested paths

do n=1,n_receivers

   rec => receiver(n)

   do m=1,rec%n_rays

      s => rec%ray(m)%source
      if (s%is_teleseismic) cycle
      path_id = rec%ray(m)%raypath_id
      nsec=s%path(path_id)%n_tf

      if (nsec > 1) then

         i1 = s%path(path_id)%sequence(2*nsec)
         i2 = s%path(path_id)%sequence(2*nsec-1)
         top_id = min(i1,i2)
         bot_id = max(i1,i2)
         if ( (rec%r > interpolate_interface(rec%lat,rec%long,intrface(top_id))+pgrid%tolerance) .or. &
              (rec%r < interpolate_interface(rec%lat,rec%long,intrface(bot_id))-pgrid%tolerance)) then
            print *,'receiver',n,' does not lie in final time field of path',m
            stop
         endif

      else

         if (s%on_interface) then

            top_id=max(1,s%topint_id-1)
            bot_id=min(n_interfaces,s%botint_id+1)
            if ( (rec%r > interpolate_interface(rec%lat,rec%long,intrface(top_id))+pgrid%tolerance) .or. &
                 (rec%r < interpolate_interface(rec%lat,rec%long,intrface(bot_id))-pgrid%tolerance)) then
               print *,'receiver',n,' does not lie in final time field of path',m
               stop
            endif

         else

            top_id=s%region_id
            bot_id=s%region_id+1
            if ( (rec%r > interpolate_interface(rec%lat,rec%long,intrface(top_id))+pgrid%tolerance) .or. &
                 (rec%r < interpolate_interface(rec%lat,rec%long,intrface(bot_id))-pgrid%tolerance)) then
               print *,'receiver',n,' does not lie in final time field of path',m
               stop
            endif

         endif

      endif

   end do

end do

print *,'finished reading receivers'



! test which paths are actually used, and pre-count the number of pp-phases arriving at this receiver

allocate(npp_receiver(n_receivers),ipp_receiver(n_receivers))

npp_receiver = 0    ! will contain the number of pp phases arriving at this receiver
ipp_receiver = 1    ! a counter for pp-phases

do n=1,n_receivers
   do m=1,receiver(n)%n_rays

      ss => receiver(n)%ray(m)%source
      path_id = receiver(n)%ray(m)%raypath_id

      ss%path(path_id)%used=.true.

      if (ss%path(path_id)%refstep /= 0) then  ! we have a pp type phase
         npp_receiver(n)=npp_receiver(n)+1
      endif

   end do
end do


! add receivers of pp-type phases to the source list, and construct their
! path sequences by inverting the tail end of the original pp-type sequence  

do n=1,n_receivers

   if (npp_receiver(n) > 0) allocate(receiver(n)%path_equivalent(receiver(n)%n_rays))

   do m=1,receiver(n)%n_rays

      ss => receiver(n)%ray(m)%source
      path_id = receiver(n)%ray(m)%raypath_id


      if (ss%path(path_id)%refstep /= 0) then  ! we have a pp type phase


      ! if the receiver does not yet have a source equivalent,create one and point to it

         if (receiver(n)%source_equivalent == 0) then

            n_sources_ppinc = n_sources_ppinc + 1
            receiver(n)%source_equivalent = n_sources_ppinc
            s => source(receiver(n)%source_equivalent)


            s%r=receiver(n)%r
            s%lat=receiver(n)%lat
            s%long=receiver(n)%long
            s%coslat=cos(s%lat)

            s%is_local = .true.

            call initialize_source(s,pgrid)

            if (.not. s%on_interface) print *,'virtual receiver source', &
                 n_sources_ppinc,' is located in region ',s%region_id

            if (s%on_interface) print *,'virtual receiver source', &
                 n_sources_ppinc,' is located on interfaces ',&
                 s%topint_id,s%botint_id

            receiver(n)%path_equivalent(m)= 1


         ! allocate storage for path info

            s%n_paths = npp_receiver(n)
            allocate(s%path(s%n_paths))
            do i=1,s%n_paths ; call path_defaults(s%path(i)) ; end do

         else   ! point to the already existing source equivalent

            s => source(receiver(n)%source_equivalent) 
            receiver(n)%path_equivalent(m) = ipp_receiver(n)

         endif


       ! set some attributes of the path

         recpath_id=ipp_receiver(n)

         s%path(recpath_id)%n_tf = ss%path(path_id)%n_tf - ss%path(path_id)%refstep 
         s%path(recpath_id)%id   = ipp_receiver(n)
         s%path(recpath_id)%valid   = .true.
         s%path(recpath_id)%used    = .true.
         s%path(recpath_id)%refstep = s%path(recpath_id)%n_tf + 1
         s%path(recpath_id)%fitting_interface = ss%path(path_id)%fitting_interface

         allocate(s%path(recpath_id)%sequence(2*s%path(recpath_id)%n_tf))
         allocate(s%path(recpath_id)%tf_sequence(s%path(recpath_id)%n_tf))
         allocate(s%path(recpath_id)%vtype_sequence(s%path(recpath_id)%n_tf))


       ! construct the path for the receiver equivalent source by inverting the orginal source path

       ! go step by step

         do i=1,2*s%path(recpath_id)%n_tf,2

            s%path(recpath_id)%sequence(i)   =  ss%path(path_id)%sequence(2*ss%path(path_id)%n_tf-i+1 )
            s%path(recpath_id)%sequence(i+1) =  ss%path(path_id)%sequence(2*ss%path(path_id)%n_tf-i )
            j=(i+1)/2
            s%path(recpath_id)%vtype_sequence(j) =   &
                 ss%path(path_id)%vtype_sequence(ss%path(path_id)%n_tf-j+1 )

          ! if it is a turning step exchange the interfaces in this step 

            if ( i > 1 ) then

               if ( s%path(recpath_id)%sequence(i+1) == s%path(recpath_id)%sequence(i-1)) then
                  j= s%path(recpath_id)%sequence(i+1)
                  s%path(recpath_id)%sequence(i+1) = s%path(recpath_id)%sequence(i)
                  s%path(recpath_id)%sequence(i) = j

               endif

            endif

         end do 

         s%path(recpath_id)%sequence(1)=0

         ipp_receiver(n)=ipp_receiver(n)+1

      endif

   end do

end do

deallocate(npp_receiver,ipp_receiver)


! now remove the unwanted sections from the original pp-type paths from the original sources

do n=1,n_sources
   s=> source(n)
   do m=1,s%n_paths
      if (s%path(m)%refstep /= 0) s%path(m)%n_tf=s%path(m)%refstep - 1
   end do
end do

! finished dealing with late reflections


! allocate the required timefields for all sources

do n=1,n_sources_ppinc

   i=n_vtypes*n_regions        ! max time fields from initial propagation through 
                               ! regions overlapping source grid
   do m=1,source(n)%n_paths
      i=i+source(n)%path(m)%n_tf
   end do
   allocate(source(n)%time_field(i))
   do j=1,i ; call time_field_defaults(source(n)%time_field(j)); source(n)%time_field(j)%id = j ; end do

end do

! if it is requested that some grids of arrival times be saved,
! read the information about which paths are to be svaed from the
! file gridsave.in

if (save_timefields_mode) then

   open(19,file='gridsave.in',iostat=io)
   if (io /= 0) then
      print *, 'save time fields mode specified, but there is a problem with the input file'
      print *, 'check that the file gridsave.in is present and not empty'
      stop 'the file gridsave.in does not appear to be present'
   endif

   save_grid_counter=0

   do
      read(19,*,iostat=io) n,nsave
      if (io<0) exit
      if (nsave>0) then
         allocate(iarr(nsave),jarr(nsave))
         read(19,*) iarr(1:nsave)    
         do m=1,nsave
            source(n)%path(iarr(m))%gridsave=.true.
            save_grid_counter=save_grid_counter+1
         end do
         read(19,*) jarr(1:nsave)    
         do m=1,nsave
            source(n)%path(iarr(m))%first_tf_to_save=jarr(m)
         end do
         deallocate(iarr,jarr)
      endif
   end do

   close(19)

   open(19,file='arrtimes.dat')
   write(19,*) pgrid%nr,pgrid%nlat,pgrid%nlong
   write(19,*) pgrid%dr0,pgrid%dlat0/deg_to_rad,pgrid%dlong0/deg_to_rad
   write(19,*) pgrid%r0,pgrid%lat0/deg_to_rad,pgrid%long0/deg_to_rad
   write(19,*) save_grid_counter
       
endif


! check if info for all sources is correct

! print *,'n_sources_ppinc',n_sources_ppinc

!do n=1,n_sources_ppinc
!   do m=1,source(n)%n_paths
!
!      print *
!      print *,'source',n,'  path', m, 'id', source(n)%path(m)%id
!      print *,'ntfinit',source(n)%n_tf_init
!      print *,'ntf =',source(n)%path(m)%n_tf,'refstep=',source(n)%path(m)%refstep
!      print '(a10,20i5)','path',source(n)%path(m)%sequence(1:2*source(n)%path(m)%n_tf)
!      print *,'fitting interface',source(n)%path(m)%fitting_interface
!
!   end do
! end do

! stop 'pptest stop'

! finished initializing!!!!!!!!!!!!!!!!!!!!
print *,'finsihed initializing sources and receivers'
print *

call system_clock(t2,count_rate,count_max)


if (display_mode) then

! create OpenDx input files showing the geometry of the problem

   do i=1,n_interfaces
      call display_interface(intersection(i))
   end do

   call display_sources
   call display_receivers

endif

call system_clock(t3,count_rate,count_max)
!--------------------------------------------------
! now calculate the required time fields for each source

global_source_counter = 0

! in no_pp_mode is true, the ray tracing will be done immediately after the
! fast marching for each source. It is called no_pp_mode because in this mode
! it is not possible to do late reflection, since this requires the time 
! fields from different sources at the same time

if (no_pp_mode) then

   call initialize_inversion(do_frechet_derivatives)

   open(11,file='arrivals.dat')
   open(21,file='frechet.dat')
   if (save_rays_mode) open(31,file='rays.dat')
   if (display_mode) open(41,file='raypos')
   if (display_mode) open(42,file='raytype')
   raypoint_counter=0

endif



do ns=1,n_sources_ppinc

   s => source(ns)

   if (parallel_mode) then
      if (mod(s%id-1,nproc)/=myid) cycle
   endif

! first evaluate the initial time fields of the sequences starting at the source

   if (s%is_local) then

      print *,'starting to initialize local source',s%id

! do the first sweep from the source through the regions overlapping source grid
! this is a big subroutine that does the entire grid refinement procedure around the source
! the call to this subroutine takes a significant part of the computation time

      call initialize_source_regions(s)

   endif


   if (s%is_teleseismic) then

!      print *,'starting to initialize teleseismic source',s%id,s%teleseismic_id
         
     call initialize_teleseismic_source(s)

   endif

   print *,'# time fields from initialization',s%n_time_fields


! we now have the first time fields for every possible sequence


! start the sequence of sweeps required for the paths

pathloop: do n=1,s%n_paths

   vtype = s%path(n)%vtype_sequence(1)     ! get the velocity type (P or S)
   print *

 ! identify the first time field of the path, created during initialisation above

   if (s%is_local) then

      if (s%on_interface) then

         if (s%path(n)%sequence(2) == s%topint_id - 1) prev_tf = s%first_tf_up(vtype)
         if (s%path(n)%sequence(2) == s%botint_id + 1) prev_tf = s%first_tf_down(vtype)

      else

         prev_tf = s%first_tf_up(vtype)

      endif

   endif

   if (s%is_teleseismic) prev_tf = s%first_tf_up(vtype)


   s%path(n)%tf_sequence(1) = prev_tf
   print *
   print '(a5,i4,a35,i4)', &
        'path',s%path(n)%id,' leg  1  using source time field',prev_tf  



! now go through the specified sequence of the path

   do m=3,2*s%path(n)%n_tf,2

      step=(m+1)/2
      vtype = s%path(n)%vtype_sequence(step)

!      print *
!      print *,'path=',n,'step=',step,'prev_tf=',prev_tf,'s%nt=',s%n_time_fields
!      print *,'children of prev_tf',s%time_field(prev_tf)%next_tf(1:4)

      istart => intersection(s%path(n)%sequence(m))      ! intersection at which the sweep starts
      inext  => intersection(s%path(n)%sequence(m+1))    ! intersection at which the sweep ends
      if (istart%id > inext%id) then
         reg    => istart%regabo
      else
         reg    => istart%regbel
      endif

      print '(a5,i4,a5,i4,a16,i4,a5,i4,a16,i4)', 'path',s%path(n)%id, &
           ' leg',step,' from interface',istart%id,'to',inext%id,'through region',reg%id


      ! set the child type of the next time field (up,down,turning,non-turning)

      if (step > 2) then   ! if the previous time field is not a source time field

         if ( istart%id == s%time_field(prev_tf)%inonstart%id ) then

            ! if the next required timefield is not derived from a turning ray
         
            if (istart%id > inext%id) ctype=1+(vtype-1)*4              ! propagate upwards
            if (istart%id < inext%id) ctype=2+(vtype-1)*4              ! propagate downwards

         else
         
            ! if the next required timefield is derived from a turning ray
       
            ! check if turning rays are indeed present
            if (.not.s%time_field(prev_tf)%turning_rays_present) then
               print '(a5,i4,a5,i4,a44)', &
                    'path',n,' leg',step,'no turning rays present, non-existing path'
               s%path(n)%valid=.false.
               print *,'step=',step,s%path(n)%sequence(1),prev_tf,s%time_field(prev_tf)%inonstart%id
               cycle pathloop
            endif

            if (istart%id > inext%id) ctype=3+(vtype-1)*4              ! propagate upwards
            if (istart%id < inext%id) ctype=4+(vtype-1)*4              ! propagate downwards

         endif

      else   ! only if the previous step is a source time field

         print *,'previous step is a source time field'

         
          if (istart%id == s%time_field(prev_tf)%reg%itop%id ) then

            if (istart%id > inext%id) ctype=1+(vtype-1)*4              ! propagate upwards
            if (istart%id < inext%id) ctype=2+(vtype-1)*4              ! propagate downwards

         else

            if (istart%id > inext%id) ctype=3+(vtype-1)*4              ! propagate upwards
            if (istart%id < inext%id) ctype=4+(vtype-1)*4              ! propagate downwards
               
         endif

      endif

!      print *,'ctype =',ctype,vtype

      if (s%time_field(prev_tf)%next_tf(ctype) == 0) then  ! if the required timefield does not exist

!         print *,'the requested field does not exist'

         ! transfer the starting times to the starting intersection
         ! here we use the fact that regional (and thus timefield) nodes are always in the order
         ! regular grid nodes, top intersection nodes, bottom intersection nodes

         if (istart%id == s%time_field(prev_tf)%reg%ibot%id) then
            i=s%time_field(prev_tf)%reg%nnode-istart%nnode+1
            j=s%time_field(prev_tf)%reg%nnode
            istart%arrivaltime=s%time_field(prev_tf)%arrivaltime(i:j)
            istart%time_gradient(1,:)=s%time_field(prev_tf)%time_gradient(1,i:j)
            istart%time_gradient(2,:)=s%time_field(prev_tf)%time_gradient(2,i:j)
            istart%time_gradient(3,:)=s%time_field(prev_tf)%time_gradient(3,i:j)
         else
            i=s%time_field(prev_tf)%reg%ngnode+1
            j=s%time_field(prev_tf)%reg%ngnode+istart%nnode
            istart%arrivaltime=s%time_field(prev_tf)%arrivaltime(i:j)
            istart%time_gradient(1,:)=s%time_field(prev_tf)%time_gradient(1,i:j)
            istart%time_gradient(2,:)=s%time_field(prev_tf)%time_gradient(2,i:j)
            istart%time_gradient(3,:)=s%time_field(prev_tf)%time_gradient(3,i:j)
         endif



         ! modify the time gradients at the interface for the next leg of the path

         if (reg%id == s%time_field(prev_tf)%reg%id) then

            ! if we go back into the same region as the previous time field, it is a reflection
            ! convert the direction of the gradient, and set the arrival time at the 
            ! intersection points where reflection is impossible (due to wave type conversion) 
            ! to huge_time so that they do not act as a source

 !           print *,'calling reflect-gradient at interface',istart%id

            call reflect_gradient(istart,s%time_field(prev_tf),vtype)

         else

            ! if the next region is another region, it is a refraction
            ! convert the direction of the gradient, and set the arrival time at the 
            ! intersection points where total reflection occurs to huge_time 
            ! so that they do not act as a source

            direction = s%time_field(prev_tf)%reg%id - reg%id

!            print *,'calling refract_gradient with direction',direction

            call refract_gradient(istart,s%time_field(prev_tf)%reg,vtype,direction)   

            
         endif


     ! we are finally set up, now do the actual fast marching across the region 
     ! generating a new time field
     !-------------------------------------------------------------------------

         call sweep_region_from_interface(reg,istart,vtype,s)


     ! attach info to the generated time field

         ! time field preceding the current one
         s%time_field(s%n_time_fields)%prev_tf = prev_tf

         ! identify the current time field as the child of the previous time field
         s%time_field(prev_tf)%next_tf(ctype) = s%n_time_fields

         ! store the index of the time field in the sequence of timefields for this path
         s%path(n)%tf_sequence(step) = s%n_time_fields

         ! pointers to start and non-start interfaces
         s%time_field(s%n_time_fields)%istart => istart
         s%time_field(s%n_time_fields)%inonstart =>  inext  
         s%time_field(s%n_time_fields)%reg =>  reg  
         s%time_field(s%n_time_fields)%vtype = vtype

         prev_tf = s%n_time_fields
         print *,'created new time field',s%n_time_fields

      else  ! the required timefield already exists

         ! step to the next (existing) time field
         prev_tf = s%time_field(prev_tf)%next_tf(ctype)
         print *,'used existing time field',prev_tf

         ! store the index of the time field in the sequence of timefields for this path
         s%path(n)%tf_sequence(step) = prev_tf
         
      endif

   end do   ! crossings in path

end do pathloop  ! path loop


do n=1,s%n_paths
   print '(a5,i5,a12,10i5)','path',s%path(n)%id,' timefields',s%path(n)%tf_sequence(1:s%path(n)%n_tf)
end do

! save the arrival times on the main grid for this path if requested

if (save_timefields_mode) then
   do n=1,s%n_paths
       if (s%path(n)%gridsave.and.s%path(n)%valid) then
         call write_arrivaltime_grid(s,s%path(n))
      endif
   end do
endif


if (file_mode) then

! store the timefields of this source on file
   global_source_counter =  global_source_counter + 1
   nfile =  global_source_counter + 1000
   s%nfile=nfile
   open(nfile,form='unformatted')

   do n=1,s%n_time_fields
      write(nfile) s%time_field(n)%arrivaltime
      write(nfile) s%time_field(n)%time_gradient
      deallocate(s%time_field(n)%arrivaltime,s%time_field(n)%time_gradient)
   end do

   close(nfile)

endif


if (no_pp_mode) then

   print *,'starting the ray tracing for source',s%id

   do n=1,n_receivers

      do m=1,receiver(n)%n_rays

         ray => receiver(n)%ray(m)
         if (ray%source%id /= s%id) cycle 
         path_id = ray%raypath_id

         if (s%path(path_id)%valid) then ! the original path was recognised as valid 
                                         ! during the timefield calculations

            if (s%path(path_id)%refstep == 0) then    !  the standard case of no reflection fitting

               call trace_ray_from_receiver(receiver(n),s,ray)

               if (ray%valid) then   ! valid ray path found

                  print '(a12,i4,a10,i4,a15,i4,a4,f10.4,2l5)','traced ray',m,'to source',s%id,&
                       ' from receiver',n,'  t=', ray%receiver_time,ray%diffracted,ray%headwave

                  k=0
                  write(11,'(4i6,f15.6,2l5)') n,ray%source%id,m,k,ray%receiver_time, &
                       ray%diffracted,ray%headwave

                  if (display_mode) call store_ray(ray)

               else   ! ray tracing found that this ray path does not exist

                  print '(a12,i4,a10,i4,a15,i4,a13)','ray',m,'to source',ray%source_id,&
                       ' from receiver',n,' is invalid'

                  k=0
                  t_arrival=-1.0_dp               
                  write(11,'(4i6,f15.6,2l5)') n,ray%source%id,m,k,t_arrival,ray%diffracted,ray%headwave

               endif

            endif
         
            if (s%path(path_id)%refstep /= 0) then  !  if reflection fitting is required..

               stop 'trying to preform reflection fit while no_pp_mode enabled'

               call trace_reflectionfit(n,m)

            endif


   ! do the frechet derivatives for this ray if required

            if (do_frechet_derivatives) then

               if (receiver(n)%ray(m)%is_multiray) then

                  stop 'inconsistency: multiray in no_pp_mode'

               else

                  ray => receiver(n)%ray(m)
                  if (ray%valid) then
                     print *, 'getting partials for rec',n,'ray',m
                     call ray_partials(ray)
                  else
                     print *,'no valid ray path for rec',n,'ray',m
                  endif

               endif

               call write_frechet_derivatives(n,m)

            endif


         else           ! the original path was recognised as invalid during the timefield calculations

            print *,'ray',m,' to receiver ',n, &
                 ' : requested path was recognised as invalid during the timefield calculations'

            ray%valid = .false.
            k=0
            t_arrival=-1.0_dp               
            write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,k,t_arrival,ray%diffracted,ray%headwave

         endif

         if (save_rays_mode) call write_valid_rays(n,m)

         call clean_ray(n,m)

      end do    ! loop over rays/paths

   end do    !loop over receivers


endif

! deallocate source specific stuff that is not needed any more

do n=1,s%n_time_fields
   if (associated(s%time_field(n)%received_turning_ray)) &
        deallocate(s%time_field(n)%received_turning_ray)
end do

if (s%is_teleseismic) then
   do n=1,n_regions
      reg=>region(n)
      if (reg%n_init > 0) deallocate(reg%init_id,reg%init_type,reg%init_arrivaltime, &
           reg%init_time_gradient)
   end do
   reg%n_init=0
endif

if (no_pp_mode) then
   do n=1,s%n_time_fields
      deallocate(s%time_field(n)%arrivaltime,s%time_field(n)%time_gradient)
   end do
endif

if (no_pp_mode) then
   print *,'finished  fast marching and ray tracing for source',s%id   
else
   print *,'finished the fast marching for source',s%id
endif

print *,'******************************************************************************'
print *

end do ! loop over sources


if (no_pp_mode) then

   close(11)
   close(21)
   if (save_rays_mode) close(31)
   if (display_mode) close(41)
   if (display_mode) close(42)

   if (display_mode) call display_stored_rays

endif

!! finished the time fields
!******************************************************************************************************

  call system_clock(t4,count_rate,count_max)


! -------------- now do the ray tracing if not in no_pp_mode  --------------------------------------------------------------

if (n_receivers > 0 .and. (.not.no_pp_mode)) then

   call initialize_inversion(do_frechet_derivatives)

   open(11,file='arrivals.dat')
   open(21,file='frechet.dat')
   if (save_rays_mode) open(31,file='rays.dat')
   if (display_mode) open(41,file='raypos')
   if (display_mode) open(42,file='raytype')
   raypoint_counter=0

   print *,'starting the ray tracing'

   do n=1,n_receivers

      do m=1,receiver(n)%n_rays

         ray => receiver(n)%ray(m)
         s   => ray%source
         path_id = ray%raypath_id

         if (s%path(path_id)%valid) then ! the original path was recognised as valid 
                                         ! during the timefield calculations

            if (s%path(path_id)%refstep == 0) then    !  the standard case of no reflection fitting

               if (file_mode) call load_source_timefields(s)
               call trace_ray_from_receiver(receiver(n),s,ray)

               if (ray%valid) then   ! valid ray path found

                  print '(a12,i4,a10,i4,a15,i4,a4,f10.4,2l5)','traced ray',m,'to source',ray%source_id,&
                       ' from receiver',n,'  t=', ray%receiver_time,ray%diffracted,ray%headwave

                  k=0
                  write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,k,ray%receiver_time,ray%diffracted,ray%headwave

                  if (display_mode) call store_ray(ray)

               else   ! ray tracing found that this ray path does not exist

                  print '(a12,i4,a10,i4,a15,i4,a13)','ray',m,'to source',ray%source_id,&
                       ' from receiver',n,' is invalid'

                  k=0
                  t_arrival=-1.0_dp               
                  write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,k,t_arrival,ray%diffracted,ray%headwave

               endif

            endif
         
            if (s%path(path_id)%refstep /= 0) then  !  if reflection fitting is required..

               if (file_mode) call load_source_timefields(receiver(n)%ray(m)%source)
               if (file_mode) call load_source_timefields(source(receiver(n)%source_equivalent))

               call trace_reflectionfit(n,m)

            endif


   ! do the frechet derivatives for this ray if required

            if (do_frechet_derivatives) then

               if (receiver(n)%ray(m)%is_multiray) then

                  do k=1,receiver(n)%ray(m)%n_subrays
                     ray => receiver(n)%ray(m)%subray(k)
                     if (ray%valid) then
                        print *, 'getting partials for rec',n,'ray',m,'subray',k
                        call ray_partials(ray)
                     else
                        print *, 'no valid ray path for rec',n,'ray',m,'subray',k
                     endif
                  end do

               else

                  ray => receiver(n)%ray(m)
                  if (ray%valid) then
                     print *, 'getting partials for rec',n,'ray',m
                     call ray_partials(ray)
                  else
                     print *,'no valid ray path for rec',n,'ray',m
                  endif

               endif

               call write_frechet_derivatives(n,m)

            endif


            if (file_mode) call clean_source_timefields(receiver(n)%ray(m)%source)
            if (file_mode .and. s%path(path_id)%refstep /= 0)  &
                 call clean_source_timefields(source(receiver(n)%source_equivalent))


         else           ! the original path was recognised as invalid during the timefield calculations

            print *,'ray',m,' to receiver ',n, &
                 ' : requested path was recognised as invalid during the timefield calculations'

            ray%valid = .false.
            k=0
            t_arrival=-1.0_dp               
            write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,k,t_arrival,ray%diffracted,ray%headwave

         endif

         if (save_rays_mode) call write_valid_rays(n,m)

         call clean_ray(n,m)

      end do    ! loop over rays/paths

   end do    !loop over receivers

   close(11)
   close(21)
   if (save_rays_mode) close(31)
   if (display_mode) close(41)
   if (display_mode) close(42)

   call system_clock(t5,count_rate,count_max)

   print *,'init time       :',dble(t2-t1)/dble(count_rate),' sec'

   print *,'propagation time:',dble(t4-t3)/dble(count_rate),' sec'

   print *,'ray tracing time:',dble(t5-t4)/dble(count_rate),' sec'

   if (display_mode) call display_stored_rays


endif  ! n_receivers > 0 and not in no_pp_mode


!-----------------------------------------------------------------------------------------------------------

end program fm3d

3dfmlib.f90
!*************************************************************
! This subroutine reads the initial parametrization of the velocity
! fields from which the values on the propagation grid are to be interpolated
! from file into the appropriate structures

subroutine initialize_velocity_grids
use mod_3dfm
implicit none

integer :: n,m,i,j,k

open(10,file='vgrids.in')

! read the velocity mode from the file

! read the number of regions in the input velocity structure
read (10,*) n_vgrids,n_vtypes

! allocate space for these regions
allocate(vgrid(n_vgrids,n_vtypes))
do m=1,n_vtypes
   do n=1,n_vgrids 
      call vgrid_defaults(vgrid(n,m)) 
   end do
end do
! read the grid properties and velocity values to be interpolated for each region

do m=1,n_vtypes

   do n=1,n_vgrids

      ! grid parameters
      read(10,*) vgrid(n,m)%nr,vgrid(n,m)%nlat,vgrid(n,m)%nlong
      read(10,*) vgrid(n,m)%dr0,vgrid(n,m)%dlat0,vgrid(n,m)%dlong0
      read(10,*) vgrid(n,m)%r0,vgrid(n,m)%lat0,vgrid(n,m)%long0
      print *,'Velocity Grid parameters\n'
      print *,vgrid(n,m)%nr,vgrid(n,m)%nlat,vgrid(n,m)%nlong
      print *,vgrid(n,m)%dr0,vgrid(n,m)%dlat0,vgrid(n,m)%dlong0
      print *,vgrid(n,m)%r0,vgrid(n,m)%lat0,vgrid(n,m)%long0


! initialize the grid

      allocate(vgrid(n,m)%r(vgrid(n,m)%nr),vgrid(n,m)%lat(vgrid(n,m)%nlat), &
           vgrid(n,m)%long(vgrid(n,m)%nlong))

      do i=1,vgrid(n,m)%nr
         vgrid(n,m)%r(i)=vgrid(n,m)%r0 + (i-1)*vgrid(n,m)%dr0
      end do

      do i=1,vgrid(n,m)%nlat
         vgrid(n,m)%lat(i)=vgrid(n,m)%lat0 + (i-1)*vgrid(n,m)%dlat0
      end do

      do i=1,vgrid(n,m)%nlong
         vgrid(n,m)%long(i)=vgrid(n,m)%long0 + (i-1)*vgrid(n,m)%dlong0
      end do

! read in the velocity values on the interpolation grid

      allocate(vgrid(n,m)%velocity(vgrid(n,m)%nr,vgrid(n,m)%nlat,vgrid(n,m)%nlong))

      do i=1,vgrid(n,m)%nr
         do j=1,vgrid(n,m)%nlat
            do k=1,vgrid(n,m)%nlong
               read (10,*) vgrid(n,m)%velocity(i,j,k)
            end do
         end do
      end do

      if (count(vgrid(n,m)%velocity > 20.0_dp) > 0 ) &
           print *,'*** WARNING *** : velocity grid contains values larger than 20 km/sec'
      if (count(vgrid(n,m)%velocity < 1.0_dp) > 0 ) &
           print *,'*** WARNING *** : velocity grid contains values less than 1 km/sec'

      vgrid(n,m)%nnode=vgrid(n,m)%nr*vgrid(n,m)%nlat*vgrid(n,m)%nlong

 ! allocate and initialize the activity flag

      allocate(vgrid(n,m)%active(vgrid(n,m)%nr,vgrid(n,m)%nlat,vgrid(n,m)%nlong))
      vgrid(n,m)%active = .false.


   end do  ! loop over interpolation regions

end do  ! vtypes

close(10)

end subroutine initialize_velocity_grids


!***********************************************************************************************
! reads the parameters of the propagation grid from file and initializes the grid
subroutine initialize_propagation_grid
use mod_3dfm
implicit none

integer :: i,j,k
real(kind=dp) :: deg_to_rad

open(10,file='propgrid.in')

allocate(pgrid)
call pgrid_defaults(pgrid)

! grid parameters
read(10,*) pgrid%nr,pgrid%nlat,pgrid%nlong
read(10,*) pgrid%dr0,pgrid%dlat0,pgrid%dlong0
read(10,*) pgrid%r0,pgrid%lat0,pgrid%long0

deg_to_rad=acos(-1.0_dp)/180.0_dp

pgrid%dlat0=pgrid%dlat0*deg_to_rad
pgrid%dlong0=pgrid%dlong0*deg_to_rad
pgrid%lat0=pgrid%lat0*deg_to_rad
pgrid%long0=pgrid%long0*deg_to_rad

pgrid%r0 =  earth_radius + pgrid%r0 - dble(pgrid%nr-1)*pgrid%dr0

pgrid%tolerance=interface_tolerance*pgrid%dr0

pgrid%rmax = pgrid%r0 + (pgrid%nr-1)*pgrid%dr0
pgrid%latmax = pgrid%lat0 + (pgrid%nlat-1)*pgrid%dlat0
pgrid%longmax = pgrid%long0 + (pgrid%nlong-1)*pgrid%dlong0

read(10,*) refinement_factor,ncell_to_be_refined

! initialize the grid

allocate(pgrid%r(pgrid%nr),pgrid%lat(pgrid%nlat),pgrid%coslat(pgrid%nlat),pgrid%long(pgrid%nlong))

do i=1,pgrid%nr
   pgrid%r(i)=pgrid%r0 + (i-1)*pgrid%dr0
end do

do i=1,pgrid%nlat
   pgrid%lat(i)=pgrid%lat0 + (i-1)*pgrid%dlat0
end do

do i=1,pgrid%nlat
   pgrid%coslat(i)=cos(pgrid%lat(i))
end do


do i=1,pgrid%nlong
   pgrid%long(i)=pgrid%long0 + (i-1)*pgrid%dlong0
end do


close(10)


allocate(pgrid%rnode_id(pgrid%nr,pgrid%nlat,pgrid%nlong))
allocate(pgrid%node_region(pgrid%nr,pgrid%nlat,pgrid%nlong))
pgrid%node_region = 0
allocate(pgrid%ccind_from_3dc(pgrid%nr,pgrid%nlat,pgrid%nlong))
do k=1,pgrid%nlong
   do j=1,pgrid%nlat
      do i=1,pgrid%nr
         nullify(pgrid%ccind_from_3dc(i,j,k)%p)
      end do
   end do
end do

pgrid%is_main_grid = .true.

end subroutine initialize_propagation_grid



!*************************************************************
! This subroutine reads the initial parametrization of the interface positions
! from which the values on the propagation grid are to be interpolated
! from file into the appropriate structures

subroutine initialize_interfaces
use mod_3dfm
implicit none

integer :: n,i,j
integer :: nlat,nlong
real(kind=dp) :: dlat0,dlong0,lat0,long0,h,hb

open(10,file='interfaces.in')


! read the number of interfaces
read (10,*) n_interfaces

! allocate space for these interfaces (and the associated intersections and regions for future use)
allocate(intrface(n_interfaces))
do n=1,n_interfaces ; call interface_defaults(intrface(n)) ; intrface(n)%id=n ; end do


! read the grid properties and radius values to be interpolated for the internal interfaces

! grid parameters
   read(10,*) nlat,nlong
   read(10,*) dlat0,dlong0
   read(10,*) lat0,long0


do n=1,n_interfaces

   intrface(n)%nlat = nlat
   intrface(n)%nlong = nlong
   intrface(n)%dlat0 = dlat0
   intrface(n)%dlong0 = dlong0
   intrface(n)%lat0 = lat0
   intrface(n)%long0 = long0


! initialize the grid

   allocate(intrface(n)%lat(intrface(n)%nlat),intrface(n)%long(intrface(n)%nlong))

   do i=1,intrface(n)%nlat
      intrface(n)%lat(i)=intrface(n)%lat0 + (i-1)*intrface(n)%dlat0
   end do

   do i=1,intrface(n)%nlong
      intrface(n)%long(i)=intrface(n)%long0 + (i-1)*intrface(n)%dlong0
   end do


! read in the radius values on the interpolation grid

   allocate(intrface(n)%r(intrface(n)%nlat,intrface(n)%nlong))

   do i=1,intrface(n)%nlat
      do j=1,intrface(n)%nlong
         read (10,*) intrface(n)%r(i,j)
      end do
   end do

   intrface(n)%nnode=intrface(n)%nlat*intrface(n)%nlong

end do  ! loop over interfaces

close(10)

! check top and bottom interfaces are not outside the propagation grid

     if (count(intrface(1)%r > pgrid%r(pgrid%nr)) > 0) stop ' ERROR: surface above propagation grid'
     if (count(intrface(n_interfaces)%r < pgrid%r(1)) > 0) stop ' ERROR: bottom below propagation grid'


! correct for intersecting interfaces, higher takes priority

  do n=2,n_interfaces

     do j=1,intrface(n)%nlong
        do i=1,intrface(n)%nlat

    ! higher has priority, EXCEPT bottom interface

           if (n < n_interfaces) then

              hb=intrface(n_interfaces)%r(i,j) 

              if (intrface(n)%r(i,j) < hb) then
                 intrface(n)%r(i,j) = hb
                 intrface(n)%pinched = .true.
                 intrface(n_interfaces)%pinched = .true.
              endif

           endif

    ! check if interface above is crossed

           h=intrface(n-1)%r(i,j)

           if (intrface(n)%r(i,j) > h) then
              intrface(n)%r(i,j) = h
              intrface(n)%pinched = .true.
              intrface(n-1)%pinched = .true.
           endif

        end do
     end do

  end do


end subroutine initialize_interfaces


!*****************************************************
! this function returns the position (radius) of an interface at a horizontal position
! using 2D bicubic spline interpolation

function interpolate_interface(lat,long,iface)
use mod_3dfm
implicit none

real(kind=dp) :: interpolate_interface
real(kind=dp) :: lat,long
type(Tinterface)  :: iface

integer        :: i,j,ilat,ilong

real(kind=dp) :: u,v,bu(4),bv(4),value

ilat=floor((lat-iface%lat0)/iface%dlat0)+1
ilong=floor((long-iface%long0)/iface%dlong0)+1

if (ilong < 2 .or. ilong > (iface%nlong-2)) then 
   print *,'interpolate_interface : interpolation outside range :ilong'
   print *,iface%id,ilong,ilat
   print *,long,lat
   print *,'if this happens during initialization your interface parameter '
   print *,'grid may not cover the entire propagation grid'
   stop
endif

if (ilat < 2 .or. ilat > (iface%nlat-2)) then
   print *,'interpolate_interface : interpolation outside range:ilat',iface%id
   print *,ilat,iface%nlat-2, lat,iface%lat0+(ilat-1)*iface%dlat0,iface%lat0+(ilat)*iface%dlat0
   print *,'if this happens during initialization your interface parameter '
   print *,'grid may not cover the entire propagation grid'
   stop
endif

u=(lat-iface%lat(ilat))/iface%dlat0
v=(long-iface%long(ilong))/iface%dlong0

bu(1)=(1.0_dp-u)**3/6.0_dp
bu(2)=(4.0_dp-6.0_dp*u**2+3.0_dp*u**3)/6.0_dp
bu(3)=(1.0_dp+3.0*u+3.0_dp*u**2-3.0_dp*u**3)/6.0_dp
bu(4)=u**3/6.0_dp
bv(1)=(1.0_dp-v)**3/6.0_dp
bv(2)=(4.0_dp-6.0_dp*v**2+3.0_dp*v**3)/6.0_dp
bv(3)=(1.0_dp+3.0*v+3.0_dp*v**2-3.0_dp*v**3)/6.0_dp
bv(4)=v**3/6.0_dp


value=0.0_dp
do j=1,4
   do i=1,4
      value=value+bu(i)*bv(j)*iface%r(ilat+i-2,ilong+j-2)
   end do
end do

interpolate_interface=value

end function interpolate_interface
!*****************************************************

! this subroutine returns the upward normal of an interface at a horizontal position
! using 2D bicubic spline interpolation

subroutine interface_normal(lat,long,iface,norm_r,norm_lat,norm_long,h)
use mod_3dfm
implicit none

real(kind=dp) :: lat,long,norm_r,norm_lat,norm_long
type(Tinterface)  :: iface

integer        :: i,j,ilat,ilong

real(kind=dp) :: u,v,bu(4),bv(4),h,dhdu,dhdv,bpu(4),bpv(4),dhdlat,dhdlong,norm

ilat=floor((lat-iface%lat0)/iface%dlat0)+1
ilong=floor((long-iface%long0)/iface%dlong0)+1

if (ilong < 2 .or. ilong > (iface%nlong-2)) then 
   print *,'interface_normal : interpolation outside range :ilong'
   print *,iface%id,ilong,ilat
   print *,long,lat
   stop
endif

if (ilat < 2 .or. ilat > (iface%nlat-2)) then
   print *,'interface_normal : interpolation outside range:ilat',iface%id
   print *,ilat,iface%nlat-2, lat,iface%lat0+(ilat-1)*iface%dlat0,iface%lat0+(ilat)*iface%dlat0
   stop
endif

u=(lat-iface%lat(ilat))/iface%dlat0
v=(long-iface%long(ilong))/iface%dlong0

bu(1)=(1.0_dp-u)**3/6.0_dp
bu(2)=(4.0_dp-6.0_dp*u**2+3.0_dp*u**3)/6.0_dp
bu(3)=(1.0_dp+3.0*u+3.0_dp*u**2-3.0_dp*u**3)/6.0_dp
bu(4)=u**3/6.0_dp
bv(1)=(1.0_dp-v)**3/6.0_dp
bv(2)=(4.0_dp-6.0_dp*v**2+3.0_dp*v**3)/6.0_dp
bv(3)=(1.0_dp+3.0*v+3.0_dp*v**2-3.0_dp*v**3)/6.0_dp
bv(4)=v**3/6.0_dp

bpu(1)=-0.5_dp*(1.0_dp-u)**2
bpu(2)=-2.0_dp*u+1.5_dp*u**2
bpu(3)=0.5_dp+u-1.5_dp*u**2
bpu(4)=0.5_dp*u**2
bpv(1)=-0.5_dp*(1.0_dp-v)**2
bpv(2)=-2.0_dp*v+1.5_dp*v**2
bpv(3)=0.5_dp+v-1.5_dp*v**2
bpv(4)=0.5_dp*v**2


h=0.0_dp
dhdu=0.0_dp
dhdv=0.0_dp
do j=1,4
   do i=1,4
      h   =h   +bu(i)*bv(j)*iface%r(ilat+i-2,ilong+j-2)
      dhdu=dhdu+bpu(i)*bv(j)*iface%r(ilat+i-2,ilong+j-2)
      dhdv=dhdv+bu(i)*bpv(j)*iface%r(ilat+i-2,ilong+j-2)
   end do
end do

dhdlat=dhdu/(h*iface%dlat0)
dhdlong=dhdv/(h*cos(lat)*iface%dlong0)

norm=sqrt(1.0_dp+dhdlat**2+dhdlong**2)

norm_r=1.0_dp/norm
norm_lat=-dhdlat/norm
norm_long=-dhdlong/norm

end subroutine interface_normal


!*****************************************************
! this function returns the velocity (propagation speed) at a 3D position
! using 3D bicubic spline interpolation

function interpolate_velocity(r,lat,long,gridv)
use mod_3dfm
implicit none


real(kind=dp) :: r,lat,long,interpolate_velocity
type(Tvelocity_grid)        :: gridv

integer        :: i,j,k,ir,ilat,ilong

real(kind=dp) :: u,v,w,bu(4),bv(4),bw(4),value

ir=floor((r-gridv%r0)/gridv%dr0)+1
ilat=floor((lat-gridv%lat0)/gridv%dlat0)+1
ilong=floor((long-gridv%long0)/gridv%dlong0)+1

!print *,ir,r,gridv%r0,gridv%dr0

if (ir < 2 .or. ir > (gridv%nr-2)) then
   print *,r,ir,gridv%nr
   print *,gridv%r(gridv%nr),gridv%r(gridv%nr-1),gridv%r(1),gridv%r(2)
   print *, 'interpolate_velocity : interpolation outside range ir'
   print *,'if this happens during initialization your velocity parameter '
   print *,'grid may not cover the entire propagation grid'
   stop
endif

if (ilong < 2 .or. ilong > (gridv%nlong-2))  then
   print *, 'interpolate_velocity : interpolation outside range ilong'
   print *,'if this happens during initialization your velocity parameter '
   print *,'grid may not cover the entire propagation grid'
   stop
endif

if (ilat < 2 .or. ilat > (gridv%nlat-2))  then
   print *, 'interpolate_velocity : interpolation outside range ilat'
   print *,'if this happens during initialization your velocity parameter '
   print *,'grid may not cover the entire propagation grid'
   stop
endif

u=(r-gridv%r(ir))/gridv%dr0
v=(lat-gridv%lat(ilat))/gridv%dlat0
w=(long-gridv%long(ilong))/gridv%dlong0

bu(1)=(1.0_dp-u)**3/6.0_dp
bu(2)=(4.0_dp-6.0_dp*u**2+3.0_dp*u**3)/6.0_dp
bu(3)=(1.0_dp+3.0*u+3.0_dp*u**2-3.0_dp*u**3)/6.0_dp
bu(4)=u**3/6.0_dp
bv(1)=(1.0_dp-v)**3/6.0_dp
bv(2)=(4.0_dp-6.0_dp*v**2+3.0_dp*v**3)/6.0_dp
bv(3)=(1.0_dp+3.0*v+3.0_dp*v**2-3.0_dp*v**3)/6.0_dp
bv(4)=v**3/6.0_dp
bw(1)=(1.0_dp-w)**3/6.0_dp
bw(2)=(4.0_dp-6.0_dp*w**2+3.0_dp*w**3)/6.0_dp
bw(3)=(1.0_dp+3.0*w+3.0_dp*w**2-3.0_dp*w**3)/6.0_dp
bw(4)=w**3/6.0_dp


value=0.0_dp
do k=1,4
   do j=1,4
      do i=1,4
         value=value+bu(i)*bv(j)*bw(k)*gridv%velocity(ir+i-2,ilat+j-2,ilong+k-2)
      end do
   end do
end do


interpolate_velocity=value   !(5.0_dp*6271.0_dp)/r 


end function interpolate_velocity
!*****************************************************

! this subroutine returns the gradient of the propagation speed at a 3d position
! using 3D bicubic spline interpolation

subroutine velocity_gradient(r,lat,long,gridv,dvdr,dvdlat,dvdlong,vel)
use mod_3dfm
implicit none

real(kind=dp) :: r,lat,long,dvdr,dvdlat,dvdlong,vel
type(Tvelocity_grid)  :: gridv

integer        :: i,j,k,ir,ilat,ilong

real(kind=dp) :: u,v,w,bu(4),bv(4),bw(4),bpu(4),bpv(4),bpw(4),dvdu,dvdv,dvdw,vv

ir=floor((r-gridv%r0)/gridv%dr0)+1
ilat=floor((lat-gridv%lat0)/gridv%dlat0)+1
ilong=floor((long-gridv%long0)/gridv%dlong0)+1

if (ir < 2 .or. ir > (gridv%nr-2)) then
   print *,r,ir
   print *,gridv%r(gridv%nr),gridv%r(gridv%nr-1),gridv%r(1),gridv%r(2)
   stop ' velocity_gradient: interpolation outside range ir'
endif
if (ilong < 2 .or. ilong > (gridv%nlong-2)) then 
   print *,'velocity_gradient : interpolation outside range :ilong'
   print *,ilong,ilat
   print *,long,lat
   stop
endif

if (ilat < 2 .or. ilat > (gridv%nlat-2)) then
   print *,' velocity_gradient: interpolation outside range'
   print *,ilat
   stop
endif

u=(r-gridv%r(ir))/gridv%dr0
v=(lat-gridv%lat(ilat))/gridv%dlat0
w=(long-gridv%long(ilong))/gridv%dlong0

bu(1)=(1.0_dp-u)**3/6.0_dp
bu(2)=(4.0_dp-6.0_dp*u**2+3.0_dp*u**3)/6.0_dp
bu(3)=(1.0_dp+3.0*u+3.0_dp*u**2-3.0_dp*u**3)/6.0_dp
bu(4)=u**3/6.0_dp
bv(1)=(1.0_dp-v)**3/6.0_dp
bv(2)=(4.0_dp-6.0_dp*v**2+3.0_dp*v**3)/6.0_dp
bv(3)=(1.0_dp+3.0*v+3.0_dp*v**2-3.0_dp*v**3)/6.0_dp
bv(4)=v**3/6.0_dp
bw(1)=(1.0_dp-w)**3/6.0_dp
bw(2)=(4.0_dp-6.0_dp*w**2+3.0_dp*w**3)/6.0_dp
bw(3)=(1.0_dp+3.0*w+3.0_dp*w**2-3.0_dp*w**3)/6.0_dp
bw(4)=w**3/6.0_dp


bpu(1)=-0.5_dp*(1.0_dp-u)**2
bpu(2)=-2.0_dp*u+1.5_dp*u**2
bpu(3)=0.5_dp+u-1.5_dp*u**2
bpu(4)=0.5_dp*u**2
bpv(1)=-0.5_dp*(1.0_dp-v)**2
bpv(2)=-2.0_dp*v+1.5_dp*v**2
bpv(3)=0.5_dp+v-1.5_dp*v**2
bpv(4)=0.5_dp*v**2
bpw(1)=-0.5_dp*(1.0_dp-w)**2
bpw(2)=-2.0_dp*w+1.5_dp*w**2
bpw(3)=0.5_dp+w-1.5_dp*w**2
bpw(4)=0.5_dp*w**2

vv=0.0_dp
dvdu=0.0_dp
dvdv=0.0_dp
dvdw=0.0_dp
do j=1,4
   do i=1,4
      do k=1,4
         vv  = vv   + bu(i)*bv(j)*bw(k)*gridv%velocity(ir+i-2,ilat+j-2,ilong+k-2)
         dvdu= dvdu + bpu(i)*bv(j)*bw(k)*gridv%velocity(ir+i-2,ilat+j-2,ilong+k-2)
         dvdv= dvdv + bu(i)*bpv(j)*bw(k)*gridv%velocity(ir+i-2,ilat+j-2,ilong+k-2)
         dvdw= dvdw + bu(i)*bv(j)*bpw(k)*gridv%velocity(ir+i-2,ilat+j-2,ilong+k-2)
      end do
   end do
end do

vel    = vv
dvdr   = dvdu
dvdlat = dvdv/r
dvdlong= dvdw/(r*cos(lat))

return

end subroutine velocity_gradient


!*****************************************************
!****************************************************************************************************


! ****************************************************************************************
! this subroutine takes a grid and an interface as input, and produces
! the intersection between the two as output
! the intersection is a collection of nodes that lie on the intersections
! between the grid connections and the interface
! in addition, each intersection contains pointers to the grid cells cut by the interface
! and arrays containing pointers from interface nodes to the cut grid cells they are
! associated with, and the inverse, pointers from the cut grid cells to the interface nodes
! associated with them

subroutine find_intersection(isec,iface,grid)
use mod_3dfm_nointerfaces
implicit none

! subroutine arguments
type(Tintersection)              :: isec
type(Tinterface)                 :: iface
type(Tpropagation_grid),target   :: grid

! local arrays
real(kind=dp), dimension(:,:,:),allocatable            :: rdiff
real(kind=dp), dimension(:,:)  ,allocatable            :: r_interface
type(Tinteger_coordinates), dimension(:)  ,allocatable :: ijk_isec
integer, dimension(:)  ,allocatable                    :: type_isec 

integer,dimension(:,:,:), allocatable                  :: nextranodes 
integer,dimension(:,:,:,:),allocatable                 :: extranodes
integer,dimension(:),allocatable                       :: counter

! local variables
integer                                                :: i,j,k,n,m,nodecount
real(kind=dp)                                          :: rint
logical                                                :: diag = .false.

! external functions
real(kind=dp)                                          :: interpolate_interface


! id of interface on which the intersection is based

isec%iface_id = iface%id

! pointer to grid on which this intersection is defined

isec%grid => grid

! set flag indicating pinched intersection

isec%pinched=iface%pinched

! first allocate local arrays

allocate(rdiff(grid%nr,grid%nlat,grid%nlong),r_interface(grid%nlat,grid%nlong)) 
nodecount=grid%nr*grid%nlat*grid%nlong
allocate(ijk_isec(nodecount),type_isec(nodecount))
allocate(nextranodes(grid%nr+1,grid%nlat+1,grid%nlong+1))
allocate(extranodes(10,grid%nr+1,grid%nlat+1,grid%nlong+1))


! evaluate height above the interface of the propagation grid nodes
 
do k=1,grid%nlong
   do j=1,grid%nlat
      r_interface(j,k)=interpolate_interface(grid%lat(j),grid%long(k),iface)
      do i=1,grid%nr
         rdiff(i,j,k)=grid%r(i)-r_interface(j,k)
      end do

   end do
end do


!--------------------------
! find the intersection nodes , the cut cells and the region borders
! note that we have allocated an extra layer of cells around the real grid
! so that we don't have to test whether the cut cell is in the grid all the time

! irg_abo(j,k) and irg_bel(j,k) contain the r-index of the first regular grid node above 
! and below the interface

allocate(isec%irg_abo(grid%nlat,grid%nlong),isec%irg_bel(grid%nlat,grid%nlong))


! initialize irg_abo and irg_bel with defaults

do k=1,grid%nlong
   do j=1,grid%nlat
      if (r_interface(j,k) > grid%r(grid%nr)) then
         isec%irg_bel(j,k) = grid%nr + 1
         isec%irg_abo(j,k) = grid%nr + 1
      endif
      if (r_interface(j,k) < grid%r(1)) then
         isec%irg_bel(j,k) = 0
         isec%irg_abo(j,k) = 0
      endif
   end do
end do




! go over all regular grid nodes and test if there is an intersection between the node and
! the next in r , lat or long

nodecount=0
nextranodes=0

do k=1,grid%nlong
   do j=1,grid%nlat
      do i=1,grid%nr

         if (abs(rdiff(i,j,k)) <= grid%tolerance)  then ! node coincides with a grid node

            ! we have a new interface node           
            nodecount=nodecount+1 

            ! store the base integer coordinates and the offset type (in r, lat or long)
            ijk_isec(nodecount)%ir=i ;ijk_isec(nodecount)%ilat=j 
            ijk_isec(nodecount)%ilong=k ;type_isec(nodecount)=0

            ! in this case the node is part of 8 adjacent cells and we want to build up a list of
            ! nodes belonging to each cell

            ! first add one to the intersection node count of the affected grid cells
            nextranodes(i+1,j+1,k+1)=nextranodes(i+1,j+1,k+1)+1
            nextranodes(i+1,j,k+1)  =nextranodes(i+1,j,k+1)+1
            nextranodes(i+1,j+1,k)  =nextranodes(i+1,j+1,k)+1
            nextranodes(i+1,j,k)    =nextranodes(i+1,j,k)+1

            nextranodes(i,j+1,k+1)=nextranodes(i,j+1,k+1)+1
            nextranodes(i,j,k+1)  =nextranodes(i,j,k+1)+1
            nextranodes(i,j+1,k)  =nextranodes(i,j+1,k)+1
            nextranodes(i,j,k)    =nextranodes(i,j,k)+1
            
      ! then store the number of the new intersection node in the list of the affected grid cells
            extranodes(nextranodes(i+1,j+1,k+1),i+1,j+1,k+1) =nodecount
            extranodes(nextranodes(i+1,j,k+1)  ,i+1,j,k+1)   =nodecount
            extranodes(nextranodes(i+1,j+1,k)  ,i+1,j+1,k)   =nodecount
            extranodes(nextranodes(i+1,j,k)    ,i+1,j,k)     =nodecount

            extranodes(nextranodes(i,j+1,k+1),i,j+1,k+1) =nodecount
            extranodes(nextranodes(i,j,k+1)  ,i,j,k+1)   =nodecount
            extranodes(nextranodes(i,j+1,k)  ,i,j+1,k)   =nodecount
            extranodes(nextranodes(i,j,k)    ,i,j,k)     =nodecount



            ! store region boundary
            isec%irg_abo(j,k)=i+1
            isec%irg_bel(j,k)=i-1


            if(diag) write(13,'(a8,4i5)') 'it=0',nodecount,ijk_isec(nodecount)%ir, &
                 ijk_isec(nodecount)%ilat,ijk_isec(nodecount)%ilong
            

         else

            if (i<grid%nr)    then   ! if next in r is not outside the grid

               if (sign(1.0_dp,rdiff(i,j,k))*rdiff(i+1,j,k) < -grid%tolerance) then   

                ! interface intersects grid between grid node and the next in r

                  ! we have a new interface node
                  nodecount=nodecount+1 

                  ! store the base integer coordinates and the offset type (in r, lat or long)

                  ijk_isec(nodecount)%ir=i ;ijk_isec(nodecount)%ilat=j 
                  ijk_isec(nodecount)%ilong=k ;type_isec(nodecount)=1

                  ! the new interface node is part of 4 grid cells,and we want to build up a list of
                  ! nodes belonging to each cell

                  ! first add one to the intersection node count of the affected grid cells
                  nextranodes(i+1,j+1,k+1)=nextranodes(i+1,j+1,k+1)+1
                  nextranodes(i+1,j,k+1)  =nextranodes(i+1,j,k+1)+1
                  nextranodes(i+1,j+1,k)  =nextranodes(i+1,j+1,k)+1
                  nextranodes(i+1,j,k)    =nextranodes(i+1,j,k)+1
                  
         ! then store the number of the new intersection node in the list of the affected grid cells
                  extranodes(nextranodes(i+1,j+1,k+1),i+1,j+1,k+1) =nodecount
                  extranodes(nextranodes(i+1,j,k+1)  ,i+1,j,k+1)   =nodecount
                  extranodes(nextranodes(i+1,j+1,k)  ,i+1,j+1,k)   =nodecount
                  extranodes(nextranodes(i+1,j,k)    ,i+1,j,k)     =nodecount

                  ! store region boundary
                  isec%irg_abo(j,k)=i+1
                  isec%irg_bel(j,k)=i


                  if(diag) write(13,'(a8,4i5)') 'it=1',nodecount,ijk_isec(nodecount)%ir, &
                       ijk_isec(nodecount)%ilat,ijk_isec(nodecount)%ilong

               endif

            endif

            if (j<grid%nlat)  then     ! if next in lat is not outside the grid

               if (sign(1.0_dp,rdiff(i,j,k))*rdiff(i,j+1,k) < -grid%tolerance) then  

             ! node between grid node and the next in lat
             
                  nodecount=nodecount+1 
                  ijk_isec(nodecount)%ir=i ;ijk_isec(nodecount)%ilat=j 
                  ijk_isec(nodecount)%ilong=k ;type_isec(nodecount)=2

                  nextranodes(i+1,j+1,k+1)=nextranodes(i+1,j+1,k+1)+1
                  nextranodes(i,j+1,k+1)=nextranodes(i,j+1,k+1)+1
                  nextranodes(i+1,j+1,k)=nextranodes(i+1,j+1,k)+1
                  nextranodes(i,j+1,k)=nextranodes(i,j+1,k)+1
                  extranodes(nextranodes(i+1,j+1,k+1),i+1,j+1,k+1) =nodecount
                  extranodes(nextranodes(i,j+1,k+1)  ,i,j+1,k+1)   =nodecount
                  extranodes(nextranodes(i+1,j+1,k)  ,i+1,j+1,k)   =nodecount
                  extranodes(nextranodes(i,j+1,k)    ,i,j+1,k)     =nodecount

                  if(diag) write(13,'(a8,4i5)') 'it=2',nodecount,ijk_isec(nodecount)%ir, &
                       ijk_isec(nodecount)%ilat,ijk_isec(nodecount)%ilong

               endif

            endif

            if (k<grid%nlong) then    ! if next in long is not outside the grid

               if (sign(1.0_dp,rdiff(i,j,k))*rdiff(i,j,k+1) < -grid%tolerance) then  

                  ! node between grid node and the next in long
                  
                  nodecount=nodecount+1 
                  ijk_isec(nodecount)%ir=i ;ijk_isec(nodecount)%ilat=j 
                  ijk_isec(nodecount)%ilong=k ;type_isec(nodecount)=3

                  nextranodes(i+1,j+1,k+1)=nextranodes(i+1,j+1,k+1)+1
                  nextranodes(i+1,j,k+1)=nextranodes(i+1,j,k+1)+1
                  nextranodes(i,j+1,k+1)=nextranodes(i,j+1,k+1)+1
                  nextranodes(i,j,k+1)=nextranodes(i,j,k+1)+1
                  extranodes(nextranodes(i+1,j+1,k+1),i+1,j+1,k+1) =nodecount
                  extranodes(nextranodes(i+1,j,k+1)  ,i+1,j,k+1)   =nodecount
                  extranodes(nextranodes(i,j+1,k+1)  ,i,j+1,k+1)   =nodecount
                  extranodes(nextranodes(i,j,k+1)    ,i,j,k+1)     =nodecount

                  if(diag) write(13,'(a8,4i5)') 'it=3',nodecount,ijk_isec(nodecount)%ir, &
                       ijk_isec(nodecount)%ilat,ijk_isec(nodecount)%ilong

               endif

            endif

         endif
      end do
   end do
end do

print *,'intersection ',isec%id, ':',nodecount,' nodes found'

if (nodecount == 0) then
   isec%nnode=nodecount
!   print *,'no intersection points found,intersection',isec%id,' not created'
   deallocate(rdiff,r_interface,ijk_isec,type_isec,extranodes,nextranodes)
   deallocate(isec%irg_abo,isec%irg_bel)
   return
endif

!print *,isec%id,'irgbel,irgabo',isec%irg_bel(11,71),isec%irg_abo(11,71)


! allocate the arrays required to store the intersection

isec%nnode=nodecount      ! the number of nodes in the intersection

! the # of cells cut by the intersection
isec%n_ccells=count(nextranodes(2:grid%nr,2:grid%nlat,2:grid%nlong)>0) 

allocate(isec%r(nodecount),isec%lat(nodecount),isec%long(nodecount),isec%coslat(nodecount))
allocate(isec%intype(nodecount))
allocate(isec%ccell_from_inode(8,isec%nnode))
allocate(isec%ccells(isec%n_ccells))
allocate(isec%n_inodes(isec%n_ccells))
allocate(isec%inodes(maxval(nextranodes),isec%n_ccells))
allocate(isec%time_gradient(3,nodecount))
allocate(isec%normal(3,nodecount))


! store the 3d-indices of the cut cells and all interface nodes that belong to them

isec%n_inodes=0
n=0
do k=2,grid%nlong
   do j=2,grid%nlat
      do i=2,grid%nr

         if (nextranodes(i,j,k)>0) then  ! this cell is a cut cell

            n=n+1

            isec%ccells(n)%ir=i-1        ! integer coordinates of this cut cell
            isec%ccells(n)%ilat=j-1
            isec%ccells(n)%ilong=k-1

            ! store intersection nodes of this cut cell
            isec%n_inodes(n)=nextranodes(i,j,k)  ! store intersection nodes of this cut cell
            isec%inodes(1:isec%n_inodes(n),n)=extranodes(1:isec%n_inodes(n),i,j,k)

            ! construct the pointer from the corresponding main grid cell to this cut cell
            if (.not.associated(grid%ccind_from_3dc(i-1,j-1,k-1)%p)) then
               allocate(grid%ccind_from_3dc(i-1,j-1,k-1)%p(n_interfaces))
               grid%ccind_from_3dc(i-1,j-1,k-1)%p = 0
            endif
            grid%ccind_from_3dc(i-1,j-1,k-1)%p(isec%iface_id)=n

         endif

      end do
   end do
end do


! construct the inverse pointers from inodes to the cut cell they are in

allocate(counter(nodecount))
counter=0
isec%ccell_from_inode=0

do n=1,isec%n_ccells
   do m=1,isec%n_inodes(n)
      counter(isec%inodes(m,n))=counter(isec%inodes(m,n))+1
      if(counter(isec%inodes(m,n))>8) stop 'find_intersection:  counter > 8'
      isec%ccell_from_inode(counter(isec%inodes(m,n)),isec%inodes(m,n))=n
   end do
end do

deallocate(counter)


! find the coordinates of the intersection nodes

do n=1,isec%nnode

   select case (type_isec(n))

   case (0)          ! node coincides with a grid node

      isec%r(n)    = grid%r(ijk_isec(n)%ir)
      isec%lat(n)  = grid%lat(ijk_isec(n)%ilat)
      isec%coslat(n)  = grid%coslat(ijk_isec(n)%ilat)
      isec%long(n) = grid%long(ijk_isec(n)%ilong)
      isec%intype(n)= type_isec(n)

   case (1)          ! node between grid node and the next in r

      isec%r(n)    = r_interface(ijk_isec(n)%ilat,ijk_isec(n)%ilong)
      isec%lat(n)  = grid%lat(ijk_isec(n)%ilat)
      isec%coslat(n)  = grid%coslat(ijk_isec(n)%ilat)
      isec%long(n) = grid%long(ijk_isec(n)%ilong)
      isec%intype(n)= type_isec(n)

   case (2)          ! node between grid node and the next in lat

      isec%r(n)    = grid%r(ijk_isec(n)%ir)
      isec%lat(n)  = find_zero(ijk_isec(n),type_isec(n))
      isec%coslat(n)=cos(isec%lat(n))
      isec%long(n) = grid%long(ijk_isec(n)%ilong)
      isec%intype(n)= type_isec(n)

   case (3)          ! node between grid node and the next in long

      isec%r(n)    = grid%r(ijk_isec(n)%ir)
      isec%lat(n)  = grid%lat(ijk_isec(n)%ilat)
      isec%coslat(n)  = grid%coslat(ijk_isec(n)%ilat)
      isec%long(n) = find_zero(ijk_isec(n),type_isec(n))
      isec%intype(n)= type_isec(n)


   case default 
      stop 'illegal intersection type in subroutine find_intersection'

   end select

end do

! deallocate all the local temporary arrays

deallocate(rdiff,r_interface,ijk_isec,type_isec,extranodes,nextranodes)


! calculate the surface normals

do n=1,isec%nnode
   call interface_normal(isec%lat(n),isec%long(n),iface,isec%normal(1,n), &
        isec%normal(2,n),isec%normal(3,n),rint)
end do


contains

  function find_zero(ijk,itype)
  ! bisection root finding subroutine from NR

    use mod_3dfm
    implicit none

    real(kind=dp) :: find_zero

    ! function arguments
    type(Tinteger_coordinates)     :: ijk     ! base integer coordinates of the root
    integer     :: itype   ! root is offset from base coordinates in lat (itype=2) or long(itype=3)

    real(kind=dp) :: f,fmid,x1,x2,xmid,rtbis,dx 
    integer :: nn



    select case (itype)

        case(2)    ! intersection node lies in positive lat direction
           
           x1=grid%lat(ijk%ilat)
           x2=grid%lat(ijk%ilat+1)
           f=rdiff(ijk%ir,ijk%ilat,ijk%ilong)
           fmid=rdiff(ijk%ir,ijk%ilat+1,ijk%ilong)
           if (f*fmid >= 0.0_dp) stop 'root not bracketed in find_intersection'
           if (f < 0.0_dp) then
              rtbis=x1
              dx=x2-x1
           else
              rtbis=x2
              dx=x1-x2
           endif
           do nn=1,10
              dx=dx*0.5_dp
              xmid=rtbis+dx
              fmid=grid%r(ijk%ir)-interpolate_interface(xmid,grid%long(ijk%ilong),iface)
              if (fmid <= 0.0_dp) rtbis=xmid
           end do


        case(3)    ! intersection node lies in positive long direction

           x1=grid%long(ijk%ilong)
           x2=grid%long(ijk%ilong+1)
           f=rdiff(ijk%ir,ijk%ilat,ijk%ilong)
           fmid=rdiff(ijk%ir,ijk%ilat,ijk%ilong+1)
           if (f*fmid >= 0.0_dp) stop 'root not bracketed in find_intersection'
           if (f < 0.0_dp) then
              rtbis=x1
              dx=x2-x1
           else
              rtbis=x2
              dx=x1-x2
           endif
           do nn=1,10
              dx=dx*0.5_dp
              xmid=rtbis+dx
              fmid=grid%r(ijk%ir)-interpolate_interface(grid%lat(ijk%ilat),xmid,iface)
              if (fmid <= 0.0_dp) rtbis=xmid
           end do


        case default
           stop 'illegal intersection type in function find_zero'

      end select

      find_zero=rtbis

   end function find_zero

end subroutine find_intersection

!***********************************************************************************************
! This subroutine constructs the list of nodes (regular, top/bottom intersection) 
! belonging to a region

subroutine define_region(reg,itop,ibot,grid)
use mod_3dfm_nointerfaces
implicit none

type(Tregion),target     :: reg
type(Tintersection)      :: itop,ibot
type(Tpropagation_grid),target :: grid

integer :: m,i,j,k,istart,iend


if (.not.associated(itop%r) .and. itop%nnode /= 0) stop 'define regions: intersections not allocated'

! set the velocity grid to be used for this region

reg%ivgrid=itop%iface_id

! pointer to the grid on which this region is defined

reg%grid => grid


! assign the proper region # to the nodes of the main propagation grid  
! itop or ibot may be above or below the grid and have no intersection nodes

if (itop%nnode > 0 .and. ibot%nnode > 0) then      
   do k=1,grid%nlong
      do j=1,grid%nlat
         istart=max(1,ibot%irg_abo(j,k))
         iend=min(grid%nr,itop%irg_bel(j,k))
         if (iend >= istart) then
            do i=istart,iend
               grid%node_region(i,j,k) = reg%id
            end do
         endif
      end do
   end do
endif
if (itop%nnode > 0 .and. ibot%nnode == 0) then      
   do k=1,grid%nlong
      do j=1,grid%nlat
         do i=1,min(grid%nr,itop%irg_bel(j,k))
            grid%node_region(i,j,k) = reg%id
         end do
      end do
   end do
endif
if (itop%nnode == 0 .and. ibot%nnode > 0) then      
   do k=1,grid%nlong
      do j=1,grid%nlat
         do i=max(1,ibot%irg_abo(j,k)),grid%nr
            grid%node_region(i,j,k) = reg%id
         end do
      end do
   end do
endif
if (itop%nnode == 0 .and. ibot%nnode == 0) then      
   do k=1,grid%nlong
      do j=1,grid%nlat
         do i=1,grid%nr
            grid%node_region(i,j,k) = reg%id
         end do
      end do
   end do
endif


! print *,'def reg: grid node region id assigned'

! derive and store useful information about the regions

! register the region with bounding intersections
   itop%regbel => reg 
   ibot%regabo => reg 


! # of grid nodes in this region
   reg%ngnode=count(grid%node_region == reg%id)  ! regular grid nodes only
   reg%nnode=reg%ngnode + itop%nnode + ibot%nnode ! total

!   print *,'reg def: nnode = ',reg%nnode


! make a 1-D list of all nodes (grid + intersection) nodes in this region

!  array reg%node contains pointers from the list to the grid/interface nodes

   allocate(reg%node(reg%nnode))


! these arrays contain pointers from intersection arrays to the regional node lists above and below

   if (itop%nnode > 0) allocate(itop%rbel_node_id(itop%nnode))
   if (ibot%nnode > 0) allocate(ibot%rabo_node_id(ibot%nnode))



! storage for regional node coordinates
   allocate(reg%r(reg%nnode))
   allocate(reg%lat(reg%nnode))
   allocate(reg%coslat(reg%nnode))
   allocate(reg%long(reg%nnode))

 !  print *,'reg def : pos arrrays allocated'

! start constructing the 1-D pointer arrays
   m=0

! add the regular grid points
   do k=1,grid%nlong
      do j=1,grid%nlat
         do i=1,grid%nr
            if (grid%node_region(i,j,k) == reg%id) then
               m=m+1
               reg%node(m)%i1=i 
               reg%node(m)%i2=j 
               reg%node(m)%i3=k 
         ! grid%rnode_id contains the index of the node in its regional node list
               grid%rnode_id(i,j,k)=m  
               reg%r(m) = grid%r(i)
               reg%lat(m) = grid%lat(j)
               reg%coslat(m) = grid%coslat(j)
               reg%long(m) = grid%long(k)
            endif
         end do
      end do
   end do

! print *,'reg def regular nodes initialized'

! add the top bounding intersection nodes if they exist 
   if (itop%nnode > 0) then
      do i=1,itop%nnode
         m=m+1
         reg%node(m)%i1=0 
         reg%node(m)%i2=itop%id 
         reg%node(m)%i3=i
         itop%rbel_node_id(i) = m
         reg%r(m) = itop%r(i)
         reg%lat(m) = itop%lat(i)
         reg%coslat(m) = itop%coslat(i)
         reg%long(m) = itop%long(i)
      end do
   endif

! print *,'reg def top nodes initialized'

! add the bottom bounding intersection nodes if they exist
   if (ibot%nnode > 0) then
      do i=1,ibot%nnode
         m=m+1
         reg%node(m)%i1=0 
         reg%node(m)%i2=ibot%id 
         reg%node(m)%i3=i
         ibot%rabo_node_id(i) = m
         reg%r(m) = ibot%r(i)
         reg%lat(m) = ibot%lat(i)
         reg%coslat(m) = ibot%coslat(i)
         reg%long(m) = ibot%long(i)
      end do
   endif

   if (m /= reg%nnode) stop 'define region: node count mismatch'

!   print *,'reg def bot nodes initialized'

return

end subroutine define_region

!*************************************************************************************************
! calculates the velocities on the nodes of a propagation grid with cubic spline interpolation
subroutine velocities_to_grid(grid)
use mod_3dfm
implicit none

type(Tpropagation_grid)  :: grid

integer  :: i,j,k,ivg,m

real(kind=dp)  :: interpolate_velocity


allocate(grid%velocity(grid%nr,grid%nlat,grid%nlong,n_vtypes))

! velocities on regular grid

if (grid%is_main_grid) then

do m=1,n_vtypes
   do k=1,grid%nlong
      do j=1,grid%nlat
         do i=1,grid%nr
            if (grid%node_region(i,j,k)>0) then
               ivg = region(grid%node_region(i,j,k))%ivgrid
               grid%velocity(i,j,k,m) = interpolate_velocity(grid%r(i),grid%lat(j), &
                    grid%long(k),vgrid(ivg,m))
            endif
         end do
      end do
   end do
end do

endif

if (grid%is_source_grid) then

do m=1,n_vtypes
   do k=1,grid%nlong
      do j=1,grid%nlat
         do i=1,grid%nr
            if (grid%node_region(i,j,k)>0) then
               ivg = sregion(grid%node_region(i,j,k))%ivgrid
               grid%velocity(i,j,k,m) = interpolate_velocity(grid%r(i),grid%lat(j), &
                    grid%long(k),vgrid(ivg,m))
            endif
         end do
      end do
   end do
end do

endif

return

end subroutine velocities_to_grid

!*************************************************************************************************
! calculates the velocities on the nodes of an intersection with cubic spline interpolation
subroutine velocities_to_intersection(isec)
use mod_3dfm
implicit none

type(Tintersection)  :: isec
integer  :: i,j,n,m,vg_abo,vg_bel,isig_abo,isig_bel,iface_id

real(kind=dp)  :: interpolate_velocity,interpolate_interface,h,vel


! velocities on intersection nodes

   if (associated(isec%regabo)) then   

! if there is a defined region above the intersection,get velocities on the top

      allocate(isec%vel_top(isec%nnode,n_vtypes)) 

      vg_abo = isec%regabo%ivgrid

      do m=1,n_vtypes
         do i=1,isec%nnode
            isec%vel_top(i,m)=interpolate_velocity(isec%r(i),isec%lat(i),isec%long(i),vgrid(vg_abo,m))
         end do
      end do

   endif

   if (associated(isec%regbel)) then   

! if there is a defined region below the intersection,get velocities on the bottom

      allocate(isec%vel_bot(isec%nnode,n_vtypes))

      vg_bel = isec%regbel%ivgrid

      do m=1,n_vtypes
         do i=1,isec%nnode
            isec%vel_bot(i,m)=interpolate_velocity(isec%r(i),isec%lat(i),isec%long(i),vgrid(vg_bel,m))
         end do
      end do
   endif



! adjust velocities at the intersection nodes where the interfaces are pinched


   if (isec%pinched) then

      m=0
      iface_id = isec%iface_id


  ! find the interfaces significantly above or below the current one at the node positions

nodeloop:      do i=1,isec%nnode

      ! above
         isig_abo = 0         
         if (iface_id > 1) then
            do j=iface_id-1,1,-1
               h = interpolate_interface(isec%lat(i),isec%long(i),intrface(j))
               if (h-isec%r(i) > pgrid%tolerance) then
                  isig_abo = j
                  exit
               end if
            end do
         end if

      ! below
         isig_bel = n_interfaces+1         
         if (iface_id < n_interfaces) then
            do j=iface_id+1,n_interfaces
               h = interpolate_interface(isec%lat(i),isec%long(i),intrface(j))
               if (isec%r(i)-h > pgrid%tolerance) then
                  isig_bel = j
                  exit
               end if
            end do
         end if


      ! go to next node if interfaces not pinched at this position

         if (isig_abo == iface_id-1 .and. isig_bel == iface_id+1) cycle nodeloop



      ! now adjust the velocity values at the node 


         if (isig_abo == 0) then      ! the special case that the interface coincides with the surface

            ! set top and bottom velocity to that from the first significant region below
            ! this ensures the wave passes without modification

            vg_bel=region(isig_bel-1)%ivgrid

            do n=1,n_vtypes
               vel = interpolate_velocity(isec%r(i),isec%lat(i),isec%long(i),vgrid(vg_bel,n))

               if (iface_id > 1) then    ! vel_top of intersection(1) does not exist
                  isec%vel_top(i,n)=vel
               endif
               isec%vel_bot(i,n)=vel

            end do

            m=m+1
            cycle nodeloop

         endif


         if (isig_bel == n_interfaces+1) then      

       ! the special case that the interface coincides with the bottom

           ! set top and bottom velocity to that from the first significant region above
           ! this ensures the wave passes without modification

            vg_abo=region(isig_abo)%ivgrid

            do n=1,n_vtypes
               vel = interpolate_velocity(isec%r(i),isec%lat(i),isec%long(i),vgrid(vg_abo,n))

               if (iface_id < n_interfaces) then 

           ! vel_bot of intersection(n_interfaces) does not exist

                  isec%vel_bot(i,n)=vel
               endif
               isec%vel_top(i,n)=vel

           end do

            m=m+1
            cycle nodeloop

         endif

       ! the normal case: make only the lowest of the locally pinched interfaces significant

         vg_abo=region(isig_abo)%ivgrid
         do n=1,n_vtypes
            vel = interpolate_velocity(isec%r(i),isec%lat(i),isec%long(i),vgrid(vg_abo,n))
            isec%vel_top(i,n) = vel 
            if (isig_bel /= iface_id+1) isec%vel_bot(i,n) = vel
         end do
         m=m+1

      end do nodeloop ! loop over intersection nodes

      print '(a50,2i7)','corrected velocities at pinched nodes of iface ',isec%iface_id,m

   endif  ! if interface is pinched


   return

end subroutine velocities_to_intersection

!*************************************************************************************************
! copies the velocities stored in grid and intersections to the corresponding nodes of a region
subroutine velocities_to_region(reg,grid)
use mod_3dfm
implicit none

type(Tpropagation_grid)  :: grid
type(Tregion)  :: reg
integer  :: n,m

! assign velocities in regions

   allocate(reg%velocity(reg%nnode,n_vtypes))

   do n=1,n_vtypes
      do m=1,reg%nnode

         if (reg%node(m)%i1 /= 0) then

            reg%velocity(m,n)= grid%velocity(reg%node(m)%i1,reg%node(m)%i2,reg%node(m)%i3,n)

         else

            if (reg%node(m)%i2 == reg%itop%id) then
               reg%velocity(m,n)= reg%itop%vel_bot(reg%node(m)%i3,n)
            else
               if (reg%node(m)%i2 /= reg%ibot%id)  stop 'transfer velocities: ordering problem'
               reg%velocity(m,n)= reg%ibot%vel_top(reg%node(m)%i3,n) 
            endif

         endif
      end do
   end do


! tag the nodes on the velocity grid that actually influence the velocity field in this region
! since the velocity grids are square and in the simplest case cover the entire propagation grid
! this can be a small subset only

   call tag_active_vgrid_nodes(reg)


return

end subroutine velocities_to_region



!*****************************************************************************************************
!----------------------------------------------------------------------------------------------
! finds out where a source lies relative to propagation grid and intersections
! and stores this information inside the derived type Tsource associated with it

subroutine initialize_source(s,grid)
  use mod_3dfm
  implicit none

  integer       :: n,i,j,k
  type(Tsource) :: s 
  type(Tpropagation_grid) :: grid
  real(kind=dp) :: dist,dist_below,interpolate_interface,h
  logical,dimension(:),allocatable       :: source_on_interface


!  print *,'entering initialize_source'

! test for source in grid

  if (s%lat > pgrid%lat(pgrid%nlat)+pgrid%tolerance/s%r ) stop 'ERROR: source position beyond maximum lat'
  if (s%lat < pgrid%lat(1) -pgrid%tolerance/s%r) stop 'ERROR:source position beyond minimum lat'
  if (s%long > pgrid%long(pgrid%nlong)+pgrid%tolerance/s%r ) stop 'ERROR:source position beyond maximum long'
  if (s%long < pgrid%long(1)-pgrid%tolerance/s%r ) stop 'ERROR:source position beyond minimum long'
  if (s%r > interpolate_interface(s%lat,s%long,intrface(1))+pgrid%tolerance) &
       stop 'ERROR:source above surface'
  if (s%r < interpolate_interface(s%lat,s%long,intrface(n_interfaces))-pgrid%tolerance ) &
       stop 'ERROR:source below lowest interface'

!  print*, 'after test',s%r,interpolate_interface(s%lat,s%long,intrface(n_interfaces))-pgrid%tolerance

! determine grid cell in which the source is located

  s%ir    =  floor((s%r - grid%r0)/grid%dr0 + 1)
  s%ilat  =  floor((s%lat - grid%lat0)/grid%dlat0 + 1)
  s%ilong =  floor((s%long - grid%long0)/grid%dlong0 + 1)

! correct if source lies exactly on grid boundary 
 
  s%ir = min(s%ir,grid%nr-1)
  s%ilat = min(s%ilat,grid%nlat-1)
  s%ilong = min(s%ilong,grid%nlong-1)

  s%on_grid=.false.
  s%on_interface=.false.
  s%n_cnode = 0

! allocate arrays that will contain the indices of the source time fields in the time field array
  allocate(s%first_tf_up(n_vtypes),s%first_tf_down(n_vtypes))


! test where the source lies


! first test if the source lies exactly on any interface

  allocate(source_on_interface(n_interfaces))
  source_on_interface=.false.

  dist_below=1.0_dp
  do n=n_interfaces,1,-1

     dist = s%r-interpolate_interface(s%lat,s%long,intrface(n))

     source_on_interface(n) = abs(dist) < 2.0_dp*grid%tolerance 
     if (source_on_interface(n))   print *,'source lies exactly on interface',n

     if (dist < 0.0_dp .and. dist_below > 0.0_dp) s%region_id = n

     dist_below=dist

  enddo

  s%on_interface = count(source_on_interface(1:n_interfaces)) > 0
  s%on_pinched_interface = count(source_on_interface(1:n_interfaces)) > 1


  if (s%on_interface) then

     s%topint_id = 0 
     s%botint_id = n_interfaces+1

     do i=1,n_interfaces

        h = interpolate_interface(s%lat,s%long,intrface(i))
        if (h-s%r > pgrid%tolerance) s%topint_id = i

     end do

     do i=n_interfaces,1,-1

        h = interpolate_interface(s%lat,s%long,intrface(i))
        if (s%r-h > pgrid%tolerance) s%botint_id = i

     end do

     s%topint_id=s%topint_id+1
     s%botint_id=s%botint_id-1

     s%topreg_id=s%topint_id-1
     s%botreg_id=s%botint_id

     if (s%topint_id == 1 .or. s%botint_id == n_interfaces) then
        s%n_tf_init = 1
     else
        s%n_tf_init = 2
     endif

!     print *,'istest',s%topint_id,s%botint_id,s%topreg_id,s%botreg_id
!     stop
  endif

  deallocate(source_on_interface)


! test if the source lies exactly on a regular grid node or not

  if (.not. s%on_interface) then

     do i=0,1
        do j=0,1
           do k=0,1

              dist=sqrt((s%r-grid%r(s%ir+i))**2 + (s%lat-grid%lat(s%ilat+j))**2 +  &
                   (s%long-grid%long(s%ilong+k))**2 ) 
              if (dist < grid%tolerance) s%on_grid=.true.

           end do
        end do
     end do
 
     s%n_tf_init = 1
             
  endif  ! if not on interface



  return

end subroutine initialize_source

!----------------------------------------------------------------------------------------
! this is the basic unit used to find a fast marching solution in a region
! arrival times and gradients are copied from the starting intersection
! to the region, node status is set appropriately and then the region is
! passed on to the fast marching routine propagate
! also tests for and flags turning rays, i.e. nodes on the starting interface 
! that receive a new time during fast marching


subroutine sweep_region_from_interface(reg,istart_in,vtype,s)
  use mod_3dfm_nointerfaces
  implicit none

  integer                                :: vtype
  integer                                :: i,j,i1,i2,i3,n_turning
  type(Tintersection),pointer            :: itop,ibot,istart
  type(Tintersection),target             :: istart_in
  type(Tregion)                          :: reg
  type(Tsource)                          :: s
  logical,dimension(:),allocatable       :: received_turning_ray

!----------------------------------------------------------
! prepare for propagation in a region


  if (reg%id > n_regions .or. reg%id < 1) stop 'invalid region requested'
  if (istart_in%iface_id > n_interfaces .or. istart_in%iface_id < 1) &
       stop ' invalid starting intersection'
  if (.not.(istart_in%id == reg%itop%id .or. istart_in%id == reg%ibot%id ))  &
       stop ' requested starting intersection and region incompatible'


! local names for  starting, top and bottom intersections

  istart => istart_in
  itop => reg%itop
  ibot => reg%ibot

  allocate(reg%arrivaltime(reg%nnode),reg%time_gradient(3,reg%nnode))
  allocate(reg%node_status(reg%nnode))
  allocate(received_turning_ray(istart%nnode))
  received_turning_ray=.false.

  reg%node_status=-1

  if (.not. associated(istart%arrivaltime))  &
       stop 'requested starting interface does not have arrival times'

  allocate(istart%starttime(istart%nnode))


! set start time, time gradient and node status narrow band in regional 
! nodes that belong to the starting interface

  do i=1,istart%nnode
     istart%starttime(i)=istart%arrivaltime(i)
     if (istart%id == ibot%id) then
        j=istart%rabo_node_id(i)
     else
        j=istart%rbel_node_id(i)
     endif
     reg%arrivaltime(j)= istart%arrivaltime(i)
     reg%time_gradient(1:3,j)=istart%time_gradient(1:3,i)
     reg%node_status(j)=1
  end do

 ! print *,'region',reg%id,' starting times set'


!-------------------------------------------------------------------------------------
! do the fast marching sweep

  call propagate(reg,vtype)

  print *,'propagation through region',reg%id,' finished'


! transfer regional travel times to interfaces and regular grid
! and check for turning rays

  n_turning=0
  if (itop%nnode > 0 .and. .not. associated(itop%arrivaltime)) allocate(itop%arrivaltime(itop%nnode))
  if (ibot%nnode > 0 .and. .not. associated(ibot%arrivaltime)) allocate(ibot%arrivaltime(ibot%nnode))

  do i=1,reg%nnode

     i1 = reg%node(i)%i1 ; i2 = reg%node(i)%i2 ; i3 = reg%node(i)%i3

     if (i1 /= 0) then

!        grid%arrivaltime(i1,i2,i3)=reg%arrivaltime(i)

     else

        if (i2 /= itop%id .and. i2 /= ibot%id)  &
             stop 'intersection id mismatch in sweep_region_from_interface' 

        if (i2 == itop%id) then
           itop%arrivaltime(i3) = reg%arrivaltime(i)
           itop%time_gradient(1:3,i3) = reg%time_gradient(1:3,i)        
        else
           ibot%arrivaltime(i3) = reg%arrivaltime(i)
           ibot%time_gradient(1:3,i3) = reg%time_gradient(1:3,i)
        endif

        if (i2 == istart%id) then
           if (istart%arrivaltime(i3) < istart%starttime(i3)) then 
              n_turning=n_turning+1
              received_turning_ray(i3) = .true.
           endif
        endif

     endif

  end do
      
  if (n_turning > 0) print *,n_turning,' starting intersection nodes received a turning ray'


  ! transfer time field to the array of saved time fields

  s%n_time_fields=s%n_time_fields+1
  allocate(s%time_field(s%n_time_fields)%arrivaltime(reg%nnode))
  s%time_field(s%n_time_fields)%arrivaltime=reg%arrivaltime
  allocate(s%time_field(s%n_time_fields)%time_gradient(3,reg%nnode))
  s%time_field(s%n_time_fields)%time_gradient=reg%time_gradient


  ! if  turning rays were present, identify the intersection nodes that received them

  if (n_turning > 0) then
     s%time_field(s%n_time_fields)%turning_rays_present = .true.
     allocate(s%time_field(s%n_time_fields)%received_turning_ray(istart%nnode))
     s%time_field(s%n_time_fields)%received_turning_ray = received_turning_ray
  endif

  deallocate(reg%arrivaltime,reg%time_gradient,reg%node_status,istart%starttime)
  deallocate(received_turning_ray)

end subroutine sweep_region_from_interface

!----------------------------------------------------------------------------------------
! A simple version of the routine above used only for regions on the refined source grid
subroutine sweep_sregion_from_interface(reg,istart_in,vtype)
  use mod_3dfm_nointerfaces
  implicit none

  integer                                :: vtype
  integer                                :: i,j
  type(Tintersection),pointer            :: itop,ibot,istart
  type(Tintersection),target             :: istart_in
  type(Tregion)                          :: reg

!----------------------------------------------------------
! prepare for propagation in a region


  if (reg%id > n_regions .or. reg%id < 1) stop 'invalid region requested'
  if (istart_in%iface_id > n_interfaces .or. istart_in%iface_id < 1) &
       stop ' invalid starting intersection'
  if (.not.(istart_in%id == reg%itop%id .or. istart_in%id == reg%ibot%id ))  &
       stop ' requested starting intersection and region incompatible'


! local names for  starting, top and bottom intersections

  istart => istart_in
  itop => reg%itop
  ibot => reg%ibot

  allocate(reg%arrivaltime(reg%nnode),reg%time_gradient(3,reg%nnode))
  allocate(reg%node_status(reg%nnode))
  reg%node_status=-1

  if (.not. associated(istart%arrivaltime))  &
       stop 'requested starting interface does not have arrival times'

  allocate(istart%starttime(istart%nnode))


! set start time, time gradient and node status narrow band in regional 
! nodes that belong to the starting interface
  do i=1,istart%nnode
     istart%starttime(i)=istart%arrivaltime(i)
     if (istart%id == ibot%id) then
        j=istart%rabo_node_id(i)
     else
        j=istart%rbel_node_id(i)
     endif
     reg%arrivaltime(j)= istart%arrivaltime(i)
     reg%time_gradient(1:3,j)=istart%time_gradient(1:3,i)
     reg%node_status(j)=1
  end do

 ! print *,'sregion',reg%id,' starting times set'


!-------------------------------------------------------------------------------------
! do the fast marching sweep

  call propagate(reg,vtype)

  print *,'propagation through sregion',reg%id,' finished'

end subroutine sweep_sregion_from_interface


!----------------------------------------------------------------------------------------

subroutine sweep_region_from_source(reg,s,vtype)
  use mod_3dfm
  implicit none

  type(Tsource)   :: s
  type(Tpropagation_grid),pointer :: grid
  integer         :: vtype
  integer         :: m,j,i1,i2,i3
  type(Tintersection),pointer :: itop,ibot
  type(Tregion)   :: reg

  real(kind=dp)   :: dist,vel_av,vel_source,interpolate_velocity,dtr


  grid => reg%grid

! validate the arguments

  if (.not.s%on_interface) then
     if (reg%id /= s%region_id) stop 'region and source inconsistent'
  endif

  dtr=acos(-1.0_dp)/180.0_dp

! some local variable names to make expressions shorter

  itop => reg%itop
  ibot => reg%ibot

  vel_source=interpolate_velocity(s%r,s%lat,s%long,vgrid(reg%ivgrid,vtype))

! allocate space for regional arrays that do not need to be saved for later

  allocate(reg%arrivaltime(reg%nnode),reg%time_gradient(3,reg%nnode),reg%node_status(reg%nnode))

! set all regional points to FAR
  reg%node_status=-1
  reg%arrivaltime=huge_time


! set start time and node status NARROW BAND at the source

  if (s%n_cnode == 1) then   ! source lies exactly on a node, only 1 starting point with time 0

     i1=s%cnode(1)%i1  ;i2=s%cnode(1)%i2  ; i3=s%cnode(1)%i3  

     if (i1 /= 0) then  ! node is a regular grid point

        j= grid%rnode_id(i1,i2,i3)
        reg%arrivaltime(j)= 0.0_dp
        reg%time_gradient(1:3,j)= 0.0_dp
        reg%node_status(j)=1 
       
!        print *,'time 0 assigned to',j,i1,i2,i3
!        print *,reg%r(j),reg%lat(j)/dtr,reg%long(j)/dtr


     else                          ! node is an interface node

        if (i2 == ibot%id) then
           j=ibot%rabo_node_id(i3)
        else
           j=itop%rbel_node_id(i3)
        endif
        reg%arrivaltime(j)= 0.0_dp
        reg%time_gradient(1:3,j)= 0.0_dp
        reg%node_status(j)= 1

        print *,'source node',j,'set to zero time'
!        print *,reg%r(j),reg%lat(j)/dtr,reg%long(j)/dtr

     endif

  else                     ! source does not lie on a node, several starting points

     do m=1,s%n_cnode

     i1=s%cnode(m)%i1  ;i2=s%cnode(m)%i2  ; i3=s%cnode(m)%i3  

     if (i1 /= 0) then  ! node is a regular grid point

        j= grid%rnode_id(i1,i2,i3)
        vel_av=0.5_dp*(grid%velocity(i1,i2,i3,vtype)+vel_source)
        dist=sqrt((s%r-grid%r(i1))**2 + (s%r*(s%lat-grid%lat(i2)))**2 + &
             (s%r*s%coslat*(s%long-grid%long(i3)))**2 ) 
        reg%arrivaltime(j)= dist/vel_av
        reg%time_gradient(1,j)= (grid%r(i1)-s%r)/(dist*grid%velocity(i1,i2,i3,vtype))
        reg%time_gradient(2,j)= s%r*(grid%lat(i2)-s%lat)/(dist*grid%velocity(i1,i2,i3,vtype))
        reg%time_gradient(3,j)= s%r*s%coslat*(grid%long(i3)-s%long)/ &
             (dist*grid%velocity(i1,i2,i3,vtype))

        reg%node_status(j)=1   
     
     else                          ! node is an interface node

        if (i2 == ibot%id) then   ! node lies on the bottom interface of the region

           j=ibot%rabo_node_id(i3)
           dist=sqrt((s%r-ibot%r(i3))**2 + (s%r*(s%lat-ibot%lat(i3)))**2 + &
                (s%r*s%coslat*(s%long-ibot%long(i3)))**2 )
           reg%arrivaltime(j)= dist/ibot%vel_top(i3,vtype)
           reg%time_gradient(1,j)= (ibot%r(i3)-s%r)/(dist*ibot%vel_top(i3,vtype))
           reg%time_gradient(2,j)= s%r*(ibot%lat(i3)-s%lat)/(dist*ibot%vel_top(i3,vtype))
           reg%time_gradient(3,j)= s%r*s%coslat*(ibot%long(i3)-s%long)/(dist*ibot%vel_top(i3,vtype))

        else                             ! node lies on the top interface of the region

           j=itop%rbel_node_id(i3)
           dist=sqrt((s%r-itop%r(i3))**2 + (s%r*(s%lat-itop%lat(i3)))**2 + &
                (s%r*s%coslat*(s%long-itop%long(i3)))**2 )
           reg%arrivaltime(j)= dist/itop%vel_bot(i3,vtype)
           reg%time_gradient(1,j)= (itop%r(i3)-s%r)/(dist*itop%vel_bot(i3,vtype))
           reg%time_gradient(2,j)= s%r*(itop%lat(i3)-s%lat)/(dist*itop%vel_bot(i3,vtype))
           reg%time_gradient(3,j)= s%r*s%coslat*(itop%long(i3)-s%long)/(dist*itop%vel_bot(i3,vtype))

        endif

        reg%node_status(j)=1

     endif

     end do

  endif

!  print *,'region',reg%id,' starting times set'

!  print *,'calling propagate in sweep region from source'


!-------------------------------------------------------------------------------------
! do the fast marching sweep

  call propagate(reg,vtype)

  print *,'propagation through region',reg%id,' finished'

  return

end subroutine sweep_region_from_source


!*****************************************************************************************************
!----------------------------------------------------------------------------------------------------

! This subroutine initializes the source in the refined source grid and its intersections

subroutine initialize_refined_source(s,sc,grid,reg,itop,ibot)

  use mod_3dfm_nointerfaces
  implicit none

  integer       :: n,m,i,j,k,is,js,ks,iface
  type(Tsource) :: s,sc
  type(Tpropagation_grid) :: grid
  type(Tintersection),target   :: itop,ibot
  type(Tregion)         :: reg
  type(Tintersection),pointer   :: isec
  real(kind=dp) :: dist



! copy position from coarse source
  s%r=sc%r
  s%lat=sc%lat
  s%long=sc%long
  s%coslat=sc%coslat
  s%region_id=reg%id

! determine grid cell in which the source is located

  s%ir    =  floor((s%r - grid%r0)/grid%dr0 + 1)
  s%ilat  =  floor((s%lat - grid%lat0)/grid%dlat0 + 1)
  s%ilong =  floor((s%long - grid%long0)/grid%dlong0 + 1)

! correct if source lies exactly on grid boundary 
 
  s%ir = min(s%ir,grid%nr-1)
  s%ilat = min(s%ilat,grid%nlat-1)
  s%ilong = min(s%ilong,grid%nlong-1)

  s%on_grid=.false.
  s%on_interface=sc%on_interface
  s%n_cnode = 0

  allocate(s%first_tf_up(n_vtypes),s%first_tf_down(n_vtypes))

!-----------------------------------------


! test where the source lies

  if (s%on_interface) then

     if (sc%topint_id == ibot%id ) then
        s%interface_id = ibot%id
        isec => ibot
     else if (sc%botint_id == itop%id) then
        s%interface_id = itop%id
        isec => itop
     else
        stop 'illegal source interface in initialize refined source'
     endif


 !     print *,'isec%id',n,isec%id,isec%iface_id

     iface=isec%iface_id


        ! test if the source lies exactly on an intersection node

     do m=1,isec%nnode

        dist =sqrt((s%r-isec%r(m))**2 +(s%r*(s%lat-isec%lat(m)))**2 + &
             (s%r*s%coslat*(s%long-isec%long(m)))**2 )

        if (abs(dist) < grid%tolerance) then  ! source lies exactly on an interface node

           ! assign time value 0 to the source node

           s%n_cnode = s%n_cnode + 1
           s%cnode(s%n_cnode)%i1=0
           s%cnode(s%n_cnode)%i2=isec%id
           s%cnode(s%n_cnode)%i3=m
           s%on_grid=.true.
           print *,'source lies on node ',m,' of intersection ',isec%id

           exit

        endif
     end do


     if (.not. s%on_grid) then   ! source lies on interface but not on an interface node

        print *,'source lies exactly on interface',iface,'but not on a node'

        !  find the nodes of intersection that are part of the cell containing the source

        m= grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p(iface)
        do i=1,isec%n_inodes(m)
              
           k=isec%inodes(i,m)
           s%n_cnode = s%n_cnode + 1
           s%cnode(s%n_cnode)%i1=0
           s%cnode(s%n_cnode)%i2=isec%id
           s%cnode(s%n_cnode)%i3=k

        end do

        ! the regular nodes of the cut cell containing the source


        do i=0,1
           is=s%ir+i
           do j=0,1
              js=s%ilat+j
              do k=0,1
                 ks=s%ilong+k
                 if (is > 0 .and. is <= grid%nr .and. &
                      js > 0 .and. js <= grid%nlat .and. &
                      ks > 0 .and. ks <= grid%nlong .and. &
                      grid%node_region(is,js,ks) == s%region_id) then
                    s%n_cnode = s%n_cnode + 1
                    s%cnode(s%n_cnode)%i1=is
                    s%cnode(s%n_cnode)%i2=js
                    s%cnode(s%n_cnode)%i3=ks
                 endif
              end do
           end do
        end do


     end if

  endif

! test if the source lies exactly on a regular grid node or not

  if (.not. s%on_interface) then

     do i=0,1
        do j=0,1
           do k=0,1

!              print *,s%r,s%lat,s%long
!              print *,s%ir+i,s%ilat+j,s%ilong+k
!              print *,grid%r(s%ir+i),grid%lat(s%ilat+j),grid%long(s%ilong+k)
!              print *

              dist=sqrt((s%r-grid%r(s%ir+i))**2 + (s%r*(s%lat-grid%lat(s%ilat+j)))**2 &
                   + (s%r*s%coslat*(s%long-grid%long(s%ilong+k)))**2 ) 
!              print *,i,j,k,dist
              if (dist < grid%tolerance) then



                 s%on_grid=.true.
                 s%n_cnode = s%n_cnode + 1
                 s%cnode(s%n_cnode)%i1=s%ir+i
                 s%cnode(s%n_cnode)%i2=s%ilat+j
                 s%cnode(s%n_cnode)%i3=s%ilong+k
                 print *,'source on grid but not interface'

              endif

           end do
        end do
     end do
 
             
     if (.not. s%on_grid) then

        print *,'source does not lie on grid'

! test if the cell in which the source resides is cut

        if (associated(grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p)) then

           print *,'source lies in a cut cell'

        ! if so, make a list of the nodes in the cut cell

        ! first the intersection nodes

           do n=1,2

              if (n == 1) isec => itop
              if (n == 2) isec => ibot

              if (grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p(isec%iface_id) /= 0) then

                 m= grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p(isec%iface_id)
                 do i=1,isec%n_inodes(m)
              
                    k=isec%inodes(i,m)
                    s%n_cnode = s%n_cnode + 1
                    s%cnode(s%n_cnode)%i1=0
                    s%cnode(s%n_cnode)%i2=isec%id
                    s%cnode(s%n_cnode)%i3=k
                 
                 end do
              endif
           end do

        !then the regular nodes

           do i=-1,2
              is=s%ir+i
              do j=-1,2
                 js=s%ilat+j
                 do k=-1,2
                    ks=s%ilong+k
                    if (is > 0 .and. is <= grid%nr .and. &
                         js > 0 .and. js <= grid%nlat .and. &
                         ks > 0 .and. ks <= grid%nlong .and. &
                         grid%node_region(is,js,ks) == s%region_id) then
                       s%n_cnode = s%n_cnode + 1
                       s%cnode(s%n_cnode)%i1=is
                       s%cnode(s%n_cnode)%i2=js
                       s%cnode(s%n_cnode)%i3=ks
!                       print '(6i5,i8)',i,j,k,s%cnode(s%n_cnode)%i1,is,js,ks,grid%rnode_id(is,js,ks)
                    endif
                 end do
              end do
           end do

!           do i=0,1
!              do j=0,1
!                 do k=0,1
!                    if (grid%node_region(s%ir+i,s%ilat+j,s%ilong+k) == s%region_id) then
!                       s%n_cnode = s%n_cnode + 1
!                       s%cnode%i1=s%ir+i
!                       s%cnode%i2=s%ilat+j
!                       s%cnode%i3=s%ilong+k
!                    endif
!                 end do
!              end do
!           end do

        else   ! the cell in which the source lies is not cut

        ! make a list of the nodes in the regular cell

           print *,'source lies in an uncut cell'

           do i=-1,2
              is=s%ir+i
              do j=-1,2
                 js=s%ilat+j
                 do k=-1,2
                    ks=s%ilong+k
                    if (is > 0 .and. is <= grid%nr .and. &
                         js > 0 .and. js <= grid%nlat .and. &
                         ks > 0 .and. ks <= grid%nlong .and. &
                         grid%node_region(is,js,ks) == s%region_id) then
                       s%n_cnode = s%n_cnode + 1
                       s%cnode(s%n_cnode)%i1=is
                       s%cnode(s%n_cnode)%i2=js
                       s%cnode(s%n_cnode)%i3=ks
!                       print '(6i5,i8)',i,j,k,s%cnode(s%n_cnode)%i1,is,js,ks,grid%rnode_id(is,js,ks)
                    endif
                 end do
              end do
           end do
        endif    ! cut cell or else

     endif  ! if not on grid

  endif  ! if not on interface

!  print *,'source nodes'
!  do n=1,s%n_cnode
!     print '(4i5)',n,s%cnode(n)%i1,s%cnode(n)%i2,s%cnode(n)%i3
!  end do

  return

end subroutine initialize_refined_source



!-----------------------------------------------------------------------------------------------------

! This is an alternative refined source initialization (currently not used)
! Initializes a larger region region around the source with simple analytical estimates

subroutine initialize_refined_source2(s,sc,grid,reg,itop,ibot)

  use mod_3dfm_nointerfaces
  implicit none

  integer       :: n,m,i,j,k,is,js,ks,source_node
  type(Tsource) :: s,sc
  type(Tpropagation_grid) :: grid
  type(Tintersection),target   :: itop,ibot
  type(Tregion)         :: reg
  type(Tintersection),pointer   :: isec
  real(kind=dp) :: dist,interpolate_interface



! copy position from coarse source
  s%r=sc%r
  s%lat=sc%lat
  s%long=sc%long
  s%region_id=reg%id

! determine grid cell in which the source is located

  s%ir    =  floor((s%r - grid%r0)/grid%dr0 + 1)
  s%ilat  =  floor((s%lat - grid%lat0)/grid%dlat0 + 1)
  s%ilong =  floor((s%long - grid%long0)/grid%dlong0 + 1)

!  if (s%ir < 1 .or. s%ir>grid%nr-1) stop 'source outside propagation grid'
!  if (s%ilat < 1 .or. s%ilat>grid%nlat-1) stop 'source outside propagation grid'
!  if (s%ilong < 1 .or. s%ilong>grid%nlong-1) stop 'source outside propagation grid'


  s%on_grid=.false.
  s%on_interface=.false.
  s%n_cnode = 0


!-----------------------------------------


! test where the source lies

  do n=1,2

     if (n == 1) isec => itop
     if (n == 2) isec => ibot

!     print *,'isec%id',n,isec%id,isec%iface_id


     if (isec%nnode /= 0) then  ! intersection corresponds to a real interface 

!     print *,n,s%r,interpolate_interface(s%lat,s%long,intrface(iface))

        dist = s%r-interpolate_interface(s%lat,s%long,intrface(isec%iface_id))


! first test if the source lies exactly on the interface

        if (abs(dist) < grid%tolerance) then  ! source lies on interface

           s%interface_id = isec%iface_id
           s%on_interface=.true.
           print *,'source lies exactly on interface',isec%id

        ! test if the source lies exactly on an intersection node

           do m=1,isec%nnode

              dist =sqrt((s%r-isec%r(m))**2 +(s%r*(s%lat-isec%lat(m)))**2 + &
                   (s%r*s%coslat*(s%long-isec%long(m)))**2 )

              if (dist < grid%tolerance) then  ! source lies exactly on an interface node

              ! assign time value 0 to the source node

                 s%n_cnode = s%n_cnode + 1
                 s%cnode(s%n_cnode)%i1=0
                 s%cnode(s%n_cnode)%i2=isec%id
                 s%cnode(s%n_cnode)%i3=m
                 s%on_grid=.true.

                 if (n == 1) source_node=isec%rbel_node_id(m)
                 if (n == 2) source_node=isec%rabo_node_id(m)

                 call get_source_neighbours(source_node,s,reg,grid)                 

                 print *,'source lies on node ',m,' of intersection ',isec%id
                 exit

              endif
           end do

           if (.not. s%on_grid) then   ! source lies on interface but not on an interface node

           !  find the nodes of intersection that are part of the cell containing the source

              m= grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p(n)
              do i=1,isec%n_inodes(m)
              
                 k=isec%inodes(i,m)
                 s%n_cnode = s%n_cnode + 1
                 s%cnode(s%n_cnode)%i1=0
                 s%cnode(s%n_cnode)%i2=isec%id
                 s%cnode(s%n_cnode)%i3=k

              end do
           end if

           exit ! we don't need to test other interfaces if the source lies exactly on the present one

        endif    ! test for source on interface

     endif  ! test for real intersection

  end do



! test if the source lies exactly on a regular grid node or not

  if (.not. s%on_interface) then

     do i=0,1
        do j=0,1
           do k=0,1

!              print *,s%r,s%lat,s%long
!              print *,grid%r(s%ir+i),grid%r(s%ilat++j),grid%long(s%ilong+k)
!              print *

              dist=sqrt((s%r-grid%r(s%ir+i))**2 + (s%r*(s%lat-grid%lat(s%ilat+j)))**2 &
                   + (s%r*s%coslat*(s%long-grid%long(s%ilong+k)))**2 ) 
!              print *,i,j,k,dist
              if (dist < grid%tolerance) then

                 s%on_grid=.true.
                 is=s%ir+i ; js=s%ilat+j ; ks=s%ilong+k

                 source_node= grid%rnode_id(is,js,ks)
                 write(22,*) 'snode =',is,js,ks
                 call get_source_neighbours(source_node,s,reg,grid)
                 exit

              endif

           end do
        end do
     end do
 
             
     if (.not. s%on_grid) then

! test if the cell in which the source resides is cut

        if (associated(grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p)) then

        ! if so, make a list of the nodes in the cut cell

        ! first the intersection nodes

           do n=1,2

              if (n == 1) isec => itop
              if (n == 2) isec => ibot

              if (grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p(isec%iface_id) /= 0) then

                 m= grid%ccind_from_3dc(s%ir,s%ilat,s%ilong)%p(isec%iface_id)
                 do i=1,isec%n_inodes(m)
              
                    k=isec%inodes(i,m)
                    s%n_cnode = s%n_cnode + 1
                    s%cnode(s%n_cnode)%i1=0
                    s%cnode(s%n_cnode)%i2=isec%id
                    s%cnode(s%n_cnode)%i3=k
                 
                 end do
              endif
           end do

        !then the regular nodes

           do i=0,1
              do j=0,1
                 do k=0,1
                    if (grid%node_region(s%ir+i,s%ilat+j,s%ilong+k) == s%region_id) then
                       s%n_cnode = s%n_cnode + 1
                       s%cnode%i1=s%ir+i
                       s%cnode%i2=s%ilat+j
                       s%cnode%i3=s%ilong+k
                    endif
                 end do
              end do
           end do

        else   ! the cell in which the source lies is not cut

        ! make a list of the nodes in the regular cell

           do i=0,1
              do j=0,1
                 do k=0,1

                    s%n_cnode = s%n_cnode + 1
                    s%cnode(s%n_cnode)%i1=s%ir+i
                    s%cnode(s%n_cnode)%i2=s%ilat+j
                    s%cnode(s%n_cnode)%i3=s%ilong+k

                 end do
              end do
           end do
        endif    ! cut cell or else

     endif  ! if not on grid

  endif  ! if not on interface


  do n=1,s%n_cnode
     write(22,*) n,s%cnode(n)%i1,s%cnode(n)%i2,s%cnode(n)%i3
  end do

  return

end subroutine initialize_refined_source2

!***********************************************************************************************




!***********************************************************************************************

!***********************************************************************************************
! This subroutine takes a source as input, constructs a refined grid around the source 
! sweeps through the refined grid, transfers the
! arrivaltimes from the refined grid regions to the regions in the main grid and does a sweep 
! through these regions in the upward and downward directions.

subroutine initialize_source_regions(s)
use mod_3dfm
implicit none

type(Tsource)       :: s  ! the source with its properties on the main grid
type(Tsource)       :: ss ! the source with its properties on the refined grid

type(Tintersection),pointer :: itop,ibot    ! top and bottom intersections of refined grid
type(Tintersection),pointer :: itopc,ibotc  ! top and bottom intersections of main grid
type(Tregion),pointer       :: reg          ! the source region in the main grid
type(Tregion),pointer       :: sreg         ! the refined source region 


! local stuff
integer :: n,m,i,j,k,i1,i2,i3,prev_tf,nstart
integer :: nrmax,nrmin,nlatmax,nlatmin,nlongmax,nlongmin,vtype
real(kind=dp)  :: rmin,rmax,latmin,latmax,longmin,longmax,t_short
integer        ::t1,t2

call system_clock(t1)

! first construct the refined grid around the source

allocate(sgrid)
call pgrid_defaults(sgrid)

sgrid%is_source_grid = .true.

! limits of volume of main grid to be refined, taking main grid boundaries into account

i =  nint((s%r - pgrid%r0)/pgrid%dr0 + 1)
j =  nint((s%lat - pgrid%lat0)/pgrid%dlat0 + 1)
k =  nint((s%long - pgrid%long0)/pgrid%dlong0 + 1)

nrmax = min(pgrid%nr , i + ncell_to_be_refined)
nrmin = max( 1 , i - ncell_to_be_refined)
nlatmax = min(pgrid%nlat , j + ncell_to_be_refined)
nlatmin = max( 1 , j - ncell_to_be_refined)
nlongmax = min(pgrid%nlong , k + ncell_to_be_refined)
nlongmin = max( 1 , k - ncell_to_be_refined)

! origin of the refined source grid in propagation grid index coordinates

sgrid%index_r0 = nrmin
sgrid%index_lat0 = nlatmin
sgrid%index_long0 = nlongmin

! calculate the refined source grid parameters

sgrid%nr      = (nrmax - nrmin)*refinement_factor + 1
sgrid%nlat    = (nlatmax - nlatmin)*refinement_factor + 1
sgrid%nlong   = (nlongmax - nlongmin)*refinement_factor + 1
sgrid%dr0     = pgrid%dr0/refinement_factor
sgrid%dlat0   = pgrid%dlat0/refinement_factor
sgrid%dlong0  = pgrid%dlong0/refinement_factor
sgrid%r0      = pgrid%r(s%ir) - (s%ir - nrmin)*pgrid%dr0
sgrid%lat0    = pgrid%lat(s%ilat) - (s%ilat - nlatmin)*pgrid%dlat0
sgrid%long0   = pgrid%long(s%ilong) - (s%ilong - nlongmin)*pgrid%dlong0


! initialize the grid coordinates

allocate(sgrid%r(sgrid%nr),sgrid%lat(sgrid%nlat),sgrid%long(sgrid%nlong),sgrid%coslat(sgrid%nlat))

do i=1,sgrid%nr
   sgrid%r(i)=sgrid%r0 + (i-1)*sgrid%dr0
end do

do i=1,sgrid%nlat
   sgrid%lat(i)=sgrid%lat0 + (i-1)*sgrid%dlat0
end do

do i=1,sgrid%nlat
   sgrid%coslat(i)=cos(sgrid%lat(i))
end do


do i=1,sgrid%nlong
   sgrid%long(i)=sgrid%long0 + (i-1)*sgrid%dlong0
end do

sgrid%tolerance = refinement_factor*interface_tolerance*sgrid%dr0


! allocate storage for the refined source grid, intersections and regions

allocate(sgrid%rnode_id(sgrid%nr,sgrid%nlat,sgrid%nlong))
sgrid%rnode_id=0
allocate(sgrid%node_region(sgrid%nr,sgrid%nlat,sgrid%nlong))
sgrid%node_region=0
allocate(sgrid%ccind_from_3dc(sgrid%nr,sgrid%nlat,sgrid%nlong))
do k=1,sgrid%nlong
   do j=1,sgrid%nlat
      do i=1,sgrid%nr
         nullify(sgrid%ccind_from_3dc(i,j,k)%p)
      end do
   end do
end do

allocate(sgrid%arrivaltime(sgrid%nr,sgrid%nlat,sgrid%nlong))
sgrid%arrivaltime = huge_time

n_sintersections=n_intersections
allocate(sintersection(n_sintersections))
do i=1,n_sintersections ; call intersection_defaults(sintersection(i)) ; sintersection(i)%id=i ; end do
n_sregions=n_regions
allocate(sregion(n_sregions))
do i=1,n_sregions ; call region_defaults(sregion(i)) ; sregion(i)%id=i ; end do

do n=1,n_sintersections
   call find_intersection(sintersection(n),intrface(n),sgrid)
   if (sintersection(n)%nnode > 0) then
      allocate(sintersection(n)%arrivaltime(sintersection(n)%nnode))
      allocate(sintersection(n)%time_gradient(3,sintersection(n)%nnode))
   endif
end do

print *,'intersections found'


! set a flag at nodes on the grid that are completely regular, i.e. none of the connected cells
! has irregular nodes

call tag_regular_nodes(sgrid)
print *,'nodes of refined source grid tagged'



do n=1,n_sregions

     sregion(n)%id=n

! pointers to the intersections that are the boundaries of the region
     sregion(n)%itop => sintersection(n)
     sregion(n)%ibot => sintersection(n+1)
     sregion(n)%ivgrid=region(n)%ivgrid

     if (sregion(n)%itop%nnode == 0 .and.sregion(n)%ibot%nnode == 0 .and. &
          sregion(n)%id /= s%region_id ) then

        print *,'sregion',n,'does not exist in refined source grid'
        sregion(n)%nnode = 0
        cycle

     endif

   call define_region(sregion(n),sintersection(n),sintersection(n+1),sgrid)

   print *,'refined source region',n,' defined, nnode =',sregion(n)%nnode

end do


! transfer the velocity values to all refined grid nodes

call velocities_to_grid(sgrid)

do n=1,n_sintersections
   if (sintersection(n)%nnode /= 0) call velocities_to_intersection(sintersection(n))
end do
do n=1,n_sregions
  if (sregion(n)%nnode > 0)  call velocities_to_region(sregion(n),sgrid)
end do

print *,'refined velocities transferred'



!!------------------------------------------------------------------------------------------------
! determine the paths up and down from the source region covering the main regions overlapping 
! with the refined source region
!------------------------------------------------------------------------------------------------

! set default values for the intialization paths

do i=1,2
   do j=1,2
      call path_defaults(s%init_path(i,j))
      s%init_path(i,j)%used = .false.
   end do
end do

! first construct the path up if it exists

do vtype=1,n_vtypes

! ! first construct the path up if it exists ( init_path(1,vtype) )

if (.not.(s%on_interface .and. s%topint_id == 1)) then

!   print *,'constructing source sequence up'

   if (s%on_interface) then
      nstart = s%topreg_id - 1
   else
      nstart = s%region_id - 1
   endif

   s%init_path(1,vtype)%n_tf = 1
   do n=nstart,1,-1     ! count the number of timefields up
      if (sregion(n)%nnode > 0) then
         s%init_path(1,vtype)%n_tf = s%init_path(1,vtype)%n_tf + 1
      endif
   end do

   allocate(s%init_path(1,vtype)%sequence(2*s%init_path(1,vtype)%n_tf))
   allocate(s%init_path(1,vtype)%tf_sequence(s%init_path(1,vtype)%n_tf))

   s%init_path(1,vtype)%id  = 0
   s%init_path(1,vtype)%valid= .true.
   s%init_path(1,vtype)%used = .true.
   s%init_path(1,vtype)%refstep = 0
   s%init_path(1,vtype)%fitting_interface = 0

   s%init_path(1,vtype)%sequence(1) = 0
   if (s%on_interface) then
      s%init_path(1,vtype)%sequence(2) = region(s%topreg_id)%itop%iface_id
   else
      s%init_path(1,vtype)%sequence(2) = region(s%region_id)%itop%iface_id
   endif

! construct the sequence of the first path
   if (s%init_path(1,vtype)%n_tf > 1) then   
      do i=2,s%init_path(1,vtype)%n_tf
         s%init_path(1,vtype)%sequence(2*i-1) = s%init_path(1,vtype)%sequence(2*i-2)
         s%init_path(1,vtype)%sequence(2*i) = s%init_path(1,vtype)%sequence(2*i-1) - 1
      end do
   endif

endif


! ! then construct the path down if it exists ( init_path(2,vtype) )

if (.not.(s%on_interface .and. s%botint_id == n_interfaces)) then

!   print *,'constructing source sequence down'

   if (s%on_interface) then
      nstart = s%botreg_id + 1
   else
      nstart = s%region_id + 1
   endif

   s%init_path(2,vtype)%n_tf = 1
   do n=nstart,n_sregions     ! count the number of timefields down
      if (sregion(n)%nnode > 0) then
         s%init_path(2,vtype)%n_tf = s%init_path(2,vtype)%n_tf + 1
      endif
   end do

   allocate(s%init_path(2,vtype)%sequence(2*s%init_path(2,vtype)%n_tf))
   allocate(s%init_path(2,vtype)%tf_sequence(s%init_path(2,vtype)%n_tf))

   s%init_path(2,vtype)%id  = 0
   s%init_path(2,vtype)%valid= .true.
   s%init_path(2,vtype)%used = .true.
   s%init_path(2,vtype)%refstep = 0
   s%init_path(2,vtype)%fitting_interface = 0

   s%init_path(2,vtype)%sequence(1) = 0
   if (s%on_interface) then
      s%init_path(2,vtype)%sequence(2) = region(s%botreg_id)%ibot%iface_id
   else
      s%init_path(2,vtype)%sequence(2) = region(s%region_id)%ibot%iface_id
   endif

! construct the sequence of the second path
   if (s%init_path(2,vtype)%n_tf > 1) then   
      do i=2,s%init_path(2,vtype)%n_tf
         s%init_path(2,vtype)%sequence(2*i-1) = s%init_path(2,vtype)%sequence(2*i-2)
         s%init_path(2,vtype)%sequence(2*i) = s%init_path(2,vtype)%sequence(2*i-1) + 1
      end do
   endif

endif

end do   ! vtypes

! paths through the regions overlapping with the source region have been defined


! do the initialization for the number of velocity types in the problem

do vtype=1,n_vtypes

print *,'--------------------------------------------------'
print *,'starting source initialization for vtype =',vtype

!--------------------------------------------------------
! Do the propagation through the source regions first
   
! if the source lies on an interface initialize both regions above and below, else only source region

do n=1,s%n_tf_init

   if (s%on_interface) then   ! initialize the regions above and below the interface

      if (s%n_tf_init == 2) then
         if (n==1) sreg => sregion(s%topreg_id)
         if (n==2) sreg => sregion(s%botreg_id)
      else
         if (s%topreg_id > 0) sreg => sregion(s%topreg_id)
         if (s%botreg_id < n_regions) sreg => sregion(s%botreg_id)        
      endif

   else                       ! only he one region in which the source lies
      sreg => sregion(s%region_id)
   endif

   itop => sintersection(sreg%id) ! top intersection of the refined source region
   ibot => sintersection(sreg%id+1)   ! bottom intersection of the refined source region


   ! initialize the source in the refined source grid. 

   call initialize_refined_source(ss,s,sgrid,sreg,itop,ibot)

   print *,'refined source initialized in sregion',sreg%id


   ! sweep the refined source grid

   call sweep_region_from_source(sreg,ss,vtype)

   print *,'refined source region',sreg%id,' swept'


   ! transfer regional travel times to interfaces and regular grid

!  if (itop%nnode > 0 .and. .not. associated(itop%arrivaltime)) allocate(itop%arrivaltime(itop%nnode))
!  if (ibot%nnode > 0 .and. .not. associated(ibot%arrivaltime)) allocate(ibot%arrivaltime(ibot%nnode))


   do i=1,sreg%nnode

      i1 = sreg%node(i)%i1 ; i2 = sreg%node(i)%i2 ; i3 = sreg%node(i)%i3

      if (i1 == 0) then

         sintersection(i2)%arrivaltime(i3) = sreg%arrivaltime(i)
         sintersection(i2)%time_gradient(1:3,i3) = sreg%time_gradient(1:3,i)

      else

         sgrid%arrivaltime(i1,i2,i3)= sreg%arrivaltime(i)

      endif

   end do

enddo



!--------------------------------------------------------------------------------
! now sweep up and down through the remaining regions of the refined source grid

! sweep up

if (.not.s%on_interface .or. (s%on_interface .and. s%topreg_id > 1)) then

   if (s%on_interface) then
      nstart = s%topreg_id - 1
   else
      nstart = s%region_id - 1
   endif

   do n=nstart,1,-1

      sreg => sregion(n)

      if (sreg%nnode > 0) then

         call refract_gradient(sintersection(n+1),sregion(sreg%id+1),vtype,1)

         call sweep_sregion_from_interface(sreg,sintersection(n+1),vtype)

         print *,'refined region',n,' swept'

! transfer regional travel times to interfaces 

      
         do i=1,sreg%nnode

            i1 = sreg%node(i)%i1 ; i2 = sreg%node(i)%i2 ; i3 = sreg%node(i)%i3

            if (i1 == 0) then

               sintersection(i2)%arrivaltime(i3) = sreg%arrivaltime(i)
               sintersection(i2)%time_gradient(1:3,i3) = sreg%time_gradient(1:3,i)
            else

               sgrid%arrivaltime(i1,i2,i3)= sreg%arrivaltime(i)

            endif

         end do

      endif

   end do

endif


! sweep down

if (.not.s%on_interface .or. (s%on_interface .and. s%botreg_id < n_regions)) then

   if (s%on_interface) then
      nstart = s%botreg_id + 1
   else
      nstart = s%region_id + 1
   endif

   do n=nstart,n_sregions

      sreg => sregion(n)

      if (sreg%nnode > 0) then

         print *,'refracting at interface',n
         call refract_gradient(sintersection(n),sregion(sreg%id-1),vtype,-1)

         print *,'starting sweep from interface',n,'into region',sreg%id
         call sweep_sregion_from_interface(sreg,sintersection(n),vtype)

         print *,'refined region',n,' swept'

         do i=1,sreg%nnode

            i1 = sreg%node(i)%i1 ; i2 = sreg%node(i)%i2 ; i3 = sreg%node(i)%i3

            if (i1 == 0) then

!               write(1,*) i,i1,i2,i3 

               sintersection(i2)%arrivaltime(i3) = sreg%arrivaltime(i)
               sintersection(i2)%time_gradient(1:3,i3) = sreg%time_gradient(1:3,i)


            else

!               write(1,*) i,i1,i2,i3

               sgrid%arrivaltime(i1,i2,i3)= sreg%arrivaltime(i)

            endif

         end do
         
      endif

   end do

endif

call system_clock(t2)
print *,'refined grid for source',s%id,'took', real(t2-t1)/10000.,' sec'
print *

!-------------------------------------------------------------------------------------------------
! transfer the results to the region(s) in which the source resides on the main propagation grid

! first find shortest traveltime on refined region boundary

t_short=huge_time

rmin=pgrid%r(1)+pgrid%tolerance
rmax=pgrid%r(pgrid%nr)-pgrid%tolerance
latmin=pgrid%lat(1)+pgrid%tolerance/pgrid%dr0
latmax=pgrid%lat(pgrid%nlat)-pgrid%tolerance/pgrid%dr0
longmin=pgrid%long(1)+pgrid%tolerance/pgrid%dr0
longmax=pgrid%long(pgrid%nlong)-pgrid%tolerance/pgrid%dr0

if (sgrid%r(1) > rmin)  &
     t_short = min(t_short,minval(sgrid%arrivaltime(1,1:sgrid%nlat,1:sgrid%nlong)))
if (sgrid%r(sgrid%nr) < rmax) &
     t_short = min(t_short,minval(sgrid%arrivaltime(sgrid%nr,1:sgrid%nlat,1:sgrid%nlong)))
if (sgrid%lat(1) > latmin) &
     t_short = min(t_short,minval(sgrid%arrivaltime(1:sgrid%nr,1,1:sgrid%nlong)))
if (sgrid%lat(sgrid%nlat) < latmax) &
     t_short = min(t_short,minval(sgrid%arrivaltime(1:sgrid%nr,sgrid%nlat,1:sgrid%nlong)))
if (sgrid%long(1) > longmin) &
     t_short = min(t_short,minval(sgrid%arrivaltime(1:sgrid%nr,1:sgrid%nlat,1)))
if (sgrid%long(sgrid%nlong) < longmax) &
     t_short = min(t_short,minval(sgrid%arrivaltime(1:sgrid%nr,1:sgrid%nlat,sgrid%nlong)))

print *,'shortest time on refined grid boundary is',t_short



! transfer the refined source region(s) to the main propagation grid source region(s)
! if the source lies on an interface initialize both regions above and below, else only source region

do n=1,2

   if (s%init_path(n,vtype)%used) then

      if (s%on_interface) then   ! initialize the regions above and below the interface

         if (n==1) then ; sreg => sregion(s%topreg_id) ; reg => region(s%topreg_id) ; endif
         if (n==2) then ; sreg => sregion(s%botreg_id) ; reg => region(s%botreg_id) ; endif

      else                       ! only the one region in which the source lies

         sreg => sregion(s%region_id) ; reg => region(s%region_id)

      endif

      itopc => intersection(reg%id)      ! the top intersection of the normal source region
      ibotc => intersection(reg%id+1)    ! the bottom intersection of the normal source region
      allocate(reg%arrivaltime(reg%nnode),reg%time_gradient(3,reg%nnode),reg%node_status(reg%nnode))
      reg%node_status=-1
      reg%arrivaltime = huge_time


      print *,'before transefr refined region'

      call transfer_refined_region(sreg,reg,t_short)

      print *,'transferred times from fine grid to coarse for region',reg%id


   ! create a narrow band around the nodes transferred from the fine grid

      call create_narrow_band(reg,vtype)

      print *,'narrow band created in region',reg%id,'alive/nb',count(reg%node_status == 0),&
           count(reg%node_status == 1)

! do the fast marching sweep across the main grid region containing the source

      call propagate(reg,vtype)

      print *,'propagation through region',reg%id,' finished'

! transfer regional travel times to interfaces 

      if (itopc%nnode > 0 .and. .not. associated(itopc%arrivaltime)) then
         allocate(itopc%arrivaltime(itopc%nnode))
         allocate(itopc%time_gradient(3,itopc%nnode))
      endif
      if (ibotc%nnode > 0 .and. .not. associated(ibotc%arrivaltime)) then
         allocate(ibotc%arrivaltime(ibotc%nnode))
         allocate(ibotc%time_gradient(3,ibotc%nnode))
      endif

      do i=1,reg%nnode

         i1 = reg%node(i)%i1 ; i2 = reg%node(i)%i2 ; i3 = reg%node(i)%i3

         if (i1 == 0) then

            intersection(i2)%arrivaltime(i3) = reg%arrivaltime(i)
            intersection(i2)%time_gradient(1:3,i3) = reg%time_gradient(1:3,i)

         endif

      end do
      
  ! transfer time field to the array of saved time fields

      s%n_time_fields=s%n_time_fields+1
      allocate(s%time_field(s%n_time_fields)%arrivaltime(reg%nnode))
      s%time_field(s%n_time_fields)%arrivaltime=reg%arrivaltime
      allocate(s%time_field(s%n_time_fields)%time_gradient(3,reg%nnode))
      s%time_field(s%n_time_fields)%time_gradient=reg%time_gradient

      print *,'results written to timefield',s%n_time_fields

  ! pointer to source region

      s%time_field(s%n_time_fields)%reg =>  reg 

      s%time_field(s%n_time_fields)%vtype =  vtype

      s%init_path(n,vtype)%tf_sequence(1) = s%n_time_fields

      if (s%on_interface) then
         if (reg%id == s%topreg_id) then
            s%time_field(s%n_time_fields)%istart => intersection(s%topint_id)
            s%time_field(s%n_time_fields)%inonstart =>  intersection(s%topint_id - 1)
         else
            s%time_field(s%n_time_fields)%istart => intersection(s%botint_id)
            s%time_field(s%n_time_fields)%inonstart =>  intersection(s%botint_id + 1)
         endif
      else
         s%time_field(s%n_time_fields)%istart => ibotc
         s%time_field(s%n_time_fields)%inonstart =>  itopc
      endif

  ! register the first (source) time fields with the associated source

      if (n == 1) s%first_tf_up(vtype) = s%n_time_fields
      if (n == 2) s%first_tf_down(vtype) = s%n_time_fields


! deallocate  everything that is no longer required

      deallocate(reg%arrivaltime,reg%time_gradient,reg%node_status)


  ! if the source does not lie on an interface, the first timefields up and down are the same
  ! we do not need to do the direction down (loop index n = 2)
      if (.not.s%on_interface) then
         if (n==2) stop 'error n==2 in source region init'
         s%init_path(2,vtype)%tf_sequence(1) = s%n_time_fields         
         s%first_tf_down(vtype) = s%n_time_fields
         exit
      endif


   endif

end do

! timefields 1(/2) describing the source region(s) on the main grid ha(s)(ve) been established
print *,'timefields in main grid regions connected to the source established, vtype=',vtype



!--------------------------------------------------------------------------------------------
! now construct the actual time field sequences on the main grid up (path 1) or down (path 2)
! through the regions overlapping with the refined source grid


! the upward path 1

! only if there are regions above the source region

if (s%init_path(1,vtype)%used .and. s%init_path(1,vtype)%n_tf > 1) then 

   prev_tf = s%first_tf_up(vtype)

   print *,'starting upward sweep in main grid regions overlapping source region'

   do n=2,s%init_path(1,vtype)%n_tf


! transfer the refined source region to the main propagation grid region

      if (s%on_interface) then
         reg => region(s%topreg_id-n+1)
      else
         reg =>  region(s%region_id-n+1)         ! the region in the main propagation grid
      endif

      itopc => reg%itop                       ! the top intersection of the region
      ibotc => reg%ibot                       ! the bottom intersection of the region
      allocate(reg%arrivaltime(reg%nnode),reg%time_gradient(3,reg%nnode),reg%node_status(reg%nnode))
      reg%node_status=-1
      reg%arrivaltime = huge_time

      call transfer_refined_region(sregion(reg%id),reg,t_short)

      print *,'transferred times from fine grid to coarse for region',reg%id


 ! create narrow band around transferred nodes

      call create_narrow_band(reg,vtype)
      print *,'narrow band created in region',reg%id,'alive/nb',count(reg%node_status == 0),&
           count(reg%node_status == 1)


 ! add the intersection nodes to the narrow band

      if (.not. associated(ibotc%arrivaltime)) &
           stop 'requested starting interface does not have arrival times'

      allocate(ibotc%starttime(ibotc%nnode)) 
! start time will be compared to arrival time after propagation to find turning rays


! set start time, time gradient and node status narrow band in regional 
! nodes that belong to the starting interface

      call refract_gradient(ibotc,region(ibotc%iface_id),vtype,1)

      do i=1,ibotc%nnode

         j=ibotc%rabo_node_id(i)
         if (reg%node_status(j) /= 0) then     

    ! only for nodes that did not receive a value from the refined grid

            reg%arrivaltime(j)= ibotc%arrivaltime(i)
            reg%time_gradient(1:3,j)=ibotc%time_gradient(1:3,i)
            reg%node_status(j)=1

         endif
         ibotc%starttime(i)=reg%arrivaltime(j)
      end do

!      print *,'region',reg%id,' starting times set'


! do the fast marching sweep across the main grid region 

      call propagate(reg,vtype)

      print *,'propagation through region',reg%id,' finished'


! transfer regional travel times to interfaces 

      if (itopc%nnode > 0 .and. .not. associated(itopc%arrivaltime)) then
         allocate(itopc%arrivaltime(itopc%nnode))
         allocate(itopc%time_gradient(3,itopc%nnode))
      endif
      if (ibotc%nnode > 0 .and. .not. associated(ibotc%arrivaltime)) then
         allocate(ibotc%arrivaltime(ibotc%nnode))
         allocate(ibotc%time_gradient(3,ibotc%nnode))
      endif

      do i=1,reg%nnode

         i1 = reg%node(i)%i1 ; i2 = reg%node(i)%i2 ; i3 = reg%node(i)%i3

         if (i1 == 0) then

            intersection(i2)%arrivaltime(i3) = reg%arrivaltime(i)
            intersection(i2)%time_gradient(1:3,i3) = reg%time_gradient(1:3,i)

         endif

      end do
      
  ! transfer time field to the array of saved time fields

      s%n_time_fields=s%n_time_fields+1
      allocate(s%time_field(s%n_time_fields)%arrivaltime(reg%nnode))
      s%time_field(s%n_time_fields)%arrivaltime=reg%arrivaltime
      allocate(s%time_field(s%n_time_fields)%time_gradient(3,reg%nnode))
      s%time_field(s%n_time_fields)%time_gradient=reg%time_gradient

      print *,'results written to timefield',s%n_time_fields

     ! attach info to the generated time field

      ! time field preceding the current one
      s%time_field(s%n_time_fields)%prev_tf = prev_tf

      ! identify the current time field as the child of the previous time field
      if (prev_tf <= s%n_tf_init) then
         s%time_field(prev_tf)%next_tf(1+(vtype-1)*4) = s%n_time_fields
      else
         s%time_field(prev_tf)%next_tf(1+(vtype-1)*4) = s%n_time_fields
      endif

      ! store the index of the time field in the sequence of timefields for this path
      s%init_path(1,vtype)%tf_sequence(n) = s%n_time_fields

      ! pointers to start and non-start interfaces
      s%time_field(s%n_time_fields)%istart => ibotc
      s%time_field(s%n_time_fields)%inonstart =>  itopc  
      s%time_field(s%n_time_fields)%reg =>  reg
  
      s%time_field(s%n_time_fields)%vtype =  vtype

      prev_tf = s%n_time_fields

     ! check for turning rays
      if (count(ibotc%arrivaltime /= ibotc%starttime) > 0) then
         print *,'turning rays were present on interface',ibotc%iface_id
         s%time_field(s%n_time_fields)%turning_rays_present=.true.
         allocate(s%time_field(s%n_time_fields)%received_turning_ray(ibotc%nnode)) 
         s%time_field(s%n_time_fields)%received_turning_ray = ibotc%arrivaltime /= ibotc%starttime
      endif

! deallocate  everything that is no longer required

      deallocate(reg%arrivaltime,reg%time_gradient,reg%node_status,ibotc%starttime)

   end do  ! end of loop over steps of upward path


endif   ! upward path exist



! the downward path 2

! only if the region below the source region exists

if (s%init_path(2,vtype)%used .and. s%init_path(2,vtype)%n_tf > 1) then

   prev_tf =s%first_tf_down(vtype)

   print *,'starting downward sweep in main grid regions overlapping source region'

   do n=2,s%init_path(2,vtype)%n_tf

      print *,'timefield',n,' of the downward sequence started'

! transfer the refined source region to the main propagation grid region

      if (s%on_interface) then
         reg => region(s%botreg_id+n-1)
      else
         reg =>  region(s%region_id+n-1)         ! the region in the main propagation grid
      endif

      itopc => reg%itop                       ! the top intersection of the region
      ibotc => reg%ibot                       ! the bottom intersection of the region
      allocate(reg%arrivaltime(reg%nnode),reg%time_gradient(3,reg%nnode),reg%node_status(reg%nnode))
      reg%node_status=-1
      reg%arrivaltime = huge_time

      call transfer_refined_region(sregion(reg%id),reg,t_short)

      print *,'transferred times from fine grid to coarse for region',reg%id


 ! create narrow band around transferred nodes

      call create_narrow_band(reg,vtype)
      print *,'narrow band created in region',reg%id,'alive/nb',count(reg%node_status == 0),&
           count(reg%node_status == 1)


 ! add the intersection nodes to the narrow band

      if (.not. associated(itopc%arrivaltime)) &
           stop 'requested starting interface does not have arrival times'

      allocate(itopc%starttime(itopc%nnode))


! set start time, time gradient and node status narrow band in regional 
! nodes that belong to the starting interface

      call refract_gradient(itopc,region(itopc%iface_id-1),vtype,-1)

      do i=1,itopc%nnode

         j=itopc%rbel_node_id(i)
         if (reg%node_status(j) /= 0) then     

       ! only for nodes that did not receive a value from the refined grid

            reg%arrivaltime(j)= itopc%arrivaltime(i)
            reg%time_gradient(1:3,j)=itopc%time_gradient(1:3,i)
            reg%node_status(j)=1

         endif
         itopc%starttime(i)=reg%arrivaltime(j)
      end do


! do the fast marching sweep across the main grid region 

      call propagate(reg,vtype)

      print *,'propagation through region',reg%id,' finished'

! transfer regional travel times to interfaces 

      if (itopc%nnode > 0 .and. .not. associated(itopc%arrivaltime)) then
         allocate(itopc%arrivaltime(itopc%nnode))
         allocate(itopc%time_gradient(3,itopc%nnode))
      endif
      if (ibotc%nnode > 0 .and. .not. associated(ibotc%arrivaltime)) then
         allocate(ibotc%arrivaltime(ibotc%nnode))
         allocate(ibotc%time_gradient(3,ibotc%nnode))
      endif

      do i=1,reg%nnode

         i1 = reg%node(i)%i1 ; i2 = reg%node(i)%i2 ; i3 = reg%node(i)%i3

         if (i1 == 0) then

            intersection(i2)%arrivaltime(i3) = reg%arrivaltime(i)
            intersection(i2)%time_gradient(1:3,i3) = reg%time_gradient(1:3,i)

         endif

      end do
      
  ! transfer time field to the array of saved time fields

      s%n_time_fields=s%n_time_fields+1
      allocate(s%time_field(s%n_time_fields)%arrivaltime(reg%nnode))
      s%time_field(s%n_time_fields)%arrivaltime=reg%arrivaltime
      allocate(s%time_field(s%n_time_fields)%time_gradient(3,reg%nnode))
      s%time_field(s%n_time_fields)%time_gradient=reg%time_gradient

      print *,'results written to timefield',s%n_time_fields

     ! attach info to the generated time field

      ! time field preceding the current one
      s%time_field(s%n_time_fields)%prev_tf = prev_tf

      ! identify the current time field as the child of the previous time field
      if (prev_tf <= s%n_tf_init) then
         s%time_field(prev_tf)%next_tf(4+(vtype-1)*4) = s%n_time_fields
!         print *,'timefield',s%n_time_fields,'is ctype 4 of time field',prev_tf
      else
         s%time_field(prev_tf)%next_tf(2+(vtype-1)*4) = s%n_time_fields
!         print *,'timefield',s%n_time_fields,'is ctype 2 of time field',prev_tf
      endif

      ! store the index of the time field in the sequence of timefields for this path
      s%init_path(2,vtype)%tf_sequence(n) = s%n_time_fields

      ! pointers to start and non-start interfaces
      s%time_field(s%n_time_fields)%istart => itopc
      s%time_field(s%n_time_fields)%inonstart =>  ibotc  
      s%time_field(s%n_time_fields)%reg =>  reg  

      s%time_field(s%n_time_fields)%vtype =  vtype

      prev_tf = s%n_time_fields

     ! check for turning rays
      if (count(itopc%arrivaltime /= itopc%starttime) > 0) then
         print *,'turning rays were present on interface',itopc%iface_id
         s%time_field(s%n_time_fields)%turning_rays_present=.true.
         allocate(s%time_field(s%n_time_fields)%received_turning_ray(itopc%nnode)) 
         s%time_field(s%n_time_fields)%received_turning_ray = itopc%arrivaltime /= itopc%starttime
      endif


! deallocate  everything that is no longer required

      deallocate(reg%arrivaltime,reg%time_gradient,reg%node_status,itopc%starttime)

!      print *,'timefield',n,' of the downward sequence finished'

   end do  ! end of loop over steps of downward path

endif   ! downward path exist

end do  ! vtypes

! print *,'finished inside source_reg_init'

do m=1,n_vtypes
do n=1,2
   if (s%init_path(n,m)%used) then
      if (n==1) print '(a45,i3,a1,11i5)','timefields on upward init path for vtype',m,':', &
           s%init_path(n,m)%tf_sequence(1:s%init_path(n,m)%n_tf)
      if (n==2) print '(a45,i3,a1,11i5)','timefields on downward init path for vtype',m,':', &
           s%init_path(n,m)%tf_sequence(1:s%init_path(n,m)%n_tf)
   else
      if (n==1) print *,'there was no upward initialization path for vtype',m 
      if (n==2) print *,'there was no downward initialization path for vtype',m
   endif
end do
end do


! deallocate all storage related to the refined source grid

call clean_grid(sgrid)
deallocate(sgrid)
do n=1,n_sintersections
   call clean_intersection(sintersection(n))
end do
deallocate(sintersection)
do n=1,n_sregions
   call clean_region(sregion(n))
end do
deallocate(sregion)

! stop ' temp stop in initsource'

return

end subroutine initialize_source_regions

!***************************************************************

subroutine clean_grid(grid)

use mod_3dfm
implicit none

integer :: i,j,k

type(Tpropagation_grid) :: grid

if (associated(grid%r)) deallocate(sgrid%r)
if (associated(grid%lat)) deallocate(grid%lat)
if (associated(grid%long)) deallocate(grid%long)
if (associated(grid%arrivaltime)) deallocate(grid%arrivaltime)
if (associated(grid%time_gradient)) deallocate(grid%time_gradient)
if (associated(grid%rnode_id)) deallocate(grid%rnode_id)
if (associated(grid%node_region)) deallocate(grid%node_region)
if (associated(grid%velocity)) deallocate(grid%velocity)
if (associated(grid%coslat)) deallocate(grid%coslat)
if (associated(grid%fully_regular)) deallocate(grid%fully_regular)

if (associated(grid%ccind_from_3dc)) then
   do i=1,grid%nr
      do j=1,grid%nlat
         do k=1,grid%nlong
            if (associated(grid%ccind_from_3dc(i,j,k)%p)) deallocate(grid%ccind_from_3dc(i,j,k)%p)
         end do
      end do
   end do
   deallocate(grid%ccind_from_3dc)
endif

return

end subroutine clean_grid

!*******************************************************************

subroutine clean_intersection(isec)

use mod_3dfm
implicit none

type(Tintersection) :: isec

if (associated(isec%r)) deallocate(isec%r) 
if (associated(isec%lat)) deallocate(isec%lat)
if (associated(isec%coslat)) deallocate(isec%coslat)
if (associated(isec%long)) deallocate(isec%long)
if (associated(isec%normal)) deallocate(isec%normal)
if (associated(isec%vel_top)) deallocate(isec%vel_top)
if (associated(isec%vel_bot)) deallocate(isec%vel_bot)
if (associated(isec%arrivaltime)) deallocate(isec%arrivaltime)
if (associated(isec%starttime)) deallocate(isec%starttime)
if (associated(isec%time_gradient)) deallocate(isec%time_gradient)
if (associated(isec%intype)) deallocate(isec%intype)
if (associated(isec%ccells)) deallocate(isec%ccells)
if (associated(isec%n_inodes)) deallocate(isec%n_inodes)
if (associated(isec%inodes)) deallocate(isec%inodes)
if (associated(isec%ccell_from_inode)) deallocate(isec%ccell_from_inode)
if (associated(isec%rabo_node_id)) deallocate(isec%rabo_node_id)
if (associated(isec%rbel_node_id)) deallocate(isec%rbel_node_id)
if (associated(isec%irg_abo)) deallocate(isec%irg_abo)
if (associated(isec%irg_bel)) deallocate(isec%irg_bel)


end subroutine clean_intersection

!*******************************************************************

subroutine clean_region(reg)

use mod_3dfm
implicit none
type(Tregion)  :: reg


if (associated(reg%node)) deallocate(reg%node) 
if (associated(reg%node_status)) deallocate(reg%node_status)
if (associated(reg%arrivaltime)) deallocate(reg%arrivaltime)
if (associated(reg%time_gradient)) deallocate(reg%time_gradient)
if (associated(reg%velocity)) deallocate(reg%velocity)
if (associated(reg%r)) deallocate(reg%r)
if (associated(reg%lat)) deallocate(reg%lat)
if (associated(reg%coslat)) deallocate(reg%coslat)
if (associated(reg%long)) deallocate(reg%long)
if (associated(reg%init_id)) deallocate(reg%init_id)
if (associated(reg%init_arrivaltime)) deallocate(reg%init_arrivaltime)
if (associated(reg%init_time_gradient)) deallocate(reg%init_time_gradient)

end subroutine clean_region


!-----------------------------------------------------------------------------------------
! this logical function returns true if centernode has neighbours that are not alive
! used to evaluate status of nodes that received an arrival time from the refined source grid

  function non_alive_neighbours(centernode,reg,grid)
    use mod_3dfm
    implicit none

    logical        :: non_alive_neighbours
    type(Tregion)  :: reg
    type(Tpropagation_grid) :: grid
    integer        :: i1,i2,i3           ! identify a node (i,j,k or 0,interface,inode)
    integer        :: n,m,i,j,k,ii,jj,kk,icell ! local variables
    integer        :: centernode,n_concell
    type(Tinteger_coordinates)  :: concell(8)
    type(Tintersection),pointer  :: isec


    non_alive_neighbours = .false.

! store the identifiers of the centernode node in local variables. Reminder:
! for regular grid nodes i1,i2,i3  correspond to ir,ilat,ilong of the node
! for intersection nodes i1,i2,i3 are 0, intersection #, node # in intersection 

   i1 = reg%node(centernode)%i1 ; i2 = reg%node(centernode)%i2 ; i3 = reg%node(centernode)%i3


! make a list of grid cells of which the new alive node is part

   if (i1 /= 0) then   ! centernode is a regular grid node. Use spatial correlation of grid

      n_concell = 0
      do i=0,1
         if ((i1-i>0) .and. (i1-i<grid%nr)) then
            do j=0,1
               if ((i2-j>0) .and. (i2-j<grid%nlat)) then
                  do k=0,1
                     if ((i3-k>0) .and. (i3-k<grid%nlong)) then
                        n_concell=n_concell+1
                        concell(n_concell)%ir=i1-i  
                        concell(n_concell)%ilat=i2-j 
                        concell(n_concell)%ilong=i3-k
                     endif
                  end do
               end if
            end do
         end if
      end do

   else          ! centernode is an intersection node. Use the connections that have been found before

      if (i2 == reg%itop%id) then ; isec => reg%itop ; else ; isec => reg%ibot ; endif

      n_concell = 0
      do n=1,8
         icell = isec%ccell_from_inode(n,i3)
         if (icell > 0) then
            n_concell=n_concell+1
            concell(n_concell)%ir= isec%ccells(icell)%ir
            concell(n_concell)%ilat=isec%ccells(icell)%ilat 
            concell(n_concell)%ilong=isec%ccells(icell)%ilong
         endif
      enddo

   endif

   do n=1,n_concell

! find the intersection nodes of connected cell n

! explanation: each intersection has a 1-D list of cells cut by the interface, 
! and a list of intersection nodes that are part of each cut cell. Each regular 
! grid cell has a pointer ccind_from_3dc(i,j,k)%p (Cut Cell INDex FROM 3D Coordinates)
! associated with it, where p is a pointer to an integer array with as many elements 
! as there are intersections. If a cell is cut by interface n, the pointer 
! ccind_from_3dc(i,j,k)%p is allocated, and the variable ccind_from_3dc(i,j,k)%p(n) 
! contains the array index of cell (i,j,k) in the 1D cut cell list of intersection n 

! test whether the cell is cut by any interface, in that case the pointer to the local list of
! interfaces cutting the cell has been allocated

      if (associated(grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p)) then

! if so, check if the cell is cut by the top intersection

    ! icell is the index of the current connected cell in the list of cells cut by interface reg%itop.
                  
        icell=grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%itop%iface_id)

     ! if icell == 0 the cell is not cut be the top interface
    
         if(icell /= 0) then
 
            isec => reg%itop
            do jj=1,isec%n_inodes(icell)

!   m is the node number in the regional node list of node  jj in the list of inteface 
!   nodes that are part of cut cell icell
                m=isec%rbel_node_id(isec%inodes(jj,icell))

               if ( reg%node_status(m) < 0 ) then  
                  non_alive_neighbours = .true.
                  return
               endif

            end do

         end if


! then check the bottom intersection

        icell=grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%ibot%iface_id)
         if (icell /= 0) then

            isec => reg%ibot
            do jj=1,isec%n_inodes(icell)

               m=isec%rabo_node_id(isec%inodes(jj,icell))
!                  print *,'its regional node number is ',m,node_is_counted(m),centernode

               if ( reg%node_status(m) < 0 ) then  
                  non_alive_neighbours = .true.
                  return
               endif

            end do

         endif

      endif


! then find the regular grid nodes of connected cell n


      do i=0,1
         ii=concell(n)%ir+i
         do j=0,1
            jj=concell(n)%ilat+j
            do k=0,1
               kk=concell(n)%ilong+k

               ! reduced connectivity for regular nodes if cell is completely regular

               if ((i1 == 0 .or. (i1 /= 0 .and. abs(i1-ii)+abs(i2-jj)+abs(i3-kk) == 1))) then 
                  
                 if (grid%node_region(ii,jj,kk) == reg%id) then  
               ! node has to belong to the current region

                     m=grid%rnode_id(ii,jj,kk)

                     if ( reg%node_status(m) < 0 ) then  
                        non_alive_neighbours = .true.
                        return
                     endif

                  endif

               endif

            end do
         end do
      end do


   end do  ! loop over connected cells

   return

  end function non_alive_neighbours
!-----------------------------------------------------------------------------------------
! creates the list of neighbours (connected nodes) of a source

  subroutine get_source_neighbours(centernode,s,reg,grid)
    use mod_3dfm
    implicit none

    type(Tsource)  :: s
    type(Tregion)  :: reg
    type(Tpropagation_grid) :: grid
    integer        :: i1,i2,i3           ! identify a node (i,j,k or 0,interface,inode)
    integer        :: n,m,i,j,k,ii,jj,kk,icell ! local variables
    integer        :: centernode,n_concell,mstore(30)
    type(Tinteger_coordinates)  :: concell(8)
    type(Tintersection),pointer  :: isec

    

! store the identifiers of the centernode node in local variables. Reminder:
! for regular grid nodes i1,i2,i3  correspond to ir,ilat,ilong of the node
! for intersection nodes i1,i2,i3 are 0, intersection #, node # in intersection 

   i1 = reg%node(centernode)%i1 ; i2 = reg%node(centernode)%i2 ; i3 = reg%node(centernode)%i3


! make a list of grid cells of which the new alive node is part

   if (i1 /= 0) then   ! centernode is a regular grid node. Use spatial correlation of grid

      n_concell = 0
      do i=0,1
         if ((i1-i>0) .and. (i1-i<grid%nr)) then
            do j=0,1
               if ((i2-j>0) .and. (i2-j<grid%nlat)) then
                  do k=0,1
                     if ((i3-k>0) .and. (i3-k<grid%nlong)) then
                        n_concell=n_concell+1
                        concell(n_concell)%ir=i1-i  
                        concell(n_concell)%ilat=i2-j 
                        concell(n_concell)%ilong=i3-k
                     endif
                  end do
               end if
            end do
         end if
      end do

   else          ! centernode is an intersection node. Use the connections that have been found before

      if (i2 == reg%itop%id) then ; isec => reg%itop ; else ; isec => reg%ibot ; endif

      n_concell = 0
      do n=1,8
         icell = isec%ccell_from_inode(n,i3)
         if (icell > 0) then
            n_concell=n_concell+1
            concell(n_concell)%ir= isec%ccells(icell)%ir
            concell(n_concell)%ilat=isec%ccells(icell)%ilat 
            concell(n_concell)%ilong=isec%ccells(icell)%ilong
         endif
      enddo

   endif

   do n=1,n_concell

! find the intersection nodes of connected cell n


! test whether the cell is cut by any interface, in that case the pointer to the local list of
! interfaces cutting the cell has been allocated

      if (associated(grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p)) then

! if so, check if the cell is cut by the top intersection

  ! icell is the index of the current connected cell in the list of cells cut by interface reg%itop.
                  
        icell=grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%itop%iface_id)

     ! if icell == 0 the cell is not cut be the top interface
    
         if(icell /= 0) then
 
            isec => reg%itop
            do jj=1,isec%n_inodes(icell)

!   m is the node number in the regional node list of node  jj in the list of 
!  inteface nodes that are part of cut cell icell
                m=isec%rbel_node_id(isec%inodes(jj,icell))

               if ( s%n_cnode == 0 ) then  
                  s%n_cnode=1
                  mstore(s%n_cnode) = m
                  s%cnode(s%n_cnode)%i1=reg%node(m)%i1
                  s%cnode(s%n_cnode)%i2=reg%node(m)%i2
                  s%cnode(s%n_cnode)%i3=reg%node(m)%i3
               else
                  if ( count(mstore(1:s%n_cnode) == m) == 0 ) then
                     s%n_cnode=s%n_cnode+1
                     mstore(s%n_cnode) = m
                     s%cnode(s%n_cnode)%i1=reg%node(m)%i1
                     s%cnode(s%n_cnode)%i2=reg%node(m)%i2
                     s%cnode(s%n_cnode)%i3=reg%node(m)%i3
                  endif
               endif

            end do

         end if


! then check the bottom intersection

        icell=grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%ibot%iface_id)

         if (icell /= 0) then

            isec => reg%ibot
            do jj=1,isec%n_inodes(icell)

               m=isec%rabo_node_id(isec%inodes(jj,icell))

               if ( s%n_cnode == 0 ) then  
                  s%n_cnode=1
                  mstore(s%n_cnode) = m
                  s%cnode(s%n_cnode)%i1=reg%node(m)%i1
                  s%cnode(s%n_cnode)%i2=reg%node(m)%i2
                  s%cnode(s%n_cnode)%i3=reg%node(m)%i3
               else
                  if ( count(mstore(1:s%n_cnode) == m) == 0 ) then
                     s%n_cnode=s%n_cnode+1
                     mstore(s%n_cnode) = m
                     s%cnode(s%n_cnode)%i1=reg%node(m)%i1
                     s%cnode(s%n_cnode)%i2=reg%node(m)%i2
                     s%cnode(s%n_cnode)%i3=reg%node(m)%i3
                  endif
               endif

            end do

         endif

      endif


! then find the regular grid nodes of connected cell n

!      write(22,*) 'regular nodes of cell',n

      do i=0,1
         ii=concell(n)%ir+i
         do j=0,1
            jj=concell(n)%ilat+j
            do k=0,1
               kk=concell(n)%ilong+k

                 if (grid%node_region(ii,jj,kk) == reg%id) then  
                 ! node has to belong to the current region

                     m=grid%rnode_id(ii,jj,kk)

                     if ( s%n_cnode == 0 ) then  
                        s%n_cnode=1
                        mstore(s%n_cnode) = m
                        s%cnode(s%n_cnode)%i1=reg%node(m)%i1
                        s%cnode(s%n_cnode)%i2=reg%node(m)%i2
                        s%cnode(s%n_cnode)%i3=reg%node(m)%i3
                     else
                        if ( count(mstore(1:s%n_cnode) == m) == 0 ) then
                           s%n_cnode=s%n_cnode+1
                           mstore(s%n_cnode) = m
                           s%cnode(s%n_cnode)%i1=reg%node(m)%i1
                           s%cnode(s%n_cnode)%i2=reg%node(m)%i2
                           s%cnode(s%n_cnode)%i3=reg%node(m)%i3
                        endif
                     endif


                  endif


            end do
         end do
      end do


   end do  ! loop over connected cells

   return

  end subroutine get_source_neighbours

!****************************************************************************************************************
! this subroutine modifies the time gradient at an intersection to represent the reflected wave
! tgrad(reflected)= tgrad - 2*(tgrad.n)*n, where n is the normal to the interface (no conversion)
! or in the same way as a refraction if a wave type conversion takes place

  subroutine reflect_gradient(isec,tf_prev,vtype)
    use mod_3dfm
    implicit none

    type(Tintersection)  :: isec
    type(Ttime_field)    :: tf_prev
    integer  :: n,vtype,n_no_ref
    real(kind=dp),dimension(:),pointer :: vel
    real(kind=dp)        :: direction,det,grad_perp,grad_perp_reflected,grad_par(3)

    if (tf_prev%vtype == vtype) then  ! simple reflection ( no wave type conversion)

       do n=1,isec%nnode

          isec%time_gradient(1:3,n)=isec%time_gradient(1:3,n)- &
               2.0_dp*dot_product(isec%normal(1:3,n),isec%time_gradient(1:3,n))*isec%normal(1:3,n)

       end do

    else          ! reflection with wave type conversion

       if (isec%id == tf_prev%reg%itop%id) then   ! reflection off the top interface
          vel => isec%vel_bot(:,vtype)
          direction = -1.0_dp
       else                                       ! reflection off the bottom interface
          vel => isec%vel_top(:,vtype)
          direction = 1.0_dp
       endif
       n_no_ref=0


       do n=1,isec%nnode

          grad_perp=dot_product(isec%normal(1:3,n),isec%time_gradient(1:3,n))
          grad_par=isec%time_gradient(1:3,n)-grad_perp*isec%normal(1:3,n)

          det=1.d0/(vel(n)**2) - sum(grad_par**2)

          if (det > 0.0_dp) then
             
             ! the refracted ray exists

             grad_perp_reflected=sqrt(det)
             isec%time_gradient(1:3,n)=grad_par + &
                  sign(grad_perp_reflected,direction)*isec%normal(1:3,n)

          else

             ! no converted reflection possible

             isec%arrivaltime(n)=huge_time
             isec%time_gradient(1:3,n)=0.0_dp
             n_no_ref= n_no_ref+1

          endif

       end do

       if (n_no_ref > 0) print *,n_no_ref,' intersection points could not reflect the converted wave'
    
    endif




    return

  end subroutine reflect_gradient


!****************************************************************************************************
! this subroutine modifies the time gradient at an intersection to represent the refracted wave
! tgrad(parallel) =  tgrad - (tgrad.n)*n, where n is the normal to the interface 
! then solve |tgrad(perpendicular)|^2 + |tgrad(parallel)|^2 = s^2 where s is the slowness 
! in the region refracted into.
! also checks for total reflection and supresses such points as a source of refracted waves

  subroutine refract_gradient(isec,reg,vtype,direction)
    use mod_3dfm
    implicit none

    type(Tintersection)  :: isec
    type(Tregion)        ::reg

    integer  :: n,n_tot_ref,direction,vtype
    real(kind=dp)         :: grad_par(3),grad_perp,grad_perp_refracted,det,direction_real
    real(kind=dp),dimension(:),pointer  :: vel


    if (isec%id == reg%itop%id) vel => isec%vel_top(:,vtype)
    if (isec%id == reg%ibot%id) vel => isec%vel_bot(:,vtype)
    n_tot_ref=0
    direction_real = dble(direction) 

    do n=1,isec%nnode

       grad_perp=dot_product(isec%normal(1:3,n),isec%time_gradient(1:3,n))
       grad_par=isec%time_gradient(1:3,n)-grad_perp*isec%normal(1:3,n)

       det=1.d0/(vel(n)**2) - sum(grad_par**2)

       if (det > 0.0_dp) then
             
          ! the refracted ray exists

          grad_perp_refracted=sqrt(det)
          isec%time_gradient(1:3,n)=grad_par + &
               sign(grad_perp_refracted,direction_real)*isec%normal(1:3,n)

       else

          ! total reflection

          isec%arrivaltime(n)=huge_time
          isec%time_gradient(1:3,n)=0.0_dp
          n_tot_ref= n_tot_ref+1

       endif

    end do

    if (n_tot_ref > 0) print *,'total reflection occurred at ',n_tot_ref,' intersection nodes'

    return

  end subroutine refract_gradient

!*****************************************************************************************************
! this subroutine applies a correction to the time gradient at a node as derived from
! the direction of the wavefront used for the final update (for irregular updates only).
! we assume that the change in velocity between the location where the wavefront  used for 
! the update is evaluated and the velocity
! at the update node occurs as a step at a plane with a normal given by the norm of the
! velocity gradient. The corrected time gradient is the incoming time gradient refracted at this
! "interface". Experiments show that this significantly improves the constancy of the ray parameter
! along a ray in 1-D velocity distribution tests


  subroutine refract_locally(r,lat,long,gridv,tgrad)
    use mod_3dfm
    implicit none

    real(kind=dp)         :: grad_par(3),grad_perp,grad_perp_refracted,det,r,lat,long,vel
    real(kind=dp),dimension(3)  :: tgrad,vgrad,normal
    type(Tvelocity_grid)  :: gridv


    call velocity_gradient(r,lat,long,gridv,vgrad(1),vgrad(2),vgrad(3),vel)

    if (sum(vgrad**2) < 1.0e-30_dp) return  ! avoid division by zero if gradient = 0

    normal=vgrad/sqrt(sum(vgrad**2))

    grad_perp=dot_product(normal,tgrad)

    grad_par=tgrad-grad_perp*normal
    
    det=1.d0/(vel**2) - sum(grad_par**2)

    if (det > 0.0_dp) then
             
       ! the refracted ray exists

       grad_perp_refracted=sqrt(det)
       tgrad=grad_par + sign(grad_perp_refracted,grad_perp)*normal

    endif


    return

  end subroutine refract_locally
!*********************************************************************************
!--------------------------------------------------------------------------------\
!*******************************************************************************************
! This subroutine creates a triangualtion of a 2-D set of points
! It calls malcolm Sambridge's triangualtion routines that are a part
! of the natural neighbours package

      subroutine triangulate(t)

      use mod_3dfm

! all default floats are set to double

      implicit double precision (a-h,o-z)

! argument definition

      type(Ttriangulation) :: t

! local array definition

      double precision,dimension(:,:),allocatable::points,centres
      integer,dimension(:,:),allocatable::neighbour,SPfromTR
      integer,dimension(:),allocatable::vis_tlist,vis_elist,add_tlist,hulltriangles
      integer,dimension(:),allocatable::nnn,nnlist,ntrilist
      logical,dimension(:),allocatable::lt_work,ln_work
      double precision,dimension(:,:),allocatable::work_d1,work_d2,work_d3,work_d4
      double precision,dimension(:,:),allocatable::work_d5,work_d6,work_d7
      real,dimension(:,:),allocatable::work_r1
      real,dimension(:),allocatable::work_r2
      integer,dimension(:),allocatable::work_i1,work_i3
      integer,dimension(:,:),allocatable::work_i2


! other variables

      double precision eps
      integer dmode
      logical clockwise

!---------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! start the triangulation 
!------------------------------------------------------------------------------------


! nl is the initial number of points

      nl=t%npoints

! nlmax is the maximum number of points

      nlmax=nl

! ntmax is the maximum number of delaunay triangles built on the points

      ntmax=nlmax*3

! nhmax is the maximum number of points on the convex hull formed
! by the lagrangian cloud

      nhmax=nlmax

! npmax is nlmax

      npmax=nlmax

! nnpnmax is the maximum number of neighbours per point

      nnpnmax=50

! nvmax and nmax are the size of working arrays

      nvmax=ntmax
      nmax=3*ntmax+npmax

! eps is a small number

      eps=tiny(eps)

! allocates memory for the delaunay triangulation

     allocate (points(2,nlmax),centres(3,ntmax))
     allocate (neighbour(3,ntmax),SPfromTR(3,ntmax))
     allocate (vis_tlist(nvmax),vis_elist(nvmax),add_tlist(nvmax))
     allocate (hulltriangles(nhmax),nnn(npmax+1),nnlist(nmax),ntrilist(nmax))
     allocate (lt_work(ntmax),ln_work(nlmax))
     allocate (work_d1(2,nnpnmax),work_d2(2,nnpnmax),work_d3(2,nnpnmax))
     allocate (work_d4(2,nnpnmax),work_d5(2,nnpnmax),work_d6(2,nnpnmax))
     allocate (work_d7(2,nnpnmax),work_r1(nnpnmax,2),work_r2(nnpnmax))
     allocate (work_i1(nnpnmax),work_i2(2,nnpnmax),work_i3(nnpnmax))





!     print *,'before triangulation'
!--------------------------------------------------------------------------------------


! initialise variables for nn_setup

         dmode=-2
         nmode=0
         clockwise=.true.
         nohalt_hull=0
         loc=1
         SPfromTR=0
         neighbour=0
         points=0.d0
         field=0.d0
         centres=0.d0
         hulltriangles=0
         nnn=0
         nnlist=0
         ntrilist=0
         vis_tlist=0
         vis_elist=0
         add_tlist=0
         lt_work=.false.
         ln_work=.false.

         do i=1,nl
            points(1,i)=t%points(1,i)
            points(2,i)=t%points(2,i)
         enddo

! nn_setup (calculates delaunay triangles and other info)

         call nn2d_setup &
              (nl,ntmax,nhmax,npmax,nnpnmax,nmax, &
              points,dmode,nmode,clockwise,t%points(1,:),nt,SPfromTR, &
              centres,neighbour,nh,hulltriangles,nohalt_hull, &
              loc,nnn,nnlist,ntrilist, &
              eps,nvmax,vis_tlist,vis_elist,add_tlist, &
              lt_work,ln_work)



!  
         t%ntriangles = nt
         allocate(t%points_from_triangle(3,nt),t%triangle_neighbours(3,nt))
         do i=1,nt
            t%points_from_triangle(1:3,i)=SPfromTR(1:3,i)
            t%triangle_neighbours(1:3,i)=neighbour(1:3,i)
         end do


       ! the inverse connectivity

         allocate(t%n_triangles_from_point(t%npoints))
         t%n_triangles_from_point = 0
         do i=1,nt
            do j=1,3
               t%n_triangles_from_point(SPfromTR(j,i)) = t%n_triangles_from_point(SPfromTR(j,i)) + 1
            end do
         end do

         nmax = maxval(t%n_triangles_from_point)
         allocate(t%triangles_from_point(nmax,t%npoints))
         t%n_triangles_from_point = 0
         do i=1,nt
            do j=1,3
               ipoint = SPfromTR(j,i)
               t%n_triangles_from_point(ipoint) = t%n_triangles_from_point(ipoint) + 1
               t%triangles_from_point(t%n_triangles_from_point(ipoint),ipoint) = i
            end do
         end do




! deallocate everything


      deallocate (points,centres)
      deallocate (neighbour,SPfromTR)
      deallocate (vis_tlist,vis_elist,add_tlist)
      deallocate (hulltriangles,nnn,nnlist,ntrilist)
      deallocate (lt_work,ln_work)
      deallocate (work_d1,work_d2,work_d3,work_d4,work_d5,work_d6,work_d7)
      deallocate (work_r1,work_r2,work_i1,work_i2,work_i3)


         return

       end subroutine triangulate

!****************************************************************************

subroutine clean_triangulation(tri)
  use mod_3dfm

  type(Ttriangulation) ::tri

  if (associated(tri%points_from_triangle)) deallocate(tri%points_from_triangle)
  if (associated(tri%triangle_neighbours)) deallocate(tri%triangle_neighbours)
  if (associated(tri%points)) deallocate(tri%points)
  if (associated(tri%triangles_from_point)) deallocate(tri%triangles_from_point)
  if (associated(tri%n_triangles_from_point)) deallocate(tri%n_triangles_from_point)

  return
end subroutine clean_triangulation

!**************************************************************************************
! This subroutine sets a flag at all nodes of the propagation grid that are not
! connected to any cut cell. This is used to speed up the fast marching calculations
! 
subroutine tag_regular_nodes(grid)
  use mod_3dfm

  type(Tpropagation_grid)  :: grid

  integer  :: i,j,k
  logical  :: not_regular

  if (.not. associated(grid%ccind_from_3dc)) &
       stop 'nodes cannot be tagged when ccind_from_3dc not initialized'

  allocate(grid%fully_regular(grid%nr,grid%nlat,grid%nlong))

  grid%fully_regular = .false.

  do k=2, grid%nlong-1
     do j=2,grid%nlat-1
        do i=2,grid%nr-1

           not_regular = associated(grid%ccind_from_3dc(i-1,j-1,k-1)%p ) .or. &
                associated(grid%ccind_from_3dc(i  ,j-1,k-1)%p ) .or. &
                associated(grid%ccind_from_3dc(i-1,j  ,k-1)%p ) .or. &
                associated(grid%ccind_from_3dc(i  ,j  ,k-1)%p ) .or. &
                associated(grid%ccind_from_3dc(i-1,j-1,k  )%p ) .or. &
                associated(grid%ccind_from_3dc(i  ,j-1,k  )%p ) .or. &
                associated(grid%ccind_from_3dc(i-1,j  ,k  )%p ) .or. &
                associated(grid%ccind_from_3dc(i  ,j  ,k  )%p )

           grid%fully_regular(i,j,k) = .not. not_regular

        end do
     end do
  end do

end subroutine tag_regular_nodes

!********************************************************************************
! This subroutine transfers the nodes from a region in the refined source grid to 
! the corresponding region on the main propagation grid and sets the node status of the
! nodes in the propagation grid region that received a value to 0 (final value reached)

subroutine transfer_refined_region(sreg,reg,t_short)

use mod_3dfm
implicit none

type(Tregion) :: sreg,reg
real(kind=dp) :: t_short

type(Tintersection),pointer :: isec

integer :: k,m,n,i1,i2,i3,j1,j2,j3,icell,kk,intype
real(kind=dp)  :: xir,xilat,xilong,dxir,dxilat,dxilong,dist
integer  :: nrmin,nlatmin,nlongmin
!logical  :: pdiag = .false.


nrmin = sgrid%index_r0
nlatmin = sgrid%index_lat0
nlongmin = sgrid%index_long0

tloop : do n=1,sreg%nnode

!   pdiag = (reg%id == 3) .and. (n == 359764)

!   if (pdiag) print *,'trying to transfer xnode'

   if (sreg%arrivaltime(n) <= t_short) then  ! transfer only times less than the first boundary hit

!      if (pdiag) print *,'passed t-short'

   i1=sreg%node(n)%i1  ; i2=sreg%node(n)%i2   ; i3= sreg%node(n)%i3
   

   if (i1 /= 0) then           ! the regional node is a regular node

!      if (pdiag) print *,'xnode is regular'

      if (mod(i1-1,refinement_factor) == 0 .and. mod(i2-1,refinement_factor) == 0 &
           .and. mod(i3-1,refinement_factor) == 0) then

         j1 = nrmin    +   (i1-1)/refinement_factor
         j2 = nlatmin  +   (i2-1)/refinement_factor
         j3 = nlongmin +   (i3-1)/refinement_factor

!         print *,n,'regular node',j1,j2,j3
!         print *,pgrid%rnode_id(j1,j2,j3), reg%nnode

         reg%arrivaltime(pgrid%rnode_id(j1,j2,j3)) = sreg%arrivaltime(n)
         reg%time_gradient(1:3,pgrid%rnode_id(j1,j2,j3)) = sreg%time_gradient(1:3,n)
         reg%node_status(pgrid%rnode_id(j1,j2,j3)) = 0

!         print *,n,'regular node'

      endif

   else   ! the regional node is an intersection node

!      if (pdiag) print *,'irregular sreg node',i1,i2,i3,sintersection(i2)%intype(i3)

! choose the corresponding intersection on the main grid

      if (i2 == sreg%itop%id) then
         isec => reg%itop
      else
         isec => reg%ibot
      endif

! coarse grid index coordinates

      xir      = (sintersection(i2)%r(i3)-pgrid%r0)/pgrid%dr0
      xilat    = (sintersection(i2)%lat(i3)-pgrid%lat0)/pgrid%dlat0
      xilong   = (sintersection(i2)%long(i3)-pgrid%long0)/pgrid%dlong0

!      if (pdiag) print '(a10,3f15.6)','gi cs',xir,xilat,xilong

! distance from nearest coarse grid coordinate plane

      dxir    = abs(xir - anint(xir))
      dxilat  = abs(xilat - anint(xilat))
      dxilong = abs(xilong - anint(xilong))

!      if (pdiag) print '(a10,3f15.6)','npdist',dxir,dxilat,dxilong


! coordinates of the coarse grid cell containing the refined region node

      j1= min(int(xir)+1,pgrid%nr-1)
      j2= min(int(xilat)+1,pgrid%nlat-1)
      j3= min(int(xilong)+1,pgrid%nlong-1)

!      if (pdiag) print *,'cgcell',j1,j2,j3


! convert type of inode on refined grid to type of inode on coarse grid

      if (sintersection(i2)%intype(i3) /= 0) then

         intype = sintersection(i2)%intype(i3)

      else

         if (dxir<interface_tolerance.and.dxilat<interface_tolerance.and.&
              dxilong<interface_tolerance) then
            intype = 0
         else if(dxilat < interface_tolerance .and. dxilong < interface_tolerance) then
            intype = 1
         else if(dxir < interface_tolerance .and. dxilong < interface_tolerance) then
            intype = 2
         else if(dxir < interface_tolerance .and. dxilat < interface_tolerance) then
            intype = 3
         else
            cycle tloop
         endif

      endif

! test whether the refined regional node lies on r,lat or long connection

      select case (intype)

         case(0)

            ! if all coordinates coincide exactly with the coarse grid lines, there is 
            ! a corresponding coarse grid intersection node

            
            if (dxir < interface_tolerance .and. dxilat < interface_tolerance &
                 .and. dxilong < interface_tolerance ) then 

!               if (pdiag) print *,n,'type 0'

               ! get the cut cell index of the coarse grid cell
               if (associated(pgrid%ccind_from_3dc(j1,j2,j3)%p)) then

                  icell = pgrid%ccind_from_3dc(j1,j2,j3)%p(isec%iface_id)

       ! test the nodes of this cut cell for coincidence with the refined inode under consideration
                  do m=1,isec%n_inodes(icell)
                  
                     k=isec%inodes(m,icell)  ! the intersection node #

                     if (i2 == sreg%itop%id) then 
                        kk=isec%rbel_node_id(k)
                     else
                        kk=isec%rabo_node_id(k)
                     endif


                     if (isec%intype(k) == 0) then  
              ! only consider nodes of the same type as the refined inode

                        dist =sqrt((sintersection(i2)%r(i3)-isec%r(k))**2 &
                             + (isec%r(k)*(sintersection(i2)%lat(i3)-isec%lat(k)))**2 &
                             + (isec%r(k)*isec%coslat(k)*(sintersection(i2)%long(i3) &
                             -isec%long(k)))**2)

                        if (dist < pgrid%tolerance) then
                           reg%arrivaltime(kk) = sreg%arrivaltime(n)
                           reg%time_gradient(1:3,kk) = sreg%time_gradient(1:3,n)
                           reg%node_status(kk) = 0
                           cycle tloop
                        endif

                     endif

                  end do

               endif

            endif


         case(1)    ! inode lies on an r-connection

            ! if other coordinates coincide exactly with the coarse grid lines, there is 
            ! a corresponding coarse grid intersection node

            
            if (dxilat < interface_tolerance .and. dxilong < interface_tolerance) then 

!               if (pdiag) print *,n,'type 1'

               ! get the cut cell index of the coarse grid cell
               if (associated(pgrid%ccind_from_3dc(j1,j2,j3)%p)) then

                  icell = pgrid%ccind_from_3dc(j1,j2,j3)%p(isec%iface_id)

       ! test the nodes of this cut cell for coincidence with the refined inode under consideration
                  do m=1,isec%n_inodes(icell)
                  
                     k=isec%inodes(m,icell)  ! the intersection node #

                     if (i2 == sreg%itop%id) then 
                        kk=isec%rbel_node_id(k)
                     else
                        kk=isec%rabo_node_id(k)
                     endif


                     if (isec%intype(k) == 1) then  
                ! only consider nodes of the same type as the refined inode

                        dist =sqrt((sintersection(i2)%r(i3)-isec%r(k))**2 &
                             +(isec%r(k)*(sintersection(i2)%lat(i3)-isec%lat(k)))**2 &
                             + (isec%r(k)*isec%coslat(k)*(sintersection(i2)%long(i3) &
                             -isec%long(k)))**2)


                        if (dist < pgrid%tolerance) then
                           reg%arrivaltime(kk) = sreg%arrivaltime(n)
                           reg%time_gradient(1:3,kk) = sreg%time_gradient(1:3,n)
                           reg%node_status(kk) = 0
                           cycle tloop
                        endif

                     endif

                  end do

               endif

            endif


         case(2) ! inode lies on a lat-connection

            if (dxir < interface_tolerance .and. dxilong < interface_tolerance) then

!               if (pdiag) print *,n,'type 2'
               if (associated(pgrid%ccind_from_3dc(j1,j2,j3)%p)) then

                  icell = pgrid%ccind_from_3dc(j1,j2,j3)%p(isec%iface_id)

                  do m=1,isec%n_inodes(icell)
                  
                     k=isec%inodes(m,icell)

                     if (i2 == sreg%itop%id) then 
                        kk=isec%rbel_node_id(k)
                     else
                        kk=isec%rabo_node_id(k)
                     endif

                     if (isec%intype(k) == 2) then

                        dist =sqrt((sintersection(i2)%r(i3)-isec%r(k))**2 &
                             +(isec%r(k)*(sintersection(i2)%lat(i3)-isec%lat(k)))**2 &
                             + (isec%r(k)*isec%coslat(k)*(sintersection(i2)%long(i3) &
                             -isec%long(k)))**2)

                        if (dist < pgrid%tolerance) then
                           reg%arrivaltime(kk) = sreg%arrivaltime(n)
                           reg%time_gradient(1:3,kk) = sreg%time_gradient(1:3,n)
                           reg%node_status(kk) = 0
                           cycle tloop
                        endif

                     endif

                  end do

               endif

            endif


         case(3)    ! inode lies on a long-connection

            if (dxilat < interface_tolerance .and. dxir < interface_tolerance) then

!               if (pdiag) print *,n,'type 3'

               if (associated(pgrid%ccind_from_3dc(j1,j2,j3)%p)) then

                  icell = pgrid%ccind_from_3dc(j1,j2,j3)%p(isec%iface_id)

                  do m=1,isec%n_inodes(icell)
                  
                     k=isec%inodes(m,icell)


                     if (i2 == sreg%itop%id) then 
                        kk=isec%rbel_node_id(k)
                     else
                        kk=isec%rabo_node_id(k)
                     endif

                     if (isec%intype(k) == 3) then

                        dist =sqrt((sintersection(i2)%r(i3)-isec%r(k))**2 &
                             +(isec%r(k)*(sintersection(i2)%lat(i3)-isec%lat(k)))**2 &
                             + (isec%r(k)*isec%coslat(k)*(sintersection(i2)%long(i3) &
                             -isec%long(k)))**2)

                        if (dist < pgrid%tolerance) then
                           reg%arrivaltime(kk) = sreg%arrivaltime(n)
                           reg%time_gradient(1:3,kk) = sreg%time_gradient(1:3,n)
                           reg%node_status(kk) = 0
                           cycle tloop
                        endif

                     endif

                  end do

               endif

            endif

      end select


   endif   ! refined node is an inode

   endif  ! arrivaltime < t_short

end do tloop

end subroutine transfer_refined_region

!*****************************************************************************************************
!*****************************************************************************************************
subroutine write_valid_rays(n,m)

  use mod_3dfm
  type(Tray),pointer                   :: ray
  integer                              :: i,j,k,m,n,zero

  zero=0

  if (receiver(n)%ray(m)%is_multiray) then

     do k=1,receiver(n)%ray(m)%n_subrays

        ray => receiver(n)%ray(m)%subray(k)

        if (ray%valid) then

           write(31,'(5i6)') n,ray%source_id,m,k,ray%nsections

           do i=1,ray%nsections

              write(31,'(2i6,2l5)') ray%section(i)%npoints,ray%section(i)%reg%id, &
                   ray%section(i)%diffracted,ray%section(i)%headwave

              do j=ray%section(i)%npoints,1,-1
                 write(31,'(3f17.8)') ray%section(i)%point(1:3,j)
              end do

           end do

        else

           write(31,'(5i6)') n,ray%source_id,m,k,zero

        endif

     end do

  else

     k = 0

     ray => receiver(n)%ray(m)

     if (ray%valid) then

        write(31,'(5i6)') n,ray%source_id,m,zero,ray%nsections

        do i=1,ray%nsections

           write(31,'(2i6,2l5)') ray%section(i)%npoints,ray%section(i)%reg%id, &
                ray%section(i)%diffracted,ray%section(i)%headwave

           do j=ray%section(i)%npoints,1,-1
              write(31,'(3f17.8)') ray%section(i)%point(1:3,j)
           end do

        end do

     else

        write(31,'(5i6)') n,ray%source_id,m,zero,zero

     endif

  endif


end subroutine write_valid_rays


!*****************************************************************************************************
!*****************************************************************************************************
subroutine clean_ray(n,m)

  use mod_3dfm

  integer                              :: i,k,m,n
  type(Tray),pointer                   :: ray

  if (receiver(n)%ray(m)%is_multiray) then

     do k=1,receiver(n)%ray(m)%n_subrays

        ray => receiver(n)%ray(m)%subray(k)

        if (associated(ray%pdev)) deallocate(ray%pdev)
        if (associated(ray%pdev_indx)) deallocate(ray%pdev_indx)

        do i=1,ray%nsections
           if (associated(ray%section(i)%point)) deallocate(ray%section(i)%point)           
        end do

     end do

  else

     ray => receiver(n)%ray(m)

     if (associated(ray%pdev)) deallocate(ray%pdev)
     if (associated(ray%pdev_indx)) deallocate(ray%pdev_indx)

     do i=1,ray%nsections
        if (associated(ray%section(i)%point)) deallocate(ray%section(i)%point)
     end do

  endif

end subroutine clean_ray

!*****************************************************************************************************
subroutine write_frechet_derivatives(n,m)

  use mod_3dfm

  integer                              :: i,k,m,n,zero
  type(Tray),pointer                   :: ray

  zero=0

  if (receiver(n)%ray(m)%is_multiray) then

     do k=1,receiver(n)%ray(m)%n_subrays

        ray => receiver(n)%ray(m)%subray(k)

        if (ray%valid) then

           write(21,'(5i6)') n,ray%source_id,m,k,ray%n_pdev

           do i=1,ray%n_pdev
              write(21,'(i10,f17.8)') ray%pdev_indx(i),ray%pdev(i)
           end do

        else

           write(21,'(5i6)') n,ray%source_id,m,k,zero

        endif

     end do

  else

     k = 0

     ray => receiver(n)%ray(m)

     if (ray%valid) then

        write(21,'(5i6)') n,ray%source_id,m,k,ray%n_pdev

        do i=1,ray%n_pdev
           write(21,'(i10,e17.8)') ray%pdev_indx(i),ray%pdev(i)
        end do

     else

        write(21,'(5i6)') n,ray%source_id,m,k,zero

     endif

  endif

end subroutine write_frechet_derivatives
!
!**************************************************************
!

subroutine load_source_timefields(s)

  use mod_3dfm
  implicit none

  type(Tsource)  :: s
  integer   :: n,nfile

  nfile=s%nfile
  open(nfile,form='unformatted')

   do n=1,s%n_time_fields

      allocate(s%time_field(n)%arrivaltime(s%time_field(n)%reg%nnode))
      allocate(s%time_field(n)%time_gradient(3,s%time_field(n)%reg%nnode))
      read(nfile) s%time_field(n)%arrivaltime
      read(nfile) s%time_field(n)%time_gradient

   end do

   close(nfile)

end subroutine load_source_timefields

subroutine clean_source_timefields(s)

  use mod_3dfm
  implicit none

  type(Tsource)  :: s
  integer  :: n

  do n=1,s%n_time_fields

     deallocate(s%time_field(n)%arrivaltime,s%time_field(n)%time_gradient)

  end do

end subroutine clean_source_timefields

!
!*******************************************************************************
!
 subroutine write_arrivaltime_grid(src,path)

  use mod_3dfm
  implicit none

  type(Tpath) :: path
  type(Tsource) :: src
  integer       :: n,i,i1,i2,i3,ir,ilat,ilong,ntype
  type(Ttime_field),pointer :: tf
  type(Tregion),pointer :: reg
  real(kind=dp) :: r,lat,long,dxr,dxlat,dxlong


  if (associated(pgrid%arrivaltime)) deallocate(pgrid%arrivaltime)

  allocate(pgrid%arrivaltime(pgrid%nr,pgrid%nlat,pgrid%nlong))
  pgrid%arrivaltime=huge_time

  do n=path%first_tf_to_save, path%n_tf

     tf=> src%time_field(path%tf_sequence(n))
     reg => tf%reg

     do i=1, reg%nnode

        i1=reg%node(i)%i1 ; i2=reg%node(i)%i2 ; i3=reg%node(i)%i3 

        if (i1 /= 0) then
           pgrid%arrivaltime(i1,i2,i3)=min(tf%arrivaltime(i),pgrid%arrivaltime(i1,i2,i3))
        endif
        
        if (i1==0) then

           r=reg%r(i)-pgrid%r0 ; lat=reg%lat(i)-pgrid%lat0 ; long=reg%long(i)-pgrid%long0

           ir=nint(r/pgrid%dr0)+1
           ilat=nint(lat/pgrid%dlat0)+1
           ilong=nint(long/pgrid%dlong0)+1

           dxr=abs(r-pgrid%dr0*(ir-1))
           dxlat=reg%r(i)*abs(lat-pgrid%dlat0*(ilat-1))          
           dxlong=reg%r(i)*cos(reg%lat(i))*abs(long-pgrid%dlong0*(ilong-1)) 

           ntype=intersection(i2)%intype(i3)
           select case(ntype)
              case(0)
                 pgrid%arrivaltime(ir,ilat,ilong)= &
                      min(tf%arrivaltime(i),pgrid%arrivaltime(ir,ilat,ilong))

              case(1)
                 if (dxr<=pgrid%tolerance) pgrid%arrivaltime(ir,ilat,ilong)= &
                      min(tf%arrivaltime(i),pgrid%arrivaltime(ir,ilat,ilong))

              case(2)
                 if (dxlat<=pgrid%tolerance) pgrid%arrivaltime(ir,ilat,ilong)= &
                      min(tf%arrivaltime(i),pgrid%arrivaltime(ir,ilat,ilong))

              case(3)
                 if (dxlong<=pgrid%tolerance) pgrid%arrivaltime(ir,ilat,ilong)= &
                      min(tf%arrivaltime(i),pgrid%arrivaltime(ir,ilat,ilong))

           end select

        endif

     end do

  end do

  where (pgrid%arrivaltime>1.e10) pgrid%arrivaltime=-1.0


  write(19,*) src%id,path%id,path%first_tf_to_save
  do i3=1,pgrid%nlong
     do i2=1,pgrid%nlat
        do i1=1,pgrid%nr
           write(19,'(f12.5)') pgrid%arrivaltime(i1,i2,i3)
        end do
     end do
  end do


  deallocate(pgrid%arrivaltime)


end subroutine write_arrivaltime_grid
!
ellip.f
c	ELLIPTICITY CORRECTIONS - Polynomial interpolation
c
      subroutine ellip()
      real sc0
      real sc1
      real sc2
c SQUARE ROOT OF 3. TIMES .5
      real s3
c EPICENTAL CO-LATITUDE - RADIANS
      real ecolat
c EPICENTRAL DISTANCE - RADIANS
      real edist
c AZIMUTH FROM EPICENTER TO RECEIVER - RADIANS
      real az
c EPICENTRAL DEPTH (km)
      real edepth
c ELLIPTICITY CORRECTION -OUTPUT
      real tcor
      character*(*) phase
c TAU's
      real t0, t1, t2
c	CONSTANTS FOR POLYNOMIAL
      real t0con(8,10), t1con(8,10), t2con(8,10)
      integer ii, j
      real adepth
      data (t0con(1,j),j=1,10) /-0.01711,-1.7791,0.,0.,0.,-0.9630,-13.
     &2326,13.7390,0.,0./
      data (t0con(2,j),j=1,10) /-0.08291,-2.1455,2.4538,-0.7907,0.,2.
     &0258,-12.9357,2.1287,5.2668,-0.9229/
      data (t0con(3,j),j=1,10) /-1.5022,-0.0943,1.9655,-1.1661,.1393,3.
     &4920,-9.9051,-0.3875,5.3581,-0.0686/
      data (t0con(4,j),j=1,10) /2.9971,-2.9549,0.4082,0.,0.,28.1650,9.
     &2160,-17.9030,-5.2995,3.2029/
      data (t0con(5,j),j=1,10) /3.6775,-2.2221,0.,0.,0.,-1.3127,-6.2476,
     &1.6684,0.,0./
      data (t0con(6,j),j=1,10) /-10.6238,15.4993,-7.4840,1.0673,0.,3.
     &2763,-6.4596,-0.4923,0.,0./
      data (t0con(7,j),j=1,10) /-0.01332,-3.2777,-1.2243,7.5246,0.,-3.
     &4856,-10.3187,43.4834,-70.5341,-50.2287/
      data (t0con(8,j),j=1,10) /-0.07859,-4.0924,4.6116, -1.4760,0.,2.
     &9104,-17.8661, 4.6262,7.1486,-1.9154/
      data (t1con(1,j),j=1,10) /.0040,-0.7841,6.0441,-17.5535,0.,-0.
     &2549,2.0519,-19.0605,-37.8235,54.5110/
      data (t1con(2,j),j=1,10) /-.0048, .0839,-2.2705,2.4137,-0.5957,-2.
     &4241,-4.2792,1.9728,3.5644,-0.5285/
      data (t1con(3,j),j=1,10) /.0033,-1.3485,0.1735,1.1583,-0.4162,-0.
     &1096,0.2576,-0.5978,0.1888,0.1600/
      data (t1con(4,j),j=1,10) /2.6249,-.0025,-0.2086,-0.0184,0.,-1.
     &5077,0.9904,0.3513,0.,0./
      data (t1con(5,j),j=1,10) /3.4213,-0.9359,0.,0.,0.,0.,0.,0.,0.,0./
      data (t1con(6,j),j=1,10) /-8.0633,8.0238,-1.7407,0.,0.,0.,0.,0.,0.
     &,0./
      data (t1con(7,j),j=1,10) /0.0109,-1.2300,8.9145,-27.5847,0.,-0.
     &6951,5.6201,-33.0908,-83.8233,102.4333/
      data (t1con(8,j),j=1,10) /-0.0311,0.1896,-4.0694,4.2599,-1.0387,-
     &3.9368,-8.4379,2.6814,6.9535,-0.6086/
      data (t2con(1,j),j=1,10) /0.0107,0.0275,-0.6912,0.0347,0.1157,-0.
     &1836,0.,0.0296,0.,0./
      data (t2con(2,j),j=1,10) /0.0107,0.0275,-0.6912,0.0347,0.1157,-0.
     &1836,0.,0.0296,0.,0./
      data (t2con(3,j),j=1,10) /0.0005,-0.01231,-1.0156,0.4396,0.,0.,0.,
     &0.,0.,0./
      data (t2con(4,j),j=1,10) /-3.5838,2.1474,-0.3548,0.,0.,-1.3369,-5.
     &4889,0.6809,1.5096,-0.0763/
      data (t2con(5,j),j=1,10) /-2.9912,1.0313,0.,0.,0.,0.,0.,0.,0.,0./
      data (t2con(6,j),j=1,10) /3.2814,-7.1224,3.5418,-0.5115,0.,0.,0.,
     &0.,0.,0./
      data (t2con(7,j),j=1,10) /0.00025,0.1685,-2.2435,3.3433,0.,-0.
     &0503,0.5353,1.5362,-14.3118,-3.2938/
      data (t2con(8,j),j=1,10) /0.0843,-0.2917,-0.6767,-0.2934,0.2779,-
     &0.4336,0.0306,0.07113,0.,0./
      data s3 /.8660254/
      save sc0, sc1, sc2
c
c	INITIAL CALL TO SET UP CONSTANTS
      entry ellref(ecolat)
c
      sc0 = .25*(1.+3.*cos(2.*ecolat))
      sc1 = s3 * sin(2.*ecolat)
      sc2 = s3 * sin(ecolat)*sin(ecolat)
      return
c
c	CALLED ONCE FOR EACH STATION - RETURNS ELLIPTICITY CORRECTION IN tcor
      entry ellcor(edist, azim, edepth, phase,tcor)
      adepth = edepth/6371.
c DETERMINE INDEX FOR POLYNOMIAL
      if(.not.(phase .eq. 'P'))goto 23000
         if(.not.(edist .lt. (15.*(3.14159265/180.))))goto 23002
            ii = 1
            goto 23003
c        else
23002       continue
            ii = 2
23003    continue
         goto 23001
c     else
23000    continue
         if(.not.(phase .eq. 'PcP'))goto 23004
            ii = 3
            goto 23005
c        else
23004       continue
            if(.not.(phase .eq. 'PKPab'))goto 23006
               ii = 4
               goto 23007
c           else
23006          continue
               if(.not.(phase .eq. 'PKPbc'))goto 23008
                  ii = 5
                  goto 23009
c              else
23008             continue
                  if(.not.(phase .eq. 'PKIKP'))goto 23010
                     ii = 6
                     goto 23011
c                 else
23010                continue
                     if(.not.(phase .eq. 'S'))goto 23012
                        if(.not.(edist .lt. (15.*(3.14159265/180.))))
     &                    goto 23014
                           ii = 7
                           goto 23015
c                       else
23014                      continue
                           ii = 8
23015                   continue
                        goto 23013
c                    else
23012                   continue
                        tcor = 0.
                        return
23013                continue
23011             continue
23009          continue
23007       continue
23005    continue
23001 continue
c	COMPUTE TAU's
c      write(6,*) ii
      t0 = t0con(ii,1) + edist*(t0con(ii,2) + edist*(t0con(ii,3) + 
     &edist*(t0con(ii,4) + edist*t0con(ii,5)))) +adepth*(t0con(ii,6) + 
     &adepth*t0con(ii,7)) + adepth*edist*(t0con(ii,8) + t0con(ii,9)*
     &adepth + t0con(ii,10)*edist)
      t1 = t1con(ii,1) + edist*(t1con(ii,2) + edist*(t1con(ii,3) + 
     &edist*(t1con(ii,4) + edist*t1con(ii,5)))) +adepth*(t1con(ii,6) + 
     &adepth*t1con(ii,7)) + adepth*edist*(t1con(ii,8) + t1con(ii,9)*
     &adepth + t1con(ii,10)*edist)
      t2 = t2con(ii,1) + edist*(t2con(ii,2) + edist*(t2con(ii,3) + 
     &edist*(t2con(ii,4) + edist*t2con(ii,5)))) +adepth*(t2con(ii,6) + 
     &adepth*t2con(ii,7)) + adepth*edist*(t2con(ii,8) + t2con(ii,9)*
     &adepth + t2con(ii,10)*edist)
      tcor = sc0 * t0 + sc1 * cos(azim) * t1 + sc2 * cos(2.*azim) * t2
      return
      end


C ELLIPTICITY CORRECTIONS FOR AK135 MODEL (full set of phases)
C
       SUBROUTINE kellip()
C==========================================================================
C                                                                         
C    Ellipticity correction for any given phase using
C    Dziewonski & Gilbert representation
C                                                   
C      The ellipticity corrections are found by linear interpolation       
C    in terms of values calculated for the ak135 model for a wide 
C    range of phases to match the output of the iasp software 
C
C     first call:  ellref(ecolat) 
C                        - to set up source dependent constants
C     2nd call  :  ellcor(phase,edist,depth,ecolat,azim,tcor,abrt) 
C                        - to calculate correction for a station                                                                                     C                                                                         
C    Parameters: 
C    character  
C          phase : a  string specifying the PHASE,   -e.g P, ScP etc.  
C                                                        
C    real 
C          edist  :  epicentral distance to station (in degrees)     
C          edepth :  depth of event         
C          ecolat :  epicentral co-latitude of source (in radians) 
C          azim   :  azimuth from source to station (in radians)
C                                
C          tcor   :  time correction for path to allow for ellipticity
C 
C    logical 
C          abrt   :  a logical variable -usally set to .FALSE.  
C                    which is set to .TRUE. if a phase for      
C                    which no data is available is chosen       
C                                                                         
C==========================================================================
C   B.L.N. Kennett RSES,ANU        May 1995, August 1996                 
C   (based on earlier routine by D.J. Brown)
C   with input from W. Spakman, Utrecht
C=========================================================================
      character *(*) phase
      character*8 phcod(57)
      integer phind(57),phspn(57),phnch(57)
      real edist,edepth,ecolat,azim,
     ^     sc0,sc1,sc2,s3,tcor,
     ^     tau0, a0,b0,h0,d0,e0,f0,g0,
     ^     tau1, a1,b1,h1,d1,e1,f1,g1,
     ^     tau2, a2,b2,h2,d2,e2,f2,g2
      real dpth(6),delta(50)
      real t0(50,6),t1(50,6),t2(50,6)
      integer Ne,Nd
      logical abrt
      data phcod/
     & "Pup   ","P     ","Pdiff ","PKPab ","PKPbc ","PKPdf ",
     & "PKiKP ","pP    ","pPKPab","pPKPbc","pPKPdf","pPKiKP",
     & "sP    ","sPKPab","sPKPbc","sPKPdf","sPKiKP","PcP   ",
     & "ScP   ","SKPab ","SKPbc ","SKPdf ","SKiKP ","PKKPab",
     & "PKKPbc","PKKPdf","SKKPab","SKKPbc","SKKPdf","PP    ",
     & "P'P'  ","Sup   ","S     ","Sdiff ","SKSac ","SKSdf ",
     & "pS    ","pSKSac","pSKSdf","sS    ","sSKSac","sSKSdf",
     & "ScS   ","PcS   ","PKSab ","PKSbc ","PKSdf ","PKKSab",
     & "PKKSbc","PKKSdf","SKKSac","SKKSdf","SS    ","S'S'  ",
     & "SP    ","PS    ","PnS   "/
      data phind/
     &        1,      14,      91,     136,     165,     178,
     &      235,     364,     433,     462,     475,     532,
     &      661,     742,     771,     784,     841,     970,
     &     1047,    1100,    1113,    1134,    1195,    1316,
     &     1337,    1382,    1507,    1516,    1573,    1702,
     &     1827,    1932,    1945,    2022,    2067,    2132,
     &     2197,    2234,    2295,    2356,    2425,    2490,
     &     2551,    2628,    2681,    2694,    2711,    2772,
     &     2781,    2838,    2967,    3140,    3273,    3398,
     &     3587,    3656,    3697/
      data phspn/
     &        3,      19,      11,       7,       3,      14,
     &       32,      17,       7,       3,      14,      32,
     &       20,       7,       3,      14,      32,      19,
     &       13,       3,       5,      15,      30,       5,
     &       11,      31,       2,      14,      32,      31,
     &       26,       3,      19,      11,      16,      16,
     &        9,      15,      15,      17,      16,      15,
     &       19,      13,       3,       4,      15,       2,
     &       14,      32,      43,      33,      31,      47,
     &       17,      10,       6/ 
      data phnch/
     &        3,       1,       5,       5,       5,       5,
     &        5,       2,       6,       6,       6,       6,
     &        2,       6,       6,       6,       6,       3,
     &        3,       5,       5,       5,       5,       6,
     &        6,       6,       6,       6,       6,       2,
     &        4,       3,       1,       5,       5,       5,
     &        2,       6,       6,       2,       6,       6,
     &        3,       3,       5,       5,       5,       6,
     &        6,       6,       6,       6,       2,       4,
     &        2,       2,       3/ 
      data dpth/ 0.0, 100.0, 200.0, 300.0, 500.0, 700.0 /
      save sc0,sc1,sc2
c...
c     In addition to the phase names listed above a number of phase aliases
c     are available in the routine phase_alias, e.g. Pn --> P etc
c     The input phase code is first checked against the phcod array
c     and next against the phase aliases.
c<sc>
c	                     initial call to set up source dependent constants
      entry kellref(ecolat)
c                                            
      s3 = sqrt(3.0)/2.0
      sc0 = 0.25*(1.0+3.0*cos(2.0*ecolat))
      sc1 = s3*sin(2.0*ecolat)
      sc2 = s3*sin(ecolat)*sin(ecolat)
      return
c<sc>
c<ec>                                           phase identification
      entry kellcor(phase,edist,edepth,ecolat,azim,tcor,abrt)
*      write(6,*) "phase,edist,edepth,ecolat,azim"
*      write(6,*)  phase,edist,edepth,ecolat,azim
      Nd = 6
      NUMPH = 57
      deldst = 5.0
      abrt = .FALSE.
c                                             check on the length of phase
      l=len(phase)
      if(l.lt.8) then
       stop 
     >    'character variable `phase` should have at least length 8'
      endif

c                                             select phase
      ip = -1
      nc=min(lnblk(phase),8)
      do 10 i=1,NUMPH
        if(nc.ne.phnch(i)) goto 10
        if (phase(1:nc) .eq. phcod(i)(1:nc)) then
          ip = i
          go to 11
        endif
 10   continue
 11   continue

      if(ip.eq.-1) then
c                                             check phase aliases
        call phase_alias(phase,edist,ip)
      endif
      Ne = phspn(ip)
*      write(6,*) "ip:",ip
c                                              phase not found
      if(ip.lt.0) then
        write(6,*) phase,'  is not available'
        abrt = .true.
        return
      endif
c                                               special case of upgoing waves
*
c
c                                                acquire phase information
       nr = phind(ip)
*       write(6,*) "nrec:",nr
       read(21,61,rec=nr) phcod(ip),np,d1,d2
*       write(6,*) "phcode,np,d1,d2: ", phcod(ip),np,d1,d2
       nr = nr+1
       if(np.ne.Ne) write(6,*) "HELP! - index wrong"
       do 15 i=1,np
         read(21,62,rec=nr) delta(i)
         nr = nr+1
         read(21,63,rec=nr) (t0(i,m),m=1,6)
         nr = nr+1
         read(21,63,rec=nr) (t1(i,m),m=1,6)
         nr = nr+1
         read(21,63,rec=nr) (t2(i,m),m=1,6)
         nr = nr+1
 15    continue         
 61    format(a8,i10,2f10.0)
 62    format(f10.0)
 63    format(6f10.4)
c                                  distance index
       idist = 1 + int( (edist-d1)/ deldst )
       if(edist.lt.d1) idist =1
       if(edist.gt.d2) idist= np-1
c                                  depth index
       do 25 j = 1,Nd-1
         if ((dpth(j).le.edepth).and.(dpth(j+1).ge.edepth))then
            jdepth = j
            goto 26
         endif
 25    continue
 26    continue
*       write(6,*) "idist, jdepth;",idist,jdepth
c
*                      need to allow for zero entries (where phase
*                      description strongly depth dependent)
c tau0
         a0 = t0(idist,jdepth)
         b0 = t0(idist,jdepth+1)
         h0 = t0(idist+1,jdepth+1)
         d0 = t0(idist+1,jdepth)
         e0 = a0 + 
     ^       (d0-a0)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f0 = b0 + 
     ^       (h0-b0)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g0 = e0 + 
     ^       (f0-e0)*(edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau0 = g0
c tau1
         a1 = t1(idist,jdepth)
         b1 = t1(idist,jdepth+1)
         h1 = t1(idist+1,jdepth+1)
         d1 = t1(idist+1,jdepth)
         e1 = a1 + 
     ^       (d1-a1)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f1 = b1 + 
     ^       (h1-b1)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g1 = e1 + 
     ^       (f1-e1)*(edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau1 = g1
c tau2
         a2 = t2(idist,jdepth)
         b2 = t2(idist,jdepth+1)
         h2 = t2(idist+1,jdepth+1)
         d2 = t2(idist+1,jdepth)
         e2 = a2 + 
     ^       (d2-a2)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f2 = b2 + 
     ^       (h2-b2)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g2 = e2 + 
     ^       (f2-e2)*(edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau2 = g2
c
*         write(6,*) "tau0,tau1,tau2:",tau0,tau1,tau2
         caz = cos(azim)
         cbz = cos(2.0*azim)
*         write(6,*) "azim,caz,cbz",azim,caz,cbz    
c
         tcor = sc0*tau0 + sc1*cos(azim)*tau1 + sc2*cos(2.0*azim)*tau2
c
      return
c<ec>
      end
      subroutine phase_alias(phase,delta,ip)

c-    check for alternative phase names
c     input phase, delta
c     output ip (index of phcod)

      character*(*) phase
      if(phase(1:3).eq.'Pg ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sg ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPg ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPg ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSg ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSg ') then
c       phase='sS      '
        ip=40
c
      elseif(phase(1:3).eq.'Pb ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sb ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPb ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPb ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSb ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSb ') then
c       phase='sS      '
c
      elseif(phase(1:3).eq.'Pn ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sn ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPn ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPn ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSn ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSn ') then
c       phase='sS      '
        ip=40
      else if(phase(1:4).eq.'SPn ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SPb ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SPg ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SnP ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'PSn ') then
c       phase='PS      '
        ip=56
      else if(phase(1:5).eq.'PnPn ') then
c       phase='PP      '
        ip=30
      else if(phase(1:5).eq.'SnSn ') then
c       phase='SS      '
        ip=53
c                                       upgoing P, S
      else if(phase(1:2).eq.'p ') then
c       phase='Pup     '
        ip=1  
      else if(phase(1:2).eq.'s ') then
c       phase='Sup     '
        ip=32 
c                                        
      else if(delta.le.100.0.and.phase.eq.'pPdiff  ') then
c       phase='pP      '
        ip=8
      else if(delta.le.100.0.and.phase.eq.'sPdiff  ') then
c       phase='sP      '
        ip=13
      else if(delta.le.100.0.and.phase.eq.'pSdiff  ') then
c       phase='pS      '
        ip=37
      else if(delta.le.100.0.and.phase.eq.'sSdiff  ') then
c       phase='sS      '
        ip=40
      else if(delta.ge.100.0.and.phase.eq.'pPdiff  ') then
c       phase='Pdiff   '
        ip=8
      else if(delta.ge.100.0.and.phase.eq.'sPdiff  ') then
c       phase='Pdiff   '
        ip=13
      else if(delta.ge.100.0.and.phase.eq.'pSdiff  ') then
c       phase='Sdiff   '
        ip=37
      else if(delta.ge.100.0.and.phase.eq.'sSdiff  ') then
c       phase='Sdiff    '
        ip=40
c            
      else if(delta.le.165.0.and.phase.eq.'PKPdiff ') then
c       phase='PKPbc '
        ip=5
      else if(delta.le.165.0.and.phase.eq.'pPKPdiff') then
c       phase='pPKPbc '
        ip=10
      else if(delta.le.165.0.and.phase.eq.'sPKPdiff') then
c       phase='sPKPbc '
        ip=15
c                             
      else if(phase(1:4).eq."P'P'") then
c       phase="P'P'    "
        ip =31
      else if(phase(1:4).eq."S'S'") then
c       phase="P'P'    "
        ip =54
c                            diffractions (approx)
      else if(delta.gt.100.0.and.phase.eq.'pPdiff  ') then
c       phase='Pdiff   '
        ip=8
      else if(delta.gt.100.0.and.phase.eq.'sPdiff  ') then
c       phase='Pdiff   '
        ip=13
      else if(delta.gt.100.0.and.phase.eq.'pSdiff  ') then
c       phase='Sdiff   '
        ip=37
      else if(delta.gt.100.0.and.phase.eq.'sSdiff  ') then
c       phase='Sdiff    '
c
      else
        ip=-1
      endif
      return
      end
c-----
      integer function lnblk(s)
      character *(*) s
      l = len(s)
      do i=l,1,-1
        if (s(i:i).gt.' ') then
          lnblk=i
          return
        endif
      end do   
      lnblk=0
      return
      end
frechet.f90
!************************************************************************************************
! this subroutine reads from file for which variables the inverson is to be done
! and establishes a global index for each variable
!
! the order is first the velocity grids and then interface grids, with the order inside
! each velocity or interface grid varying in the order r fastest, lat medium, long slowest
! e.g. n_loc = ir + (ilat-1)*nlat + (ilong-1)*nlat*nlong for velocities
! source position is at the end

subroutine initialize_inversion(do_frechet_derivatives)
  use mod_3dfm
  implicit none

  integer     :: ngrids,nifaces,nsources,i,j,igrid,iiface,isource,group,last_index
  logical     :: do_frechet_derivatives

  integer     :: input_list(1000),type_list(100),vtype

  n_inv_parms = 0
  n_inv_active= 0
  group = 0

  open(1,file='frechet.in')

  read(1,*) i
  if (i == 0) then
     do_frechet_derivatives = .false.
     close(1)
     return
  else
     do_frechet_derivatives = .true.
  endif



  read (1,*) ngrids
  n_inv_vgrid=ngrids

  if (n_inv_vgrid > 0) then

     read(1,*) input_list(1:ngrids)
     read(1,*) type_list(1:ngrids)

     allocate(vgrids_to_be_inv(ngrids))

     do i=1,ngrids
        igrid = input_list(i)
        vtype = type_list(i)
        vgrids_to_be_inv(i)%igrid=igrid
        vgrids_to_be_inv(i)%vtype=vtype
        vgrid(igrid,vtype)%to_be_inverted = .true.
        group=group+1
        if (group == 1 ) then
           vgrid(igrid,vtype)%start_index=1
        else
           vgrid(igrid,vtype)%start_index=last_index  
        endif
        last_index = vgrid(igrid,vtype)%start_index + vgrid(igrid,vtype)%nnode
        n_inv_parms = n_inv_parms+vgrid(igrid,vtype)%nnode
        n_inv_active= n_inv_active+count(vgrid(igrid,vtype)%active)
     end do

  endif

  read (1,*) nifaces
  n_inv_iface=nifaces

  if (n_inv_iface > 0) then

     read(1,*) input_list(1:nifaces)

     allocate(ifaces_to_be_inv(nifaces))

     do i=1,nifaces
!        read(1,*) iiface
        iiface=input_list(i)
        ifaces_to_be_inv(i)=iiface
        intrface(iiface)%to_be_inverted = .true.
        group=group+1
        if (group == 1 ) then
           intrface(iiface)%start_index=1
        else
           intrface(iiface)%start_index=last_index
        endif
        last_index = intrface(iiface)%start_index + intrface(iiface)%nnode
        n_inv_parms = n_inv_parms + intrface(iiface)%nnode
        n_inv_active= n_inv_active + intrface(iiface)%nnode
     end do

  endif

  read(1,*) nsources
  if (nsources > 0) then ; locate_source=.true.; else ; locate_source=.false. ; endif

  if (locate_source) then

     n_inv_source=nsources

     read(1,*) input_list(1:nsources)

     allocate(sources_to_be_inv(nsources))
     do i=1,nsources
        isource = input_list(i)
        sources_to_be_inv(i)=isource
        source(isource)%to_be_inverted = .true.
        group=group+1
        if (group == 1) then
           source(isource)%start_index=1
        else
           source(isource)%start_index=last_index 
        endif
        last_index = source(isource)%start_index + 4
        n_inv_parms =n_inv_parms+4
        n_inv_active=n_inv_active+4
        j=4
     end do
  end if

  close(1)

  return


end subroutine initialize_inversion


!*****************************************************

! this subroutine takes a ray as input and produces the partial derivatives of the
! arrival time of this ray with respect to the velocity grids it encounters
! partial derivatives are stored in compact row storage as part of the ray derived type

subroutine ray_partials(ray)
  use mod_3dfm
  implicit none

  real(kind=dp)                                     :: r,lat,long
  type(Tray)                                        :: ray
  type(Tray_section),pointer                        :: raysec,raysec_prev
  type(Tvelocity_grid), pointer                     :: gridv
  type(Tintersection),pointer                       :: isec
  type(Tinterface),pointer                          :: iface
  integer, dimension(:),pointer                     :: regid

  real(kind=dp),dimension(:),allocatable             :: dtdpar,r_interface
  logical,dimension(:),allocatable                   :: iface_pinched

  integer                            :: start_index,glob_inv_index

  integer        :: i,j,k,ii,jj,kk,ir,ilat,ilong,n,m,nr,nlat,nlatnr,icell,ninodes,n_pinched_interfaces
  integer        :: iftop,ifbot,vtype
  real(kind=dp) :: u,v,w,bu(4),bv(4),bw(4),vel,weight,dl,pinch_modifier
  real(kind=dp) :: dpos(3),gradt_in(3),normal(3),h,gw(20),norm,dist,distmin
  real(kind=dp) :: gradt_in_perp,gradt_out_perp,gradt_in_par,geo_factor,vel_out,det,vel_in
  real(kind=dp) :: interpolate_velocity,interpolate_interface
  logical       :: rec_region_done,do_interface_partials

  ! allocate a temporary array to store the row in the inversion matrix corresponding
  ! to this ray. At the end the zeros will be removed by conversion to CRS

  allocate(dtdpar(n_inv_parms))
  dtdpar=0.0_dp

 ! some local arrays to assist with the identification of pinched interfaces
  allocate(r_interface(n_interfaces),iface_pinched(n_interfaces))

 ! raysec_prev is a pointer to the previous ray section in the sense of ray path integration

  nullify(raysec_prev)

 ! these are used to suppress interface derivative calculation in the receiver section
  rec_region_done = .false.
  do_interface_partials = .false.

  ! loop over the sections of this ray starting at receiver

  do m=ray%nsections,1,-1

!     print *,'ray_partial: section',m

     raysec => ray%section(m)

   ! if the ray section contains only two points (start and end) it is considered to lie
   ! between pinched interfaces, and does not influence the arrival time

     if (raysec%npoints <= 2) then

        cycle    ! go to the next section in the sense of path integration, the previous in the time field sequence

     else

   ! for the section starting at the receiver, there are no partial derivatives at the starting interface
   ! also set the previous section to the first significant (not pinched) section on the ray

        if (.not. rec_region_done) then
           raysec_prev => raysec
           rec_region_done = .true.
        else
           do_interface_partials = .true.
        endif

!        print *,'more than 2 points'

     ! first the partial derivatives with respect to the velocity grid

        if (raysec%tf%vtype == 0) then
           print *,raysec%tf%vtype
           print *,raysec%tf%id,raysec%tf%reg%id,raysec%tf%istart%id
           stop 'illegal vtype in ray_partials'
        endif

        vtype = raysec%tf%vtype
        gridv => vgrid(raysec%reg%ivgrid,vtype)   ! the velocity grid that applies to this section


        if (gridv%to_be_inverted) then

!           print *,'vgrid partials for section',m,'velocity grid',raysec%reg%ivgrid

           start_index = gridv%start_index     ! the starting index of the velocity grid parameters in the global list

           nr= gridv%nr
           nlat = gridv%nlat
           nlatnr=nlat*nr


           do n=2,raysec%npoints

              r    =  (raysec%point(1,n)+raysec%point(1,n-1))*0.5_dp
              lat  =  (raysec%point(2,n)+raysec%point(2,n-1))*0.5_dp
              long =  (raysec%point(3,n)+raysec%point(3,n-1))*0.5_dp


              dl=sqrt(  (raysec%point(1,n)-raysec%point(1,n-1))**2 + &
                   (r*(raysec%point(2,n)-raysec%point(2,n-1)))**2 + &
                   (r*cos(lat)*(raysec%point(3,n)-raysec%point(3,n-1)))**2 )



              ir=floor((r-gridv%r0)/gridv%dr0)+1
              ilat=floor((lat-gridv%lat0)/gridv%dlat0)+1
              ilong=floor((long-gridv%long0)/gridv%dlong0)+1

              u=(r-gridv%r(ir))/gridv%dr0
              v=(lat-gridv%lat(ilat))/gridv%dlat0
              w=(long-gridv%long(ilong))/gridv%dlong0

              bu(1)=(1.0_dp-u)**3/6.0_dp
              bu(2)=(4.0_dp-6.0_dp*u**2+3.0_dp*u**3)/6.0_dp
              bu(3)=(1.0_dp+3.0*u+3.0_dp*u**2-3.0_dp*u**3)/6.0_dp
              bu(4)=u**3/6.0_dp
              bv(1)=(1.0_dp-v)**3/6.0_dp
              bv(2)=(4.0_dp-6.0_dp*v**2+3.0_dp*v**3)/6.0_dp
              bv(3)=(1.0_dp+3.0*v+3.0_dp*v**2-3.0_dp*v**3)/6.0_dp
              bv(4)=v**3/6.0_dp
              bw(1)=(1.0_dp-w)**3/6.0_dp
              bw(2)=(4.0_dp-6.0_dp*w**2+3.0_dp*w**3)/6.0_dp
              bw(3)=(1.0_dp+3.0*w+3.0_dp*w**2-3.0_dp*w**3)/6.0_dp
              bw(4)=w**3/6.0_dp


              vel=0.0_dp

              do k=1,4
                 kk=ilong+k-2
                 do j=1,4
                    jj=ilat+j-2
                    do i=1,4
                       ii=ir+i-2

                       weight=bu(i)*bv(j)*bw(k)
                       vel = vel + weight*gridv%velocity(ii,jj,kk)

                    end do
                 end do
              end do

              do k=1,4
                 kk=ilong+k-2
                 do j=1,4
                    jj=ilat+j-2
                    do i=1,4
                       ii=ir+i-2

                       weight=bu(i)*bv(j)*bw(k)

                       glob_inv_index=start_index + nlatnr*(kk-1) + nr*(jj-1) + (ii-1)

                       dtdpar(glob_inv_index) = dtdpar(glob_inv_index) - dl*weight/vel**2

                    end do
                 end do
              end do


           end do  ! ray section points loop

!        print *,'vgrid partials for section',m,'finished'

        endif    ! vgrid is to be inverted


     ! now the partial derivatives with respect to interface position

        if (do_interface_partials) then   

!           print *,'considering interface partials for section',m

           isec => raysec%istart
           iface => intrface(isec%iface_id)

           if (iface%to_be_inverted) then

!           print *,'interface partials for iface',iface%id,raysec_prev%reg%id,raysec%reg%id


        ! the first point of the ray section is the intersection point

           r    =  raysec%point(1,1)
           lat  =  raysec%point(2,1)
           long =  raysec%point(3,1)



       ! test if interfaces are pinched at the position the ray hits 

           iface_pinched = .false.  ! elements of this array will be set to .true. if pinched at this position

           if (raysec%reg%id /= raysec_prev%reg%id) then

           ! if it is a refraction, we can use the current and previous region to identify
           ! pinched interfaces

              n_pinched_interfaces = abs(raysec_prev%reg%id - raysec%reg%id)
              iftop = min(raysec_prev%reg%id,raysec%reg%id) + 1
              ifbot = max(raysec_prev%reg%id,raysec%reg%id)
              iface_pinched(iftop:ifbot) = .true. 

           else

            ! if it is a reflection, can't use path to identify pinched interface
            ! must do it directly

              do i=1,n_interfaces
                 r_interface(i)=interpolate_interface(lat,long,intrface(i))
              end do

              n_pinched_interfaces = 0

              do i=1,n_interfaces

                 if (abs(r_interface(i)-r_interface(iface%id)) < pgrid%tolerance) then

                    n_pinched_interfaces=n_pinched_interfaces+1
                    iface_pinched(i) = .true.

                 endif

              end do             

           endif


       ! the frechet derivatives will be divided by this modifier to account for pinched interfaces

           pinch_modifier = dble(count(iface_pinched(1:n_interfaces) .and. intrface(1:n_interfaces)%to_be_inverted))

!           print *,'pinch modifier is',pinch_modifier


              start_index = iface%start_index ! starting index of the intrface parameters in the global list
              nlat=iface%nlat

           ! first interpolate the value of the incoming time gradient at the interface

           ! determine if interface is at the top or bottom of the region

              if (isec%id == raysec%reg%itop%id) then
                 regid => isec%rbel_node_id
              else
                 if (isec%id == raysec%reg%ibot%id) then
                    regid => isec%rabo_node_id
                 else
                    stop 'error a in ray partials'
                 endif
              endif


           ! get the cut cell

              ir=min(pgrid%nr-1,floor((r-pgrid%r0)/pgrid%dr0)+1)
              ilat=min(pgrid%nlat-1,floor((lat-pgrid%lat0)/pgrid%dlat0)+1)
              ilong=min(pgrid%nlong-1,floor((long-pgrid%long0)/pgrid%dlong0)+1)


              if (associated(pgrid%ccind_from_3dc(ir,ilat,ilong)%p)) then

            ! if the first point lies in a cut cell (always except problematic cases)

                 icell = pgrid%ccind_from_3dc(ir,ilat,ilong)%p(isec%iface_id)

                 if (icell == 0)  stop 'error 2 in ray partials'


           ! interpolate based on inodes of the cut cell

                 ninodes=isec%n_inodes(icell)
                 do j=1,ninodes
                    gw(j)=1.d0/(sqrt((r-isec%r(isec%inodes(j,icell)))**2 + &
                         (r*(lat-isec%lat(isec%inodes(j,icell))))**2 &
                         +(r*cos(lat)*(long-isec%long(isec%inodes(j,icell)))**2)) + 0.05_dp*pgrid%dr0)
                 end do

                 gw(1:ninodes) = gw(1:ninodes)/sum(gw(1:ninodes))

                 gradt_in=0.0_dp
                 do j=1,ninodes
                    gradt_in = gradt_in + gw(j)*raysec%tf%time_gradient(1:3,regid(isec%inodes(j,icell)))
                 end do

              else

            ! if the first point does not lie in a cut cell (interfaces crossing back and
            ! forth along a coordinate plane on a scale less than the grid spacing)
            ! take the gradient at the closest point

                 distmin=10.0*earth_radius**2
                 do j=1,raysec%reg%nnode
                    dist=(r-raysec%reg%r(j))**2 + (r*(lat-raysec%reg%lat(j)))**2 + &
                         (r*cos(lat)*(long-raysec%reg%long(j)))**2
                    if (dist < distmin) then
                       distmin=dist
                       jj=j
                    endif
                 end do

                 gradt_in=raysec%tf%time_gradient(1:3,jj)

              endif

!         ensure the interpolated time gradient has the correct norm

              vel_in = interpolate_velocity(r,lat,long,vgrid(raysec%reg%ivgrid,vtype))
              norm = sqrt(sum(gradt_in**2))
              gradt_in = gradt_in/(vel_in*norm)

!           write(24,*) 'gin',gradt_in

           ! we have the incoming time gradient, now the surface normal

              call interface_normal(lat,long,intrface(isec%iface_id),normal(1),normal(2),normal(3),h)

!           write(24,*) 'nor',normal

           ! the component of the incoming time gradient normal to the surface 

              gradt_in_perp = dot_product(gradt_in,normal)


           ! the component of the outgoing time gradient normal to the surface 
        
              if (raysec%reg%id == raysec_prev%reg%id)  then   ! reflection

                 gradt_out_perp = -gradt_in_perp

              else                                                  ! refraction

                 gradt_in_par=sqrt(sum((gradt_in-gradt_in_perp*normal)**2))
              
                 vel_out=interpolate_velocity(r,lat,long,vgrid(raysec_prev%reg%ivgrid,raysec_prev%tf%vtype))

                 det = 1.0_dp/vel_out**2 - gradt_in_par**2
                 if (det < 0.0_dp) then

                    print *, 'warning 3 in ray partials : violating total reflection'
                    print *, 'section',m
                    print *, 'region_in',raysec%tf%reg%id,'region_out',raysec_prev%tf%reg%id
                    print *, 'vel_in',vel_in,'vel_out',vel_out
                    print *, det,1.0_dp/vel_out**2,gradt_in_par**2
                    print *, gradt_in
                    print *, gradt_in_perp*normal
                    print *, 1.0/sqrt(sum(gradt_in**2)),1.0/gradt_in_perp,1.0/sqrt(gradt_in_perp**2+gradt_in_par**2)
                    print *, 1.0/gradt_in_par,vel_out
!                    stop 'error 3 in ray partials : violating total reflection'
                    det = 0.0_dp
                 endif

                 gradt_out_perp = sign(sqrt(det),gradt_in_perp)

              endif


           ! finally we can calculate the geometrical factor

              geo_factor=(gradt_in_perp-gradt_out_perp)*normal(1)


           ! now calculate the partial derivatives

              ilat=floor((lat-iface%lat0)/iface%dlat0)+1
              ilong=floor((long-iface%long0)/iface%dlong0)+1

              v=(lat-iface%lat(ilat))/iface%dlat0
              w=(long-iface%long(ilong))/iface%dlong0

              bv(1)=(1.0_dp-v)**3/6.0_dp
              bv(2)=(4.0_dp-6.0_dp*v**2+3.0_dp*v**3)/6.0_dp
              bv(3)=(1.0_dp+3.0*v+3.0_dp*v**2-3.0_dp*v**3)/6.0_dp
              bv(4)=v**3/6.0_dp
              bw(1)=(1.0_dp-w)**3/6.0_dp
              bw(2)=(4.0_dp-6.0_dp*w**2+3.0_dp*w**3)/6.0_dp
              bw(3)=(1.0_dp+3.0*w+3.0_dp*w**2-3.0_dp*w**3)/6.0_dp
              bw(4)=w**3/6.0_dp


              do k=1,4
                 kk=ilong+k-2
                 do j=1,4
                    jj=ilat+j-2

                    weight=bv(j)*bw(k)

                    do n=1,n_interfaces

                       if (iface_pinched(n) .and. intrface(n)%to_be_inverted) then

                          glob_inv_index=intrface(n)%start_index + nlat*(kk-1) + (jj-1)

                          dtdpar(glob_inv_index) = dtdpar(glob_inv_index) + weight*geo_factor/pinch_modifier

                       endif

                    end do

                 end do
              end do
   
           endif   ! iface to be inverted

        endif  ! section is not the first significant (first section considered starts at receiver, not an interface)

        raysec_prev => raysec  ! points to the last significant (non-pinched) section along the ray

     endif   ! interfaces not pinched at position ray passes through

  end do ! ray sections loop

!  print *,'finished ray section loop'

  ! the partial derivatives with respect to source position and time (source is last point on the ray)

  if (locate_source) then

     if (ray%source%to_be_inverted) then

        n=ray%section(1)%npoints

        dpos(1)= ray%section(1)%point(1,n)-raysec%point(1,n-1)
        dpos(2)= ray%section(1)%point(1,n)*(ray%section(1)%point(2,n)-raysec%point(2,n-1))
        dpos(3)= ray%section(1)%point(1,n)*cos(ray%section(1)%point(2,n))*(ray%section(1)%point(3,n)-raysec%point(3,n-1))

        dpos=dpos/sqrt(sum(dpos**2))

        vel=interpolate_velocity(ray%section(1)%point(1,n),ray%section(1)%point(2,n), &
             ray%section(1)%point(3,n),vgrid(ray%section(1)%reg%ivgrid,ray%section(1)%tf%vtype))


        glob_inv_index=ray%source%start_index

        dtdpar(glob_inv_index:glob_inv_index+2)=dpos/vel

        dtdpar(glob_inv_index+3) = -1.0_dp

     endif

  endif


! store the partial derivatives for this ray (row in inversion matrix)) in CRS

  ray%n_pdev = count(dtdpar /= 0.0_dp)
  allocate(ray%pdev(ray%n_pdev),ray%pdev_indx(ray%n_pdev))
  m=0
  do n=1,n_inv_parms
     if (dtdpar(n) /= 0.0_dp) then
        m=m+1
        ray%pdev(m) = dtdpar(n)
        ray%pdev_indx(m)=n
     endif
  end do

  deallocate(dtdpar)

  return

end subroutine ray_partials

!*************************************************************************************************

subroutine decode_global_index(glob_index,i1,i2,i3,p_type,p_index,p_subindex)
  use mod_3dfm
  implicit none

  integer                       ::  glob_index,i1,i2,i3,p_type,p_index,p_subindex,n_groups
  integer                       ::  start(31)
  integer                       ::  i,n,m,group,nr,nlat,loc_index

  p_subindex = 0

  n_groups=n_inv_vgrid+n_inv_iface

  if (locate_source) n_groups=n_groups+n_inv_source

  if (n_groups > 30) stop 'too many groups in decode_global_index'

  m=0
  do i=1,n_inv_vgrid
     m=m+1
     start(m)=vgrid(vgrids_to_be_inv(i)%igrid,vgrids_to_be_inv(i)%vtype)%start_index
  end do

  do i=1,n_inv_iface
     m=m+1
     start(m)=intrface(ifaces_to_be_inv(i))%start_index
  end do

  if (locate_source) then
     do i=1,n_inv_source
        m=m+1
        start(m)=source(sources_to_be_inv(i))%start_index
     end do
  endif

  start(n_groups+1)=n_inv_parms+1

  do n=1,n_groups
     if (glob_index < start(n+1)) then
        group=n
        exit
     endif 
  end do

  if (group <= n_inv_vgrid) then

!     write(24,*) 'group',group,start(1:n_groups)

     p_type=1
     p_index=vgrids_to_be_inv(group)%igrid
     p_subindex=vgrids_to_be_inv(group)%vtype
     
     nr=vgrid(p_index,p_subindex)%nr
     nlat=vgrid(p_index,p_subindex)%nlat

     loc_index=glob_index-start(group)
     i3= loc_index/(nlat*nr) +1                      
     i2=(loc_index-(i3-1)*nlat*nr)/nr + 1            
     i1= loc_index-(i3-1)*nlat*nr-(i2-1)*nr + 1    


     return

  endif

  if (group <= (n_inv_vgrid+n_inv_iface) ) then

!     write(24,*) 'group',group,start(1:n_groups)

     p_type=2
     p_index=ifaces_to_be_inv(group-n_inv_vgrid)
     
     nlat=intrface(p_index)%nlat

     loc_index=glob_index-start(group)
     i3= loc_index/nlat + 1                          
     i2= loc_index-(i3-1)*nlat + 1                   
     i1=0

     return

  endif

  if (group <= n_groups ) then

     p_type=3
     p_index=sources_to_be_inv(group-n_inv_vgrid-n_inv_iface)

     loc_index=glob_index-start(group)
     i1=loc_index+1 ; i2=0 ; i3=0

     return

  endif

  stop 'something wrong in decode_global_index'

end subroutine decode_global_index

!*********************************************************************************************
! this subroutine sets a flag at the nodes of the velocity grids belonging to a region
! to true if the velocity grid node actually influences the velocity field in the region

subroutine tag_active_vgrid_nodes(reg)

  use mod_3dfm
  implicit none

  type(Tregion)  :: reg
  type(Tvelocity_grid),pointer :: gridv

  integer :: ir,ilat,ilong,n,vtype


  do vtype=1,n_vtypes

     gridv => vgrid(reg%ivgrid,vtype) 

     do n=1,reg%nnode

        ir=floor((reg%r(n)-gridv%r0)/gridv%dr0)+1
        ilat=floor((reg%lat(n)-gridv%lat0)/gridv%dlat0)+1
        ilong=floor((reg%long(n)-gridv%long0)/gridv%dlong0)+1

        gridv%active(ir-1:ir+2,ilat-1:ilat+2,ilong-1:ilong+2) = .true.

     end do  ! regional node loop

  end do  ! vtypes

end subroutine tag_active_vgrid_nodes
ftest.f90
program ftest


  a=10.0
  n=7899
  write(n,'(f12.3)') a

end program ftest



idefs_ak135.f90
subroutine get_n_interfaces(nif)
  integer :: nif
  nif=6
  return
end subroutine get_n_interfaces



function iface1(lat,long)
  real(kind=8)  :: iface1,lat,longev
  iface1=0.0d0
  return
end function iface1

function iface2(lat,long)
  real(kind=8)  :: iface2,lat,long
  iface2=20.0d0
  return
end function iface2

function iface3(lat,long)
  real(kind=8)  :: iface3,lat,long
  iface3=35.0d0
  return
end function iface3

function iface4(lat,long)
  real(kind=8)  :: iface4,lat,long
  iface4=410.0d0
  return
end function iface4

function iface5(lat,long)
  real(kind=8)  :: iface5,lat,long
  iface5=660.0d0
  return
end function iface5


! bottom of the propagation grid

function iface6(lat,long)
  real(kind=8)  :: iface6,lat,long
  iface6=1.0d3
  return
end function iface6


!function iface6(lat,long)
!  real(kind=8)  :: iface6,lat,long
!  iface6=2.8915d3
!  return
!end function iface6

! not used below this line

function iface7(lat,long)
  real(kind=8)  :: iface7,lat,long
  iface7=1.0d3
  return
end function iface7

function iface8(lat,long)
  real(kind=8)  :: iface8,lat,long
  iface8=1.0d3
  return
end function iface8

function iface9(lat,long)
  real(kind=8)  :: iface9,lat,long
  iface9=1.0d3
  return
end function iface9

function iface10(lat,long)
  real(kind=8)  :: iface10,lat,long
  iface10=1.0d3
  return
end function iface10

function iface11(lat,long)
  real(kind=8)  :: iface11,lat,long
  iface11=1.0d3
  return
end function iface11
idefs_cmplx.f90
subroutine get_n_interfaces(nif)
  integer :: nif
  nif=6
  return
end subroutine get_n_interfaces


function iface1(lat,long)
  real(kind=8)  :: iface1,lat,long
  iface1=0.0d0
  return
end function iface1

function iface2(lat,long)
  real(kind=8)  :: iface2,lat,long,x,y,z,a,b,h1,h2
  real(kind=8)  :: deg_to_rad,twopi,wlat,wlong,xx,yy
  deg_to_rad=acos(-1.0d0)/180.0d0
  twopi=2.0d0*acos(-1.0d0)
  wlat=4.0d0*deg_to_rad
  wlong=4.0d0*deg_to_rad
  x=lat*deg_to_rad/1.0d0
  y=long*deg_to_rad/1.0d0
  xx=(lat+0.02542)/1.0504
  yy=(long+0.02542)/1.0504
  z=50.0d0
  a=z/15.
  b=z/10.
  h1=z*0.05d0
  h2=z*0.15d0-a*(1.0-cos(twopi*x/(wlat))*cos(twopi*y/(wlong)))-(b*(xx+yy))
  iface2=min(h1,h2)
  return
end function iface2

function iface3(lat,long)
  real(kind=8)  :: iface3,lat,long,x,y,z,a,b,h1,h2
  real(kind=8)  :: deg_to_rad,twopi,wlat,wlong,xx,yy
  deg_to_rad=acos(-1.0d0)/180.0d0
  twopi=2.0d0*acos(-1.0d0)
  wlat=4.0d0*deg_to_rad
  wlong=4.0d0*deg_to_rad
  x=lat*deg_to_rad/1.0d0
  y=long*deg_to_rad/1.0d0
  xx=(lat+0.02542)/1.0504
  yy=(long+0.02542)/1.0504
  z=50.0d0
  a=z/10.
  b=z/10.
  h1=z*0.20d0-a*(1.0 - cos(twopi*x/(wlat))*cos(twopi*y/(wlong)))-(b*(xx+yy))
  h2=z*0.1d0
  iface3=min(h1,h2)
  return
end function iface3

function iface4(lat,long)
  real(kind=8)  :: iface4,lat,long,x,y,z,a,b
  real(kind=8)  :: deg_to_rad,twopi,wlat,wlong,xx,yy
  deg_to_rad=acos(-1.0d0)/180.0d0
  twopi=2.0d0*acos(-1.0d0)
  wlat=4.0d0*deg_to_rad
  wlong=4.0d0*deg_to_rad
  x=lat*deg_to_rad/1.0d0
  y=long*deg_to_rad/1.0d0
  xx=(lat+0.02542)/1.0504
  yy=(long+0.02542)/1.0504
  z=50.0d0
  a=z/3.
  b=z/5.
  iface4=z*0.08d0+a*(1.0-cos(twopi*x/wlat)*cos(twopi*y/wlong))+b*(xx**2+yy**2)
  return
end function iface4

function iface5(lat,long)
  real(kind=8)  :: iface5,lat,long,x,y,z,a,b
  real(kind=8)  :: deg_to_rad,twopi,wlat,wlong,xx,yy
  deg_to_rad=acos(-1.0d0)/180.0d0
  twopi=2.0d0*acos(-1.0d0)
  wlat=4.0d0*deg_to_rad
  wlong=4.0d0*deg_to_rad
  x=lat*deg_to_rad/1.0d0
  y=long*deg_to_rad/1.0d0
  xx=(lat+0.02542)/1.0504
  yy=(long+0.02542)/1.0504
  z=50.0d0
  a=z/3.
  b=z/5.
  iface5=z*0.28d0+a*(1.0-cos(twopi*x/wlat)*cos(twopi*y/wlong))+b*(xx**2+yy**2)
  return
end function iface5

function iface6(lat,long)
  real(kind=8)  :: iface6,lat,long
  iface6=5.0d1
  return
end function iface6


! not used below this line

function iface7(lat,long)
  real(kind=8)  :: iface7,lat,long
  iface7=1.0d3
  return
end function iface7

function iface8(lat,long)
  real(kind=8)  :: iface8,lat,long
  iface8=1.0d3
  return
end function iface8

function iface9(lat,long)
  real(kind=8)  :: iface9,lat,long
  iface9=1.0d3
  return
end function iface9

function iface10(lat,long)
  real(kind=8)  :: iface10,lat,long
  iface10=1.0d3
  return
end function iface10

function iface11(lat,long)
  real(kind=8)  :: iface11,lat,long
  iface11=1.0d3
  return
end function iface11
idefs_simple.f90
subroutine get_n_interfaces(nif)
  integer :: nif
  nif=4
  return
end subroutine get_n_interfaces


function iface1(lat,long)
  real(kind=8)  :: iface1,lat,long
  iface1=0.0d0
  return
end function iface1

function iface2(lat,long)
  real(kind=8)  :: iface2,lat,long,x,y,z,b
  x=lat/20.0
  y=long/20.0
  z=1000.0
  b=z/1.0
  iface2=z*0.65 - b*x 
  return
end function iface2

function iface3(lat,long)
  real(kind=8)  :: iface3,lat,long,x,y,z,b
  x=lat/20.0
  y=long/20.0
  z=1000.0
  b=z/5.
  iface3=z*0.45 + b*x 
  return
end function iface3

function iface4(lat,long)
  use mod_3dfm
  real(kind=8)  :: iface4,lat,long
  iface4=1.0d3
  return
end function iface4




! not used below this line

function iface5(lat,long)
  real(kind=8)  :: iface5,lat,long
  iface5=1.0d3
  return
end function iface5


function iface6(lat,long)
  real(kind=8)  :: iface6,lat,long
  iface6=1.0d3
  return
end function iface6


function iface7(lat,long)
  real(kind=8)  :: iface7,lat,long
  iface7=1.0d3
  return
end function iface7

function iface8(lat,long)
  real(kind=8)  :: iface8,lat,long
  iface8=1.0d3
  return
end function iface8

function iface9(lat,long)
  real(kind=8)  :: iface9,lat,long
  iface9=1.0d3
  return
end function iface9

function iface10(lat,long)
  real(kind=8)  :: iface10,lat,long
  iface10=1.0d3
  return
end function iface10

function iface11(lat,long)
  real(kind=8)  :: iface11,lat,long
  iface11=1.0d3
  return
end function iface11
libsun.f
      subroutine warn(msg)
      character*(*) msg
      write(*,100) msg
 100  format(1x,a)
      return
      end
      subroutine tnoua(ia,nc)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine tnoua writes the first nc charcters of the
c   character string ia to the standard 
c   output without the trailing newline (allowing user input on the 
c   same line).  Programmed on 17 September 1980 by R. Buland.
c
      save
      character*(*) ia
      write(*,100)ia(1:nc)
 100  format(a,$)
      return
      end
      subroutine dasign(lu,mode,ia,len)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine dasign opens (connects) logical unit lu to the disk file
c   named by the character string ia with mode mode.  If iabs(mode) = 1,
c   then open the file for reading.  If iabs(mode) = 2, then open the
c   file for writing.  If iabs(mode) = 3, then open a scratch file for
c   writing.  If mode > 0, then the file is formatted.  If mode < 0,
c   then the file is unformatted.  All files opened by dasign are
c   assumed to be direct access.  Programmed on 3 December 1979 by
c   R. Buland.
c
      save
      character*(*) ia
      logical exst
c
      if(mode.ge.0) nf=1
      if(mode.lt.0) nf=2
      ns=iabs(mode)
      if(ns.le.0.or.ns.gt.3) ns=3
      go to (1,2),nf
 1    go to (11,12,13),ns
 11   open(lu,file=ia,status='old',form='formatted',
     1 access='direct',recl=len)
      return
 12   inquire(file=ia,exist=exst)
      if(exst) go to 11
 13   open(lu,file=ia,status='new',form='formatted',
     1 access='direct',recl=len)
      return
 2    go to (21,22,23),ns
 21   open(lu,file=ia,status='old',form='unformatted',access='direct',
     1 recl=len)
      return
 22   inquire(file=ia,exist=exst)
      if(exst) go to 21
 23   open(lu,file=ia,status='new',form='unformatted',access='direct',
     1 recl=len)
      return
      end
      subroutine vexit(ierr)
c      call exit(ierr)
      if (ierr .NE.0) print *,'exit called with error code',ierr
      stop
      end
libtau.f
      subroutine tabin(in,modnam)
      include 'ttlim.inc'
      character*(*) modnam
c     logical log
      character*8 phcd,phdif(6)
      character cdum*20
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
c
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
c      data tauc,xc/jtsm*0d0,jxsm*0d0/
c

      nin=in
      phdif(1)='P'
      phdif(2)='S'
      phdif(3)='pP'
      phdif(4)='sP'
      phdif(5)='pS'
      phdif(6)='sS'
c++
      call asnag1(nin,-1,1,'Enter model name:',modnam)

c++ 
      read(nin) nasgr,nl,len2,xn,pn,tn,mt,nseg,nbrn,ku,km,fcs,nafl,
     1 indx,kndx
      read(nin) pm,zm,ndex
      read(nin) pu,pux
      read(nin) phcd,px,xt,jndx
      read(nin) pt,taut
      read(nin) coef
      call retrns(nin)

c
      nb=index(modnam,' ')-1
      if(nb.le.0) nb=len(modnam)
      cdum=modnam(1:nb)//'.tbl'

      call dasign(nin,-1,cdum,nasgr)
c
      do 11 nph=1,2
 11   pu(ku(nph)+1,nph)=pm(1,nph)
c
c     write(10,*)'nasgr nl len2',nasgr,nl,len2
c     write(10,*)'nseg nbrn mt ku km',nseg,nbrn,mt,ku,km
c     write(10,*)'xn pn tn',xn,pn,tn
c     write(10,200)(i,(ndex(i,j),pm(i,j),zm(i,j),j=1,2),i=1,mt(2))
c200  format(/(1x,i3,i7,2f12.6,i7,2f12.6))
c     write(10,201)(i,(pu(i,j),j=1,2),i=1,ku(2)+1)
c201  format(/(1x,i3,2f12.6))
c     write(10,201)(i,(pux(i,j),j=1,2),i=1,km(2))
c     write(10,202)(i,(nafl(i,j),j=1,3),(indx(i,j),j=1,2),(kndx(i,j),
c    1 j=1,2),(fcs(i,j),j=1,3),i=1,nseg)
c202  format(/(1x,i3,7i5,3f5.0))
c     cn=180./3.1415927
c     write(10,203)(i,(jndx(i,j),j=1,2),(px(i,j),j=1,2),(cn*xt(i,j),
c    1 j=1,2),phcd(i),i=1,nbrn)
c203  format(/(1x,i3,2i5,2f12.6,2f12.2,2x,a))
c     write(10,204)(i,pt(i),taut(i),(coef(j,i),j=1,5),i=1,jout)
c204  format(/(1x,i5,0p2f12.6,1p5d10.2))
c
      tn=1./tn
      dn=3.1415927/(180.*pn*xn)
      odep=-1.
      ki=0
      msrc(1)=0
      msrc(2)=0
      k=1
      do 3 i=1,nbrn
      jidx(i)=jndx(i,2)
      do 4 j=1,2
 4    dbrn(i,j)=-1d0
 8    if(jndx(i,2).le.indx(k,2)) go to 7
      k=k+1
      go to 8
 7    if(nafl(k,2).gt.0) go to 9
      ind=nafl(k,1)
      l=0
      do 10 j=jndx(i,1),jndx(i,2)
      l=l+1
 10   tp(l,ind)=pt(j)
 9    if(nafl(k,1).gt.0.and.(phcd(i)(1:1).eq.'P'.or.
     1 phcd(i)(1:1).eq.'S')) go to 3
      do 5 j=1,6
      if(phcd(i).eq.phdif(j)) go to 6
 5    continue
      go to 3
 6    dbrn(i,1)=1d0
      phdif(j)=' '
 3    continue
c     write(10,205)(i,phcd(i),(dbrn(i,j),j=1,2),jidx(i),i=1,nbrn)
c205  format(/(1x,i5,2x,a,2f8.2,i5))
c     write(10,206)(i,(tp(i,j),j=1,2),i=1,jbrnu)
c206  format(/(1x,i5,2f12.6))
      return
      end
      subroutine depset(dep,usrc)
      save 
      include 'ttlim.inc'
      logical dop,dos,segmsk,prnt
      character*8 phcd
      real usrc(2)
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      common/prtflc/segmsk(jseg),prnt(2)
c      data segmsk,prnt/jseg*.true.,2*.false./
c

      if(amax1(dep,.011).ne.odep) go to 1
      dop=.false.
      dos=.false.
      do 2 i=1,nseg
      if(.not.segmsk(i).or.iidx(i).gt.0) go to 2
      if(iabs(nafl(i,1)).le.1) dop=.true.
      if(iabs(nafl(i,1)).ge.2) dos=.true.
 2    continue
      if(.not.dop.and..not.dos) return
      go to 3
c
 1    nph0=0
      int0(1)=0
      int0(2)=0
      mbr1=nbrn+1
      mbr2=0
      dop=.false.
      dos=.false.
      do 4 i=1,nseg
      if(.not.segmsk(i)) go to 4
      if(iabs(nafl(i,1)).le.1) dop=.true.
      if(iabs(nafl(i,1)).ge.2) dos=.true.
 4    continue
      do 5 i=1,nseg
      if(nafl(i,2).gt.0.or.odep.lt.0.) go to 5
      ind=nafl(i,1)
      k=0
      do 15 j=indx(i,1),indx(i,2)
      k=k+1
 15   pt(j)=tp(k,ind)
 5    iidx(i)=-1
      do 6 i=1,nbrn
 6    jndx(i,2)=-1
      if(ki.le.0) go to 7
      do 8 i=1,ki
      j=kk(i)
 8    pt(j)=pk(i)
      ki=0
c   Sample the model at the source depth.
 7    odep=amax1(dep,.011)
      rdep=dep
      if(rdep.lt..011) rdep=0.
      zs=amin1(alog(amax1(1.-rdep*xn,1e-30)),0.)
      hn=1./(pn*(1.-rdep*xn))
      if(prnt(1).or.prnt(2)) write(10,100)dep
 100  format(/1x,'Depth =',f7.2/)
c

 3    if(nph0.gt.1) go to 12
      if(dop) call depcor(1)
      if(dos) call depcor(2)
      go to 14
 12   if(dos) call depcor(2)
      if(dop) call depcor(1)
c
c   Interpolate all tau branches.
c
 14   j=1
      do 9 i=1,nseg
      if(.not.segmsk(i)) go to 9
      nph=iabs(nafl(i,1))
c     print *,'i iidx nph msrc nafl =',i,iidx(i),nph,msrc(nph),nafl(i,1)
      if(iidx(i).gt.0.or.(msrc(nph).le.0.and.nafl(i,1).gt.0)) go to 9
      iidx(i)=1
      if(nafl(i,2).le.0) int=nafl(i,1)
      if(nafl(i,2).gt.0.and.nafl(i,2).eq.iabs(nafl(i,1)))
     1  int=nafl(i,2)+2
      if(nafl(i,2).gt.0.and.nafl(i,2).ne.iabs(nafl(i,1)))
     1  int=iabs(nafl(i,1))+4
      if(nafl(i,2).gt.0.and.nafl(i,2).ne.nafl(i,3)) int=nafl(i,2)+6
 11   if(jndx(j,1).ge.indx(i,1)) go to 10
      j=j+1
      go to 11
 10   idel(j,3)=nafl(i,1)
c      print *,'spfit:  j int =',j,int 
      call spfit(j,int)
      mbr1=min0(mbr1,j)
      mbr2=max0(mbr2,j)
      if(j.ge.nbrn) go to 9
      j=j+1
c     print *,'j jidx indx jndx =',j,jidx(j),indx(i,2),jndx(j,2)
      if(jidx(j).le.indx(i,2).and.jndx(j,2).gt.0) go to 10
 9    continue
c     write(10,*)'mbr1 mbr2',mbr1,mbr2
c     write(10,*)'msrc isrc odep zs us',msrc,isrc,odep,sngl(zs),
c    1 sngl(us(1)),sngl(us(2))
c     write(10,200)ki,(i,iidx(i),kk(i),pk(i),i=1,nseg)
c200  format(/10x,i5/(1x,3i5,f12.6))
      usrc(1)=us(1)/pn
      usrc(2)=us(2)/pn

      return
      end
c
c----------------------------------------------------------------
c
      subroutine depcor(nph)
      save
      include 'ttlim.inc'
      character*8 phcd
      logical noend,noext,segmsk,prnt
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp,ua,taua
      double precision tup(jrec),umod,zmod,tauus1(2),tauus2(2),xus1(2),
     1 xus2(2),ttau,tx,sgn,umin,dtol,u0,u1,z0,z1,fac,du
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      common/pdec/ua(5,2),taua(5,2),deplim,ka
      common/prtflc/segmsk(jseg),prnt(2)
      equivalence (tauc,tup)
c      data tol,dtol,deplim,ka,lpower/.01,1d-6,1.1,4,7/ 
      data tol,dtol,lpower/.01,1d-6,7/
c
c     write(10,*)'depcor:  nph nph0',nph,nph0
c      print *,'entering depcor'
      if(nph.eq.nph0) go to 1
      nph0=nph
      us(nph)=umod(zs,isrc,nph)
c   If we are in a high slowness zone, find the slowness of the lid.
      umin=us(nph)
      ks=isrc(nph)
c     write(10,*)'ks us',ks,sngl(umin)
      do 2 i=1,ks
      if(pm(i,nph).gt.umin) go to 2
      umin=pm(i,nph)
 2    continue
c   Find where the source slowness falls in the ray parameter array.
      n1=ku(nph)+1
      do 3 i=2,n1
      if(pu(i,nph).gt.umin) go to 4
 3    continue
      k2=n1
      if(pu(n1,nph).eq.umin) go to 50
      write(6,*)'Source slowness too large.'
      call abort
 4    k2=i
c50   write(10,*)'k2 umin',k2,sngl(umin)
c
c   Read in the appropriate depth correction values.
c
 50   noext=.false.
      sgn=1d0
      if(msrc(nph).eq.0) msrc(nph)=1
c   See if the source depth coincides with a model sample
      ztol=xn*tol/(1.-xn*odep)
      if(dabs(zs-zm(ks+1,nph)).gt.ztol) go to 5
      ks=ks+1
      go to 6
 5    if(dabs(zs-zm(ks,nph)).gt.ztol) go to 7
c   If so flag the fact and make sure that the right integrals are
c   available.
 6    noext=.true.
      if(msrc(nph).eq.ks) go to 8
      call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
      go to 11
c   If it is necessary to interpolate, see if appropriate integrals
c   have already been read in.
 7    if(msrc(nph).ne.ks+1) go to 9
      ks=ks+1
      sgn=-1d0
      go to 8
 9    if(msrc(nph).eq.ks) go to 8
c   If not, read in integrals for the model depth nearest the source
c   depth.
      if(dabs(zm(ks,nph)-zs).le.dabs(zm(ks+1,nph)-zs)) go to 10
      ks=ks+1
      sgn=-1d0
 10   call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
c   Move the depth correction values to a less temporary area.
 11   do 31 i=1,ku(nph)
 31   tauu(i,nph)=tup(i)
      k=ku(nph)
      do 12 i=1,km(nph)
      k=k+1
      xc(i)=tup(k)
 12   xu(i,nph)=tup(k)
c     write(10,*)'bkin',ks,sngl(sgn),sngl(tauu(1,nph)),sngl(xu(1,nph))
c
c   Fiddle pointers.
c
 8    msrc(nph)=ks
c     write(10,*)'msrc sgn',msrc(nph),sngl(sgn)
      noend=.false.
      if(dabs(umin-pu(k2-1,nph)).le.dtol*umin) k2=k2-1
      if(dabs(umin-pu(k2,nph)).le.dtol*umin) noend=.true.
      if(msrc(nph).le.1.and.noext) msrc(nph)=0
      k1=k2-1
      if(noend) k1=k2
c     write(10,*)'noend noext k2 k1',noend,noext,k2,k1
      if(noext) go to 14
c
c   Correct the integrals for the depth interval [zm(msrc),zs].
c
      ms=msrc(nph)
      if(sgn)15,16,16
 16   u0=pm(ms,nph)
      z0=zm(ms,nph)
      u1=us(nph)
      z1=zs
      go to 17
 15   u0=us(nph)
      z0=zs
      u1=pm(ms,nph)
      z1=zm(ms,nph)
 17   mu=1
c     write(10,*)'u0 z0',sngl(u0),sngl(z0)
c     write(10,*)'u1 z1',sngl(u1),sngl(z1)
      do 18 k=1,k1
      call tauint(pu(k,nph),u0,u1,z0,z1,ttau,tx)
      tauc(k)=tauu(k,nph)+sgn*ttau
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 18
      xc(mu)=xu(mu,nph)+sgn*tx
c     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 18   continue
      go to 39
c   If there is no correction, copy the depth corrections to working
c   storage.
 14   mu=1
      do 40 k=1,k1
      tauc(k)=tauu(k,nph)
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 40
      xc(mu)=xu(mu,nph)
c     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 40   continue
c
c   Calculate integrals for the ray bottoming at the source depth.
c
 39   xus1(nph)=0d0
      xus2(nph)=0d0
      mu=mu-1
      if(dabs(umin-us(nph)).gt.dtol.and.dabs(umin-pux(mu,nph)).le.dtol)
     1  mu=mu-1
c   This loop may be skipped only for surface focus as range is not
c   available for all ray parameters.
      if(msrc(nph).le.0) go to 1
      is=isrc(nph)
      tauus2(nph)=0d0
      if(dabs(pux(mu,nph)-umin).gt.dtol.or.dabs(us(nph)-umin).gt.dtol)
     1  go to 48
c   If we happen to be right at a discontinuity, range is available.
      tauus1(nph)=tauc(k1)
      xus1(nph)=xc(mu)
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph)),'  *'
      go to 33
c   Integrate from the surface to the source.
 48   tauus1(nph)=0d0
      j=1
      if(is.lt.2) go to 42
      do 19 i=2,is
      call tauint(umin,pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph),ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
 19   j=i
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph))
 42   if(dabs(zm(is,nph)-zs).le.dtol) go to 33
c   Unless the source is right on a sample slowness, one more partial
c   integral is needed.
      call tauint(umin,pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph))
 33   if(pm(is+1,nph).lt.umin) go to 41
c   If we are in a high slowness zone, we will also need to integrate
c   down to the turning point of the shallowest down-going ray.
      u1=us(nph)
      z1=zs
      do 35 i=is+1,mt(nph)
      u0=u1
      z0=z1
      u1=pm(i,nph)
      z1=zm(i,nph)
      if(u1.lt.umin) go to 36
      call tauint(umin,u0,u1,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
 35   xus2(nph)=xus2(nph)+tx
c36   write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)),
c    1 sngl(xus2(nph))
 36   z1=zmod(umin,i-1,nph)
      if(dabs(z0-z1).le.dtol) go to 41
c   Unless the turning point is right on a sample slowness, one more
c   partial integral is needed.
      call tauint(umin,u0,umin,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
      xus2(nph)=xus2(nph)+tx
c     write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)),
c    1 sngl(xus2(nph))
c
c   Take care of converted phases.
c
 41   iph=mod(nph,2)+1
      xus1(iph)=0d0
      xus2(iph)=0d0
      tauus1(iph)=0d0
      tauus2(iph)=0d0
      go to (59,61),nph
 61   if(umin.gt.pu(ku(1)+1,1)) go to 53
c
c   If we are doing an S-wave depth correction, we may need range and
c   tau for the P-wave which turns at the S-wave source slowness.  This
c   would bd needed for sPg and SPg when the source is in the deep mantle.
c
      do 44 j=1,nbrn
      if((phcd(j)(1:2).ne.'sP'.and.phcd(j)(1:2).ne.'SP').or.
     1 px(j,2).le.0d0) go to 44
c     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1),
c    1 px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 44   continue
      go to 53
c
c   If we are doing an P-wave depth correction, we may need range and
c   tau for the S-wave which turns at the P-wave source slowness.  This
c   would be needed for pS and PS.
c
 59   do 60 j=1,nbrn
      if((phcd(j)(1:2).ne.'pS'.and.phcd(j)(1:2).ne.'PS').or.
     1 px(j,2).le.0d0) go to 60
c     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1),
c    1 px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 60   continue
      go to 53
c
c   Do the integral.
 45   j=1
c     write(10,*)'Depcor:  do pS or sP integral - iph =',iph
      do 46 i=2,mt(iph)
      if(umin.ge.pm(i,iph)) go to 47
      call tauint(umin,pm(j,iph),pm(i,iph),zm(j,iph),zm(i,iph),ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
 46   j=i
 47   z1=zmod(umin,j,iph)
      if(dabs(zm(j,iph)-z1).le.dtol) go to 53
c   Unless the turning point is right on a sample slowness, one more
c   partial integral is needed.
      call tauint(umin,pm(j,iph),umin,zm(j,iph),z1,ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
c     write(10,*)'is ks tauusp xusp',j,ks,sngl(tauus1(iph)),
c    1 sngl(xus1(iph))
c
 53   ua(1,nph)=-1d0
c     if(odep.ge.deplim.or.odep.le..1) go to 43
      if(odep.ge.deplim) go to 43
      do 57 i=1,nseg
      if(.not.segmsk(i)) go to 57
      if(nafl(i,1).eq.nph.and.nafl(i,2).eq.0.and.iidx(i).le.0) go to 58
 57   continue
      go to 43
c
c   If the source is very shallow, we will need to insert some extra
c   ray parameter samples into the up-going branches.
c
 58   du=amin1(1e-5+(odep-.4)*2e-5,1e-5)
c     write(10,*)'Add:  nph is ka odep du us =',nph,is,ka,odep,
c    1 sngl(du),sngl(us(nph))
      lp=lpower
      k=0
      do 56 l=ka,1,-1
      k=k+1
      ua(k,nph)=us(nph)-(l**lp)*du
      lp=lp-1
      taua(k,nph)=0d0
      j=1
      if(is.lt.2) go to 54
      do 55 i=2,is
      call tauint(ua(k,nph),pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph),
     1 ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
 55   j=i
c     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 54   if(dabs(zm(is,nph)-zs).le.dtol) go to 56
c   Unless the source is right on a sample slowness, one more partial
c   integral is needed.
      call tauint(ua(k,nph),pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
c     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 56   continue
      go to 43
c
c   Construct tau for all branches.
c

 1    mu=mu+1
 43   j=1
c     write(10,*)'mu',mu
c     write(10,*)'killer loop:'
      do 20 i=1,nseg
      if(.not.segmsk(i)) go to 20
c     write(10,*)'i iidx nafl nph',i,iidx(i),nafl(i,1),nph
      if(iidx(i).gt.0.or.iabs(nafl(i,1)).ne.nph.or.(msrc(nph).le.0.and.
     1 nafl(i,1).gt.0)) go to 20
c
      iph=nafl(i,2)
      kph=nafl(i,3)
c   Handle up-going P and S.
      if(iph.le.0) iph=nph
      if(kph.le.0) kph=nph
      sgn=isign(1,nafl(i,1))
      i1=indx(i,1)
      i2=indx(i,2)
c     write(10,*)'i1 i2 sgn iph',i1,i2,sngl(sgn),iph
      m=1
      do 21 k=i1,i2
      if(pt(k).gt.umin) go to 22
 23   if(dabs(pt(k)-pu(m,nph)).le.dtol) go to 2115
      m=m+1
      go to 23
 2115 tau(1,k)=taut(k)+sgn*tauc(m)
 21   continue 

      k=i2
c     write(10,*)'k m',k,m
      go to 24
c22   write(10,*)'k m',k,m
 22   if(dabs(pt(k-1)-umin).le.dtol) k=k-1
      ki=ki+1
      kk(ki)=k
      pk(ki)=pt(k)
      pt(k)=umin
      fac=fcs(i,1)
c     write(10,*)'ki fac',ki,sngl(fac)
      tau(1,k)=fac*(tauus1(iph)+tauus2(iph)+tauus1(kph)+tauus2(kph))+
     1 sgn*tauus1(nph)
c     write(10,*)'&&&&& nph iph kph tauus1 tauus2 tau =',
c    1 nph,iph,kph,sngl(tauus1(1)),sngl(tauus1(2)),sngl(tauus2(1)),
c    2 sngl(tauus2(2)),sngl(tau(1,k))
 24   m=1
 26   if(jndx(j,1).ge.indx(i,1)) go to 25
      j=j+1
      go to 26
 25   jndx(j,2)=min0(jidx(j),k)
      if(jndx(j,1).lt.jndx(j,2)) go to 37
      jndx(j,2)=-1
      go to 20
c37   write(10,*)'j jndx jidx',j,jndx(j,1),jndx(j,2),jidx(j),' ',
c    1 phcd(j)
 37   do 30 l=1,2
 28   if(dabs(pux(m,nph)-px(j,l)).le.dtol) go to 27
      if(m.ge.mu) go to 29
      m=m+1
      go to 28
 27   xbrn(j,l)=xt(j,l)+sgn*xc(m)
c     write(10,*)'x up:  j l m  ',j,l,m
      go to 30
 29   xbrn(j,l)=fac*(xus1(iph)+xus2(iph)+xus1(kph)+xus2(kph))+
     1 sgn*xus1(nph)
c     write(10,*)'x up:  j l end',j,l
c     write(10,*)'&&&&& nph iph kph xus1 xus2 xbrn =',
c    1 nph,iph,kph,sngl(xus1(1)),sngl(xus1(2)),sngl(xus2(1)),
c    2 sngl(xus2(2)),sngl(xbrn(j,l))
 30   continue
      if(j.ge.nbrn) go to 20
      j=j+1
      if(jndx(j,1).le.k) go to 25
 20   continue
      return
      end

      double precision function umod(zs,isrc,nph)
      save 
      include 'ttlim.inc'
      character*31 msg
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn
      double precision zs,uend,dtol,zmod
      dimension isrc(2)
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      data dtol/1d-6/
c
      m1=mt(nph)
      do 1 i=2,m1
      if(zm(i,nph).le.zs) go to 2
 1    continue
      dep=(1d0-dexp(zs))/xn
      write(msg,100)dep
      write(6,100)dep
 100  format('Source depth (',f6.1,') too deep.')
      write(6,*)msg
      call abort
 2    if(dabs(zs-zm(i,nph)).le.dtol.and.dabs(zm(i,nph)-zm(i+1,nph)).le.
     1 dtol) go to 3
      j=i-1
      isrc(nph)=j
      umod=pm(j,nph)+(pm(i,nph)-pm(j,nph))*(dexp(zs-zm(j,nph))-1d0)/
     1 (dexp(zm(i,nph)-zm(j,nph))-1d0)
      return
 3    isrc(nph)=i
      umod=pm(i+1,nph)
      return
c
      entry zmod(uend,js,nph)
      i=js+1
      zmod=zm(js,nph)+dlog(dmax1((uend-pm(js,nph))*(dexp(zm(i,nph)-
     1 zm(js,nph))-1d0)/(pm(i,nph)-pm(js,nph))+1d0,1d-30))
      return
      end

      subroutine bkin(lu,nrec,len,buf)
c
c $$$$$ calls no other routines $$$$$
c
c   Bkin reads a block of len double precision words into array buf(len)
c   from record nrec of the direct access unformatted file connected to
c   logical unit lu.
c
      save
      double precision buf(len),tmp
c
      if(nrec.le.0) go to 1
      read(lu,rec=nrec)buf
      tmp=buf(1)
      return
c   If the record doesn't exist, zero fill the buffer.
 1    do 2 i=1,len
 2    buf(i)=0d0
      return
      end

      subroutine tauint(ptk,ptj,pti,zj,zi,tau,x)
      save
c
c $$$$$ calls warn $$$$$
c
c   Tauint evaluates the intercept (tau) and distance (x) integrals  for
c   the spherical earth assuming that slowness is linear between radii
c   for which the model is known.  The partial integrals are performed
c   for ray slowness ptk between model radii with slownesses ptj and pti
c   with equivalent flat earth depths zj and zi respectively.  The partial
c   integrals are returned in tau and x.  Note that ptk, ptj, pti, zj, zi,
c   tau, and x are all double precision.
c
      character*71 msg
      double precision ptk,ptj,pti,zj,zi,tau,x
      double precision xx,b,sqk,sqi,sqj,sqb
c
      if(dabs(zj-zi).le.1d-9) go to 13
      if(dabs(ptj-pti).gt.1d-9) go to 10
      if(dabs(ptk-pti).le.1d-9) go to 13
      b=dabs(zj-zi)
      sqj=dsqrt(dabs(ptj*ptj-ptk*ptk))
      tau=b*sqj
      x=b*ptk/sqj
      go to 4
 10   if(ptk.gt.1d-9.or.pti.gt.1d-9) go to 1
c   Handle the straight through ray.
      tau=ptj
      x=1.5707963267948966d0
      go to 4
 1    b=ptj-(pti-ptj)/(dexp(zi-zj)-1d0)
      if(ptk.gt.1d-9) go to 2
      tau=-(pti-ptj+b*dlog(pti/ptj)-b*dlog(dmax1((ptj-b)*pti/
     1 ((pti-b)*ptj),1d-30)))
      x=0d0
      go to 4
 2    if(ptk.eq.pti) go to 3
      if(ptk.eq.ptj) go to 11
      sqk=ptk*ptk
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(sqb.gt.1d-30) go to 5
      xx=0d0
      x=ptk*(dsqrt(dabs((pti+b)/(pti-b)))-dsqrt(dabs((ptj+b)/
     1 (ptj-b))))/b
      go to 6
 5    if(b*b.lt.sqk) go to 7
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*
     1 (sqb*sqj+b*ptj-sqk)),1d-30))
      x=ptk*xx/sqb
      go to 6
 7    xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptk*dabs(pti-b)),1d0),-1d0))-
     1 dasin(dmax1(dmin1((b*ptj-sqk)/(ptk*dabs(ptj-b)),1d0),-1d0))
      x=-ptk*xx/sqb
 6    tau=-(sqi-sqj+b*dlog((pti+sqi)/(ptj+sqj))-sqb*xx)
      go to 4
 3    sqk=pti*pti
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 8
      xx=dlog(dmax1((ptj-b)*(b*pti-sqk)/((pti-b)*(sqb*sqj+b*ptj-sqk)),
     1 1d-30))
      x=pti*xx/sqb
      go to 9
 8    xx=dsign(1.5707963267948966d0,b-pti)-dasin(dmax1(dmin1((b*ptj-
     1 sqk)/(pti*dabs(ptj-b)),1d0),-1d0))
      x=-pti*xx/sqb
 9    tau=-(b*dlog(pti/(ptj+sqj))-sqj-sqb*xx)
      go to 4
 11   sqk=ptj*ptj
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 12
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*(b*ptj-sqk)),
     1 1d-30))
      x=ptj*xx/sqb
      go to 14
 12   xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptj*dabs(pti-b)),1d0),-1d0))-
     1 dsign(1.5707963267948966d0,b-ptj)
      x=-ptj*xx/sqb
 14   tau=-(b*dlog((pti+sqi)/ptj)+sqi-sqb*xx)
c
c   Handle various error conditions.
c
 4    if(x.ge.-1d-10) go to 15
      write(msg,100)ptk,ptj,pti,tau,x
 100  format('Bad range: ',1p5d12.4)
      call warn(msg)
 15   if(tau.ge.-1d-10) go to 16
      write(msg,101)ptk,ptj,pti,tau,x
 101  format('Bad tau: ',1p5d12.4)
      call warn(msg(1:69))
 16   return
c   Trap null integrals and handle them properly.
 13   tau=0d0
      x=0d0
      return
      end


      subroutine spfit(jb,int)
      save 
      include 'ttlim.inc'
      character*3 disc
      character*8 phcd
      logical newgrd,makgrd,segmsk,prnt
c     logical log
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
      double precision pmn,dmn,dmx,hm,shm,thm,p0,p1,tau0,tau1,x0,x1,pe,
     1 pe0,spe0,scpe0,pe1,spe1,scpe1,dpe,dtau,dbrnch,cn,x180,x360,dtol,
     2 ptol,xmin
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      common/prtflc/segmsk(jseg),prnt(2)
      data dbrnch,cn,x180,x360,xmin,dtol,ptol/2.5307274d0,57.295779d0,
     1 3.1415927d0,6.2831853d0,3.92403d-3,1d-6,2d-6/
c
      if(prnt(1)) write(10,102)
      i1=jndx(jb,1)
      i2=jndx(jb,2)
c     write(10,*)'Spfit:  jb i1 i2 pt =',jb,i1,i2,sngl(pt(i1)),
c    1 sngl(pt(i2))
      if(i2-i1.gt.1.or.dabs(pt(i2)-pt(i1)).gt.ptol) go to 14
      jndx(jb,2)=-1
      return
 14   newgrd=.false.
      makgrd=.false.
      if(dabs(px(jb,2)-pt(i2)).gt.dtol) newgrd=.true.
c     write(10,*)'Spfit:  px newgrd =',sngl(px(jb,2)),newgrd
      if(.not.newgrd) go to 10
      k=mod(int-1,2)+1
      if(int.ne.int0(k)) makgrd=.true.
c     write(10,*)'Spfit:  int k int0 makgrd =',int,k,int0(k),makgrd
      if(int.gt.2) go to 12
c     call query('Enter xmin:',log)
c     read *,xmin
c     xmin=xmin*xn
      xmin=xn*amin1(amax1(2.*odep,2.),25.)
c     write(10,*)'Spfit:  xmin =',xmin,xmin/xn
      call pdecu(i1,i2,xbrn(jb,1),xbrn(jb,2),xmin,int,i2)
      jndx(jb,2)=i2
 12   nn=i2-i1+1
      if(makgrd) call tauspl(1,nn,pt(i1),tcoef(1,1,k))
c     write(10,301,iostat=ios)jb,k,nn,int,newgrd,makgrd,
c    1 xbrn(jb,1),xbrn(jb,2),(i,pt(i-1+i1),tau(1,i-1+i1),
c    2 (tcoef(j,i,k),j=1,5),i=1,nn)
c301  format(/1x,4i3,2l3,2f12.8/(1x,i5,0p2f12.8,1p5d10.2))
      call fitspl(1,nn,tau(1,i1),xbrn(jb,1),xbrn(jb,2),tcoef(1,1,k))
      int0(k)=int
      go to 11
 10   call fitspl(i1,i2,tau,xbrn(jb,1),xbrn(jb,2),coef)
 11   pmn=pt(i1)
      dmn=xbrn(jb,1)
      dmx=dmn
      mxcnt=0
      mncnt=0
c     call appx(i1,i2,xbrn(jb,1),xbrn(jb,2))
c     write(10,300)(i,pt(i),(tau(j,i),j=1,3),i=i1,i2)
c300  format(/(1x,i5,4f12.6))
      pe=pt(i2)
      p1=pt(i1)
      tau1=tau(1,i1)
      x1=tau(2,i1)
      pe1=pe-p1
      spe1=dsqrt(dabs(pe1))
      scpe1=pe1*spe1
      j=i1
      is=i1+1
      do 2 i=is,i2
      p0=p1
      p1=pt(i)
      tau0=tau1
      tau1=tau(1,i)
      x0=x1
      x1=tau(2,i)
      dpe=p0-p1
      dtau=tau1-tau0
      pe0=pe1
      pe1=pe-p1
      spe0=spe1
      spe1=dsqrt(dabs(pe1))
      scpe0=scpe1
      scpe1=pe1*spe1
      tau(4,j)=(2d0*dtau-dpe*(x1+x0))/(.5d0*(scpe1-scpe0)-1.5d0*spe1*
     1 spe0*(spe1-spe0))
      tau(3,j)=(dtau-dpe*x0-(scpe1+.5d0*scpe0-1.5d0*pe1*spe0)*tau(4,j))/
     1 (dpe*dpe)
      tau(2,j)=(dtau-(pe1*pe1-pe0*pe0)*tau(3,j)-(scpe1-scpe0)*tau(4,j))/
     1 dpe
      tau(1,j)=tau0-scpe0*tau(4,j)-pe0*(pe0*tau(3,j)+tau(2,j))
      xlim(1,j)=dmin1(x0,x1)
      xlim(2,j)=dmax1(x0,x1)
      if(xlim(1,j).ge.dmn) go to 5
      dmn=xlim(1,j)
      pmn=pt(j)
      if(x1.lt.x0) pmn=pt(i)
 5    disc=' '
      if(dabs(tau(3,j)).le.1d-30) go to 4
      shm=-.375d0*tau(4,j)/tau(3,j)
      hm=shm*shm
      if(shm.le.0d0.or.(hm.le.pe1.or.hm.ge.pe0)) go to 4
 7    thm=tau(2,j)+shm*(2d0*shm*tau(3,j)+1.5d0*tau(4,j))
      xlim(1,j)=dmin1(xlim(1,j),thm)
      xlim(2,j)=dmax1(xlim(2,j),thm)
      if(thm.ge.dmn) go to 6
      dmn=thm
      pmn=pe-hm
 6    disc='max'
      if(tau(4,j).lt.0d0) disc='min'
      if(disc.eq.'max') mxcnt=mxcnt+1
      if(disc.eq.'min') mncnt=mncnt+1
 4    if(prnt(1)) write(10,100,iostat=ios)disc,j,pt(j),
     1 (tau(k,j),k=1,4),(cn*xlim(k,j),k=1,2)
 100  format(1x,a,i5,f10.6,1p4e10.2,0p2f7.2)
      dmx=dmax1(dmx,xlim(2,j))
 2    j=i
c     if(prnt(1)) write(10,100,iostat=ios)'   ',j,pt(j)
      xbrn(jb,1)=dmn
      xbrn(jb,2)=dmx
      xbrn(jb,3)=pmn
      idel(jb,1)=1
      idel(jb,2)=1
      if(xbrn(jb,1).gt.x180) idel(jb,1)=2
      if(xbrn(jb,2).gt.x180) idel(jb,2)=2
      if(xbrn(jb,1).gt.x360) idel(jb,1)=3
      if(xbrn(jb,2).gt.x360) idel(jb,2)=3
      if(int.gt.2) go to 1
      phcd(jb)=phcd(jb)(1:1)
      i=jb
      do 8 j=1,nbrn
      i=mod(i,nbrn)+1
      if(phcd(i)(1:1).eq.phcd(jb).and.phcd(i)(2:2).ne.'P'.and.
     1 (pe.ge.px(i,1).and.pe.le.px(i,2))) go to 9
 8    continue
      go to 1
 9    phcd(jb)=phcd(i)
      if(dabs(pt(i2)-pt(jndx(i,1))).le.dtol) phcd(jb)=phcd(i-1)
 1    if(prnt(1).and.prnt(2)) write(10,102)
 102  format()
      if(dbrn(jb,1).le.0d0) go to 3
      dbrn(jb,1)=dmx
      dbrn(jb,2)=dbrnch
      if(prnt(2)) write(10,101,iostat=ios)phcd(jb),
     1 (jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),
     2 (cn*dbrn(jb,k),k=1,2),(idel(jb,k),k=1,3),int,newgrd,makgrd
 101  format(1x,a,2i5,2f8.2,f8.4,2f8.2,4i3,2l2)
      go to 15
 3    if(prnt(2)) write(10,103,iostat=ios)phcd(jb),
     1 (jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),
     2 (idel(jb,k),k=1,3),int,newgrd,makgrd
 103  format(1x,a,2i5,2f8.2,f8.4,16x,4i3,2l2)
 15   if(mxcnt.gt.mncnt.or.mncnt.gt.mxcnt+1)
     1 call warn('Bad interpolation on '//phcd(jb))
      return
      end
c
c-----------------------------------------------------------------
c
      subroutine pdecu(i1,i2,x0,x1,xmin,int,len)
      save 
      include 'ttlim.inc'
      double precision us,pt,tau,xlim,xbrn,dbrn,ua,taua
      double precision x0,x1,xmin,dx,dx2,sgn,rnd,xm,axm,x,h1,h2,hh,xs
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/pdec/ua(5,2),taua(5,2),deplim,ka
c
c     write(10,*)'Pdecu:  ua =',sngl(ua(1,int))
      if(ua(1,int).le.0d0) go to 17
c     write(10,*)'Pdecu:  fill in new grid'
      k=i1+1
      do 18 i=1,ka
      pt(k)=ua(i,int)
      tau(1,k)=taua(i,int)
 18   k=k+1
      pt(k)=pt(i2)
      tau(1,k)=tau(1,i2)
      go to 19
c
 17   is=i1+1
      ie=i2-1
      xs=x1
      do 11 i=ie,i1,-1
      x=xs
      if(i.ne.i1) go to 12
      xs=x0
      go to 14
 12   h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      xs=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 14   if(dabs(x-xs).le.xmin) go to 15
 11   continue
      len=i2
      return
 15   ie=i
      if(dabs(x-xs).gt..75d0*xmin.or.ie.eq.i2) go to 16
      xs=x
      ie=ie+1
 16   n=max0(idint(dabs(xs-x0)/xmin+.8d0),1)
      dx=(xs-x0)/n
      dx2=dabs(.5d0*dx)
      sgn=dsign(1d0,dx)
      rnd=0d0
      if(sgn.gt.0d0) rnd=1d0
      xm=x0+dx
      k=i1
      m=is
      axm=1d10
      do 1 i=is,ie
      if(i.lt.ie) go to 8
      x=xs
      go to 5
 8    h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      x=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 5    if(sgn*(x-xm).le.dx2) go to 2
      if(k.lt.m) go to 3
      do 4 j=m,k
 4    pt(j)=-1d0
 3    m=k+2
      k=i-1
      axm=1d10
 7    xm=xm+dx*idint((x-xm-dx2)/dx+rnd)
 2    if(dabs(x-xm).ge.axm) go to 1
      axm=dabs(x-xm)
      k=i-1
 1    continue
      if(k.lt.m) go to 9
      do 6 j=m,k
 6    pt(j)=-1d0
 9    k=i1
      do 10 i=is,i2
      if(pt(i).lt.0d0) go to 10
      k=k+1
      pt(k)=pt(i)
      tau(1,k)=tau(1,i)
 10   continue
 19   len=k
c     write(10,300)(i,pt(i),tau(1,i),i=i1,len)
c300  format(/(1x,i5,0pf12.6,1pd15.4))
      return
      end
      subroutine tauspl(i1,i2,pt,coef)
c
c $$$$$ calls only library routines $$$$$
c
c   Given ray parameter grid pt;i (pt sub i), i=i1,i1+1,...,i2, tauspl
c   determines the i2-i1+3 basis functions for interpolation I such
c   that:
c
c      tau(p) = a;1,i + Dp * a;2,i + Dp**2 * a;3,i + Dp**(3/2) * a;4,i
c
c   where Dp = pt;n - p, pt;i <= p < pt;i+1, and the a;j,i's are
c   interpolation coefficients.  Rather than returning the coefficients,
c   a;j,i, which necessarily depend on tau(pt;i), i=i1,i1+1,...,i2 and
c   x(pt;i) (= -d tau(p)/d p | pt;i), i=i1,i2, tauspl returns the
c   contribution of each basis function and its derivitive at each
c   sample.  Each basis function is non-zero at three grid points,
c   therefore, each grid point will have contributions (function values
c   and derivitives) from three basis functions.  Due to the basis
c   function normalization, one of the function values will always be
c   one and is not returned in array coef with the other values.
c   Rewritten on 23 December 1983 by R. Buland.
c
      save
      double precision pt(i2),coef(5,i2)
      double precision del(5),sdel(5),deli(5),d3h(4),d1h(4),dih(4),
     1 d(4),ali,alr,b3h,b1h,bih,th0p,th2p,th3p,th2m
c
      n2=i2-i1-1
      if(n2.le.-1) return
      is=i1+1
c
c   To achieve the requisite stability, proceed by constructing basis
c   functions G;i, i=0,1,...,n+1.  G;i will be non-zero only on the
c   interval [p;i-2,p;i+2] and will be continuous with continuous first
c   and second derivitives.  G;i(p;i-2) and G;i(p;i+2) are constrained
c   to be zero with zero first and second derivitives.  G;i(p;i) is
c   normalized to unity.
c
c   Set up temporary variables appropriate for G;-1.  Note that to get
c   started, the ray parameter grid is extrapolated to yeild p;i, i=-2,
c   -1,0,1,...,n.
      del(2)=pt(i2)-pt(i1)+3d0*(pt(is)-pt(i1))
      sdel(2)=dsqrt(dabs(del(2)))
      deli(2)=1d0/sdel(2)
      m=2
      do 1 k=3,5
      del(k)=pt(i2)-pt(i1)+(5-k)*(pt(is)-pt(i1))
      sdel(k)=dsqrt(dabs(del(k)))
      deli(k)=1d0/sdel(k)
      d3h(m)=del(k)*sdel(k)-del(m)*sdel(m)
      d1h(m)=sdel(k)-sdel(m)
      dih(m)=deli(k)-deli(m)
 1    m=k
      l=i1-1
      if(n2.le.0) go to 10
c   Loop over G;i, i=0,1,...,n-3.
      do 2 i=1,n2
      m=1
c   Update temporary variables for G;i-1.
      do 3 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 3
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 3    m=k
      l=l+1
      del(5)=pt(i2)-pt(l+1)
      sdel(5)=dsqrt(dabs(del(5)))
      deli(5)=1d0/sdel(5)
      d3h(4)=del(5)*sdel(5)-del(4)*sdel(4)
      d1h(4)=sdel(5)-sdel(4)
      dih(4)=deli(5)-deli(4)
c   Construct G;i-1.
      ali=1d0/(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*dih(1)*del(3))*
     1 del(3))
      alr=ali*(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*del(3)*
     1 deli(2)-sdel(3))*del(3))
      b3h=d3h(2)+alr*d3h(1)
      b1h=d1h(2)+alr*d1h(1)
      bih=dih(2)+alr*dih(1)
      th0p=d1h(1)*b3h-d3h(1)*b1h
      th2p=d1h(3)*b3h-d3h(3)*b1h
      th3p=d1h(4)*b3h-d3h(4)*b1h
      th2m=dih(3)*b3h-d3h(3)*bih
c   The d;i's completely define G;i-1.
      d(4)=ali*((dih(1)*b3h-d3h(1)*bih)*th2p-th2m*th0p)/((dih(4)*b3h-
     1 d3h(4)*bih)*th2p-th2m*th3p)
      d(3)=(th0p*ali-th3p*d(4))/th2p
      d(2)=(d3h(1)*ali-d3h(3)*d(3)-d3h(4)*d(4))/b3h
      d(1)=alr*d(2)-ali
c   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
c   G;i-1(p;i-1) need not be constructed as it is normalized to unity.
      coef(1,l)=(.125d0*del(5)*sdel(5)-(.75d0*sdel(5)+.375d0*deli(5)*
     1 del(4)-sdel(4))*del(4))*d(4)
      if(i.ge.3) coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+
     1 .375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)
c   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
      coef(3,l)=-.75d0*(sdel(5)+deli(5)*del(4)-2d0*sdel(4))*d(4)
      if(i.ge.2) coef(4,l-1)=-.75d0*((sdel(2)+deli(2)*del(3)-
     1 2d0*sdel(3))*d(2)-(d1h(1)+dih(1)*del(3))*d(1))
      if(i.ge.3) coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-
     1 2d0*sdel(2))*d(1)
 2    continue
c   Loop over G;i, i=n-2,n-1,n,n+1.  These cases must be handled
c   seperately because of the singularities in the second derivitive
c   at p;n.
 10   do 4 j=1,4
      m=1
c   Update temporary variables for G;i-1.
      do 5 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 5
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 5    m=k
      l=l+1
      del(5)=0d0
      sdel(5)=0d0
      deli(5)=0d0
c   Construction of the d;i's is different for each case.  In cases
c   G;i, i=n-1,n,n+1, G;i is truncated at p;n to avoid patching across
c   the singularity in the second derivitive.
      if(j.lt.4) go to 6
c   For G;n+1 constrain G;n+1(p;n) to be .25.
      d(1)=2d0/(del(1)*sdel(1))
      go to 9
c   For G;i, i=n-2,n-1,n, the condition dG;i(p)/dp|p;i = 0 has been
c   substituted for the second derivitive continuity condition that
c   can no longer be satisfied.
 6    alr=(sdel(2)+deli(2)*del(3)-2d0*sdel(3))/(d1h(1)+dih(1)*del(3))
      d(2)=1d0/(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*deli(2)*
     1 del(3)-sdel(3))*del(3)-(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*
     2 dih(1)*del(3))*del(3))*alr)
      d(1)=alr*d(2)
      if(j-2)8,7,9
c   For G;n-1 constrain G;n-1(p;n) to be .25.
 7    d(3)=(2d0+d3h(2)*d(2)+d3h(1)*d(1))/(del(3)*sdel(3))
      go to 9
c   No additional constraints are required for G;n-2.
 8    d(3)=-((d3h(2)-d1h(2)*del(4))*d(2)+(d3h(1)-d1h(1)*del(4))*
     1 d(1))/(d3h(3)-d1h(3)*del(4))
      d(4)=(d3h(3)*d(3)+d3h(2)*d(2)+d3h(1)*d(1))/(del(4)*sdel(4))
c   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
 9    if(j.le.2) coef(1,l)=(.125d0*del(3)*sdel(3)-(.75d0*sdel(3)+.375d0*
     1 deli(3)*del(4)-sdel(4))*del(4))*d(3)-(.125d0*d3h(2)-(.75d0*
     2 d1h(2)+.375d0*dih(2)*del(4))*del(4))*d(2)-(.125d0*d3h(1)-(.75d0*
     3 d1h(1)+.375d0*dih(1)*del(4))*del(4))*d(1)
      if(l-i1.gt.1) coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+
     1 .375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)
c   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
      if(j.le.2) coef(3,l)=-.75d0*((sdel(3)+deli(3)*del(4)-
     1 2d0*sdel(4))*d(3)-(d1h(2)+dih(2)*del(4))*d(2)-(d1h(1)+
     2 dih(1)*del(4))*d(1))
      if(j.le.3.and.l-i1.gt.0) coef(4,l-1)=0d0
      if(l-i1.gt.1) coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-
     1 2d0*sdel(2))*d(1)
 4    continue
      return
      end
      subroutine fitspl(i1,i2,tau,x1,xn,coef)
c
c $$$$$ calls only library routines $$$$$
c
c   Given ray parameter grid p;i (p sub i), i=1,2,...,n, corresponding
c   tau;i values, and x;1 and x;n (x;i = -dtau/dp|p;i); tauspl finds
c   interpolation I such that:  tau(p) = a;1,i + Dp * a;2,i + Dp**2 *
c   a;3,i + Dp**(3/2) * a;4,i where Dp = p;n - p and p;i <= p < p;i+1.
c   Interpolation I has the following properties:  1) x;1, x;n, and
c   tau;i, i=1,2,...,n are fit exactly, 2) the first and second
c   derivitives with respect to p are continuous everywhere, and
c   3) because of the paramaterization d**2 tau/dp**2|p;n is infinite.
c   Thus, interpolation I models the asymptotic behavior of tau(p)
c   when tau(p;n) is a branch end due to a discontinuity in the
c   velocity model.  Note that array a must be dimensioned at least
c   a(4,n) though the interpolation coefficients will be returned in
c   the first n-1 columns.  The remaining column is used as scratch
c   space and returned as all zeros.  Programmed on 16 August 1982 by
c   R. Buland.
c
      save 
      double precision tau(4,i2),x1,xn,coef(5,i2),a(2,100),ap(3),
     1 b(100),alr,g1,gn
c
      if(i2-i1)13,1,2
 1    tau(2,i1)=x1
 13   return
 2    n=0
      do 3 i=i1,i2
      n=n+1
      b(n)=tau(1,i)
      do 3 j=1,2
 3    a(j,n)=coef(j,i)
      do 4 j=1,3
 4    ap(j)=coef(j+2,i2)
      n1=n-1
c
c   Arrays ap(*,1), a, and ap(*,2) comprise n+2 x n+2 penta-diagonal
c   matrix A.  Let x1, tau, and xn comprise corresponding n+2 vector b.
c   Then, A * g = b, may be solved for n+2 vector g such that
c   interpolation I is given by I(p) = sum(i=0,n+1) g;i * G;i(p).
c
c   Eliminate the lower triangular portion of A to form A'.  A
c   corresponding transformation applied to vector b is stored in
c   a(4,*).
      alr=a(1,1)/coef(3,i1)
      a(1,1)=1d0-coef(4,i1)*alr
      a(2,1)=a(2,1)-coef(5,i1)*alr
      b(1)=b(1)-x1*alr
      j=1
      do 5 i=2,n
      alr=a(1,i)/a(1,j)
      a(1,i)=1d0-a(2,j)*alr
      b(i)=b(i)-b(j)*alr
 5    j=i
      alr=ap(1)/a(1,n1)
      ap(2)=ap(2)-a(2,n1)*alr
      gn=xn-b(n1)*alr
      alr=ap(2)/a(1,n)
c   Back solve the upper triangular portion of A' for coefficients g;i.
c   When finished, storage g(2), a(4,*), g(5) will comprise vector g.
      gn=(gn-b(n)*alr)/(ap(3)-a(2,n)*alr)
      b(n)=(b(n)-gn*a(2,n))/a(1,n)
      j=n
      do 6 i=n1,1,-1
      b(i)=(b(i)-b(j)*a(2,i))/a(1,i)
 6    j=i
      g1=(x1-coef(4,i1)*b(1)-coef(5,i1)*b(2))/coef(3,i1)
c
      tau(2,i1)=x1
      is=i1+1
      ie=i2-1
      j=1
      do 7 i=is,ie
      j=j+1
 7    tau(2,i)=coef(3,i)*b(j-1)+coef(4,i)*b(j)+coef(5,i)*b(j+1)
      tau(2,i2)=xn
      return
      end
      subroutine trtm(delta,max,n,tt,dtdd,dtdh,dddp,phnm)
      save 
      include 'ttlim.inc'
      character*(*) phnm(max)
      character*8 ctmp(60)
      dimension tt(max),dtdd(max),dtdh(max),dddp(max),tmp(60,4),
     1 iptr(60)
      double precision us,pt,tau,xlim,xbrn,dbrn
      double precision x(3),cn,dtol,pi,pi2
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      data cn,dtol,atol,pi,pi2/.017453292519943296d0,1d-6,.005,
     1 3.1415926535897932d0,6.2831853071795865d0/
c
c      if (delta < 1.0) then
c      do i=1,jout
c         write(32,*) i,tau(1,i)
c      end do
c      stop
c      endif
      n=0
      if(mbr2.le.0) return
      x(1)=dmod(dabs(cn*delta),pi2)
      if(x(1).gt.pi) x(1)=pi2-x(1)
      x(2)=pi2-x(1)
      x(3)=x(1)+pi2
      if(dabs(x(1)).gt.dtol) go to 9
      x(1)=dtol
      x(3)=-10d0
 9    if(dabs(x(1)-pi).gt.dtol) go to 7
      x(1)=pi-dtol
      x(2)=-10d0
 7    do 1 j=mbr1,mbr2
 1    if(jndx(j,2).gt.0) call findtt(j,x,max,n,tmp,tmp(1,2),tmp(1,3),
     1 tmp(1,4),ctmp)
      if(n-1)3,4,5
 4    iptr(1)=1
      go to 6
 5    call r4sort(n,tmp,iptr)
 6    k=0
      do 2 i=1,n
      j=iptr(i)
      if(k.le.0) go to 8
      if(phnm(k).eq.ctmp(j).and.abs(tt(k)-tmp(j,1)).le.atol) go to 2
 8    k=k+1
      tt(k)=tmp(j,1)
      dtdd(k)=tmp(j,2)
      dtdh(k)=tmp(j,3)
      dddp(k)=tmp(j,4)
      phnm(k)=ctmp(j)
 2    continue
      n=k
 3    return
      end
      subroutine findtt(jb,x0,max,n,tt,dtdd,dtdh,dddp,phnm)
      save 
      include 'ttlim.inc'
      character*(*) phnm(max)
      character*8 phcd
      character*67 msg
      dimension tt(max),dtdd(max),dtdh(max),dddp(max)
      double precision us,pt,tau,xlim,xbrn,dbrn
      double precision x,x0(3),p0,p1,arg,dp,dps,delp,tol,ps,deps
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/pcdc/phcd(jbrn)
      data tol/3d-6/,deps/1d-10/
c
      nph=iabs(idel(jb,3))
      hsgn=isign(1,idel(jb,3))*hn
      dsgn=(-1.)**idel(jb,1)*dn
      dpn=-1./tn
      do 10 ij=idel(jb,1),idel(jb,2)
      x=x0(ij)
      dsgn=-dsgn
      if(x.lt.xbrn(jb,1).or.x.gt.xbrn(jb,2)) go to 12
      j=jndx(jb,1)
      is=j+1
      ie=jndx(jb,2)
      do 1 i=is,ie
      if(x.le.xlim(1,j).or.x.gt.xlim(2,j)) go to 8
      le=n
      p0=pt(ie)-pt(j)
      p1=pt(ie)-pt(i)
      delp=dmax1(tol*(pt(i)-pt(j)),1d-3)
      if(dabs(tau(3,j)).gt.1d-30) go to 2
      dps=(x-tau(2,j))/(1.5d0*tau(4,j))
      dp=dsign(dps*dps,dps)
      dp0=dp
      if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 9
      if(n.ge.max) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*(tau(1,j)+dp*(tau(2,j)+dps*tau(4,j))+ps*x)
      dtdd(n)=dsgn*ps
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*.75d0*tau(4,j)/dmax1(dabs(dps),deps)
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 8
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
      go to 8
 2    do 4 jj=1,2
      go to (5,6),jj
 5    arg=9d0*tau(4,j)*tau(4,j)+32d0*tau(3,j)*(x-tau(2,j))
      if(arg.ge.0d0) go to 3
      write(msg,100)arg
 100  format('Bad sqrt argument:',1pd11.2,'.')
      call warn(msg(1:30))
 3    dps=-(3d0*tau(4,j)+dsign(dsqrt(dabs(arg)),tau(4,j)))/(8d0*
     1 tau(3,j))
      dp=dsign(dps*dps,dps)
      dp0=dp
      go to 7
 6    dps=(tau(2,j)-x)/(2d0*tau(3,j)*dps)
      dp=dsign(dps*dps,dps)
 7    if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 4
      if(n.ge.max) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+ps*x)
      dtdd(n)=dsgn*ps
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*(2d0*tau(3,j)+.75d0*tau(4,j)/dmax1(dabs(dps),deps))
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 4
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
 4    continue
 9    if(n.gt.le) go to 8
      write(msg,101)phcd(jb),x,dp0,dp,p1,p0
 101  format('Failed to find phase:  ',a,f8.1,4f7.4)
      call warn(msg)
 8    j=i
 1    continue
c
 12   if(x.lt.dbrn(jb,1).or.x.gt.dbrn(jb,2)) go to 10
      if(n.ge.max) go to 13
      j=jndx(jb,1)
      i=jndx(jb,2)
      dp=pt(i)-pt(j)
      dps=dsqrt(dabs(dp))
      n=n+1
      tt(n)=tn*(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+
     1 pt(j)*x)
      dtdd(n)=dsgn*sngl(pt(j))
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-pt(j)*pt(j))))
      dddp(n)=dpn*(2d0*tau(3,j)+.75d0*tau(4,j)/dmax1(dps,deps))
      ln=index(phcd(jb),' ')-1
      if(ln.le.0) ln=len(phcd(jb))
      phnm(n)=phcd(jb)(1:ln)//'diff'
 10   continue
      return
 13   write(msg,102)max
 102  format('More than',i3,' arrivals found.')
      call warn(msg(1:28))
      return
      end
      subroutine asnag1(lu,mode,n,ia,ib)
c
c $$$$$ calls assign, iargc, and getarg $$$$$
c
c   Asnag1 assigns logical unit lu to a direct access disk file
c   with mode "mode" and record length "len".  See dasign for 
c   details.  The n th argument is used as the model name.  If there 
c   is no n th argument and ib is non-blank, it is taken to be the 
c   model name.  If ib is blank, the user is prompted for the
c   model name using the character string in variable ia as the
c   prompt.  Programmed on 8 October 1980 by R. Buland.
c
      save
      logical log
      character*(*) ia,ib
      character cdum*30
c
c      if(iargc().lt.n) go to 1
c      call getarg(n,ib)
c      go to 2
c
c 1    if(ib.ne.' ') go to 2
c      call query(ia,log)
c      read(*,100)ib
c 100  format(a)
c
c 2    nb=index(ib,' ')-1
c      if(nb.le.0) nb=len(ib)
c      cdum=ib(1:nb)//'.hed'
      cdum='ak135.hed'
      call assign(lu,mode,cdum)
      return
      end
      subroutine assign(lu,mode,ia)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine assign opens (connects) logical unit lu to the disk file
c   named by the character string ia with mode mode.  If iabs(mode) = 1,
c   then open the file for reading.  If iabs(mode) = 2, then open the
c   file for writing.  If iabs(mode) = 3, then open a scratch file for
c   writing.  If mode > 0, then the file is formatted.  If mode < 0,
c   then the file is unformatted.  All files opened by assign are
c   assumed to be sequential.  Programmed on 3 December 1979 by
c   R. Buland.
c
      save
      character*(*) ia
      logical exst
c
      if(mode.ge.0) nf=1
      if(mode.lt.0) nf=2
      ns=iabs(mode)
      if(ns.le.0.or.ns.gt.3) ns=3
      go to (1,2),nf
 1    go to (11,12,13),ns
 11   open(lu,file=ia,status='old',form='formatted')
      rewind lu
      return
 12   inquire(file=ia,exist=exst)
      if(exst) go to 11
 13   open(lu,file=ia,status='new',form='formatted')
      return
 2    go to (21,22,23),ns
 21   open(lu,file=ia,status='old',form='unformatted')
      rewind lu
      return
 22   inquire(file=ia,exist=exst)
      if(exst) go to 21
 23   open(lu,file=ia,status='new',form='unformatted')
      return
      end
      subroutine retrns(lu)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine retrns closes (disconnects) logical unit lu from the
c   calling program.  Programmed on 3 December 1979 by R. Buland.
c
      save
      close(unit=lu)
      return
      end
      subroutine query(ia,log)
c
c $$$$$ calls tnoua $$$$$
c
c   Subroutine query scans character string ia (up to 78 characters) for
c   a question mark or a colon.  It prints the string up to and
c   including the flag character plus two blanks with no newline on the
c   standard output.  If the flag was a question mark, query reads the
c   users response.  If the response is 'y' or 'yes', log is set to
c   true.  If the response is 'n' or 'no', log is set to false.  Any
c   other response causes the question to be repeated.  If the flag was
c   a colon, query simply returns allowing user input on the same line.
c   If there is no question mark or colon, the last non-blank character
c   is treated as if it were a colon.  If the string is null or all
c   blank, query prints an astrisk and returns.  Programmed on 3
c   December 1979 by R. Buland.
c
      save
      logical log
      character*(*) ia
      character*81 ib
      character*4 ans
      nn=len(ia)
      log=.true.
      ifl=1
      k=0
c   Scan ia for flag characters or end-of-string.
      do 1 i=1,nn
      ib(i:i)=ia(i:i)
      if(ib(i:i).eq.':') go to 7
      if(ib(i:i).eq.'?') go to 3
      if(ib(i:i).eq.'\0') go to 5
      if(ib(i:i).ne.' ') k=i
 1    continue
c   If we fell off the end of the string, branch if there were any non-
c   blank characters.
 5    if(k.gt.0) go to 6
c   Handle a null or all blank string.
      i=1
      ib(i:i)='*'
      go to 4
c   Handle a string with no question mark or colon but at least one
c   non-blank character.
 6    i=k
c   Append two blanks and print the string.
 7    i=i+2
      ib(i-1:i-1)=' '
      ib(i:i)=' '
c   Tnoua prints the first i characters of ib without a newline.
 4    call tnoua(ib,i)
      if(ifl.gt.0) return
c   If the string was a yes-no question read the response.
      read 102,ans
 102  format(a4)
      call uctolc(ans,-1)
c   If the response is yes log is already set properly.
      if(ans.eq.'y   '.or.ans.eq.'yes ') return
c   If the response is no set log to false.  Otherwise repeat the
c   question.
      if(ans.ne.'n   '.and.ans.ne.'no  ') go to 4
      log=.false.
      return
 3    ifl=-ifl
      go to 7
      end
      subroutine uctolc(ia,ifl)
c
c $$$$$ calls only library routines $$$$$
c
c   Subroutine uctolc converts alphabetic characters in string ia from
c   upper case to lower case.  If ifl < 0 all characters are converted.
c   Otherwise characters enclosed by single quotes are left unchanged.
c   Programmed on 21 January by R. Buland.  Calling sequence changed
c   on 11 December 1985 by R. Buland.
c
      character*(*) ia
      data nfl/1/
      if(ifl.lt.0) nfl=1
c   Scan the string.
      n=len(ia)
      do 1 i=1,n
      if(ifl.lt.0) go to 2
c   Look for single quotes.
      if(ia(i:i).eq.'''') nfl=-nfl
c   If we are in a quoted string skip the conversion.
      if(nfl.lt.0) go to 1
c   Do the conversion.
 2    ib=ichar(ia(i:i))
      if(ib.lt.65.or.ib.gt.90) go to 1
      ia(i:i)=char(ib+32)
 1    continue
      return
      end
      subroutine r4sort(n,rkey,iptr)
c
c $$$$$ calls no other routine $$$$$
c
c   R4sort sorts the n elements of array rkey so that rkey(i), 
c   i = 1, 2, 3, ..., n are in asending order.  R4sort is a trivial
c   modification of ACM algorithm 347:  "An efficient algorithm for
c   sorting with minimal storage" by R. C. Singleton.  Array rkey is
c   sorted in place in order n*alog2(n) operations.  Coded on
c   8 March 1979 by R. Buland.  Modified to handle real*4 data on
c   27 September 1983 by R. Buland.
c
      save
      dimension rkey(n),iptr(n),il(10),iu(10)
c   Note:  il and iu implement a stack containing the upper and
c   lower limits of subsequences to be sorted independently.  A
c   depth of k allows for n<=2**(k+1)-1.
      if(n.le.0) return
      do 1 i=1,n
 1    iptr(i)=i
      if(n.le.1) return
      r=.375
      m=1
      i=1
      j=n
c
c   The first section interchanges low element i, middle element ij,
c   and high element j so they are in order.
c
 5    if(i.ge.j) go to 70
 10   k=i
c   Use a floating point modification, r, of Singleton's bisection
c   strategy (suggested by R. Peto in his verification of the
c   algorithm for the ACM).
      if(r.gt..58984375) go to 11
      r=r+.0390625
      go to 12
 11   r=r-.21875
 12   ij=i+(j-i)*r
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 20
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 20   l=j
      if(rkey(iptr(j)).ge.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(j)
      iptr(j)=it
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 39   tmpkey=rkey(iptr(ij))
      go to 40
c
c   The second section continues this process.  K counts up from i and
c   l down from j.  Each time the k element is bigger than the ij
c   and the l element is less than the ij, then interchange the
c   k and l elements.  This continues until k and l meet.
c
 30   it=iptr(l)
      iptr(l)=iptr(k)
      iptr(k)=it
 40   l=l-1
      if(rkey(iptr(l)).gt.tmpkey) go to 40
 50   k=k+1
      if(rkey(iptr(k)).lt.tmpkey) go to 50
      if(k.le.l) go to 30
c
c   The third section considers the intervals i to l and k to j.  The
c   larger interval is saved on the stack (il and iu) and the smaller
c   is remapped into i and j for another shot at section one.
c
      if(l-i.le.j-k) go to 60
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 80
 60   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 80
c
c   The fourth section pops elements off the stack (into i and j).  If
c   necessary control is transfered back to section one for more
c   interchange sorting.  If not we fall through to section five.  Note
c   that the algorighm exits when the stack is empty.
c
 70   m=m-1
      if(m.eq.0) return
      i=il(m)
      j=iu(m)
 80   if(j-i.ge.11) go to 10
      if(i.eq.1) go to 5
      i=i-1
c
c   The fifth section is the end game.  Final sorting is accomplished
c   (within each subsequence popped off the stack) by rippling out
c   of order elements down to their proper positions.
c
 90   i=i+1
      if(i.eq.j) go to 70
      if(rkey(iptr(i)).le.rkey(iptr(i+1))) go to 90
      k=i
      kk=k+1
      ib=iptr(kk)
 100  iptr(kk)=iptr(k)
      kk=k
      k=k-1
      if(rkey(ib).lt.rkey(iptr(k))) go to 100
      iptr(kk)=ib
      go to 90
      end
      function iupcor(phnm,dtdd,xcor,tcor)
      save
      include 'ttlim.inc'
      character*(*) phnm
      character*8 phcd
      double precision us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,xu,
     1 px,xt,taut,coef,tauc,xc,tcoef,tp
      double precision x,dp,dps,ps,cn
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      data oldep,jp,js/-1.,2*0/,cn/57.295779d0/
c
      iupcor=1
c     print *,'oldep odep',oldep,odep
      if(oldep.eq.odep) go to 1
      oldep=odep
c   Find the upgoing P branch.
c     print *,'mbr1 mbr2',mbr1,mbr2
      do 2 jp=mbr1,mbr2
c     print *,'jp phcd xbrn',jp,'  ',phcd(jp),xbrn(jp,1)
      if((phcd(jp).eq.'Pg'.or.phcd(jp).eq.'Pb'.or.phcd(jp).eq.'Pn'.or.
     1 phcd(jp).eq.'P').and.xbrn(jp,1).le.0d0) go to 3
 2    continue
      jp=0
c   Find the upgoing S branch.
 3    do 4 js=mbr1,mbr2
c     print *,'js phcd xbrn',js,'  ',phcd(js),xbrn(js,1)
      if((phcd(js).eq.'Sg'.or.phcd(js).eq.'Sb'.or.phcd(js).eq.'Sn'.or.
     1 phcd(js).eq.'S').and.xbrn(js,1).le.0d0) go to 1
 4    continue
      js=0
c
c1    print *,'jp js',jp,js
 1    if(phnm.ne.'P'.and.phnm.ne.'p') go to 5
      jb=jp
      if(jb)14,14,6
c
 5    if(phnm.ne.'S'.and.phnm.ne.'s') go to 13
      jb=js
      if(jb)14,14,6
c
 6    is=jndx(jb,1)+1
      ie=jndx(jb,2)
      ps=abs(dtdd)/dn
c     print *,'jb is ie dtdd dn ps',jb,is,ie,dtdd,dn,ps
      if(ps.lt.pt(is-1).or.ps.gt.pt(ie)) go to 13
      do 7 i=is,ie
c     print *,'i pt',i,pt(i)
      if(ps.le.pt(i)) go to 8
 7    continue
      go to 13
c
 8    j=i-1
      dp=pt(ie)-ps
      dps=dsqrt(dabs(dp))
      x=tau(2,j)+2d0*dp*tau(3,j)+1.5d0*dps*tau(4,j)
c     print *,'j pt dp dps x',j,pt(ie),dp,dps,x
      tcor=tn*(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+ps*x)
      xcor=cn*x
c     print *,'iupcor xcor tcor',iupcor,xcor,tcor
      return
c
 13   iupcor=-1
 14   xcor=0.
      tcor=0.
c     print *,'iupcor xcor tcor',iupcor,xcor,tcor
      return
      end
      subroutine brnset(nn,pcntl,prflg)
c
c   Brnset takes character array pcntl(nn) as a list of nn tokens to be
c   used to select desired generic branches.  Prflg(3) is the old
c   prnt(2) debug print flags in the first two elements plus a new print
c   flag which controls a branch selected summary from brnset.  Note that
c   the original two flags controlled a list of all tau interpolations
c   and a branch range summary respectively.  The original summary output
c   still goes to logical unit 10 (ttim1.lis) while the new output goes
c   to the standard output (so the caller can see what happened).  Each
c   token of pcntl may be either a generic branch name (e.g., P, PcP,
c   PKP, etc.) or a keyword (defined in the data statement for cmdcd
c   below) which translates to more than one generic branch names.  Note
c   that generic branch names and keywords may be mixed.  The keywords
c   'all' (for all branches) and 'query' (for an interactive token input
c   query mode) are also available.
c
      save
      parameter(ncmd=4,lcmd=16)
      include 'ttlim.inc'
      logical prflg(3),segmsk,prnt,fnd,all
      character*(*) pcntl(nn)
      character*8 phcd,segcd(jbrn),cmdcd(ncmd),cmdlst(lcmd),phtmp,
     1 phlst(jseg)
      double precision zs,pk,pu,pux,tauu,xu,px,xt,taut,coef,tauc,xc,
     1 tcoef,tp
      dimension nsgpt(jbrn),ncmpt(2,ncmd)
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
c   Segmsk is a logical array that actually implements the branch
c   editing in depset and depcor.
      common/prtflc/segmsk(jseg),prnt(2)
c
c   The keywords do the following:
c      P      gives P-up, P, Pdiff, PKP, and PKiKP
c      P+     gives P-up, P, Pdiff, PKP, PKiKP, PcP, pP, pPdiff, pPKP,
c             pPKiKP, sP, sPdiff, sPKP, and sPKiKP
c      S+     gives S-up, S, Sdiff, SKS, sS, sSdiff, sSKS, pS, pSdiff,
c             and pSKS
c      basic  gives P+ and S+ as well as ScP, SKP, PKKP, SKKP, PP, and
c             P'P'
c   Note that generic S gives S-up, Sdiff, and SKS already and so
c   doesn't require a keyword.
c
      data cmdcd/'P','P+','basic','S+'/
      data cmdlst/'P','PKiKP','PcP','pP','pPKiKP','sP','sPKiKP','ScP',
     1 'SKP','PKKP','SKKP','PP','S','ScS','sS','pS'/
      data ncmpt/1,2,1,7,1,13,13,16/
c
c   Take care of the print flags.
      prnt(1)=prflg(1)
      prnt(2)=prflg(2)
      if(prnt(1)) prnt(2)=.true.
c   Copy the token list into local storage.
      no=min0(nn,jseg)
      do 23 i=1,no
 23   phlst(i)=pcntl(i)
c   See if we are in query mode.
      if(no.gt.1.or.(phlst(1).ne.'query'.and.phlst(1).ne.'QUERY'))
     1 go to 1
c
c   In query mode, get the tokens interactively into local storage.
c
 22   print *,'Enter desired branch control list at the prompts:'
      no=0
 21   call query(' ',fnd)
      if(no.ge.jseg) go to 1
      no=no+1
      read 100,phlst(no)
 100  format(a)
c   Terminate the list of tokens with a blank entry.
      if(phlst(no).ne.' ') go to 21
      no=no-1
      if(no.gt.0) go to 1
c   If the first token is blank, help the user out.
      print *,'You must enter some branch control information!'
      print *,'     possibilities are:'
      print *,'          all'
      print 101,cmdcd
 101  format(11x,a)
      print *,'          or any generic phase name'
      go to 22
c
c   An 'all' keyword is easy as this is already the default.
 1    all=.false.
      if(no.eq.1.and.(phlst(1).eq.'all'.or.phlst(1).eq.'ALL'))
     1 all=.true.
      if(all.and..not.prflg(3)) return
c
c   Make one or two generic branch names for each segment.  For example,
c   the P segment will have the names P and PKP, the PcP segment will
c   have the name PcP, etc.
c
      kseg=0
      j=0
c   Loop over the segments.
      do 2 i=1,nseg
      if(.not.all) segmsk(i)=.false.
c   For each segment, loop over associated branches.
 9    j=j+1
      phtmp=phcd(j)
c   Turn the specific branch name into a generic name by stripping out
c   the crustal branch and core phase branch identifiers.
      do 3 l=2,8
 6    if(phtmp(l:l).eq.' ') go to 4
      if(phtmp(l:l).ne.'g'.and.phtmp(l:l).ne.'b'.and.phtmp(l:l).ne.'n')
     1 go to 5
      if(l.lt.8) phtmp(l:)=phtmp(l+1:)
      if(l.ge.8) phtmp(l:)=' '
      go to 6
 5    if(l.ge.8) go to 3
      if(phtmp(l:l+1).ne.'ab'.and.phtmp(l:l+1).ne.'ac'.and.
     1 phtmp(l:l+1).ne.'df') go to 3
      phtmp(l:)=' '
      go to 4
 3    continue
c4    print *,'j phcd phtmp =',j,' ',phcd(j),' ',phtmp
c
c   Make sure generic names are unique within a segment.
 4    if(kseg.lt.1) go to 7
      if(phtmp.eq.segcd(kseg)) go to 8
 7    kseg=kseg+1
      segcd(kseg)=phtmp
      nsgpt(kseg)=i
c     if(prflg(3)) print *,'kseg nsgpt segcd =',kseg,nsgpt(kseg),' ',
c    1 segcd(kseg)
 8    if(jidx(j).lt.indx(i,2)) go to 9
 2    continue
      if(all) go to 24
c
c   Interpret the tokens in terms of the generic branch names.
c
      do 10 i=1,no
c   Try for a keyword first.
      do 11 j=1,ncmd
      if(phlst(i).eq.cmdcd(j)) go to 12
 11   continue
c
c   If the token isn't a keyword, see if it is a generic branch name.
      fnd=.false.
      do 14 k=1,kseg
      if(phlst(i).ne.segcd(k)) go to 14
      fnd=.true.
      l=nsgpt(k)
      segmsk(l)=.true.
c     print *,'Brnset:  phase found - i k l segcd =',i,k,l,' ',
c    1 segcd(k)
 14   continue
c   If no matching entry is found, warn the caller.
      if(.not.fnd) print *,'Brnset:  phase ',phlst(i),' not found.'
      go to 10
c
c   If the token is a keyword, find the matching generic branch names.
 12   j1=ncmpt(1,j)
      j2=ncmpt(2,j)
      do 15 j=j1,j2
      do 15 k=1,kseg
      if(cmdlst(j).ne.segcd(k)) go to 15
      l=nsgpt(k)
      segmsk(l)=.true.
c     print *,'Brnset:  cmdlst found - j k l segcd =',j,k,l,' ',
c    1 segcd(k)
 15   continue
 10   continue
c
c   Make the caller a list of the generic branch names selected.
c
 24   if(.not.prflg(3)) return
      fnd=.false.
      j2=0
c   Loop over segments.
      do 16 i=1,nseg
      if(.not.segmsk(i)) go to 16
c   If selected, find the associated generic branch names.
      j2=j2+1
      do 17 j1=j2,kseg
      if(nsgpt(j1).eq.i) go to 18
 17   continue
      print *,'Brnset:  Segment pointer (',i,') missing?'
      go to 16
 18   do 19 j2=j1,kseg
      if(nsgpt(j2).ne.i) go to 20
 19   continue
      j2=kseg+1
c   Print the result.
 20   j2=j2-1
c     if(.not.fnd) print *,'Brnset:  the following phases have '//
c    1 'been selected -'
      fnd=.true.
c     print 102,i,(segcd(j),j=j1,j2)
c102  format(10x,i5,5(2x,a))
 16   continue
      return
      end
c
      block data

      include "ttlim.inc"

      logical segmsk,prnt
      common/prtflc/ segmsk(jseg),prnt(2)
      data segmsk,prnt/jseg*.true.,2*.false./ 

      double precision zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
c
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      data tauc,xc/jtsm*0d0,jxsm*0d0/

      integer ka
      real deplim
      double precision ua,taua
      common/pdec/ua(5,2),taua(5,2),deplim,ka
      data deplim,ka/1.1,4/ 

      end
logtest.f90
program logtest

  logical l1,l2,l3,l4,l5

  open(1,file='mode_set.in')
  read(1,*) l1
  read(1,*) l2  
  read(1,*) l3
  read(1,*) l4
  read(1,*) l6

  close(1)

  print *,l1,l2,l3,l4,l5

end program logtest
make_aklayers.f90
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



make_akvelP.f90
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



make_akvelPS.f90
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
! This program generates velocities specific to the ak135 reference model. 
! The 210 km S-velocity discontinuity and discontinuities inside the core are 
! not treated as explicit discontinuities. Only use in conjunction with the
! program make_aklayers that defines the position of the velocity discontuities 
! treated explicitly


implicit none 
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

real(kind=dp)  :: rmin, latmin,longmin

real(kind=dp)  :: stretch = 1.01_dp
REAL(KIND=dp)      :: earth_radius = 6371.0_dp
integer        :: nv0 = 6
integer        :: n_vtypes = 2
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



make_cmplayers.f90
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



make_cmpvelocities.f90
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



make_interfaces.f90
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



make_vgrids.f90
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



matchref.f90
subroutine match_reflections(interface_id,tf1,tf2,extrema,n_extrema,s,sr)

! this subroutine finds stationary points of a travel time field on an intersection

  use mod_3dfm
  implicit none

  integer                                    :: interface_id
  type(Ttime_field)                          :: tf1,tf2
  real(kind=dp)                              :: extrema(3,*)
  integer                                    :: n_extrema
  type(Tintersection),pointer                :: isec    
  type(Tsource)                              :: s,sr
  real(kind=dp),dimension(:,:),allocatable   :: tgsum,tg1,tg2
  real(kind=dp),dimension(:),allocatable     :: tsum,tgnorm
  logical,dimension(:),allocatable           :: node_invalid,node_is_locmin,node_is_tmin
  integer                                    :: n,i,nnode,ii
  integer                                    :: ind_start,ind_end


  type(Ttriangulation)                       :: tri
  real(kind=dp)                              :: norm(3),vec(3)
  real(kind=dp)                              :: dotp
  real(kind=dp)                              :: deg_to_rad
  real(kind=dp)                              :: locmin_ratio 
  integer                                    :: iv,ivv,ip,ip1,ipp1,it,it1,itt,itt1
  

  locmin_ratio = 0.05_dp


! some validation of the input parameters

  if (tf2%reg%id /= tf1%reg%id ) stop 'regions of specified time fields unmatched in match_reflection'
  if (interface_id /= tf1%reg%itop%iface_id .and. interface_id /= tf1%reg%ibot%iface_id)  &
       stop 'illegal interface in match_reflection'


  deg_to_rad = acos(-1.0_dp)/180.0_dp


! 
! get the gradients for the fitting interface from the regional time gradient array

  isec => intersection(interface_id)
  nnode = isec%nnode
  allocate(tgsum(3,nnode),tg1(3,nnode),tg2(3,nnode),node_invalid(nnode),tsum(nnode),tgnorm(nnode),&
       node_is_locmin(nnode),node_is_tmin(nnode))

  if (interface_id == tf1%reg%itop%iface_id) then

     ind_start = tf1%reg%ngnode + 1
     ind_end = ind_start + nnode - 1

  else

     ind_end = tf1%reg%nnode
     ind_start = ind_end - nnode + 1

  endif

  tg1=tf1%time_gradient(1:3,ind_start:ind_end)
  tg2=tf2%time_gradient(1:3,ind_start:ind_end)
  tgsum=tg1+tg2
  tsum = tf1%arrivaltime(ind_start:ind_end) + tf2%arrivaltime(ind_start:ind_end)

! construct the norm of the parallel time gradient

  do i=1,nnode

!     call interface_normal(isec%lat(i),isec%long(i),intrface(interface_id),norm(1),norm(2),norm(3),h)
     norm = isec%normal(1:3,i)
     dotp=dot_product(norm,tgsum(1:3,i))
     vec(1:3) = tgsum(1:3,i) - dotp*norm(1:3)
     tgnorm(i)=sqrt(sum(vec**2))

  end do


! set nodes too close to the source or receiver as invalid since no valid time gradient may be available

  do n=1,nnode

     node_invalid(n) = (sqrt( (s%r-isec%r(n))**2 + (isec%r(n)*(s%lat-isec%lat(n)))**2 + &
          (isec%r(n)*(s%long-isec%long(n)))**2) < isec%r(n)*pgrid%dlat0) .or. &
          (sqrt( (sr%r-isec%r(n))**2 + (isec%r(n)*(sr%lat-isec%lat(n)))**2 + &
          (isec%r(n)*(sr%long-isec%long(n)))**2) < isec%r(n)*pgrid%dlat0) 

  end do


! now do a Delaunay triangulation of the intersection nodes

  call triangulation_defaults(tri)
  tri%npoints=nnode
  allocate(tri%points(2,nnode))
  tri%points(1,:) = isec%lat
  tri%points(2,:) = isec%long

  call triangulate(tri)


! find the modes that are local minima of tgnorm

  node_is_locmin=.true.

tgminloop:  do ip=1,nnode    ! all nodes

     do it=1,tri%n_triangles_from_point(ip)   ! triangles connected to zero level node

        it1=tri%triangles_from_point(it,ip)

     ! don't consider the node of a conncted triangle if on the edge (may be artifcial minimum)

        if (count(tri%triangle_neighbours(1:3,it1) == 0) /= 0 ) then
           node_is_locmin(ip) = .false.
           cycle tgminloop
        endif

        do iv=1,3    ! first level nodes of triangles connected to zero level node

           ip1 = tri%points_from_triangle(iv,it1)

           if (tgnorm(ip1) < tgnorm(ip)) then
              node_is_locmin(ip) = .false.
              cycle tgminloop
           endif

           do itt=1,tri%n_triangles_from_point(ip1) ! triangles connected to first level nodes

              itt1=tri%triangles_from_point(itt,ip1)

              do ivv=1,3    ! second level nodes of triangles connected to first level nodes

                 ipp1 = tri%points_from_triangle(ivv,itt1)

                 if (tgnorm(ipp1) < tgnorm(ip)) then
                    node_is_locmin(ip) = .false.
                    cycle tgminloop
                 endif

              end do

           end do

        end do

     end do

  end do  tgminloop



! find the modes that are local minima of tsum

  node_is_tmin=.true.

minloop:  do ip=1,nnode    ! all nodes

     do it=1,tri%n_triangles_from_point(ip)   ! triangles connected to zero level node

        it1=tri%triangles_from_point(it,ip)

     ! don't consider the node of a conncted triangle if on the edge (may be artifcial minimum)

        if (count(tri%triangle_neighbours(1:3,it1) == 0) /= 0 ) then
           node_is_locmin(ip) = .false.
           cycle minloop
        endif

        do iv=1,3    ! first level nodes of triangles connected to zero level node

           ip1 = tri%points_from_triangle(iv,it1)

           if (tsum(ip1) < tsum(ip)) then
              node_is_tmin(ip) = .false.
              cycle minloop
           endif

           do itt=1,tri%n_triangles_from_point(ip1) ! triangles connected to first level nodes

              itt1=tri%triangles_from_point(itt,ip1)

              do ivv=1,3    ! second level nodes of triangles connected to first level nodes

                 ipp1 = tri%points_from_triangle(ivv,itt1)

                 if (tsum(ipp1) < tsum(ip)) then
                    node_is_tmin(ip) = .false.
                    cycle minloop
                 endif

              end do

           end do

        end do

     end do

  end do  minloop

! select further to take only nodes with a sufficiently small value at the minimum and that are not
! right next to the source or receiver

do n=1,nnode
  node_is_locmin(n) = node_is_locmin(n) .and. (tgnorm(n) < locmin_ratio*(2.0_dp*sqrt(sum(tg1(1:3,n)**2)))) &
       .and. .not.node_invalid(n)
end do

do n=1,nnode
  node_is_tmin(n) = node_is_tmin(n) .and. .not.node_invalid(n)
end do

! copy the reflection points to the output array

   n_extrema=count(node_is_locmin)

   if (n_extrema > 0) then

      ii=0
      do n=1,nnode
         if(node_is_locmin(n)) then
            ii=ii+1
            extrema(1,ii)= isec%r(n)
            extrema(2,ii)= isec%lat(n)
            extrema(3,ii)= isec%long(n)

         endif
      end do

   else

      return

   endif

  deallocate(tgsum,tg1,tg2,node_invalid,tsum,tgnorm,node_is_locmin,node_is_tmin)
  call clean_triangulation(tri)

  return

end subroutine match_reflections

!***************************************************************************************************
subroutine trace_reflectionfit(n,m)

! this subroutine generates the ray paths for a phase that includes a late reflection search

  use mod_3dfm

  integer   :: i,k,m,n,ifit,nsec,npsec,iref,recpath_id,path_id,nsec_comp

  type(Tsource),pointer                  :: s,ss
  type(Treceiver)                        :: receiver_at_reflection
  type(Tray),pointer                     :: ray,subray
  type(Tray)                             :: ray_to_src,ray_to_rec
  type(Tray_section),pointer             :: raysec,subsec
  type(Ttime_field),pointer              :: tf1,tf2
  real(kind=dp)                          :: reflections(3,10),deg_to_rad,t_arrival
  integer                                :: n_reflections

  deg_to_rad = acos(-1.0_dp)/180.0_dp

  ray => receiver(n)%ray(m)
  s   => ray%source
  path_id = ray%raypath_id

  ss   => source(receiver(n)%source_equivalent)
  recpath_id = receiver(n)%path_equivalent(m)
  tf1  => s%time_field(s%path(path_id)%tf_sequence(s%path(path_id)%n_tf))
  tf2  => ss%time_field(ss%path(recpath_id)%tf_sequence(ss%path(recpath_id)%n_tf))
  ifit =  s%path(path_id)%fitting_interface

!  print *,'calling match-reflections with'
!  print *,'source 1',s%id
!  print *,'path 1',path_id
!  print *,'tf 1',tf1%id
!  print *,'source 2',ss%id
!  print *,'path 2',recpath_id
!  print *,'tf 2',tf2%id

  call ray_defaults(ray_to_rec)
  call ray_defaults(ray_to_src)

  call match_reflections(ifit,tf1,tf2,reflections,n_reflections,s,ss)

  if (n_reflections > 0) print '(a5,i3,a12,i3,a3,i5,a20)', &
       'ray',m,'from source',s%id, ':',n_reflections,'reflection points'

  if (n_reflections == 0) then
     print '(a12,i4,a10,i4,a15,i4,a23)','ray',m,'to source',ray%source_id,&
          ' from receiver',n,' : no reflections found' 
     ray%valid = .false.

     k=0
     t_arrival=-1.0_dp               
     write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,k,t_arrival,ray%diffracted,ray%headwave

     return   ! go to next ray
  endif


! if there are reflections, the array of subrays within the ray associated with the path
! is allocated. The subrays will contain raypaths from all the reflection points found.

  if (n_reflections > 0) then
     ray%is_multiray = .true.
     ray%n_subrays = n_reflections
     allocate(ray%subray(n_reflections))
     do i=1,n_reflections ; call ray_defaults(ray%subray(i)) ; end do
  endif

!  define the path of the two rays from reflection point to source and receiver

  ray_to_src%source     => s
  ray_to_src%raypath_id = path_id
  ray_to_rec%source     => ss
  ray_to_rec%raypath_id = recpath_id


  do iref=1,n_reflections

     subray => ray%subray(iref)

!  create a receiver at the reflection point

     receiver_at_reflection%r    = reflections(1,iref)
     receiver_at_reflection%lat  = reflections(2,iref)
     receiver_at_reflection%long = reflections(3,iref)

     ray_to_src%valid = .true.
     ray_to_rec%valid = .true.

!  trace ray from reflection point to source


     call trace_ray_from_receiver(receiver_at_reflection,s,ray_to_src)
     if (.not.ray_to_src%valid) then   ! no valid ray path found
        print *,'subray from reflection point',iref,' to source was not valid'
        ray%subray(iref)%valid = .false.

        k=0
        t_arrival=-1.0_dp               
        write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,iref,t_arrival,subray%diffracted,subray%headwave

        cycle
     endif
   
!  trace ray from reflection point to receiver

     call trace_ray_from_receiver(receiver_at_reflection,ss,ray_to_rec)
     if (.not.ray_to_rec%valid) then   ! no valid ray path found
        print *,'subray from reflection point',iref,' to receiver was not valid'
        ray%subray(iref)%valid = .false.

        k=0
        t_arrival=-1.0_dp               
        write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,iref,t_arrival,subray%diffracted,subray%headwave

        cycle
     endif

  
! assemble total ray from the two pieces

  ! first assign the attributes of the composite ray

     subray%valid = .true.
     subray%source => s
     subray%raypath_id = path_id
     subray%nsections=ray_to_src%nsections+ray_to_rec%nsections
     subray%source_id = s%id
     subray%receiver_time = ray_to_src%receiver_time + ray_to_rec%receiver_time 
     subray%receiver_time_gradient(1:3) = 0.0_dp
     subray%diffracted = ray_to_src%diffracted .or. ray_to_rec%diffracted 
     subray%headwave = ray_to_src%headwave .or. ray_to_rec%headwave

     allocate(subray%section(subray%nsections))
     do i=1,subray%nsections ; call ray_section_defaults(subray%section(i)) ; end do

  ! transfer the ray sections from the reflection point to receiver part
  ! here the direction of travel has to be inverted

     do nsec=ray_to_rec%nsections,1,-1

        nsec_comp = ray_to_src%nsections + ray_to_rec%nsections - nsec + 1
        npsec = ray_to_rec%section(nsec)%npoints

        allocate(subray%section(nsec_comp)%point(3,npsec))

        raysec => ray_to_rec%section(nsec)
        subsec => subray%section(nsec_comp)

        subsec%ray => subray
        subsec%reg => raysec%reg 
        subsec%istart => raysec%iend
        subsec%iend => raysec%istart
        subsec%tf => raysec%tf
        subsec%source => raysec%source
        subsec%npoints  =  raysec%npoints      
        subsec%place_in_sequence  = nsec_comp                   
        subsec%diffracted = raysec%diffracted
        subsec%headwave = raysec%headwave

        do i=1,subsec%npoints
           subsec%point(1:3,i) = raysec%point( 1:3 , raysec%npoints - i + 1)
        end do

        deallocate(raysec%point)

     end do

  ! transfer the ray sections  from the reflection point to source part
  ! here the direction of travel is the same 

     do nsec=ray_to_src%nsections,1,-1

        nsec_comp = nsec 
        npsec = ray_to_src%section(nsec)%npoints

        allocate(subray%section(nsec_comp)%point(3,npsec))

        raysec => ray_to_src%section(nsec)
        subsec => subray%section(nsec_comp)

        subsec%ray => subray
        subsec%reg => raysec%reg 
        subsec%istart => raysec%istart
        subsec%iend => raysec%iend
        subsec%tf => raysec%tf
        subsec%source => raysec%source
        subsec%npoints  =  raysec%npoints      
        subsec%place_in_sequence  = nsec_comp                 
        subsec%diffracted = raysec%diffracted
        subsec%headwave = raysec%headwave

        do i=1,subsec%npoints
           subsec%point(1:3,i) = raysec%point(1:3,i)
        end do

        deallocate(raysec%point)

     end do

   ! have to set the starting interface of the ray section from the reflection point to the source
   ! to the reflection interface 

     subray%section(ray_to_src%nsections)%istart => intersection(ifit)


     deallocate(ray_to_src%section,ray_to_rec%section)

     print '(a12,i4,a10,i4,a15,i4,a8,i3,a4,f10.5,2l3,3f10.4)','traced ray',m,'to source',ray%source_id, &
          ' from receiver',n,'subray ',iref,'  t=', subray%receiver_time, subray%diffracted,subray%headwave, &
          reflections(1,iref), reflections(2:3,iref)/deg_to_rad

     write(11,'(4i6,f15.6,2l5)') n,ray%source_id,m,iref,subray%receiver_time,subray%diffracted,subray%headwave

     if (display_mode) call store_ray(subray)

  end do  ! loop over reflection points

  return

end subroutine trace_reflectionfit

!***************************************************************************************************
function interpolate_triangle(pos,field,vpos)

! linear interpolation in a triangle, the weight of each vertex is proportional to
! the distance from the point along the normal of the opposite side

  use mod_3dfm

  real(kind=dp)  :: pos(2),field(3),vpos(2,3),w(3),x(2),e(2),d(2),h
  real(kind=dp)  :: interpolate_triangle


  x = pos - vpos(1:2,2)
  e = vpos(1:2,3) - vpos(1:2,2)
  e = e/sqrt(sum(e**2))
  d = vpos(1:2,1) - vpos(1:2,2)
  h = sqrt(sum((d - dot_product(d,e)*e)**2))
  w(1)= sqrt(sum((x - dot_product(x,e)*e)**2))/h

  x = pos - vpos(1:2,3)
  e = vpos(1:2,1) - vpos(1:2,3)
  e = e/sqrt(sum(e**2))
  d = vpos(1:2,2) - vpos(1:2,3)
  h = sqrt(sum((d - dot_product(d,e)*e)**2))
  w(2)= sqrt(sum((x - dot_product(x,e)*e)**2))/h


  x = pos - vpos(1:2,1)
  e = vpos(1:2,2) - vpos(1:2,1)
  e = e/sqrt(sum(e**2))
  d = vpos(1:2,3) - vpos(1:2,1)
  h = sqrt(sum((d - dot_product(d,e)*e)**2))
  w(3)= sqrt(sum((x - dot_product(x,e)*e)**2))/h

  w=w/sum(w(1:3))

  h = dot_product(w,field)

  interpolate_triangle=dot_product(w,field)

  return

end function interpolate_triangle
mod_3dfm.f90
module type_definitions

  INTEGER, PARAMETER     :: sp = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

! declarations of derived types
!------------------------------------------------------

type Tpointer_to_integer_array

   integer,dimension(:),pointer   :: p 
   
end type Tpointer_to_integer_array


type Tpropagation_grid    ! a regular grid used for wave front propagation with fast marching

    integer       :: nr,nlong,nlat     ! # of grid cells in each direction
    REAL(KIND=dp) :: dr0,dlat0,dlong0  ! grid step sizes
    REAL(KIND=dp) :: r0,lat0,long0     ! position of grid origin
    REAL(KIND=dp) :: rmax,latmax,longmax     ! position of grid corner opposite origin
    REAL(KIND=dp) :: tolerance         ! tolerance for deciding if a position coincides exactly with the grid or not
    logical       :: is_main_grid 
    logical       :: is_source_grid 

    integer       :: nnode             ! total # of nodes

    integer       :: index_r0,index_lat0,index_long0  ! indices of the origin of a refined source grid in the main grid



    REAL(KIND=dp), DIMENSION (:), pointer :: r 
    REAL(KIND=dp), DIMENSION (:), pointer :: lat 
    REAL(KIND=dp), DIMENSION (:), pointer :: long 
    REAL(KIND=dp), DIMENSION (:), pointer :: coslat 


    REAL(KIND=dp), DIMENSION (:,:,:,:), pointer     :: velocity 
    REAL(KIND=dp), DIMENSION (:,:,:), pointer     :: arrivaltime 
    REAL(KIND=dp), DIMENSION (:,:,:,:), pointer   :: time_gradient 

    integer, DIMENSION (:,:,:), pointer :: node_region           ! identifies region of which this node is part
    integer, DIMENSION (:,:,:), pointer :: rnode_id              ! index of this node in the regional node list of
                                                                          ! the region to which it belongs

    logical, dimension(:,:,:),  pointer :: fully_regular       ! is node surrounded by regular cells only


    ! Here it gets complicated: ccind_from_3dc stands for cut cell index from 3D coordinates.
    ! If a cell i,j,k of the grid is cut by an interface, and therefore has intersection nodes,
    ! the pointer ccind_from_3dc(i,j,k) is allocated, and points to an integer array ccind_from_3dc(i,j,k)%p
    ! of length n_intersections. If the cell is cut by intersection n, ccind_from_3dc(i,j,k)%p(n) contains
    ! the index of the cut cell i,j,k in the cut cell list of intersection n. For an intersection m that does not
    ! cut the cell  ccind_from_3dc(i,j,k)%p(m) is zero. In this way we can identify the intersection nodes that
    ! are part of a cell from its position in the grid.
    ! 

    type(Tpointer_to_integer_array), DIMENSION (:,:,:), pointer :: ccind_from_3dc 

end type Tpropagation_grid

!-----------------------------------------------------------------------------------

type Tvelocity_grid   ! regular grid defining the velocity as a function of position
                      ! Note: usually each velocity grid is defined over the entire propagation
                      ! grid, but only some nodes actually influence the corresponding region  

    integer       :: nr,nlong,nlat     ! # of grid cells in each direction
    REAL(KIND=dp) :: dr0,dlat0,dlong0  ! grid step sizes
    REAL(KIND=dp) :: r0,lat0,long0     ! position of grid origin


    integer       :: nnode             ! total # of nodes

    integer       :: start_index       ! for use in inversion if the grid is a velocity grid
    logical       :: to_be_inverted ! for use in inversion if the grid is a velocity grid


    REAL(KIND=dp), DIMENSION (:), pointer :: r 
    REAL(KIND=dp), DIMENSION (:), pointer :: lat 
    REAL(KIND=dp), DIMENSION (:), pointer :: long 


    REAL(KIND=dp), DIMENSION (:,:,:), pointer     :: velocity 
    logical,dimension(:,:,:),pointer              :: active   ! set to true for the nodes that actually
                                                              ! influence the region to which the grid belongs 


end type Tvelocity_grid


!-------------------------------

type Tpath   ! contains information about the sequence of interface interactions defining a path in
             ! multi-stage fast marching

    integer                            :: id                ! index of path in user-defined list
    integer                            :: n_tf              ! # of time fields on each path
    integer,dimension(:) ,pointer      :: sequence    ! the sequence of interfaces on each path
    integer,dimension(:) ,pointer      :: tf_sequence ! the sequence of time fields on each path
                                                               ! number refers to index in time_field array in source
    integer,dimension(:) ,pointer      :: vtype_sequence    ! the sequence of velocity types on each path

    logical                            :: valid       ! flag for valid path, required sequence may not exist
    logical                            :: used         ! flag set if path is actually used
    logical                            :: gridsave     ! flag set if a grid of arrival times is to be saved
                                                       ! for this path
    integer                            :: first_tf_to_save


    integer                            :: refstep          ! =0 if no late reflections, step indicating refl otherwise
    integer                            :: fitting_interface  ! interface at which fitting is performed

end type Tpath

!-------------------------------------------------------------------------------------------------------


type Tbackpointer    ! three integer numbers that uniquely identify each node
                     ! for regular grid nodes these are the grid coordinates in r,lat and long
                     ! for interface nodes i1=0, i2= id of interface, i3 = number of node in list of interface
                     ! nodes

    integer :: i1
    integer :: i2
    integer :: i3

end type Tbackpointer

!-------------------------------
type Tsource

    integer                  :: id                  ! id # of source
    REAL(KIND=dp)            :: r,lat,long,coslat   ! position
    integer                  :: ir,ilat,ilong       ! main grid cell containing source
    type(Tbackpointer)       :: cnode(100)          ! nodes connected to the source
    integer                  :: n_cnode             ! # of nodes connected to the source

    logical                  :: on_grid,on_interface,on_pinched_interface
    integer                  :: region_id,interface_id
    integer                  :: topreg_id,botreg_id,topint_id,botint_id

    logical                  :: is_local,is_teleseismic
    integer                  :: teleseismic_id
    character(len=8)         :: teleseismic_phase

    integer                  :: nfile

    integer                  :: n_tf_init        ! # of timefields in which source lies (1/2)
                                                    ! 1 if in a region, 2 if on an interface that is not top or bottom

    integer,dimension(:),pointer   :: first_tf_up   ! indices of source time fields
    integer,dimension(:),pointer   :: first_tf_down ! indices of source time fields


! these are the paths (sequence of reflections/refractions) originating from this source that have to be computed
    integer                              :: n_paths       ! # of paths for this source
    type(Tpath),dimension(:),pointer     :: path     ! the actual paths
    type(Tpath),dimension(2,2)           :: init_path

! these are the time fields containing regional traveltimes for all computed paths
    integer                                 :: n_time_fields      ! # of time fields for this source
    type(Ttime_field),dimension(:),pointer  :: time_field    ! list of time fields


! parameters having to do with the calculation of Frechet derivatives
    integer                  :: start_index         ! location of source parametrs in inversion parameter array
    logical                  :: to_be_inverted      ! flag showing whether the the source parameters are fixed or free

     
end type Tsource


!-------------------------------!

type Treceiver

    integer                                    :: id                    ! identifies the receiver
    REAL(KIND=dp)                              :: r,lat,long,coslat     ! position
    REAL(KIND=dp), dimension(:),pointer        :: arrivaltime           ! the arrival time at the receiver (not used)
    integer                                    :: n_rays                ! number of paths ending at this receiver
    type(Tray),dimension(:),pointer            :: ray                   ! rays of the paths ending at this receiver


  ! the variables below are only used if the path contains a reflection fit
  ! they identify the sequence of time fields from the receiver to the intermediate fitting interface

    integer                                    :: source_equivalent  ! index in source list of this receiver
    integer,dimension(:),pointer               :: path_equivalent       ! corresponding path in the source_equivalent
                                                                        ! path list

    
end type Treceiver

!-------------------------------!

type Tinterface
! this type defines the position of an interface
!  note: type Tinterface defines of the position of the interface (used as input by the cubic spline interpolation)
!        type Tintersection contains the actual nodes of the interface and everything related to them

! grid definitions
    integer :: nlat,nlong,id                        ! # of points in lat,long
    REAL(KIND=dp) :: dlat0,dlong0                   ! size of intervals
    REAL(KIND=dp) :: lat0,long0                     ! position of grid origin
    logical       :: pinched                        ! true if the interface touches another


! parametrs of the inversion
    integer       :: nnode                          ! # of parameters describing the interface position    
    integer       :: start_index                    ! start index of the interface parameters in the global parametr list
    logical       :: to_be_inverted        ! is the position of this interface is to be inverted for?


    REAL(KIND=dp), DIMENSION (:), pointer :: lat      
    REAL(KIND=dp), DIMENSION (:), pointer :: long 

    REAL(KIND=dp), DIMENSION (:,:), pointer :: r    ! the actual radius values for the interface at the nodes

end type Tinterface

!-------------------------------!

type Tinteger_coordinates

    integer :: ir,ilat,ilong

end type Tinteger_coordinates

!-------------------------------

type Tintersection
! this type contains information about where an interface intersects with a grid
!  note: type Tinterface defines of the position of the interface (used as input by the cubic spline interpolation)
!        type Tintersection contains the actual nodes of the interface and everything related to them

   integer :: nnode,n_ccells                         ! # of intersection nodes, # of grid cells cut by the interface
   integer :: id                                 ! intersection id
   integer :: iface_id                            ! corresponding interface id
   type(Tpropagation_grid),pointer :: grid  ! pointer to the grid on which the intersection is defined        

! flags
   logical :: pinched           ! true if the the intersection touches another

   REAL(KIND=dp), DIMENSION (:), pointer              :: r             ! positions of the intersection nodes
   REAL(KIND=dp), DIMENSION (:), pointer              :: lat 
   REAL(KIND=dp), DIMENSION (:), pointer              :: long 
   REAL(KIND=dp), DIMENSION (:), pointer              :: coslat 
   REAL(KIND=dp), DIMENSION (:,:), pointer            :: normal        ! vector giving the normal of the interface
                                                                                ! at every intersection node

   REAL(KIND=dp), DIMENSION (:,:), pointer              :: vel_top       ! velocities on both side of the interface
   REAL(KIND=dp), DIMENSION (:,:), pointer              :: vel_bot 

   REAL(KIND=dp), DIMENSION (:), pointer              :: arrivaltime   ! time a front arrives at an intersection node
   REAL(KIND=dp), DIMENSION (:), pointer              :: starttime     ! time a front starts from an intersection node
   REAL(KIND=dp), DIMENSION (:,:), pointer            :: time_gradient ! time gradient at an interface node
                                                 ! the time gradient is normal to the wave front and has size 1/velocity

   integer,dimension(:),pointer           :: intype      ! type of inode (0:on grid,1,2,3 on r,lat long connection) 


! the variables below allow us to identify the grid cells that an intersection node is part of
! and the inverse, the intersection nodes that are part of a given grid cell

   type(Tinteger_coordinates), DIMENSION (:), pointer :: ccells           ! list of cells cut by this interface
   integer,dimension(:), pointer                      :: n_inodes         ! number of inodes in cut cell
   integer,dimension(:,:),pointer                     :: inodes           ! actual inodes in cut cell
   integer, DIMENSION (:,:), pointer                  :: ccell_from_inode ! pointer from inode to cut cells it is in


! these variables connect the intersection with the regions (abov and below) that is is part of

   type(Tregion),pointer                                  :: regabo       ! the region above the intersection
   type(Tregion),pointer                                  :: regbel      ! the region above the intersection

   integer, DIMENSION (:), pointer          :: rabo_node_id   !index of each intersection node in the node 
                                                                       ! list of the region above the intersection

   integer, DIMENSION (:), pointer          :: rbel_node_id   !index of each intersection node in the node 
                                                                       ! list of the region below the intersection

   integer, DIMENSION (:,:), pointer          :: irg_abo   !first regular grid node above the interface at (j,k)
   integer, DIMENSION (:,:), pointer          :: irg_bel  !first regular grid node below the interface at (j,k)


end type Tintersection

!-------------------------------

type Tregion
! this type contains a region of the propagation grid between interfaces
! the region is a volume of space between two interfaces, and contains nodes of the regular grid 
! and both bounding intersections. Fast marching proceeds on a region by region basis

   integer      :: id                           ! region identification
   integer      :: ivgrid                       ! velocity definition grid associated with this region
   type(Tintersection),pointer :: itop,ibot        ! intersections at the top and bottom of this region
   type(Tpropagation_grid),pointer :: grid  ! pointer to the grid on which the region is defined   
   integer      :: ngnode                         ! number of propagation grid nodes in this region


! arrays below define a 1-D list of regular + intersection nodes constituting this region
! this 1-D array is used for the fast marching 

   integer                                     :: nnode   ! total number of nodes in this region including boundary nodes
   type(Tbackpointer),dimension(:),pointer     :: node    ! points back from 1-D list of nodes in the 
                                                                   ! region to gridnode or intnode it corresponds with
   integer, dimension(:), pointer              :: node_status  ! fast marching status
   real(kind=dp),dimension(:),pointer          :: arrivaltime  ! 
   real(kind=dp),dimension(:,:),pointer        :: time_gradient  !
   real(kind=dp),dimension(:,:),pointer        :: velocity     !
   real(kind=dp),dimension(:),pointer          :: r    !
   real(kind=dp),dimension(:),pointer          :: lat   !
   real(kind=dp),dimension(:),pointer          :: long     !
   REAL(KIND=dp),DIMENSION(:),pointer          :: coslat 

! nodes that have been initialized from a teleseismic source
   integer                                     :: n_init       ! number of initialized nodes
   integer,dimension(:),pointer                :: init_id      ! the list of initialized nodes
   integer,dimension(:),pointer                :: init_type      ! the list of initialized nodes
   real(kind=dp),dimension(:),pointer          :: init_arrivaltime  ! arrival time of the init nodes
   real(kind=dp),dimension(:,:),pointer        :: init_time_gradient ! incoming gradient at the init node

end type Tregion

!-------------------------------

type Ttime_field
! the type Ttimefield contains a fast marching solution for the arrival times
! and the time gradients in a region


   integer                                    :: id                             ! identifies the time field
   integer                                    :: vtype                 ! the type of velocity used for this tf
   integer                                    :: nfile          ! number of file in which data will be stored   
   type(Tregion), pointer                     :: reg                   ! the region this time field refers to
   type(Tintersection),pointer                :: istart                ! the intersection from which the FMM started 
   type(Tintersection),pointer                :: inonstart            ! the other intersection of the region
   real(kind=dp),dimension(:),pointer         :: arrivaltime  !  
   real(kind=dp),dimension(:,:),pointer       :: time_gradient  !  

   logical                                    :: turning_rays_present  ! indicates if turning rays hit the starting
                                                                                ! intersection during FMM
   logical,dimension(:),pointer               :: received_turning_ray  ! if turning rays were present, this array
                                                       ! indicates which nodes of the starting intersection received them


 ! these pointers allow the time fields to be organised in a tree structure that allows for easy re-use
 ! of timefields already calculated by other path sequences

   integer                                    :: prev_tf      ! pointer to the previous time field on the path/ray
   integer                                    :: next_tf(8)   ! pointer to the next time fields, possibilities:
                                                              ! (1) up from inonstart (vtype = 1)
                                                              ! (2) down from inonstart 1
                                                              ! (3) up from istart (only possible if turning rays exist) 1
                                                              ! (4) down from istart (only possible if turning rays exist) 1
                                                              ! (1) up from inonstart (vtype = 2)
                                                              ! (2) down from inonstart 2
                                                              ! (3) up from istart (only possible if turning rays exist) 2
                                                              ! (4) down from istart (only possible if turning rays exist) 2
end type Ttime_field

!-------------------------------

type Tray_section
! a ray section is the part of a ray between two interfaces, except the first one which starts at source, and the
! last one which ends at a receiver

   type(Tray),pointer                         :: ray      ! the ray that this section is part of
   type(Tregion), pointer                     :: reg      ! the region in which the ray section lies 
   type(Tintersection),pointer                :: istart   ! the intersection at which the integration to find
                                                          ! the ray starts. Note that the rays are found integrating
                                                          ! backward in time!
   type(Tintersection),pointer                :: iend     ! the intersection at which the ray integration ends
   type(Ttime_field),pointer                  :: tf       ! the time field used for finding the ray
   type(Tsource),pointer                      :: source   ! the source of the ray

   integer                                    :: npoints           ! # of points on this ray section
   real(kind=dp),dimension(:,:),pointer       :: point  ! ! the actual positions of the points on the ray section

   integer                                    :: place_in_sequence ! the position of the section in collection of ray 
                                                                   ! sections defining the ray
   logical                                    :: diffracted 
   logical                                    :: headwave 

end type Tray_section

!-------------------------------

type Tray
! a ray is a collection of ray sections that define a given ray/path/phase

   integer                                      :: nsections      ! # of sections on the ray
   type(Tray_section),dimension(:),pointer      :: section   ! the actual ray sections
   type(Tsource),pointer                        :: source    ! the source of this ray
   integer                                      :: raypath_id      ! number of ray in list of rays from this source
   integer                                      :: source_id 
   real(kind=dp)                                :: receiver_time      ! time the ray arrives at the receiver
   real(kind=dp),dimension(3)                   :: receiver_time_gradient ! time gradient (direction) of ray at the receiver
   logical                                      :: diffracted 
   logical                                      :: headwave 
   logical                                      :: is_multiray  ! true if reflection fit is required
   integer                                      :: n_subrays          ! # of reflections found in fit
   type(Tray),dimension(:),pointer              :: subray              ! contains the rays in case of a the reflection fit 

! variables relating to the inversion process
   integer                                      :: n_pdev             ! # of non-zero partial derivatives based on this ray
   integer,dimension(:),pointer                 :: pdev_indx          ! list of inversion parameters for which the 
                                                                      ! partial derivative of the arrival time is non-zero
   real(kind=dp),dimension(:),pointer           :: pdev               ! partial derivative of the arrival time with respect
                                                                      ! to the corresponding inversion parameter in the 
                                                                      ! array pdev_indx

   logical                                      :: valid     ! set to false if this ray does not exist             

end type Tray


type Ttriangulation
! this type contains a 2-dimensional triangulation of a point set

   integer                                         :: npoints    ! number of points
   real(kind=dp),dimension(:,:),pointer            :: points  ! the coordinates of the points

   integer                                         :: ntriangles  ! # of triangles
   integer,dimension(:,:),pointer                  :: points_from_triangle  ! the points that are part of a given
                                                                                     ! triangle
   integer,dimension(:,:),pointer                  :: triangle_neighbours   ! the triangles that neighbour a
                                                                                     ! triangle

   integer,dimension(:),pointer                    :: n_triangles_from_point  ! # of triangles connected to a
                                                                                       ! given point
   integer,dimension(:,:),pointer                  :: triangles_from_point    ! indices of triangles connected to a
                                                                                       ! given point

end type Ttriangulation


type Tgrid_identifier

   integer   :: igrid,vtype

end type Tgrid_identifier


!**************************************************************************************************************************

contains

! the subroutines below can be called to give variables inside instances of the derived types defined
! above default values when allocated. Initialization inside the derived type definition
! is a Fortran 95 feature, and we had to remove it to ensure fortran 90 compatibility.


  subroutine pgrid_defaults(grid)

    type(Tpropagation_grid)  :: grid

    grid%is_main_grid = .false.
    grid%is_source_grid  = .false.
    grid%nnode = 0      

    nullify(grid%r) 
    nullify(grid%lat) 
    nullify(grid%long) 
    nullify(grid%coslat) 
    nullify(grid%velocity) 
    nullify(grid%arrivaltime) 
    nullify(grid%time_gradient) 
    nullify(grid%node_region)   
    nullify(grid%rnode_id)      
    nullify(grid%fully_regular) 
    nullify(grid%ccind_from_3dc) 

  end subroutine pgrid_defaults


  subroutine vgrid_defaults(grid)

    type(Tvelocity_grid)  :: grid

    grid%to_be_inverted  = .false.
    grid%nnode = 0      
    nullify(grid%r) 
    nullify(grid%lat) 
    nullify(grid%long) 
    nullify(grid%velocity)
    nullify(grid%active)

  end subroutine vgrid_defaults


  subroutine source_defaults(source)
    type(Tsource)  :: source

    source%on_grid = .false.
    source%on_interface = .false.
    source%on_pinched_interface = .false.
    source%region_id = 0
    source%interface_id = 0
    source%topreg_id = 0
    source%botreg_id = 0
    source%topint_id = 0
    source%botint_id = 0
    source%is_local = .false.
    source%is_teleseismic = .false.
    source%teleseismic_id = 0
    source%teleseismic_phase = ' '
    source%n_tf_init = 1 
    nullify(source%first_tf_up)
    nullify(source%first_tf_down)
    source%n_paths = 0 
    source%to_be_inverted  = .false.
    nullify(source%path)  
    source%n_time_fields = 0   
    nullify(source%time_field)

  end subroutine source_defaults


  subroutine path_defaults(path)
    type(Tpath)  :: path

    path%id = 0               
    path%n_tf = 0             
    nullify(path%sequence)
    nullify(path%tf_sequence)
    nullify(path%vtype_sequence)
    path%valid = .true.  
    path%used  = .false. 
    path%gridsave  = .false.       
    path%refstep = 0    
    path%fitting_interface = 0

  end subroutine path_defaults


  subroutine receiver_defaults(rec)
    type(Treceiver) :: rec

    rec%id  = 0              
    nullify(rec%arrivaltime)   
    rec%n_rays  = 0      
    nullify(rec%ray)                 
    rec%source_equivalent = 0
    nullify(rec%path_equivalent)
                                
  end subroutine receiver_defaults


  subroutine intersection_defaults(isec)
    type(Tintersection)  :: isec
    isec%nnode = 0
    isec%n_ccells = 0 
    isec%id = 0       
    isec%iface_id = 0 
    nullify(isec%grid)
    isec%pinched = .false.  
    nullify(isec%r)
    nullify(isec%lat)
    nullify(isec%long)
    nullify(isec%coslat)
    nullify(isec%normal)
    nullify(isec%vel_top)
    nullify(isec%vel_bot)
    nullify(isec%arrivaltime)
    nullify(isec%starttime)
    nullify(isec%time_gradient)
    nullify(isec%intype)
    nullify(isec%ccells)
    nullify(isec%n_inodes)
    nullify(isec%inodes)
    nullify(isec%ccell_from_inode)
    nullify(isec%regabo)
    nullify(isec%regbel)
    nullify(isec%rabo_node_id)
    nullify(isec%rbel_node_id)
    nullify(isec%irg_abo)
    nullify(isec%irg_bel)

  end subroutine intersection_defaults


  subroutine interface_defaults(iface)
    type(Tinterface)  :: iface

    iface%pinched = .false.       
    iface%nnode = 0                 
    iface%to_be_inverted = .false. 
    nullify(iface%lat)
    nullify(iface%long)
    nullify(iface%r)

  end subroutine interface_defaults


  subroutine region_defaults(reg)
    type(Tregion)  :: reg

    reg%id = 0          
    reg%ivgrid = 0      
    nullify(reg%grid)
    reg%ngnode   = 0 
    reg%nnode  = 0 
    reg%n_init = 0
    nullify(reg%node)
    nullify(reg%node_status)
    nullify(reg%arrivaltime)
    nullify(reg%time_gradient)
    nullify(reg%velocity)
    nullify(reg%r) 
    nullify(reg%lat)
    nullify(reg%long)
    nullify(reg%coslat)
    nullify(reg%init_id)
    nullify(reg%init_arrivaltime)
    nullify(reg%init_time_gradient)

  end subroutine region_defaults


  subroutine time_field_defaults(tf)
    type(Ttime_field)  :: tf

    tf%id = 0  
    tf%vtype = 0
    nullify(tf%reg)
    nullify(tf%istart)
    nullify(tf%inonstart)
    nullify(tf%arrivaltime)
    nullify(tf%time_gradient)
    tf%turning_rays_present = .false.
    nullify(tf%received_turning_ray)
    tf%prev_tf = 0                   
    tf%next_tf(1:8) = 0              

  end subroutine time_field_defaults


  subroutine ray_defaults(ray)
    type(Tray)  :: ray

    ray%nsections = 0  
    nullify(ray%section)
    nullify(ray%source)
    ray%raypath_id = 0 
    ray%source_id = 0
    ray%diffracted = .false.
    ray%headwave = .false.
    ray%is_multiray = .false.
    ray%n_subrays = 0        
    nullify(ray%subray)
    ray%n_pdev = 0           
    nullify(ray%pdev_indx)
    nullify(ray%pdev)
    ray%valid = .true.    

  end subroutine ray_defaults


  subroutine ray_section_defaults(raysec)
    type(Tray_section) :: raysec

    nullify(raysec%ray)
    nullify(raysec%reg)
    nullify(raysec%istart)
    nullify(raysec%iend)
    nullify(raysec%tf)
    nullify(raysec%source)
    raysec%npoints = 0       
    nullify(raysec%point)
    raysec%place_in_sequence = 0
    raysec%diffracted = .false.
    raysec%headwave = .false.

  end subroutine ray_section_defaults


  subroutine triangulation_defaults(tri)
    type(Ttriangulation) :: tri

    tri%npoints = 0  
    nullify(tri%points)
    tri%ntriangles = 0
    nullify(tri%points_from_triangle)
    nullify(tri%triangle_neighbours)
    nullify(tri%n_triangles_from_point)
    nullify(tri%triangles_from_point)

  end subroutine triangulation_defaults

end module type_definitions

!*************************************************************************************************************************


module global_variables

  use type_definitions

! global parameters

  REAL(KIND=dp), PARAMETER   :: interface_tolerance = 0.005_dp !intersection nodes snap to regular grid if closer than this
                                                    ! fraction of radial grid cell size (to avoid tiny triangles)

  REAL(KIND=dp), PARAMETER   :: huge_time = 1.0e20_dp  ! default value for a time that is larger than any realistic value

  REAL(KIND=dp), PARAMETER   :: earth_radius = 6371.0_dp

  integer   :: refinement_factor          ! reduction in size for refined source grid
  integer   :: ncell_to_be_refined        ! extent of refined grid around the source in main grid cells 
                                           ! refinement parameters are overwritten by values in propgrid.in 
  integer   :: global_source_counter
  integer   :: raypoint_counter

  logical   :: file_mode
  logical   :: no_pp_mode
  logical   :: parallel_mode
  logical   :: display_mode
  logical   :: save_rays_mode
  logical   :: save_timefields_mode

! below the definition of variables that will be accessible to every subprogram

!.........................................................................
! this array of type Tinterface defines  the location of the interfaces on the propagation grids
! by B-spline interpolation 
  integer                                                  :: n_interfaces
  type(Tinterface),dimension(:),pointer                    :: intrface 


!.........................................................................
! these are the velocity grids used to define the velocity on the propagation grids
! by B-spline interpolation
  integer                                                  :: n_vgrids
  integer                                                  :: n_vtypes
  type(Tvelocity_grid),dimension(:,:),pointer          :: vgrid 



!.........................................................................
! the main propagation grid
  type(Tpropagation_grid),pointer                         :: pgrid 

! these are the intersections associated with the main propgation grid
  integer                                                  :: n_intersections
  type(Tintersection),dimension(:),pointer                 :: intersection 

! these are the regions of the main propagation grid
  integer                                                  :: n_regions
  type(Tregion),dimension(:),pointer                       :: region 


!.........................................................................
! the fine grid around the source
  type(Tpropagation_grid),pointer                         :: sgrid 

! these are the intersections associated with the refined source grid
  integer                                                  :: n_sintersections
  type(Tintersection),dimension(:),pointer                 :: sintersection 

! these are the regions of the refined source grid
  integer                                                  :: n_sregions
  type(Tregion),dimension(:),pointer                       :: sregion 

!.........................................................................
! the receivers
  integer                                                  :: n_receivers
  type(Treceiver),dimension(:),pointer                     :: receiver 

! the sources
  integer                                                  :: n_sources   ! the number of sources defined in the input
  integer                                                  :: n_sources_ppinc ! n_sources + # of virtual sources at receiver
                                                                            ! positions required for reflection matching
  type(Tsource),dimension(:),pointer                       :: source

!.........................................................................
! parameters of the inversion 
  integer                                                  :: n_inv_parms        ! total # of inversion parameters
  integer                                                  :: n_inv_active       ! total # of active inversion parameters
  integer                                                  :: n_inv_vgrid        ! # of velocity grids to be solved for
  type(Tgrid_identifier) ,dimension(:),pointer             :: vgrids_to_be_inv   ! list of velocity grids to be solved for
  integer                                                  :: n_inv_iface        ! # of interfaces to be solved for
  integer,dimension(:),pointer                             :: ifaces_to_be_inv   ! list of interfaces to be solved for

  logical                                                  :: locate_source      ! solve for source position and time or not
  integer                                                  :: n_inv_source       ! # of sources to be solved for
  integer,dimension(:),pointer                             :: sources_to_be_inv  ! list of sources to be solved for

end module global_variables

!*************************************************************************************************************


module interface_definitions
! explicit interfaces for subroutines that have pointer/target arguments

interface
   subroutine propagate(regin,vtype)
     use type_definitions
     type(Tregion),target    :: regin
     integer                 :: vtype
   end subroutine
end interface

interface
   subroutine trace_ray_from_receiver(rec,s,ray)
     use type_definitions
     type(Tsource),target          :: s         
     type(Treceiver)               :: rec       
     type(Tray),target             :: ray       
   end subroutine
end interface

interface
   subroutine find_intersection(isec,iface,grid)
     use type_definitions
     type(Tintersection)               :: isec
     type(Tinterface)                  :: iface
     type(Tpropagation_grid),target    :: grid
   end subroutine
end interface

interface
   subroutine define_region(reg,itop,ibot,grid)
     use type_definitions
     type(Tregion),target     :: reg
     type(Tintersection)      :: itop,ibot
     type(Tpropagation_grid),target :: grid
   end subroutine
end interface

interface
   subroutine sweep_region_from_interface(reg,istart_in,vtype,s)
     use type_definitions
     type(Tregion)                          :: reg
     type(Tintersection),target             :: istart_in
     integer                                :: vtype
     type(Tsource)                          :: s
   end subroutine
end interface

interface
   subroutine sweep_sregion_from_interface(reg,istart_in,vtype)
     use type_definitions
     type(Tregion)                          :: reg
     type(Tintersection),target             :: istart_in
     integer                                :: vtype
   end subroutine
end interface

interface
   subroutine initialize_refined_source(s,sc,grid,reg,itop,ibot)
     use type_definitions
     type(Tsource) :: s
     type(Tsource) :: sc
     type(Tpropagation_grid) :: grid
     type(Tregion)         :: reg
     type(Tintersection),target   :: itop,ibot
   end subroutine
end interface


interface
   subroutine initialize_refined_source2(s,sc,grid,reg,itop,ibot)
     use type_definitions
     type(Tsource) :: s
     type(Tsource) :: sc
     type(Tpropagation_grid) :: grid
     type(Tregion)         :: reg
     type(Tintersection),target   :: itop,ibot
   end subroutine
end interface

end module interface_definitions

!**********************************************************************************************************

module mod_3dfm

  use type_definitions
  use global_variables
  use interface_definitions

end module mod_3dfm

module mod_3dfm_nointerfaces

  use type_definitions
  use global_variables

end module mod_3dfm_nointerfaces
nn_subsf.f
c------------------------------------------------------------------------
c
c	nn2d_setup - Performs all setup procedures for natural neighbour 
c		     interpolation routine nn2D.
c
c	Input:
c	        np			number of nodes
c	        nt_max			maximum number of triangles
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		np_max			maximum number of nodes 
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~20 in calling program)
c		nmax			maximum sum of the number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c		points(2,np)		array of node co-ordinates
c	        dmode			Delaunay calculation mode (integer)
c	        nmode			NN setup mode
c	        clockwise		logical for the vertices input order  
c               data(np)                data values at each node
c		nnn			Integer work array used by build_nv
c		nnlist			Integer work array used by build_nv
c		ntwork		        Integer work array used by build_nv
c		nohalt_hull		determines error response in routine
c					calculate_hulltriangles
c		eps 			tolerance parameter used by delaun
c					(see delaun for details)
c		vis_tlist		Integer work array used by delaun
c		vis_elist		Integer work array used by delaun
c		add_tlist		Integer work array used by delaun
c		nv_max			size of delaun work arrays vis_tlist
c					vis_elist, & add_tlist (passed to 
c					delaun for error checking)
c
c	Output:
c	        nt			number of triangles
c	        vertices(3,nt)		array of triangle vertices 
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        nh			number of triangles with an edge
c					on the convex hull
c		hulltriangles(nh)	array of triangles with an edge 
c					on the convex hull
c		loc			an initial guess triangle for point 
c					location routine `triloc' used by nn2D
c					(set somewhere near the centre)
c
c	Operation modes:
c
c		 The setup routine will perform different tasks depending
c		 on the input parameters dmode and nmode (see table below).
c		 Depending on the modes used some work arrays may be set 
c		 to size 1 to save memory. The "Memory Savings" column in the
c		 table below shows the dimension statement that may
c		 be used in the calling program if the routine is ONLY EVER 
c		 CALLED IN THE CORRESPONDING MODE. 
c
c		 PARAMETERS	ACTION	 		MEMORY SAVINGS
c
c		 nmode = 1	Delaunay only 		real*8 centres(3,1)
c							integer hulltriangles(1)
c		 nmode = 0	Delaunay + nn setup	
c		 nmode = -1	nn setup only		
c
c		 dmode > 0	Delaunay read in from   integer vis_tlist(1)
c				logical unit dmode.	integer vis_elist(1) 
c							integer add_tlist(1) 
c		 dmode = 0      Qhull used		Same as dmode > 0.
c		 dmode = -1	Delaun + X-sort         integer nnn(1)
c							integer nnlist(1)
c							integer ntwork(1)
c		 dmode = -2	Delaun + no sort        Same as dmode=-1
c
c		 dmode = 0 & nmode=1			integer neighbour(3,1)
c		 
c		 A call with nmode = -1 can only be made after a call 
c		 with nmode = 1.
c
c	Comments:
c
c		 If the arrays are used then they should be dimensioned
c		 in the calling program in the following way:
c
c		 real*8		points(2,np_max)
c		 real*8		centres(3,nt_max)
c		 integer	vertices(3,nt_max)
c		 integer	neighbour(3,nt_max)
c		 integer	hulltriangles(nh_max)
c		 
c		 Note: nh_max can usually be dimensioned much less than nt_max
c		 because hulltriangles, stores in a compact form, all 
c		 triangles with an edge on the convex hull. Except for
c		 very irregular point distributions nh << nt. If nh is 
c		 determined to be > nh_max then an error is reported and the
c		 program is halted (unless nohalt parameter is set to 1). 
c		 The array hulltriangles is only used by nn2Do see routine 
c		 calculate_hulltriangles. If nh_max = 1 then hulltriangles 
c		 is not calculated.
c
c		 The initial guess triangle 'loc' is set in nn_setup but
c		 at each call it will be set to the triangle found during
c		 the previous call to nn2D. The user may modify its value
c		 if the input point (x,y) is known to be in, or near, a
c		 particular triangle.
c
c		 If dmode > 0 the the deluanay tessellation is read in from
c		 logical unit `dmode' instead of being calculated internally 
c		 This can be useful if qhullf fails because
c		 of precision errors. The Deluanay may be determined
c		 externally to this program using a double precision version
c		 or another algorithm, e.g. Fortune's sweepline method.
c
c		 If 50 > dmode > 0 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ZERO, 
c		 which is the standard output format of codes qhull 
c		 (quickhull method) and voronoi (sweepline method).
c		 If clockwise = .true. (.false.) then the vertices are assumed 
c		 to be in clockwise (anti-clockwise) order. Note program
c		 qhull outputs vertices in anti-clockwise order while 
c		 voronoi in clockwise order. The internal format is 
c		 anti-clockwise and nodes numbered from ONE.
c
c		 If dmode => 50 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ONE,
c		 which is the output format of program del (using delaun). 
c
c		 Three other work arrays are produced as a `by product'
c		 of the routine build_nv which calculates the neighbour
c		 array. These must be dimensioned in the calling program in 
c		 the following way (unless delaun is used for calculating the
c		 Delaunay because it already determines the neighbour array)
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntwork(nmax) : dummy work array
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 Each of these are useful lists that describe features of
c		 the Voronoi diagram. Both nnlist and ntwork are stored in
c		 a compact format to avoid zeros. They are only used 
c		 in the setup routine and the memory may be freed once
c		 initialization is completed.
c		 
c
c		 Calls are made to: qhullf, ccentres, build_nv and 
c				    calculate_hulltriangles, delaun.
c		 
c					M. Sambridge, RSES, April 1994.
c  					        (Last modified 10/4/96)
c
c------------------------------------------------------------------------
c
	Subroutine nn2d_setup
     &             (np,nt_max,nh_max,np_max,nnpn_max,nmax,
     &              points,dmode,nmode,clockwise,data,nt,vertices,
     &              centres,neighbour,nh,hulltriangles,nohalt_hull,
     &              loc,nnn,nnlist,ntwork,
     &              eps,nv_max,vis_tlist,vis_elist,add_tlist,
     &		    lt_work,ln_work)

	real*8		points(2,*)
	real*8		centres(3,*)
        real*8          data(*)
	real*8		eps
	integer		vertices(3,*)
	integer		neighbour(3,*)
        integer		hulltriangles(*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntwork(*)
        integer         vis_tlist(*)
        integer         vis_elist(*)
        integer         add_tlist(*)
        integer         dmode,nmode
	logical*1	lt_work(*)
	logical*1	ln_work(*)
	logical		nnwrite
	logical		clockwise
	logical		ldummy
        logical         timing

        common/nnswitches/nnwrite,lud

        common/timing/timing,t_loc,t_int,t_setup

        if(timing)a = cputime(t1,t2)

        if(nmode.eq.1.or.nmode.eq.0)then

           if(dmode.eq.0)then
c                                       calculate Delaunay using qhull 
 
              call qhullf(np,2,2,nt_max,0,points,nt,vertices)

           else if(dmode.eq.-1.or.dmode.eq.-2)then

c					sort the points in ascending x order
c					and rearrange data points also
	      if(dmode.eq.-1)then

c		  write(*,*)' X sort in progress'
	          call hpsort_d(np,1,points,data)
c		  write(*,*)' X sort done'
 		  write(*,*)' Input points have been sorted in order '
 		  write(*,*)' of first co-ordinate'

	      end if

c                                       calculate Delaunay using delaun 
 
              call delaun (points,np,neighbour,vertices,nt,nt_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     0,ldummy,0,0,0)

	   else
c                                       read in Delaunay vertices

              nt = 0
              i1 = 1
              i2 = 2
              if(clockwise)then
                 i1 = 2
                 i2 = 1
              end if
              read(dmode,*)
  1           read(dmode,*,end=3,err=2)
     &        vertices(i1,nt+1),vertices(i2,nt+1),vertices(3,nt+1)
              nt = nt + 1
              if(nt.ge.nt_max)then
                 write(*,*) 'Error in nn_setup: too many triangles'
                 write(*,*) 'Remedy: increase size of parameter nt_max'
                 write(*,*) '        in calling program.'
                 stop 
              end if
              go to 1
  2           write(*,*)
     &        'Error in nn_setup: read error in Delaunay input file'
              stop
  3           continue
     
           end if
 
c					adjust array vertices to
c					range from nodes 1 to np
c
	   if(dmode.ge.0.and.dmode.lt.50)then
	      do 5 i = 1,nt
	         vertices(1,i) = vertices(1,i) + 1
	         vertices(2,i) = vertices(2,i) + 1
	         vertices(3,i) = vertices(3,i) + 1
 5            continue
	   end if

	end if
c
c					Perform set up for nn interpolation
c
        if(nmode.eq.0.or.nmode.eq.-1)then

c					set initial guess for 
c					triangle location procedure
           loc = nt/2

c                                       Calculate Circumcentres

           call ccentres(points,vertices,nt,centres)

c                                       Build neighbour matrix
c					(if not already built)

           if(dmode.ge.0)then
              call build_nv
     &        (np,vertices,nt,np_max,nmax,
     &         neighbour,nnn,nnlist,ntwork)
           end if

c					calculate hulltriangles

           if(nh_max.gt.1) call calculate_hulltriangles
     &     (neighbour,nt,nh_max,hulltriangles,nh,nohalt_hull)

c					initialize logical work arrays 
           do i=1,nt
              lt_work(i) = .false.
           end do
           do i=1,np
              ln_work(i) = .false.
           end do

	end if

        if(timing)then
           a = cputime(t1,t2)
           t_setup = t_setup + t1
        end if

	return
	end
c
c------------------------------------------------------------------------
c

c------------------------------------------------------------------------
c
c	Heapsort - Sorts an array into ascending order.
c		   Modified from numerical recipes 2, to use a 
c		   two dimensional double precision array.
c
c		   Sorting is done on index ka (ka=1 or 2).
c
c------------------------------------------------------------------------
c
c
      SUBROUTINE hpsort_d(n,ka,ra,da)
      INTEGER n
c     REAL ra(n)
      REAL*8 ra(2,n)
      REAL*8 da(n)
      INTEGER i,ir,j,l,ka
c     REAL rra
      REAL*8 rra(2)
      REAL*8 dda
      kb = 1
      if(ka.eq.1)kb = 2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra(ka)=ra(ka,l)
          rra(kb)=ra(kb,l)
	  dda = da(l)
        else
          rra(ka)=ra(ka,ir)
          rra(kb)=ra(kb,ir)
	  dda = da(ir)
          ra(ka,ir)=ra(ka,1)
          ra(kb,ir)=ra(kb,1)
          da(ir)=da(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(ka,1)=rra(ka)
            ra(kb,1)=rra(kb)
            da(1)=dda
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(ka,j).lt.ra(ka,j+1))j=j+1
          endif
          if(rra(ka).lt.ra(ka,j))then
            ra(ka,i)=ra(ka,j)
            ra(kb,i)=ra(kb,j)
            da(i)=da(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(ka,i)=rra(ka)
        ra(kb,i)=rra(kb)
        da(i)=dda
      goto 10
      END
c
c------------------------------------------------------------------------
c
c	pangle - pseudo angle routine 
c
c		 returns a number between 0 and 360 which is NOT the
c		 angle made by the line from p1 to p2 with the horizontal
c		 but which has the same order properties as that angle,
c		 i.e. has the same order of angles as arctan dy/dx.
c	 	 This function involves only simple products and quotients.
c
c		 From Sedgewick (1990) `Algorithms in C' (Addison Wesley)
c
c						M. Sambridge 1996.
c
c------------------------------------------------------------------------
c
	Function pangle(p1,p2)

	real*8		p1(2)
	real*8		p2(2)

	dx = p2(1) - p1(1)
        ax = abs(dx)
	dy = p2(2) - p1(2)
        ay = abs(dy)
        t = 0.
        a = ax+ay
	if(a.ne.0.)then
           t = dy/a
        end if
        if(dx.lt.0.)then
          t = 2-t
        else if(dy.lt.0.)then 
          t = 4+t
        end if
        pangle = t*90

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum - calculates circum-centre of three points
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c	Comments:
c
c		 Solves 3x3 linear system of equations.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum(pa,pb,pc,centre)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_d - calculates circum-centre of three points and 
c	           derivatives of circum-centre with respect to 
c		   co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_d(pa,pb,pc,centre,vx,vy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2)
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2

c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_dd - calculates circum-centre of three points and 
c		    1st and 2nd derivatives of circum-centre with 
c		    respect to co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_dd(pa,pb,pc,centre,vx,vy,vxx,vyy,vxy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2),vxx(2),vyy(2),vxy(2)
c
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2
c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2
c						Second derivatives
	   f11 = 2*vx(1) - 1.
	   f22 = 2*vy(2) - 1.
	   f12 = vy(1) + vx(2)

           vxx(1) = f11*denum1
           vxx(2) = f11*denum2
           vyy(1) = f22*denum1
           vyy(2) = f22*denum2
           vxy(1) = f12*denum1
           vxy(2) = f12*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c						heapsort modified to 
c						sort two arrays  
c
c------------------------------------------------------------------------
c
      SUBROUTINE hpsort_two(n,ra,rb)
      INTEGER n
      REAL ra(n)
      REAL*8 rb(2,n)
      INTEGER i,ir,j,l
      REAL rra
      REAL*8 rrb1,rrb2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
          rrb1=rb(1,l)
          rrb2=rb(2,l)
        else
          rra=ra(ir)
          rrb1=rb(1,ir)
          rrb2=rb(2,ir)
          ra(ir)=ra(1)
          rb(1,ir)=rb(1,1)
          rb(2,ir)=rb(2,1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            rb(1,1)=rrb1
            rb(2,1)=rrb2
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            rb(1,i)=rb(1,j)
            rb(2,i)=rb(2,j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra
        rb(1,i)=rrb1
        rb(2,i)=rrb2
      goto 10
      END
c
c------------------------------------------------------------------------
c
c					Numerical Recipes routine index
c
c------------------------------------------------------------------------
c
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
c
c------------------------------------------------------------------------
c
c	second_v_area - calculates the area of a second-order Voronoi 
c			cell using an un-ordered list of vertices and a
c		        a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		n			number of vertices
c
c	Output:
c		area			area of polygon X 2
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Comments:
c		 The vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area(x,n,p,a,area)

      real*8	x(2)
      real*8	p(2,*)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
      area = abs(area)
      
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_d - calculates the area of a second-order Voronoi 
c			  cell using an un-ordered list of vertices and a
c		          a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first two vertices 
c					of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx,
c					      df(2) = df/dy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_d(x,n,p,dp,a,area,df)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	df(2)
      real*4    theta
      real*4    a(*)
c
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
      end if

      area = abs(area)
       
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_dd - calculates the area of a second-order Voronoi 
c			   cell using an un-ordered list of vertices and a
c		           a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first 
c					two vertices of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx
c					      df(2) = df/dy
c		ddf(3)			second derivatives of area X 2
c					      ddf(1) = d2f/dxx
c					      ddf(2) = d2f/dyy
c					      ddf(3) = d2f/dxy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two and xplot routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_dd(x,n,p,dp,ddp,a,area,df,ddf)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	ddp(6,2)
      real*8	df(2),ddf(3)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
        a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

c						calculate 2nd derivatives
     
      ddf(1) =   ddp(1,1)*d222n + ddp(2,1)*d1n12 
     &         + ddp(1,2)*d2321 + ddp(2,2)*d1113
     &         + 2.*(dp(1,1)*dp(2,2) - dp(2,1)*dp(1,2))
      ddf(2) =   ddp(3,1)*d222n + ddp(4,1)*d1n12 
     &         + ddp(3,2)*d2321 + ddp(4,2)*d1113
     &         + 2.*(dp(3,1)*dp(4,2) - dp(4,1)*dp(3,2))
      ddf(3) =   ddp(5,1)*d222n + ddp(6,1)*d1n12 
     &         + ddp(5,2)*d2321 + ddp(6,2)*d1113
     &         + dp(1,1)*dp(4,2) - dp(2,1)*dp(3,2)
     &         + dp(3,1)*dp(2,2) - dp(4,1)*dp(1,2)

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
        ddf(1) = -ddf(1)
        ddf(2) = -ddf(2)
        ddf(3) = -ddf(3)
      end if

      area = abs(area)
       
      return
      end
c
c------------------------------------------------------------------------
c
c	delaun - calculates delaunay triangulation incrementally 
c	 	 for a set of points in 2-D using a variation of
c		 Lawson's algorithm.
c
c	Input:
c		points(2,np)		array of node co-ordinates
c		num			number of nodes to be used
c               vis_tlist(nv_max)       List of triangles visible from 
c					current point.
c               vis_elist(nv_max)       List of edges visible from 
c					current point.
c               add_tlist(nv_max)       work array used by routine addpoint
c               eps                     distance from an interface for a
c                                       a point to be considered on an
c                                       interface (real*8). Prevents zero
c                                       area triangles resulting from rounding
c                                       error when nodes are co-linear.
c		nv_max			size of work arrays
c		mode			(=0,1,2,3) operation mode (see below)
c		inactive(np)		logical array. If mode=1 then the i-th
c					node is ignored if active(i) = .true.
c		nfirst			If mode=3 then nfirst is the first
c					node added to an existing triangulation
c		itstart			If mode=3 then itstart is a first
c					guess triangle containing node first node
c		subset(np)		logical array. If mode=2 then only
c					the nodes (subset(i),i=1,num) are used.
c
c	Output:
c               v(3,*)           	array of triangle vertices
c               numtri                  number of triangles in current
c                                       triangulation.
c               e(3,*)                  adjacency matrix of neighbouring
c                                       triangles. e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2) 
c					in triangle j, stored counterclockwise 
c					about j.  
c                                       (This is the `opposite' definition)
c
c	Comments:
c
c       This routine calculates the Delaunay triangulation of a set of nodes 
c	using a variation of Lawson's method.  Each node is added sequentially 
c	and the Delaunay triangulation is updated. If the new node is inside 
c	the convex hull of the existing triangulation then the standard Lawson 
c	method is used. If it is outside then the list of triangle edges 
c	which are visible from the new point is calculated using routine 
c	visiblelist and each of these is used as the start of the swapping 
c	routine addpoint.
c
c	Four different operation modes are allowed.
c
c	MODE = 0:
c	The `standard' mode. All nodes from 1 to num are included. The arrays 
c	`subset' and `inactive' are unused and may be set to dummy variables
c       (saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 1:
c	All nodes from 1 to num are included except those for which
c	inactive(i) is set to true. The array `subset' is unused and may
c	be set to a dummy variable (saving memory). The variables nfirst 
c	and itstart are also unused.
c
c	MODE = 2:
c	Only nodes from subset(1) to subset(num) are included. The array
c	`inactive' is unused and may be set to a dummy variable
c	(saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 3:
c	Used to add nodes from nfirst to num to an existing triangulation.
c	Nodes for which inactive(i) is set to true are ignored.
c	The array `subset' is unused and may be set to a dummy variable.
c
c       The performance may be sensitive to the order in which the nodes are
c       added so these can be sorted before calling this routine if desired.
c
c	This routine was converted to use the `opposite' definition of
c	the adjacency matrix on 30/1/96.
c
c
c	Calls are made to triloc_del,visiblelist,insert_point,addpoint.
c
c					         M. Sambridge, Dec. 1994.
c					Modified by J. Braun, Sept. 1995.
c					(last change 30/1/96: multiple modes,
c					 and uses opposite definition of
c					 adjacency matrix)
c
c------------------------------------------------------------------------
c
	Subroutine delaun (points,num,e,v,numtri,numtri_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     mode,inactive,nfirst,itstart,subset)

	real*8		points(2,*)
	real*8		x,y
        real*8          eps,del1,del2,del
 	integer		vis_tlist(*),vis_elist(*),add_tlist(*)
 	integer		v(3,*)
 	integer		e(3,*)
 	integer		subset(*)
	integer		t,p,ccw
        logical         out
	logical		newpoint
	logical		inactive(*)

        if (mode.eq.0.or.mode.eq.1.or.mode.eq.2) then

c					We are calculating Delaunay 
c					of all input points or a
c					subset of all input points

c					find first two active nodes
           if(mode.eq.0)then
              i1=1 
              i2=2 
              nodestart=3
           else if(mode.eq.1)then
              i1=0
              i2=0
              do i=1,num
                 if (i2.ne.0) goto 2222
                 if (i1.ne.0.and..not.inactive(i)) i2=i
                 if (i1.eq.0.and..not.inactive(i)) i1=i
              end do
 2222         continue
              nodestart = i2+1
           else if(mode.eq.2)then
              i1 = subset(1)
              i2 = subset(2)
              nodestart = 3
           end if
c                                       Find three non-colinear points
c                                       to form the first triangle 
           v(1,1) = i1
           v(2,1) = i2
           do 10 j=nodestart,num
              i = j
              if(mode.eq.2)then
                  i = subset(j)
              else if(mode.eq.1.and.inactive(i))then
                  go to 10
              end if
              istart=i
	      del1 = (points(2,i1)-points(2,i))
     &              *(points(1,i2)-points(1,i))
	      del2 = (points(1,i1)-points(1,i))
     &              *(points(2,i2)-points(2,i))
              del = del1-del2
              if(dabs(del).gt.eps) goto 11111
 10        continue
           stop 'all input data are in a line...'
11111      v(3,1) = istart

c					Initialize adjacency matrix
 	   e(1,1) = 0
 	   e(2,1) = 0
 	   e(3,1) = 0
c					Ensure initial triangle 
c					is in ccw order
c					
 	   if(ccw(points(1,v(1,1)),
     &            points(1,v(2,1)),
     &            points(1,v(3,1)),k).eq.-1)then
                  itemp = v(1,1)
                  v(1,1) = v(2,1)
                  v(2,1) = itemp
c	          write(*,*)' initial triangle was cw'
 	    end if
	
c					Initialize variables
 	    numtri = 1
	    t = 1

        else if (mode.eq.3) then
c					We are adding nodes to an 
c					existing triangulation
c					Perform initialization
           nodestart=nfirst
           t = itstart
           istart = 0
           if(t.le.0.or.t.gt.numtri)t=1

        end if 
c					Incrementally update the 
c					Delaunay triangulation

 	do 100 j=nodestart,num


           p = j
           if(mode.eq.1.or.mode.eq.3)then
             if (inactive(j)) goto 100
           else if(mode.eq.2)then
	     p = subset(j)
           end if
             
           if(p.eq.istart)go to 100

	   x = points(1,p)
	   y = points(2,p)

c					locate triangle 
c					containing current node

	   call triloc_del(x,y,points,v,e,t,eps,out,ipos,iface)


	   if(out)then

c					point is outside of convex hull, 
c					so find list of edges that are 
c					visible from current point 

	      call visiblelist(points,e,v,x,y,t,ipos,eps,
     &                         vis_tlist,vis_elist,nvis)

c					for each visible edge 
c					start swapping algorithm

              newpoint = .true.

	      if(nvis.gt.nv_max)then
                 write(*,*)' Error in subroutine delaun:'
                 write(*,*)' Too many visible triangles
     &                       from current point'
                 write(*,*)' Remedy: increase size of parameter nv_max'
                 write(*,*)'         in calling program'
                 write(*,*)'         Number of visible triangles '
                 write(*,*)'         for this point             =',nvis
                 write(*,*)'       Current value of nv_max    =',nv_max
                 stop
              end if

	      do 60 i=1,nvis
                 t = vis_tlist(i)
                 ipos = vis_elist(i)
                 jpos = mod(vis_elist(i),3)+1
c	         write(6,*)' visible t =',t,' node',v(ipos,t),v(jpos,t)
 	         call addpoint
     &           (points,e,v,p,t,ipos,numtri,newpoint,add_tlist) 
                 newpoint = .false.
 60           continue

	   else
 
c	      write(6,*)' point located in triangle',t

c					add node to inside of convex hull
c					using swapping algorithm

  	      call insertpoint(points,e,v,p,t,numtri,iface) 

           end if

           if (numtri.gt.numtri_max) then
              write (*,*) 'Error in subroutine delaun:'
              write (*,*) 'Too many triangles'
              write(*,*)' Remedy: increase size of parameter numtri_max'
              write(*,*)'         in calling program'
              write(*,*)'         Number of triangles '
              write(*,*)'         for this point             =',numtri
              write(*,*)'         Current value of numtri_max    =',
     &                  numtrimax_max
              stop
           endif

 100    continue

	return
	end
c------------------------------------------------------------------------
c
c	Subroutine visiblelist - calculates all sides of triangles 
c		                 visible from the point (x,y), which is
c			         outside of the convex hull.
c			
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		x,y			Co-ordinates of test point p
c		t			Index of any triangle on hull 
c					that is visible from point p.
c					(Usually given by routine Triloc_del.)
c		tpos			Position of edge in triangle t
c					(using Sloan's adjacency convention)
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c		nvis			Number of triangles visible from point p
c		vis_tlist		List of triangles visible from p
c		vis_elist		List of edges visible from p
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c		 This routine was converted from using Sloan's definition of
c		 the adjacency matrix to the `opposite' definition on 30/1/96.
c
c	Calls routine visible.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c		
	Subroutine visiblelist
     &             (points,neighbour,vertices,x,y,t,tpos,eps,
     &              vis_tlist,vis_elist,nvis)

	real*8		points(2,*)
	real*8		x,y
	real*8		eps
 	integer		vertices(3,*)
 	integer		neighbour(3,*)
 	integer		vis_tlist(*),vis_elist(*)
 	integer		t,tpos,pos,t1,t2,edg,tnew
	logical		visible
	logical		special
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/

	nvis = 1
	vis_tlist(1) = t
	vis_elist(1) = tpos
    	inode = vertices(tpos,t)
cd      write(6,100)t,inode,vertices(mod(tpos,3)+1,t)
    	pos = c1(tpos)
      	jnode = vertices(pos,t)
    	t1 = neighbour(c2(pos),t)
        special = .false.
	if(t1.eq.0)then
	  t1 = t
          tnew = 0
          special = .true.
	end if

  5     continue
        if(.not.special)then
           pos = edg(t1,jnode,vertices)
           tnew = neighbour(c2(pos),t1)
cd         write(6,*)' tnew =',tnew,' t1',t1,' jnode',jnode
        end if
        special = .false.
        if(tnew.eq.0)then
  6        continue
  	   if(visible(x,y,points,vertices,t1,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t1
	      vis_elist(nvis) = pos
cd            write(6,100)t1,jnode,vertices(mod(pos,3)+1,t1)
	   else
cd	      write(6,200)t1,jnode,vertices(mod(pos,3)+1,t1)
              go to 10
	   end if
           pos = c1(pos)
	   jnode = vertices(pos,t1)
    	   tnew = neighbour(c2(pos),t1)
	   if(tnew.eq.0) go to 6
           t1 = tnew
	   go to 5 
        else
cd	   write(6,300)t1,jnode,vertices(mod(pos,3)+1,t1)
	   t1 = tnew
	   go to 5 
	end if

  10	jnode = inode
    	pos = c2(tpos)
    	t2 = neighbour(c2(pos),t)
        special = .false.
	if(t2.eq.0)then
	  t2 = t
          tnew = 0
          special = .true.
	end if

  15    continue
        if(.not.special)then
           pos = c2(edg(t2,jnode,vertices))
           tnew = neighbour(c2(pos),t2)
        end if
        special = .false.
        if(tnew.eq.0)then
  16       continue
  	   if(visible(x,y,points,vertices,t2,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t2
	      vis_elist(nvis) = pos
cd            write(6,100)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
	   else
cd	      write(6,200)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
              go to 20
	   end if
	   jnode = vertices(pos,t2)
    	   pos = c2(pos)
    	   tnew = neighbour(c2(pos),t2)
	   if(tnew.eq.0)go to 16
	   t2 = tnew
	   go to 15
        else
cd	   write(6,300)t2,jnode,vertices(mod(pos,3)+1,t2)
	   t2 = tnew
	   go to 15
	end if

 20	continue
	      
c 100    format
c     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is visible')
c 200    format
c     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not visible')
c 300    format
c     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not on convex hull')

	return
	end
c
c------------------------------------------------------------------------
c
c	Function visible - determines whether the triangle t is visible
c		           from the point p on edge tpos.
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        t			Triangle to be tested
c	        tpos			Edge to be tested in triangle t 
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        visible			Logical: = true if edge is visible  
c	                                         = false if edge is not visible
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c	Calls no other routines.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Function visible(x,y,points,vertices,t,tpos,eps)

	real*8		points(2,*)
	real*8		del1,del2
	real*8		x,y
 	integer		vertices(3,*)
 	integer		t,tpos
	logical		visible
	real*8		eps,del
	integer		c1(3)
        save            c1
	data		c1/2,3,1/

        j = c1(tpos)
c						test edge tpos in triangle t
        i1 = vertices(tpos,t)
        i2 = vertices(j,t)
        del1 = (points(2,i1)-y)*(points(1,i2)-x)
        del2 = (points(1,i1)-x)*(points(2,i2)-y)
	del = del1-del2
        if(del.gt.eps)then
           visible = .true.
	else
           visible = .false.
	end if

	return
	end
c
c------------------------------------------------------------------------
c
c	addpoint - inserts a point into an existing delaunay triangulation 
c		   when point is outside of triangulation (but attached to
c		   triangle t) using the stacking procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle on convex hull visible 
c					from input point 
c		numtri			number of triangles in current
c					triangulation.
c		tpos			position of start node in triangle t
c		tri			list of triangles visible from point p
c		newpoint		logical = true if t is the first
c					triangle on the hull visible from p
c
c	Output:
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	When newpoint = false then there are multiple triangles
c	from the new point to the convex hull, and addpoint must
c	be called once for each attached triangle. In this case
c	the initialization of the adjacency list includes the
c	neighouring triangles already processed by addpoint, i.e.
c	those from point p to the hull.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine addpoint (points,e,v,p,t,tpos,numtri,newpoint,tri) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t,tpos,ccw
	logical		swap
	integer		tri(*)
	logical		newpoint
        save 		ip,lp
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/

	if(newpoint)then
	   ip = 0
	   lp = 0
        end if

c			Add new node to existing triangulation

c					create new triangle
        numtri = numtri + 1
        v1 = v(tpos,t)
        v2 = v(c1(tpos),t)
        if(ccw(points(1,v1),
     &         points(1,v2),
     &         points(1,p),k).eq.-1)then
               itemp = v1
               v1 = v2
               v2 = itemp
        end if
        v(1,numtri) = p
        v(2,numtri) = v1
        v(3,numtri) = v2

c					initialize adjacency list including
c					neighbouring triangles attached
c					from the point to the hull.

        e(c2(1),numtri) = 0
        e(c2(2),numtri) = t
        e(c2(3),numtri) = 0

c				
        if(.not.newpoint)then
           do 10 j=1,lp
              k = tri(j)
              if(v(2,k).eq.v1)then
c                write(6,*)' v1 match with node 2'
c                write(6,*)' current triangle',numtri,' new',k
c                write(6,*)' nodes:',v(1,k),v(2,k),v(3,k) 
c                write(6,*)' e mat:',e(c2(1),k),e(c2(2),k),e(c2(3),k) 
                 e(c2(1),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v1)then
                 e(c2(1),numtri) = k
                 e(c2(3),k) = numtri
              end if
              if(v(2,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(3),k) = numtri
              end if
 10        continue
        end if

c
c					initialize stack

 	call stackinit

c                                       update adjacency list
c                                       for triangle on old boundary
        e(c2(tpos),t) = numtri


c					add new triangle on stack

	call push(numtri)

c					loop while stack is not empty

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl) 
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = l
	   end if

	   if(b.ne.0)then
              call push(r)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = r
           end if
	   if(c.ne.0) e(edg(c,l,e),c)=r

	else

c					record triangle attached to p
	   ip = ip + 1
           tri(ip) = l

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
        call stackflush()

	lp = ip

c       write(6,*)' Number of triangles attached to last point',ip
c1	write(6,*)(tri(i),i=1,ip)
c	write(6,*)' triangles attached to last point on hull'
c       do 100 i=1,ip
c          it=tri(i)
c          do 101 k=1,3
c             l=mod(k,3)+1
c             if(e(k,it).eq.0)then
c                write(6,*)' t',it,' edge ',v(k,it),v(l,it)
c             end if
c 101      continue 
c 100   continue 
c

	return
	end
c
c------------------------------------------------------------------------
c
c	insertpoint - inserts a point into an existing delaunay triangulation 
c		      (when new point is inside triangle t) using the stacking 
c		      procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle containing input point
c		numtri			number of triangles in current
c					triangulation.
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Output:
c
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The new point is assumed to be inside the convex hull of the
c	existing triangulation.
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine insertpoint (points,e,v,p,t,numtri,iface) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t
	logical		swap
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/


c					add new node to existing triangulation
        if(iface.eq.0)then
	   a = e(c2(1),t)
	   b = e(c2(2),t)
	   c = e(c2(3),t)
	   v1 = v(1,t)
	   v2 = v(2,t)
	   v3 = v(3,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+2	
	   e(c2(2),t) = a	
	   e(c2(3),t) = numtri+1	
        

c					create new triangles
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v2
	   v(3,numtri)=v3
	   e(c2(1),numtri)=t
	   e(c2(2),numtri)=b
	   e(c2(3),numtri)=numtri+1
	   numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=numtri-1
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        else
           j = iface
           k = c1(j)
           i = c1(k)
	   a = e(c2(i),t)
	   b = e(c2(j),t)
	   c = e(c2(k),t)
	   v1 = v(i,t)
	   v2 = v(j,t)
	   v3 = v(k,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+1
	   e(c2(2),t) = a	
	   e(c2(3),t) = 0 	

c					create new triangle
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=0
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        end if
  
c
c					initialize stack

 	call stackinit

c					add new triangles on stack
c					and update adjacency list

	if(a.ne.0)call push(t)

	if(b.ne.0)then
          e(edg(b,t,e),b)=numtri-1
          call push(numtri-1)
        end if

	if(c.ne.0)then
          e(edg(c,t,e),c)=numtri
          call push(numtri)
        end if
c					loop while stack is not empty

        if(a.eq.0.and.b.eq.0.and.c.eq.0)go to 100

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl)
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   end if
	   if(b.ne.0) call push(r)
	   if(c.ne.0) e(edg(c,l,e),c)=r

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
 100    continue
        call stackflush()

	return
	end
c
c------------------------------------------------------------------------
c
c	Function edg - finds edge in triangle l which is adjacent 
c		       to triangle k.
c
c		       (From Sloan 1987)
c
c------------------------------------------------------------------------
c
	Function edg(l,k,e)
c
	integer		l,k,i,e(3,*),edg
c
	do 10 i=1,3
	   if(e(i,l).eq.k)then
              edg = i
              return
           end if
 10     continue

	write(*,*)' ***Error in function edg***'
	write(*,*)' ***Triangles not adjacent***'
	write(*,*)' triangle = ',l,' looking for triangle',k

	stop
	end
c
c------------------------------------------------------------------------
c
c	logical function swap - checks to see if point p lies 
c			        inside circumcircle about points p1,p2,p3
c				using the algorithm of Cline and Renka
c				(see Sloan 1987).
c
c------------------------------------------------------------------------
c
	Function swap(p1,p2,p3,p)

	logical		swap

	real*8		p(2),p1(2),p2(2),p3(2)
	real*8		x13,y13,x23,y23,x1p,y1p,x2p,y2p
	real*8		cosa,cosb,sina,sinb

	x13=p1(1)-p3(1)
	y13=p1(2)-p3(2)
	x23=p2(1)-p3(1)
	y23=p2(2)-p3(2)
	x1p=p1(1)-p(1)
	y1p=p1(2)-p(2)
	x2p=p2(1)-p(1)
	y2p=p2(2)-p(2)

	cosa = x13*x23 + y13*y23
	cosb = x2p*x1p + y1p*y2p

	if((cosa.ge.0.d0).and.(cosb.ge.0.d0))then
            swap = .false.
        else if((cosa.lt.0.d0).and.(cosb.lt.0.d0))then
            swap = .true.
        else
            sina=x13*y23-x23*y13
            sinb=x2p*y1p-x1p*y2p
	    if((sina*cosb+sinb*cosa).lt.0.d0)then
                swap = .true.
            else
                swap = .false.
            end if
        end if

	return
	end
c
c------------------------------------------------------------------------
c
c	Triloc_del - locates the triangle containing point x,y
c
c	Input:
c		x,y			co-ordinates of input points	
c		points(2,np)		array of node co-ordinates	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c	        loc			first guess of triangle containing
c					(x, y).
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c	        out			=true if (x,y) is outside of
c					the convex hull, otherwise = false. 
c	        k			index of face through which the
c					algorithm last passed (used by
c					routine visbilelist if out = .true.)
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Comments:
c		 If (x,y) is outside convex hull loc is a visible triangle
c		 on the hull, out is set to .true., and k is set to the
c		 index of the face of triangle loc visible from the input point
c		 (used as a starting point by the routine visiblelist)
c
c		 This version also returns the parameter iface. 
c		 If iface .ne. 0 then the input point is on the face of 
c		 triangle t between nodes iface and mod(iface,3)+1 and 
c		 it is also on the convex hull.
c
c		 A point is assumed to be on the edge (or its extension)
c		 between two nodes if it is inside the triangle at a 
c		 distance >= eps.
c
c		 Can be extended to higher dimensions using a similar
c		 stepping mechanism but without angular test.
c
c	         This routine was converted from using Sloan's definition of
c	         the adjacency matrix to the `opposite' definition on 30/1/96.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine Triloc_del
     &                       (x,y,points,vertices,neighbour,loc,eps,
     &                        out,k,iface)
c
	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		p1,p2
        logical		out
	real*8		x,y,del1,del2,del
	real*8		eps
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/
        logical         new

	out = .false.
        new = .true.
        ic = 0

 10     continue
c					point is outside convex hull
        if(out)return
        iface = 0

        do 20 i=1,3
	   j = c1(i)
c	   k = c1(j)
c					definition of adjacency matrix

c					use Sloan's 
c					definition of adjacency matrix
	   k = i

           p1 = vertices(i,loc)
           p2 = vertices(j,loc)
	   del1 = (points(2,p1)-y)*(points(1,p2)-x)
	   del2 = (points(1,p1)-x)*(points(2,p2)-y)
           del = del1-del2
 	   if(dabs(del).le.eps)then
              iface = i
	   else if(del.gt.0.d0)then
	      if(neighbour(c2(k),loc).eq.0)then
                 out = .true.
	      else
	         loc = neighbour(c2(k),loc)
	      end if
              if(.not.new.and.loc.eq.loc1)then
                 write(*,100) 
                 write(*,*)' Current triangle:',loc, 
     &           ' last three:',loc1,loc2,loc3
                 write(*,*)' New point      x:',x,' y:',y
                 write(*,*)' Triangle ',loc,
     &           ' v:',(vertices(j,loc),j=1,3),
     &           ' n:',(neighbour(c2(j),loc),j=1,3)
c                write(*,*)' del',del,' del1',del1,' del2',del2
                 write(*,101) 
                 stop
              end if
              if(new)then
                ic = ic + 1
                if(ic.eq.3)new = .false.
              end if
              loc1 = loc2
              loc2 = loc3
              loc3 = loc
	      go to 10
	   end if
 20     continue
	
c						check if input point is
c						on the convex hull
c
        if(neighbour(c2(iface),loc).ne.0)iface = 0

c       if(iface.ne.0)then
c          j = mod(iface,3)+1
c          jj = vertices(iface,loc)
c          kk = vertices(j,loc)
c          write(*,*)' point on triangle between nodes ',
c    &               jj,' and',kk
c          write(*,*)' point is on the convex hull'
c       end if

 100    format(/'Error in subroutine Triloc_del:',//
     &  ' Infinite loop detected in walking triangle algorithm',/,
     &  ' Probably due to rounding error creating a flat triangle'/)

 101    format(/1x,'Remedy: '/
     &  ' Either increase size of parameter eps in calling routine '/
     &  ' or re-order input points by running program nn_hull '/)

	return
	end
c
c------------------------------------------------------------------------
c
c	Function ccw - used to test the orientation of three points
c
c		Input : 	points p1,p2 and p3 (vectors 2x1)
c				(e.g. p(1,2) = x co-ordinate of p2)
c
c		Output: 	ccw,I
c
c		ccw    k
c	  	 1     0	:The direction p1,p2,p3 is ccw (+ve)   
c	  	-1     0	:The direction p1,p2,p3 is  cw (-ve)   
c	  	 1     1	:p1,p2,p3 are colinear & p2 in middle  
c	  	-1     1	:p1,p2,p3 are colinear & p1 in middle
c	  	 0     1	:p1,p2,p3 are colinear & p3 in middle 
c
c
c				Calls no other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
      Integer Function ccw(p1,p2,p3,k)
c     
      real*8		p1(2),p2(2),p3(2)
      real*8		dx1,dx2,dy1,dy2,a,b
c     integer		ccw
c
      dx1 = p2(1) - p1(1)
      dx2 = p3(1) - p1(1)
      dy1 = p2(2) - p1(2)
      dy2 = p3(2) - p1(2)
      a = dx1*dy2
      b = dy1*dx2
      if (a.gt.b)then
         k = 0
         ccw = 1
      else if(a.lt.b)then
         k = 0
         ccw = -1
      else if(dx1*dx2.lt.0.0.or.dy1*dy2.lt.0.0)then
         k = 1
         ccw = -1 
      else if((dx1*dx1+dy1*dy1).lt.(dx2*dx2+dy2*dy2))then
         k = 1
         ccw = 1 
      else
         k = 1
         ccw = 0 
      end if
      return
      end
c
c------------------------------------------------------------------------
c
c	function theta - returns a real number between 0 and 360
c		         which has the same ordering as the angle
c			 between the line (a,b) and the horizontal.
c
c------------------------------------------------------------------------
c
	function theta(a,b)

	real*8		a(2),b(2)
	real*4		theta

	dx = b(1) - a(1)  
	ax = abs(dx)
	dy = b(2) - a(2)  
	ay = abs(dy)

	theta = 0.
	d = ax+ay
	if(d.ne.0)theta=dy/d

	if(dx.lt.0.0)then
           theta = 2.-theta
	else if(dy.lt.0.0)then
           theta = 4.+theta
	end if
	theta = theta*90.

	return
	end
c
c-----------------------------------------------------------------------------
c
c	Subroutine qhullf (np,i,j,nt_max,k,points,nt,vertices)
c
c-----------------------------------------------------------------------------
c
	Subroutine qhullf (np,i,j,nt_max,k,points,nt,vertices)
c
	real*8		points(2,*)
	integer		vertices(3,*)

	write(*,*)' '
	write(*,*)' Error in subroutine nn2d_setup'
	write(*,*)' qhull is not installed'
	write(*,*)' Delaunay triangulation must be either'
	write(*,*)' calculated with routine delaun (dmode=-1 or -2)'
	write(*,*)' or read in from a file (dmode>0; logical unit=dmode)'
	write(*,*)' '

	stop
	end
c
c------------------------------------------------------------------------
c
c       plot_c - dummy routine
c
c
c------------------------------------------------------------------------
c
        Subroutine plot_c(xs,ys,xs2,ys2)

	real*8		xs2,ys2
c                                               do nothing
        return
        end

c
c------------------------------------------------------------------------
c
c       plot_tc - dummy routine
c
c------------------------------------------------------------------------
c
        Subroutine plot_tc(n,points,vertices,centres)

        real*8          points(2,*)
        real*8          centres(3,*)
        integer         vertices(3,*)

c                                               do nothing
        return
        end
 
c
c ----------------------------------------------------------------------------
c
c       cputime - calls system dependent routine to calculate cputime
c		  since last call.
c
c       Calls dtime.
c						M. Sambridge, June 2001
c
c ----------------------------------------------------------------------------
c
        Function cputime(t1,t2)
        real*4 t1,t2
        real*4 tarray(2)
        
c        cputime = dtime(tarray)
c        t1 = tarray(1)
c        t2 = tarray(2)
        cputime=1.0

        return
        end
c------------------------------------------------------------------------
c
c	Ccentres - calculates centres of all Delaunay circumcircles
c
c
c	Input:
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c
c	Output:
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c
c	Comments:
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine ccentres(points,vertices,nt,centres)
c
	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x1,x2,x3,y1,y2,y3,x,y
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	integer		vertices(3,*)
c						Find centres of all
c						Delaunay Circumcircles
	do 5 i= 1,nt

	   x1 = points(1,vertices(1,i))
	   x2 = points(1,vertices(2,i))
	   x3 = points(1,vertices(3,i))
	   y1 = points(2,vertices(1,i))
	   y2 = points(2,vertices(2,i))
	   y3 = points(2,vertices(3,i))

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1
	   x = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1*0.5d0
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0*dy2m1)/
     &         (denom)

	   y = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1)*0.5d0)/ 
     &         (denom)

	   centres(1,i) = x
	   centres(2,i) = y
           x1 = x - x1
           y1 = y - y1
	   centres(3,i) = x1*x1 + y1*y1

 5	continue

	return
	end
c
c------------------------------------------------------------------------
c
c	build_nv - Builds neighbour array for Delaunay triangulation in 2-D.
c
c	Input:	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c		np_max			maximum number of nodes
c		nmax			maximum total number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c
c	Output:
c		neighbour(3,nt)		array of neighbouring triangles
c
c	Comments:
c		 Assumes input list of vertices in anticlockwise sequence
c		 and produces an anticlockwise list of neighbour triangles.
c		 The value of neighbour(i,j) is the index of the neighbouring
c		 triangle opposite node i in triangle j.
c
c		 Three temporary work arrays are used and must be dimensioned
c		 in the calling program in the following way:
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntwork(nmax) : dummy array (NOTE NOT triangles
c					  attached to each node)
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c					(using ideas by J.Braun)
c
c------------------------------------------------------------------------
c
	Subroutine build_nv
     &             (np,vertices,nt,np_max,nmax,
     &              neighbour,nnn,nnlist,ntwork)
c
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntwork(*)
	logical		nnwrite

	common/nnswitches/nnwrite,lud

	if(nnwrite)write(*,*)' Building neighbour v ...'
c
c					initialize neighbour list
	do 5 i = 1,3
	   do 4 j = 1,nt
	      neighbour(i,j) = 0
 4         continue
 5      continue
c					initialize work arrays
        do 6 i = 1,nmax
	   nnlist(i) = 0
	   ntwork(i) = 0
 6      continue

        do 7 i = 1,np
	   nnn(i) = 0
 7      continue

	do 10 it = 1,nt
	   i1 = vertices(1,it)
	   i2 = vertices(2,it)
	   i3 = vertices(3,it)
	   nnn(i1) = nnn(i1) + 1
	   nnn(i2) = nnn(i2) + 1
	   nnn(i3) = nnn(i3) + 1
 10     continue

c					turn nnn into a running sum
	itemp = nnn(1)+1
	nnn(1) = 1
	do 20 j = 2,np+1
	   itemp2  = itemp 
	   itemp   = itemp + nnn(j)+1
	   nnn(j) = itemp2 + 1
 20     continue
c       write(*,*)' size of array =',nnn(np+1)-1
c       write(*,*)' 3nt+np        =',3*nt+np

	if(nnn(np+1).ge.nmax)then
           write(*,*)'Error: array sizes too small in subroutine '
     &               ,'build_nv'
           write(*,*)'       maximum number of neighbours for all nodes'
           write(*,*)'       is too small: current value =',nmax
           write(*,*)'       Increase size of parameter nmax'
           write(*,*)'       to at least',nnn(np+1)
           write(*,*)'       This will be satisfied if nmax is set'
           write(*,*)'       to 3*nt_max+np_max in calling program' 
	   stop
	end if

	do 25 it = 1,nt
	   i1 = vertices(1,it) 
	   i2 = vertices(2,it) 
	   i3 = vertices(3,it) 
c						compare neighbours i1 i2
c						(remove go to ?)
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   jt = 0
	   do 30 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i2
	         ntwork(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 31
	      else if(nnlist(j).eq.i2.and.ntwork(j).ne.it)then
                 jt = ntwork(j)
	         go to 31
	      end if
  30       continue
  31       continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i2)
	      j2 = nnn(i2+1) - 1 
	      do 32 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntwork(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 33
	         end if
  32          continue
	   end if
  33       continue

	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i2
	      neighbour(3,it) = jt
	      k1 = vertices(1,jt)
	      k2 = vertices(2,jt)
	      k3 = vertices(3,jt)
	      if(k1.ne.i1.and.k1.ne.i2)then 
	         neighbour(1,jt) = it
	      else if(k2.ne.i1.and.k2.ne.i2)then 
	         neighbour(2,jt) = it
	      else
	         neighbour(3,jt) = it
	      end if
           end if
c						compare neighbours i1 i3
	   jt = 0
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   do 130 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntwork(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 131
	      else if(nnlist(j).eq.i3.and.ntwork(j).ne.it)then
                 jt = ntwork(j)
	         go to 131
	      end if
  130      continue
  131      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 132 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntwork(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 133
	         end if
  132         continue
	      end if
  133      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i3
	     neighbour(2,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i1.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i1.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if
c						compare neighbours i2 i3
	   jt = 0
	   j1 = nnn(i2)
	   j2 = nnn(i2+1) - 1 
	   do 230 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntwork(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 231
 	      else if(nnlist(j).eq.i3.and.ntwork(j).ne.it)then
                 jt = ntwork(j)
 	         go to 231
	      end if
  230      continue
  231      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 232 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i2
	            ntwork(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 233
	         end if
  232         continue
	      end if
  233      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i2 and i3
	     neighbour(1,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i2.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i2.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if

 25     continue

	if(nnwrite)write(*,*)' built neighbour v'

	return
	end
c------------------------------------------------------------------------
c
c       Calculate_hulltriangles - finds all triangles with a face on
c                                 the convex hull by searching through
c                                 the entries in the array neighbour.
c
c       Input:
c               neighbour(3,nt)         array of neighbouring tetrahedra
c               nt                      number of tetrahedra
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		nohalt			determines error response
c       Output:
c		hulltriangles(nh)	array of triangles with an edge
c					on the convex hull
c               nh                      number of tetrahedra with an edge
c                                       on the convex hull
c       Comments:
c
c                This routine fills up the array hulltriangles which
c                is only used by routine nn2Do, i.e the `pseudo-extension' 
c		 Watson's nn-interpolation method to points outside of the 
c		 convex hull. If nnext is set to false then hulltriangles
c		 is never used and the array can be set to size 1.
c
c		 If nohalt = 0 then the routine will stop with an error
c		 message if nh > nh_max. If nohalt .ne. 0 and nh > nh_max
c		 then it will return nh = -1. 
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, May 1995.
c
c------------------------------------------------------------------------
c
	Subroutine calculate_hulltriangles
     &             (neighbour,nt,nh_max,hulltriangles,nh,nohalt)
c
	integer		neighbour(3,*)
	integer		hulltriangles(*)

c                                               store list of triangles
c                                               which have an edge on the
c                                               convex hull.
c                                               (used by routine nn2D)
        nh = 1
        do 100 j = 1,nt
           if(neighbour(1,j).eq.0.or.
     &        neighbour(2,j).eq.0.or.
     &        neighbour(3,j).eq.0)then
              hulltriangles(nh) = j
              nh = nh + 1
              if(nh.gt.nh_max.and.nohalt.eq.0)then
                  write(*,*)' Error array storing outward facing '
                  write(*,*)' triangles on convex hull is too small.'
                  write(*,*)' Increase size of parameter nh_max'
                  stop
              else if(nh.gt.nh_max.and.nohalt.ne.0)then
                  nh = -1
                  return
              end if
           end if
 100    continue
        nh = nh -1
 
	return
	end
c
propagate.f90
! module mod_prop contains variables that relate to the fast marching through a
! region (grid between two interfaces including points on the interface)
!
module mod_prop
use mod_3dfm_nointerfaces

type(Tregion),pointer                              :: reg     ! the argument of the subroutine
type(Tpropagation_grid),pointer                    :: grid    ! grid on which the region is defined
real(kind=dp),dimension(:),pointer                 :: velocity
type(Tvelocity_grid),pointer                       :: velgrid

integer, dimension(:),pointer              :: node_from_tree_ind ! pointers from tree position to the associated node
integer, dimension(:),pointer              :: node_is_counted    ! array to keep track of whether a node is already
                                                                 ! counted when constructing a unique neighbour list
logical,dimension(:),pointer               :: suspect_time_gradient ! flag indicating that time gradient should be re-evaluated


real(kind=dp),dimension(:),pointer         :: cosla,sinla,coslo,sinlo ! precomputed sines and cosines of the regional nodes

type(Tinteger_coordinates),dimension(8)    :: concell            ! a list of cells connected to a given node
integer,dimension(200)                     :: connode            ! a list of nodes connected to a given node
integer,dimension(8)                       :: concell_nnode      ! the # of connected nodes in a concell
integer,dimension(8)                       :: concell_nbase      ! start of nodes from a cell in the list
logical,dimension(8)                       :: concell_is_regular ! true if cell contains no intersection nodes
integer,dimension(200)                     :: trialnode          ! the unique list of nodes connected to a given node


integer                                    :: ntree              ! # of nodes in the tree
integer                                    :: newalive           ! index of node most recently turned alive
integer                                    :: n_connode,n_concell! # of connected nodes or cells
integer                                    :: n_trial            ! number of neighbours to a node

integer                                    :: regular_order = 2
integer                                    :: n_det_negative


real(kind=dp)                              :: oldtime
integer,dimension(:),pointer               :: in
logical                                    :: diag = .false.
logical                                    :: diat = .false.
logical                                    :: fullcon = .false.
integer                                    :: testnode = 66566

end module mod_prop


!***********************************************************************************************
! the following subroutine does a fast marching sweep across the region that is
! its argument. Before calling this routine the regional nodes that have a starting time
! must have their values, and their node_status must be 0 for nodes whose values are already fixed,
! or 1 for nodes in the narrow band, the rest of the nodes must have node_status -1.
!
subroutine propagate(regin,vtype)

use mod_prop
implicit none

type(Tregion),target                       :: regin
integer                                    :: vtype
integer                                    :: n,nt ! local variables
logical                                    :: first_step

reg => regin
grid => reg%grid
velocity => reg%velocity(:,vtype)
velgrid => vgrid(reg%ivgrid,vtype)

n_det_negative = 0

!print *,'entering propagate'
!print *,'region',reg%id,reg%ivgrid
!print *,'top', reg%itop%id,reg%itop%iface_id
!print *,'bot', reg%ibot%id,reg%ibot%iface_id

if (reg%id == 99) then
   diat=.true.
   testnode=pgrid%rnode_id(37,9,5)
   write (24,*) 'testnode set to',testnode
   write (24,'(a20,3f12.5)') 'test coordinates',reg%r(testnode),reg%lat(testnode),reg%long(testnode)
   write (24,'(3i8)') reg%node(testnode)%i1,reg%node(testnode)%i2,reg%node(testnode)%i3
else
   diat=.false.
endif

! the time gradient derived during regular updates when fewer than 3 nodes are used can be very inaccurate
! set this flag in that case and correct afterwards

allocate(suspect_time_gradient(reg%nnode))
suspect_time_gradient = .false.

! precalculate required sines and cosines

allocate(cosla(reg%nnode),sinla(reg%nnode),coslo(reg%nnode),sinlo(reg%nnode))

cosla=cos(reg%lat)
sinla=sin(reg%lat)
coslo=cos(reg%long)
sinlo=sin(reg%long)


! set up the initial tree

allocate(node_from_tree_ind(reg%nnode),node_is_counted(reg%nnode))
node_is_counted=0
ntree = 0

do n=1,reg%nnode
   if (reg%node_status(n) == 1) call add_node_to_tree(n) 
end do


!start main propagation loop

first_step=.true.
do while (ntree > 0)

! take the node with the smallest time from the narrowband, set it to alive and adjust the tree

   newalive=node_from_tree_ind(1)
   reg%node_status(newalive) = 0
   call remove_root_from_tree

   if (ntree == 0 .and. .not.first_step) exit  ! stop when the tree structure is empty, but not in the first step since there may
                                               ! be only a single starting point (the source)

   first_step=.false.

   if (diag) print *,'root removed',newalive
   if (diag) print '(i5,3f12.5,3i5)',newalive,reg%r(newalive),reg%lat(newalive),reg%long(newalive),reg%node(newalive)

! find the neighbours of newalive, a list of node numbers is returned in the array connode 

   call find_connected_nodes(newalive)

   if (diag) print *,'neighbour list found',n_connode

   if (diag) then
      do n=1,n_connode
         print '(2i5,3f12.5,3i5)',n,connode(n),reg%r(connode(n)),reg%lat(connode(n)),reg%long(connode(n)),&
              reg%node(connode(n))
      end do

     print *
     print *,'arrivaltime at newalive is',reg%arrivaltime(newalive)
     print *

     stop
   endif

! save the list of neighbours, because the array connode will be used again when finding neighbours of the
! neighbours of newalive that have to be updated 

   n_trial=n_connode
   trialnode(1:n_trial)=connode(1:n_trial)


!  loop over neighbours and take appropopriate action

   do nt=1,n_trial

! if the neighbour is already in the narrow band, update it

      if (reg%node_status(trialnode(nt)) > 0) then

         oldtime=reg%arrivaltime(trialnode(nt))
         call update_time(trialnode(nt))

         if (diag) print '(i5,a6,f12.5,a6,f12.5)',nt,'nb ',reg%arrivaltime(trialnode(nt)),'was',oldtime

         if ( reg%arrivaltime(trialnode(nt)) < oldtime ) call update_tree(trialnode(nt))

      endif

! if the neighbour is far, add it to the narrow band

      if (reg%node_status(trialnode(nt)) == -1) then

         reg%arrivaltime(trialnode(nt)) = huge_time
         oldtime=reg%arrivaltime(trialnode(nt))
         call update_time(trialnode(nt))

         if (diag) print '(i5,a6,f12.5)',nt,'far',reg%arrivaltime(trialnode(nt))

         if ( reg%arrivaltime(trialnode(nt)) < oldtime ) call add_node_to_tree(trialnode(nt))

      endif
   end do   ! loop over neighbours of newalive

   if (diag) print *


end do  ! main propagation loop


deallocate(cosla,sinla,coslo,sinlo)
deallocate(node_from_tree_ind,node_is_counted)


if (n_det_negative > 0) print *,'warning!!!! determinant in regular update was negative ',n_det_negative,' times'


! correct suspect time gradients

   do n=1,reg%nnode
      if (suspect_time_gradient(n)) then
         call fit_gradient_at_node(reg,n,reg%time_gradient(1,n),reg%time_gradient(2,n),reg%time_gradient(3,n))
      endif
   end do

deallocate(suspect_time_gradient)

!-----------------------------------------------------------------

end subroutine propagate


!-----------------------------------------------------------------------------------------
! when returning,this intrinsic subroutine has filled the arrays giving the connected cells
! and the list of unique neighbours (connected nodes) of the node centernode

  subroutine find_connected_nodes(centernode)
    use mod_prop
    implicit none
    integer                                    :: i1,i2,i3           ! identify a node (i,j,k or 0,interface,inode)
    integer                                    :: n,m,i,j,k,ii,jj,kk,icell ! local variables
    integer   :: centernode
    type(Tintersection),pointer                :: isec
    

! store the identifiers of the centernode node in local variables. Reminder:
! for regular grid nodes i1,i2,i3  correspond to ir,ilat,ilong of the node
! for intersection nodes i1,i2,i3 are 0, intersection #, node # in intersection 

   i1 = reg%node(centernode)%i1 ; i2 = reg%node(centernode)%i2 ; i3 = reg%node(centernode)%i3

!   print *,'finding neighbours of',i1,i2,i3

   if (i1 /= 0) then
      if (grid%fully_regular(i1,i2,i3)) then

         n_connode = 6
         connode(1)=grid%rnode_id(i1-1,i2,i3)
         connode(2)=grid%rnode_id(i1+1,i2,i3)
         connode(3)=grid%rnode_id(i1,i2-1,i3)
         connode(4)=grid%rnode_id(i1,i2+1,i3)
         connode(5)=grid%rnode_id(i1,i2,i3-1)
         connode(6)=grid%rnode_id(i1,i2,i3+1)

         return

      endif
   endif

! make a list of grid cells of which the new alive node is part

   if (i1 /= 0) then   ! centernode is a regular grid node. Use spatial correlation of grid

      n_concell = 0
      do i=0,1
         if ((i1-i>0) .and. (i1-i<grid%nr)) then
            do j=0,1
               if ((i2-j>0) .and. (i2-j<grid%nlat)) then
                  do k=0,1
                     if ((i3-k>0) .and. (i3-k<grid%nlong)) then
                        n_concell=n_concell+1
                        concell(n_concell)%ir=i1-i  
                        concell(n_concell)%ilat=i2-j 
                        concell(n_concell)%ilong=i3-k
                     endif
                  end do
               end if
            end do
         end if
      end do

   else          ! centernode is an intersection node. Use the connections that have been found before

      if (i2 == reg%itop%id) then ; isec => reg%itop ; else ; isec => reg%ibot ; endif
      n_concell = 0
      do n=1,8
         icell = isec%ccell_from_inode(n,i3)
         if (icell > 0) then
            n_concell=n_concell+1
            concell(n_concell)%ir= isec%ccells(icell)%ir
            concell(n_concell)%ilat=isec%ccells(icell)%ilat 
            concell(n_concell)%ilong=isec%ccells(icell)%ilong
         endif
      enddo

   endif

!print *,'connected cells'
!do n=1,n_concell
!   print *,n,concell(n)
!end do


! make a list of nodes in the connected cells

! make sure the centernode itself is not counted as one of the neighbours

   node_is_counted(centernode)=centernode

   n_connode=0

   do n=1,n_concell

      ! if centernode is a regular node (i1 /= 0) , the cell is potentially regular

      concell_is_regular(n) = (i1 /= 0)

! find the intersection nodes of connected cell n

!   explanation: each intersection has a 1-D list of cells cut by the interface, and a list of intersection
!   nodes that are part of each cut cell. Each regular grid cell has a pointer ccind_from_3dc(i,j,k)%p
!   (Cut Cell INDex FROM 3D Coordinates)
!   associated with it, where p is a pointer to an integer array with as many elements as there are intersections.
!   If a cell is cut by interface n, the pointer ccind_from_3dc(i,j,k)%p is allocated, and the variable
!   ccind_from_3dc(i,j,k)%p(n) contains the array index of cell (i,j,k) in the 1D cut cell list of intersection n 

! test whether the cell is cut by any interface, in that case the pointer to the local list of
! interfaces cutting the cell has been allocated

      if (associated(grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p)) then

! if so, check if the cell is cut by the top intersection

     ! icell is the index of the current connected cell in the list of cells cut by interface reg%itop.
                  
         icell = grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%itop%iface_id)

     ! if icell == 0 the cell is not cut be the top interface
    
         if(icell /= 0) then
 
            concell_is_regular(n) = .false.

            isec => reg%itop
            do jj=1,isec%n_inodes(icell)

!   m is the node number in the regional node list of node  jj in the list of inteface nodes that are part of cut cell icell

                m=isec%rbel_node_id(isec%inodes(jj,icell))

               if ( node_is_counted(m) /= centernode ) then   ! if node is not yet counted as a neighbour of centernode

                  n_connode=n_connode+1
                  connode(n_connode)=m
                  node_is_counted(m) = centernode

               endif

            end do

         end if


! then check the bottom intersection

         icell = grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%ibot%iface_id)
         if (icell /= 0) then

            concell_is_regular(n) = .false.

            isec => reg%ibot
            do jj=1,isec%n_inodes(icell)

               m=isec%rabo_node_id(isec%inodes(jj,icell))

               if ( node_is_counted(m) /= centernode ) then 
                  n_connode=n_connode+1
                  connode(n_connode)=m
                  node_is_counted(m) = centernode

               endif

            end do

         endif

      endif


! then find the regular grid nodes of connected cell n

      do i=0,1
         ii=concell(n)%ir+i
         do j=0,1
            jj=concell(n)%ilat+j
            do k=0,1
               kk=concell(n)%ilong+k

               ! reduced connectivity for regular nodes if cell is completely regular
               if (((i1 /= 0 .and. abs(i1-ii)+abs(i2-jj)+abs(i3-kk) == 1)).or.(.not.concell_is_regular(n))) then

                 if (grid%node_region(ii,jj,kk) == reg%id) then  ! node has to belong to the current region

                     m=grid%rnode_id(ii,jj,kk)
                     if ( node_is_counted(m) /= centernode ) then 
                        n_connode=n_connode+1
                        connode(n_connode)=m
                        node_is_counted(m) = centernode
                     endif

                  endif

               endif

            end do
         end do
      end do

   end do  ! loop over connected cells


   return

  end subroutine find_connected_nodes


!-----------------------------------------------------------------------------------------
! when returning,this intrinsic subroutine has filled the arrays giving the connected cells
! of the node centernode and the ALIVE nodes that are part of these connected cells (non-unique,
! in the sense that nodes that are part of more than one connected cell appear more than one 
! time

  subroutine find_connected_cells(centernode)
    use mod_prop
    implicit none
    integer                                    :: i1,i2,i3           ! identify a node (i,j,k or 0,interface,inode)
    integer                                    :: n,m,i,j,k,ii,jj,kk,icell ! local variables
    integer  :: centernode
    type(Tintersection),pointer                :: isec

! store the identifiers of the centernode node in local variables. Reminder:
! for regular grid nodes i1,i2,i3  correspond to ir,ilat,ilong of the node
! for intersection nodes i1,i2,i3 are 0, intersection #, node # in intersection 

   i1 = reg%node(centernode)%i1 ; i2 = reg%node(centernode)%i2 ; i3 = reg%node(centernode)%i3


! make a list of grid cells of which the new node is part

   if (i1 /= 0) then   ! centernode is a regular grid node. Use spatial correlation of grid

      n_concell = 0
      do i=0,1
         if ((i1-i>0) .and. (i1-i<grid%nr)) then
            do j=0,1
               if ((i2-j>0) .and. (i2-j<grid%nlat)) then
                  do k=0,1
                     if ((i3-k>0) .and. (i3-k<grid%nlong)) then
                        n_concell=n_concell+1
                        concell(n_concell)%ir=i1-i  
                        concell(n_concell)%ilat=i2-j 
                        concell(n_concell)%ilong=i3-k
                     endif
                  end do
               end if
            end do
         end if
      end do

   else          ! centernode is an intersection node. Use the connections that have been found before

      if (i2 == reg%itop%id) then ; isec => reg%itop ; else ; isec => reg%ibot ; endif
      n_concell = 0
      do n=1,8
         icell = isec%ccell_from_inode(n,i3)
         if (icell > 0) then
            n_concell=n_concell+1
            concell(n_concell)%ir= isec%ccells(icell)%ir
            concell(n_concell)%ilat=isec%ccells(icell)%ilat 
            concell(n_concell)%ilong=isec%ccells(icell)%ilong
         endif
      enddo

   endif

! make a list of nodes in the connected cells

   n_connode=0

   do n=1,n_concell

     ! if centernode is a regular node (i1 /= 0) , the cell is potentially regular

     concell_is_regular(n) = (i1 /= 0)

! the intersection nodes of these cells

! test whether the cell is cut by any interface, in that case the pointer to the local list of
! interfaces cutting the cell has been allocated

      if (associated(grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p)) then

! if so, check the top intersection

     ! icell is the index of the current connected cell in the list of cells cut by interface reg%itop
     ! if the current connected cell is not cut by this interface icell == 0

         icell = grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%itop%iface_id)
         if(icell /= 0) then

            concell_is_regular(n)=.false.

            isec => reg%itop
            do jj=1,isec%n_inodes(icell)

               m=isec%rbel_node_id(isec%inodes(jj,icell))

               ! take only alive nodes (note that this excludes centernode itself)
               if ( reg%node_status(m) == 0 ) then
                  n_connode=n_connode+1
                  connode(n_connode)=m
               endif

            end do

         end if


! then do the same for the bottom intersection

         icell = grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%ibot%iface_id)
         if (icell /= 0) then

            concell_is_regular(n)=.false.


            isec => reg%ibot

            if (icell < 1) print *, 'icell < 1',icell
            if (icell > isec%n_ccells) then
               print *,'icell > n_ccells',icell,isec%n_ccells
               print *,'cell',concell(n)%ir,concell(n)%ilat,concell(n)%ilong
               print *,'interface',ii
               print *,grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(1:4)
            endif

            do jj=1,isec%n_inodes(icell)

               m=isec%rabo_node_id(isec%inodes(jj,icell))
               if ( reg%node_status(m) == 0 ) then
                  n_connode=n_connode+1
                  connode(n_connode)=m
               endif

            end do

         endif

      endif


! the regular grid nodes of these cells

      do i=0,1
         ii=concell(n)%ir+i
         do j=0,1
            jj=concell(n)%ilat+j
            do k=0,1
               kk=concell(n)%ilong+k

               ! reduced connectivity for regular nodes if cell is completely regular
               if ((i1 /= 0 .and. abs(i1-ii)+abs(i2-jj)+abs(i3-kk) == 1).or.(.not.concell_is_regular(n))) then

                  if (grid%node_region(ii,jj,kk) == reg%id) then  ! node has to belong to the current region

                     m=grid%rnode_id(ii,jj,kk)
                     if ( reg%node_status(m) == 0 ) then  ! node has to be alive(note that this excludes centernode itself)
                        n_connode=n_connode+1
                        connode(n_connode)=m
                     endif

                  endif

               endif

            end do
         end do
      end do

      concell_nnode(n)=n_connode  

   end do  ! loop over connected cells


! store the starting index and length of the list of alive nodes for each cell, since updating is
! done by cell
 
   concell_nbase(1)=1
   do n=1,n_concell

      concell_nnode(n)=concell_nnode(n)-concell_nbase(n)+1
      if (n < n_concell) concell_nbase(n+1)=concell_nbase(n)+concell_nnode(n)

   end do

   return

  end subroutine find_connected_cells

!-------------------------------------------------------------------------
!
!
  subroutine update_time(update_node)
    use mod_prop
    implicit none

    integer       :: update_node
    integer       :: n,i,j,k            ! local variables
    integer       :: jl,kl, i1,i2,i3,id
    integer       :: anode(100),nanode  ! list of nodes used for the update


   i1 = reg%node(update_node)%i1 ; i2 = reg%node(update_node)%i2 ; i3 = reg%node(update_node)%i3

   if (i1 /= 0) then

      if (grid%fully_regular(i1,i2,i3)) then

!         print *,'node is fully regular'     ! deldiag

         do k=-1,1,2
            do j=-1,1,2
               do i=-1,1,2

                  nanode=0

                  id = grid%rnode_id(i1+i,i2,i3)
                  if (reg%node_status(id) == 0) then
                     nanode=nanode+1
                     anode(nanode)= id
                  endif
                  id = grid%rnode_id(i1,i2+j,i3)
                  if (reg%node_status(id) == 0) then
                     nanode=nanode+1
                     anode(nanode)= id
                  endif
                  id = grid%rnode_id(i1,i2,i3+k)
                  if (reg%node_status(id) == 0) then
                     nanode=nanode+1
                     anode(nanode)= id
                  endif

                  if (nanode > 0)  then

                  ! put newalive at position 1 if it is there
                     do n=1,nanode
                        if (anode(n) == newalive) then
                           anode(n)=anode(1)
                           anode(1)=newalive
                        endif
                     end do

                     if (anode(1) == newalive) call regular_update2(update_node,anode(1:nanode),nanode)

                  endif

               end do
            end do
         end do

         return

      endif

   endif

!
! first find the cells to which the node to be updated is connected, and the list of
! ALIVE nodes in these cells

    call find_connected_cells(update_node)

! now loop over the connected cells, and if possible use the alive nodes
! in each cell to find a new arrival time at the node to be updated
 
    do n=1,n_concell


       if (concell_nnode(n) > 0) then   ! if there are any alive nodes in this cell


! copy alive nodes in this connected cell to local array
       nanode=concell_nnode(n)

       if (nanode > 30) then
          print *,'nanode =',nanode

          do i=1,nanode

             print'(4i5,4f12.5)',i,reg%node(anode(i))%i1,reg%node(anode(i))%i2,reg%node(anode(i))%i3, &
                  reg%r(anode(i)),reg%lat(anode(i)),reg%long(anode(i))

          enddo
          stop 'subroutine update_time: nanode >30'
       endif

       anode(1:nanode)=connode(concell_nbase(n):concell_nbase(n)+nanode-1)


! put newalive at position 1 if it is there
       do i=1,nanode
          if (anode(i) == newalive) then
             j=anode(1)
             anode(1)=anode(i)
             anode(i)=j
          endif
       end do

       if (diat) then
          if (update_node == testnode) then
             write(24,*) '************************'
             write(24,'(2i8,a15,3i5)') update_node,newalive,'connected cell',concell(n)
             write(24,'(10i12)')anode(1:nanode)
             write(24,'(10i12)')reg%node_status(anode(1:nanode))
             write(24,'(10i12)')reg%node(anode(1:nanode))%i1
             write(24,'(10i12)')reg%node(anode(1:nanode))%i2
             write(24,'(10i12)')reg%node(anode(1:nanode))%i3
             write(24,'(10f12.5)')reg%r(anode(1:nanode))
             write(24,'(10f12.5)')reg%lat(anode(1:nanode))
             write(24,'(10f12.5)')reg%long(anode(1:nanode))
             write(24,'(10f12.5)')reg%arrivaltime(anode(1:nanode))
          endif
       endif


    ! test if newalive is among the nodes that will be used for the update
    ! otherwise the update has already been done and does not have to be repeated

          if (anode(1) == newalive) then 

             if (concell_is_regular(n)) then

                ! do a regular fast marching update, second order if possible, in this octant

                call regular_update2(update_node,anode(1:nanode),nanode)

             else

             ! do an irregular update in this octant

                select case (nanode)
          
                case(1)     ! there is 1 alive node in the cell

                   call time_from_1_node(update_node,anode(1))

                case(2)     ! there are 2 alive nodes in the cell

                   call time_from_1_node(update_node,anode(1))
                   call time_from_2_nodes(update_node,anode(1),anode(2))

                case(3)     ! there are 3 alive nodes in the cell

                   call time_from_1_node(update_node,anode(1))
                   call time_from_2_nodes(update_node,anode(1),anode(2))
                   call time_from_2_nodes(update_node,anode(1),anode(3))
                   call time_from_3_nodes(update_node,anode(1),anode(2),anode(3))

                case(4:)     ! there are more than 3 alive nodes in the cell

                   call time_from_1_node(update_node,anode(1))

                   do jl=2,nanode
                      call time_from_2_nodes(update_node,anode(1),anode(jl))
                   end do

                   do jl=2,nanode-1
                      do kl=jl+1,nanode
                         call time_from_3_nodes(update_node,anode(1),anode(jl),anode(kl))
                      end do
                   end do


                end select

             endif   ! regular or irregular

          endif ! if newalive is among the nodes used for the update

       endif   ! if there are any alive nodes in this cell


    end do ! loop over connected cells


  end subroutine update_time

!***********************************************************************************************
subroutine time_from_1_node(unode,node1)
  use mod_prop
  implicit none

  integer           :: unode,node1
  real(kind=dp)     :: atime,ttime,dist,xu(3),x1(3) 
  real(kind=dp)     :: a(3,3),b(3)

!  print *,'updating from 1 node'     ! deldiag


! for all irregular updates the coordinates of the nodes are first converted to Cartesian
! afterwards the time gradient( wave vector of updating wave) is converted back to spherical

  xu(1)=reg%r(unode)*cosla(unode)*coslo(unode)
  xu(2)=reg%r(unode)*cosla(unode)*sinlo(unode)
  xu(3)=reg%r(unode)*sinla(unode)

  x1(1)=reg%r(node1)*cosla(node1)*coslo(node1)
  x1(2)=reg%r(node1)*cosla(node1)*sinlo(node1)
  x1(3)=reg%r(node1)*sinla(node1)

  dist=sqrt(sum((xu-x1)**2))

  if (dist < 0.01_dp*pgrid%tolerance) then
     reg%arrivaltime(unode) = reg%arrivaltime(node1)
     reg%time_gradient(1:3,unode) =  reg%time_gradient(1:3,node1) 
     return
   endif


  ttime= dist/(0.5*(velocity(unode)+velocity(node1)))
  atime= reg%arrivaltime(node1) + ttime

  if ( atime < reg%arrivaltime(unode)) then

     reg%arrivaltime(unode) = atime

! convert time gradient back to spherical

     b=(xu-x1)/(dist*velocity(unode))
     a(1,1)=cosla(unode)*coslo(unode) ; a(2,1)=-sinla(unode)*coslo(unode) ; a(3,1)=-sinlo(unode)
     a(1,2)=cosla(unode)*sinlo(unode) ; a(2,2)=-sinla(unode)*sinlo(unode) ; a(3,2)= coslo(unode)
     a(1,3)=sinla(unode)              ; a(2,3)= cosla(unode)              ; a(3,3)= 0.0_dp
     reg%time_gradient(1:3,unode) = matmul(a,b)

     call refract_locally(reg%r(unode),reg%lat(unode),reg%long(unode),velgrid, &
          reg%time_gradient(1:3,unode))

     if (diat.and.unode == testnode) write (24,*) 'arrival adjusted'

  endif

  if (diat.and.unode == testnode) write(24,'(a20,2i8,2f12.5)') 'updating from 1 node',newalive,node1,&
       reg%arrivaltime(unode),reg%arrivaltime(node1)

  return

end subroutine time_from_1_node


!**********************************************************************************************
subroutine time_from_2_nodes(unode,anode1,anode2)
  use mod_prop
  implicit none

  integer           :: unode,node1,node2,anode1,anode2
  real(kind=dp)     :: atime  

  real(kind=dp)     :: wn1(3),wn2(3)
  real(kind=dp)     :: a(3,3),b(3),p(3),q(3),u1,u2,det,ae,be,ce
  real(kind=dp)     :: n1(3),n2(3),r1(3),r2(3),n0(3),xu(3),x1(3),x2(3)

  integer           :: i,vscount,indx(3),status


! we take a geometric approach to solving for the wavefront normal n. The 3 equations
! to be fulfilled are:
! t2 = t1+n.(x2-x1)/v12        --> a(1,1)*n1 + a(1,2)*n2 + a(1,3)*n3  = b(1)
!  0 = n . [(x2-x) x (x1-x)]   --> a(2,1)*n1 + a(2,2)*n2 + a(2,3)*n3  = b(2)
! |n| = 0                      --> n1**2 +   n2**2 +  n3**2   = 1
!
! The second equation comes from the condition that the wavefront normal n lies in the plane
! containing the 2 known points and the point at which the arrivaltime is evaluated.
! The method of solution first calculates the line defining the intersection of
! the two planes defined by the two linear equations, then the 2 intersections of this line
! with  the sphere defined by the normalization constraint


  if (diat.and.unode == testnode) write(24,*) 'updating from 2 nodes'

  node1=anode1 ; node2=anode2 

! order from latest to earliest arrival time

  if (reg%arrivaltime(node1) < reg%arrivaltime(node2)) then ; i=node1 ; node1=node2 ; node2=i ; endif


! convert nodes to Cartesian

  xu(1)=reg%r(unode)*cosla(unode)*coslo(unode)
  xu(2)=reg%r(unode)*cosla(unode)*sinlo(unode)
  xu(3)=reg%r(unode)*sinla(unode)

  x1(1)=reg%r(node1)*cosla(node1)*coslo(node1)
  x1(2)=reg%r(node1)*cosla(node1)*sinlo(node1)
  x1(3)=reg%r(node1)*sinla(node1)

  x2(1)=reg%r(node2)*cosla(node2)*coslo(node2)
  x2(2)=reg%r(node2)*cosla(node2)*sinlo(node2)
  x2(3)=reg%r(node2)*sinla(node2)


! first define some local variables (r1,r2 difference vectors between node1,node2 and unode)

  r1=x1-xu
  r2=x2-xu


! first check if nodes coincide, if so copy times and gradients and return

  if (sqrt(sum(r1**2)) < 0.01_dp*pgrid%tolerance) then
     reg%arrivaltime(unode) = reg%arrivaltime(node1)
     reg%time_gradient(1:3,unode) =  reg%time_gradient(1:3,node1) 
     return
  endif
  if (sqrt(sum(r2**2)) < 0.01_dp*pgrid%tolerance) then
     reg%arrivaltime(unode) = reg%arrivaltime(node2)
     reg%time_gradient(1:3,unode) =  reg%time_gradient(1:3,node2) 
     return
  endif
  if (sqrt(sum((x2-x1)**2)) < 0.01_dp*pgrid%tolerance) then
     call time_from_1_node(unode,node1)
     return
  endif

 ! construct the two linear equations
  b(1)=reg%arrivaltime(node1)-reg%arrivaltime(node2)
  b(2)=0.0_dp

!  a(1,1:3)=(x1-x2)/velocity(node1)
  a(1,1:3)=(x1-x2)/(0.5_dp*(velocity(node1)+velocity(node2)))


  a(2,1)=r1(2)*r2(3)-r1(3)*r2(2)
  a(2,2)=r1(3)*r2(1)-r1(1)*r2(3)
  a(2,3)=r1(1)*r2(2)-r1(2)*r2(1)

  n0=a(2,:)    ! save the normal to the plane containing the 3 points (2 alive+update) for future use


! the direction of the line that is the intersection between the two planes defined by the
! two linear equations is the cross product  q = a(1,1:3) x a(2,1:3)

  q(1)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  q(2)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  q(3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)


! a third equation is added that defines the plane that is normal to the direction of the 
! intersection line and goes through (0,0,0)

  b(3)=0.0_dp
  a(3,1)=q(1) ; a(3,2)=q(2) ; a(3,3)=q(3)


! now find the intersection point of the intersection line with the plane that is normal to it
! and goes through (0,0,0). This gives the minimum length support vector

  call ludcmp(a,indx,status)
  if (status == 1) then

     print *,a(1,1:3),b(1)
     print *,a(2,1:3),b(2)
     print *,a(3,1:3),b(3)

     stop 'singular matrix in time_from_2_nodes'
  endif
  call lubksb(a,indx,b)
  p=b

! we now have the intersection line defined as p+uq , where u is a scalar variable
! solve the quadratic normalization constraint |n| = 1 for u

   ae=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
   be=2.0_dp*(p(1)*q(1)+p(2)*q(2)+p(3)*q(3))
   ce=p(1)*p(1)+p(2)*p(2)+p(3)*p(3)-1.0_dp

   det=be*be-4.0_dp*ae*ce

   if (det < 0.0_dp) return   ! there is no valid solution if the line does not intersect the sphere


   u1=0.5_dp*(-be+sqrt(det))/ae
   u2=0.5_dp*(-be-sqrt(det))/ae

! the two solutions for the wave vector

   wn1=p+u1*q
   wn2=p+u2*q

  if (diat.and.unode == testnode) write(24,*)' wn1',wn1
  if (diat.and.unode == testnode) write(24,*)' wn2',wn2

! now test if either of these two solutions is valid (obeys the causality constraint)
! -wn  must lie within the wedge defined by the planes containing x1-x0 or x2-x0
! and the normal to the plane of the three points x0,x1,x2, which is just the vector a(2,:)
! calculated above (saved in n0 before calling ludcmp because a is modified there)


! calculate the inward pointing normals of the two planes defining the wedge

   n1(1)=n0(2)*r1(3)-n0(3)*r1(2)
   n1(2)=n0(3)*r1(1)-n0(1)*r1(3)
   n1(3)=n0(1)*r1(2)-n0(2)*r1(1)
   n1=n1/sqrt(sum(n1**2))
   if (n1(1)*r2(1)+n1(2)*r2(2)+n1(3)*r2(3) < 0.0_dp) n1=-n1   ! define the sign of the normal so the other alive point is
                                                              ! on the positive side 
   n2(1)=n0(2)*r2(3)-n0(3)*r2(2)
   n2(2)=n0(3)*r2(1)-n0(1)*r2(3)
   n2(3)=n0(1)*r2(2)-n0(2)*r2(1)
   n2=n2/sqrt(sum(n2**2))
   if (n2(1)*r1(1)+n2(2)*r1(2)+n2(3)*r1(3) < 0.0_dp) n2=-n2   ! define the sign of the normal so the other alive point is
                                                              ! on the positive side 


  if (diat.and.unode == testnode) write(24,*)' n1',n1
  if (diat.and.unode == testnode) write(24,*)' n2',n2

  if (diat.and.unode == testnode) write(24,*)'dp1 ',dot_product(wn1,n1) , dot_product(wn1,n2)
  if (diat.and.unode == testnode) write(24,*)'dp2 ',dot_product(wn2,n1) , dot_product(wn2,n2)

   vscount=0
   if (dot_product(wn1,n1) <= 0.0_dp .and. dot_product(wn1,n2) <= 0.0_dp ) then

      atime= reg%arrivaltime(node1)-dot_product(wn1,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
      if (atime >= reg%arrivaltime(node1)) then
         if ( atime < reg%arrivaltime(unode)) then
            reg%arrivaltime(unode) = atime

            b=wn1/(sqrt(sum(wn1**2))*velocity(node1))
            a(1,1)=cosla(unode)*coslo(unode) ; a(2,1)=-sinla(unode)*coslo(unode) ; a(3,1)=-sinlo(unode)
            a(1,2)=cosla(unode)*sinlo(unode) ; a(2,2)=-sinla(unode)*sinlo(unode) ; a(3,2)= coslo(unode)
            a(1,3)=sinla(unode)              ; a(2,3)= cosla(unode)              ; a(3,3)= 0.0_dp
            reg%time_gradient(1:3,unode) = matmul(a,b)

            x1=reg%time_gradient(1:3,unode)

            call refract_locally(reg%r(unode),reg%lat(unode),reg%long(unode),velgrid, &
                 reg%time_gradient(1:3,unode))

            if (diat.and.unode == testnode) then
               write (24,*) 'arrival adjusted'
               write(24,'(a5,3f12.6)') 'wn1',wn1
               write(24,'(a5,3f12.6)') 'b  ',b
               write(24,'(a5,3f12.6)') 'tg ',x1
               write(24,'(a5,3f12.6)') 'tgr',reg%time_gradient(1:3,unode)

            endif
         endif
         vscount=vscount+1
   if (diat.and.unode == testnode) write (24,'(a20,3i8,2f12.5)') 'update from 2 node',newalive,node1,node2,&
        reg%arrivaltime(unode),reg%arrivaltime(node1)
      endif

   end if

   if (dot_product(wn2,n1) <= 0.0_dp .and. dot_product(wn2,n2) <= 0.0_dp ) then

      atime= reg%arrivaltime(node1)-dot_product(wn2,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
      if (atime >= reg%arrivaltime(node1)) then
         if ( atime < reg%arrivaltime(unode)) then
            reg%arrivaltime(unode) = atime

            b=wn2/(sqrt(sum(wn2**2))*velocity(node1))
            a(1,1)=cosla(unode)*coslo(unode) ; a(2,1)=-sinla(unode)*coslo(unode) ; a(3,1)=-sinlo(unode)
            a(1,2)=cosla(unode)*sinlo(unode) ; a(2,2)=-sinla(unode)*sinlo(unode) ; a(3,2)= coslo(unode)
            a(1,3)=sinla(unode)              ; a(2,3)= cosla(unode)              ; a(3,3)= 0.0_dp
            reg%time_gradient(1:3,unode) = matmul(a,b)

            call refract_locally(reg%r(unode),reg%lat(unode),reg%long(unode),velgrid, &
                 reg%time_gradient(1:3,unode))

            if (diat.and.unode == testnode) write (24,*) 'arrival adjusted'
         endif
         vscount=vscount+1
   if (diat.and.unode == testnode) write (24,'(a20,3i8,2f12.5)') 'update from 2 node',newalive,node1,node2, &
        reg%arrivaltime(unode),reg%arrivaltime(node1)
      endif
   end if

   if (vscount == 2 .and. u1 /= u2 ) then
      if (dot_product(wn1,wn2) < 0.99999_dp) then
      print *,' r1',r1
      print *,' wn1',wn1
      print *,' wn2',wn2
      print *,'dp1 ',dot_product(wn1,n1) , dot_product(wn1,n2)
      print *,'dp2 ',dot_product(wn2,n1) , dot_product(wn2,n2)
      stop 'time_from_2_nodes: two different valid solutions found'
      endif
   endif
   return

end subroutine time_from_2_nodes

!**********************************************************************************************
subroutine time_from_3_nodes(unode,anode1,anode2,anode3)
  use mod_prop
  implicit none

  integer           :: unode,node1,node2,node3,anode1,anode2,anode3
  real(kind=dp)     :: atime  


  real(kind=dp)     :: wn1(3),wn2(3)
  real(kind=dp)     :: a(3,3),b(3),p(3),q(3),u1,u2,det,ae,be,ce,dist_from_plane,cangle(3)
  real(kind=dp)     :: n1(3),n2(3),n3(3),r1(3),r2(3),r3(3),xu(3),x1(3),x2(3),x3(3)

  integer           :: i,vscount,indx(3),status,min_loc(1)

  logical           :: wn1_inplane,wn2_inplane


! we take a geometric approach to solving for the wavefront normal n . The 3 equations
! to be fulfilled are:
! t2=t1+n.(x2-x1)/v12   --> a(1,1)*n1 + a(1,2)*n2 + a(1,3)*n3  = b(1)
! t3=t1+n.(x3-x1)/v13   --> a(2,1)*n1 + a(2,2)*n2 + a(2,3)*n3  = b(2)
! |n| = 0               --> n1**2 +   n2**2 +  n3**2   = 1
! the method of solution first calculates the line defining the intersection of
! the two planes defined by the two linear equations, then the intersections of this line
! with  the sphere defined by the normalization constraint


  if (diat.and.unode == testnode) write(24,*) 'updating from 3 nodes',node1,node2,node3

  node1=anode1 ; node2=anode2 ; node3=anode3

! order from latest to earliest arrival time

  if (reg%arrivaltime(node1) < reg%arrivaltime(node2)) then ; i=node1 ; node1=node2 ; node2=i ; endif
  if (reg%arrivaltime(node2) < reg%arrivaltime(node3)) then ; i=node2 ; node2=node3 ; node3=i ; endif
  if (reg%arrivaltime(node1) < reg%arrivaltime(node2)) then ; i=node1 ; node1=node2 ; node2=i ; endif


! construct the two linear equations
  b(1)=reg%arrivaltime(node1)-reg%arrivaltime(node2)
  b(2)=reg%arrivaltime(node1)-reg%arrivaltime(node3)

! convert to Cartesian

  xu(1)=reg%r(unode)*cosla(unode)*coslo(unode)
  xu(2)=reg%r(unode)*cosla(unode)*sinlo(unode)
  xu(3)=reg%r(unode)*sinla(unode)

  x1(1)=reg%r(node1)*cosla(node1)*coslo(node1)
  x1(2)=reg%r(node1)*cosla(node1)*sinlo(node1)
  x1(3)=reg%r(node1)*sinla(node1)

  x2(1)=reg%r(node2)*cosla(node2)*coslo(node2)
  x2(2)=reg%r(node2)*cosla(node2)*sinlo(node2)
  x2(3)=reg%r(node2)*sinla(node2)

  x3(1)=reg%r(node3)*cosla(node3)*coslo(node3)
  x3(2)=reg%r(node3)*cosla(node3)*sinlo(node3)
  x3(3)=reg%r(node3)*sinla(node3)


   r1=x1-xu
   r2=x2-xu
   r3=x3-xu


! first check if nodes coincide, if so copy times and gradients and return

  if (sqrt(sum(r1**2)) < 0.01_dp*pgrid%tolerance) then
      atime= reg%arrivaltime(node1)
      reg%arrivaltime(unode) = reg%arrivaltime(node1)
      reg%time_gradient(1:3,unode) =  reg%time_gradient(1:3,node1) 
      return
   endif
  if (sqrt(sum(r2**2)) < 0.01_dp*pgrid%tolerance) then
     reg%arrivaltime(unode) = reg%arrivaltime(node2)
     reg%time_gradient(1:3,unode) =  reg%time_gradient(1:3,node2) 
     return
   endif     
  if (sqrt(sum(r3**2)) < 0.01_dp*pgrid%tolerance) then
     reg%arrivaltime(unode) = reg%arrivaltime(node3)
     reg%time_gradient(1:3,unode) =  reg%time_gradient(1:3,node3) 
     return
   endif    
  if (sqrt(sum((x2-x1)**2)) < 0.1_dp*pgrid%tolerance) then
     call time_from_2_nodes(unode,node1,node3)
     return
  endif
  if (sqrt(sum((x2-x3)**2)) < 0.1_dp*pgrid%tolerance) then
     call time_from_2_nodes(unode,node1,node2)
     return
  endif
  if (sqrt(sum((x3-x1)**2)) < 0.1_dp*pgrid%tolerance) then
     call time_from_2_nodes(unode,node1,node2)
     return
  endif

!  a(1,1:3)=(x1-x2)/velocity(node1)
!  a(2,1:3)=(x1-x3)/velocity(node1)
  a(1,1:3)=(x1-x2)/(0.5_dp*(velocity(node1)+velocity(node2)))
  a(2,1:3)=(x1-x3)/(0.5_dp*(velocity(node1)+velocity(node3)))

! the direction of the line that is the intersection between the two planes defined by the
! two linear equations is the cross product  q = a(1:3) x b(1:3)

  q(1)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  q(2)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  q(3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)



! a third equation is added that defines he plane that is normal to the direction of the 
! intersection line and goes through (0,0,0)

  b(3)=0.0_dp
  a(3,1)=q(1) ; a(3,2)=q(2) ; a(3,3)=q(3)


! now find the intersection point of the intersection line with the plane

  call ludcmp(a,indx,status)
  if (status == 1) stop 'singular matrix in time_from_3_nodes'
  call lubksb(a,indx,b)
  p=b


! we now have the intersection line defined as p+uq , where u is a scalar variable
! solve the quadratic normalization constraint |n| = 1 for u

   ae=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
   be=2.0_dp*(p(1)*q(1)+p(2)*q(2)+p(3)*q(3))
   ce=p(1)*p(1)+p(2)*p(2)+p(3)*p(3)-1.0_dp

   det=be*be-4.0_dp*ae*ce

   if (det < 0.0_dp) return  ! if the line does not intersect the sphere there is no valid solution

   u1=0.5_dp*(-be+sqrt(det))/ae
   u2=0.5_dp*(-be-sqrt(det))/ae


! the two solutions for the wavefront normal

   wn1=p+u1*q
   wn2=p+u2*q

   if (diat.and.unode == testnode) write(24,*) ' wn1',wn1
   if (diat.and.unode == testnode) write(24,*) ' wn2',wn2

! now test if either of these two solutions is valid (obeys the causality constraint)
! -wn  must lie within the sector defined by the vectors x1-x0,x2-x0,x3-x0



   if (diat.and.unode == testnode) write(24,*) 'r1 ',r1
   if (diat.and.unode == testnode) write(24,*) 'r2 ',r2
   if (diat.and.unode == testnode) write(24,*) 'r3 ',r3

! calculate the inward pointing normals to the three planes defining the sector

   n1(1)=r1(2)*r2(3)-r1(3)*r2(2)
   n1(2)=r1(3)*r2(1)-r1(1)*r2(3)
   n1(3)=r1(1)*r2(2)-r1(2)*r2(1)
   n1=n1/sqrt(sum(n1**2))
 

! check for coplanarity of update and alive points, since then the causality criterion 
! becomes degenerate and breaks down
    
   dist_from_plane= (n1(1)*r3(1)+n1(2)*r3(2)+n1(3)*r3(3))/sqrt(sum(r3**2))

   if (abs(dist_from_plane) < 0.001_dp) then ! update node coplanar with alive points used for update

!      if (reg%node_status(unode) == 1) return


      wn1_inplane=abs(dot_product(wn1,n1)) < 0.001_dp
      wn2_inplane=abs(dot_product(wn2,n1)) < 0.001_dp

   ! quick test whether any of the two wave front normals lies in the nodal plane
   ! this will reject most cases with little computation

      if (.not.(wn1_inplane.or.wn2_inplane)) then

         if (diat.and.unode == testnode) write(24,*) 'coplanar reject1 ',unode,node1,node2,node3
         return

      else

         ! find the two nodes with the largest angle between them, and apply 2 point update with them

         cangle(1) = dot_product(r2,r3)/sqrt(sum(r2**2)*sum(r3**2))
         cangle(2) = dot_product(r3,r1)/sqrt(sum(r3**2)*sum(r1**2))
         cangle(3) = dot_product(r1,r2)/sqrt(sum(r1**2)*sum(r2**2))

         min_loc=minloc(cangle)

         select case(min_loc(1))

            case(1)
               n2(1)=n1(2)*r3(3)-n1(3)*r3(2)
               n2(2)=n1(3)*r3(1)-n1(1)*r3(3)
               n2(3)=n1(1)*r3(2)-n1(2)*r3(1)
               n2=n2/sqrt(sum(n2**2))
               if (n2(1)*r1(1)+n2(2)*r1(2)+n2(3)*r1(3) < 0.0_dp) n2=-n2

               n3(1)=n1(2)*r2(3)-n1(3)*r2(2)
               n3(2)=n1(3)*r2(1)-n1(1)*r2(3)
               n3(3)=n1(1)*r2(2)-n1(2)*r2(1)
               n3=n3/sqrt(sum(n3**2))
               if (n3(1)*r1(1)+n3(2)*r1(2)+n3(3)*r1(3) < 0.0_dp) n3=-n3

            case(2)
               n2(1)=n1(2)*r3(3)-n1(3)*r3(2)
               n2(2)=n1(3)*r3(1)-n1(1)*r3(3)
               n2(3)=n1(1)*r3(2)-n1(2)*r3(1)
               n2=n2/sqrt(sum(n2**2))
               if (n2(1)*r2(1)+n2(2)*r2(2)+n2(3)*r2(3) < 0.0_dp) n2=-n2

               n3(1)=n1(2)*r1(3)-n1(3)*r1(2)
               n3(2)=n1(3)*r1(1)-n1(1)*r1(3)
               n3(3)=n1(1)*r1(2)-n1(2)*r1(1)
               n3=n3/sqrt(sum(n3**2))
               if (n3(1)*r2(1)+n3(2)*r2(2)+n3(3)*r2(3) < 0.0_dp) n3=-n3

            case(3)
               n2(1)=n1(2)*r1(3)-n1(3)*r1(2)
               n2(2)=n1(3)*r1(1)-n1(1)*r1(3)
               n2(3)=n1(1)*r1(2)-n1(2)*r1(1)
               n2=n2/sqrt(sum(n2**2))
               if (n2(1)*r3(1)+n2(2)*r3(2)+n2(3)*r3(3) < 0.0_dp) n2=-n2

               n3(1)=n1(2)*r2(3)-n1(3)*r2(2)
               n3(2)=n1(3)*r2(1)-n1(1)*r2(3)
               n3(3)=n1(1)*r2(2)-n1(2)*r2(1)
               n3=n3/sqrt(sum(n3**2))
               if (n3(1)*r3(1)+n3(2)*r3(2)+n3(3)*r3(3) < 0.0_dp) n3=-n3

         end select

         vscount=0
         if (dot_product(wn1,n2) <= 0.0_dp .and. dot_product(wn1,n3) <= 0.0_dp .and. wn1_inplane) then
            
            atime= reg%arrivaltime(node1)-dot_product(wn1,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
            if (atime >= reg%arrivaltime(node1)) then
               if ( atime < reg%arrivaltime(unode)) then

                  reg%arrivaltime(unode) = atime

                  b=wn1/(sqrt(sum(wn1**2))*velocity(node1))
                  a(1,1)=cosla(unode)*coslo(unode) ; a(2,1)=-sinla(unode)*coslo(unode) ; a(3,1)=-sinlo(unode)
                  a(1,2)=cosla(unode)*sinlo(unode) ; a(2,2)=-sinla(unode)*sinlo(unode) ; a(3,2)= coslo(unode)
                  a(1,3)=sinla(unode)              ; a(2,3)= cosla(unode)              ; a(3,3)= 0.0_dp
                  reg%time_gradient(1:3,unode) = matmul(a,b)

                  call refract_locally(reg%r(unode),reg%lat(unode),reg%long(unode),velgrid, &
                       reg%time_gradient(1:3,unode))

                  if (diat.and.unode == testnode) write (24,*) 'arrival adjusted'
               endif

               vscount=vscount+1
               if(diat.and.unode == testnode) write(24,'(a15,4i8,2f12.5)') '3 node updatec',newalive,node1,node2,node3,&
                    reg%arrivaltime(unode),reg%arrivaltime(node1)
            endif

         end if

         if (dot_product(wn2,n2) <= 0.0_dp .and. dot_product(wn2,n3) <= 0.0_dp .and. wn2_inplane) then

            atime= reg%arrivaltime(node1)-dot_product(wn2,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
            if (atime >= reg%arrivaltime(node1)) then
               if ( atime < reg%arrivaltime(unode)) then

                  reg%arrivaltime(unode) = atime

                  b=wn2/(sqrt(sum(wn2**2))*velocity(node1))
                  a(1,1)=cosla(unode)*coslo(unode) ; a(2,1)=-sinla(unode)*coslo(unode) ; a(3,1)=-sinlo(unode)
                  a(1,2)=cosla(unode)*sinlo(unode) ; a(2,2)=-sinla(unode)*sinlo(unode) ; a(3,2)= coslo(unode)
                  a(1,3)=sinla(unode)              ; a(2,3)= cosla(unode)              ; a(3,3)= 0.0_dp
                  reg%time_gradient(1:3,unode) = matmul(a,b)

                  call refract_locally(reg%r(unode),reg%lat(unode),reg%long(unode),velgrid, &
                       reg%time_gradient(1:3,unode))

                  if (diat.and.unode == testnode) write (24,*) 'arrival adjusted'
               endif
               vscount=vscount+1
               if(diat.and.unode == testnode) write(24,'(a15,4i8,2f12.5)') '3 node updatec',newalive,node1,node2,node3,&
                    reg%arrivaltime(unode),reg%arrivaltime(node1)
            endif
         end if

         if (vscount == 2 .and. u1 /= u2 ) then
            if (dot_product(wn1,wn2) < 0.999_dp) then
               print *,'minloc(1)',min_loc(1)
               print *,' r1',r1
               print *,' r2',r2
               print *,' r3',r3
               print *,' wn1',wn1
               print *,' wn2',wn2
               print *,'dpt ',dot_product(wn1,n1),dot_product(wn2,n1)
               print *,'dp1 ',dot_product(wn1,n2) , dot_product(wn1,n3)
               print *,'dp2 ',dot_product(wn2,n2) , dot_product(wn2,n3)
               stop 'time_from_3_nodes: two different valid solutions found in coplanar case'
            endif
         endif

         return

      endif  ! coplanar case passing first rejection test

   endif  ! coplanar case

   if (dist_from_plane < 0.0_dp) n1=-n1                       ! define the sign of the normal so the other alive point is
                                                              ! on the positive side 

   n2(1)=r1(2)*r3(3)-r1(3)*r3(2)
   n2(2)=r1(3)*r3(1)-r1(1)*r3(3)
   n2(3)=r1(1)*r3(2)-r1(2)*r3(1)
   n2=n2/sqrt(sum(n2**2))
   if (n2(1)*r2(1)+n2(2)*r2(2)+n2(3)*r2(3) < 0.0_dp) n2=-n2   ! define the sign of the normal so the other alive point is
                                                              ! on the positive side 

   n3(1)=r2(2)*r3(3)-r2(3)*r3(2)
   n3(2)=r2(3)*r3(1)-r2(1)*r3(3)
   n3(3)=r2(1)*r3(2)-r2(2)*r3(1)
   n3=n3/sqrt(sum(n3**2))
   if (n3(1)*r1(1)+n3(2)*r1(2)+n3(3)*r1(3) < 0.0_dp) n3=-n3   ! define the sign of the normal so the other alive point is
                                                              ! on the positive side 

   if (diat.and.unode == testnode) write(24,*) 'n1 ',n1
   if (diat.and.unode == testnode) write(24,*) 'n2 ',n2
   if (diat.and.unode == testnode) write(24,*) 'n3 ',n3
   if (diat.and.unode == testnode) write(24,*) 'dp1 ',dot_product(wn1,n1) , dot_product(wn1,n2), dot_product(wn1,n3)
   if (diat.and.unode == testnode) write(24,*) 'dp2 ',dot_product(wn2,n1) , dot_product(wn2,n2), dot_product(wn2,n3)
   if (diat.and.unode == testnode) write(24,*) reg%arrivaltime(node1),dot_product(wn1,r1) , dot_product(wn2,r1)
   if (diat.and.unode == testnode) write(24,*) reg%arrivaltime(unode)
   if (diat.and.unode == testnode) write(24,*) 


! the solution is only acceptable if it lies on the correct side of all 3 planes

   vscount=0
   if (dot_product(wn1,n1) <= 0.0_dp .and. dot_product(wn1,n2) <= 0.0_dp .and. dot_product(wn1,n3) <= 0.0_dp) then

      atime= reg%arrivaltime(node1)-dot_product(wn1,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
      if (atime >= reg%arrivaltime(node1)) then
         if ( atime < reg%arrivaltime(unode)) then

            reg%arrivaltime(unode) = atime

            b=wn1/(sqrt(sum(wn1**2))*velocity(node1))
            a(1,1)=cosla(unode)*coslo(unode) ; a(2,1)=-sinla(unode)*coslo(unode) ; a(3,1)=-sinlo(unode)
            a(1,2)=cosla(unode)*sinlo(unode) ; a(2,2)=-sinla(unode)*sinlo(unode) ; a(3,2)= coslo(unode)
            a(1,3)=sinla(unode)              ; a(2,3)= cosla(unode)              ; a(3,3)= 0.0_dp
            reg%time_gradient(1:3,unode) = matmul(a,b)

            call refract_locally(reg%r(unode),reg%lat(unode),reg%long(unode),velgrid, &
                 reg%time_gradient(1:3,unode))

            if (diat.and.unode == testnode) write (24,'(a20,4f14.6)') 'arrival adjusted',atime,reg%time_gradient(1:3,unode)
         endif
         vscount=vscount+1
         if(diat.and.unode == testnode) write(24,'(a15,4i8,2f12.5)') '3 node update',newalive,node1,node2,node3, &
              reg%arrivaltime(unode),reg%arrivaltime(node1)

         if (diat.and.unode == testnode)  then
            write(24,*) 'updated from wn1'
            write(24,*) node1,node2,node3,reg%arrivaltime(unode),atime
            write(24,*) reg%node(unode)
            write(24,*) reg%node(node1)
            write(24,*) reg%node(node2)
            write(24,*) reg%node(node3)
            write(24,*)
         endif
      endif
 
   end if

   if (dot_product(wn2,n1) <= 0.0_dp .and. dot_product(wn2,n2) <= 0.0_dp .and. dot_product(wn2,n3) <= 0.0_dp) then

      atime= reg%arrivaltime(node1)-dot_product(wn2,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
      if (atime >= reg%arrivaltime(node1)) then
         if ( atime < reg%arrivaltime(unode)) then

            reg%arrivaltime(unode) = atime

            b=wn2/(sqrt(sum(wn2**2))*velocity(node1))
            a(1,1)=cosla(unode)*coslo(unode) ; a(2,1)=-sinla(unode)*coslo(unode) ; a(3,1)=-sinlo(unode)
            a(1,2)=cosla(unode)*sinlo(unode) ; a(2,2)=-sinla(unode)*sinlo(unode) ; a(3,2)= coslo(unode)
            a(1,3)=sinla(unode)              ; a(2,3)= cosla(unode)              ; a(3,3)= 0.0_dp
            reg%time_gradient(1:3,unode) = matmul(a,b)

            call refract_locally(reg%r(unode),reg%lat(unode),reg%long(unode),velgrid, &
                 reg%time_gradient(1:3,unode))

            if (diat.and.unode == testnode) write (24,'(a20,4f14.6)') 'arrival adjusted',atime,reg%time_gradient(1:3,unode)
         endif
         vscount=vscount+1
         if(diat.and.unode == testnode) write(24,'(a15,4i8,2f12.5)') '3 node update',newalive,node1,node2,node3,&
              reg%arrivaltime(unode),reg%arrivaltime(node1)

         if (diat.and.unode == testnode) then
            write(24,*) 'updated from wn2'
            write(24,*) node1,node2,node3,reg%arrivaltime(unode),atime
            write(24,*) reg%node(unode)
            write(24,*) reg%node(node1)
            write(24,*) reg%node(node2)
            write(24,*) reg%node(node3)
           write(24,*)
         endif
      endif

   end if

   if (vscount == 2 .and. u1 /= u2 ) then
 
      if(dot_product(wn1,wn2) < 0.99999_dp) then
      print '(a10,5f12.5,3i5)','unode :',reg%r(unode),reg%lat(unode),reg%long(unode),&
           velocity(unode),reg%arrivaltime(unode),reg%node(unode)
      print '(a10,5f12.5,3i5)','node1 :',reg%r(node1),reg%lat(node1),reg%long(node1),&
           velocity(node1),reg%arrivaltime(node1),reg%node(node1)
      print '(a10,5f12.5,3i5)','node2 :',reg%r(node2),reg%lat(node2),reg%long(node2),&
           velocity(node2),reg%arrivaltime(node2),reg%node(node2)
      print '(a10,5f12.5,3i5)','node3 :',reg%r(node3),reg%lat(node3),reg%long(node3),&
           velocity(node3),reg%arrivaltime(node3),reg%node(node3)

      print *,u1,u2,det

      print *,' r1',r1
      print *,' r2',r2
      print *,' r3',r3

      print *,' n1',n1
      print *,' n2',n2
      print *,' n3',n3

      print *,' wn1',wn1
      print *,' wn2',wn2
      print *,'dp1 ',dot_product(wn1,n1) , dot_product(wn1,n2), dot_product(wn1,n3)
      print *,'dp2 ',dot_product(wn2,n1) , dot_product(wn2,n2), dot_product(wn2,n3)

      stop 'time_from_3_nodes: two different valid solutions found'
      endif
   endif

   return

end subroutine time_from_3_nodes

!**********************************************************************************************************
subroutine regular_update2(unode,anode,n_anode)
  use mod_prop
  implicit none

  integer           :: unode,n_anode,anode(n_anode)
  integer           :: iu,ju,ku,ia,ja,ka,ias,jas,kas,di,dj,dk,secnode,ioff,joff,koff
  real(kind=dp)     :: a,b,c,t1,t2,d2,det,t1r,t2r,t1lat,t2lat,t1long,t2long,velr
  integer           :: n,i,j,k
  logical           :: r_info,lat_info,long_info
  real(kind=dp)     :: deg_to_rad,tmax

  deg_to_rad=acos(-1.0_dp)/180.0_dp

  if (n_anode > 3) stop 'too many nodes for regular update'

  if (diat.and.unode == testnode) then
     write(24,*) 'regular update with',n_anode,'connected nodes'
     write(24,'(a8,3i6,f12.6)') 'node :',reg%node(unode)%i1,reg%node(unode)%i2,reg%node(unode)%i3,reg%arrivaltime(unode)
     do n=1,n_anode
     write(24,'(a8,3i6,f12.6)') 'anode :',reg%node(anode(n))%i1,reg%node(anode(n))%i2,reg%node(anode(n))%i3,&
          reg%arrivaltime(anode(n))
     end do
  endif



! integer coordinates of the node to be updated

  iu=reg%node(unode)%i1 ; ju=reg%node(unode)%i2 ; ku=reg%node(unode)%i3 


! loop over the available nodes to accumulate the coefficients of the quadratic equation for the new time 

  a=0.0_dp ; b=0.0_dp ; c=0.0_dp

  t2r=0.0_dp ; t2lat=0.0_dp ; t2long=0.0_dp

  r_info=.false. ;lat_info=.false. ;long_info=.false.  

  tmax=0.0

  do n=1,n_anode

     tmax=max(tmax,reg%arrivaltime(anode(n)))

    ! integer coordinates of the alive node used in the update
 
     ia=reg%node(anode(n))%i1 ; ja=reg%node(anode(n))%i2 ; ka=reg%node(anode(n))%i3 


     if (ia == 0) stop 'interface node in regular update'

!!!
     if (ia /= iu) then       ! the alive node is offset in r

        r_info=.true.
        ioff=ia
        di = ia - iu
        ias= iu + 2*di

        velr = 0.5_dp*(velocity(unode)+velocity(anode(n)))

        ! test if the next node in the r direction is available for a second order update

        if (ias >= 1 .and. ias <= grid%nr) then         ! if it lies inside the grid 

           if (grid%node_region(ias,ja,ka) == reg%id) then        ! if it lies inside the region

              secnode=grid%rnode_id(ias,ja,ka)             ! get the # of the next node in the regional node list

                 ! if the next node is alive and causal

                 if (reg%node_status(secnode) == 0 .and. reg%arrivaltime(anode(n)) >= reg%arrivaltime(secnode)) then

                 ! second order unsafe at high wavefront curvature

                 if (dot_product(reg%time_gradient(1:3,anode(n)),reg%time_gradient(1:3,secnode)) >= & 
                      0.71_dp/(velocity(anode(n))*velocity(secnode)) ) then

                   ! .and. &
                   ! abs(velocity(anode(n))-velocity(secnode)) < 10.2_dp .and. &
                   ! abs(velocity(unode)-velocity(anode(n))) < 10.2_dp ) then

                 ! use second order stencil

                 t1r = reg%arrivaltime(anode(n))
                 t2r = reg%arrivaltime(secnode)
!                 d2 = 1.0_dp/grid%dr0**2
!                 d2 = 0.25_dp*((velocity(unode)+velocity(anode(n)))/grid%dr0)**2
                 d2 = (velr/grid%dr0)**2
                 a = a + 2.25_dp*d2
                 b = b + (-6.0_dp*t1r + 1.5_dp*t2r)*d2
                 c = c + (4.0_dp*t1r*t1r - 2.0_dp*t1r*t2r + 0.25_dp*t2r*t2r)*d2 

                 if (diat.and.unode == testnode) write(24,*) 'anode',n,'offset in r, second order'
                 if (diat.and.unode == testnode) write(24,*) a,b,c,t1,t2

                 cycle   ! done for this neighbouring alive point

                 endif

              endif   ! alive test

           endif    ! region test

        endif    ! grid test

       ! next node is not available, use first order stencil  

        t1r = reg%arrivaltime(anode(n))
!        d2 = 1.0_dp/grid%dr0**2
!        d2 = 0.25_dp*((velocity(unode)+velocity(anode(n)))/grid%dr0)**2
        d2 = (velr/grid%dr0)**2
        a = a + d2
        b = b + (-2.0_dp*t1r)*d2
        c = c + (t1r*t1r)*d2 

        if (diat.and.unode == testnode) write(24,*) 'anode',n,'offset in r, first order'
        if (diat.and.unode == testnode) write(24,*) a,b,c,t1
      
     else

        if (ja /= ju) then    ! the alive node is offset in lat

           lat_info=.true.
           joff=ja
           dj = ja - ju
           jas= ju + 2*dj


           ! test if the next node in the lat direction is available for a second order update

           if (jas >= 1 .and. jas <= grid%nlat) then             ! if it lies inside the grid 

              if (grid%node_region(ia,jas,ka) == reg%id) then            ! if it lies inside the region

                 secnode=grid%rnode_id(ia,jas,ka)             ! get the # of the next node in the regional node list

                 ! if the next node is alive and causal

                 if (reg%node_status(secnode) == 0 .and. reg%arrivaltime(anode(n)) >= reg%arrivaltime(secnode)) then

                 ! second order unsafe at high wavefront curvature

                    if (dot_product(reg%time_gradient(1:3,anode(n)),reg%time_gradient(1:3,secnode)) >= &
                      0.71_dp/(velocity(anode(n))*velocity(secnode))) then
!                         0.71_dp/velocity(anode(n))**2) then

                    ! use second order stencil

                    t1lat = reg%arrivaltime(anode(n))
                    t2lat = reg%arrivaltime(secnode)
!                    d2 = 1.0_dp/(grid%r(iu)*grid%dlat0)**2
                    d2 = 0.25_dp*((velocity(unode)+velocity(anode(n)))/(grid%r(iu)*grid%dlat0))**2
                    a = a + 2.25_dp*d2
                    b = b + (-6.0_dp*t1lat + 1.5_dp*t2lat)*d2
                    c = c + (4.0_dp*t1lat*t1lat - 2.0_dp*t1lat*t2lat + 0.25_dp*t2lat*t2lat)*d2 

                    if (diat.and.unode == testnode) write(24,*) 'anode',n,'offset in lat, second order'
                    if (diat.and.unode == testnode) write(24,*) a,b,c,t1,t2

                    cycle   ! done for this neighbouring alive point

                    endif

                 endif   ! alive test

              endif    ! region test

           endif    ! grid test

           ! next node is not available, use first order stencil  

           t1lat = reg%arrivaltime(anode(n))
!           d2 = 1.0_dp/(grid%r(iu)*grid%dlat0)**2
           d2 = 0.25_dp*((velocity(unode)+velocity(anode(n)))/(grid%r(iu)*grid%dlat0))**2
           a = a + d2
           b = b + (-2.0_dp*t1lat)*d2
           c = c + (t1lat*t1lat)*d2    
  
           if (diat.and.unode == testnode) write(24,*) 'anode',n,'offset in lat, first order'
           if (diat.and.unode == testnode) write(24,*) a,b,c,t1

        else                  ! the alive node is offset in long

           long_info=.true.
           koff=ka
           dk = ka - ku
           kas= ku + 2*dk


           ! test if the next node in the long direction is available for a second order update

           if (kas >= 1 .and. kas <= grid%nlong) then            ! if it lies inside the grid 

              if (grid%node_region(ia,ja,kas) == reg%id) then            ! if it lies inside the region

                 secnode=grid%rnode_id(ia,ja,kas)             ! get the # of the next node in the regional node list

                 ! if the next node is alive and causal

                 if (reg%node_status(secnode) == 0 .and. reg%arrivaltime(anode(n)) >= reg%arrivaltime(secnode)) then


                 ! second order unsafe at high wavefront curvature

                    if (dot_product(reg%time_gradient(1:3,anode(n)),reg%time_gradient(1:3,secnode)) >= &
                      0.71_dp/(velocity(anode(n))*velocity(secnode))) then
!                         0.71_dp/velocity(anode(n))**2) then

                    ! use second order stencil

                    t1long = reg%arrivaltime(anode(n))
                    t2long = reg%arrivaltime(secnode)
!                    d2 = 1.0_dp/(grid%r(iu)*grid%coslat(ju)*grid%dlong0)**2
                    d2 = 0.25_dp*((velocity(unode)+velocity(anode(n)))/(grid%coslat(ju)*grid%r(iu)*grid%dlong0))**2
                    a = a + 2.25_dp*d2
                    b = b + (-6.0_dp*t1long + 1.5_dp*t2long)*d2
                    c = c + (4.0_dp*t1long*t1long - 2.0_dp*t1long*t2long + 0.25_dp*t2long*t2long)*d2 

                    if (diat.and.unode == testnode) write(24,*) 'anode',n,'offset in long, second order'
                    if (diat.and.unode == testnode) write(24,*) a,b,c,t1,t2

                    cycle   ! done for this neighbouring alive point

                    endif

                 endif   ! alive test

              endif    ! region test

           endif    ! grid test

           ! next node is not available, use first order stencil  

           t1long = reg%arrivaltime(anode(n))
!           d2 = 1.0_dp/(grid%r(iu)*grid%coslat(ju)*grid%dlong0)**2
           d2 = 0.25_dp*((velocity(unode)+velocity(anode(n)))/(grid%coslat(ju)*grid%r(iu)*grid%dlong0))**2
           a = a + d2
           b = b + (-2.0_dp*t1long)*d2
           c = c + (t1long*t1long)*d2      

           if (diat.and.unode == testnode) write(24,*) 'anode',n,'offset in long, first order'
           if (diat.and.unode == testnode) write(24,*) a,b,c,t1

        endif

     endif

  end do  ! end loop over alive points in octant

! move the RHS of the Eikonal equation (s**2) to the LHS

  c = c - 1.0_dp      !/velocity(unode)**2 !- 1.0_dp

  if (diat.and.unode == testnode) write(24,*) 'final abc'
  if (diat.and.unode == testnode) write(24,*) a,b,c

  det= b*b - 4.0_dp*a*c

  if (det < 0.0_dp) then

!     if ( 0.5_dp*sqrt(abs(det))/abs(a) < 0.05_dp) then 
     if ( 0.5_dp*sqrt(abs(det))/abs(a) < 1.0e20_dp) then 

        det=0.0_dp
        n_det_negative = n_det_negative + 1

     else

        print *,grid%dr0,grid%r(iu)*grid%dlat0,grid%r(iu)*grid%coslat(ju)*grid%dlong0
        i=0 ; j=0 ; k=0
        if (r_info) i=i+1
        if (t2r /= 0.0_dp) i=i+1
        if (lat_info) j=j+1
        if (t2lat /= 0.0_dp) j=j+1
        if (long_info) k=k+1
        if (t2long /= 0.0_dp) k=k+1
        print '(i8,6i5)',unode,reg%node(unode),i,j,k
        do n=1,n_anode
           print '(i8,3i5,4f9.3)',anode(n),reg%node(anode(n)),reg%arrivaltime(anode(n)),reg%time_gradient(1:3,anode(n))
           print *,pgrid%r(reg%node(anode(n))%i1), &
                pgrid%lat(reg%node(anode(n))%i2)/deg_to_rad, &
                pgrid%long(reg%node(anode(n))%i3)/deg_to_rad
           print *
        enddo
        print *,a,b,c
        print *,'det time=',0.5_dp*sqrt(abs(det))/abs(a)
        print *,'time with det=0',0.5_dp*(-b)/a
        stop ' negative determinant too large'

     endif

  endif

  t1 = 0.5_dp*(-b + sqrt(det))/a

! hard enforcement of causality
 
  t1=max(t1,tmax)

  if (t1 < reg%arrivaltime(unode)) then  ! node will be updated

     reg%arrivaltime(unode) = t1  ! the arrival time

     reg%time_gradient(1:3,unode)=0.0_dp   ! and the time gradient that led to this update

     if (r_info) then
        if (t2r == 0.0_dp) then
           reg%time_gradient(1,unode)=(t1-t1r)/(grid%r(iu)-grid%r(ioff))
        else
           reg%time_gradient(1,unode)=0.5_dp*(3.0_dp*t1-4.0_dp*t1r+t2r)/(grid%r(iu)-grid%r(ioff))
        endif
     endif

     if (lat_info) then
        if (t2lat == 0.0_dp) then
           reg%time_gradient(2,unode)=(t1-t1lat)/(reg%r(unode)*(grid%lat(ju)-grid%lat(joff)))
        else
           reg%time_gradient(2,unode)=0.5_dp*(3.0_dp*t1-4.0_dp*t1lat+t2lat)/(reg%r(unode)*(grid%lat(ju)-grid%lat(joff)))
        endif
     endif

     if (long_info) then
        if (t2long == 0.0_dp) then
           reg%time_gradient(3,unode)=(t1-t1long)/(reg%r(unode)*cosla(unode)*(grid%long(ku)-grid%long(koff)))
        else
           reg%time_gradient(3,unode)=0.5_dp*(3.0_dp*t1-4.0_dp*t1long+t2long)/(reg%r(unode)*cosla(unode)*&
                (grid%long(ku)-grid%long(koff)))
        endif
     endif

     ! flag the time gradient as suspect if it does not use information in all three dimensions
     suspect_time_gradient(unode)= .not.(r_info .and. lat_info .and. long_info)

  endif

  if (diat.and.unode == testnode) write(24,*) 'newt, current estimate',reg%arrivaltime(unode),t1

  return

end subroutine regular_update2

!**********************************************************************************************************
! This subroutine creates a narrow band around the set of alive points transferred from the
! refined source grid to the main source region

 
subroutine create_narrow_band(regin,vtype)
  use mod_prop
  implicit none

  type(Tregion),target                       :: regin
  integer                                    :: vtype
  integer                                    :: n,m,nt ! local variables


!  print *,'entering create narrow band'

  reg =>regin
  grid => reg%grid
  velocity => reg%velocity(:,vtype)
  velgrid => vgrid(reg%ivgrid,vtype)


  allocate(node_is_counted(reg%nnode))
! precalculate required sines and cosines

  allocate(cosla(reg%nnode),sinla(reg%nnode),coslo(reg%nnode),sinlo(reg%nnode),suspect_time_gradient(reg%nnode))

  cosla=cos(reg%lat)
  sinla=sin(reg%lat)
  coslo=cos(reg%long)
  sinlo=sin(reg%long)

  do n=1,reg%nnode

     if (reg%node_status(n) == 0) then

        call find_connected_nodes(n)

! save the list of neighbours, because the array connode will be used again when finding neighbours of the
! neighbours that have to be updated 

        n_trial=n_connode
        trialnode(1:n_trial)=connode(1:n_trial)


!  loop over neighbours and take appropopriate action


        do m=1,n_trial

           nt=trialnode(m)
           if (reg%node_status(trialnode(m)) == -1) then

              reg%node_status(trialnode(m)) = 1
              
              reg%arrivaltime(trialnode(m))= huge_time
              
              call update_nb_time(trialnode(m))

           endif

        end do

     endif

  end do

  deallocate(cosla,sinla,coslo,sinlo,node_is_counted,suspect_time_gradient)

!  print *,'leaving create_narrowband'

  return

end subroutine create_narrow_band
!-------------------------------------------------------------------------
!
! This subroutine is used only by the subroutine create_narrow_band 
! It updates the time at a point, and is only different from the normal 
! routine update_time in not using the special properties of the most recent 
! alive point (newalive)

  subroutine update_nb_time(update_node)
    use mod_prop
    implicit none

    integer       :: update_node
    integer       :: n            ! local variables
    integer       :: il,jl,kl
    integer       :: anode(20),nanode
!
! first find the cells to which the node to be updated is connected, and the list of
! ALIVE nodes in these cells

    call find_connected_cells(update_node)

! now loop over the connected cells, and if possible use the alive nodes
! in each cell to find a new arrival time at the node to be updated
 
    do n=1,n_concell


       if (concell_nnode(n) > 0) then   ! if there are any alive nodes in this cell


! copy alive nodes in this connected cell to local array
       nanode=concell_nnode(n)
       if (nanode > 20) stop 'subroutine update_time: nanode >20'
       anode(1:nanode)=connode(concell_nbase(n):concell_nbase(n)+nanode-1)


             if (concell_is_regular(n)) then

                ! do a regular fast marching update, second order if possible, in this octant

                call regular_update2(update_node,anode(1:nanode),nanode)

             else

             ! do an irregular update in this octant

                select case (nanode)
          
                case(1)     ! there is 1 alive node in the cell

                   call time_from_1_node(update_node,anode(1))

                case(2)     ! there are 2 alive nodes in the cell

                   do jl=1,2
                      call time_from_1_node(update_node,anode(jl))
                   end do

                   call time_from_2_nodes(update_node,anode(1),anode(2))

                case(3)     ! there are 3 alive nodes in the cell

                   do jl=1,3
                      call time_from_1_node(update_node,anode(jl))
                   end do

                   do jl=1,2
                      do kl=jl+1,3
                         call time_from_2_nodes(update_node,anode(jl),anode(kl))
                      end do
                   end do

                   call time_from_3_nodes(update_node,anode(1),anode(2),anode(3))

                case(4:)     ! there are more than 3 alive nodes in the cell

                   do jl=1,nanode
                      call time_from_1_node(update_node,anode(jl))
                   end do 

                   do jl=2,nanode-1
                      do kl=jl+1,nanode
                      call time_from_2_nodes(update_node,anode(jl),anode(kl))
                      end do
                   end do

                   do il=2,nanode-2
                      do jl=il+1,nanode-1
                         do kl=jl+1,nanode
                         
                            call time_from_3_nodes(update_node,anode(il),anode(jl),anode(kl))
                         end do
                      end do
                   end do


                end select

             endif   ! regular or irregular

       endif   ! if there are any alive nodes in this cell


    end do ! loop over connected cells

    return

  end subroutine update_nb_time


!*********************************************************************************************************
!-------------------------------------------------------------------------
! below are the subroutines used to manipulate the binary tree
! NOTE!!!!!! the naming of the directions in the binary tree used for sorting 
! the narrow band has been reversed from the definition used by Sethian. This means:
! the root (smallest value) is at the BOTTOM, not the top
! UP is in the direction of increasing values, and increasing index in the tree array
! DOWN is in the direction of decreasing values, and decreasing index in the tree array
!
!-----------------------------------------------------------------
! adds a node to the binary tree at its correct position

  subroutine add_node_to_tree(nn)
    use mod_prop
    implicit none

    integer :: nn, parent_index, child_index, temp
!
! First, increase the size of the tree by one.
!
    ntree=ntree+1
!
! Put new value at top of tree
!
    reg%node_status(nn)=ntree
    node_from_tree_ind(ntree)=nn
!
! Now filter the new value down to its correct position
!
    child_index = ntree
    parent_index = ntree/2

    do while (parent_index > 0)

       if (reg%arrivaltime(nn) < reg%arrivaltime(node_from_tree_ind(parent_index))) then

          reg%node_status(nn)=parent_index
          reg%node_status(node_from_tree_ind(parent_index)) = child_index
          temp=node_from_tree_ind(child_index)
          node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
          node_from_tree_ind(parent_index) = temp
          child_index=parent_index
          parent_index=child_index/2

       else

          parent_index = 0

       endif

    end do
    return

  end subroutine add_node_to_tree

!-----------------------------------------------------------------
! removes the root from the tree and adjusts the tree structure above

  subroutine remove_root_from_tree
    use mod_prop
    implicit none

    integer  :: parent_index,child_index,temp


! if only one node left in the tree, nothing is left to be done

    if(ntree == 1) then ; ntree=0 ; return ; endif


! put the node with largest travel time at the root position

    reg%node_status(node_from_tree_ind(ntree))=1
    node_from_tree_ind(1)=node_from_tree_ind(ntree)


! Reduce size of tree by one

    ntree=ntree-1


! Now move the new root up to its correct position

    parent_index=1
    child_index=2*parent_index

    do while (child_index < ntree)

! Check which of the two children is smallest - use the smallest

       if(reg%arrivaltime(node_from_tree_ind(child_index)) > reg%arrivaltime(node_from_tree_ind(child_index+1))) &
            child_index=child_index+1


!  Check whether the child is smaller than the parent; if so, then swap,
!  if not, then we are done

       if(reg%arrivaltime(node_from_tree_ind(child_index)) < reg%arrivaltime(node_from_tree_ind(parent_index))) then

          reg%node_status(node_from_tree_ind(parent_index)) = child_index
          reg%node_status(node_from_tree_ind(child_index)) = parent_index
          temp = node_from_tree_ind(child_index)
          node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
          node_from_tree_ind(parent_index) = temp
          parent_index=child_index
          child_index=2*parent_index

       else

          child_index=ntree+1

       endif

    end do

    if (child_index == ntree) then

      if(reg%arrivaltime(node_from_tree_ind(child_index)) < reg%arrivaltime(node_from_tree_ind(parent_index))) then

          reg%node_status(node_from_tree_ind(parent_index)) = child_index
          reg%node_status(node_from_tree_ind(child_index)) = parent_index
          temp = node_from_tree_ind(child_index)
          node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
          node_from_tree_ind(parent_index) = temp

       endif
    endif

    return

  end subroutine remove_root_from_tree


!-----------------------------------------------------------------
! updates the tree structure after the value of the traveltime at a node has been
! updated (which should always be a decrease)

  subroutine update_tree(nn)
    use mod_prop
    implicit none

    integer  :: nn,temp,parent_index,child_index

!
! Filter the updated value down to its correct position
!
    child_index=reg%node_status(nn)
    parent_index=child_index/2

    do while (parent_index > 0)
       if (reg%arrivaltime(nn) < reg%arrivaltime(node_from_tree_ind(parent_index))) then

          reg%node_status(nn)=parent_index
          reg%node_status(node_from_tree_ind(parent_index)) = child_index
          temp=node_from_tree_ind(child_index)
          node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
          node_from_tree_ind(parent_index) = temp
          child_index=parent_index
          parent_index=child_index/2

       else

          parent_index = 0

       endif

    end do
    return

  end subroutine update_tree

!----------------------------------------------------------------------------------------
!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP(A,INDX,status)
   use mod_3dfm
 implicit none
 real(kind=dp), PARAMETER              :: TINY=1.5D-16
 REAL(kind=dp)                         :: AMAX,DUM, SUM, A(3,3),VV(3)
 INTEGER                               :: status, INDX(3),i,j,k,imax

 status=0

 DO I=1,3
   AMAX=0._dp
   DO J=1,3
     IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     status = 1
     RETURN
   END IF
   VV(I) = 1.0_dp / AMAX
 END DO ! i loop

 DO J=1,3
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.0_dp
   DO I=J,3
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*ABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,3
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.3) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,3
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

 RETURN
 END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine LUBKSB(A,INDX,B)
   use mod_3dfm
 real(kind=dp)         ::SUM, A(3,3),B(3)
 INTEGER              ::INDX(3),ii,ll,i,j

 II = 0

 DO I=1,3
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.0_dp) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=3,1,-1
   SUM = B(I)
   IF(I < 3) THEN
     DO J=I+1,3
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine LUBKSB

!*************************************************************************************************************

! this subroutine is used to find the gradient of the arrival time at a node
! if the gradient stored during the fast marching is suspect, which is always 
! when it was derived from a regular update from less than 3 nodes
! If all cells connected to the node are regular, just use finite differences on the grid.
! Otherwise, the local gradient is derived from a least squares fit of a constant gradient
! field to the values of the arrival times at the nodes connected to the one at which
! the gradient is evaluated

subroutine fit_gradient_at_node(reg,centernode,dtdr,dtdlat,dtdlong)
use mod_3dfm
implicit none

type(Tregion)                :: reg

integer,parameter            :: maxfitnodes = 50
integer                      :: n,m,i,j,k,i1,i2,i3,n_cnode,nminus,nplus
integer                      :: ii,jj,kk,nfit,icell
integer                      :: centernode,n_concell,mstore(maxfitnodes)
type(Tinteger_coordinates)   :: concell(8)
logical                      :: concell_irregular(8)
real(kind=dp)                :: r(maxfitnodes),lat(maxfitnodes),long(maxfitnodes),atime(maxfitnodes)
real(kind=dp)                :: dtdr,dtdlat,dtdlong,wmax,wmin,dx,dy,dz,dw

real(kind=dp),allocatable ::v(:,:),w(:),sol(:)
real(kind=dp),allocatable ::a(:,:),b(:)

type(Tintersection),pointer      :: isec
type(Tpropagation_grid),pointer :: grid


grid => reg%grid


! first make list of connecetd cells


! store the identifiers of the centernode node in local variables. Reminder:
! for regular grid nodes i1,i2,i3  correspond to ir,ilat,ilong of the node
! for intersection nodes i1,i2,i3 are 0, intersection #, node # in intersection 

   i1 = reg%node(centernode)%i1 ; i2 = reg%node(centernode)%i2 ; i3 = reg%node(centernode)%i3


! make a list of grid cells of which the node is part

   if (i1 /= 0) then   ! centernode is a regular grid node. Use spatial correlation of grid

      n_concell = 0
      do i=0,1
         if ((i1-i>0) .and. (i1-i<grid%nr)) then
            do j=0,1
               if ((i2-j>0) .and. (i2-j<grid%nlat)) then
                  do k=0,1
                     if ((i3-k>0) .and. (i3-k<grid%nlong)) then
                        n_concell=n_concell+1
                        concell(n_concell)%ir=i1-i  
                        concell(n_concell)%ilat=i2-j 
                        concell(n_concell)%ilong=i3-k
                     endif
                  end do
               end if
            end do
         end if
      end do

   else          ! centernode is an intersection node. Use the connections that have been found before

      if (i2 == reg%itop%id) then ; isec => reg%itop ; else ; isec => reg%ibot ; endif
      n_concell = 0
      do n=1,8
         icell = isec%ccell_from_inode(n,i3)
         if (icell > 0) then
            n_concell=n_concell+1
            concell(n_concell)%ir= isec%ccells(icell)%ir
            concell(n_concell)%ilat=isec%ccells(icell)%ilat 
            concell(n_concell)%ilong=isec%ccells(icell)%ilong
         endif
      enddo

   endif

   ! test whether any connected cell is cut by an interface, in that case the pointer to the local list of
   ! interfaces cutting the cell has been allocated

   do n=1,n_concell
      concell_irregular(n) = associated(grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p)
   end do


   ! if all connected cells are regular, gradient is easy

   if (count(concell_irregular(1:n_concell)) == 0) then

      if (i1 == 0) stop 'fit_gradient_at_node : interface node in regular gradient evaluation'

      k = max(1,i1-1)
      if ( k<1 .or. k>grid%nr .or. i2<1 .or. i2>grid%nlat .or. i3<1 .or. i3>grid%nlong ) then

         print *,k,i1,i2,i3
         stop

      endif

      nminus = grid%rnode_id(max(1,i1-1),i2,i3)
      nplus  = grid%rnode_id(min(grid%nr,i1+1),i2,i3)
      dtdr = (reg%arrivaltime(nplus)-reg%arrivaltime(nminus))/(reg%r(nplus)-reg%r(nminus))

      nminus = grid%rnode_id(i1,max(1,i2-1),i3)
      nplus  = grid%rnode_id(i1,min(grid%nlat,i2+1),i3)
      dtdlat = (reg%arrivaltime(nplus)-reg%arrivaltime(nminus))/(reg%r(centernode)*(reg%lat(nplus)-reg%lat(nminus)))

      nminus = grid%rnode_id(i1,i2,max(1,i3-1))
      nplus  = grid%rnode_id(i1,i2,min(grid%nlong,i3+1))
      dtdlong = (reg%arrivaltime(nplus)-reg%arrivaltime(nminus))/&
           (reg%r(centernode)*reg%coslat(centernode)*(reg%long(nplus)-reg%long(nminus)))

      return

   endif



! some connected cells are irregular, use fitting mode for gradient

   n_cnode=0

   do n=1,n_concell

! find the intersection nodes of connected cell n

!   explanation: each intersection has a 1-D list of cells cut by the interface, and a list of intersection
!   nodes that are part of each cut cell. Each regular grid cell has a pointer ccind_from_3dc(i,j,k)%p
!   (Cut Cell INDex FROM 3D Coordinates)
!   associated with it, where p is a pointer to an integer array with as many elements as there are intersections.
!   If a cell is cut by interface n, the pointer ccind_from_3dc(i,j,k)%p is allocated, and the variable
!   ccind_from_3dc(i,j,k)%p(n) contains the array index of cell (i,j,k) in the 1D cut cell list of intersection n 

      if (concell_irregular(n)) then

! if so, check if the cell is cut by the top intersection

     ! icell is the index of the current connected cell in the list of cells cut by interface reg%itop.
                  
         icell = grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%itop%iface_id)

     ! if icell == 0 the cell is not cut be the top interface
    
         if(icell /= 0) then
 
            isec => reg%itop
            do jj=1,isec%n_inodes(icell)

!   m is the node number in the regional node list of node  jj in the list of 
!   interface nodes that are part of cut cell icell

               m=isec%rbel_node_id(isec%inodes(jj,icell))
               if (m == centernode) cycle

               if ( n_cnode == 0 ) then  
                  n_cnode=1
                  mstore(n_cnode) = m
                  r(n_cnode)     = reg%r(m)
                  lat(n_cnode)   = reg%lat(m)
                  long(n_cnode)  = reg%long(m)
                  atime(n_cnode) = reg%arrivaltime(m)
               else
                  if ( count(mstore(1:n_cnode) == m) == 0 ) then
                     n_cnode=n_cnode+1
                     mstore(n_cnode) = m
                     r(n_cnode)     = reg%r(m)
                     lat(n_cnode)   = reg%lat(m)
                     long(n_cnode)  = reg%long(m)
                     atime(n_cnode) = reg%arrivaltime(m)
                  endif
               endif

            end do

         end if


! then check the bottom intersection

         icell = grid%ccind_from_3dc(concell(n)%ir,concell(n)%ilat,concell(n)%ilong)%p(reg%ibot%iface_id)
         if (icell /= 0) then

            isec => reg%ibot
            do jj=1,isec%n_inodes(icell)

               m=isec%rabo_node_id(isec%inodes(jj,icell))
               if (m == centernode) cycle

               if ( n_cnode == 0 ) then  
                  n_cnode=1
                  mstore(n_cnode) = m
                  r(n_cnode)     = reg%r(m)
                  lat(n_cnode)   = reg%lat(m)
                  long(n_cnode)  = reg%long(m)
                  atime(n_cnode) = reg%arrivaltime(m)
               else
                  if ( count(mstore(1:n_cnode) == m) == 0 ) then
                     n_cnode=n_cnode+1
                     mstore(n_cnode) = m
                     r(n_cnode)     = reg%r(m)
                     lat(n_cnode)   = reg%lat(m)
                     long(n_cnode)  = reg%long(m)
                     atime(n_cnode) = reg%arrivaltime(m)
                  endif
               endif

            end do

         endif

      endif


! then find the regular grid nodes of connected cell n

      do i=0,1
         ii=concell(n)%ir+i
         do j=0,1
            jj=concell(n)%ilat+j
            do k=0,1
               kk=concell(n)%ilong+k

                 if (grid%node_region(ii,jj,kk) == reg%id) then  ! node has to belong to the current region

                     m=grid%rnode_id(ii,jj,kk)
                     if (m == centernode) cycle

                     if ( n_cnode == 0 ) then  
                        n_cnode=1
                        mstore(n_cnode) = m
                        r(n_cnode)     = reg%r(m)
                        lat(n_cnode)   = reg%lat(m)
                        long(n_cnode)  = reg%long(m)
                        atime(n_cnode) = reg%arrivaltime(m)
                     else
                        if ( count(mstore(1:n_cnode) == m) == 0 ) then
                           n_cnode=n_cnode+1
                           mstore(n_cnode) = m
                           r(n_cnode)     = reg%r(m)
                           lat(n_cnode)   = reg%lat(m)
                           long(n_cnode)  = reg%long(m)
                           atime(n_cnode) = reg%arrivaltime(m)
                        endif
                     endif


                  endif


            end do
         end do
      end do

   end do  ! loop over connected cells


! now do the actual least squares fit using singular value decomposition

      nfit=3

      allocate(a(n_cnode,nfit),b(n_cnode),v(nfit,nfit),w(nfit),sol(nfit))

      do i=1,n_cnode

         dx=r(i)-reg%r(centernode)
         dy=reg%r(centernode)*(lat(i)-reg%lat(centernode))
         dz=reg%r(centernode)*reg%coslat(centernode)*(long(i)-reg%long(centernode))
         dw=1.0_dp/(sqrt(dx*dx+dy*dy+dz*dz) + grid%tolerance)

         a(i,1)=dx*dw
         a(i,2)=dy*dw
         a(i,3)=dz*dw
         b(i)  =(atime(i)- reg%arrivaltime(centernode))*dw

      end do

      call svdcmp(a,n_cnode,nfit,n_cnode,nfit,w,v)

!	find maximum singular value

      wmax=maxval(abs(w(1:nfit)))

!	define "small"

      wmin=wmax*1.e-6_dp
 
!	zero the "small" singular values
      where (abs(w(1:nfit)) < wmin) w(1:nfit)=0.0d0

      call svbksb(a,w,v,n_cnode,nfit,n_cnode,nfit,b,sol)
   
      dtdr=sol(1)
      dtdlat=sol(2)
      dtdlong=sol(3)

      deallocate(a,b,v,w,sol)

      return


 end subroutine fit_gradient_at_node

!**********************************************************************************************************






rays.f90
subroutine trace_ray_from_receiver(rec,s,ray)
  use mod_3dfm_nointerfaces
  implicit none

  type(Tsource),target         :: s                    ! the source of the ray
  type(Treceiver)              :: rec                  ! the receiver at which it arrives
  type(Tray),target            :: ray                  ! the ray to be traced

  real(kind=dp)                :: r_in,lat_in,long_in       ! starting points for each ray section
  real(kind=dp)                :: r_out,lat_out,long_out    ! end points of each ray section
  integer                      :: n,m,nsections,ntf,n_shift               ! counter, number of sections on the rays

  logical                      :: ldum=.false.
  logical                      :: outside
  real(kind=dp)                :: interpolate_interface


  nsections=s%path(ray%raypath_id)%n_tf
  ray%nsections=nsections
  allocate(ray%section(nsections))
  do n=1,nsections ; call ray_section_defaults(ray%section(n)) ; end do

  r_in = rec%r ; lat_in=rec%lat ; long_in=rec%long

! first fill in some info on all the ray sections

  do n=nsections,1,-1

     ray%section(n)%ray => ray
     ray%section(n)%source => s
     ray%section(n)%place_in_sequence = n

     ntf = s%path(ray%raypath_id)%tf_sequence(n)
     if (ntf > 0) then
        ray%section(n)%tf => s%time_field(ntf)
     else
        print *,'tracing ray with 1 timefield'
        if (nsections /= 1 ) stop 'error 1 in trace ray from receiver'
        if (rec%r >= interpolate_interface(rec%lat,rec%long,intrface(s%topint_id))) ntf=1
        if (rec%r <  interpolate_interface(rec%lat,rec%long,intrface(s%botint_id))) ntf=2
        ray%section(n)%tf => s%time_field(ntf)
     endif

     ray%section(n)%reg => ray%section(n)%tf%reg
     if (n < nsections) ray%section(n)%istart => intersection(s%path(ray%raypath_id)%sequence(2*n+1))
     if (n > 1) ray%section(n)%iend => intersection(s%path(ray%raypath_id)%sequence(2*n-1))
     if (n == 1 .and. s%is_teleseismic) ray%section(n)%iend => intersection(n_intersections)
     if (n < nsections) ray%section(n)%istart => ray%section(n+1)%iend

  end do


! now do the actual integration section by section

  do n=nsections,1,-1

     if (n == nsections) then  ! start by storing the time and direction of arrival of this ray at the receiver

        call interpolate_arrivaltime(ray%section(n)%tf,rec%r,rec%lat,rec%long,ray%receiver_time)

        call interpolate_time_gradient(ray%section(n)%tf,rec%r,rec%lat,rec%long,&
             ray%receiver_time_gradient(1),ray%receiver_time_gradient(2),ray%receiver_time_gradient(3), &
             ldum,outside,ray%section(n))
           
     endif


     call trace_ray_section(ray%section(n),r_in,lat_in,long_in,r_out,lat_out,long_out,outside)

     if (outside) exit  ! exit if the ray path has gone through the sides of the box

     if (.not.ray%valid) exit

     ! make end points the starting points of the next section unless it was the last section
     if (n > 1) then
        r_in = r_out ; lat_in=lat_out ; long_in=long_out
     endif

  end do

  if (outside.and.s%is_local) stop 'error in ray tracing:ray outside grid for local source'

  if (outside.and.n>1) then
     
   ! if the ray from a teleseismic source has exited the box through the sides
   ! the info on the ray has to be modified, unless it was the last section....  

     print *,'ray exited grid before bottom region was reached'
     do m=n,nsections
        n_shift=m-n+1
        ray%section(n_shift)=ray%section(m)
        ray%section(n_shift)%place_in_sequence = n_shift
     end do
     ray%nsections=nsections-n+1

  endif


  return

end subroutine trace_ray_from_receiver

!***************************************************************************************************************

subroutine trace_ray_section(raysec,r_in,lat_in,long_in,r_out,lat_out,long_out,outside)
  use mod_3dfm
  implicit none

  type(Tray_section)                          :: raysec          ! the ray section
  type(Ttime_field),pointer                   :: tf              ! the time field to be used
  type(Ttime_field),pointer                   :: tf_next         ! the next time field as viewed 
                                                                 ! from the direction of tracing
  type(Tpropagation_grid),pointer            :: grid            ! the grid on which the time field is based
  type(Tinterface),pointer                    :: iftop,ifbot     ! the top, bottom interfaces of the time field 
  type(Tinterface),pointer                    :: ifend   ! the interface at which the ray section ends
  integer                                     :: ifstart_id,ifend_id
  type(Tsource),pointer                       :: s               ! the source of the ray
  real(kind=dp)                               :: r_in,lat_in,long_in   ! ray section start
  real(kind=dp)                               :: r_out,lat_out,long_out,lat,long   ! ray section end
  real(kind=dp)                               :: interpolate_interface ! external function

  real(kind=dp),allocatable,dimension(:,:)    :: pt     ! temporary storage for ray section
  real(kind=dp)                               :: pp(3),grad0(3),grad1(3),r_interface
  real(kind=dp)                               :: grad2(3),gtp_in,norm(3),h
  real(kind=dp)                               :: dlmax,dl ! max, actual step size of ray integration
  real(kind=dp)                               :: v_average = 5.0_dp   ! used to estimate max step size
  real(kind=dp)                               :: p_tolerance          ! stopping criterion for distance to interface
  real(kind=dp)                               :: s_tolerance          ! stopping criterion for distance to source
  real(kind=dp)                               :: x1,x2,f,fmid,xmid,dx,rtbis
  integer                                     :: maxpoints = 10000
  integer                                     :: n,nn,iter,seqnum,start_direction,i
  logical                                     :: section_is_refraction,outside

  logical  :: verbose

!  print *,'entering trace ray section'

  ! define some local pointers

  tf => raysec%tf
  iftop => intrface(tf%reg%itop%iface_id)
  ifbot => intrface(tf%reg%ibot%iface_id)
  grid => pgrid
  s => raysec%source

  seqnum=raysec%place_in_sequence

! a ray section with others before and after
  if (seqnum < raysec%ray%nsections .and.seqnum > 1 ) then 
     ifstart_id = raysec%istart%id
     ifend_id = raysec%iend%id
     ifend => intrface(ifend_id)
     tf_next => raysec%ray%section(seqnum-1)%tf
     start_direction = raysec%ray%section(seqnum+1)%reg%id - tf%reg%id
     section_is_refraction = raysec%ray%section(seqnum-1)%reg%id /= raysec%ray%section(seqnum)%reg%id
  endif

! a ray section that starts at the receiver with others after
  if (seqnum == raysec%ray%nsections .and. seqnum > 1) then 
     ifstart_id = 0
     ifend_id = raysec%iend%id
     ifend => intrface(ifend_id)
     tf_next => raysec%ray%section(seqnum-1)%tf
     start_direction = 0
     section_is_refraction = raysec%ray%section(seqnum-1)%reg%id /= raysec%ray%section(seqnum)%reg%id
  endif

! the last ray section to the source
  if (raysec%place_in_sequence == 1 ) then 

     if (s%is_local) then
        ifstart_id = 0
        ifend_id = 0
        nullify(ifend)
        nullify(tf_next)
        if (seqnum+1 <= raysec%ray%nsections) then
           start_direction = raysec%ray%section(seqnum+1)%reg%id - tf%reg%id
        else
           start_direction=0
        endif
     endif

     if (s%is_teleseismic) then
        ifstart_id = 0 ! raysec%istart%id
        ifend_id = raysec%iend%id
        ifend => intrface(ifend_id)
        nullify(tf_next) 
        if (seqnum+1 <= raysec%ray%nsections) then
           start_direction = raysec%ray%section(seqnum+1)%reg%id - tf%reg%id
        else
           start_direction=0
        endif
        section_is_refraction=.false.  ! must be set so that the next time field is not tested
                                       ! for total reflection
     endif

  endif

  ! test whether time field exists
  if (.not.associated(tf%time_gradient)) stop 'trace_ray_section: no timefield gradient'


  ! allocate local array for temporary storage of this ray section
  allocate(pt(3,maxpoints))

  dlmax=0.1_dp*grid%dr0*v_average   ! the maximum step size for the integration
  p_tolerance =5.e-4_dp*grid%dr0    ! maximum allowed difference between predictor and 
                                    ! predictor+corrector step in integration 
  s_tolerance = 1.0_dp*grid%dr0     ! distance to the source at which the integration stops

!  print *
!  print *
!  print '(a25,i5,3f15.7)','starting ray section',raysec%place_in_sequence,r_in,lat_in,long_in
!  print *,'into region',tf%reg%id, 'with top/bot iface',iftop%id,ifbot%id  !tf%reg%itop%id,tf%reg%ibot%id
!  print *,'using time field',tf%id
!  print *,'start/end ifaces',ifstart_id,ifend_id
!  print *,'s_tolerance',s_tolerance


! set starting point
  n=1
  pt(1,1)=r_in ; pt(2,1)=lat_in ; pt(3,1)=long_in 



ploop:  do

      ! find the next point 

           n=n+1

           if (n > maxpoints) then
              print *,'more than maxpoints in raysection'
              print '(i6,3f12.5)',n-1,pt(1:3,n-1)
              print '(a6,3f12.5)','source',s%r,s%lat,s%long
              stop
           endif

           dl=dlmax

           if (n ==0) then ; verbose =.true. ; else; verbose=.false. ; endif
            call interpolate_time_gradient(tf,pt(1,n-1),pt(2,n-1),pt(3,n-1),grad0(1),grad0(2),grad0(3) &
                ,verbose,outside,raysec)
            if (outside) then ; n=n-1 ; exit ploop ; endif
           if (.not.raysec%ray%valid) then; deallocate(pt) ; return ; endif
           verbose=.false.

           ! first guess at the position of the next point
           ! note: time gradients have units time/length, convert to angle for angular coordinates

           pp(1)=pt(1,n-1)-dl*grad0(1)
           pp(2)=pt(2,n-1)-dl*grad0(2)/pp(1)          
           pp(3)=pt(3,n-1)-dl*grad0(3)/(pp(1)*cos(pp(2)))


           ! now test if we can get into the next requested region at all

           if (n == 2 .and. start_direction /= 0 ) then ! on the first step and when the starting direction is defined

              call interface_normal(lat_in,long_in,intrface(raysec%istart%id),norm(1),norm(2),norm(3),h) 

              ! if -grad(T) points straight back through the interface at which we start, quit this section

              if (dot_product(grad0,start_direction*norm) > 0.0_dp) then

                 print *,'can not get into region',tf%reg%id,' from region',raysec%ray%section(seqnum+1)%reg%id
!                 print '(i5,6f12.5)',raysec%place_in_sequence,grad0,start_direction*norm

                 raysec%ray%valid=.false.
                 deallocate(pt)
                 return

              endif

           endif


           ! test if an interface is crossed

           r_interface=interpolate_interface(pp(2),pp(3),iftop)
           if (pp(1) >= r_interface) then
              
              if (iftop%id == ifend_id) then

                 if (section_is_refraction) then  

                    ! when the ray is refracted into the region check for total reflection (headwave),
                    ! if so keep going along the interface until the ray can get through

                    call interface_normal(pp(2),pp(3),intrface(ifend_id),norm(1),norm(2),norm(3),h)
                    call interpolate_time_gradient(tf,r_interface,pp(2),pp(3),grad1(1),grad1(2),grad1(3) &
                         ,verbose,outside,raysec)   
                    call interpolate_time_gradient(tf_next,r_interface,pp(2),pp(3),grad2(1),grad2(2),grad2(3) &
                         ,verbose,outside,raysec)

                    if (outside) then ; n=n-1 ; exit ploop ; endif
                    if (.not.raysec%ray%valid) then; deallocate(pt) ; return ; endif

                    gtp_in =sqrt(sum((grad2 - dot_product(norm,grad2)*norm)**2))


                    if (gtp_in <= sqrt(sum(grad1**2))) then  ! the ray can pass through

                       exit ploop

                    else   ! total reflection at this point, keep going along the interface

                       pt(2,n)=pp(2)          
                       pt(3,n)=pp(3)
                       pt(1,n)=r_interface
                       raysec%headwave =.true.
                       raysec%ray%headwave =.true.     

                       ! special return condition if source is reached
                       if(s%is_local .and. raysec%place_in_sequence == 1) then
                          dx= sqrt((pt(1,n)-s%r)**2+(s%r*(pt(2,n)-s%lat))**2+(s%r*s%coslat*(pt(3,n)-s%long))**2)
                          if (dx < s_tolerance) then
                             pt(1,n)=s%r ;pt(2,n)=s%lat ;pt(3,n)=s%long 
                             raysec%npoints=n
                             allocate(raysec%point(3,n))
                             raysec%point(1:3,1:n)=pt(1:3,1:n)
                             deallocate(pt)
                             return
                          endif
                       endif

                       cycle ploop

                    endif

                 else   ! section is a reflection, we can stop this section here

                    exit ploop

                 endif

              else  ! if the ray hits the non-end interface, the requested phase is a diffraction
                    ! just keep going along the interface and hopefully at some point the ray will detach again
                 
                 pt(2,n)=pp(2)          
                 pt(3,n)=pp(3)
                 pt(1,n)=r_interface
                 raysec%diffracted =.true.
                 raysec%ray%diffracted =.true.  

                 ! special return condition if source is reached
                 if(s%is_local .and. raysec%place_in_sequence == 1) then
                    dx= sqrt((pt(1,n)-s%r)**2+(s%r*(pt(2,n)-s%lat))**2+(s%r*s%coslat*(pt(3,n)-s%long))**2)
                    if (dx < s_tolerance) then
                       pt(1,n)=s%r ;pt(2,n)=s%lat ;pt(3,n)=s%long 
                       raysec%npoints=n
                       allocate(raysec%point(3,n))
                       raysec%point(1:3,1:n)=pt(1:3,1:n)
                       deallocate(pt)
                       return
                    endif
                 endif

                 cycle ploop
                 
              endif

           endif  ! if top interface crossed


           r_interface=interpolate_interface(pp(2),pp(3),ifbot)
           if (pp(1) <= r_interface) then

              if (ifbot%id == ifend_id) then

                 if (section_is_refraction) then  

                    ! when the ray is refracted into the region check for total reflection (headwave),
                    ! if so keep going along the interface until the ray can get through

                    call interface_normal(pp(2),pp(3),intrface(ifend_id),norm(1),norm(2),norm(3),h)
                    call interpolate_time_gradient(tf,r_interface,pp(2),pp(3),grad1(1),grad1(2),grad1(3) &
                         ,verbose,outside,raysec)   
                    call interpolate_time_gradient(tf_next,r_interface,pp(2),pp(3),grad2(1),grad2(2),grad2(3) &
                         ,verbose,outside,raysec) 

                    if (outside) then ; n=n-1 ; exit ploop ; endif
                    if (.not.raysec%ray%valid) then; deallocate(pt) ; return ; endif

                    gtp_in =sqrt(sum((grad2 - dot_product(norm,grad2)*norm)**2))

                    if (gtp_in <= sqrt(sum(grad1**2))) then  ! the ray can pass through

!                       print '(a30,4f12.5)','refracting through top 1  ',r_interface,pp(1:3)

                       exit ploop

                    else   ! total reflection at this point, keep going along the interface

                       pt(2,n)=pp(2)          
                       pt(3,n)=pp(3)
                       pt(1,n)=r_interface
                       raysec%headwave =.true.
                       raysec%ray%headwave =.true. 
    
                       ! special return condition if source is reached
                       if(s%is_local .and. raysec%place_in_sequence == 1) then
                          dx= sqrt((pt(1,n)-s%r)**2+(s%r*(pt(2,n)-s%lat))**2+(s%r*s%coslat*(pt(3,n)-s%long))**2)
                          if (dx < s_tolerance) then
                             pt(1,n)=s%r ;pt(2,n)=s%lat ;pt(3,n)=s%long 
                             raysec%npoints=n
                             allocate(raysec%point(3,n))
                             raysec%point(1:3,1:n)=pt(1:3,1:n)
                             deallocate(pt)
                             return
                          endif
                       endif

                       cycle ploop

                    endif

                 else   ! section is a reflection, we can stop this section here

                    exit ploop

                 endif

              else  ! if the ray hits the non-end interface, the requested phase is a diffraction
                    ! just keep going along the interface and hopefully at some point the ray will detach again

                 pt(2,n)=pp(2)          
                 pt(3,n)=pp(3)
                 pt(1,n)=r_interface
                 raysec%diffracted =.true.
                 raysec%ray%diffracted =.true.

                 ! special return condition if source is reached
                 if(s%is_local .and. raysec%place_in_sequence == 1) then
                    dx= sqrt((pt(1,n)-s%r)**2+(s%r*(pt(2,n)-s%lat))**2+(s%r*s%coslat*(pt(3,n)-s%long))**2)
                    if (dx < s_tolerance) then
                       pt(1,n)=s%r ;pt(2,n)=s%lat ;pt(3,n)=s%long 
                       raysec%npoints=n
                       allocate(raysec%point(3,n))
                       raysec%point(1:3,1:n)=pt(1:3,1:n)
                       deallocate(pt)
                       return
                    endif
                 endif

                 cycle ploop

              endif
           endif

           ! iterate the step until result is sufficiently accurate
           iter = 1
itloop   : do 

              pp(1)=pt(1,n-1)-dl*grad0(1)
              pp(2)=pt(2,n-1)-dl*grad0(2)/pp(1)          
              pp(3)=pt(3,n-1)-dl*grad0(3)/(pp(1)*cos(pp(2)))


              ! calculate gradient at provisional next position

              call interpolate_time_gradient(tf,pp(1),pp(2),pp(3),grad1(1),grad1(2),grad1(3),verbose,outside,raysec)
              if (outside) then ; n=n-1 ; exit ploop ; endif
              if (.not.raysec%ray%valid) then; deallocate(pt) ; return ; endif

              ! repeat the step with the average of the gradients at the aoriginal and provisional points

              pt(1,n)=pt(1,n-1)-dl*(grad0(1)+grad1(1))*0.5_dp
              pt(2,n)=pt(2,n-1)-dl*(grad0(2)+grad1(2))/(pt(1,n)+pt(1,n-1))          
              pt(3,n)=pt(3,n-1)-dl*(grad0(3)+grad1(3))/((pt(1,n)+pt(1,n-1))*cos(pt(2,n)))

              ! test again if an interface is crossed during the iteration, if so exit

              r_interface=interpolate_interface(pt(2,n),pt(3,n),iftop)
              if (pt(1,n) >= r_interface) then

                 if (iftop%id == ifend_id) then

                    if (section_is_refraction) then  

                       ! when the ray is refracted into the region check for total reflection (headwave),
                       ! if so keep going along the interface until the ray can get through

                       call interface_normal(pt(2,n),pt(3,n),intrface(ifend_id),norm(1),norm(2),norm(3),h)
                       call interpolate_time_gradient(tf,r_interface,pt(2,n),pt(3,n),grad1(1),grad1(2),grad1(3) &
                            ,verbose,outside,raysec)   
                       call interpolate_time_gradient(tf_next,r_interface,pt(2,n),pt(3,n),grad2(1),grad2(2),grad2(3) &
                            ,verbose,outside,raysec) 

                       if (outside) then ; n=n-1 ; exit ploop ; endif
                       if (.not.raysec%ray%valid) then; deallocate(pt) ; return ; endif

                       gtp_in =sqrt(sum((grad2 - dot_product(norm,grad2)*norm)**2))

                       if (gtp_in <= sqrt(sum(grad1**2))) then ! the ray can get through here

!                          print '(a30,4f12.5)','refracting through top 2  ',r_interface,pt(1:3,n)
                          pp=pt(1:3,n)
                          exit ploop

                       else  ! total reflection, keep going along the interface

                          pt(1,n)=r_interface
                          raysec%headwave =.true.
                          raysec%ray%headwave =.true.   
                          cycle ploop

                       endif

                    else  ! in case the section is a reflection

                       pp=pt(1:3,n)
                       exit ploop

                    endif

                 else  ! it is the non-end interface, keep going along it

                    pt(1,n)=r_interface
                    raysec%diffracted =.true.
                    raysec%ray%diffracted =.true.
                    cycle ploop
                 
                 endif

              endif

              r_interface=interpolate_interface(pt(2,n),pt(3,n),ifbot)
              if (pt(1,n) <= r_interface) then

                 if (ifbot%id == ifend_id) then

                    if (section_is_refraction) then  

                       ! when the ray is refracted into the region check for total reflection (headwave),
                       ! if so keep going along the interface until the ray can get through

                       call interface_normal(pt(2,n),pt(3,n),intrface(ifend_id),norm(1),norm(2),norm(3),h)
                       call interpolate_time_gradient(tf,r_interface,pt(2,n),pt(3,n),grad1(1),grad1(2),grad1(3) &
                            ,verbose,outside,raysec)   
                       call interpolate_time_gradient(tf_next,r_interface,pt(2,n),pt(3,n),grad2(1),grad2(2),grad2(3) &
                            ,verbose,outside,raysec) 

                       if (outside) then ; n=n-1 ; exit ploop ; endif
                       if (.not.raysec%ray%valid) then; deallocate(pt) ; return ; endif

                       gtp_in =sqrt(sum((grad2 - dot_product(norm,grad2)*norm)**2))

                       if (gtp_in <= sqrt(sum(grad1**2))) then ! the ray can get through here

                          pp=pt(1:3,n)
                          exit ploop

                       else  ! total reflection, keep going along the interface

                          pt(1,n)=r_interface
                          raysec%headwave =.true.
                          raysec%ray%headwave =.true.   
                          cycle ploop

                       endif

                    else  ! in case the section is a reflection

                       pp=pt(1:3,n)
                       exit ploop

                    endif

                 else

                    pt(1,n)=r_interface
                    raysec%diffracted =.true.
                    raysec%ray%diffracted =.true.
                    cycle ploop
                 
                 endif

              endif
       
              ! if no interface is crossed, test for convergence, if sufficient exit
              if (sqrt(sum((pt(1:3,n)-pp(1:3))**2)) < p_tolerance) exit itloop

              ! safety stop for no convergence
              iter = iter+1
              if (iter > 10) stop 'trace_ray_section: no convergence after max iterations'

              ! reduce step size if initial size not accurate enough
              dl=0.5_dp*dl

           end do itloop


    ! special return condition if source is reached
           if(s%is_local .and. raysec%place_in_sequence == 1) then
              dx= sqrt((pt(1,n)-s%r)**2+(s%r*(pt(2,n)-s%lat))**2+(s%r*s%coslat*(pt(3,n)-s%long))**2)
              if (dx < s_tolerance) then
                 pt(1,n)=s%r ;pt(2,n)=s%lat ;pt(3,n)=s%long 
                 raysec%npoints=n
                 allocate(raysec%point(3,n))
                 raysec%point(1:3,1:n)=pt(1:3,1:n)
                 deallocate(pt)
                 return
              endif
           endif

        end do ploop

        if(s%is_local .and. raysec%place_in_sequence == 1)  stop ' unable to find local source in ray section 1'   

        if (outside) then

           raysec%npoints=n
           allocate(raysec%point(3,n))
           raysec%point(1:3,1:n)=pt(1:3,1:n)
           deallocate(pt)
           return

        endif

   ! find the intersection of the ray with the end interface

   x1=0.0_dp
   x2=dl
   f=pt(1,n-1) - interpolate_interface(pt(2,n-1),pt(3,n-1),ifend)
   fmid=pp(1)-r_interface


! check that interface is actually crossed, if not error stop

   if (f*fmid > 0.01_dp*pgrid%tolerance) then

      print *,raysec%place_in_sequence
!      print '(a10,4f14.7)','parms',dx,p_tolerance,s_tolerance
      print '(a10,4f14.7)','rt1',pt(1:3,n-1),f
      print '(a10,4f14.7)','rt2',pp(1:3),fmid
      print '(a10,4f14.7)','grad0',grad0(1:3),dl
      print *,'top if',interpolate_interface(pt(2,n-1),pt(3,n-1),intrface(tf%reg%itop%id))
      print *,'bot if',interpolate_interface(pt(2,n-1),pt(3,n-1),intrface(tf%reg%ibot%id))
      print *,'004 if',interpolate_interface(pt(2,n-1),pt(3,n-1),intrface(4))
!      do i=1,tf%reg%nnode
!         write(29,'(i5,4f12.6)') i,tf%arrivaltime(i),tf%time_gradient(1:3,i)
!         if (sqrt(sum(tf%time_gradient(1:3,i)**2)) < 0.01_dp) &
!              print'(i5,4f12.6,3i5)', i,tf%arrivaltime(i),tf%time_gradient(1:3,i),&
!              tf%reg%node(i)%i1,tf%reg%node(i)%i3,tf%reg%node(i)%i3
!      end do
      stop 'root not bracketed in trace_ray_section'
   endif


! catch pinched interfaces (end point = starting point)

   if (abs(f*fmid) <= 0.01_dp*pgrid%tolerance) then

      if ( ((raysec%iend%id == ifbot%id .and. fmid<0.0_dp) .or. &
         (raysec%iend%id == iftop%id .and. fmid>0.0_dp))) then


         pt(1:3,n) = pt(1:3,n-1)

         pt(1,n)=interpolate_interface(pt(2,n),pt(3,n),ifend)
!         print *,'no end point iteration necessary'
!         print '(a20,3f15.7)','last section point',pt(1:3,n)

         r_out=pt(1,n) ;lat_out=pt(2,n) ;long_out=pt(3,n) 

!         write(1,'(i5,6f12.5)') n,pt(1:3,n),grad0

         raysec%npoints=n
         allocate(raysec%point(3,n))
         raysec%point(1:3,1:n)=pt(1:3,1:n)

         deallocate(pt)

         return


      else

         print *,raysec%place_in_sequence
!         print '(a10,4f14.7)','parms',dx,p_tolerance,s_tolerance
         print '(a10,4f14.7)','rt1',pt(1:3,n-1),f
         print '(a10,4f14.7)','rt2',pp(1:3),fmid
         stop 'error traversing pinched interface'

      endif
   endif


! normal case of root bracketed

   if (f*fmid < -0.01_dp*pgrid%tolerance) then

      if (f < 0.0_dp) then
         rtbis=x1
         dx=x2-x1
      else
         rtbis=x2
         dx=x1-x2
      endif
      do nn=1,10
         dx=dx*0.5_dp
         xmid=rtbis+dx
         lat = pt(2,n-1)-xmid*grad0(2)/pt(1,n-1)
         long= pt(3,n-1)-xmid*grad0(3)/(pt(1,n-1)*cos(pt(2,n-1)))
         fmid=(pt(1,n-1)-xmid*grad0(1))-interpolate_interface(lat,long,ifend)
         if (fmid <= 0.0_dp) rtbis=xmid
      end do

      pt(1,n)=pt(1,n-1)-rtbis*grad0(1)
      pt(2,n)=pt(2,n-1)-rtbis*grad0(2)/pt(1,n-1)          
      pt(3,n)=pt(3,n-1)-rtbis*grad0(3)/(pt(1,n-1)*cos(pt(2,n-1)))

      pt(1,n)=interpolate_interface(pt(2,n),pt(3,n),ifend)

!      print '(a20,3f15.7)','last section point',pt(1:3,n)

   endif

   r_out=pt(1,n) ;lat_out=pt(2,n) ;long_out=pt(3,n) 

!   write(1,'(i5,6f12.5)') n,pt(1:3,n),grad0


   ! finally copy the ray section points from temporary array pt into the ray section
  raysec%npoints=n
  allocate(raysec%point(3,n))
  raysec%point(1:3,1:n)=pt(1:3,1:n)

  deallocate(pt)

  return

end subroutine trace_ray_section
!*************************************************************************************************************
!**********************************************************************************************************
! This routine returns the time gradient at an arbitray position within a time field
! In a regular cell it uses 
! In a cut cell this modified, could probably be improved for cut cells
subroutine interpolate_time_gradient(tf,r,lat,long,dtdr,dtdlat,dtdlong,verbose,outside,raysec)
  use mod_3dfm
  implicit none

  real(kind=dp)                              :: r,lat,long,dtdr,dtdlat,dtdlong
  type(Ttime_field)                          :: tf
  type(Tregion),pointer                      :: reg
  type(Tpropagation_grid),pointer            :: grid  
  type(Tray_section)                         :: raysec
  real(kind=dp)                              :: w(30),w1,w2,w3,tgrad(3,30)
  real(kind=dp)                              :: dxplus,dxmin,xbase,dx,avnorm(3),dplane,d(3)
  integer                                    :: nnode,i,j,k,m,ii,jj,kk,icell,ir,ilat,ilong,inode,n_inode,inode1
  logical    :: verbose,outside
  integer    :: locnode(20)
  real(kind=dp)                               :: norm_r,norm_lat,norm_long,h,nor,nolat,nolong
  real(kind=dp)                               :: interpolate_interface

  
  reg=> tf%reg
  grid => reg%grid

  ! find cell in which position resides

  ir    =  floor((r - grid%r0)/grid%dr0 + 1)
  ilat  =  floor((lat - grid%lat0)/grid%dlat0 + 1)
  ilong =  floor((long - grid%long0)/grid%dlong0 + 1)

  if (ir == grid%nr) ir=ir-1
  if (ir == 0) ir=ir+1

  outside = .false.
  if ((ir <1 .or. ir > grid%nr-1).or.(ilat <1 .or. ilat > grid%nlat-1).or.(ilong <1 .or. ilong > grid%nlong-1)) then
     print *,' outside grid in raytracing'
!     print *,ir,ilat,ilong
!     print *,grid%nr-1,grid%nlat-1,grid%nlong-1
!     print *,r,lat,long
     outside = .true.
     return
  endif



  nnode=0
  n_inode=0
  avnorm=0.0_dp
  locnode=0

  ! find the nodes that belong to this cell


  ! first check if the cell is cut by an interface

  if (associated(pgrid%ccind_from_3dc(ir,ilat,ilong)%p)) then

     ! check if the cell is cut by the top intersection

     ! icell is the index of the current connected cell in the list of cells cut by interface reg%itop.
                  
     icell = pgrid%ccind_from_3dc(ir,ilat,ilong)%p(reg%itop%iface_id)

     ! if icell == 0 the cell is not cut by the top interface
    
     if(icell /= 0) then
 
        ii=reg%itop%id
        do jj=1,intersection(ii)%n_inodes(icell)

!   m is the node number in the regional node list of node  jj in the list of inteface nodes that are part of cut cell icell

           inode =(intersection(ii)%inodes(jj,icell))
           m=intersection(ii)%rbel_node_id(inode)

           n_inode=n_inode+1
           if (n_inode == 1) inode1=m
           nnode=nnode+1
           locnode(nnode)=m
           tgrad(1:3,nnode)= tf%time_gradient(1:3,m) 


           select case (intersection(ii)%intype(inode))

              case(0)
                 w1=1.0_dp - abs(r-reg%r(m))/grid%dr0
                 w2=1.0_dp - abs(lat-reg%lat(m))/grid%dlat0
                 w3=1.0_dp - abs(long-reg%long(m))/grid%dlong0

              case(1)
                 xbase=grid%r0+grid%dr0*floor((reg%r(m)-grid%r0)/grid%dr0)
                 dxplus=(xbase+grid%dr0)-reg%r(m)
                 dxmin=reg%r(m)-xbase
                 dx=r-reg%r(m)
                 if (dx >= 0.0_dp)  w1=1.0_dp - abs(dx)/dxplus
                 if (dx < 0.0_dp)   w1=1.0_dp - abs(dx)/dxmin
                 w2=1.0_dp - abs(lat-reg%lat(m))/grid%dlat0
                 w3=1.0_dp - abs(long-reg%long(m))/grid%dlong0

              case(2)
                 xbase=grid%lat0+grid%dlat0*floor((reg%lat(m)-grid%lat0)/grid%dlat0)
                 dxplus=(xbase+grid%dlat0)-reg%lat(m)
                 dxmin=reg%lat(m)-xbase
                 dx=lat-reg%lat(m)
                 w1=1.0_dp - abs(r-reg%r(m))/grid%dr0
                 if (dx >= 0.0_dp)  w2=1.0_dp - abs(dx)/dxplus
                 if (dx < 0.0_dp)   w2=1.0_dp - abs(dx)/dxmin
                 w3=1.0_dp - abs(long-reg%long(m))/grid%dlong0

              case(3)
                 xbase=grid%long0+grid%dlong0*floor((reg%long(m)-grid%long0)/grid%dlong0)
                 dxplus=(xbase+grid%dlong0)-reg%long(m)
                 dxmin=reg%long(m)-xbase
                 dx=long-reg%long(m)
                 w1=1.0_dp - abs(r-reg%r(m))/grid%dr0
                 w2=1.0_dp - abs(lat-reg%lat(m))/grid%dlat0
                 if (dx >= 0.0_dp)  w3=1.0_dp - abs(dx)/dxplus
                 if (dx < 0.0_dp)   w3=1.0_dp - abs(dx)/dxmin

           end select


           avnorm=avnorm+intersection(ii)%normal(1:3,inode)
           w(nnode)=w1*w2*w3


!           if (intersection(ii)%intype(inode) /= 0) print *,'xbase',xbase,dxplus,dxmin

!  write(33,'(3i5,10f12.5)') ii,inode,intersection(ii)%intype(inode),reg%r(m),reg%lat(m),&
!           reg%long(m),xbase,dx,dxplus,dxmin,w1,w2,w3

        end do

     end if


! then check the bottom intersection

     icell = pgrid%ccind_from_3dc(ir,ilat,ilong)%p(reg%ibot%iface_id)
     if (icell /= 0) then

        ii=reg%ibot%id
        do jj=1,intersection(ii)%n_inodes(icell)
!   m is the node number in the regional node list of node  jj in the list of inteface nodes that are part of cut cell icell

           inode =(intersection(ii)%inodes(jj,icell))
           m=intersection(ii)%rabo_node_id(inode)

           n_inode=n_inode+1
           if (n_inode == 1) inode1=m
           nnode=nnode+1
           locnode(nnode)=m
           tgrad(1:3,nnode)= tf%time_gradient(1:3,m) 


           select case (intersection(ii)%intype(inode))

              case(0)
                 w1=1.0_dp - abs(r-reg%r(m))/grid%dr0
                 w2=1.0_dp - abs(lat-reg%lat(m))/grid%dlat0
                 w3=1.0_dp - abs(long-reg%long(m))/grid%dlong0

              case(1)
                 xbase=grid%r0+grid%dr0*floor((reg%r(m)-grid%r0)/grid%dr0)
                 dxplus=(xbase+grid%dr0)-reg%r(m)
                 dxmin=reg%r(m)-xbase
                 dx=r-reg%r(m)
                 if (dx >= 0.0_dp)  w1=1.0_dp - abs(dx)/dxplus
                 if (dx < 0.0_dp)   w1=1.0_dp - abs(dx)/dxmin
                 w2=1.0_dp - abs(lat-reg%lat(m))/grid%dlat0
                 w3=1.0_dp - abs(long-reg%long(m))/grid%dlong0

              case(2)
                 xbase=grid%lat0+grid%dlat0*floor((reg%lat(m)-grid%lat0)/grid%dlat0)
                 dxplus=(xbase+grid%dlat0)-reg%lat(m)
                 dxmin=reg%lat(m)-xbase
                 dx=lat-reg%lat(m)
                 w1=1.0_dp - abs(r-reg%r(m))/grid%dr0
                 if (dx >= 0.0_dp)  w2=1.0_dp - abs(dx)/dxplus
                 if (dx < 0.0_dp)   w2=1.0_dp - abs(dx)/dxmin
                 w3=1.0_dp - abs(long-reg%long(m))/grid%dlong0

              case(3)
                 xbase=grid%long0+grid%dlong0*floor((reg%long(m)-grid%long0)/grid%dlong0)
                 dxplus=(xbase+grid%dlong0)-reg%long(m)
                 dxmin=reg%long(m)-xbase
                 dx=long-reg%long(m)
                 w1=1.0_dp - abs(r-reg%r(m))/grid%dr0
                 w2=1.0_dp - abs(lat-reg%lat(m))/grid%dlat0
                 if (dx >= 0.0_dp)  w3=1.0_dp - abs(dx)/dxplus
                 if (dx < 0.0_dp)   w3=1.0_dp - abs(dx)/dxmin

           end select

           avnorm=avnorm+intersection(ii)%normal(1:3,inode)
           w(nnode)=w1*w2*w3

!  write(33,'(3i5,10f12.5)') ii,inode,intersection(ii)%intype(inode),reg%r(m),reg%lat(m),&
!           reg%long(m),xbase,dx,dxplus,dxmin,w1,w2,w3
        end do

     endif

  endif  ! if cell is cut


! then find the regular grid nodes of connected cell n

  do i=0,1
     ii=ir+i
     do j=0,1
        jj=ilat+j
        do k=0,1
           kk=ilong+k

           if (pgrid%node_region(ii,jj,kk) == reg%id) then  ! node has to belong to the current region

              m=pgrid%rnode_id(ii,jj,kk)

              nnode=nnode+1
              locnode(nnode)=m
              tgrad(1:3,nnode)= tf%time_gradient(1:3,m)  
              w1=1.0_dp - abs(r-reg%r(m))/grid%dr0
              w2=1.0_dp - abs(lat-reg%lat(m))/grid%dlat0
              w3=1.0_dp - abs(long-reg%long(m))/grid%dlong0
              w(nnode)=w1*w2*w3 

           endif


        end do
     end do
  end do

  if (n_inode > 500) then
  avnorm=avnorm/n_inode
     do i=1,n_inode
        if (intersection(reg%node(inode1)%i2)%intype(reg%node(inode1)%i3) /= 0) then
           d(1) = r-reg%r(inode1)
           d(2) = lat-reg%lat(inode1)*r
           d(3) = long-reg%long(inode1)*r*cos(lat)
           dplane= abs(dot_product(d,avnorm))/grid%dr0
           w(i)=w(i)/(dplane+0.5_dp)
        endif
     end do
  end if

  w(1:nnode)=w(1:nnode)/sum(w(1:nnode))

!  write(1,*) 'itgrad:'
!  do n=1,nnode
!     write(1,'(2i5,4f12.5,3i5)') n,n_inode,w(n),tgrad(1:3,n),ir,ilat,ilong
!  end do


  dtdr    = sum(tgrad(1,1:nnode)*w(1:nnode))  
  dtdlat  = sum(tgrad(2,1:nnode)*w(1:nnode))
  dtdlong = sum(tgrad(3,1:nnode)*w(1:nnode))

!  if (tf%id == 2) write(15,'(3i5,6f12.6)') ir,ilat,ilong,r,lat,long,dtdr,dtdlat,dtdlong

  if (sqrt(dtdr**2+dtdlat**2+dtdlong**2) < 0.01_dp) then
     print *
     print *,'----- WARNING -------'
     print *,'an interpolated time gradient was unphysically small '
     print *,'this occurred for ray',raysec%ray%raypath_id,' from source',raysec%source%id
     print *,'this sometimes happens when the specified source/path/receiver combination'
     print *,'is not physically consistent. The ray will be declared invalid'
     print *,'If you think the ray should exist please report this message as a bug'
     print *
     raysec%ray%valid = .false.

     return

     print *,'nodes:',nnode
     print *,'cell:',ir,ilat,ilong
     print *,'itop,ibot',tf%reg%itop%id,tf%reg%ibot%id
     if (.not.associated(pgrid%ccind_from_3dc(ir,ilat,ilong)%p)) then
        print *,'gradient error in uncut cell'
        print '(4f15.7)',r,lat,long,interpolate_interface(lat,long,intrface(4))
        do i=0,1
           ii=ir+i
           do j=0,1
              jj=ilat+j
              do k=0,1
                 kk=ilong+k

                 m=pgrid%rnode_id(ii,jj,kk)
                 print '(5i5,3f15.7)',m,ii,jj,kk,pgrid%node_region(ii,jj,kk),pgrid%r(ii),&
                      pgrid%lat(jj),pgrid%long(kk)

              end do
           end do
        end do
        stop
     else
        print *,'ccind',pgrid%ccind_from_3dc(ir,ilat,ilong)%p(1:n_interfaces)
     endif


     call interface_normal(lat,long,intrface(tf%inonstart%id),norm_r,norm_lat,norm_long,h)

     do i=1,nnode  
        m=locnode(i)
        call interface_normal(tf%reg%lat(m),tf%reg%long(m),intrface(tf%inonstart%iface_id),nor,nolat,nolong,h)
 

        print '(i5,4f12.6,3i5,2f12.6)',i,w(i),tgrad(1:3,i),tf%reg%node(m)%i1, &
             tf%reg%node(m)%i2,tf%reg%node(m)%i3, &
             (norm_r*tgrad(1,i)+norm_lat*tgrad(2,i)+norm_long*tgrad(3,i))/sqrt(sum(tgrad(1:3,i)**2)), &
             (nor*tgrad(1,i)+nolat*tgrad(2,i)+nolong*tgrad(3,i))/sqrt(sum(tgrad(1:3,i)**2))
!        print '(3f12.6)',nor,nolat,nolong
!        print *
     end do

     stop 'interpolate gradient error in cut cell'

  endif

  return

  end subroutine interpolate_time_gradient

!*********************************************************************************
!**********************************************************************************************************
! This subroutine estimates the arrival time at an arbitrary location within a time field
! It look for the closest regional node, and estimates the time to be
! t = t(node) + (x_node-x).grad(T)_node

subroutine interpolate_arrivaltime(tf,r,lat,long,atime)
  use mod_3dfm
  implicit none

  real(kind=dp)                              :: r,lat,long
  type(Ttime_field)                          :: tf
  type(Tregion),pointer                      :: reg
  type(Tintersection),pointer                :: isec
  type(Tpropagation_grid),pointer            :: grid  
  real(kind=dp)                              :: dist,dr(3),distmin,atime

  integer                                    :: node,i,j,k,m,ii,jj,kk,icell,ir,ilat,ilong,inode

  reg=> tf%reg
  grid => pgrid

  ! find cell in which position resides

  ir    =  floor((r - grid%r0)/grid%dr0 + 1)
  ilat  =  floor((lat - grid%lat0)/grid%dlat0 + 1)
  ilong =  floor((long - grid%long0)/grid%dlong0 + 1)

  if (ir == grid%nr) ir=ir-1
  if (ir == 0) ir=ir+1


  if ((ir <1 .or. ir > grid%nr-1).or.(ilat <1 .or. ilat > grid%nlat-1).or.(ilong <1 .or. ilong > grid%nlong-1)) then
     print *,' outside grid in raytracing'
     print *,ir,ilat,ilong
     print *,grid%nr-1,grid%nlat-1,grid%nlong-1
     print *,r,lat,long
  endif



  distmin=1.e100_dp
  node=0

  ! find the nodes that belong to this cell


  ! first check if the cell is cut by an interface

  if (associated(pgrid%ccind_from_3dc(ir,ilat,ilong)%p)) then

     ! check if the cell is cut by the top intersection

     ! icell is the index of the current connected cell in the list of cells cut by interface reg%itop.
                  
     icell = pgrid%ccind_from_3dc(ir,ilat,ilong)%p(reg%itop%iface_id)

     ! if icell == 0 the cell is not cut be the top interface
    
     if(icell /= 0) then
 
        isec => reg%itop
        do jj=1,isec%n_inodes(icell)

!   m is the node number in the regional node list of node  jj in the list of inteface nodes that are part of cut cell icell

           inode =(isec%inodes(jj,icell))
           m=isec%rbel_node_id(inode)

           dr(1)= (r-reg%r(m))
           dr(2)= (lat-reg%lat(m))*r
           dr(3)= (long-reg%long(m))*r*cos(lat)

           dist=sum(dr**2)
           if (dist < distmin) then
              node=m
              distmin=dist
           endif

        end do

     end if


! then check the bottom intersection

     icell = pgrid%ccind_from_3dc(ir,ilat,ilong)%p(reg%ibot%iface_id)
     if (icell /= 0) then

        isec => reg%ibot
        do jj=1,isec%n_inodes(icell)
!   m is the node number in the regional node list of node  jj in the list of inteface nodes that are part of cut cell icell

           inode =(isec%inodes(jj,icell))
           m=isec%rabo_node_id(inode)

           dr(1)= (r-reg%r(m))
           dr(2)= (lat-reg%lat(m))*r
           dr(3)= (long-reg%long(m))*r*cos(lat)

           dist=sum(dr**2)
           if (dist < distmin) then
              node=m
              distmin=dist
           endif

        end do

     endif

  endif  ! if cell is cut


! then find the regular grid nodes of connected cell n

  do i=0,1
     ii=ir+i
     do j=0,1
        jj=ilat+j
        do k=0,1
           kk=ilong+k

           if (pgrid%node_region(ii,jj,kk) == reg%id) then  ! node has to belong to the current region

              m=pgrid%rnode_id(ii,jj,kk)

              dr(1)= (r-reg%r(m))
              dr(2)= (lat-reg%lat(m))*r
              dr(3)= (long-reg%long(m))*r*cos(lat)

              dist=sum(dr**2)
              if (dist < distmin) then
                 node=m
                 distmin=dist
              endif


           endif


        end do
     end do
  end do

!  print *,'closest node is',reg%node(node)%i1,reg%node(node)%i2,reg%node(node)%i3
!  print *, 'distance', distmin, tf%arrivaltime(node)

  dr(1)= (r-reg%r(node))
  dr(2)= (lat-reg%lat(node))*r
  dr(3)= (long-reg%long(node))*r*cos(lat)

  atime=tf%arrivaltime(node)+dot_product(tf%time_gradient(1:3,node),dr)

  return

  end subroutine interpolate_arrivaltime

!*********************************************************************************
recarray.f90
program recarray

  nrec=119
  write(1,*),nrec 

  do i=1,nrec
     
     r=0.0
     xlat=1.0+i*0.2
     xlong=xlat
     write(1,'(3f6.2)')r,xlat,xlong 
     write(1,'(a1)')'5'
     write(1,'(a10)')'1 1 1 1 1'
     write(1,'(a10)')'1 2 3 4 5'

  end do

end program recarray
setbrn.f
      program setbrn
c
c Herewith the new version of setbran.f with separated table and header
c and correct index assignments for direct access (hopefully)
c
      include 'limits.inc'
      save
      character*8 code,phcd
      character*20 modnam
      double precision zm,pm,pb,pu,taup,xp,taul,px,xt,xl,pux,pt,taut,
     1 coef,xa
      double precision tmp(nsl1,2),xm(nsl1,2),deg,dtol,zmax,zoc,zic,
     1 z0
      dimension ndx2(nsr0,2)
      common/umodc/zm(nsr0,2),pm(nsr0,2),ndex(nsr0,2),mt(2)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/brnc/phcd(jbrn)
      common/msqc/pux(jbrn,2),km(2),midx(jbrn,2)
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data nin,nout,xmin,dtol/1,2,200.,1d-6/
      deg=180d0/3.1415927d0
c
c     write(6,*) "rec length for dasign:"
c     read(5,*) ndasr
c
      call assign(nin,-1,'remodl.hed')
      read(nin)ndasr,modnam,zmax,zoc,zic,kb,(pb(i),i=1,kb(2)),
     1 mt,lt,lbb,lcb,xn,pn,tn
      read(nin)((lbrk(i,nph),i=1,lbb(nph)),(code(i,nph),
     1 i=1,lcb(nph)),(zm(i,nph),pm(i,nph),ndex(i,nph),i=1,mt(nph)),
     2 (lvz(i,nph),taul(i,nph),xl(i,nph),i=1,lt(nph)),nph=1,2)
      call retrns(nin)
      print *,'ndasr =',ndasr,'  modnam = ',modnam
c
      call dasign(nin,-1,'remodl.tbl',ndasr)
      nrec=0
      do 1 nph=1,2
      n1=kb(nph)
      ind=0
      do 2 k=1,n1
 2    xm(k,nph)=0d0
 3    nrec=nrec+1
      read(nin,rec=nrec)z0,n,(tmp(k,1),k=1,n),(tmp(k,2),k=1,n)
      if(ind.gt.0) go to 4
      if(dabs(z0-zoc).le.dtol) go to 4
      j=1
      do 5 k=2,n
      xm(k,nph)=dmax1(xm(k,nph),dabs(tmp(j,2)-tmp(k,2)))
 5    j=k
      if(n+1.eq.n1) xm(n1,nph)=tmp(n,2)
      go to 3
 4    ind=ind+1
      do 6 k=1,n
      taup(k,ind,nph)=tmp(k,1)
 6    xp(k,ind,nph)=tmp(k,2)
      if(ind.lt.3) go to 3
 1    continue
c
      xmin=xn*xmin
c
c^      call assign(10,2,'setbrn1.lis')
c^      write(10,*)'kb mt lt lbb lcb',kb,mt,lt,lbb,lcb
c^      write(10,*)'xn pn tn xmin',xn,pn,tn,xmin
      cn=1./xn
c^      write(10,209)(i,(lbrk(i,j),code(i,j),j=1,2),i=1,lbb(1))
 209  format(/(1x,2i5,2x,a,i5,2x,a))
c^      write(10,210)(i,lbrk(i,2),code(i,2),i=lbb(1)+1,lbb(2))
 210  format(1x,i5,15x,i5,2x,a)
c^      write(10,200,iostat=ios)(i,(zm(i,j),pm(i,j),ndex(i,j),j=1,2),
c^     1 i=1,mt(1))
 200  format(/(1x,i5,2f12.6,i5,2x,2f12.6,i5))
c^      write(10,201,iostat=ios)(i,zm(i,2),pm(i,2),ndex(i,2),
c^     1 i=mt(1)+1,mt(2))
 201  format(1x,i5,31x,2f12.6,i5)
c^      write(10,217)((nph,i,lvz(i,nph),taul(i,nph),deg*xl(i,nph),
c^     1 i=1,lt(nph)),nph=1,2)
 217  format(/(1x,3i5,f12.6,f12.2))
c^      write(10,202)(i,pb(i),cn*xm(i,1),cn*xm(i,2),i=1,kb(1))
 202  format(/(5x,i5,f12.6,2f12.2))
c^      write(10,203)(i,pb(i),cn*xm(i,2),i=kb(1)+1,kb(2))
 203  format(5x,i5,f12.6,12x,f12.2)
c^      call retrns(10)
c^      call assign(10,2,'setbrn2.lis')
c
      do 8 nph=1,2
      n1=kb(nph)
      do 9 i=2,n1
      xm(i,nph)=xm(i-1,nph)+xm(i,nph)
      pu(i,nph)=pb(i)
 9    kuse(i,nph)=-1
      do 8 ind=3,2,-1
      jnd=ind-1
      do 8 i=1,n1
      taup(i,ind,nph)=taup(i,ind,nph)-taup(i,jnd,nph)
 8    xp(i,ind,nph)=xp(i,ind,nph)-xp(i,jnd,nph)
      do 10 nph=1,2
 10   call pdecx(kb(nph),nph,xm,2.)
c
c^      write(10,*)'ku',ku
c^      write(10,204)(i,(pu(i,nph),cn*xm(i,nph),cn*(xm(i+1,nph)-
c^     1 xm(i,nph)),nph=1,2),i=1,ku(1))
 204  format(/(5x,i5,2(f12.6,2f12.2)))
c^      write(10,205)(i,pu(i,2),cn*xm(i,2),cn*(xm(i+1,2)-xm(i,2)),
c^     1 i=ku(1)+1,ku(2))
 205  format(5x,i5,36x,f12.6,2f12.2)
c^      do 207 nph=1,2
c^ 207  write(10,206)(i,pb(i),(taup(i,j,nph),j=1,3),(deg*xp(i,j,nph),
c^     1 j=1,3),i=1,kb(nph))
c^ 206  format(/(1x,i5,4f10.6,3f10.2))
c
      call layout
c     write(10,214)(i,pb(i),(kuse(i,j),j=1,2),i=1,kb(2))
c214  format(/(5x,i5,f12.6,2i5))
      do 11 nph=1,2
      n1=kb(nph)
      k=0
      do 12 i=1,n1
      if(kuse(i,nph).lt.0) go to 12
      k=k+1
      pu(k,nph)=pb(i)
 12   continue
 11   ku(nph)=k
      call kseq
      call mseq
c^      write(10,215)(i,(pu(i,j),j=1,2),i=1,ku(1))
 215  format(/(5x,i5,2f12.6))
c^      write(10,216)(i,pu(i,2),i=ku(1)+1,ku(2))
 216  format(5x,i5,12x,f12.6)
c^      write(10,208)(i,(nafl(i,j),j=1,3),(indx(i,j),j=1,2),(kndx(i,j),
c^     1 j=1,2),(fcs(i,j),j=1,3),i=1,nseg)
 208  format(/(1x,8i6,3f6.1))
c^      write(10,211)(i,(jndx(i,j),j=1,2),(mndx(i,j),j=1,2),(px(i,j),
c^     1 j=1,2),(deg*xt(i,j),j=1,2),phcd(i),i=1,nbrn)
 211  format(/(1x,i3,4i5,2f12.6,2f10.2,2x,a))
c^      write(10,218)(i,(midx(i,j),j=1,2),(pux(i,j),j=1,2),
c^     1 i=1,max0(km(1),km(2)))
 218  format(/(1x,i3,2i5,2f12.6))
c^      write(10,212,iostat=ios)(i,pt(i),taut(i),deg*xa(i),
c^     1 cn*(xa(i)-xa(i+1)),(coef(j,i),j=1,5),i=1,nl)
 212  format(/(1x,i4,0p2f12.6,2f10.2,1p5d10.2))
c^      call retrns(10)
c
c^      call assign(10,2,'setbrn3.lis')
      do 20 nph=1,2
      mt(nph)=mt(nph)-3
      ku(nph)=ku(nph)-1
 20   km(nph)=km(nph)-1
c     icor=33  -  originally 32 records used as header in setbrn
c                 and 2 records used as header in remodl.
      icor=3
      do 14 nph=1,2
      m1=mt(nph)
      icor=icor-3
      do 14 i=2,m1
 14   ndx2(i,nph)=ndex(i,nph)+icor
      len1=ku(2)+km(2)
      len0=8*len1
      len2=5*nl
c^      write(10,*)'nseg nbrn mt ku km len len1',nseg,nbrn,mt,ku,km,len0,
     1 len1
c^      write(10,*)
      nasgr = len0
      write(6,*) 'reclength for direct access', nasgr
c++ 
c     write(6,*) 'enter model name'
c     read(5,*) cmodel
      nb=index(modnam,' ')-1
      if(nb.le.0) nb=len(modnam)
c     cnam1 = cmodel(1:nb)//'.tbl'
c     cnam2 = cmodel(1:nb)//'.hed'
      write(6,*) 'header file  :',modnam(1:nb)//'.hed'
      write(6,*) 'table file   :',modnam(1:nb)//'.tbl'
      call assign(nout,-2,modnam(1:nb)//'.hed')
c++ 
      write(nout) nasgr,nl,len2,xn,pn,tn,mt,nseg,nbrn,ku,km,fcs,nafl,
     1 indx,kndx
      write(nout) pm,zm,ndx2
      write(nout) pu,pux
      write(nout) phcd,px,xt,jndx
      write(nout) pt,taut
      write(nout) coef
      call retrns(nout)
c
      call dasign(nout,-2,modnam(1:nb)//'.tbl',nasgr)
      nrec = 0
      do 16 nph=1,2
      m1=mt(nph)
      n1=ku(nph)
      k1=km(nph)
c^      write(10,*)'nph m1 n1 k1',nph,m1,n1,k1
      do 16 m=2,m1
      if(ndex(m,nph).eq.ndex(m-1,nph)) go to 16
      read(nin,rec=ndex(m,nph))z0,n,(tmp(k,1),k=1,n),(tmp(k,2),k=1,n)
c^      write(10,*)'m nph ndex n',m,nph,ndex(m,nph),n
      k=0
      l=1
      do 17 i=1,n
      if(kuse(i,nph).lt.0) go to 17
      if(dabs(pux(l,nph)-pb(i)).gt.dtol) go to 18
      tmp(l,2)=tmp(i,2)
      l=l+1
 18   k=k+1
      tmp(k,1)=tmp(i,1)
 17   continue
c^      write(10,*)'k l nrec',k,l-1,nrec+1,ndx2(m,nph),sngl(tmp(1,1))
      if(k.ge.n1) go to 19
      k=k+1
      do 21 i=k,n1
 21   tmp(i,1)=0d0
 19   if(l.gt.k1) go to 23
      do 22 i=l,k1
 22   tmp(i,2)=0d0
 23   nrec=nrec+1
      write(nout,rec=nrec)(tmp(i,1),i=1,n1),(tmp(i,2),i=1,k1)
 16   continue
c
      call retrns(10)
      call retrns(nin)
      call retrns(nout)
      call vexit(0)
      end
c
      subroutine pdecx(n1,nph,xm,fac)
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl
      double precision xm(nsl1,2),ptol,pa,pax,plim
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      data ptol/.03/
c
      call collct(1,n1,xm(1,nph),fac*xmin)
      k=0
      plim=.7d0*pu(n1,nph)
      do 1 i=1,n1
      if(xm(i,nph).lt.0d0) go to 1
      if(pu(i,nph).lt.plim) go to 2
      if(pu(i,nph)-pu(k,nph).le.ptol) go to 2
      pa=pu(k,nph)+.75d0*(pu(i,nph)-pu(k,nph))
      pax=1d10
      m=0
      do 3 j=i1,i
      if(dabs(pu(j,nph)-pa).ge.pax) go to 3
      m=j
      pax=dabs(pu(j,nph)-pa)
 3    continue
      if(m.eq.i1.or.m.eq.i) go to 2
      k=k+1
      pu(k,nph)=pu(m,nph)
      xm(k,nph)=0d0
      kuse(m,nph)=1
 2    k=k+1
      i1=i
      pu(k,nph)=pu(i,nph)
      xm(k,nph)=xm(i,nph)
      kuse(i,nph)=1
 1    continue
      ku(nph)=k
      return
      end
c
      subroutine collct(i1,i2,x,xmn)
c
c $$$$$ calls varn $$$$$
c
      save
      double precision x(i2)
      data cn/6371./
c
      is=i1+1
      ie=i2-1
      if(ie.lt.is) return
      k1=i1
      var=0.
      m=0
      do 1 i=is,ie
      dx1=dabs(x(k1)-x(i))-xmn
      dx2=dabs(x(k1)-x(i+1))-xmn
      if(abs(dx2).ge.abs(dx1)) go to 2
      x(i)=-x(i)
      go to 1
 2    if(k1.le.i1) kb=i
      k1=i
      var=var+dx1*dx1
      m=m+1
 1    continue
      dx1=dabs(x(k1)-x(i2))-xmn
      var=var+dx1*dx1
      m=m+1
 7    if(m.le.1) return
      k1=i1
      k2=kb
      ks=kb+1
      nch=0
      do 8 i=ks,i2
      if(x(i).lt.0d0) go to 8
      k0=k1
      k1=k2
      k2=i
      var1=varn(x,k0,k1,k2,k1-1,xmn,var,m,m1)
      var2=varn(x,k0,k1,k2,k1+1,xmn,var,m,m2)
      if(amin1(var1/m1,var2/m2).ge.var/m) go to 6
      nch=nch+1
      x(k1)=-x(k1)
      if(var1/m1-var2/m2)3,4,5
 4    if(m1-m2)3,3,5
 3    k1=k1-1
      x(k1)=dabs(x(k1))
      var=var1
      m=m1
      go to 6
 5    k1=k1+1
      x(k1)=dabs(x(k1))
      var=var2
      m=m2
 6    if(k0.eq.i1) kb=k1
 8    continue
      if(nch.gt.0) go to 7
      return
      end
c
      function varn(x,k0,k1,k2,kt,xmn,var,m,mn)
c
c $$$$$ calls only library routines $$$$$
c
      save
      double precision x(k2)
c
      dx1=dabs(x(k0)-x(k1))-xmn
      dx2=dabs(x(k1)-x(k2))-xmn
      varn=var-dx1*dx1-dx2*dx2
      if(kt.le.k0.or.kt.ge.k2) go to 1
      dx1=dabs(x(k0)-dabs(x(kt)))-xmn
      dx2=dabs(dabs(x(kt))-x(k2))-xmn
      varn=varn+dx1*dx1+dx2*dx2
      mn=m
      return
 1    dx1=dabs(x(k0)-dabs(x(k2)))-xmn
      varn=varn+dx1*dx1
      mn=m-1
      return
      end
c
      subroutine layout
c
c   Layout contains the program for the desired travel-time segments 
c   implemented as calls to the mk_br entry points.  Each call does one 
c   segment (which may have many branches).
c
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision dir(3),cref(3),sref(3)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data dir,cref,sref/1d0,1d0,1d0,1d0,2d0,2d0,2d0,2d0,2d0/
c
c   Initialize variables.
      nseg=0
      nbrn=0
      nl=0
      do 1 j=1,3
      do 1 i=1,jseg
 1    fcs(i,j)=0.
      do 2 i=1,jout
 2    taut(i)=0d0
      do 3 j=1,2
      do 3 i=1,jbrn
 3    xt(i,j)=0d0
c
c   Do all of the segments.
c
c   P (up-going branch).
      print *,'Layout:  do Pup'
      call mkubr(ku(1),   +1)
c   P, Pdiff, and PKP.
      print *,'Layout:  do P and PKP'
      call mkdbr(1,lbb(1),-1,3,1,1,dir)
c   PKiKP.
      print *,'Layout:  do PKiKP'
      call mkrbr(2,       -1,2,1,1,dir)
c   pP.
      print *,'Layout:  do pP'
      call mkdbr(1,lbb(1),+1,3,1,1,dir)
c   sP.
      print *,'Layout:  do sP'
      call mkdbr(1,lbb(1),+2,3,1,1,dir)
c   pPKiKP.
      print *,'Layout:  do pPKiKP'
      call mkrbr(2,       +1,2,1,1,dir)
c   sPKiKP.
      print *,'Layout:  do sPKiKP'
      call mkrbr(2,       +2,2,1,1,dir)
c   PcP.
      print *,'Layout:  do PcP'
      call mkrbr(3,       -1,1,1,1,dir)
c   ScP.
      print *,'Layout:  do ScP'
      call mkrbr(3,       -2,1,2,1,dir)
c   SKP.
      print *,'Layout:  do SKP'
      call mkdbr(1,3,     -2,3,2,1,dir)
c   SKiKP.
      print *,'Layout:  do SKiKP'
      call mkrbr(2,       -2,2,2,1,dir)
c   PKKP.
      print *,'Layout:  do PKKP'
      call mkdbr(1,3,     -1,3,1,1,cref)
c   SKKP.
      print *,'Layout:  do SKKP'
      call mkdbr(1,3,     -2,3,2,1,cref)
c   PP and P'P'.
      print *,'Layout:  do PP, P''P'''
      call mkdbr(1,lbb(1),-1,3,1,1,sref)
c   S (up-going branch).
      print *,'Layout:  do Sup'
      call mkubr(ku(2),   +2)
c   S, Sdiff, and SKS.
      print *,'Layout:  do S and SKS'
      call mkdbr(1,lbb(2),-2,3,2,2,dir)
c   pS
      print *,'Layout:  do pS'
      call mkdbr(1,lbb(1),+1,3,2,2,dir)
c   sS
      print *,'Layout:  do sS'
      call mkdbr(1,lbb(2),+2,3,2,2,dir)
c   ScS
      print *,'Layout:  do ScS'
      call mkrbr(4,       -2,1,2,2,dir)
c   PcS
      print *,'Layout:  do PcS'
      call mkrbr(3,       -1,1,1,2,dir)
c   PKS
      print *,'Layout:  do PKS'
      call mkdbr(1,3,     -1,3,1,2,dir)
c   PKKS
      print *,'Layout:  do PKKS'
      call mkdbr(1,3,     -1,3,1,2,cref)
c   SKKS
      print *,'Layout:  do SKKS'
      call mkdbr(1,3,     -2,3,2,2,cref)
c   SS and S'S'.
      print *,'Layout:  do SS and S''S'''
      call mkdbr(1,lbb(2),-2,3,2,2,sref)
c   SP
      print *,'Layout:  do SP'
      call mkcbr(4,lbb(1),-2,1,2,1,sref)
c   PS
      print *,'Layout:  do PS'
      call mkcbr(4,lbb(1),-1,1,1,2,sref)
      return
      end
c
      subroutine mkdbr(l1,l2,isgn,lyr,nph,kph,fac)
c
c   Mkdbr sets up a simple refracted wave segment.  L1 and l2 point to the 
c   lbrk array of slowness break point pointers.  Note that the P and S 
c   break point arrays don't necessarily line up layer by layer.  This is 
c   not generally a problem as most phases need only worry about the 
c   pointer to the surface slowness for one wave type and a pointer 
c   somewhere in the core (which is constrained to be the same for both 
c   P and S).  Isgn is positive if the wave starts out going up and 
c   negative if the wave starts out going down.  Iabs(isng) is 1 if the 
c   wave starts out as a P wave and 2 if the wave starts out as an S wave.  
c   Lyr gives the number of major layers (mantle, outer core, and inner 
c   core) that the wave penetrates.  Nph and kph give the wave type (1 for 
c   P and 2 for S) on the down-going and up-going legs of the ray path 
c   respectively.  Fac is a three element array giving the number of 
c   repeats of the ray path in each major layer.  This scheme incorporates 
c   turning rays (e.g., P and S), turning rays reflected, but not 
c   converted at the surface (e.g., PP and SS), up-going rays reflected 
c   and/or converted at the surface into turning rays (e.g., pP and sP), 
c   turning rays converted during transmission through an interface (e.g., 
c   SKP and PKS), and rays which turn multiple times while reflecting from 
c   the bottom side of a layer (e.g., PKKP or SKKP).  Mkdbr does not 
c   include up-going rays (to the receiver), rays reflected from the top 
c   side of a discontinuity, or rays which are reflected and converted at 
c   the free surface.  See mkubr, mkrbr, and mkcbr respectively for 
c   routines which handle these types of rays.
c
      save
      include 'limits.inc'
      character*8 code,phcd,ks
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision fac(3)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/brnc/phcd(jbrn)
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data ks/'KKKKKKKK'/
c
c   Remember the programming as part of the final phase construction is 
c   done in depcor.
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
c   Using l1 and l2 to get the breakpoints has some shortcommings, 
c   particularly for converted phases.  It would be more general to 
c   have separate indicies for the breakpoints and the layers covered.
      if(l1.gt.1) kndx(nseg,1)=lbrk(l1-1,nph)
      kndx(nseg,2)=min0(min0(lbrk(l2,nph),lbrk(l2,kph)),
     1 lbrk(l2,iabs(isgn)))
      print *,'Mkdbr:  l1 l2 isgn lyr nph kph =',l1,l2,isgn,lyr,nph,kph
      print *,'Mkdbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 14 m=1,lyr
      fcs(nseg,m)=fac(m)
 14   xfc=amax1(xfc,fcs(nseg,m))
c
c   Set up the required slownesses, taus and distances.
c
      j=kndx(nseg,1)
      lz1=1
      lz2=1
c   Loop over the layers of interest.
      do 1 i=l1,l2
c   Be sure that the phase cuts off at the right place.
      l=min0(lbrk(i,nph),kndx(nseg,2))
c   Skip all total internal reflections.
      if(code(i,nph)(1:1).eq.'r'.or.j.ge.l) go to 1
c   Set the starting branch pointer.
      nbrn=nbrn+1
      nt=nl+1
      jndx(nbrn,1)=nt
c   Copy in the desired slownesses.
      do 2 k=j,l
      nl=nl+1
      pt(nl)=pb(k)
c   Add up the tau contributions.
      do 2 m=1,lyr
 2    taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
c   Take care of branch end pointers and slownesses.
      mndx(nbrn,1)=j
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(j)
      px(nbrn,2)=pb(l)
c   Add up distance contributions for the branch end points only.
      do 3 m=1,lyr
      xt(nbrn,1)=xt(nbrn,1)+fac(m)*(xp(j,m,nph)+xp(j,m,kph))
 3    xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
c   Take care of the contribution due to low velocity zones for the 
c   down-going leg(s).
      if(j.ne.lvz(lz1,nph)) go to 9
      do 11 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,nph)
 11   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,nph)
      taut(nt)=taut(nt)+fac(1)*taul(lz1,nph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz1,nph)
      lz1=lz1+1
c   Take care of the contributions due to low velocity zones for the 
c   up-going leg(s).
 9    if(j.ne.lvz(lz2,kph)) go to 10
      do 12 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,kph)
 12   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,kph)
      taut(nt)=taut(nt)+fac(1)*taul(lz2,kph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz2,kph)
      lz2=lz2+1
c   Decimate the slownesses if the branch is oversampled in distance.
 10   call pdect(jndx(nbrn,1),nl,j,iabs(isgn),xfc)
c   Set up the interpolation.
      call tauspl(jndx(nbrn,1),nl,pt,coef)
c   Remember the final branch end slowness value.
      jndx(nbrn,2)=nl
c
c   Take care of the branch name.  First, set up a default.
      phcd(nbrn)=code(i,nph)(2:2)//code(i,kph)(3:)
      if(idint(fac(1)+.5d0).gt.1) go to 5
      if(idint(fac(2)+.5d0).le.1) go to 4
c   Re-do the name if the ray is reflected from the underside of the 
c   core-mantle boundary.
      ind=idint(fac(2)-.5d0)
      phcd(nbrn)=code(i,nph)(2:2)//ks(1:ind)//code(i,kph)(3:)
      go to 4
c   Re-do the name if the ray is reflected from the surface.
 5    if(code(i,nph)(3:3).eq.' ') phcd(nbrn)=code(i,nph)(2:2)//
     1 code(i,kph)(2:)
      if(code(i,nph)(3:3).ne.' '.and.code(i,nph)(3:3).ne.'K')
     1  phcd(nbrn)=code(i,nph)(2:3)//code(i,kph)(2:)
      if(code(i,nph)(3:3).eq.'K') phcd(nbrn)=code(i,nph)(2:2)//''''//
     1 code(i,kph)(2:2)//''''//code(i,kph)(5:)
c   Take care .
 4    ind=max0(index(phcd(nbrn),'KSab'),index(phcd(nbrn),'S''ab'))
      if(phcd(nbrn)(1:1).eq.'S'.and.ind.gt.0) phcd(nbrn)(ind+2:ind+3)=
     1 'ac'
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
 1    j=l
      indx(nseg,2)=nl
      return
c
c   Mkubr handles up-going P and S.  L1 and isgn are as for mkdbr (except 
c   that l1 actually plays the role of l2 with the beginning break point 
c   assumed to be zero).  The other arguments are not needed.
c
      entry mkubr(l1,isgn)
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=0
      nafl(nseg,3)=0
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      l=kb(iabs(isgn))
      kndx(nseg,2)=l
      print *,'Mkubr:  l1 isgn =',l1,isgn
      print *,'Mkubr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      nbrn=nbrn+1
      jndx(nbrn,1)=nl+1
      do 6 k=1,l1
      nl=nl+1
      pt(nl)=pu(k,isgn)
 6    xa(nl)=0d0
      mndx(nbrn,1)=1
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(1)
      px(nbrn,2)=pb(l)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
      phcd(nbrn)=code(1,iabs(isgn))(2:2)
      indx(nseg,2)=nl
      return
c
c   Mkrbr handles reflected phases possibly with a conversion such as 
c   PcP, PcS, and PkiKP.  Arguments are as for mkdbr (except that l1 
c   actually plays the role of l2 with the beginning break point assumed 
c   to be zero).
c
      entry mkrbr(l1,isgn,lyr,nph,kph,fac)
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      l=min0(lbrk(l1,nph),lbrk(l1,kph))
      kndx(nseg,2)=l
      print *,'Mkrbr:  l1 isgn lyr nph kph =',l1,isgn,lyr,nph,kph
      print *,'Mkrbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 15 m=1,lyr
      fcs(nseg,m)=fac(m)
 15   xfc=amax1(xfc,fcs(nseg,m))
      if(lyr.ge.2) xfc=2.
c
      nbrn=nbrn+1
      jndx(nbrn,1)=nl+1
      do 7 k=1,l
      nl=nl+1
      pt(nl)=pb(k)
      do 7 m=1,lyr
 7    taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
      mndx(nbrn,1)=1
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(1)
      px(nbrn,2)=pb(l)
      do 8 m=1,lyr
 8    xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
      call pdect(jndx(nbrn,1),nl,1,iabs(isgn),xfc)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
      if(lyr.eq.1) phcd(nbrn)=code(l1,nph)(2:2)//'c'//code(l1,kph)(2:2)
      if(lyr.eq.2) phcd(nbrn)=code(l1,nph)(2:2)//code(l1,kph)(3:)
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
      indx(nseg,2)=nl
      return
c
c   Mkcbr handles phases reflected and converted at the surface such as 
c   PS and SP.  Arguments are as for mkdbr.
c
      entry mkcbr(l1,l2,isgn,lyr,nph,kph,fac)
      if(nph.gt.0.and.kph.gt.0.and.nph.ne.kph) go to 29
      print *,'Mkcbr:  bad call - nph kph =',nph,kph
      call vexit(1)
 29   nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      if(l1.gt.1) kndx(nseg,1)=min0(lbrk(l1,nph),lbrk(l1,kph))
      kndx(nseg,2)=min0(min0(lbrk(l2,nph),lbrk(l2,kph)),
     1 lbrk(l2,iabs(isgn)))
      print *,'Mkcbr:  l1 l2 isgn lyr nph kph =',l1,l2,isgn,lyr,nph,kph
      print *,'Mkcbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 16 m=1,lyr
      fcs(nseg,m)=fac(m)
 16   xfc=amax1(xfc,fcs(nseg,m))
c
      j=kndx(nseg,1)
      lz1=1
      lz2=1
      ik=l1
c
      print *,'Mkcbr:  start loop'
      do 17 in=l1,l2
 31   l=min0(lbrk(in,nph),kndx(nseg,2))
      if(code(in,nph)(1:1).eq.'r'.or.j.ge.l) go to 17
      l=min0(lbrk(ik,kph),kndx(nseg,2))
      if(code(ik,kph)(1:1).ne.'r'.and.j.lt.l.or.ik.ge.l2) go to 28
      j=max0(j,l)
      ik=ik+1
      go to 31
c
 28   if(lbrk(in,nph).le.lbrk(ik,kph)) go to 26
      l=min0(lbrk(ik,kph),kndx(nseg,2))
      print *,'kph:  kph ik j l code =',kph,ik,j,l,' ',code(ik,kph)
      isw=2
      go to 27
 26   l=min0(lbrk(in,nph),kndx(nseg,2))
      print *,'nph:  nph in j l code =',nph,in,j,l,' ',code(in,nph)
      isw=1
c
 27   nbrn=nbrn+1
      nt=nl+1
      jndx(nbrn,1)=nt
      do 18 k=j,l
      nl=nl+1
      pt(nl)=pb(k)
      do 18 m=1,lyr
 18   taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
      mndx(nbrn,1)=j
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(j)
      px(nbrn,2)=pb(l)
      do 19 m=1,lyr
      xt(nbrn,1)=xt(nbrn,1)+fac(m)*(xp(j,m,nph)+xp(j,m,kph))
 19   xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
      if(j.ne.lvz(lz1,nph)) go to 20
      do 21 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,nph)
 21   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,nph)
      taut(nt)=taut(nt)+fac(1)*taul(lz1,nph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz1,nph)
      lz1=lz1+1
 20   if(j.ne.lvz(lz2,kph)) go to 22
      do 23 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,kph)
 23   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,kph)
      taut(nt)=taut(nt)+fac(1)*taul(lz2,kph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz2,kph)
      lz2=lz2+1
 22   call pdect(jndx(nbrn,1),nl,j,iabs(isgn),xfc)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
c
      if(code(in,nph)(3:3).eq.' ') phcd(nbrn)=code(in,nph)(2:2)//
     1 code(ik,kph)(2:)
      if(code(in,nph)(3:3).ne.' '.and.code(in,nph)(3:3).ne.'K')
     1  phcd(nbrn)=code(in,nph)(2:3)//code(ik,kph)(2:)
      if(code(in,nph)(3:3).eq.'K') phcd(nbrn)=code(in,nph)(2:2)//''''//
     1 code(ik,kph)(2:2)//''''//code(ik,kph)(5:)
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
      print *,'phcd:  in ik phcd =',in,ik,' ',phcd(nbrn)
      if(isw.le.1) go to 17
      ik=ik+1
      j=max0(j,l)
      go to 31
 17   j=max0(j,l)
      indx(nseg,2)=nl
      return
      end
c
      subroutine pdect(i1,i2,j1,nph,fac)
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision h1,h2,hh
      dimension ib(2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
c
      xmn=fac*xmin
      isg=1
      do 1 i=1,2
      ib(i,1)=i1
 1    ib(i,2)=i2
      ii=i1+1
      ie=i2-1
      xa(i1)=xt(nbrn,1)
      do 2 i=ii,ie
      h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
 2    xa(i)=-(h2*taut(i-1)-(h2+h1)*taut(i)+h1*taut(i+1))/hh
      xa(i2)=xt(nbrn,2)
      do 3 i=ii,ie
      if((xa(i+1)-xa(i))*(xa(i)-xa(i-1)).gt.0d0) go to 3
      isg=2
      ib(1,2)=i-2
      ib(2,1)=i+2
 3    continue
      do 4 it=1,isg
 4    call collct(ib(it,1),ib(it,2),xa,xmn)
      k=i1-1
      j=j1
      do 5 i=i1,i2
      if(xa(i).lt.0d0) go to 5
      k=k+1
      pt(k)=pt(i)
      taut(k)=taut(i)
      xa(k)=xa(i)
      kuse(j,nph)=1
 5    j=j+1
      if(k.eq.nl) return
      ii=k+1
      do 6 i=ii,nl
 6    taut(i)=0d0
      nl=k
      return
      end
c
      subroutine kseq
c
c   Kseq makes a correspondence between model slownesses in array pb and 
c   the subset of the same slownesses used for sampling tau which are 
c   stored in pu (separate sets for P and S).  The net result is to 
c   translate the kndx pointers to critical slowness values (bounding 
c   branches actually implemented) from pointing into pb to pointing 
c   into pu.
c   
c
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl
      dimension kl(2),kk(jseg,2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      data kl/0,0/
c
c   Compile a sorted list of unique kndx values in the first column of 
c   kk.
c
      do 1 i=1,nseg
      nph=iabs(nafl(i,1))
      k=kl(nph)
      do 2 j=1,2
      if(k.le.0) go to 3
      do 4 m=1,k
      if(kk(m,nph,1)-kndx(i,j))4,2,5
 4    continue
 3    k=k+1
      kk(k,nph,1)=kndx(i,j)
      go to 2
 5    do 6 l=k,m,-1
 6    kk(l+1,nph,1)=kk(l,nph,1)
      k=k+1
      kk(m,nph,1)=kndx(i,j)
 2    continue
 1    kl(nph)=k
c
c   Make the correspondence between pb and pu for each kndx and save it 
c   in the second column of kk.
c
      do 7 nph=1,2
      n1=ku(nph)
      k=1
      ki=kk(k,nph,1)
      do 8 i=1,n1
      if(pu(i,nph)-pb(ki))8,9,10
 9    kk(k,nph,2)=i
      if(k.ge.kl(nph)) go to 7
      k=k+1
      ki=kk(k,nph,1)
 8    continue
 10   write(*,100)ki,pb(ki),nph
 100  format(1x,'Kseq:  pb(',i3,') =',f7.4,' not found in pu(*,',i1,
     1 ').')
      call vexit(1)
 7    continue
c
c   Replace each kndx pb pointer with the corresponding pu pointer.
c
      do 11 i=1,nseg
      nph=iabs(nafl(i,1))
      k=kl(nph)
      do 11 j=1,2
      do 12 m=1,k
      if(kk(m,nph,1)-kndx(i,j))12,11,13
 12   continue
 13   write(*,101)kndx(i,j)
 101  format(1x,'Kseq:  kndx value',i4,' not translated.')
      call vexit(1)
 11   kndx(i,j)=kk(m,nph,2)
      return
      end
c
      subroutine mseq
c         partial reordering of tables
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pux
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/msqc/pux(jbrn,2),km(2),midx(jbrn,2)
      data km/0,0/
c
      is=1
      do 1 i=1,nbrn
 8    if(jndx(i,2).le.indx(is,2)) go to 7
      is=is+1
      go to 8
 7    nph=iabs(nafl(is,1))
      k=km(nph)
      do 2 j=1,2
      if(k.le.0) go to 3
      do 4 m=1,k
      if(midx(m,nph)-mndx(i,j))4,2,5
 4    continue
 3    k=k+1
      midx(k,nph)=mndx(i,j)
      pux(k,nph)=px(i,j)
      go to 2
 5    do 6 l=k,m,-1
      midx(l+1,nph)=midx(l,nph)
 6    pux(l+1,nph)=pux(l,nph)
      k=k+1
      midx(m,nph)=mndx(i,j)
      pux(m,nph)=px(i,j)
 2    continue
 1    km(nph)=k
      return
      end
sphdist.f
*-----        Location Routines      -------------------------------- + ----
c                                                                     ydist
      subroutine ydist
     ^        (dlats, dlons, dlatr, dlonr, delta, cazim, bazim, azima)
c
c     AUTHOR:  Brian L.N. Kennett  RSES, ANU
c     DATE:    January 1985
c     PURPOSE:
c             YDIST        Calculates distance and azimuth
c                          for spheroidal earth between
c                          specified geographic source and
c                          receiver station coordinates
c
*----------------------------------------------------------------------*------*
c     PARAMETERS
c
      real    dlats, dlons, dlatr, dlonr, delta, cazim, bazim
c
c     dlats  latitude of source
c     dlons  longitude of source
c     dlatr  latitude of receiver
c     dlonr  longitude of receiver
c     delta  angular distance
c     cazim  apparent azimuth at an array
c     bazim   azimuth from epicentre to receiver
c
*----------------------------------------------------------------------*------*
c
c     implicit real*8 (a-h,o-z)
      real*8 ecc,re,ec1,pi,pib2,degr,rlats,rlons,rlatr
      real*8 rlonr,glats,glatr,sps,cps,spr,cpr,rs,rr
      real*8 trs,prs,trr,prr,AS,BS,CS,DS,ES,GS,HS,KS
      real*8 AR,BR,CR,DR,ER,GR,HR,KR
      real*8 cosdr,deltar,sindr,deltak,szs,czs,szr,czr
      real*8 e,x,y
      real azima
c                          radius on spheroid
      gra(x,y,e) = dsqrt( (1.0d0-e)**2 /
     &                   ((1.0d0-e*y)**2 + e*e*x*y ) )
      ecc = 0.003367
      re = 6378.388
      ec1 = (1.0d0-ecc)**2
      pi = 3.141592653589793
      pib2 = pi/2.0
      degr = pi/180.0
      rlats = dlats*degr
      rlons = dlons*degr
      rlatr = dlatr*degr
      rlonr = dlonr*degr
c                          geocentric coordinates
      glats = datan2 ( ec1*dsin(rlats) ,dcos(rlats) )
      glatr = datan2 ( ec1*dsin(rlatr) ,dcos(rlatr) )
      sps = dsin(glats)**2
      cps = dcos(glats)**2
      spr = dsin(glatr)**2
      cpr = dcos(glatr)**2
c                          radii at source,receiver
      rs = re*gra(sps,cps,ecc)
      rr = re*gra(spr,cpr,ecc)
c
      trs = pib2 - glats
      prs = dlons*degr
      trr = pib2 - glatr
      prr = dlonr*degr
c                          direction cosines for source
      AS = dsin(trs)*dcos(prs)
      BS = dsin(trs)*dsin(prs)
      CS = dcos(trs)
      DS = dsin(prs)
      ES = -dcos(prs)
      GS = dcos(trs)*dcos(prs)
      HS = dcos(trs)*dsin(prs)
      KS = -dsin(trs)
c                          direction cosines for receiver
      AR = dsin(trr)*dcos(prr)
      BR = dsin(trr)*dsin(prr)
      CR = dcos(trr)
      DR = dsin(prr)
      ER = -dcos(prr)
      GR = dcos(trr)*dcos(prr)
      HR = dcos(trr)*dsin(prr)
      KR = -dsin(trr)
c                          distance
      cosdr = AS*AR + BS*BR + CS*CR
      deltar = dacos(cosdr)
      sindr = dsin(deltar)
c
      deltak = deltar*0.5d0*(rr+rs)
      delta = deltar/degr
c                          azimuth
      szs = DS*AR + ES*BR
      czs = GS*AR + HS*BR + KS*CR
      szr = DR*AS + ER*BS
      czr = GR*AS + HR*BS + KR*CS
c                          azima - azimuth to source
c                          bazim - backazimuth from source
c                          cazim - apparent azimuth at an array
      if (szr.eq.0.0) then
        bazim = 0.0
        if(dlats.gt.dlatr)then
           azima = 360.0
        else
           azima = 180.0
        endif
      else
        bazim = datan2(-szs ,-czs ) /degr
        azima = datan2(-szr ,-czr ) /degr
      end if
      if( bazim .lt. 0.0) bazim = bazim + 360.0
      cazim = azima + 180.0
      if( azima .lt. 0.0) azima = azima + 360.0
c
      if( cazim.lt. 0.0) cazim = cazim + 360.0
c
      return
      end
stack.f90
module stack

  type Tnode
     integer                              :: value
     type(Tnode),pointer                  :: next
  end type Tnode

  type(Tnode),pointer                     :: head,t  
   
  type Tnode2
     integer                              :: value1,value2
     type(Tnode2),pointer                 :: next
  end type Tnode2

  type(Tnode2),pointer                    :: head2,t2  

end module stack

  subroutine stackinit
    use stack
    allocate(head)
    nullify(head%next) ; head%value = 0
  end subroutine stackinit

  subroutine push(p)
    use stack
    integer :: p
    allocate(t)
    t%next => head%next ; t%value = p
    head%next => t
    nullify(t)
  end subroutine push

  subroutine pop(x)
    use stack
    integer :: x
    t => head%next ; head%next => t%next
    x = t%value
    deallocate(t)
  end subroutine pop

  subroutine stackempty(i)
    use stack
    integer  :: i
    i=0
    if (.not.associated(head%next)) i=1
  end subroutine stackempty

  subroutine stackflush
    use stack
    deallocate(head)
  end subroutine stackflush


  subroutine stackpairinit
    use stack
    allocate(head2)
    nullify(head2%next) ; head2%value1 = 0; head2%value2 = 0
  end subroutine stackpairinit

  subroutine pushpair(p1,p2)
    use stack
    integer :: p1,p2
    allocate(t2)
    t2%next => head2%next ; t2%value1 = p1 ; t2%value2 = p2
    head2%next => t2
    nullify(t2)
  end subroutine pushpair

  subroutine poppair(x1,x2)
    use stack
    integer :: x1,x2
    t2 => head2%next ; head2%next => t2%next
    x1 = t2%value1 ; x2 = t2%value2 
    deallocate(t2)
  end subroutine poppair

  subroutine stackpairempty(i)
    use stack
    integer  :: i
    i=0
    if (.not.associated(head2%next)) i=1
  end subroutine stackpairempty

  subroutine stackpairflush
    use stack
    deallocate(head2)
  end subroutine stackpairflush

svdlib.f90
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
        implicit none
        INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)
      integer, PARAMETER  :: NMAX=500
      INTEGER  :: m,mp,n,np
      real(kind=dp) :: a(mp,np),v(np,np),w(np)
      INTEGER  :: i,its,j,jj,k,l,nm
      real(kind=dp) :: anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag

      g=0.0_dp
      scale=0.0_dp
      anorm=0.0_dp
      do i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0_dp
        s=0.0_dp
        scale=0.0_dp
        if(i.le.m)then
          do k=i,m
            scale=scale+abs(a(k,i))
          end do
          if(scale.ne.0.0_dp)then
            do k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
            end do
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do j=l,n
              s=0.0_dp
              do k=i,m
                s=s+a(k,i)*a(k,j)
              end do
              f=s/h
              do k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              end do
            end do
            do k=i,m
              a(k,i)=scale*a(k,i)
            end do
          endif
        endif
        w(i)=scale *g
        g=0.0_dp
        s=0.0_dp
        scale=0.0_dp
        if((i.le.m).and.(i.ne.n))then
          do k=l,n
            scale=scale+abs(a(i,k))
          end do
          if(scale.ne.0.0)then
            do k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
            end do
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do k=l,n
              rv1(k)=a(i,k)/h
            end do
            do j=l,m
              s=0.0
              do k=l,n
                s=s+a(j,k)*a(i,k)
              end do
              do k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
              end do
            end do
            do k=l,n
              a(i,k)=scale*a(i,k)
            end do
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      end do
      do i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0_dp)then
            do j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
            end do
            do j=l,n
              s=0.0_dp
              do k=l,n
                s=s+a(i,k)*v(k,j)
              end do
              do k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
              end do
            end do
          endif
          do j=l,n
            v(i,j)=0.0_dp
            v(j,i)=0.0_dp
          end do
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
      end do
      do  i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do j=l,n
          a(i,j)=0.0_dp
        end do
        if(g.ne.0.0_dp)then
          g=1.0_dp/g
          do j=l,n
            s=0.0_dp
            do k=l,m
              s=s+a(k,i)*a(k,j)
            end do
            f=(s/a(i,i))*g
            do k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            end do
          end do
          do j=i,m
            a(j,i)=a(j,i)*g
          end do
        else
          do j= i,m
            a(j,i)=0.0_dp
          end do
        endif
        a(i,i)=a(i,i)+1.0_dp
      end do
      do k=n,1,-1
        do its=1,30
          do l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
          end do                        
1         c=0.0_dp
          s=1.0_dp
          do i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1._dp/h
            c= (g*h)
            s=-(f*h)
            do j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
            end do
          end do
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0_dp)then
              w(k)=-z
              do j=1,n
                v(j,k)=-v(j,k)
              end do
            endif
            goto 3
          endif
          if(its.eq.50) stop 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0_dp)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0_dp
          s=1.0_dp
          do j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
            end do
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0_dp)then
              z=1.0_dp/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
            end do
          end do
          rv1(l)=0.0_dp
          rv1(k)=f
          w(k)=x
        end do

3       continue

      end do
      return
      end subroutine 

      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
        INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)
      integer, PARAMETER :: NMAX=500
      INTEGER   :: m,mp,n,np
      real(kind=dp) :: b(mp),u(mp,np),v(np,np),w(np),x(np)
      INTEGER   :: i,j,jj
      real(kind=dp) :: s,tmp(NMAX)
      do j=1,n
        s=0.
        if(w(j).ne.0.)then
          do i=1,m
            s=s+u(i,j)*b(i)
          end do
          s=s/w(j)
        endif
        tmp(j)=s
      end do
      do j=1,n
        s=0.
        do jj=1,n
          s=s+v(j,jj)*tmp(jj)
        end do
        x(j)=s
      end do
      return
      END

      FUNCTION pythag(a,b)
        INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)
      real(kind=dp) :: pythag
      real(kind=dp) :: a,b
      real(kind=dp) :: absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1._dp+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1._dp+(absa/absb)**2)
        endif
      endif
      return
      END

teleseismic.f90
!*****************************************************************************************************
!***********************************************************************************************
! This subroutine takes a teleseismic source as input,  initializes the bottom
! intersection of the grid and the sides facing the source, and does a sweep through all regions 
! in the upward direction to the surface.

subroutine initialize_teleseismic_source(s)
use mod_3dfm
implicit none

type(Tsource)       :: s  ! the source with its properties on the main grid

type(Tintersection),pointer       :: itopc,ibotc ! top and bottom intersections of main grid 
type(Tregion),pointer             :: reg          ! the source region in the main grid
logical                           :: do_northside,do_southside,do_westside,do_eastside
integer,dimension(:),allocatable  :: init_type

! local stuff
integer :: n,m,i,j,i1,i2,i3,vtype,n_init,n_tot_ref
real(kind=4)  :: latsw,latnw,latne,latse
real(kind=4)  :: lonsw,lonnw,lonne,lonse
real(kind=4)  :: slat,slon,deltas,cazim,bazim,azima
real(kind=8)  :: deg_to_rad 

deg_to_rad = acos(-1.0_dp)/180.0_dp

print *,'initializing teleseismic source',s%id,s%teleseismic_id

! allocate the arrays that will contain the indices of the source time fields in the time field list
allocate(s%first_tf_up(n_vtypes),s%first_tf_down(n_vtypes))


! determine which sides of the grid face the source

do_northside = .false. ; do_southside = .false.
do_westside = .false.  ; do_eastside = .false.

! corners of the gird on the surface and source coordinates in real*4

latsw=(pgrid%lat0)/deg_to_rad
latnw=(pgrid%lat0 + (pgrid%nlat-1)*pgrid%dlat0)/deg_to_rad
latne=(pgrid%lat0 + (pgrid%nlat-1)*pgrid%dlat0)/deg_to_rad
latse=(pgrid%lat0)/deg_to_rad

lonsw=(pgrid%long0)/deg_to_rad
lonnw=(pgrid%long0)/deg_to_rad
lonne=(pgrid%long0 + (pgrid%nlong-1)*pgrid%dlong0)/deg_to_rad
lonse=(pgrid%long0 + (pgrid%nlong-1)*pgrid%dlong0)/deg_to_rad

slat=(s%lat)/deg_to_rad
slon=(s%long)/deg_to_rad

call ydist(slat,slon,latsw,lonsw,deltas,cazim,bazim,azima)
do_westside = (sin(cazim*deg_to_rad)> 0.01)  
call ydist(slat,slon,latnw,lonnw,deltas,cazim,bazim,azima)
do_westside = (sin(cazim*deg_to_rad)> 0.01).and.do_westside

call ydist(slat,slon,latnw,lonnw,deltas,cazim,bazim,azima)
do_northside = (cos(cazim*deg_to_rad) < -0.01)  
call ydist(slat,slon,latne,lonne,deltas,cazim,bazim,azima)
do_northside = (cos(cazim*deg_to_rad) < -0.01).and.do_northside

call ydist(slat,slon,latne,lonne,deltas,cazim,bazim,azima)
do_eastside = (sin(cazim*deg_to_rad) < -0.01)  
call ydist(slat,slon,latse,lonse,deltas,cazim,bazim,azima)
do_eastside = (sin(cazim*deg_to_rad) < -0.01).and.do_eastside

call ydist(slat,slon,latse,lonse,deltas,cazim,bazim,azima)
do_southside = (cos(cazim*deg_to_rad) > 0.01)  
call ydist(slat,slon,latsw,lonsw,deltas,cazim,bazim,azima)
do_southside = (cos(cazim*deg_to_rad) > 0.01).and.do_southside

 print *,'exposed sides: N',do_northside,' S ',do_southside,' W ', &
      do_westside,' E ',do_eastside


! make a list of nodes to be initialized, region by region

do m=n_regions,1,-1

   reg => region(m)
   allocate(init_type(reg%nnode))
   init_type = 0
   n_init = 0

!   print *,'region',reg%id,reg%nnode

   do n=1,reg%nnode

      if (reg%node(n)%i1==0 .and. reg%node(n)%i2==n_intersections) then
         n_init=n_init+1
         init_type(n)=1
         cycle
      endif
         
      if (do_northside .and. reg%r(n)*abs(reg%lat(n)-pgrid%latmax)<pgrid%tolerance) then
         n_init=n_init+1
         init_type(n)=2
         cycle
      endif

      if (do_southside .and. reg%r(n)*abs(reg%lat(n)-pgrid%lat0)<pgrid%tolerance) then
         n_init=n_init+1
         init_type(n)=3
         cycle
      endif

      if (do_eastside .and. reg%r(n)*abs(reg%long(n)-pgrid%longmax)<pgrid%tolerance) then
         n_init=n_init+1
         init_type(n)=4
         cycle
      endif

      if (do_westside .and. reg%r(n)*abs(reg%long(n)-pgrid%long0)<pgrid%tolerance) then
         n_init=n_init+1
         init_type(n)=5
         cycle
      endif

   end do

   reg%n_init = n_init

!   print *,n_init,'nodes in list'
!   print *,count(init_type == 1),count(init_type == 2),count(init_type == 3),count(init_type == 4),count(init_type == 5)


   if (n_init > 0) then
      allocate(reg%init_id(reg%n_init),reg%init_type(reg%n_init),reg%init_arrivaltime(reg%n_init), &
           reg%init_time_gradient(3,reg%n_init))
      i=0
      do n=1,reg%nnode
         if (init_type(n) /= 0) then
            i=i+1
            reg%init_id(i)=n
            reg%init_type(i)=init_type(n)
         endif
      end do
   endif

   deallocate(init_type)

end do  ! loop making init node list for each region

print *,'init nodes identified'

! for all regions, get the starting times and time gradients on the bottom interface
! and relevant sides from the ttimes program

do m=1,n_regions
   if (region(m)%n_init>0) call teleseismic_initialization(region(m),s)
end do

print *

! all boundary nodes hit by the incoming front are now identified and have arrival times
! incoming time gradients  if a solution for that node exists

!----------------------------------------------------------------------------------------

! a teleseismic source is intialized by doing a path sequence going from bottom to top
! for 1 or 2 velocity types as required

s%n_tf_init = n_regions

! sweep through all regions towards the surface

vloop: do vtype=1,n_vtypes

   do n=1,s%n_tf_init

      reg => region(n_regions - n + 1)

      itopc => intersection(reg%id)      ! the top intersection of the region
      ibotc => intersection(reg%id+1)    ! the bottom intersection of the region
      allocate(reg%arrivaltime(reg%nnode),reg%time_gradient(3,reg%nnode),reg%node_status(reg%nnode))
      reg%node_status=-1
      reg%arrivaltime = huge_time
      reg%time_gradient=0.0_dp

      if (reg%id < n_regions) then  ! if the region is not the bottom region

   ! create a narrow band on the starting intersection from previous result

         call refract_gradient(ibotc,region(reg%id+1),vtype,1) 

         i=reg%nnode-ibotc%nnode+1
         j=reg%nnode
         reg%arrivaltime(i:j)=ibotc%arrivaltime
         reg%time_gradient(1:3,i:j)=ibotc%time_gradient(1:3,1:ibotc%nnode)
         reg%node_status(i:j) = 1


      endif


   ! create a narrow band from the incoming teleseismic wave
   ! ( copy arrival times from init, refract gradients, set node status)

      if (reg%n_init>0) call refract_teleseismic_front(reg,vtype,n_tot_ref)

      if (n==1 .and. n_tot_ref == reg%n_init) then

         print *,'teleseismic front incompatible with start from  bottom region'
         print *,'all boundary nodes experienced total reflection or inside hits'
         print *,'for vtype =',vtype,' No paths starting with this type are possible'
         print *
         cycle vloop

      endif


! do the fast marching sweep across the main grid region containing the source

      call propagate(reg,vtype)

      print *,'propagation through region',reg%id,' finished'

! transfer regional travel times to interfaces 

      if (itopc%nnode > 0 .and. .not. associated(itopc%arrivaltime)) then
         allocate(itopc%arrivaltime(itopc%nnode))
         allocate(itopc%time_gradient(3,itopc%nnode))
      endif
      if (ibotc%nnode > 0 .and. .not. associated(ibotc%arrivaltime)) then
         allocate(ibotc%arrivaltime(ibotc%nnode))
         allocate(ibotc%time_gradient(3,ibotc%nnode))
      endif

      do i=1,reg%nnode

         i1 = reg%node(i)%i1 ; i2 = reg%node(i)%i2 ; i3 = reg%node(i)%i3

         if (i1 == 0) then

            intersection(i2)%arrivaltime(i3) = reg%arrivaltime(i)
            intersection(i2)%time_gradient(1:3,i3) = reg%time_gradient(1:3,i)

         endif

      end do
      
!   print *,'min time on top isec',minval(itopc%arrivaltime)

  ! transfer time field to the array of saved time fields

      s%n_time_fields=s%n_time_fields+1
      allocate(s%time_field(s%n_time_fields)%arrivaltime(reg%nnode))
      s%time_field(s%n_time_fields)%arrivaltime=reg%arrivaltime
      allocate(s%time_field(s%n_time_fields)%time_gradient(3,reg%nnode))
      s%time_field(s%n_time_fields)%time_gradient=reg%time_gradient

      print *,'results written to timefield',s%n_time_fields

  
! save pointer to region,start and non-start interfaces and tree structure pointers

      s%time_field(s%n_time_fields)%reg =>  reg 

      s%time_field(s%n_time_fields)%istart => ibotc
      s%time_field(s%n_time_fields)%inonstart =>  itopc
      s%time_field(s%n_time_fields)%vtype = vtype

      if (n==1) s%first_tf_up(vtype)=s%n_time_fields

      if (n>1) then
         s%time_field(s%n_time_fields)%prev_tf = n-1 + (vtype-1)*s%n_tf_init
         s%time_field(s%n_time_fields-1)%next_tf(1+(vtype-1)*4)= s%n_time_fields
         print *,'field is ctype',1+(vtype-1)*4,'of timefield',s%n_time_fields-1
      endif

! deallocate  everything that is no longer required

!      print *,'timefield attributes saved'

      deallocate(reg%arrivaltime,reg%time_gradient,reg%node_status)

   end do

end do vloop ! vtypes


print *,'initial timefields for teleseismic source established'


! stop ' temp stop in init teleseismic source'

return

end subroutine initialize_teleseismic_source




!*****************************************************************************************************
!
subroutine  teleseismic_initialization(reg,s)

  use mod_3dfm
  implicit none
!
!     Determines ak135 traveltimes for a specified phase to
!     the boundary nodes of a region . It first solves iteratively for the   
!     angular separation at which the teleseismic ray and the ray from the
!     boundary node hit the surface with the same angle (ray parameter)
!     The time at the node is then the difference between the travel times 
!     of the two phases to this surface location 
!
!     Marthijn de Kool
!     Australian National University
!
!    
!     NOTE: Makes use of subroutines from the freeware
!           program ttimes, and uses some code from Nick Rawlinson's aktsurf.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  type(Tregion) :: reg
  type(Tsource)       :: s

  integer             :: i,j,k,n,m,phid,iter,n_suspect_init,vtype
  integer,parameter   :: maxnp=60
  real                :: tt(maxnp),dtdd(maxnp),dtdh(maxnp),dddp(maxnp)
  real                :: evlat,evlon,evdep,ndlat,ndlon,rslat,usrc(2),nddep
  real                :: deltas1,deltas2,deltas3,delta
  real                :: deltan1,deltan2,deltan3,gt_norm,f1,f2,f3,delta_min,f1_old,f2_old
  real                :: stt1,stt2,stt3,sdtdd1,sdtdh1,sdtdd2,sdtdh2,sdtdd3,sdtdh3
  real                :: ntt1,ntt2,ntt3,ndtdd1,ndtdh1,ndtdd2,ndtdh2,ndtdd3,ndtdh3
  real                :: dt1,dt2,dt3
  real                :: deltas,cazim,bazim,edist,bazr,etcor,azima
  real                :: delta_pos,delta_neg,time_error,accuracy
  real(kind=dp)       :: deg_to_rad,deg_per_km,surf_vel,det

  real,parameter      :: pi=3.1415926535
  character(len=8)    :: phcd(maxnp),phlst(maxnp),tele_phase,loc_phase,s_id,n_id
  character(len=25)   :: modnam
  logical             :: prnt(2),s1_invalid,s2_invalid

!c     tt = traveltime of phase
!c     dtdd, dtdh,dddp = partial derivatives
!c     phcd = Phase id
!c     phlst = Phase list
!c     prnt = debugging print flag
!c     modnam = Name of velocity model
!c     evlat = event latitude
!c     evlon = event longitude
!c     evdep = event depth
!c     ndlat = node latitude
!c     ndlon = node longitude
!c     rslat = station co-latitude
!c     deltas = event-station angular distance
!c     cazim,bazim,azima = azimuth information
!c     m = number of phases found
!c     edist = event-station distance
!c     bazr = adjusted bazim
!c     etcor = elliptical correction for travel time
!c     phid = Phase id of specific phase
!c     tph = traveltime of specific phase
!c     maxnp = maximum number of ak135 phases
!c     tele_phase = specific phase name for given event
!c     pi = pi
!c     er = Earth radius
!c     kmpd = Number of km per great circle degree

 
! some validation of the input

  if (.not.s%is_teleseismic) stop 'teleseismic_initialization called with local source'

  if ( count(reg%r<(earth_radius-800.0_dp)) > 0 )  then
     print *, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print *, 'sorry, teleseismic initialization not possible if nodes on the'
     print *,' bottom interface have a depth > 800 km'
     stop 'bottom interface node(s) with depth > 800 km in tele_init'
  endif

  deg_to_rad = acos(-1.0_dp)/180.0_dp
  deg_per_km = 1.d0/(earth_radius*deg_to_rad)
  gt_norm = 0.20 / deg_per_km    ! establish order of magnitude of d(Tarrival)/d(Delta) at the surface
                                   ! assuming velocity is 5 km/sec

!      Specify model type

  modnam='ak135'

!     Open travel time tables

  prnt(1)=.false.
  prnt(2)=.false.
  phlst(2)='  '
  call tabin(10,modnam)

!     get source location and phase type

  evlat = s%lat/deg_to_rad
  evlon = s%long/deg_to_rad
  evdep = earth_radius - s%r
  tele_phase = s%teleseismic_phase


! determine whether the last leg of the teleseismic phase is P or S

  i = index(tele_phase,'P',back=.true.)
  j = index(tele_phase,'S',back=.true.)
  if (i>j) then
     loc_phase = 'P'
     vtype=1
     surf_vel=5.8_dp
  else
     loc_phase = 'S'
     vtype=n_vtypes
     surf_vel=3.46_dp
  endif


!  Now loop through all the boundary nodes

  n_suspect_init = 0
  print *,'starting initialization of',reg%n_init,' boundary nodes of region',reg%id

  nodeloop: do m=1,reg%n_init

     n=reg%init_id(m)

 ! convert position of the node to input taken by ttimes
     
     ndlat = reg%lat(n)/deg_to_rad
     ndlon = reg%long(n)/deg_to_rad
     nddep = earth_radius - reg%r(n)
     if (nddep < 0.0) then
        reg%init_arrivaltime(m)=huge_time
        reg%init_time_gradient(1:3,m) = 0.d0
        cycle nodeloop
     endif


!  This is a bit of a hack fix for the case when event and grid node longitude are equal. 
!  For some reason the ellipticity corrections are wrong for this case

     if(ndlon.eq.evlon)then
        ndlon=ndlon+pgrid%dlong0/(50.0*deg_to_rad)
     endif


! get angular distance between source and node
! note that this is the SPHEROIDAL distance, not the SPHERICAL
 
     call ydist(evlat,evlon,ndlat,ndlon,deltas,cazim,bazim,azima)

! if the node lies on the surface, it is a special case. No iteration necessary

     if (abs(nddep) < pgrid%tolerance) then

        call phase_time(deltas,evdep,tele_phase,stt1,sdtdd1,sdtdh1)

        if (phid == 0) then  ! if phase does not exist, try next node
           reg%init_arrivaltime(m)=huge_time
           reg%init_time_gradient(1:3,m) = 0.d0
           cycle nodeloop
        endif

        phlst(1)=tele_phase
        call brnset(1,phlst,prnt)
        call depset(evdep,usrc)
        rslat = (90.-evlat)*0.017453292
        call ellref(rslat)
        edist=deltas*0.017453292
        bazr=bazim*0.017453292
        call ellcor(edist, bazr, evdep, tele_phase, etcor)
        stt1=stt1+etcor

        reg%init_arrivaltime(m) = dble(stt1)

        reg%init_time_gradient(2,m) = sdtdd1*deg_per_km*cos(cazim*deg_to_rad)
        reg%init_time_gradient(3,m) = sdtdd1*deg_per_km*sin(cazim*deg_to_rad)
        det = 1.0_dp/surf_vel**2 - (reg%init_time_gradient(3,m)**2 + reg%init_time_gradient(2,m)**2)
        if (det >= 0.0_dp) then
           reg%init_time_gradient(1,m) = sqrt(det)
        else
           stop 'det < 0 at surface node during teleseismic initialization'
        endif

        cycle nodeloop

     endif


 ! for all nodes not on the surface:
 ! we will use the secant method to find the angular distance of the point on the surface
 ! where the ray parameters (dtdd) of the ray from the boundary node and
 ! the teleseismic source are equal.

 ! if we can find two points bracketing the solution (function positive and negative)
 ! a solution is sure to exist, if these can not be found we assume the phase does not exist at the node

! the travel time and arrival time gradient from the source to the two initial guess points

     delta  = nddep*deg_per_km  ! the first guess at the size of the interval is the depth of the node
     deltas1 = deltas            ! point 1 lies above the node
     deltas2 = deltas + delta    
     s1_invalid = .false.
     s2_invalid = .false.

     call phase_time(deltas1,evdep,tele_phase,stt1,sdtdd1,sdtdh1)
     if (phid == 0) s1_invalid = .true.
     call phase_time(deltas2,evdep,tele_phase,stt2,sdtdd2,sdtdh2)
     if (phid == 0) s2_invalid = .true.

! if both points are invalid, a valid arrival time can not be found

     if (s1_invalid .and. s2_invalid) then
        reg%init_arrivaltime(m)=huge_time
        reg%init_time_gradient(1:3,m) = 0.d0
 !       print *,'failed because 2 invalid times'
        cycle nodeloop
     endif


 ! the travel time and arrival time gradient from the node to the two initial guess points
 ! these should always exist, so if they don't its a serious error stop
 ! note that we use a special approximation for the ray parameter if the node lies
 ! very close to the surface since for this case ttimes returns a wrong value

     deltan1 = 0.0
     deltan2=delta

     call phase_time(deltan1,nddep,loc_phase,ntt1,ndtdd1,ndtdh1)
     if (nddep <= 5.0 ) ndtdd1 = deltan1/(sqrt(deltan1**2+(nddep*deg_per_km)**2)*surf_vel*deg_per_km)
     if (phid == 0) stop 'tele_init: phase does not exist n1'

     call phase_time(deltan2,nddep,loc_phase,ntt2,ndtdd2,ndtdh2)
     if (nddep <= 5.0 ) ndtdd2 = deltan2/(sqrt(deltan2**2+(nddep*deg_per_km)**2)*surf_vel*deg_per_km)
     if (phid == 0) stop 'tele_init: phase does not exist n2'


! if only one of the two paths from the source is invalid, try if we can find better initial guesses

     if (s1_invalid .or. s2_invalid) then

    ! if more distant point p2 is valid, but not the one above the node (p1)

        if (s1_invalid) then

      ! try bringing point 1 closer to point 2
           
           deltas1 = deltas1+delta/2.0
           iter =0
           do while (iter < 6 .and. s1_invalid)

              iter = iter + 1

              s1_invalid = .false.
              call phase_time(deltas1,evdep,tele_phase,stt1,sdtdd1,sdtdh1)
              if (phid == 0) s1_invalid = .true.

         ! if the new point is valid
              if (.not.s1_invalid) then

                 f1 = sdtdd1 - ndtdd1

                 if (f1 >= 0.0) then

              ! bracketing points found, finished
                       
                    cycle

                 else
                
              ! overshot, move closer to p1

                    deltas1=deltas1-delta/2**(iter+1)
                    s1_invalid = .true.

                 endif

              else
  
              ! still not valid, move closer to p2

                 deltas1 = deltas1+delta/2**(iter+1)

              endif

           end do

        ! if iteration didn't help, write off node as invalid
           if (s1_invalid) then
              reg%init_arrivaltime(m)=huge_time
              reg%init_time_gradient(1:3,m) = 0.d0
 !             print *,'failed because node p1 invalid and search failed'
              cycle nodeloop
           endif

        endif



    ! if more distant point p2 is invalid, but not the one above the node (p1)

        if (s2_invalid) then

      ! try bringing point 2 closer to point 1
           
           deltas2 = deltas2-delta/2
           iter =0
           do while (iter < 6 .and. s2_invalid)

              iter = iter + 1

              s2_invalid = .false.
              call phase_time(deltas2,evdep,tele_phase,stt2,sdtdd2,sdtdh2)
              if (phid == 0) s2_invalid = .true.

         ! if the new point is valid
              if (.not.s2_invalid) then

                 f2 = sdtdd2 - ndtdd2

                 if (f2 <= 0.0) then

              ! bracketing points found, finished
                       
                    cycle

                 else
                
              ! overshot, move closer to p2 again

                    deltas2=deltas2+delta/2**(iter+1)
                    s2_invalid = .true.

                 endif

              else
  
              ! still not valid, move closer to p1

                 deltas2 = deltas2-delta/2**(iter+1)

              endif

           end do

        ! if iteration didn't help, write off node as invalid
           if (s2_invalid) then
              reg%init_arrivaltime(m)=huge_time
              reg%init_time_gradient(1:3,m) = 0.d0
 !             print *,'failed because node p2 invalid and search failed'
              cycle nodeloop
           endif

        endif

     endif



  ! we have two valid points

     f1 = sdtdd1 - ndtdd1
     f2 = sdtdd2 - ndtdd2

     if (f1 < 0.0) then
        
        print *,'Warning 1 in tele_init'
        print '(a5,7f12.5)','p1',deltas1,stt1,deltan1,ntt1,f1,sdtdd1,ndtdd1
        print '(a5,7f12.5)','p2',deltas2,stt2,deltan2,ntt2,f2,sdtdd2,ndtdd2
        reg%init_arrivaltime(m)=huge_time
        reg%init_time_gradient(1:3,m) = 0.d0
        cycle nodeloop
     endif

!     print *,'first two fs',f1,f2,deltan2


 ! if the root is not bracketed in the first interval, try extending it

     iter = 1
     delta_pos=deltan1
     delta_neg=deltan2

     do while (f1*f2 > 0.0)

        iter = iter + 1

        if (iter > 3) then
           reg%init_arrivaltime(m)=huge_time
           reg%init_time_gradient(1:3,m) = 0.d0
!           print *,'failed because root could not be bracketed'
           cycle nodeloop
        endif

        deltas2 = deltas2 + delta  
        deltan2 = deltan2 + delta

        call phase_time(deltas2,evdep,tele_phase,stt2,sdtdd2,sdtdh2)
        if (phid == 0) stop 'phase does not exist s2_rep'
        call phase_time(deltan2,nddep,loc_phase,ntt2,ndtdd2,ndtdh2)
        if (nddep <= 5.0 ) ndtdd2 = deltan2/(sqrt(deltan2**2+(nddep*deg_per_km)**2)*surf_vel*deg_per_km)
!        if (nddep <= 5.0 ) ndtdd2 = sin(deltan2*deg_to_rad)/(surf_vel*deg_per_km)
        if (phid == 0) stop 'phase does not exist n2_rep'

        f2 = sdtdd2 - ndtdd2

!        print *,'trying to bracket root,next f2',f2,deltan2

     end do


!  start the root finding iteration 
     deltan1=deltas1-deltas
     deltan2=deltas2-deltas
     iter = 0
     delta_min=pgrid%dlong0/(50.0*deg_to_rad)
     accuracy = pgrid%dr0/(50.0*reg%velocity(1,1))
     f3 = gt_norm ; f1_old=2.*f3 ; f2_old=2.*f3
     delta_pos=deltan1
     delta_neg=deltan2


itloop: do while (abs(f3) > 1.e-7*gt_norm  .and. abs(deltan2-deltan1) > 1.e-6*deltas) 

        iter = iter + 1

        deltan3 = deltan1 +((f1/(f1-f2))*(deltan2-deltan1))

      ! if solution deteriorates rather than improves, switch to bisection

        if ((deltan3-delta_pos)*(deltan3-delta_neg)>0.0 .or. &
             (abs(f3)>=abs(f1_old).and.abs(f3)>=abs(f2_old))) then
            deltan3=(delta_pos+delta_neg)/2.
        endif


        deltas3 = deltas + deltan3  

        call phase_time(deltas3,evdep,tele_phase,stt3,sdtdd3,sdtdh3)
        if (phid == 0) stop 'phase does not exist s3'
        s_id=phcd(phid)
        call phase_time(deltan3,nddep,loc_phase,ntt3,ndtdd3,ndtdh3)
        if (nddep <= 5.0 ) ndtdd3 = deltan3/(sqrt(deltan3**2+(nddep*deg_per_km)**2)*surf_vel*deg_per_km)
        if (phid == 0) stop 'phase does not exist n3'
        n_id=phcd(phid)

        f3 = sdtdd3 - ndtdd3

        if (abs(f3) <= 1.e-7*gt_norm ) exit itloop  ! exit if converged

        if (iter > 20) then

           dt1=stt1-ntt1 ; dt2=stt2-ntt2 ; dt3=stt3-ntt3 ;
           time_error=max(abs(dt1-dt2),abs(dt1-dt3),abs(dt2-dt3))

        ! if the solution has not fully converged but the maximum time error is acceptable
        ! still allow the solution to be used, otherwise discard and treat node as not initialized

           if (time_error < accuracy) then
              n_suspect_init = n_suspect_init+1
              exit itloop
           else
              reg%init_arrivaltime(m)=huge_time
              reg%init_time_gradient(1:3,m) = 0.d0
              cycle nodeloop
           endif

        endif

 
        f1_old=f1 ; f2_old=f2

        if (abs(f1)>abs(f2)) then
           f1 = f3 ; deltan1=deltan3 ; deltas1=deltas3 ; ntt1 = ntt3 ; stt1 = stt3 ; sdtdd1=sdtdd3 ; ndtdd1=ndtdd3
        else
           f2 = f3 ; deltan2=deltan3 ; deltas2=deltas3 ; ntt2 = ntt3 ; stt2 = stt3 ; sdtdd2=sdtdd3 ; ndtdd2=ndtdd3
        endif

        if (f1 > 0.0) then
           delta_pos = max(delta_pos,deltan1)
        else
           delta_neg = min(delta_neg,deltan1)
        endif
        if (f2 > 0.0) then
           delta_pos = max(delta_pos,deltan2)
        else
           delta_neg = min(delta_neg,deltan2)
        endif


        delta =abs(deltan2-deltan1)


     end do itloop


!    Apply elliptical corrections
!

     bazr=bazim*0.017453292

     phlst(1)=tele_phase
     call brnset(1,phlst,prnt)
     call depset(evdep,usrc)
     rslat = (90.-evlat)*0.017453292
     call ellref(rslat)
     edist=deltas3*0.017453292
     call ellcor(edist, bazr, evdep, tele_phase, etcor)
     stt3=stt3+etcor

     phlst(1)=loc_phase
     call brnset(1,phlst,prnt)
     call depset(nddep,usrc)
     rslat = (90.-ndlat)*0.017453292
     call ellref(rslat)
     edist=deltan3*0.017453292
     call ellcor(edist, bazr, nddep, loc_phase, etcor)
     ntt3=ntt3+etcor


! arrival time at the bottom interface node is the difference in arrival times 

     reg%init_arrivaltime(m) = stt3 - ntt3


 ! store the incoming (before refraction at bottom interface) time gradient in the ak135 model

     reg%init_time_gradient(1,m) = ndtdh3
     reg%init_time_gradient(2,m) = ndtdd3*deg_per_km*cos(cazim*deg_to_rad)* &
          earth_radius/reg%r(n)
     reg%init_time_gradient(3,m) = ndtdd3*deg_per_km*sin(cazim*deg_to_rad)* &
          earth_radius/reg%r(n)

!     if(reg%init_time_gradient(2,m)>0.0) then
!        print *,m,cazim,azima
!        stop 'gothcha'
!     endif

  enddo nodeloop


! print some diagnostics on how the initialization went
  
  if (n_suspect_init > 0) print *,'warning: teleseismic initialization problematic at',n_suspect_init,' nodes'
  i=count(reg%init_arrivaltime == huge_time)
  if (i>0) print *,'warning: teleseismic initialization failed at     ',i,' nodes'
  if (i==0 .and. n_suspect_init==0) print *,'teleseismic initialization successful at all nodes'

  return

  contains

    subroutine phase_time(del,dep,phase,time,raypar,tg_radial)

! this is a short intrinsic subroutine that simplifies calls to the ttimes routines
! it takes the angular separation (del), depth of the source(dep) and the phase of
! interest (phase) as input and returns the arrival time (time) and dtdd, the ray parameter
! (raypar) and tg_radial ( = dtdh, the arrival time derivative wrt source depth) as output


      real :: del,dep,time,raypar,tg_radial
      character(len=8) :: phase,n_phase,b_phase,g_phase
      integer  :: mm


! if nodes of bottom interface lie close to the surface, the local ray to the surface will be named
! Pn,Pb,Pg or Sn,Sb,Sg by the ttimes routines, check for this

      n_phase = phase ;b_phase = phase ;g_phase = phase  

      if (phase == 'P') then
         n_phase = 'Pn'
         b_phase = 'Pb'
         g_phase = 'Pg'
      endif

      if (phase == 'S') then
         n_phase = 'Sn'
         b_phase = 'Sb'
         g_phase = 'Sg'
      endif

        phlst(1)=phase
        call brnset(1,phlst,prnt)
        call depset(dep,usrc)
        call trtm(del,maxnp,mm,tt,dtdd,dtdh,dddp,phcd)


        k=1
        phid=0
        do while(k.le.mm.and.phid.eq.0)
           if( phcd(k).eq.phase .or. phcd(k).eq.n_phase &
                .or. phcd(k).eq.b_phase.or. phcd(k).eq.g_phase )then
              phid=k
           else
              k=k+1
           endif
        enddo

        if (phid > 0) then 
           time=tt(phid) 
           raypar=dtdd(phid) 
           tg_radial=dtdh(phid)
        endif

    end subroutine phase_time

end subroutine teleseismic_initialization

!*********************************
subroutine refract_teleseismic_front(reg,vtype,n_tot_ref)

! this subroutine takes the incoming teleseismic wave front normal (time gradient)
! and refracts it through the appropriate side of the grid to initialize the time
! gradient on the regional boundary. Nodes at which the incoming wave front is
! totally reflected are also removed from the initial narrow band.

  use mod_3dfm

  type(Tregion)          :: reg
  integer                :: vtype,n_tot_ref,n,m,inside_hits,late_hits
  real(kind=dp)          :: grad_perp,grad_perp_refracted,grad_par(3),det,normal(3)

  n_tot_ref = 0
  inside_hits=0
  late_hits=0

  do m=1,reg%n_init

     n=reg%init_id(m)     ! take m'th entry from list of nodes to be initialized

  ! do not use teleseismic arrival time if an internal wave arrives at the node earlier
     if (reg%init_arrivaltime(m) > reg%arrivaltime(n)) then
        late_hits=late_hits+1
        cycle
     endif

     select case (reg%init_type(m))   ! get the appropriate inward-pointing normal

        case(1)   ! bottom interface
           normal = intersection(n_intersections)%normal(1:3,reg%node(n)%i3)

        case(2)   ! north side
           normal =(/0.0_dp,-1.0_dp,0.0_dp/)

        case(3)   ! south side
           normal =(/0.0_dp,1.0_dp,0.0_dp/)

        case(4)   ! east side
           normal =(/0.0_dp,0.0_dp,-1.0_dp/)

        case(5)   ! west side
           normal =(/0.0_dp,0.0_dp,1.0_dp/)

        case default
           print *,'illegal node type in refract_teleseismic_front'

     end select

 
! decompose into parallel and perpendicular components

        grad_perp=dot_product(normal,reg%init_time_gradient(1:3,m))

     ! reject if the teleseismic wave hits the boundary from the inside
        if (grad_perp < 0.0_dp) then
           inside_hits=inside_hits+1
           write(15,'(4i5,6f12.3)') reg%init_type(m),reg%node(n)%i1, &
                reg%node(n)%i2,reg%node(n)%i3,reg%init_time_gradient(1:3,m),normal
           cycle
        endif

        grad_par=reg%init_time_gradient(1:3,m)-grad_perp*normal

 ! calculate perpendicular component on the inside of the boundary

        det=1.d0/(reg%velocity(n,vtype)**2) - sum(grad_par**2)

        if (det >= 0.0_dp) then
             
           ! the refracted ray exists

           grad_perp_refracted=sqrt(det)
           reg%time_gradient(1:3,n)=grad_par + grad_perp_refracted*normal
           reg%arrivaltime(n)=reg%init_arrivaltime(m)
           reg%node_status(n)=1

        else

           ! total reflection

           reg%arrivaltime(n)=huge_time
           reg%time_gradient(1:3,n) = 0.0_dp
           n_tot_ref = n_tot_ref + 1
           reg%node_status(n)=-1

        endif

  end do


  if (late_hits>0) print *,'region',reg%id,' :',late_hits,' nodes reached from inside before teleseismic'
  if (inside_hits>0) then
     print *,'region',reg%id,' :',inside_hits,' nodes hit by outgoing teleseismic wave '
     stop
  endif
  if (n_tot_ref > 0) then
     print *,'total reflection of the teleseismic wavefront in region',reg%id,' velocity type',vtype
     print *,'This occurred at',n_tot_ref,' out of',reg%n_init,' nodes'
  endif

  n_tot_ref = n_tot_ref + inside_hits

end subroutine refract_teleseismic_front
vdefs_ak135.f90
module akmod
  logical :: ak_initialized =.false.

  real(kind=8),target               :: depak(136),vpak(136),vsak(136),dak(136)
  real(kind=8),dimension(:),pointer :: vak

contains
 
  subroutine ak_init
    
    if (ak_initialized) return

! read in the ak135 velocity model
    open(1,file='ak135.dat')
    do i=1,136
       read(1,*) depak(i),vpak(i),vsak(i),dak(i)
    end do
    print *,'ak135 models read in'
    ak_initialized=.true.

  end subroutine ak_init

end module akmod



! the functions below return the velocity in each region

function vel1(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel1,dep,lat,long
  integer        :: vtype

  call ak_init   ! check that the model has been read in

  select case (vtype)
     case(1)          
        vak=>vpak   ! if vtype==1 vak points to the P-velocity array
     case(2)
        vak=>vsak   ! if vtype==2 vak points to the S-velocity array
     case default
        stop 'vtype can only be 1 or 2'
  end select

  vel1=vak(1)  ! for the first crustal layer the velocity is constant

  return

end function vel1

function vel2(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel2,dep,lat,long
  integer        :: vtype

  call ak_init
  select case (vtype)
     case(1)
        vak=>vpak
     case(2)
        vak=>vsak
     case default
        stop 'vtype can only be 1 or 2'
  end select

  vel2=vak(3)

  return

end function vel2


function vel3(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel3,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vak=>vpak
     case(2)
        vak=>vsak
     case default
        stop 'vtype can only be 1 or 2'
  end select


! we have to define the velocity outside the region to be constant
! and equal to that on the boundary, to ensure the proper velocity jump
! at the interface. Just interpolating could destroy the discontinuity
! if an interface node does not lie exactly on the discontinuity because
! of rounding errors

  if (dep >= depak(14)) then
     vel3=vak(14)   ! at the bottom of region 3 and below the velocity is constant
     return
  endif
  if (dep <= depak(5)) then
     vel3=vak(5)   ! at the top of region 3 and above the velocity is constant
     return
  endif

  do ii=5,14  ! inbetween do linear interpolation 
     if (dep <= depak(ii+1)) then
        vel3 = vak(ii)+(vak(ii+1)-vak(ii))*(dep-depak(ii))/(depak(ii+1)-depak(ii))
        return
     endif
  end do
  stop 'vel3 error'
end function vel3

function vel4(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel4,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vak=>vpak
     case(2)
        vak=>vsak
     case default
        stop 'vtype can only be 1 or 2'
  end select

  if (dep >= depak(20)) then
     vel4=vak(20)
     return
  endif
  if (dep <= depak(15)) then
     vel4=vak(15)
     return
  endif
  do ii=15,20
     if (dep <= depak(ii+1)) then
        vel4 = vak(ii)+(vak(ii+1)-vak(ii))*(dep-depak(ii))/(depak(ii+1)-depak(ii))
        return
     endif
  end do
  stop 'vel4 error'
end function vel4

function vel5(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel5,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vak=>vpak
     case(2)
        vak=>vsak
     case default
        stop 'vtype can only be 1 or 2'
  end select

  if (dep >= depak(67)) then
     vel5=vak(67)
     return
  endif
  if (dep <= depak(21)) then
     vel5=vak(21)
     return
  endif
  do ii=21,67
     if (dep <= depak(ii+1)) then
        vel5= vak(ii)+(vak(ii+1)-vak(ii))*(dep-depak(ii))/(depak(ii+1)-depak(ii))
        return
     endif
  end do
  stop 'vel5 error'
end function vel5

function vel6(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel6,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vak=>vpak
     case(2)
        vak=>vsak
     case default
        stop 'vtype can only be 1 or 2'
  end select

  if (dep <= depak(68)) then
     vel6=vak(68)
     return
  endif

  do ii=68,135
     if (dep <= depak(ii+1)) then
        vel6 = vak(ii)+(vak(ii+1)-vak(ii))*(dep-depak(ii))/(depak(ii+1)-depak(ii))
        return
     endif
  end do
  stop 'vel6 error'
end function vel6

!-------------------------------------
! below this functions are not used


function vel7(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel7,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel7=3.0d0
     case(2)
        vel7=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel7

function vel8(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel8,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel8=3.0d0
     case(2)
        vel8=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel8

function vel9(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel9,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel9=3.0d0
     case(2)
        vel9=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel9

function vel10(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel10,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel10=3.0d0
     case(2)
        vel10=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel10
vdefs_cmplx.f90
function vel1(dep,lat,long,vtype)

  real(kind=8)  :: vel1,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel1=2.0d0
     case(2)
        vel1=1.2d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel1


function vel2(dep,lat,long,vtype)

  real(kind=8)  :: vel2,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel2=3.0d0
     case(2)
        vel2=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel2


function vel3(dep,lat,long,vtype)

  real(kind=8)  :: vel3,dep,lat,long
  integer        :: vtype

  call random_number(rad)

  select case (vtype)
     case(1)
        vel3=4.0d0 + (dep-5.0d0)*0.05d0 + (rad-0.5d0)*0.2d0
     case(2)
        vel3=0.6*(4.0d0 + (dep-5.0d0)*0.05d0 + (rad-0.5d0)*0.2d0)
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel3


function vel4(dep,lat,long,vtype)

  real(kind=8)  :: vel4,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel4=5.0d0 + (dep-5.0d0)*0.05d0
     case(2)
        vel4=0.6*(5.0d0 + (dep-5.0d0)*0.05d0)
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel4


function vel5(dep,lat,long,vtype)

  real(kind=8)  :: vel5,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel5=8.0d0
     case(2)
        vel5=4.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel5

!--------------------------------------
!functions below this line are not used



function vel6(dep,lat,long,vtype)

  real(kind=8)  :: vel6,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel6=3.0d0
     case(2)
        vel6=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel6

function vel7(dep,lat,long,vtype)

  real(kind=8)  :: vel7,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel7=3.0d0
     case(2)
        vel7=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel7

function vel8(dep,lat,long,vtype)

  real(kind=8)  :: vel8,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel8=3.0d0
     case(2)
        vel8=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel8

function vel9(dep,lat,long,vtype)

  real(kind=8)  :: vel9,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel9=3.0d0
     case(2)
        vel9=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel9

function vel10(dep,lat,long,vtype)

  real(kind=8)  :: vel10,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel10=3.0d0
     case(2)
        vel10=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel10
vdefs_simple.f90
function vel1(dep,lat,long,vtype)

  real(kind=8)  :: vel1,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel1=3.0d0
     case(2)
        vel1=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel1


function vel2(dep,lat,long,vtype)

  real(kind=8)  :: vel2,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel2=4.0d0 + dep*0.005d0
     case(2)
        vel2=0.6*(4.0d0 + dep*0.005d0)
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel2


function vel3(dep,lat,long,vtype)

  real(kind=8)  :: vel3,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel3=5.5 + dep*0.005
     case(2)
        vel3=0.6*(5.5 + dep*0.005)
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel3


function vel4(dep,lat,long,vtype)

  real(kind=8)  :: vel4,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel4=3.0
     case(2)
        vel4=1.8
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel4


function vel5(dep,lat,long,vtype)

  real(kind=8)  :: vel5,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel5=3.0
     case(2)
        vel5=1.8
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel5


function vel6(dep,lat,long,vtype)

  real(kind=8)  :: vel6,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel6=3.0
     case(2)
        vel6=1.8
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel6


function vel7(dep,lat,long,vtype)

  real(kind=8)  :: vel7,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel7=3.0
     case(2)
        vel7=1.8
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel7


function vel8(dep,lat,long,vtype)

  real(kind=8)  :: vel8,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel8=3.0
     case(2)
        vel8=1.8
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel8


function vel9(dep,lat,long,vtype)

  real(kind=8)  :: vel9,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel9=3.0
     case(2)
        vel9=1.8
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel9


function vel10(dep,lat,long,vtype)

  real(kind=8)  :: vel10,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel10=3.0
     case(2)
        vel10=1.8
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel10
visual.f90
!*********************************************************************************
!--------------------------------------------------------------------------------\
!*******************************************************************************************


      subroutine display_interface(isec)

      use mod_3dfm

! all default floats are set to double

      implicit double precision (a-h,o-z)

! argument definition

      type(Tintersection) :: isec

! local array definition

      double precision,dimension(:,:),allocatable::points,centres
      integer,dimension(:,:),allocatable::neighbour,SPfromTR
      integer,dimension(:),allocatable::vis_tlist,vis_elist,add_tlist,hulltriangles
      integer,dimension(:),allocatable::nnn,nnlist,ntrilist
      logical,dimension(:),allocatable::lt_work,ln_work
      double precision,dimension(:,:),allocatable::work_d1,work_d2,work_d3,work_d4
      double precision,dimension(:,:),allocatable::work_d5,work_d6,work_d7
      real,dimension(:,:),allocatable::work_r1
      real,dimension(:),allocatable::work_r2
      integer,dimension(:),allocatable::work_i1,work_i3
      integer,dimension(:,:),allocatable::work_i2


! other variables

      real(kind=dp) :: x,y,z
      integer dmode
      logical clockwise


!---------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! start the triangulation 
!------------------------------------------------------------------------------------


! nl is the initial number of lagrangian points

      nl=isec%nnode

! nlmax is the maximum number of lagrangian points

      nlmax=isec%nnode*2

! ntmax is the maximum number of delaunay triangles built on the lagrangian points

      ntmax=nlmax*3

! nhmax is the maximum number of points on the convex hull formed
! by the lagrangian cloud

      nhmax=nlmax

! npmax is nlmax

      npmax=nlmax

! nnpnmax is the maximum number of neighbours per point

      nnpnmax=50

! nvmax and nmax are the size of working arrays

      nvmax=ntmax
      nmax=3*ntmax+npmax

! eps is a small number

      eps=tiny(eps)

! allocates memory for the delaunay triangulation

     allocate (points(2,nlmax),centres(3,ntmax))
     allocate (neighbour(3,ntmax),SPfromTR(3,ntmax))
     allocate (vis_tlist(nvmax),vis_elist(nvmax),add_tlist(nvmax))
     allocate (hulltriangles(nhmax),nnn(npmax+1),nnlist(nmax),ntrilist(nmax))
     allocate (lt_work(ntmax),ln_work(nlmax))
     allocate (work_d1(2,nnpnmax),work_d2(2,nnpnmax),work_d3(2,nnpnmax))
     allocate (work_d4(2,nnpnmax),work_d5(2,nnpnmax),work_d6(2,nnpnmax))
     allocate (work_d7(2,nnpnmax),work_r1(nnpnmax,2),work_r2(nnpnmax))
     allocate (work_i1(nnpnmax),work_i2(2,nnpnmax),work_i3(nnpnmax))





!     print *,'before triangulation'
!--------------------------------------------------------------------------------------


! initialise variables for nn_setup

         dmode=-2
         nmode=0
         clockwise=.true.
         nohalt_hull=0
         loc=1
         SPfromTR=0
         neighbour=0
         points=0.d0
         field=0.d0
         centres=0.d0
         hulltriangles=0
         nnn=0
         nnlist=0
         ntrilist=0
         vis_tlist=0
         vis_elist=0
         add_tlist=0
         lt_work=.false.
         ln_work=.false.

         do i=1,nl
            points(1,i)=isec%lat(i)
            points(2,i)=isec%long(i)
         enddo

! nn_setup (calculates delaunay triangles IN SURFACE COORDINATES RSPOS and other info)

         call nn2d_setup &
              (nl,ntmax,nhmax,npmax,nnpnmax,nmax, &
              points,dmode,nmode,clockwise,isec%lat(:),nt,SPfromTR, &
              centres,neighbour,nh,hulltriangles,nohalt_hull, &
              loc,nnn,nnlist,ntrilist, &
              eps,nvmax,vis_tlist,vis_elist,add_tlist, &
              lt_work,ln_work)



!  

         SPfromTR=SPfromTR-1   ! counting in DX starts at zero
         value=0.00

         open(1,file='iface'//char(48+isec%id)//'.dx')

         write (1,*) 'object 1 class array type float rank 1 shape 3 items',nl, 'data follows'
         do n=1,nl
            x=isec%r(n)*cos(isec%lat(n))*cos(isec%long(n))
            y=isec%r(n)*cos(isec%lat(n))*sin(isec%long(n))
            z=isec%r(n)*sin(isec%lat(n))
            write(1,'(3f15.6)') y,z,x
         end do
         write (1,*)

         write (1,*) 'object 2 class array type int rank 1 shape 3 items',nt, 'data follows'

         do n=1,nt
            write(1,'(3i6)') (SPfromTR(i,n),i=1,3)
         end do

         write (1,*) 'attribute "element type" string "triangles"'
         write (1,*) 'attribute "ref" string "positions"'
         write (1,*)

         write (1,*) 'object 3 class array type float rank 0 items',nt, 'data follows'
         do n=1,nt
            write(1,'(f15.6)') value
         end do
         write (1,*) 'attribute "dep" string "connections"'
         write (1,*)

         write (1,*) 'object "irregular positions irregular connections" class field'
         write (1,*) 'component "positions" value 1'
         write (1,*) 'component "connections" value 2'
         write (1,*) 'component "data" value 3'
         write (1,*) 'end'

         close(1)

!         print*,'display_interface:',nt,' triangles witten to file iface.dx'

! deallocate everything


      deallocate (points,centres)
      deallocate (neighbour,SPfromTR)
      deallocate (vis_tlist,vis_elist,add_tlist)
      deallocate (hulltriangles,nnn,nnlist,ntrilist)
      deallocate (lt_work,ln_work)
      deallocate (work_d1,work_d2,work_d3,work_d4,work_d5,work_d6,work_d7)
      deallocate (work_r1,work_r2,work_i1,work_i2,work_i3)

      return

      end subroutine display_interface

!*********************************************************************************************************
subroutine display_valid_rays

use mod_3dfm

type(Tray),dimension(:),allocatable  :: disp_ray
integer                              :: j,k,m,n

  j=0
  do n=1,n_receivers
     do m=1,receiver(n)%n_rays
        if (receiver(n)%ray(m)%valid) then
           if (receiver(n)%ray(m)%is_multiray) then
              do k=1,receiver(n)%ray(m)%n_subrays
                 if (receiver(n)%ray(m)%subray(k)%valid) then
                    j=j+1
                 endif
              end do
           else
              j=j+1
           endif
        endif
     end do
  end do


  allocate(disp_ray(j))
  j=0
  do n=1,n_receivers
     do m=1,receiver(n)%n_rays
        if (receiver(n)%ray(m)%valid) then
           if (receiver(n)%ray(m)%is_multiray) then
              do k=1,receiver(n)%ray(m)%n_subrays
                 if (receiver(n)%ray(m)%subray(k)%valid) then
                    j=j+1
                    disp_ray(j)=receiver(n)%ray(m)%subray(k)
                 endif
              end do
           else
              j=j+1
              disp_ray(j)=receiver(n)%ray(m)
           endif
        endif
     end do
  end do

  call display_ray(disp_ray,j)


  deallocate(disp_ray)

  print *,j,'rays displayed'


end subroutine display_valid_rays


!*********************************************************************************************************
subroutine display_ray(ray,nray)

use mod_3dfm

type(Tray) :: ray(nray)
integer,dimension(:),pointer :: vtype_seq
integer i,k,n,m,nsec,nray
real(kind=dp) ::c1,c2,c,x,y,z,r,lat,long


            open(1,file='ray.dx')

            c1=1.0_dp
            c2=1.1_dp
            k=0
            do m=1,nray
               nsec= ray(m)%nsections
               do i=1,nsec ; k=k+ ray(m)%section(i)%npoints ; end do 
            end do

            write (1,*) 'object 1 class array type float rank 1 shape 3 items',k, 'data follows'


            do m=1,nray
               do i=1,ray(m)%nsections
                  do n=1,ray(m)%section(i)%npoints
                     lat=ray(m)%section(i)%point(2,n)
                     long=ray(m)%section(i)%point(3,n)
                     r=ray(m)%section(i)%point(1,n)
                     x=r*cos(lat)*cos(long)
                     y=r*cos(lat)*sin(long)
                     z=r*sin(lat)
                     write(1,'(3f15.6)') y,z,x
                  end do
               end do
            end do

            write (1,*)


            write (1,*) 'object 2 class array type float rank 0 items',k, 'data follows'

            do m=1,nray
               vtype_seq => ray(m)%source%path(ray(m)%raypath_id)%vtype_sequence
               do i=1,ray(m)%nsections
                  if (vtype_seq(i) == 1) c=c1
                  if (vtype_seq(i) == 2) c=c2
                  do n=1,ray(m)%section(i)%npoints
                     write(1,'(3f15.6)') c
                  end do
               end do
            end do

            write (1,*) 'attribute "dep" string "positions"'
            write (1,*)

            write (1,*) 'object "irregular positions irregular connections" class field'
            write (1,*) 'component "positions" value 1'
            write (1,*) 'component "data" value 2'
            write (1,*) 'end'
            
            close(1)

end subroutine display_ray
!*********************************************************************************************************
subroutine display_stored_rays

  use mod_3dfm

  integer n
  real(kind=dp) ::c,x,y,z

  open(1,file='ray.dx')  
  open(41,file='raypos')
  open(42,file='raytype')

  write (1,*) 'object 1 class array type float rank 1 shape 3 items', &
       raypoint_counter, 'data follows'

  do n=1,raypoint_counter
     read(41,*) y,z,x
     write(1,'(3f15.6)') y,z,x
  end do

  write (1,*)


  write (1,*) 'object 2 class array type float rank 0 items', &
       raypoint_counter, 'data follows'

  do n=1,raypoint_counter
     read(42,*) c
     write(1,'(3f15.6)') c
  end do


  write (1,*) 'attribute "dep" string "positions"'
  write (1,*)

  write (1,*) 'object "irregular positions irregular connections" class field'
  write (1,*) 'component "positions" value 1'
  write (1,*) 'component "data" value 2'
  write (1,*) 'end'
            
  close(1)

end subroutine display_stored_rays


!*********************************************************************************************************
subroutine store_ray(ray)

use mod_3dfm

type(Tray) :: ray
integer,dimension(:),pointer :: vtype_seq
integer i,k,n,m,nsec,nray
real(kind=dp) ::c1,c2,c,x,y,z,r,lat,long


            c1=1.0_dp
            c2=1.1_dp
            k=0
            nsec= ray%nsections
            do i=1,nsec 
               raypoint_counter=raypoint_counter+ ray%section(i)%npoints 
            end do


            do i=1,ray%nsections
               do n=1,ray%section(i)%npoints
                  lat=ray%section(i)%point(2,n)
                  long=ray%section(i)%point(3,n)
                  r=ray%section(i)%point(1,n)
                  x=r*cos(lat)*cos(long)
                  y=r*cos(lat)*sin(long)
                  z=r*sin(lat)
                  write(41,'(3f15.6)') y,z,x
               end do
            end do

            vtype_seq => ray%source%path(ray%raypath_id)%vtype_sequence
            do i=1,ray%nsections
               if (vtype_seq(i) == 1) c=c1
               if (vtype_seq(i) == 2) c=c2
               do n=1,ray%section(i)%npoints
                  write(42,'(3f15.6)') c
               end do
            end do

end subroutine store_ray

!
!************************************************************************************
!

subroutine display_sources

use mod_3dfm
implicit none

integer n,nsloc
real(kind=dp) ::ttrue,x,y,z

            open(1,file='sources.dx')


            nsloc=count(source(1:n_sources)%is_local)

            write (1,*) 'object 1 class array type float rank 1 shape 3 items',nsloc, 'data follows'


            do n=1,n_sources
               if (source(n)%is_local) then
                  x=source(n)%r*cos(source(n)%lat)*cos(source(n)%long)
                  y=source(n)%r*cos(source(n)%lat)*sin(source(n)%long)
                  z=source(n)%r*sin(source(n)%lat)
                  write(1,'(3f15.6)') y,z,x
               endif
            end do

            write (1,*)


            write (1,*) 'object 2 class array type float rank 0 items',nsloc, 'data follows'

            ttrue=1.0_dp
            do n=1,nsloc
               write(1,'(3f15.6)') ttrue
            end do

            write (1,*) 'attribute "dep" string "positions"'
            write (1,*)

            write (1,*) 'object "irregular positions irregular connections" class field'
            write (1,*) 'component "positions" value 1'
            write (1,*) 'component "data" value 2'
            write (1,*) 'end'
            
            close(1)

end subroutine display_sources
!*********************************************************************************************************
subroutine display_receivers

use mod_3dfm


integer n
real(kind=dp) ::ttrue,x,y,z

            open(1,file='receivers.dx')

 
            write (1,*) 'object 1 class array type float rank 1 shape 3 items',n_receivers+1, 'data follows'

               x=receiver(1)%r*cos(receiver(1)%lat)*cos(receiver(1)%long)
               y=receiver(1)%r*cos(receiver(1)%lat)*sin(receiver(1)%long)
               z=receiver(1)%r*sin(receiver(1)%lat)
               write(1,'(3f15.6)') y,z,x

            do n=1,n_receivers
               x=receiver(n)%r*cos(receiver(n)%lat)*cos(receiver(n)%long)
               y=receiver(n)%r*cos(receiver(n)%lat)*sin(receiver(n)%long)
               z=receiver(n)%r*sin(receiver(n)%lat)
               write(1,'(3f15.6)') y,z,x
            end do

            write (1,*)


            write (1,*) 'object 2 class array type float rank 0 items',n_receivers+1, 'data follows'

            ttrue=1.0_dp
            do n=1,n_receivers+1
               write(1,'(3f15.6)') ttrue
            end do

            write (1,*) 'attribute "dep" string "positions"'
            write (1,*)

            write (1,*) 'object "irregular positions irregular connections" class field'
            write (1,*) 'component "positions" value 1'
            write (1,*) 'component "data" value 2'
            write (1,*) 'end'
            
            close(1)

end subroutine display_receivers

!*********************************************************************************************************
subroutine display_nodes(reg,display)

use mod_3dfm

type(Tregion) :: reg
integer n
real(kind=dp) ::ttrue,x,y,z
logical ::display(reg%nnode)

            open(1,file='nodes.dx')

            nnode=count(display)
 
            write (1,*) 'object 1 class array type float rank 1 shape 3 items',nnode, 'data follows'

            do n=1,reg%nnode
               if(display(n)) then
                  x=reg%r(n)*cos(reg%lat(n))*cos(reg%long(n))
                  y=reg%r(n)*cos(reg%lat(n))*sin(reg%long(n))
                  z=reg%r(n)*sin(reg%lat(n))
                  write(1,'(3f15.6)') y,z,x
               endif
            end do

            write (1,*)


            write (1,*) 'object 2 class array type float rank 0 items',nnode, 'data follows'

            ttrue=1.0_dp
            do n=1,nnode
               write(1,'(3f15.6)') ttrue
            end do

            write (1,*) 'attribute "dep" string "positions"'
            write (1,*)

            write (1,*) 'object "irregular positions irregular connections" class field'
            write (1,*) 'component "positions" value 1'
            write (1,*) 'component "data" value 2'
            write (1,*) 'end'
            
            close(1)

end subroutine display_nodes
!*********************************************************************************************************
subroutine display_vectors(reg,vectors,display)

use mod_3dfm

type(Tregion) :: reg
integer n
real(kind=dp) ::x,y,z,vectors(3,reg%nnode),xyz(3),a(3,3)
logical ::display(reg%nnode)

            open(1,file='vectors.dx')

            nnode=count(display)
 
            write (1,*) 'object 1 class array type float rank 1 shape 3 items',nnode, 'data follows'

            do n=1,reg%nnode
               if(display(n)) then
                  x=reg%r(n)*cos(reg%lat(n))*cos(reg%long(n))
                  y=reg%r(n)*cos(reg%lat(n))*sin(reg%long(n))
                  z=reg%r(n)*sin(reg%lat(n))
                  write(1,'(3f15.6)') y,z,x
               endif
            end do

            write (1,*)


            write (1,*) 'object 2 class array type float rank 1 shape 3 items',nnode, 'data follows'


            do n=1,reg%nnode
               if(display(n)) then

     a(1,1)=cos(reg%lat(n))*cos(reg%long(n)) ; a(1,2)=-sin(reg%lat(n))*cos(reg%long(n)) ; a(1,3)=-sin(reg%long(n))
     a(2,1)=cos(reg%lat(n))*sin(reg%long(n)) ; a(2,2)=-sin(reg%lat(n))*sin(reg%long(n)) ; a(2,3)= cos(reg%long(n))
     a(3,1)=sin(reg%lat(n))                  ; a(3,2)= cos(reg%lat(n))                  ; a(3,3)= 0.0_dp

                  xyz = matmul(a,vectors(1:3,n))

                  write(1,'(3f15.6)') xyz(2),xyz(3),xyz(1)
               endif

            end do

            write (1,*) 'attribute "dep" string "positions"'
            write (1,*)

            write (1,*) 'object "irregular positions irregular connections" class field'
            write (1,*) 'component "positions" value 1'
            write (1,*) 'component "data" value 2'
            write (1,*) 'end'
            
            close(1)

end subroutine display_vectors
write_plot.f90
program write_plot

  logical head,diff
  real(kind=4)slat,slong,rlat,rlong,delta,cazim,bazim,azima

  open(1,file='arrivals.dat')
  open(12,file='p35')
  open(13,file='p410')
  open(14,file='p660')
  open(15,file='r410')
  open(16,file='r660')

  slat=1.0
  slong=1.0

  do i=1,119
 
     rlat=1.0+i*0.2
     rlong=rlat
    
     call ydist(slat,slong,rlat,rlong,delta,cazim,bazim,azima)

     do j=12,16

        read(1,*) i1,i2,i3,i4,time,diff,head
        if (time==-1.0) cycle

        itype=0
        if (diff) itype=1
        if (head) itype=2

        write(j,*) delta,time,itype

     end do

  end do

  close(1)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
        
end program write_plot
