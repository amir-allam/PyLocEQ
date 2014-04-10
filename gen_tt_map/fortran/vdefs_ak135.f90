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
