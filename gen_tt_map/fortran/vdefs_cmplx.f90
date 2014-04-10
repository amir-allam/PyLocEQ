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
