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
