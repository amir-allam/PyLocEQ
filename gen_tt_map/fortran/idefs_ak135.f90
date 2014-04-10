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
