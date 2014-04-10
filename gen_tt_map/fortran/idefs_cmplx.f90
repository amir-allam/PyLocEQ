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
