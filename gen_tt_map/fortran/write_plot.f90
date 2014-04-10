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
