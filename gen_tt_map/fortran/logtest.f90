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
