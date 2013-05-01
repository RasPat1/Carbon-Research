       program angtobohr
c
c  convert coord in ang to coord in bohr for
c
      parameter(bohr=0.52917725) 
      character*80 title
      character*3 atom
c
      open (unit=9,form='formatted',status='unknown',file='ang.in') 
      open (unit=10,form='formatted',status='unknown',file='bohr.out')  
c
      read (9,*) natm
      write(10,*) natm
      read (9,'(a80)') title
      write(10,'(a80)') title
c
      do i=1,natm
      read(9,*) atom,x,y,z
      write(10,'(a3,2x,3f10.5)') atom, x/bohr,y/bohr,z/bohr
      enddo
c
      end
