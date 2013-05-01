      program fullopenc60
      
      implicit double precision (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c
c   short time approximation with out importance sampling
c
c   potential  h2 mol -  c60 open cage fullerene -  n- h2 molecules
c   99 atoms  80 carbon atoms   ********* c60 *********
c   c(1-80) h(81-94) n(95-96) 0(97-98) s(99)  ********* c60 *********
c
c  !   H2-C potential fitted 3 site leonard-jones
c  !   Minzhong, Sebastianelli,Gibbons,Bacic,Lawler,Turro
c  !   J. Chem. Phys. 130, 224306 (2009)
c  !   h2-h, h2-N, h2-O, h2-S  charmm potential shufeng
c
c     h2-h2 potential diep and johnson fit to ab-initio  cal
c     JCP vol 112, 4465 (2000)
c
c     subroutine poten(x,vpot) computes h2h2 potential for npair pairs
c     xsite(3,natmx) contains 1 configuration with 2*npair h atoms
c     pairs and the interaction of each atom of a pair with the C of
c     the nanotube
c
c   reconfiguration of weights - weight1 of kevin
c   similar to algorithm of caffarel
c
c   calculates both growth estimator of energy and average of potential
c
c   constraint version
c
c   this version does block averages and saves the potenial after each
c   time step in v1 and v2
c
c
      parameter (npmx=40,npairmx=20,ntubemx=600,nconmx=3000)
c
      parameter (pi = 3.14159265358979323846d0)
      parameter (bohr = 0.529177249d0)
      parameter (hartree2k = 315777.4464d0)
	  parameter (hartree2wn=219474.63d0)
	  parameter (deg=(two*anine*ten)/pi)
      parameter (fk2wavenumber=0.695038770d0)
c	  
      common /param/dt,hb,etrial,nconf,npart
      common /potcom/ hbm(npmx),ipair(npairmx,2),atomwt(npmx),npair
      common /concom/ tol,r2con(3*npmx),ijcon(2,3*npmx),ncon
      common /tube/ tube(3,ntubemx),natmtube,nc
      common /hchnos/ epsc,sigmac,wscale,epsh,sigmah,epsn,sigman,epso,sigmao,epss,sigmas
c      
      dimension x1(nconmx,3,npmx),x2(nconmx,3,npmx),v1(nconmx),v2(nconmx)
      dimension xsite(3,npmx)
      dimension xx(3,npmx)
      dimension irn(4),eblk(50)
      logical skip(2,npmx)
c     real*4 t1,t2,time,tsec
c     real*4 DTIME,tarray(2)
      character*10 eunit
      character*4 atom(npmx)
      character*4 chartemp(ntubemx)      
      character*80 tlabel
c
c common block variables
c /param/
c
c     H = -D*d2 /dx2 + v(x) for particle i
c     D(i)=hbsq/2m(i)=hb/M(i)=hbm(i), M=atomic weight in atomic units
c     where hb=hbar**2/(2.*mass proton) in energy units eunit
c     v(x) is sum of pair potentials
c
c     dt = time step
c     The propagator is
c     G=(1/4*pi*D*dt)**0.5*exp((x'-x))**2/4*D*dt)
c     weight=exp(-(v(x)+v(x'))/2-etrial)*dt)
c     etrial= trial energy adjusted after each call to step
c           =(etrial+eg1)/2
c 
c     nconf = number of walkers which is constant
c     step= number of steps to take
c     nblk= number of blks
c     neq = number of equilibrating blks
c     irn = random number seed
c     nwt10 = number of weights equal or greater then 10 in tstep
c
      open (unit=5,form='formatted',status='unknown',file='nano.in')
      open (unit=6,form='formatted',status='unknown',file='nano.out')
      open (unit=9,form='formatted',status='unknown',file='nano.cnf')
      open (unit=11,form='formatted',status='unknown',file='tube.in')
      open (unit=10,form='formatted',status='unknown',file='nano.jmol')
c
      read (5,'(a80)') tlabel
      read (5,'(4i4)') irn
      read (5,'(i10,f10.5,a10)') npart,hb,eunit
      npair=npart/2
      read (5,'(4i10)') nstep,nblk,nblkeq,nconf
      read (5,'(2i10)') isite
      read (5,'(2f10.5)') dt,delta
      read (5,'(f10.5)') etrial
c
      do 10 i=1,npart
      read (5,'(6x,a4,f10.5,3f10.5)') atom(i),atomwt(i),(xsite(ic,i),ic=1,3)
   10 continue
c
      read (5,'(i10,e20.10)') ncon,tol
      if (ncon.ne.0) then
         read (5,'(2i5)') (ijcon(1,i),ijcon(2,i),i=1,ncon)
         do 100 k=1,ncon
         i=ijcon(1,k)
         j=ijcon(2,k)
         r2con(k)=(xsite(1,i)-xsite(1,j))**2
     +      +(xsite(2,i)-xsite(2,j))**2
     +      +(xsite(3,i)-xsite(3,j))**2
  100    continue
         endif
c
      write (6,'(a80)') tlabel
      write (6,'(t10,''random number seed ='',t40,4i4)') irn
      write (6,'(t10,''number of particles ='',t40,i10)') npart
      write (6,'(t10,''number of pairs ='',t40,i10)') npair
      write (6,'(t10,''hbar**2/(2.*m) ='',t40,f10.5)') hb
      write (6,'(t10,''all energies are in'',t40,a10)') eunit
      write (6,'(t10,''number of configurations ='',t40,i10)') nconf
      write (6,'(t10,''no. of energies/block='',t40,i10)') nstep
      write (6,'(t10,''no. of blocks after eq.='',t40,i10)') nblk
      write (6,'(t10,''no. of blocks before eq. ='',t40,i10)') nblkeq
      write (6,'(t10,''time step ='',t40,f10.5)') dt
      write (6,'(t10,''delta ='',t40,f10.5)') delta
      write (6,'(t10,''etrial ='',t40,f15.5)') etrial
c
      write (6,'(//3x,''i'',5x,''atomic mass'',10x,''site'')')
      write (6,'(i4,2x,a4,f10.5,3f10.5)') (i,atom(i),atomwt(i)
     +   ,(xsite(ic,i),ic=1,3),i=1,npart)
c
      write (6,'(t10,''no. of constraints'',i10)') ncon
      write (6,'(t10,''constrain distances to '',e20.10)') tol
      write (6,'(t10,''pairs to constrain and distance'')')
      write (6,'(t10,2i10,f10.5)')
     +    (ijcon(1,i),ijcon(2,i),sqrt(r2con(i)),i=1,ncon)
      write (6,'(//1x)')
c
c      H atom C  H N O S atom LJ potential shufeng charmm
c
C These parameters have come directly from Shufeng 
      epsc=2.98866684d0
C Used in J. phys. Chem. A Vol 114 No 36 2010 page 9939
      sigmac=5.57469167d0
      wscale=7.5d0
c
        sigmaH=1.89315699d0
        sigmaN=3.49253373d0
        sigmaO=3.24000056d0
        sigmaS=3.74506690d0
        epsH=16.0885000d0
        epsN=33.5468415d0
        epsO=25.9852717d0
        epsS=50.3202623d0
c        
      write (6,'(/1x)')
      write (6,*) '   h2-c,h,n,o,s lj parameters shufeng charmm'
      write (6,'(/1x)') 
c      
      write (6,*) '    epsc= ',epsc, ' cm-1'
      write (6,*) '  sigmac= ',sigmac,' bohr'      
      write (6,*) ' wscale= ',wscale
      write (6,*) '    epsh= ',epsh, ' cm-1'
      write (6,*) '  sigmah ',sigmah,' bohr'
      write (6,*) '    epsn= ',epsn, ' cm-1'
      write (6,*) '  sigman ',sigman,' bohr'
      write (6,*) '    epso= ',epso, ' cm-1'
      write (6,*) '  sigmao ',sigmao,' bohr'
      write (6,*) '    epss= ',epss, ' cm-1'
      write (6,*) '  sigmas ',sigmas,' bohr'
      write (6,'(//1x)')
c
c     read in coords of nanotube
c
      read (11,*) natmtube,nc
      read (11,'(a80)') tlabel
      write (6,'(a80)') tlabel
      write (6,*) '  natmtube= ',natmtube
      write (6,*) '  ncarbon= ',nc
      write (6,'(//1x)')
c      
      do i=1,natmtube
c
      read (11,*) chartemp(i),tube(1,i),tube(2,i),tube(3,i)
c
      enddo
c      
      call setrn(irn)   
c
c load hbm and and pair table - npair h2 mol
c
      do 6 i=1,npart
    6 hbm(i)=hb/atomwt(i)
c
      iatm=1
c
      do i=1,npair
         do j=1,2
          ipair(i,j)=iatm
          iatm=iatm+1
         enddo
      enddo
c     
c  calc input potential
c
      call poten(xsite,vpot)
      write(6,*)'     vpot from sites = ',vpot, 'wavenumber'
c
       call flush(6)
c
c
c     stop 'test poten'
c
c check sites flag if on then get initial configuration from
c sites routine
c
      if (isite.eq.0) then
c
c sites flag not on so get configuration from unit 9. if unable
c then get from sites
c
         rewind 9
c
          read (9,'(i10)') n1
        if(n1.ne.nconf)then
         write(6,*)' tape9 n1 neq nconf - n1= ',n1
         stop 'tape9'
        endif
        read (9,'(3e20.10)')
     +    (((x1(k,ic,i),ic=1,3),i=1,npart),k=1,n1)
        write (6,'(/,t10,i6,'' initial configurations read from 
     +   unit 9'')') n1
c
        rewind 9      
c     
      else
c
c     load x(nconf,3,natmx) from sites and shake up the starting positions of the atoms
c
        n1=nconf
c        
        do 13 k=1,n1
c        
        do 11 i=1,npart
        do 11 ic=1,3
   11   xx(ic,i)=xsite(ic,i)+(rannyu(0)-half)*delta
c
c       call constr(xsite,xx,npart,hbm,r2con,ijcon,ncon,tol,niter)        
c
        do 12 i=1,npart
        do 12 ic=1,3
c       
   12   x1(k,ic,i)=xx(ic,i)
c   
   13   continue
c
      endif

c
      write (6,'(/,t10,i6,'' initial configurations from sites'')') n1
c
c load v1 from sites or tape9
c
      do 250 k=1,n1
c
      do 240 i=1,npart
      do 240 j=1,3
      xx(j,i)=x1(k,j,i)
  240 continue
c
      call poten(xx,vv)
      v1(k)=vv
c
  250 continue
c 
c zero out estimators
c engcum is the energy estimator
c eng2cum will accumulate squared averages for variance
c vvcum and vv2cum potential estimators
c
      engcum=zero
      eng2cum=zero
      enpcum=zero
      enp2cum=zero
      nbav=0
c
c     loop to do short time iterations
c
c     call second(t1)
c     t1=SECNDS(0.0)
c     tsec=DTIME(tarray)
c
      write (6,*)' start iterations'
c
      nstep=nstep/2
c
      write(6,*)' *************************************************'
c
      do 25 j=1,nblk
c
      eng=zero
      eng2=zero
      enp=zero
      enp2=zero
      nstep2=0
c
      do 20 i=1,nstep
c
c walk configurations on x1 and put them on x2
c
c     write (6,*)' dt gauss ',dt,gauss(0),gauss(0)
c
      call tstep(x1,x2,v1,v2,eg1,vpot1)
c
c     adjust etrial
c
      etrial=(eg1+etrial)/two
c
c     walk from 2 back to 1
c
      call tstep(x2,x1,v2,v1,eg2,vpot2)
c
c     adjust etrial
c
      etrial=(eg2+etrial)/two
c
c     calculate energies
c
      eng=eng+eg1+eg2
      eng2=eng2+eg1**2+eg2**2
      enp=enp+vpot1+vpot2
      enp2=enp2+vpot1**2+vpot2**2
      nstep2=nstep2+2
c
  20  continue
c
c write out block averages
c
      engav=eng/nstep2
      eng2av=eng2/nstep2
      sig=sqrt(abs(eng2av-(engav)**2)/nstep2)
      enpav=enp/nstep2
      enp2av=enp2/nstep2
      sigp=sqrt(abs(enp2av-(enpav)**2)/nstep2)
c      
      write(6,*)' block =',j
      write(6,'(1x,"energy by growth estimate = ",f12.5," +- "
     + ,f12.5)')engav,sig
      write(6,'(1x,"energy by potential estimate = ",f12.5," +- "
     + ,f12.5)') enpav,sigp    
c
      write(6,*)' *************************************************'
c
      call flush (6)
c
      if(j.gt.nblkeq)then
        nbav=nbav+1
        engcum=engcum+engav
        eng2cum=eng2cum+engav**2
        eblk(nbav)= engav
        enpcum=enpcum+enpav
        enp2cum=enp2cum+enpav**2
        
      endif
c
   25 continue
c
c write out final averages
c
      engcum=engcum/nbav
      eng2cum=eng2cum/nbav
      sigblk=sqrt(abs(eng2cum-(engcum)**2)/(nbav-1))
      enpcum=enpcum/nbav
      enp2cum=enp2cum/nbav
      sigpblk=sqrt(abs(enp2cum-(enpcum)**2)/(nbav-1))      
c
      write(6,*)' average over blocks '
      write(6,'(1x,"growth energy estimate = ",f12.5," +- ",f12.5)')
     + engcum,sigblk
      write(6,'(1x,"potential energy estimate = ",f12.5," +- ",f12.5)')
     + enpcum,sigpblk     
c
      write(6,'(//)')

c
c      call second(t2)
c      time=(t2-t1)/60.
c
c     tsec=DTIME(tarray)
c     time=tsec/60.0
c
      call CPU_TIME(time)
      timeout=time/60.0
c
       write(6,'(10x,''time='',f10.3,2x,''min'')') timeout
c
c  write out configurations to tape9
c
      rewind 9
      write (9,'(i10)') n1
         do kk=1,n1
          write (9,'(3e20.10)')((x1(kk,ii,jj),ii=1,3),jj=1,npart)
         enddo
c
c  write out every tenth configuration for jmol
c
       do 90 k=1,n1,10
c
      write (10,'(i5)') npart+natmtube
      write (10,*)' h2-h2 fullerene c60 open short time'
c      note !!! change to angstroms      
      write (10,'(a4,3e20.10)') (atom(i), (x1(k,ic,i)*bohr,ic=1,3),
     x       i=1,npart)
c      
      do i=1,natmtube
c
      write (10,*) chartemp(i),tube(1,i)*bohr,tube(2,i)*bohr,tube(3,i)*bohr
c
      enddo
c        
c      
  90  continue
c
c compute autocorrelation for energies after neq iterations
c
c     write(6,*)' *************************************************'
c
c     call autocor(eblk,nbav)
c
c
      write(6,'(//)')
c
      stop 'fullopenspce - short time- no importance sampling- weight1'
      end
c      
      subroutine poten(x,vpot)
      implicit double precision (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1e0,half=.5e0,third=1.d0/3.d0)
c
c     compute h atom from each pair interaction with nanotube - open fullerene
c     h atom c LJ potential - charmm parameters h-c,h,n,o,s
c
c ****  c60 open fullerene version - must change for c70 ***********
c
c     h2-h2 potential diep and johnson fit to ab-initio  cal
c     JCP vol 112, 4465 (2000)
c
c     cal_jacobi (c routine) converts cartesian to jacoby coords - fortran version
c     francesco
c     pqs containts R, theta1, theta2, tau
c
c     subroutine to compute h2h2 potential for npair pairs
c     x(3,natmx) contains 1 configuration with 2*npair h atoms
c
      parameter (pi = 3.14159265358979323846d0)
      parameter (bohr = 0.529177249d0)
      parameter (hartree2k = 315777.4464d0)
	  parameter (hartree2wn=219474.63d0)
	  parameter (fk2wavenumber=0.695038770d0)
c
      parameter (npmx=40,npairmx=20,ntubemx=600,nconmx=3000)
c      
      common /param/  npart
      common /potcom/ hbm(npmx),ipair(npairmx,2),atomwt(npmx),npair
      common /tube/ tube(3,ntubemx),natmtube,nc
      common /hchnos/ epsc,sigmac,wscale,epsh,sigmah,epsn,sigman,epso,sigmao,epss,sigmas
c
      dimension x(3,npmx),coord(3,4),cmass(4),pqs(4),hcm(3)
c
      vpot=zero
c
c     compute h2 nano interaction
c
      do i=1, npair
c
      i1=ipair(i,1)
      i2=ipair(i,2)
c
c      write (6,*) ' i i1 i2 ', i, i1,i2
c      write (6,*) '      '
c
       do ii=1,3
         hcm(ii)=(x(ii,i1)+x(ii,i2))*half
       enddo
c       
          do j=1,natmtube
c
           r1sq=zero
           r2sq=zero
           rchcmsq=zero
c
              do k=1,3
c
                 r1sq=r1sq+(tube(k,j)-x(k,i1))**2
                 r2sq=r2sq+(tube(k,j)-x(k,i2))**2
                 rchcmsq=rchcmsq+(tube(k,j)-hcm(k))**2
c
c             write (6,*) ' j k rsq1 tube(k,j)  x(k,i1) ',j,k,r1sq,tube(k,j),x(k,i1)
c
              enddo
c
           r1=sqrt(r1sq)
           r2=sqrt(r2sq)
           rchcm=sqrt(rchcmsq)
c
        if(j.lt.81) then
           v1=vlj(r1,sigmac,epsc)
           v2=vlj(r2,sigmac,epsc)
           v3=vlj(rchcm,sigmac,epsc)
c
c         write (6,*) ' j  r1  r2  v1  v2 ', j,r1,r2,v1,v2
c
           vv=v1+v2+wscale*v3
c           
        elseif (j.lt.95) then
           v1=vlj(r1,sigmah,epsh)
           v2=vlj(r2,sigmah,epsh)
           vv=v1+v2
c
        elseif (j.lt.97) then
           v1=vlj(r1,sigman,epsn)
           v2=vlj(r2,sigman,epsn)
           vv=v1+v2
c
        elseif (j.lt.99) then
           v1=vlj(r1,sigmao,epso)
           v2=vlj(r2,sigmao,epso)
           vv=v1+v2
c
       elseif (j.lt.100) then
           v1=vlj(r1,sigmas,epss)
           v2=vlj(r2,sigmas,epss)
           vv=v1+v2
c           
       endif           
c          
           vpot=vpot+vv
c
         enddo
c
c        write (6,*) ' pot h2-c pair j vv ',j,vv
c        write (6,*) '      '
c         
      enddo
c
c     write (6,*) ' vpot h2-c,h,n,o,s  total ', vpot, ' wave number '
c     write (6,*) '  '
c          
c     compute h2-h2 interaction
c 
      if (npair.ne.1) then
c      
      do i=1, npair-1
         do j=i+1,npair
c           
              i1=ipair(i,1)
              i2=ipair(i,2)
              i3=ipair(j,1)
              i4=ipair(j,2)
c
              cmass(1)=atomwt(i1)
              cmass(2)=atomwt(i2)
              cmass(3)=atomwt(i3)
              cmass(4)=atomwt(i4)
c
              do k=1,3
c
! for JULES: I changed the following four lines (your original lines commented). 
! Now coor dimensions are always (3,4)
!                 coord(1,k)=x(k,i1)
!                 coord(2,k)=x(k,i2)
!                 coord(3,k)=x(k,i3)
!                 coord(4,k)=x(k,i4)
!  change below!!!!!!!!!!!!!!!!!!!!!!!
                  coord(k,1)=x(k,i1)
                  coord(k,2)=x(k,i2)
                  coord(k,3)=x(k,i3)
                  coord(k,4)=x(k,i4)
c
              enddo
c
              call cal_jacobi(cmass,coord,pqs)
c
c             change hartee to wavenumer
c                  
                  vh2h2=h2h2pot(pqs(1),pqs(2),pqs(3),pqs(4))
                  vpot=vpot+vh2h2*hartree2wn
c
c            write(6,*)' R= ',pqs(1), ' bohr' 
c            write(6,*)' theta1= ',pqs(2), ' rad'
c            write(6,*)' theta2= ',pqs(3), ' rad'
c            write(6,*)' tau= ',pqs(4),' rad '
c            write(6,*) '           coords'
! for JULES, again here switching for coord(3,4): see coord(ic,iii)!!!!!
c           do iii=1,4
c            write(6,*) (coord(ic,iii), ic=1,3)
c           end do          
c           write(6,*)' '
c
         enddo
      enddo
c
      endif
c
      return
      end
      double precision function vlj(r,sigma,eps)
      implicit double precision (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1e0,half=.5e0,third=1.d0/3.d0)
c.      
      common /hchnos/ epsc,sigmac,wscale,epsh,sigmah,epsn,sigman,epso,sigmao,epss,sigmas
c      
c     LJ potential sig bohr  eps cm-1
c
      r6i=(sigma/r)**6
      r12i=r6i**2
      vlj=four*eps*(r12i-r6i)
      return
      end       
      double precision function gauss(x)
      implicit double precision (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c
c sample gaussian by box muller method
c samples exp(-x**2/2)
c
      parameter (npmx=40,npairmx=20,ntubemx=600,nconmx=3000)
      parameter (pi = 3.14159265358979323846d0)
      common /param/dt,hb,etrial,nconf,npart
      common /potcom/ hbm(npmx),ipair(npairmx,2),atomwt(npmx),npair
c
      gauss=cos(rannyu(0)*two*pi)*sqrt(-log(rannyu(0)))
      return
      end
      subroutine tstep(xold,xnew,vvold,vvnew,egrow,vpot)
      implicit double precision (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c
c routine to take time step dt and propagage walkers from xold to xnew
c
      parameter (npmx=40,npairmx=20,ntubemx=600,nconmx=3000)
      common /param/dt,hb,etrial,nconf,npart
      common /potcom/ hbm(npmx),ipair(npairmx,2),atomwt(npmx),npair

      common /concom/ tol,r2con(3*npmx),ijcon(2,3*npmx),ncon
      dimension xold(nconmx,3,npmx),xnew(nconmx,3,npmx)
      dimension vvold(nconmx),vvnew(nconmx)
      dimension xo(3,npmx),xn(3,npmx)
      dimension w(nconmx),iwt(nconmx)
      logical skip(2,npmx)
c
c     initialize pointer
c
      nnew=0
      vpot=zero
c
c     write(6,*)'  nold etrial ',nold,etrial
c
      do 15 i=1, nconf
c
c move walker - one configuration - dx from short time green's function
c
      do 5 m=1,npart
      do 5 ic=1,3
      xo(ic,m)=xold(i,ic,m)
      dx=sqrt(four*hbm(m)*dt)*gauss(xx)
      xn(ic,m)=xo(ic,m)+dx
   5  continue
c
c rigid constraints
c
      if (ncon.ne.0) then
         call constr(xo,xn,npart,hbm,r2con,ijcon,ncon,tol,niter,skip)
         endif
c
c     calculate w(i)
c
c     write(6,*)' xo '
c     write(6,*) ((xo(ic,jj),ic=1,3),jj=1,npart)
c     write(6,*)' xn '
c     write(6,*) ((xn(ic,jj),ic=1,3),jj=1,npart)
c
c     call poten(xo,vold)
c
      vold=vvold(i)
      call poten(xn,vnew)
c
      u=half*(vold+vnew)-etrial
c
      w(i)=exp(-u*dt)
c
c     nowt can only be 1 if nstate=1 and wt to be removed
c
      if(nowt.eq.1)then
      w(i)=zero
      endif
c     
      do 10 m=1,npart
      do 10 ic=1,3
      xold(i,ic,m)=xn(ic,m)
  10  continue
c
      vvold(i)=vnew
c
  15  continue
c
c     reconfigure weights kevin weight1
c
      call weight1(w,nconf,iwt,wtsum)
c
      do 50 i=1, nconf
c
      ipart=iwt(i)
c
c     write(6,*)' vold vnew u wt ipart ',vold,vnew,u,wt,ipart
c
      if (ipart.gt.0) then
c
        do 25 j=1,ipart
        nnew=nnew+1
c
        do 20 m=1,npart
        do 20 ic=1,3
        xnew(nnew,ic,m)=xold(i,ic,m)
  20    continue
c
        vvnew(nnew)=vvold(i)
        vpot=vpot+vvnew(nnew)
c
  25    continue
c
      endif
c
c     write(6,*)' ipart nnew ',ipart,nnew
c
  50  continue
c
      if(nconf.ne.nnew)then
         write(6,*)' nconf= ',nconf,' nnew= ',nnew
         stop 'tstep'
      endif
c
      egrow=log(nconf/wtsum)/dt+etrial
      vpot=vpot/nconf
c
      return
      end
      double precision function rannyu(x)
      implicit double precision (a-h,o-z)
      parameter(twom12=0.000244140625d0)
      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4
c
c generic statement functions
c
      ishft12(ii)=ii/4096
      mask12(ii)=mod(ii,4096)
c
c unix f77 statement functions
c
c     ishft12(ii)=rshift(ii,12)
c     mask12(ii)=and(ii,4095)
c
      i1=l1*m4+l2*m3+l3*m2+l4*m1
      i2=l2*m4+l3*m3+l4*m2
      i3=l3*m4+l4*m3
      i4=l4*m4
      l4=mask12(i4)
      i3=i3+ishft12(i4)
      l3=mask12(i3)
      i2=i2+ishft12(i3)
      l2=mask12(i2)
      l1=mask12(i1+ishft12(i2))
      rannyu=twom12*(l1+
     +       twom12*(l2+
     +       twom12*(l3+
     +       twom12*(l4))))
      return
      end
      subroutine setrn(iseed)
      implicit double precision (a-h,o-z)
      common /rnyucm/ m(4),l(4)
      integer iseed(4)
c
      do 10 i=1,4
      l(i)=mod(iseed(i),4096)
   10 continue
      l(4)=2*(l(4)/2)+1
      return
      end
      subroutine savern(iseed)
      implicit double precision (a-h,o-z)
      common /rnyucm/ m(4),l(4)
      integer iseed(4)
      do 10 i=1,4
      iseed(i)=l(i)
   10 continue
      return
      end
      block data random
      common /rnyucm/ m(4),l(4)
c
      data m / 502,1521,4071,2107/
      data l /   0,   0,   0,   1/
      end
c     subroutine autocor(eblk,nblk)
c
c autocorrelation formula from koonin computational physics
c project 8 - monte carlo solution of h2 molecule
c subroutine crltns
c
c     parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
c     parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
c     parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c
c     dimension eblk(50),c(0:15)
c
c     kmax=6
c
c     do 30 k=0,kmax 
c     ei=zero
c     eik=zero
c     esqi=zero
c     esqik=zero
c     eiek=zero
c     imax=nblk-k
c     do 20 i=1,imax
c     ei=ei+eblk(i)
c     eik=eik+eblk(i+k)
c     esqi=esqi+eblk(i)**2
c     esqik=esqik+eblk(i+k)**2
c     eiek=eiek+eblk(i)*eblk(i+k)
c 20  continue
c     fim=float(imax)
c     ei=ei/fim
c     eik=eik/fim
c     esqi=esqi/fim
c     esqik=esqik/fim
c     eiek=eiek/fim
c     c(k)=(eiek-ei*eik)/(sqrt(esqi-ei**2))/(sqrt(esqik-eik**2))
c 30  continue
c
c
c     write(6,*)'   autocorrelation  nblk  kmax',nblk,kmax
c     write(6,*)'    k               c(k) '
c     do 50 k=0,kmax
c     write(6,*) k,c(k)
c     write(8,*) k,c(k)
c 50  continue
c     return
c     end
      subroutine constr(x1,x2,npart,hbm,r2con,ijcon,ncon,tol,niter,skip)
c
c fix constraint dynamics. This is a modified version of shake
c by w.f. van gunsteren, groningen, aug. 1981
c description of shake : j. comp. phys. 23 (1977) 327.
c
      implicit double precision (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
      parameter (maxit=1000,small=1.d-6)
      dimension x1(3,npart),x2(3,npart),hbm(npart),r2con(ncon)
      dimension ijcon(2,ncon)
      logical skip(2,npart)
      dimension dx(3)
      logical ready
      niter=0
      tol2=two*tol
c     tol2=2.e0*tol
      do 10 k=1,npart
      skip(1,k)=.true.
      skip(2,k)=.false.
   10 continue
      ready=.false.
      do 20 k=1,maxit
      if (ready) return
      ready=.true.
      do 30 m=1,ncon
      i=ijcon(1,m)
      j=ijcon(2,m)
      if (.not.(skip(2,i).and.skip(2,j))) then
         toler=r2con(m)
         dx(1)=x2(1,i)-x2(1,j)
         dx(2)=x2(2,i)-x2(2,j)
         dx(3)=x2(3,i)-x2(3,j)
         diff=toler-dx(1)**2-dx(2)**2-dx(3)**2
         if (abs(diff).ge.toler*tol2) then
            dx(1)=x1(1,i)-x1(1,j)
            dx(2)=x1(2,i)-x1(2,j)
            dx(3)=x1(3,i)-x1(3,j)
            rrpr=dx(1)**2+dx(2)**2+dx(3)**2
            if (rrpr.lt.toler*small) then
               write (6,'(1x,''deviation too small in constr '')')
               write (6,*) ' i j rrpr ',i,j,rrpr
               stop
               endif
            acor=diff/(rrpr*(hbm(i)+hbm(j))*two)
            dx(1)=dx(1)*acor
            dx(2)=dx(2)*acor
            dx(3)=dx(3)*acor
            x2(1,i)=x2(1,i)+dx(1)*hbm(i)
            x2(2,i)=x2(2,i)+dx(2)*hbm(i)
            x2(3,i)=x2(3,i)+dx(3)*hbm(i)
            x2(1,j)=x2(1,j)-dx(1)*hbm(j)
            x2(2,j)=x2(2,j)-dx(2)*hbm(j)
            x2(3,j)=x2(3,j)-dx(3)*hbm(j)
            skip(1,i)=.false.
            skip(1,j)=.false.
            ready=.false.
            endif
         endif
   30 continue
      niter=niter+1
      do 40 i=1,npart
      skip(2,i)=skip(1,i)
      skip(1,i)=.true.
   40 continue
   20 continue
      write (6,'(1x,i10,'' iterations without convergence in constr'')')
     +   maxit
c
c some compilers complain about the return and some need it
c so put in stupid if statement to make it think it can return
c
      if (0.eq.0) stop
      return
      end

      subroutine weight1(wt,nwalk,iwt,wtsum)
      implicit double precision (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
c
c     from Kevin
c     test branch f77 version
c
      parameter (nmax=3000)
      dimension iwt(nmax),wt(nmax)
c
c wt(i) = input weight
c nwalk = input number of walkers
c iwt(i) = returned integer output weight
c wtsum = returned sum of input weights
c
      do i=1,nwalk
         iwt(i)=0
      enddo
c
      wtsum=zero
      do i=1,nwalk
         wtsum=wtsum+wt(i)
      enddo
      anorm=nwalk/wtsum
      do i=1,nwalk
         wt(i)=wt(i)*anorm
      enddo
c
c take care of weights less than 1 and the first 1 part
c of the weights greater than 1. 
c
      do i=1,nwalk
         if (wt(i).lt.1.0d0) then
            rn=rannyu(0)
            if (rn.lt.wt(i)) then
            iwt(i)=1
         else
            iwt(i)=0
         endif
         wt(i)=0.0d0
      else
         iwt(i)=1
         wt(i)=wt(i)-1.0d0
      endif
      enddo
c
c calculate the number of samples needed from the remaining
c weights and normalize the remaining weights so they add
c to that number
c
      iwtsum=0
      sumwt=zero
      do i=1,nwalk
        iwtsum=iwtsum+iwt(i)
        sumwt=sumwt+wt(i)
      enddo
c
      nsamp=nwalk-iwtsum
      anorm=nsamp/sumwt
      do i=1,nwalk
         wt(i)=wt(i)*anorm
      enddo
c
c get nsamp samples from remaining weights
c
      i=0
      w=0.0e0
      do k=1,nsamp
         rn=rannyu(0)
         do while (w.lt.rn)
            i=i+1
            if (i.gt.nwalk) return
            w=w+wt(i)
         enddo
         w=w-1.0d0
         iwt(i)=iwt(i)+1
      enddo
      return
      end
c***********************************************************************
      function H2H2POT(rr, theta1, theta2, phi12)
c  Calculates the value of the H2-H2 potential at a given orientation 
c  (theta1, theta2, phi12) and with the radial components v000(rr),
c  v022(rr)=v202(rr), and v224(rr).  All center-of-mass to 
c  center-of-mass distances are in bohr units.  Energy units are in
c  Hartree.
c  Based on a fit to ab initio data from JCP 113, 3480-3481 (2000) and
c  JCP 112, 4465-4473 (2000).
c***********************************************************************
      implicit none
      real*8 H2H2POT, rr, theta1, theta2, phi12,
     &       v000, v022, v224
      H2H2POT=v000(rr)+
     & v022(rr)*5.0d0*(3.0d0*(cos(theta2))**2-1.0d0)/2.0d0+
     & v022(rr)*5.0d0*(3.0d0*(cos(theta1))**2-1.0d0)/2.0d0+
     & v224(rr)*45.0d0/sqrt(1120.0d0)*
     & (2.0d0*(3.0d0*(cos(theta1))**2-1.0d0)*
     & (3.0d0*(cos(theta2))**2-1.0d0)-16.0d0*sin(theta1)*cos(theta1)*
     & sin(theta2)*cos(theta2)*cos(phi12)+(sin(theta1))**2*
     & (sin(theta2))**2*cos(2.0d0*phi12))
      return
      end
c********************** end of function potential **********************
      function v000(rr)
c  Computes the radial component v000 of the H2-H2 potential.
c  This is the isotropic contribution.  Distance and energies are in 
c  atomic units (bohr and hartree).
c***********************************************************************
      implicit none
      real*8 rr, v000
      real*8 alpha, beta, gamma, rc, c6, c8, c10, fcut
      parameter(alpha=0.655914d0, beta=1.018447d0, gamma=8.070618E-02)
      parameter(c6=13.076837d0,c8=80.700360d0,c10=3687.082967d0)
      parameter(rc=9.034308d0)
      if (rr .lt. rc) then
         fcut = exp(-(rc/rr-1.0d0)**2)
      else
         fcut = 1.0d0
      endif
      v000 = exp(alpha - beta*rr - gamma*rr*rr) -
     &     fcut*(c6/rr**6 + c8/rr**8 + c10/rr**10)
      return
      end
c***********************************************************************
      function v022(rr)
c  Computes the radial components, v022 and v202, of the H2-H2 
c  potential.  v022 = v202 due to symmetry.
c  Atomic units (bohr and hartree).
c***********************************************************************
      implicit none
      real*8 rr, v022
      real*8 alpha, beta, gamma, rc, c6, c8, c10, fcut
      parameter(alpha=-3.428227d0, beta=0.715011d0, gamma=0.100120d0)
      parameter(c6=0.288396d0,c8=8.242595d0,c10=210.984075d0)
      parameter(rc=8.422755d0)
      if (rr .lt. rc) then
         fcut = exp(-(rc/rr-1.0d0)**2)
      else
         fcut = 1.0d0
      endif
      v022 = exp(alpha - beta*rr - gamma*rr*rr) -
     &     fcut*(c6/rr**6 + c8/rr**8 + c10/rr**10)
      return
      end
c***********************************************************************
      function v224(rr)
c     This is the so-called quadrupole-quadrupole interaction term.
c***********************************************************************
      implicit none
      real*8 rr, v224, epsi
      parameter(epsi=0.135269d0)
      v224 = epsi/rr**5
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine cal_jacobi(mass,coor,pq)
c***********************************************************************
c Need new jacobi for H2-H2 pqs (define diff. in potential)
      Implicit none
      integer i,j,k
      real*8 pi
      Real*8 R,thetHH,thetFH,tau
      Real*8 rFH, rHH, ln1, ln2
      Real*8 sprod2,length2
      Real*8 sign,test
      Real*8 mass,coor,pq
      Real*8 vFH,vHH,vCC
      Real*8 n1,n2,n3
      Real*8 FHC, HHC
      Real*8 vec1,vec2,vec
      Real*8 vect1,vect2,vect3
      parameter (pi=3.1415926535897932384626d0)
c                       coor(x y z , i-th H atom)  
      dimension mass(4),coor(3,4),pq(4)
      dimension vFH(3),vHH(3),vCC(3)
      dimension n1(3),n2(3),n3(3)
      dimension FHC(3), HHC(3)
      dimension vec1(3),vec2(3),vec(3)
      dimension vect1(3),vect2(3),vect3(3)
c
c      (1)   H	
c           /			   H (3)
c	     /			    
c	(2)	F================\
c				 
c					  H (4)
c	pq(1)=R;		[au]
c	pq(2)=thetHH;		[rad]
c	pq(3)=Pi-thetFH;	[rad]
c	pq(4)=tau;		[rad]
c
c	vFH: F->H (1->0)		---   rFH
c	vFHC:FHC->H
c	vHH: H->H (2->3)		---   rHH
c	vHHC:HHC->H
c	vCC: HFC->HHC			---   R
c	n1: vCC*vFH			---   Ln1
c	n2: vCC*vHH			---   Ln2
c	thetFH: vFH|vCC
c	thetHH: vHH|vCC
c	tau	: n1|n2
c	FHC = Coordinate center of mass HF
c	HHC = Coordinate center of mass HH
c
c Center of mass coordinates:
      do j=1,3
         FHC(j)=(mass(1)*coor(j,1)+mass(2)*coor(j,2))/(mass(1)+mass(2))
         HHC(j)=(mass(3)*coor(j,3)+mass(4)*coor(j,4))/(mass(3)+mass(4))
      end do
c Calculate Vectors
      do j=1,3
         vFH(j)=coor(j,2)-coor(j,1)
         vHH(j)=coor(j,4)-coor(j,3)
         vCC(j)=HHC(j)-FHC(j)
      end do
c Get normal vectors
      call vprod2(vCC,vHH,n1)
      call vprod2(vCC,vFH,n2)
      call vprod2(n2,n1,n3)
      sign=sprod2(n3,vCC)
c Vector lengths
      rFH=length2(vFH)
      rHH=length2(vHH)
      R=length2(vCC)
      ln1=length2(n1)
      ln2=length2(n2)
c Angles
      thetFH=acos(sprod2(vFH,vCC)/(rFH*R))
      thetHH=acos(sprod2(vHH,vCC)/(rHH*R))
      test=sprod2(n1,n2)/(ln1*ln2)
      if(test.lt.-1.0) then
        test=-1.0
      end if
      if(test.gt.1.0) then
        test=1.0
      end if
      tau=acos(test)
      if(sign.lt.0.0) then
        tau=2.0*Pi-tau
      end if
c Final Quantities
      pq(1)=R
      pq(2)=thetHH
      pq(3)=Pi-thetFH
c#######################################################################
c     below IF ELSE: 
c     I added the following
c     to avoid phi=NAN when vCC parallel to vH1 or vH2 
      IF(sign.eq.0.0)THEN
        tau=0.0
        pq(4)=tau
      ELSE
        pq(4)=tau
      END IF
c########################
c     write(*,*)
c     write(*,*)'pq(1)=R= ', pq(1),' bohr'
c     write(*,*)'pq(2)=thetHH= ', pq(2),' rad'
c     write(*,*)'pq(3)=Pi-thetFH= ',pq(3),' rad'
c     write(*,*)'pq(4)=tau= ',pq(4),' rad'
      return 
      end
c***********************************************************************
      function sprod2(vec1,vec2)
c***********************************************************************
      implicit none
      integer i,j,k
      real*8 sprod2
      real*8 vec1,vec2
      dimension vec1(3),vec2(3)
      sprod2=0
      do j=1,3
      sprod2=vec1(j)*vec2(j)+sprod2
      end do
      return
      end
c***********************************************************************
      function length2(vec)
c***********************************************************************
      implicit none
      integer i,j,k
      real*8 length2
      real*8 vec
      dimension vec(3)
      length2=SQRT(vec(1)**2+vec(2)**2+vec(3)**2)
      return 
      end
c***********************************************************************
      subroutine vprod2(vect1,vect2,vect3)
c***********************************************************************
      implicit none
      integer i,j,k
      real*8 vect1,vect2,vect3
      dimension vect1(3),vect2(3),vect3(3)
      vect3(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
      vect3(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
      vect3(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)
      return
      end
C***********************************************************************


      
