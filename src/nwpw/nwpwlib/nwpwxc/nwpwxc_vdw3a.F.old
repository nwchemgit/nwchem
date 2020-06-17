C dftd3 program for computing the dispersion energy and forces from cartesian atomic coordinates
C and atomic numbers as described in
C
C S. Grimme, J. Antony, S. Ehrlich and H. Krieg
C A consistent and accurate ab initio parameterization of density functional dispersion correction
C (DFT-D) for the 94 elements H-Pu
C J. Chem. Phys, 132 (2010), 154104
C 
C if BJ-damping is used 
C S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456-1465.
C
C should be cited as well.
C
C Copyright (C) 2009 - 2011 Stefan Grimme, University of Muenster, Germany
C
C This program is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 1, or (at your option)
C any later version.

C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.

C For the GNU General Public License, see <http://www.gnu.org/licenses/>

 
      subroutine  nwpwxc_vdw3_dftd3(options,n,iz,xyz,lat,eout,g,g_lat)
      implicit none             
      character*(*) options
c number of atoms
      integer n
c cardinal numbers of elements
      integer iz(*)
c coordinates in au
      real*8 xyz(3,n)
c lattice in au
      real*8 lat(3,3)

      real*8 eout
c gradient
      real*8 g(3,n)
      real*8 g_lat(3,3)


      integer istart,iend
      integer maxat,max_elem,maxc                      
c conversion factors
      real*8 autoang,autokcal,c6conv,autoev
      parameter (maxat   =20000)
      parameter (max_elem=94)
c maximum coordination number references per element
      parameter (maxc    =5)
c coversion factors
      parameter (autoang =0.52917726d0)
      parameter (autokcal=627.509541d0)
      parameter (autoev  = 27.21138505)

c DFT-D version
      integer version
!c cardinal numbers of elements
      !integer,dimension(:), allocatable :: iz  
cc fixed atoms in geometry opt
c      logical fix(maxat)
c cut-off radii for all element pairs
      real*8 r0ab(max_elem,max_elem)
c C6 for all element pairs 
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
c how many different C6 for one element
      integer mxc(max_elem)
c C6810 
      real*8 c6,c8,c10
c coordination numbers of the atoms
      real*8,dimension(:), allocatable :: cn  
c covalent radii
      real*8 rcov(max_elem)
c atomic <r^2>/<r^4> values
      real*8 r2r4(max_elem)
c energies
      real*8 e6, e8, e10, e12, disp, e6abc        
c THE PARAMETERS OF THE METHOD (not all a "free")
      real*8 rs6, rs8, rs10, s6, s18, alp6, alp8, alp10, s42, rs18, alp
c printout option
      logical echo
c grad ?
      logical grad
c analyse results ?
      logical anal
c third-order term?
      logical noabc
c gradient calctype
      logical numgrad
c special parameters
      logical tz
c periodic boundary conditions
      logical pbc
c repetitions of the unitcell to match the rthr and c_thr
      integer rep_vdw(3),rep_cn(3)
c R^2 distance neglect threshold (important for speed in case of large systems)
      real*8 rthr,rthr2
c R^2 distance to cutoff for CN_calculation
      real*8 cn_thr
c Integer for assigning only max/min cn C6 (0=normal, 1=min, 2=max)
c local and dummy variables
      character*80 atmp,btmp,ctmp,dtmp,etmp,ftmp,func
      character*2  nwpwxc_esym 
      integer i,j,k,z,nn,iat,jat,i1,i2
      integer ida(max_elem),ipot
      real*8  x,y,dispr,displ,gdsp,dum,xx(10),dum6(86)
      real*8  dum1,dum2,dum3(3)
      logical ex,pot,fdum
      logical minc6list(max_elem),maxc6list(max_elem),minc6,maxc6

c PBE0/def2-QZVP atomic values 
!      data r2r4 /
!     .  8.0589,  3.4698, 29.0974, 14.8517, 11.8799,  7.8715,  5.5588,
!     .  4.7566,  3.8025,  3.1036, 26.1552, 17.2304, 17.7210, 12.7442,
!     .  9.5361,  8.1652,  6.7463,  5.6004, 29.2012, 22.3934, 19.0598,
!     . 16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569,
!     . 10.1428,  9.4907, 13.4606, 10.8544,  8.9386,  8.1350,  7.1251,
!     .  6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451,
!     . 11.0355, 10.1997,  9.5414,  9.0061,  8.6417,  8.9975, 14.0834,
!     . 11.8333, 10.0179,  9.3844,  8.4110,  7.5152, 32.7622, 27.5708,
!     . 23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,
!     . 15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,
!     . 14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671,  8.3549,
!     .  7.8496,  7.3278,  7.4820, 13.5124, 11.6554, 10.0959,  9.7340,
!     .  8.8584,  8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305,
!     . 15.6161, 15.1226, 16.1576 /                                       


c scale r4/r2 values of the atoms by sqrt(Z) 
c sqrt is also globally close to optimum
c together with the factor 1/2 this yield reasonable
c c8 for he, ne and ar. for larger Z, C8 becomes too large
c which effectively mimics higher R^n terms neglected due
c to stability reasons
      
!  r2r4 =sqrt(0.5*r2r4(i)*dble(i)**0.5 ) with i=elementnumber
!  the large number of digits is just to keep the results consistent
!  with older versions. They should not imply any higher accuracy than
!  the old values
      data r2r4 /
     . 2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594,
     . 3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516,
     . 6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576,
     . 4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947,
     . 6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167,
     . 5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141,
     . 6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647,
     . 4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917,
     . 6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424,
     . 5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523,
     . 5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549,
     .10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807,
     . 8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454,
     . 8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339,
     . 7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381,
     . 6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695,
     . 7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318,
     . 6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068,
     . 8.77140725,  8.65402716,  8.53923501,  8.85024712 /

c covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
c values for metals decreased by 10 %
!      data rcov/
!     .  0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67
!     ., 1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54
!     ., 1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09
!     ., 1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39
!     ., 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26
!     ., 1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57
!     ., 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53
!     ., 1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32
!     ., 1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58
!     ., 1.52, 1.53, 1.54, 1.55 /

! these new data are scaled with k2=4./3.  and converted a_0 via
! autoang=0.52917726d0
      data rcov/
     . 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
     . 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
     . 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
     . 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
     . 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
     . 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
     . 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
     . 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
     . 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
     . 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
     . 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
     . 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
     . 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
     . 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
     . 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
     . 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
     . 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
     . 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
     . 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

c k1-k3
      !include 'param'

      real*8 k1,k2,k3

c global ad hoc parameters
      parameter (k1=16.0)
      parameter (k2=4./3.) 

c reasonable choices are between 3 and 5
c this gives smoth curves with maxima around the integer values
c k3=3 give for CN=0 a slightly smaller value than computed
c for the free atom. This also yields to larger CN for atoms
c in larger molecules but with the same chem. environment
c which is physically not right
c values >5 might lead to bumps in the potential
      parameter (k3=-4.) 

      logical  inp_strtok
      external inp_strtok
      integer  ga_nodeid
      external ga_nodeid



!      write(*,'(94(F12.8,'',''))')r2r4
!      stop

c scale and convert to au
!      rcov=k2*rcov/autoang
!      write(*,'(94(F11.8,'',''))')rcov
!      stop
!      do i=1,max_elem
!         dum    =0.5*r2r4(i)*dble(i)**0.5   
c store it as sqrt because the geom. av. is taken
!         r2r4(i)=sqrt(dum)                         
!      enddo
c init
      echo=.true. 
      grad=.false.
      pot =.false.
      anal=.false.
      noabc=.true. 
      numgrad=.false.
      tz=.false.
      func=' none (read from parameter file)'
      version=0
      pbc=.false.
      !fix=.false.
      minc6=.false.
      maxc6=.false.
      minc6list=.false.
      maxc6list=.false.
      fdum=.false.
c Cutoff r^2 thresholds for the gradient in bohr^2.
c rthr influences N^2 part of the gradient.
c rthr2 influences the N^3 part of the gradient. When using
c dftd3 in combination with semi-empirical methods or FFs, and large
c (>1000 atoms) systems, rthr2 is crucial for speed:
c Recommended values are 20^2 to 25^2 bohr.

      rthr=9000.0d0   ! UR, SE
      rthr2=1600.0d0
      cn_thr=1600.0d0

c J/mol nm^6 - > au
      c6conv=1.d-3/2625.4999d0/((0.052917726d0)**6)

c set radii
c     call nwpwxc_rdab('~/.r0ab.dat',autoang,max_elem,r0ab)
      call nwpwxc_setr0ab(max_elem,autoang,r0ab)

c read C6 file by default from $HOME     
c     btmp='~/.c6ab.dat'
c     inquire(file=btmp,exist=ex)
c Muenster directory as second default
c     if(.not.ex)btmp='/usr/qc/dftd3/c6ab.dat'
c     call nwpwxc_loadc6(btmp,maxc,max_elem,c6ab,mxc)         

      

cc get coord filename
c      call getarg(1,etmp)
c      inquire(file=etmp,exist=ex)
c      if(.not.ex) call nwpwxc_printoptions       
      ex=.false.
      ipot=0

c options
      istart = 0
      do while (inp_strtok(options,' ',istart,iend))
         ftmp = options(istart:iend)
         if(index(ftmp,'-h')      .ne.0) call nwpwxc_printoptions
         if(index(ftmp,'-grad'   ).ne.0) grad=.true. 
         if(index(ftmp,'-anal'   ).ne.0) anal=.true. 
         if(index(ftmp,'-noprint').ne.0) echo=.false.
         if(index(ftmp,'-abc'    ).ne.0) noabc=.false.
         if(index(ftmp,'-pbc'    ).ne.0) pbc=.true.
         if(index(ftmp,'-num'    ).ne.0) numgrad=.true.
         if(index(ftmp,'-tz')     .ne.0) tz=.true.
         if(index(ftmp,'-old')    .ne.0) version=2
         if(index(ftmp,'-zero')   .ne.0) version=3
         if(index(ftmp,'-bj')     .ne.0) version=4
         if(index(ftmp,'-min')    .ne.0) then
           minc6=.true.
           j=0
           DO while (inp_strtok(options,' ',istart,iend))
             !call getarg(i+j+1,atmp)
             atmp = options(istart:iend)
  
             if (index(atmp,'-').eq.0.and.atmp.ne.'') then
               call nwpwxc_elem(atmp,nn)
               if (nn.gt.max_elem.or.nn.lt.1)
     .           call nwpwxc_stoprun('Could not recognize min Element')
               minc6list(nn)=.true.
               j=j+1
             else
               exit
             endif
           ENDDO
         endif
         if(index(ftmp,'-max')    .ne.0) then
           maxc6=.true.
           k=0
           DO while (inp_strtok(options,' ',istart,iend))
             !call getarg(i+k+1,atmp)
             atmp = options(istart:iend)
             if (index(atmp,'-').eq.0.and.atmp.ne.'') then
               call nwpwxc_elem(atmp,nn)
               if (nn.gt.max_elem.or.nn.lt.1)
     .           call nwpwxc_stoprun('Could not recognize max Element')
               maxc6list(nn)=.true.
               k=k+1
             else
               exit
             endif
           ENDDO
         endif
c         if(index(ftmp,'-pot')    .ne.0) then
c                                         pot=.true.
c                                         call getarg(i+1,atmp)
c                                         call nwpwxc_readl(atmp,xx,nn)
c                                         ipot=idint(xx(1))
c                                         endif
c         if(index(ftmp,'-cnthr')   .ne.0) then
c                                         call getarg(i+1,atmp)
c                                         call nwpwxc_readl(atmp,xx,nn)
c                                         rthr2=xx(1)
c                                         rthr2=rthr2**2
c                                         endif
         if(index(ftmp,'-func')  .ne.0)  then
            !call getarg(i+1,func)
            if (inp_strtok(options,' ',istart,iend)) then
               func = options(istart:iend)
            end if
            ex=.true.
         endif



c         if(index(ftmp,'-cutoff') .ne.0) then
c                                         call getarg(i+1,atmp)
c                                         call nwpwxc_readl(atmp,xx,nn)
c                                         rthr=xx(1)**2
c                                         endif
      end do

      echo = echo.and.(ga_nodeid().eq.0)




c      do i=1,iargc()
c      call getarg(i,ftmp)
c      if(index(ftmp,'-h')      .ne.0) call nwpwxc_printoptions
c      if(index(ftmp,'-grad'   ).ne.0) grad=.true. 
c      if(index(ftmp,'-anal'   ).ne.0) anal=.true. 
c      if(index(ftmp,'-noprint').ne.0) echo=.false.
c      if(index(ftmp,'-abc'    ).ne.0) noabc=.false.
c      if(index(ftmp,'-pbc'    ).ne.0) pbc=.true.
c      if(index(ftmp,'-num'    ).ne.0) numgrad=.true.
c      if(index(ftmp,'-tz')     .ne.0) tz=.true.
c      if(index(ftmp,'-old')    .ne.0) version=2
c      if(index(ftmp,'-zero')   .ne.0) version=3
c      if(index(ftmp,'-bj')     .ne.0) version=4
c      if(index(ftmp,'-min')    .ne.0) then
c        minc6=.true.
c        j=0
c        DO 
c          call getarg(i+j+1,atmp)
c          if (index(atmp,'-').eq.0.and.atmp.ne.'') then
c            call nwpwxc_elem(atmp,nn)
c            if (nn.gt.max_elem.or.nn.lt.1) 
c     .        call nwpwxc_stoprun('Could not recognize min Element')
c            minc6list(nn)=.true.
c            j=j+1
c          else
c            exit
c          endif
c        ENDDO
c      endif
c      if(index(ftmp,'-max')    .ne.0) then
c        maxc6=.true.
c        k=0
c        DO 
c          call getarg(i+k+1,atmp)
c          if (index(atmp,'-').eq.0.and.atmp.ne.'') then
c            call nwpwxc_elem(atmp,nn)
c            if (nn.gt.max_elem.or.nn.lt.1) 
c     .        call nwpwxc_stoprun('Could not recognize max Element')
c            maxc6list(nn)=.true.
c            k=k+1
c          else
c            exit
c          endif
c        ENDDO
c      endif
c      if(index(ftmp,'-pot')    .ne.0) then
c                                      pot=.true. 
c                                      call getarg(i+1,atmp)
c                                      call nwpwxc_readl(atmp,xx,nn)
c                                      ipot=idint(xx(1))
c                                      endif
c      if(index(ftmp,'-cnthr')   .ne.0) then
c                                      call getarg(i+1,atmp)
c                                      call nwpwxc_readl(atmp,xx,nn)
c                                      rthr2=xx(1)
c                                      rthr2=rthr2**2
c                                      endif
c      if(index(ftmp,'-func')  .ne.0)  then
c                                      call getarg(i+1,func)
c                                      ex=.true.
c                                      endif
c       
c
c
c      if(index(ftmp,'-cutoff') .ne.0) then
c                                      call getarg(i+1,atmp)
c                                      call nwpwxc_readl(atmp,xx,nn)
c                                      rthr=xx(1)**2
c                                      endif
cc      if(index(ftmp,'-pot')    .ne.0) then
c      enddo
c

c     Check command line input


      if(minc6.and.j.lt.1)then
       call nwpwxc_stoprun('No Element given for min/max')
      endif
      if(maxc6.and.k.lt.1)then
       call nwpwxc_stoprun('No Element given for min/max')
      endif
      do i=1,max_elem

        if (minc6list(i).and.maxc6list(i)) 
     .    call nwpwxc_stoprun('Unreasonable min/max input!')
!      if (minc6list(i)) write(*,*)'min:',i
!      if (maxc6list(i)) write(*,*)'max:',i
      enddo
c C6 hard-coded (c6ab.dat not used)
c this is alternative to nwpwxc_loadc6
      call nwpwxc_copyc6(btmp,maxc,max_elem,c6ab,mxc,
     .            minc6,minc6list,maxc6,maxc6list)   
       cn_thr=rthr2

!       write(*,*)'CN(P):',c6ab(15,15,mxc(15),1,2)
!       write(*,*)'mxc(P):',mxc(15)

c       if (pbc) then
c         call nwpwxc_pbcrdatomnumber(etmp,n)
c       else
c         call nwpwxc_rdatomnumber(etmp,n)
c       endif

!      allocations
      !allocate(xyz(3,n))
      !allocate(g(3,n))
      !allocate(iz(n))
      allocate(cn(n))

c reading coordinates and cell in VASP.5.2-format
c determing repetitions of unitcell
      if (pbc) then
              !call nwpwxc_pbcrdcoord(etmp,lat,n,xyz,iz,autoang)
              call nwpwxc_set_criteria(rthr,lat,dum3)
              rep_vdw=int(dum3)+1
              call nwpwxc_set_criteria(cn_thr,lat,dum3)
              rep_cn=int(dum3)+1
c              write(*,*)'VDW-cutoff:',sqrt(rthr)*autoang
c              write(*,*)'CN-cutoff :',sqrt(cn_thr)*autoang
c              write(*,*)'repvdw:',rep_vdw
c              write(*,*)'repcn :',rep_cn
      else !no pbc
c read coordinates, either TM or xmol file
              !call nwpwxc_rdcoord(etmp,n,xyz,iz,fix,fdum)
      endif !pbc
      if(n.lt.1)     call nwpwxc_stoprun( 'no atoms' )
      if(n.gt.maxat) call nwpwxc_stoprun( 'too many atoms' )


     
c the analytical E(3) grad is not available yet
      if(grad.and.(.not.noabc))numgrad=.true.
     
c set parameters for functionals
      if(ex) then
         call nwpwxc_setfuncpar(func,version,tz,s6,rs6,s18,rs18,alp)
      else
         call nwpwxc_rdpar(dtmp,version,s6,s18,rs6,rs18,alp)
      endif

      if(echo)then
      write(*,*)' _________________________________'
      write(*,*)'                                  '
      write(*,*)'|         DFTD3 V3.0 Rev 1        |'
      write(*,*)'| S.Grimme, University Bonn       |'
      write(*,*)'|           May  2013             |'
      write(*,*)'|   see dftd3 -h for options      |'
      write(*,*)' _________________________________'
      write(*,*)
      write(*,*)'Please cite DFT-D3 work done with this code as:'
      write(*,*)'S. Grimme, J. Antony, S. Ehrlich and H. Krieg,'
      write(*,*)'J. Chem. Phys. 132 (2010), 154104'
      write(*,*)'If used with BJ-damping cite also'
      write(*,*)'S. Grimme, S. Ehrlich and L. Goerigk,'
      write(*,*)'J. Comput. Chem. 32 (2011), 1456-1465'
      write(*,*)'For DFT-D2 the reference is'
      write(*,*)'S. Grimme, J. Comput. Chem., 27 (2006), 1787-1799'
      write(*,*)
      write(*,*)' files read :     '
      !write(*,*)trim(etmp)       
      if(.not.ex)write(*,*)trim(dtmp)       
      endif

      if(version.lt.2.or.version.gt.4)stop 'inacceptable version number'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C all calculations start here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c CNs for output
      if (pbc) then
          call nwpwxc_pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,cn_thr)
      else
          call nwpwxc_ncoord(n,rcov,iz,xyz,cn,cn_thr)
      endif

      if(version.eq.2)then
          if(echo)write(*,'(''loading DFT-D2 parameters ...'')')
          call nwpwxc_loadoldpar(autoang,max_elem,maxc,c6ab,r0ab,dum6)
c number of CNs for each element
          mxc=1
convert to au
          c6ab=c6ab*c6conv
      endif


c which atoms are present? (for printout)
      if(echo)then
      ida=0
      do i=1,n
         ida(iz(i))=ida(iz(i))+1
      enddo
      write(*,'(''C6 coefficients used:'')')
      do i=1,94
         if(ida(i).gt.0)then
            write(*,*) mxc(i),' C6 for element ',i
            do j=1,maxc
               if(c6ab(i,i,j,j,1).gt.0)then
               write(*,'(''Z='',i3,'' CN='',F6.3,5x,''C6(AA)='',F8.2)')
     .         i,c6ab(i,i,j,j,2),c6ab(i,i,j,j,1)
               endif
            enddo
         endif
      enddo
      endif
      
c output
      if (echo) then
          write(*,'(/''#               XYZ [au]  '',12x,
     .              '' R0(AA) [Ang.]''2x,
     .              ''CN'',7x,
     .              ''C6(AA)     C8(AA)   C10(AA) [au] '')
     .            ')
          x=0
          btmp=''
          do i=1,n
          z=iz(i)
          call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(i),iz(i),
     >                      cn(i),cn(i),c6)
          do j=1,n
          call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(i),iz(j),
     >                      cn(i),cn(j),dum)
          x=x+dum
          enddo
c compute C8/C10 for output
          c8 =r2r4(iz(i))**2*3.0d0*c6     
          c10=(49.0d0/40.0d0)*c8**2/c6
          dum=0.5*autoang*r0ab(z,z)
          if(version.eq.4)dum=rs6*0.5*autoang*sqrt(c8/c6)
          atmp=' '
          !if(fix(i)) then
          ! atmp='f'
          ! btmp='f'
          !endif
          write(*,'(i4,3F10.5,3x,a2,1x,a1,F7.3,2x,F7.3,3F12.1,L)')
     .    i,xyz(1:3,i),nwpwxc_esym(z),atmp,
     .    dum,cn(i),
     .    c6,c8,c10
          enddo
          write(*,'(/''molecular C6(AA) [au] = '',F12.2)')x
          if(btmp.eq.'f') then
          write(*,*)'  '
           write(*,*)'Caution: Some coordinates fixed 
     .in gradient (marked f, see above).'
          write(*,*)'  '
          endif
          if(fdum)then
          write(*,*)'Caution: Dummy atoms found and ignored.'
          endif
   
      endif

      
c testoutput of radii
c     do i=1,94
c        call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,i,i,0.d0,0.0d0,c6)
c        c8 =r2r4(i)**2*3.0d0*c6     
c        write(22,*) i, sqrt(c8/c6)
c     enddo
c     write(22,*)
c     do i=1,94
c        write(22,*) i, r0ab(i,i)  
c     enddo
c     stop

c for global ad hoc parameters see
c k3 in subroutine nwpwxc_getc6, k1 and k2 in subroutine nwpwxc_ncoord*
c fixed or dependent ones:
      rs8  = rs18       
      rs10 = rs18
      alp6 = alp
      alp8 = alp+2.
      alp10= alp8+2. 
c note: if version=4 (Becke-Johnson), a1=rs6 and a2=rs18
c       and alp* have no meaning

c*********************************************************************
c*********************************************************************
c testing code
c output of C6=f(CN)
      if(pot.and.ipot.gt.100)then
      x=0
      do i=1,100
      call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,ipot-100,ipot-100,
     .                              x,x,C6)
      write(2,*) x,c6
      x=x+0.05
      enddo
      stop
      endif
c Edisp pot curve for testing. Caution: C6 is not constant along R!
      if(pot)then
      write(*,*) 'Computing Edisp potential curve for atom ',ipot
      xyz=0
      iz(1)=ipot
      iz(2)=ipot
      n=2
      xyz(3,2)=1.0/autoang
 142  if (pbc) then
       call nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .           e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)
      else
        call nwpwxc_edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >     r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr,
     .     e6,e8,e10,e12,e6abc)
      endif
      xyz(3,2)=xyz(3,2)+0.02
      disp=-s6*e6-s18*e8
      write(42,*) xyz(3,2)*autoang,disp*autokcal
      write(43,*) xyz(3,2)        ,disp*autokcal
      if (pbc) then
          call nwpwxc_pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,cn_thr)
      else
          call nwpwxc_ncoord(n,rcov,iz,xyz,cn,cn_thr)
      endif
      call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(1),iz(2),
     >                  cn(1),cn(2),c6)
      write(2,*)xyz(3,2)*autoang,-autokcal*c6/xyz(3,2)**6
      if(xyz(3,2).lt.20) goto 142
      write(42,*)
      stop 'pot curve done'
      endif
c end testing code
c*********************************************************************
c*********************************************************************

c check if all parameters have been loaded and are resonable
      do iat=1,n-1
         do jat=iat+1,n
            if(r0ab(iz(jat),iz(iat)).lt.0.1) then
               write(*,*) iat,jat,iz(jat),iz(iat)
               call nwpwxc_stoprun( 'radius missing' )
            endif
            if (version.eq.2)then
              c6=c6ab(iz(jat),iz(iat),1,1,1)
            else 
              call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                      cn(iat),cn(jat),c6)
            endif
            if(c6.lt.1.d-6) then
               write(*,*) iat,jat,iz(jat),iz(iat),cn(iat),cn(jat)
               call nwpwxc_stoprun( 'C6 missing' )
            endif
         enddo
      enddo

c sanity check of read coordniates, based on covalnent radii. 
c Not omnipotent but better than nothing. S.E. 15.09.2011
c      call nwpwxc_checkcn(n,iz,cn,c6ab,max_elem,maxc)
      if (pbc) then
        call nwpwxc_pbccheckrcov(n,iz,rcov,xyz,lat)
      else
        call nwpwxc_checkrcov(n,iz,rcov,xyz)
      endif

cccccccccccccc
c energy call
cccccccccccccc
      if (pbc) then
        call  nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,
     >     r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .     e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)
 
      else
        call nwpwxc_edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >     r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr,
     .     e6,e8,e10,e12,e6abc)
      endif

      e6   = e6   *s6

      e6abc= e6abc*s6

      e8   = e8   *s18

      disp =-e6-e8-e6abc

c e10 has been tested once again with BJ-damping but has no good effect
c     e10  = e10  *s18
c     disp =-e6-e8-e10-e6abc

      eout = disp
c output
      if (echo) then
      if(version.lt.4)then
      write(*,'(/10x,'' DFT-D V'',i1)') version       
      else
      write(*,'(/10x,'' DFT-D V3(BJ)'')') 
      endif
      write(*,'('' DF '',a50)') func          
      write(*,'('' parameters'')') 
      if(version.eq.2)then
         write(*,'('' s6       :'',f10.4)') s6            
         write(*,'('' alpha6   :'',f10.4)') alp6        
      endif
      if(version.eq.3)then
         write(*,'('' s6       :'',f10.4)') s6            
         write(*,'('' s8       :'',f10.4)') s18           
         write(*,'('' rs6      :'',f10.4)') rs6  
         write(*,'('' rs18     :'',f10.4)') rs18          
         write(*,'('' alpha6   :'',f10.4)') alp6        
         write(*,'('' alpha8   :'',f10.4)') alp8           
         write(*,'('' k1-k3    :'',3f10.4)') k1,k2,k3     
      endif
      if(version.eq.4)then
         write(*,'('' s6       :'',f10.4)') s6            
         write(*,'('' s8       :'',f10.4)') s18           
         write(*,'('' a1       :'',f10.4)') rs6           
         write(*,'('' a2       :'',f10.4)') rs18          
         write(*,'('' k1-k3    :'',3f10.4)') k1,k2,k3     
      endif
        write(*,'('' Cutoff   :'',f10.4,'' a.u.'')') sqrt(rthr) !*autoang
        write(*,'('' CN-Cutoff:'',f10.4,'' a.u.'')')sqrt(cn_thr)!*autoang
      if (pbc) then
        write(*,'('' Rep_vdw  :'',3I3)') rep_vdw
      endif
      write(*,*)
      if (pbc) then
      write(*,'('' Edisp /kcal,au,eV:'',f11.4,X,f12.8,X,f11.7)') 
     .  disp*autokcal,disp,disp*autoev
      else
        write(*,'('' Edisp /kcal,au:'',f11.4,f12.8)') disp*autokcal,disp
      endif
      write(*,'(/'' E6    /kcal :'',f11.4)')-e6*autokcal
      if(version.gt.2)then
      write(*,'('' E8    /kcal :'',f11.4)')-e8*autokcal 
c     write(*,'('' E10   /kcal :'',f11.4)')-e10*autokcal 
      if(.not.noabc)
     .write(*,'('' E6(ABC) "   :'',2f11.6,F16.12)')-e6abc*autokcal 
      write(*,'('' % E8        :'',f6.2)') -e8/disp/0.01         
      if(.not.noabc)
     .write(*,'('' % E6(ABC)   :'',f6.2)') -e6abc/disp/0.01        
      endif
      endif

c  NOT NEEDED? ...EJB
cc this file for tmer2 read tool
c      open(unit=1,file='.EDISP')
c      write(1,*) disp        
c      close(1)

cccccccccccccccccccccccccc
c analyse Edisp pair terms
cccccccccccccccccccccccccc
      if(anal) then
       if (pbc) then

        call nwpwxc_pbcadisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,
     >           r0ab,rcov,
     .           rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal,autoang,
     .           rthr,rep_vdw,cn_thr,rep_cn,s6,s18,disp*autokcal,lat)
       else
        call nwpwxc_adisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,
     >           r0ab,rcov,
     .           rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal,
     .           autoang,rthr,cn_thr,s6,s18,disp*autokcal)
       endif !pbc
      endif !anal
 
cccccccccccccccccccccccccc
c gradient
cccccccccccccccccccccccccc
      if(grad)then
      g=0.0d0
      call cpu_time(dum1)
      if (pbc) then
      call nwpwxc_pbcgdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,
     .           rcov,s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,numgrad,
     .                 version,g,gdsp,x,g_lat,lat,rep_vdw,rep_cn,
     .                 rthr,echo,cn_thr)
 
      else
c        call nwpwxc_gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
c     .           s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,rthr,
c     .           numgrad,version,echo,g,gdsp,x,rthr2,fix)
        call nwpwxc_gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,rcov,
     .           s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,rthr,
     .           numgrad,version,echo,g,gdsp,x,rthr2)
      endif
      call cpu_time(dum2)
!      echo=.true. !Jonas
      if(echo)write(*,'('' time  '',f6.1)')dum2-dum1
!       if (pbc) call nwpwxc_stresstensor(maxc,max_elem,autoev,s6,s18,xyz,n,iz,
!     .     lat,c6ab,mxc,version,numgrad,noabc,echo,r0ab,r2r4,
!     .     rcov,rs6,rs8,rs10,alp6,alp8,alp10,rthr,rep_vdw,cn_thr,rep_cn,
!     .     g_lat)
      if (echo.and.pbc) then
        write(*,*)'Stresstensor:'
        write(*,*)g_lat


      endif

 
c check if nwpwxc_gdisp yields same energy as edisp
      if(abs((disp-gdsp)/disp).gt.1.d-3) then
         write(*,*) disp,gdsp
         call nwpwxc_stoprun('internal error')
      endif !sanitycheck
c write to energy and gradient files in TM style
      if (pbc) then
        if (echo) then
          write(*,*)'Cartesian gradient written to file dftd3_gradient.'
          write(*,*)'Cartesian cellgradient written 
     . to file dftd3_cellgradient. (a.u.)'
          g_lat=g_lat!*autoev
          call nwpwxc_pbcwregrad(n,g,g_lat)
        endif
      else
        if (echo) then
        write(*,*) 'Cartesian gradient written to file dftd3_gradient'
        call nwpwxc_outg(n,g)
        end if

      endif
      endif !grad
      !deallocate(xyz,g,iz,cn)
      deallocate(cn)

      if(echo)write(*,*) 'normal termination of dftd3'

      return
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_printoptions
      write(*,*) 'dftd3 <coord filename> [-options]'
      write(*,*) 'options:'
      write(*,*) '-func <functional name in TM style>'
      write(*,*) '-grad'
      write(*,*) '-anal (pair analysis)'
      write(*,*) '     file <fragemt> with atom numbers'
      write(*,*) '     is read for a fragement based   '
      write(*,*) '     analysis (one fragment per line,'
      write(*,*) '     atom ranges (e.g. 1-14 17-20) are allowed)'
      write(*,*) '-noprint'
      write(*,*) '-pbc (periodic boundaries; reads VASP-format)'
      write(*,*) '-abc (compute E(3))'
      write(*,*) '-cnthr (neglect threshold in Bohr for CN, default=40)'
      write(*,*) '-cutoff (neglect threshold in Bohr for E_disp, 
     . default=95)'
      write(*,*) '-old (DFT-D2)'         
      write(*,*) '-zero (DFT-D3 original zero-damping)' 
      write(*,*) '-bj   (DFT-D3 with Becke-Johnson finite-damping)' 
      write(*,*) '-tz (use special parameters for TZ-type calculations)'
      write(*,*) 'variable parameters read from ~/.dftd3par.<hostname>'
      write(*,*) 'if -func is used, -zero or -bj or -old is required!"'
      stop
      end subroutine nwpwxc_printoptions

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C set parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C>
C> ### References ###
C>
C> [1] A. Hoffmann, R. Grunzke, S. Herres-Pawlis,
C>     "Insights into the influence of dispersion correction in the
C>      theoretical treatment of guanidine-quinoline copper(I)
C>      complexes", J. Comput. Chem. (2014) <b>35</b>, pp. 1943-1950,
C>     DOI:
C>     <a href="https://doi.org/10.1002/jcc.23706">
C>     10.1002/jcc.23706</a>.
C>
      subroutine nwpwxc_setfuncpar(func,version,TZ,s6,rs6,s18,rs18,alp)
      implicit none  
      integer version
      real*8 s6,rs6,s18,alp,rs18
      character*(*) func     
      logical TZ
c double hybrid values revised according to procedure in the GMTKN30 paper

c DFT-D3 with Becke-Johnson finite-damping, variant 2 with their radii 
c SE: Alp is only used in 3-body calculations
      if(version.eq.4)then
      s6=1.0d0
      alp =14.0d0

      select case (func)
         case ("b-p")
              rs6 =0.3946
              s18 =3.2822
              rs18=4.8516
         case ("b-lyp")
              rs6 =0.4298
              s18 =2.6996
              rs18=4.2359
         case ("revpbe")
              rs6 =0.5238
              s18 =2.3550
              rs18=3.5016
         case ("rpbe")
              rs6 =0.8318
              s18 =0.1820
              rs18=4.0094
         case ("b97-d")
              rs6 =0.5545
              s18 =2.2609
              rs18=3.2297
         case ("pbe")
              rs6 =0.4289
              s18 =0.7875
              rs18=4.4407
         case ("rpw86-pbe")
              rs6 =0.4613
              s18 =1.3845
              rs18=4.5062
         case ("b3-lyp")
              rs6 =0.3981
              s18 =1.9889
              rs18=4.4211
         case ("tpss")
              rs6 =0.4535
              s18 =1.9435
              rs18=4.4752
         case ("hf")
              rs6 =0.3385
              s18 =0.9171
              rs18=2.8830
         case ("tpss0")
              rs6 =0.3768
              s18 =1.2576
              rs18=4.5865
         case ("pbe0")
              rs6 =0.4145
              s18 =1.2177
              rs18=4.8593
         case ("hse06")
              rs6 =0.383
              s18 =2.310
              rs18=5.685
         case ("revpbe38")
              rs6 =0.4309
              s18 =1.4760
              rs18=3.9446
         case ("pw6b95")
              rs6 =0.2076
              s18 =0.7257
              rs18=6.3750
         case ("b2-plyp")
              rs6 =0.3065
              s18 =0.9147
              rs18=5.0570
              s6=0.64d0
         case ("dsd-blyp")
              rs6 =0.0000
              s18 =0.2130
              rs18=6.0519
              s6=0.50d0
         case ("dsd-blyp-fc")
              rs6 =0.0009
              s18 =0.2112
              rs18=5.9807
              s6=0.50d0
         case ("bop")
              rs6 =0.4870
              s18 =3.2950
              rs18=3.5043
         case ("mpwlyp")
              rs6 =0.4831
              s18 =2.0077
              rs18=4.5323
         case ("o-lyp")
              rs6 =0.5299
              s18 =2.6205
              rs18=2.8065
         case ("pbesol")
              rs6 =0.4466
              s18 =2.9491
              rs18=6.1742
         case ("bpbe")
              rs6 =0.4567
              s18 =4.0728
              rs18=4.3908
         case ("opbe")
              rs6 =0.5512
              s18 =3.3816
              rs18=2.9444
         case ("ssb")
              rs6 =-0.0952
              s18 =-0.1744
              rs18=5.2170
         case ("revssb")
              rs6 =0.4720
              s18 =0.4389
              rs18=4.0986
         case ("otpss")
              rs6 =0.4634
              s18 =2.7495
              rs18=4.3153
         case ("b3pw91")
              rs6 =0.4312
              s18 =2.8524
              rs18=4.4693
         case ("bh-lyp")
              rs6 =0.2793
              s18 =1.0354
              rs18=4.9615
         case ("revpbe0")
              rs6 =0.4679
              s18 =1.7588
              rs18=3.7619
         case ("tpssh") !< See [1], Table 1.
              rs6 =0.4529
              s18 =2.2382
              rs18=4.6550
         case ("mpw1b95")
              rs6 =0.1955
              s18 =1.0508
              rs18=6.4177
         case ("pwb6k")
              rs6 =0.1805
              s18 =0.9383
              rs18=7.7627
         case ("b1b95")
              rs6 =0.2092
              s18 =1.4507
              rs18=5.5545
         case ("bmk")
              rs6 =0.1940
              s18 =2.0860
              rs18=5.9197
         case ("cam-b3lyp")
              rs6 =0.3708
              s18 =2.0674
              rs18=5.4743
         case ("lc-wpbe")
              rs6 =0.3919
              s18 =1.8541
              rs18=5.0897
         case ("b2gp-plyp")
              rs6 =0.0000
              s18 =0.2597
              rs18=6.3332
                s6=0.560
         case ("ptpss")
              rs6 =0.0000
              s18 =0.2804
              rs18=6.5745
                s6=0.750
         case ("pwpb95")
              rs6 =0.0000
              s18 =0.2904
              rs18=7.3141
                s6=0.820
c special HF/DFT with eBSSE correction
         case ("hf/mixed")
              rs6 =0.5607  
              s18 =3.9027  
              rs18=4.5622  
         case ("hf/sv")
              rs6 =0.4249  
              s18 =2.1849  
              rs18=4.2783   
         case ("hf/minis")
              rs6 =0.1702  
              s18 =0.9841   
              rs18=3.8506  
         case ("b3-lyp/6-31gd")
              rs6 =0.5014  
              s18 =4.0672   
              rs18=4.8409
         case ("hcth120")
              rs6=0.3563
              s18=1.0821
              rs18=4.3359
c DFTB3 (zeta=4.0)   
         case ("dftb3")
              rs6=0.7461
              s18=3.209 
              rs18=4.1906
         case DEFAULT
              call nwpwxc_stoprun( 'functional name unknown' )
      end select
      endif

c DFT-D3
      if(version.eq.3)then
      s6  =1.0d0
      alp =14.0d0
      rs18=1.0d0
c default def2-QZVP (almost basis set nwpwxc_limit)
      if(.not.TZ) then
      select case (func)
         case ("slater-dirac-exchange")
              rs6 =0.999
              s18 =-1.957
              rs18=0.697
         case ("b-lyp")
              rs6=1.094
              s18=1.682
         case ("b-p")
              rs6=1.139
              s18=1.683
         case ("b97-d")
              rs6=0.892
              s18=0.909
         case ("revpbe")
              rs6=0.923
              s18=1.010
         case ("pbe")
              rs6=1.217
              s18=0.722
         case ("pbesol")
              rs6=1.345
              s18=0.612
         case ("rpw86-pbe")
              rs6=1.224
              s18=0.901
         case ("rpbe")
              rs6=0.872
              s18=0.514
         case ("tpss")
              rs6=1.166
              s18=1.105
         case ("b3-lyp")
              rs6=1.261
              s18=1.703
         case ("pbe0")
              rs6=1.287
              s18=0.928

         case ("hse06")
              rs6=1.129
              s18=0.109
         case ("revpbe38")
              rs6=1.021
              s18=0.862
         case ("pw6b95")
              rs6=1.532
              s18=0.862
         case ("tpss0")
              rs6=1.252
              s18=1.242
         case ("b2-plyp")
              rs6=1.427
              s18=1.022
              s6=0.64
         case ("pwpb95")
              rs6=1.557
              s18=0.705
              s6=0.82
         case ("b2gp-plyp")
              rs6=1.586
              s18=0.760
              s6=0.56
         case ("ptpss")
              rs6=1.541
              s18=0.879
              s6=0.75
         case ("hf")
              rs6=1.158
              s18=1.746
         case ("mpwlyp")
              rs6=1.239
              s18=1.098
         case ("bpbe")
              rs6=1.087
              s18=2.033
         case ("bh-lyp")
              rs6=1.370
              s18=1.442
         case ("tpssh")
              rs6=1.223
              s18=1.219
         case ("pwb6k")
              rs6=1.660
              s18=0.550
         case ("b1b95")
              rs6=1.613
              s18=1.868
         case ("bop")
              rs6=0.929
              s18=1.975
         case ("o-lyp")
              rs6=0.806
              s18=1.764
         case ("o-pbe")
              rs6=0.837
              s18=2.055
         case ("ssb")
              rs6=1.215
              s18=0.663
         case ("revssb")
              rs6=1.221
              s18=0.560
         case ("otpss")
              rs6=1.128
              s18=1.494
         case ("b3pw91")
              rs6=1.176
              s18=1.775
         case ("revpbe0")
              rs6=0.949
              s18=0.792
         case ("pbe38")
              rs6=1.333
              s18=0.998
         case ("mpw1b95")
              rs6=1.605
              s18=1.118
         case ("mpwb1k")
              rs6=1.671
              s18=1.061
         case ("bmk")
              rs6=1.931
              s18=2.168
         case ("cam-b3lyp")
              rs6=1.378
              s18=1.217
         case ("lc-wpbe")
              rs6=1.355
              s18=1.279
         case ("m05")
              rs6=1.373
              s18=0.595
         case ("m052x")
              rs6=1.417
              s18=0.000
         case ("m06l")
              rs6=1.581
              s18=0.000
         case ("m06")
              rs6=1.325
              s18=0.000
         case ("m062x")
              rs6=1.619
              s18=0.000
         case ("m06hf")
              rs6=1.446
              s18=0.000
c DFTB3 (zeta=4.0)   
         case ("dftb3")
              rs6=1.235
              s18=0.673
         case ("hcth120")
              rs6=1.221
              s18=1.206
         case DEFAULT
              call nwpwxc_stoprun( 'functional name unknown' )
      end select
      else
c special TZVPP parameter
      select case (func)
         case ("b-lyp")
              rs6=1.243
              s18=2.022
         case ("b-p")
              rs6=1.221
              s18=1.838
         case ("b97-d")
              rs6=0.921
              s18=0.894
         case ("revpbe")
              rs6=0.953
              s18=0.989
         case ("pbe")
              rs6=1.277
              s18=0.777
         case ("tpss")
              rs6=1.213
              s18=1.176
         case ("b3-lyp")
              rs6=1.314
              s18=1.706
         case ("pbe0")
              rs6=1.328
              s18=0.926
         case ("pw6b95")
              rs6=1.562
              s18=0.821
         case ("tpss0")
              rs6=1.282
              s18=1.250
         case ("b2-plyp")
              rs6=1.551
              s18=1.109
              s6=0.5
         case DEFAULT
              call nwpwxc_stoprun( 'functional name unknown (TZ case)' )
      end select
      endif
      endif
c DFT-D2
      if(version.eq.2)then
      rs6=1.1d0
      s18=0.0d0
      alp=20.0d0
      select case (func)
         case ("b-lyp")
              s6=1.2  
         case ("b-p")
              s6=1.05  
         case ("b97-d")
              s6=1.25 
         case ("revpbe")
              s6=1.25 
         case ("pbe")
              s6=0.75 
         case ("tpss")
              s6=1.0  
         case ("b3-lyp")
              s6=1.05 
         case ("pbe0")
              s6=0.6   
         case ("pw6b95")
              s6=0.5   
         case ("tpss0")
              s6=0.85  
         case ("b2-plyp")
              s6=0.55 
         case ("b2gp-plyp")
              s6=0.4
         case ("dsd-blyp")
              s6=0.41
              alp=60.0d0
         case DEFAULT
              call nwpwxc_stoprun( 'functional name unknown' )
      end select

      endif

      end subroutine nwpwxc_setfuncpar

      subroutine nwpwxc_rdpar(dtmp,version,s6,s18,rs6,rs18,alp)
      implicit none
      real*8 s6,s18,rs6,rs18,alp
      integer version
      character*(*) dtmp
      character*80  ftmp
      logical ex
      real*8 xx(10)
      integer nn
c read parameter file
      call system('hostname > .tmpx')
      open(unit=43,file='.tmpx')
      read(43,'(a)')ftmp
      close(43,status='delete')
      write(dtmp,'(''~/.dftd3par.'',a)')trim(ftmp)
      inquire(file=dtmp,exist=ex)
      s6 =0
      s18=0
      rs6=0
      rs18=0
      alp =0
      if(ex)then
         open(unit=42,file=dtmp)
         read(42,'(a)',end=9)ftmp 
         call nwpwxc_readl(ftmp,xx,nn)
         if(nn.eq.6) then
                     s6  =xx(1)
                     rs6 =xx(2)
                     s18 =xx(3)
                     rs18=xx(4)
                     alp =xx(5)
            version=idint(xx(6))
         endif
  9      close(42)
      endif

      end subroutine nwpwxc_rdpar

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,rcov,
     .           rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr,
     .           e6,e8,e10,e12,e63)
      implicit none  

#include "global.fh"
#include "msgtypesf.h"

      integer n,iz(*),max_elem,maxc,version,mxc(max_elem) 
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
      real*8 rs6,rs8,rs10,alp6,alp8,alp10,rcov(max_elem)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8 e6, e8, e10, e12, e63        
      logical noabc
 
      integer iat,jat,kat
      real*8 r,r2,r6,r8,tmp,alp,dx,dy,dz,c6,c8,c10,ang,rav
      real*8 damp6,damp8,damp10,rr,thr,c9,r42,c12,r10,c14,rthr,cn_thr
      real*8 cn(n)                             
      real*8 r2ab(n*n),cc6ab(n*n),dmp(n*n),d2(3),t1,t2,t3,a1,a2
      integer*2 icomp(n*n)
      integer ij,ik,jk
      integer,external :: nwpwxc_lin

      integer taskid,np,pcount
      np     = GA_Nnodes()
      taskid = GA_Nodeid()
      pcount = 0


      e6 =0
      e8 =0
      e10=0
      e12=0
      e63=0

c     Becke-Johnson parameters
      a1=rs6
      a2=rs8

C DFT-D2
      if(version.eq.2)then

      do iat=1,n-1
         do jat=iat+1,n
            if (mod(pcount,np).eq.taskid) then
            dx=xyz(1,iat)-xyz(1,jat)
            dy=xyz(2,iat)-xyz(2,jat)
            dz=xyz(3,iat)-xyz(3,jat)
            r2=dx*dx+dy*dy+dz*dz
c           if(r2.gt.rthr) cycle
            r=sqrt(r2)
            c6=c6ab(iz(jat),iz(iat),1,1,1)
            damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
            r6=r2**3      
            e6 =e6+c6*damp6/r6
            end if
            pcount = pcount + 1
         enddo
      enddo
      if (np.gt.1) call GA_DGOP(9+MSGDBL,e6,1,'+')

      else
C DFT-D3
      call nwpwxc_ncoord(n,rcov,iz,xyz,cn,cn_thr)

      icomp=0
      do iat=1,n-1
         do jat=iat+1,n

         if (mod(pcount,np).eq.taskid) then
            dx=xyz(1,iat)-xyz(1,jat)
            dy=xyz(2,iat)-xyz(2,jat)
            dz=xyz(3,iat)-xyz(3,jat)
            r2=dx*dx+dy*dy+dz*dz
CTHR        
            if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r
c damping
            tmp=rs6*rr   
            damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
            tmp=rs8*rr     
            damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
c get C6
            call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                    cn(iat),cn(jat),c6)
            if(.not.noabc)then
            ij=nwpwxc_lin(jat,iat)
            icomp(ij)=1
c store C6 for C9, calc as sqrt
            cc6ab(ij)=sqrt(c6)
c store R^2 for abc
            r2ab(ij)=r2
c store for abc damping
            dmp(ij)=(1./rr)**(1./3.) 
            endif

            r6=r2**3      
            r8=r6*r2
c r2r4 stored in main as sqrt
            c8 =3.0d0*c6*r2r4(iz(iat))*r2r4(iz(jat))

C DFT-D3 zero-damp
            if(version.eq.3)then
               e6=e6+c6*damp6/r6
               e8=e8+c8*damp8/r8
            endif
C DFT-D3(BJ)          
            if(version.eq.4)then
c use BJ radius
               tmp=sqrt(c8/c6)              
               e6=e6+  c6/(r6+(a1*tmp+a2)**6)
               e8=e8+  c8/(r8+(a1*tmp+a2)**8)
c              c10=(49.0d0/40.0d0)*c8**2/c6
c             e10=e10+c10/(r8*r2+(a1*tmp+a2)**10)
            endif

         end if 
         pcount = pcount + 1
         enddo
      enddo
      if (np.gt.1) call GA_DGOP(9+MSGDBL,e6,1,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,e8,1,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,e10,1,'+')

      if(noabc)return

C compute non-additive third-order energy using averaged C6
      do iat=1,n-1
         do jat=iat+1,n
         if (mod(pcount,np).eq.taskid) then

            ij=nwpwxc_lin(jat,iat)
            if(icomp(ij).eq.1)then
            do kat=jat+1,n
            ik=nwpwxc_lin(kat,iat)
            jk=nwpwxc_lin(kat,jat)
            if(icomp(ik).eq.0.or.icomp(jk).eq.0) cycle
c damping func product
c           tmp=dmp(ik)*dmp(jk)*dmp(ij)
            rav=(4./3.)/(dmp(ik)*dmp(jk)*dmp(ij))
            tmp=1.d0/( 1.d0+6.d0*rav**alp8 )
c triple C6 coefficient (stored as sqrt)
            c9=cc6ab(ij)*cc6ab(ik)*cc6ab(jk)
c           write(*,*) 'C9', c9
c angular terms
c d is r^2
            d2(1)=r2ab(ij)
            d2(2)=r2ab(jk)
            d2(3)=r2ab(ik)
            t1 = (d2(1)+d2(2)-d2(3))/sqrt(d2(1)*d2(2))
            t2 = (d2(1)+d2(3)-d2(2))/sqrt(d2(1)*d2(3))
            t3 = (d2(3)+d2(2)-d2(1))/sqrt(d2(2)*d2(3))
            ang=0.375d0*t1*t2*t3+1.0d0

c C9 has negative sign
            e63=e63-tmp*c9*ang/(d2(1)*d2(2)*d2(3))**1.50d0

            enddo
            endif
         end if
         pcount = pcount + 1
         enddo
      enddo
      if (np.gt.1) call GA_DGOP(9+MSGDBL,e63,1,'+')

      endif

      end subroutine nwpwxc_edisp

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C analyse all pairs
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_adisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >                r2r4,r0ab,rcov,
     .                rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal,
     .                autoang,rthr,cn_thr,s6,s18,etot)
      implicit none  
      integer n,iz(*),max_elem,maxc,version,mxc(max_elem) 
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*),s6
      real*8 rs6,rs8,rs10,alp6,alp8,alp10,autokcal,etot,s18,autoang
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3),rcov(max_elem)
 
      integer iat,jat,i,j,k,nbin
      real*8 R0,r,r2,r6,r8,tmp,alp,dx,dy,dz,c6,c8,c10
      real*8 damp6,damp8,damp10,r42,rr,check,rthr,cn_thr,rvdw
      real*8 cn(n),i6,e6,e8,e10,edisp                   
      real*8 dist(0:15),li(0:15,2)
      real*8 xx(500),eg(10000)
      integer grplist(500,20)
      integer grpn(20),at(n)
      integer ngrp,dash
      integer nwpwxc_lin, iiii, jjjj, iii, jjj, ii, jj, ni, nj 
      integer iout(500)
      logical ex
      character*80 atmp
 
      real*8,dimension(:,:), allocatable :: ed
      allocate(ed(n,n))


c distance bins
      li(0,1)=0   
      li(0,2)=1.5 
      li(1,1)=1.5
      li(1,2)=2
      li(2,1)=2
      li(2,2)=2.3333333333
      li(3,1)=2.3333333333
      li(3,2)=2.6666666666
      li(4,1)=2.6666666666
      li(4,2)=3.0
      li(5,1)=3.0          
      li(5,2)=3.3333333333
      li(6,1)=3.3333333333
      li(6,2)=3.6666666666
      li(7,1)=3.6666666666
      li(7,2)=4.0
      li(8,1)=4.0
      li(8,2)=4.5
      li(9,1)=4.5
      li(9,2)=5.0
      li(10,1)=5.0
      li(10,2)=5.5
      li(11,1)=5.5
      li(11,2)=6.0
      li(12,1)=6.0
      li(12,2)=7.0           
      li(13,1)=7.0           
      li(13,2)=8.0           
      li(14,1)=8.0           
      li(14,2)=9.0           
      li(15,1)=9.0           
      li(15,2)=10.0          
      nbin=15

      call nwpwxc_ncoord(n,rcov,iz,xyz,cn,cn_thr)

      write(*,*)
      write(*,*)'analysis of pair-wise terms (in kcal/mol)'
      write(*,'(''pair'',2x,''atoms'',9x,''C6'',14x,''C8'',12x,
     .''E6'',7x,''E8'',7x,''Edisp'')')
      e8=0
      ed=0
      dist=0
      check=0
      do iat=1,n-1
         do jat=iat+1,n

            dx=xyz(1,iat)-xyz(1,jat)
            dy=xyz(2,iat)-xyz(2,jat)
            dz=xyz(3,iat)-xyz(3,jat)
            r2=(dx*dx+dy*dy+dz*dz)
CTHR
            if(r2.gt.rthr) cycle
            r =sqrt(r2)
            R0=r0ab(iz(jat),iz(iat))
            rr=R0/r
            r6=r2**3       

            tmp=rs6*rr   
            damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
            tmp=rs8*rr     
            damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )

            if (version.eq.2)then
              c6=c6ab(iz(jat),iz(iat),1,1,1)
              damp6=1.d0/(1.d0+exp(-alp6*(r/(rs6*R0)-1.0d0)))
              e6 =s6*autokcal*c6*damp6/r6
              e8=0.0
            else
              call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                      cn(iat),cn(jat),c6)
            endif

            if(version.eq.3)then
            e6 =s6*autokcal*c6*damp6/r6
            r8 =r6*r2
            r42=r2r4(iz(iat))*r2r4(iz(jat))
            c8 =3.0d0*c6*r42
            e8 =s18*autokcal*c8*damp8/r8
            endif

            if(version.eq.4)then
            r42=r2r4(iz(iat))*r2r4(iz(jat))
            c8 =3.0d0*c6*r42
c use BJ radius
            R0=sqrt(c8/c6)              
            rvdw=rs6*R0+rs8
            e6 =s6*autokcal*c6/(r6+rvdw**6)
            r8 =r6*r2
            e8 =s18*autokcal*c8/(r8+rvdw**8)
            endif

            edisp=-(e6+e8)
            ed(iat,jat)=edisp
            ed(jat,iat)=edisp

           write(*,'(2i4,2x,2i3,2D16.6,2F9.4,F10.5)')
     .     iat,jat,iz(iat),iz(jat),c6,c8,
     .    -e6,-e8,edisp

            check=check+edisp
            rr=r*autoang
            do i=0,nbin
               if(rr.gt.li(i,1).and.rr.le.li(i,2)) dist(i)=dist(i)+edisp
            enddo
         enddo
      enddo

      write(*,'(/''distance range (Angstroem) analysis'')')
      write(*,'( ''writing histogram data to <histo.dat>'')')
      open(unit=11,file='histo.dat')
      do i=0,nbin
         write(*,'(''R(low,high), Edisp, %tot :'',2f4.1,F12.5,F8.2)')
     .   li(i,1),li(i,2),dist(i),100.*dist(i)/etot
         write(11,*)(li(i,1)+li(i,2))*0.5,dist(i)
      enddo
      close(11)

      write(*,*) 'checksum (Edisp) ',check
      if(abs(check-etot).gt.1.d-3)stop'something is weired in adisp'

      inquire(file='fragment',exist=ex)
      if(.not.ex) return
      write(*,'(/''fragment based analysis'')')
      write(*,'( ''reading file <fragment> ...'')')
      open(unit=55,file='fragment')
      i=0
      at=0
 111  read(55,'(a)',end=222) atmp
      call nwpwxc_readfrag(atmp,iout,j)
      if(j.gt.0)then
         i=i+1
         grpn(i)=j
         do k=1,j
            grplist(k,i)=iout(k)      
            at(grplist(k,i))=at(grplist(k,i))+1
         enddo
      endif
      goto 111
 222  continue
      ngrp=i  
      k=0
      do i=1,n
         if(at(i).gt.1) stop 'something is weird in file <fragment>'
         if(at(i).eq.0)then
            k=k+1
            grplist(k,ngrp+1)=i
         endif
      enddo
      if(k.gt.0) then
         ngrp=ngrp+1
         grpn(ngrp)=k
      endif
c Implemented display of atom ranges instead of whole list of atoms
      write(*,*)'group #        atoms '
      dash=0
      do i=1,ngrp
       write(*,'(i4,3x,i4)',advance='no')i,grplist(1,i)
       do j=2,grpn(i)
        if(grplist(j,i).eq.(grplist(j-1,i)+1)) then
         if(dash.eq.0)then
          write(*,'(A1)',advance='no')'-'
          dash=1
         endif
        else
         if(dash.eq.1)then
          write(*,'(i4)',advance='no') grplist(j-1,i)
          dash=0
         endif
         write(*,'(i4)',advance='no') grplist(j,i)
        endif
       enddo 
       if(dash.eq.1)then
        write(*,'(i4)',advance='no') grplist(j-1,i)
        dash=0
       endif
      write(*,*)''
      enddo

c old display list code
c      write(*,*)'group #        atoms '
c      do i=1,ngrp      
c         write(*,'(i4,3x,100i3)')i,(grplist(j,i),j=1,grpn(i))
c      enddo

      eg=0
      iii=0
      do i=1,ngrp
         ni=grpn(i)
         iii=iii+1
         jjj=0
         do j=1,ngrp
            nj=grpn(j)
            jjj=jjj+1
            do ii=1,ni
               iiii=grplist(ii,i)
               do jj=1,nj
                  jjjj=grplist(jj,j)
                  if(jjjj.lt.iiii)cycle
                  eg(nwpwxc_lin(iii,jjj))
     >             =eg(nwpwxc_lin(iii,jjj))+ed(iiii,jjjj)
               enddo
            enddo
         enddo
      enddo

c     call nwpwxc_prmat(6,eg,ngrp,0,'intra- + inter-group dispersion energies')
      write(*,*)' group i      j     Edisp'
      k=0
      check=0
      do i=1,ngrp
      do j=1,i    
      k=k+1
      check=check+eg(k) 
      write(*,'(5x,i4,'' --'',i4,F8.2)')i,j,eg(k)
      enddo
      enddo
      write(*,*) 'checksum (Edisp) ',check

      end subroutine nwpwxc_adisp

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute gradient
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine nwpwxc_gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
c     .                 s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,rthr,
c     .                 num,version,echo,g,disp,gnorm,cn_thr,fix)
c
      subroutine nwpwxc_gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >                 r2r4,r0ab,rcov,
     .                 s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,rthr,
     .                 num,version,echo,g,disp,gnorm,cn_thr)
      implicit none  
      !include  'param'

#include "global.fh"
#include "msgtypesf.h"


      real*8 k1,k2,k3

c global ad hoc parameters
      parameter (k1=16.0)
      parameter (k2=4./3.) 

c reasonable choices are between 3 and 5
c this gives smoth curves with maxima around the integer values
c k3=3 give for CN=0 a slightly smaller value than computed
c for the free atom. This also yields to larger CN for atoms
c in larger molecules but with the same chem. environment
c which is physically not right
c values >5 might lead to bumps in the potential
      parameter (k3=-4.) 


      integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8 g(3,*),s6,s18,rcov(max_elem)
      real*8 rs6,rs8,rs10,alp10,alp8,alp6,a1,a2,r2ik        
      !logical noabc,num,echo,fix(n)
      logical noabc,num,echo
 
      integer iat,jat,i,j,kat
      real*8 R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc  
      real*8 dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
      real*8 damp6,damp8,damp10,e6,e8,e10,e12,gnorm,tmp1
      real*8 s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
      real*8 cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr,c8,cn_thr
      real*8 rthr3

      real*8 rij(3),rik(3),rjk(3),r7,r9
      real*8 rik_dist,rjk_dist
      real*8 drij(n*(n+1)/2)  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
      real*8 drik,drjk
      real*8 rcovij
      real*8 dc6,c6chk !d(C6ij)/d(r_ij)
      real*8 expterm,dcni
      real*8 dcn(n*(n+1)/2)    !dCN(iat)/d(r_ij) is equal to
                               !dCN(jat)/d(r_ij)     
      real*8 dc6_rest(n*(n+1)/2) ! saves (1/r^6*f_dmp + 3*r4r2/r^8*f_dmp) for kat loop
      integer,external :: nwpwxc_lin
      integer  linij,linik,linjk
      real*8 vec(3),vec2(3)
      real*8 dc6ij(n,n)       !dC6(iat,jat)/dCN(iat) in dc6ij(i,j)
                              !dC6(iat,jat)/cCN(jat) in dc6ij(j,i)
      logical skip(n*(n+1)/2)                        

      integer taskid,np,pcount

      np     = GA_Nnodes()
      taskid = GA_Nodeid()
      pcount = 0


      dc6ij=0.0d0



      

c this is the crucial threshold to reduce the N^3 to an
c effective N^2. 

c      write(*,*)'rthr=',rthr,'rthr2=',rthr2,'rthr3=',rthr3

      if(echo)write(*,*) 
c 2222222222222222222222222222222222222222222222222222222222222222222222222
      if(version.eq.2)then
      if(echo)write(*,*) 'doing analytical gradient O(N^2) ...'
      disp=0
      do iat=1,n-1
         do jat=iat+1,n
         if (mod(pcount,np).eq.taskid) then
            R0=r0ab(iz(jat),iz(iat))*rs6
            dx=(xyz(1,iat)-xyz(1,jat))
            dy=(xyz(2,iat)-xyz(2,jat))
            dz=(xyz(3,iat)-xyz(3,jat))
            r2  =dx*dx+dy*dy+dz*dz             
c           if(r2.gt.rthr) cycle
            r235=r2**3.5                       
            r   =sqrt(r2)
            damp6=exp(-alp6*(r/R0-1.0d0))
            damp1=1.+damp6           
            c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
            tmp1=damp6/(damp1*damp1*r235*R0)
            tmp2=6./(damp1*r*r235)
            gx1=alp6* dx*tmp1-tmp2*dx
            gx2=alp6*(-dx)*tmp1+tmp2*dx
            gy1=alp6* dy*tmp1-tmp2*dy
            gy2=alp6*(-dy)*tmp1+tmp2*dy
            gz1=alp6* dz*tmp1-tmp2*dz
            gz2=alp6*(-dz)*tmp1+tmp2*dz
            g(1,iat)=g(1,iat)-gx1*c6
            g(2,iat)=g(2,iat)-gy1*c6
            g(3,iat)=g(3,iat)-gz1*c6
            g(1,jat)=g(1,jat)-gx2*c6  
            g(2,jat)=g(2,jat)-gy2*c6      
            g(3,jat)=g(3,jat)-gz2*c6      
            disp=disp+c6*(1./damp1)/r2**3
         end if
         pcount = pcount + 1
         enddo
      enddo
      if (np.gt.1) call GA_DGOP(9+MSGDBL,disp,1,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,g,3*n,'+')
      disp=-disp
      goto 999
      endif

cNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
      if(num) then
      if(echo)write(*,*) 'doing numerical gradient O(N^3) ...'

      call nwpwxc_edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr,
     .     e6,e8,e10,e12,e6abc)
      disp=-s6*e6-s18*e8-s6*e6abc

      step=2.d-5     

      do i=1,n
      do j=1,3
      xyz(j,i)=xyz(j,i)+step        
      call nwpwxc_edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr,
     .     e6,e8,e10,e12,e6abc)
      dispr=-s6*e6-s18*e8-s6*e6abc
      xyz(j,i)=xyz(j,i)-2*step      
      call nwpwxc_edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr,
     .     e6,e8,e10,e12,e6abc)
      displ=-s6*e6-s18*e8-s6*e6abc
      g(j,i)=0.5*(dispr-displ)/step  
      xyz(j,i)=xyz(j,i)+step        
      enddo
      enddo

      else

      if(echo)write(*,*) 'doing analytical gradient O(N^3) ...'
c precompute for analytical part

c 333333333333333333333333333333333333333333333333333333333333333333333333333
c standard correction
      if (version.eq.3) then
      call nwpwxc_ncoord(n,rcov,iz,xyz,cn,cn_thr)
      s8 =s18
      s10=s18

      disp=0

      dc6ij = 0.0d0
      drij=0.0d0
      dc6_rest=0.0d0
      dcn=0.0d0
      kat=0
      skip=.true.


*     *** compute dc6ij,dc6_rest and dcn  and non dcn part of drij****
      do iat=1,n
         do jat=1,iat-1
             if (mod(pcount,np).eq.taskid) then
             linij=nwpwxc_lin(iat,jat)
             rij=xyz(:,jat)-xyz(:,iat)
             r2=sum(rij*rij)
             if (r2.le.rthr) then

                r0=r0ab(iz(jat),iz(iat))
                r42=r2r4(iz(iat))*r2r4(iz(jat))
                rcovij=rcov(iz(iat))+rcov(iz(jat))

                call nwpwxc_get_dC6_dCNij(maxc,max_elem,
     >            c6ab,mxc(iz(iat)),
     >            mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,
     >            c6,dc6ij(iat,jat),dc6ij(jat,iat))

                r=dsqrt(r2)
                r6=r2*r2*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                t6 = (r/(rs6*R0))**(-alp6)
                damp6 =1.d0/( 1.d0+6.d0*t6 )
                t8 = (r/(rs8*R0))**(-alp8)
                damp8 =1.d0/( 1.d0+6.d0*t8 )

                drij(linij)=drij(linij)-s6*(6.0/(r7)*C6*damp6)  ! d(r^(-6))/d(r_ij)
     >            -s8*(24.0/(r9)*C6*r42*damp8)

                drij(linij)=drij(linij)
     >            +s6*C6/r7*6.d0*alp6*t6*damp6*damp6     !d(f_dmp)/d(r_ij)
     >            +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8

                dc6_rest(linij)=s6/r6*damp6+3.d0*s8*r42/r8*damp8
                disp=disp-dc6_rest(linij)*c6 

                if (r2.lt.cn_thr) then
                   expterm=exp(-k1*(rcovij/r-1.d0))
                   dcn(linij)=-k1*rcovij*expterm/
     >                 (r*r*(expterm+1.d0)*(expterm+1.d0))
                end if

             end if
             end if
             pcount = pcount + 1
         end do
      end do
      if (np.gt.1) call GA_DGOP(9+MSGDBL,disp,1,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,dc6ij,n*n,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,dc6_rest,n*(n+1)/2,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,dcn,     n*(n+1)/2,'+')



*     **** compute dcn part of drij ****
      do iat=1,n
        do jat=1,iat-1

          if (mod(pcount,np).eq.taskid) then
          linij=nwpwxc_lin(iat,jat)

          dc6=(dc6ij(iat,jat)+dc6ij(jat,iat))*dcn(linij)
          drij(linij)=drij(linij)+dc6_rest(linij)*dc6 
!
!  Basically all term that depend on the coordinates of
!  3 atoms
!  This is the reason, why the gradient scales N^3            
!            
          do kat=1,jat-1
            linik=nwpwxc_lin(iat,kat)
            linjk=nwpwxc_lin(jat,kat)

            drij(linij) = drij(linij)
     >                  + dc6_rest(linik)*dc6ij(iat,kat)*dcn(linij)
     >                  + dc6_rest(linjk)*dc6ij(jat,kat)*dcn(linij)

            drij(linjk) = drij(linjk)
     >                  + dc6_rest(linik)*dc6ij(kat,iat)*dcn(linjk)
     >                  + dc6_rest(linij)*dc6ij(jat,iat)*dcn(linjk)
             
            drij(linik) = drij(linik)
     >                  + dc6_rest(linjk)*dc6ij(kat,jat)*dcn(linik) 
     >                  + dc6_rest(linij)*dc6ij(iat,jat)*dcn(linik)

          end do!kat

        end if
        pcount = pcount + 1
        enddo !jat
      enddo !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,drij,n*(n+1)/2,'+')

! After calculating all derivatives dE/dr_ij w.r.t. distances,
! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i       
      do iat=2,n
        do jat=1,iat-1
        if (mod(pcount,np).eq.taskid) then
          rij=xyz(:,jat)-xyz(:,iat)
          r=sqrt(sum(rij*rij))
          g(:,iat)=g(:,iat)+drij(nwpwxc_lin(iat,jat))*rij/r
          g(:,jat)=g(:,jat)-drij(nwpwxc_lin(iat,jat))*rij/r


        end if
        pcount = pcount + 1
        enddo !iat
      enddo !jat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,g,3*n,'+')

!      stop('Bis hier und nicht weiter!')
      endif !version 3

c BJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJ 
c Becke-Johnson finite damping 
      if (version.eq.4) then 
      a1 =rs6
      a2 =rs8
      s8 =s18

      disp=0
      call nwpwxc_ncoord(n,rcov,iz,xyz,cn,cn_thr)

      dc6ij = 0.0d0
      drij=0.0d0
      dc6_rest=0.0d0
      dcn=0.0d0
      kat=0
      skip=.true.

*     *** compute dc6ij,dc6_rest and dcn  and non dcn part of drij****
      do iat=1,n
         do jat=1,iat-1
            rij=xyz(:,jat)-xyz(:,iat)
            r2=sum(rij*rij)
            if (r2.le.rthr) then
               linij=nwpwxc_lin(iat,jat)
               r0=r0ab(iz(jat),iz(iat))
               r42=r2r4(iz(iat))*r2r4(iz(jat))
               rcovij=rcov(iz(iat))+rcov(iz(jat))

               call nwpwxc_get_dC6_dCNij(maxc,max_elem,
     >           c6ab,mxc(iz(iat)),
     >           mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,
     >           c6,dc6ij(iat,jat),dc6ij(jat,iat))

               r=dsqrt(r2)
               r4=r2*r2
               r6=r4*r2
               r7=r6*r
               r8=r6*r2
               r9=r8*r

               R0=a1*sqrt(3.0d0*r42)+a2  !c use BJ radius
               t6=(r6+R0**6)
               t8=(r8+R0**8)

               drij(linij)=drij(linij)
     >           -s6*C6*6.0d0*r4*r/(t6*t6)
     >           -s8*C6*24.0d0*r42*r6*r/(t8*t8)

               dc6_rest(linij)=s6/t6+3.d0*s8*r42/t8
               disp=disp-dc6_rest(linij)*c6  ! calculate E_disp for sanity check

               if (r2.lt.cn_thr) then
                  expterm=exp(-k1*(rcovij/r-1.d0))
                  dcn(linij)=-k1*rcovij*expterm/
     >                 (r*r*(expterm+1.d0)*(expterm+1.d0))
               end if
            end if
         end do
      end do
      if (np.gt.1) call GA_DGOP(9+MSGDBL,disp,1,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,dc6ij,n*n,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,dc6_rest,n*(n+1)/2,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,dcn,     n*(n+1)/2,'+')


*     **** compute dcn part of drij ****
      do iat=1,n
         do jat=1,iat-1
            if (mod(pcount,np).eq.taskid) then
            linij=nwpwxc_lin(iat,jat)

            dc6=(dc6ij(iat,jat)+dc6ij(jat,iat))*dcn(linij)
            drij(linij)=drij(linij) + dc6_rest(linij)*dc6            !d(C6(ij))/d(r_ij)

!  Basically all term that depend on the coordinates of
!  3 atoms
!  This is the reason, why the gradient scales N^3            
!            
            do kat=1,jat-1
               linik=nwpwxc_lin(iat,kat)
               linjk=nwpwxc_lin(jat,kat)

               drij(linij) = drij(linij)
     >                     + dc6_rest(linik)*dc6ij(iat,kat)*dcn(linij)
     >                     + dc6_rest(linjk)*dc6ij(jat,kat)*dcn(linij)

               drij(linjk) = drij(linjk)
     >                     + dc6_rest(linik)*dc6ij(kat,iat)*dcn(linjk)
     >                     + dc6_rest(linij)*dc6ij(jat,iat)*dcn(linjk)

               drij(linik) = drij(linik)
     >                     + dc6_rest(linjk)*dc6ij(kat,jat)*dcn(linik) 
     >                     + dc6_rest(linij)*dc6ij(iat,jat)*dcn(linik)

            end do!kat
            end if
            pcount = pcount + 1
         end do !jat
      end do !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,drij,n*(n+1)/2,'+')

! After calculating all derivatives dE/dr_ij w.r.t. distances,
! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i       
      do iat=2,n
        do jat=1,iat-1
          if (mod(pcount,np).eq.taskid) then
          rij=xyz(:,jat)-xyz(:,iat)
          r=sqrt(sum(rij*rij))
          g(:,iat)=g(:,iat)+drij(nwpwxc_lin(iat,jat))*rij/r
          g(:,jat)=g(:,jat)-drij(nwpwxc_lin(iat,jat))*rij/r

        end if
        pcount = pcount + 1

        enddo !iat
      enddo !jat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,g,3*n,'+')


      endif !version=4 (BJ)


      endif

 999  continue
      gnorm=sum(abs(g(1:3,1:n)))
      if(echo)then
      write(*,*)
      write(*,*)'|G|=',gnorm
      endif

c      do i=1,n                 !fixed atoms have no gradient
c       if(fix(i))g(:,i)=0.0
c      enddo

       

      end subroutine nwpwxc_gdisp


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      The   N E W   gradC6 routine    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      subroutine nwpwxc_get_dC6_dCNij(maxc,max_elem,c6ab,mxci,mxcj,
     >           cni,cnj,
     .           izi,izj,iat,jat,c6check,dc6i,dc6j)

      IMPLICIT NONE
      !include  'param'

      real*8 k1,k2,k3

c global ad hoc parameters
      parameter (k1=16.0)
      parameter (k2=4./3.) 

c reasonable choices are between 3 and 5
c this gives smoth curves with maxima around the integer values
c k3=3 give for CN=0 a slightly smaller value than computed
c for the free atom. This also yields to larger CN for atoms
c in larger molecules but with the same chem. environment
c which is physically not right
c values >5 might lead to bumps in the potential
      parameter (k3=-4.) 


      integer maxc,max_elem
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      integer mxci,mxcj   !mxc(iz(iat))
      real*8 cni,cnj
      integer iat,jat,izi,izj
      real*8  dc6i,dc6j,c6check


      integer i,j,a,b
      real*8 zaehler,nenner,dzaehler_i,dnenner_i,dzaehler_j,dnenner_j
      real*8 expterm,cn_refi,cn_refj,c6ref,r
      real*8 c6mem,r_save



      c6mem=-1.d99
      r_save=9999.0
      zaehler=0.0d0
      nenner=0.0d0

      dzaehler_i=0.d0
      dnenner_i=0.d0
      dzaehler_j=0.d0
      dnenner_j=0.d0


      DO a=1,mxci
        DO b=1,mxcj
          c6ref=c6ab(izi,izj,a,b,1)
          if (c6ref.gt.0) then
!            c6mem=c6ref
            cn_refi=c6ab(izi,izj,a,b,2)
            cn_refj=c6ab(izi,izj,a,b,3)
            r=(cn_refi-cni)*(cn_refi-cni)+(cn_refj-cnj)*(cn_refj-cnj)
            if (r.lt.r_save) then
               r_save=r
               c6mem=c6ref
            endif
            expterm=exp(k3*r)
            zaehler=zaehler+c6ref*expterm
            nenner=nenner+expterm
            dzaehler_i=dzaehler_i+c6ref*expterm*
     .             2.d0*k3*(cni-cn_refi)
            dnenner_i=dnenner_i+expterm*
     .             2.d0*k3*(cni-cn_refi)

            dzaehler_j=dzaehler_j+c6ref*expterm*
     .             2.d0*k3*(cnj-cn_refj)
            dnenner_j=dnenner_j+expterm*
     .             2.d0*k3*(cnj-cn_refj)
          endif
        ENDDO !b
      ENDDO !a

      if (nenner.gt.1.0d-99) then
        c6check=zaehler/nenner
        dc6i=((dzaehler_i*nenner)-(dnenner_i*zaehler))
     .    /(nenner*nenner)
        dc6j=((dzaehler_j*nenner)-(dnenner_j*zaehler))
     .    /(nenner*nenner)
      else
        c6check=c6mem
        dc6i=0.0d0
        dc6j=0.0d0
      endif
      end subroutine nwpwxc_get_dC6_dCNij



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C interpolate c6  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
      implicit none
      integer maxc,max_elem
      integer iat,jat,i,j,mxc(max_elem)
      real*8  nci,ncj,c6,c6mem
      real*8  c6ab(max_elem,max_elem,maxc,maxc,3)
c the exponential is sensitive to numerics
c when nci or ncj is much larger than cn1/cn2
      real*8  cn1,cn2,r,rsum,csum,tmp,tmp1   
      real*8  r_save
      !include 'param'

      real*8 k1,k2,k3

c global ad hoc parameters
      parameter (k1=16.0)
      parameter (k2=4./3.) 

c reasonable choices are between 3 and 5
c this gives smoth curves with maxima around the integer values
c k3=3 give for CN=0 a slightly smaller value than computed
c for the free atom. This also yields to larger CN for atoms
c in larger molecules but with the same chem. environment
c which is physically not right
c values >5 might lead to bumps in the potential
      parameter (k3=-4.) 


      c6mem=-1.d+99
      rsum=0.0
      csum=0.0
      c6  =0.0
      r_save=1.0d99
      do i=1,mxc(iat)
      do j=1,mxc(jat)
         c6=c6ab(iat,jat,i,j,1)
         if(c6.gt.0)then
!            c6mem=c6
            cn1=c6ab(iat,jat,i,j,2)
            cn2=c6ab(iat,jat,i,j,3)
c distance
            r=(cn1-nci)**2+(cn2-ncj)**2
            if (r.lt.r_save) then
               r_save=r
               c6mem=c6
            endif
            tmp1=exp(k3*r)
            rsum=rsum+tmp1     
            csum=csum+tmp1*c6
         endif
      enddo
      enddo

      if(rsum.gt.1.0d-99)then
         c6=csum/rsum
      else
         c6=c6mem
      endif

      end subroutine nwpwxc_getc6

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute coordination numbers by adding an inverse damping function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_ncoord(natoms,rcov,iz,xyz,cn,cn_thr)
      implicit none  
      !include 'param'

#include "global.fh"
#include "msgtypesf.h"

      real*8 k1,k2,k3

c global ad hoc parameters
      parameter (k1=16.0)
      parameter (k2=4./3.) 

c reasonable choices are between 3 and 5
c this gives smoth curves with maxima around the integer values
c k3=3 give for CN=0 a slightly smaller value than computed
c for the free atom. This also yields to larger CN for atoms
c in larger molecules but with the same chem. environment
c which is physically not right
c values >5 might lead to bumps in the potential
      parameter (k3=-4.) 

      integer iz(*),natoms,i,max_elem
      real*8 xyz(3,*),cn(*),rcov(94),input
      real*8 cn_thr

      integer iat    
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2
      integer taskid,np,pcount

      call dcopy(natoms,0.0d0,0,cn,1)
      np     = GA_Nnodes()
      taskid = GA_Nodeid()
      pcount = 0


      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
      if (mod(pcount,np).eq.taskid) then
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz 
            if (r2.gt.cn_thr) cycle 
            r=sqrt(r2)                  
c covalent distance in Bohr
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            xn=xn+damp
         endif
      end if
      pcount = pcount + 1
      enddo
c     if (iz(i).eq.19) then
c        write(*,*) "Input CN of Kalium"
c        read(*,*),input
c         cn(i)=input
c      else
         cn(i)=xn
c      endif  
      enddo

      if (np.gt.1) call GA_DGOP(9+MSGDBL,cn,natoms,'+')

      end subroutine nwpwxc_ncoord

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C load C6 coefficients from file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_loadc6(fname,maxc,max_elem,c6ab,maxci)
      implicit none
      integer maxc,max_elem,maxci(max_elem)
      real*8  c6ab(max_elem,max_elem,maxc,maxc,3)
      character*(*) fname
      character*1  atmp 
      character*80 btmp 

      real*8  x,y,f,cn1,cn2,cmax,xx(10)
      integer iat,jat,i,n,l,j,k,il,iadr,jadr,nn

      c6ab=-1
      maxci=0

c read file
      open(unit=1,file=fname)
      read(1,'(a)')btmp
 10   read(1,*,end=11) y,iat,jat,cn1,cn2  
      call nwpwxc_limit(iat,jat,iadr,jadr)
      maxci(iat)=max(maxci(iat),iadr)
      maxci(jat)=max(maxci(jat),jadr)
      c6ab(iat,jat,iadr,jadr,1)=y  
      c6ab(iat,jat,iadr,jadr,2)=cn1
      c6ab(iat,jat,iadr,jadr,3)=cn2

      c6ab(jat,iat,jadr,iadr,1)=y  
      c6ab(jat,iat,jadr,iadr,2)=cn2
      c6ab(jat,iat,jadr,iadr,3)=cn1
c     endif
      goto 10
 11   continue
      close(1)

      end subroutine nwpwxc_loadc6

      integer function nwpwxc_lin(i1,i2)
      integer i1,i2,idum1,idum2
      idum1=max(i1,i2)
      idum2=min(i1,i2)
      nwpwxc_lin=idum2+idum1*(idum1-1)/2
      return
      end function nwpwxc_lin

      subroutine nwpwxc_limit(iat,jat,iadr,jadr)
      implicit none
      integer iat,jat,iadr,jadr,i
      iadr=1
      jadr=1
      i=100
 10   if(iat.gt.100) then
         iat=iat-100
         iadr=iadr+1
         goto 10
      endif

      i=100
 20   if(jat.gt.100) then
         jat=jat-100
         jadr=jadr+1
         goto 20
      endif

      end subroutine nwpwxc_limit

c $Id$
