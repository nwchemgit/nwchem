      program libcheck
* $Id: libcheck.F,v 1.1 1997-12-10 01:00:56 d3e129 Exp $
      implicit none
#include "inp.fh"
      character*5 btype
      character*40 name
      character*100 field
      integer pmax, cmax
      parameter (pmax=30,cmax=30)
      double precision ee(pmax),cc(pmax,cmax)
      integer rex(pmax)
      logical oecp
      integer ip, ig, nprim, ng, nbasis, nwrong
c
      open(unit=7,file='library',status='old',
     &    access='sequential',form='formatted',err=90000)
c
      call inp_init(7,6)
      nbasis = 0
      nwrong = 0
00001 if (.not.inp_read()) then
        goto 70000
      else
*        call inp_outrec()
      endif
      call inp_set_field(0) 
      if (.not.inp_a(btype)) goto 80000
      if (.not.inp_a(name)) goto 80000
      nbasis = nbasis + 1
      if (.not.inp_read()) then
        goto 80000
      else
*        call inp_outrec()
      endif
      if (.not.inp_a(field)) goto 80000
      if (field(1:3).eq.'end') goto 00001
      if (btype(1:3).eq.'ecp') then
        oecp = .true.
      else
        oecp = .false.
      endif
      if (oecp) then
* read nelec line
        if (.not.inp_read()) then
          goto 80000
        else
*          call inp_outrec()
        endif
      endif
      if (.not.oecp.and.btype(1:5).ne.'basis') then
        write(6,*)' unrecognized btype',btype
        goto 80000
      endif
      do ip = 1,pmax
        ee(ip)  = 0.0d00
        rex(ip) = 0
        do ig = 1,cmax
          cc(ip,ig) = 0.0d00
        enddo
      enddo
00002 continue
      nprim = 0
00003 continue
      if (.not.inp_read()) then
        goto 80000
      else
*        call inp_outrec()
      endif
      call inp_set_field(0) 
      if (.not.inp_a(field)) goto 80000
      if (field(1:3).eq.'end') then
        call checkem(nwrong,ng,nprim,rex,ee,pmax,cc,oecp,name)
        goto 00001
      endif
      call inp_set_field(0) 
      if (oecp) then
        if (.not.inp_i(rex(nprim+1))) then
*rak:          if (.not.inp_read()) then
*rak:            goto 80000
*rak:          else
*rak:            call inp_outrec()
*rak:          endif
          call checkem(nwrong,ng,nprim,rex,ee,pmax,cc,oecp,name)
          goto 00002
        endif
        if (.not.inp_f(ee(nprim+1))) goto 80000
        ng = inp_n_field() - 2
        do ig = 1,ng
          if (.not.inp_f(cc(nprim+1,ig))) goto 80000
        enddo
        nprim = nprim + 1
        goto 00003
      else
        if (.not.inp_f(ee(nprim+1))) then
*rak:          if (.not.inp_read()) then
*rak:            goto 80000
*rak:          else
*rak:            call inp_outrec()
*rak:          endif
          call checkem(nwrong,ng,nprim,rex,ee,pmax,cc,oecp,name)
          goto 00002
        endif
        ng = inp_n_field() - 1
        do ig = 1,ng
          if (.not.inp_f(cc(nprim+1,ig))) goto 80000
        enddo
        nprim = nprim + 1
        goto 00003
      endif
70000 continue
      write(6,*)' number of basis sets/ecps checked :',nbasis
      write(6,*)' number that have a problem        :',nwrong
      stop ' no run errors '
80000 continue
*      call inp_outrec()
      write(6,*)'end of file or error reading'
      stop ' done with error'
c
90000 write(6,*)' error opening library '
      stop ' done with error'
      end
      subroutine checkem(nwrong,ng,nprim,rex,ee,ldc,cc,oecp,name)
      implicit none
      character*(*) name
      integer nwrong
      integer ng
      integer nprim
      integer ldc
      integer rex(nprim)
      double precision ee(nprim), cc(ldc,ng)
      logical oecp, okay
      integer ig, ip
      integer nzero
      double precision thresh
c
      thresh = 1.0d-15
      okay = .true.
c
      do ip = 1,nprim
        if (oecp)
     &      okay = okay.and.
     &      ((rex(ip).eq.0).or.(rex(ip).eq.1).or.(rex(ip).eq.2))
        okay = okay.and.(ee(ip).gt.0.0d00)
      enddo
      do ig = 1,ng
        nzero = 0
        do ip = 1,nprim
          if (abs(cc(ip,ig)).lt.thresh)nzero=nzero+1
        enddo
        okay = okay.and.nzero.lt.nprim
      enddo
      if (okay) goto 00001
      nwrong = nwrong + 1
      if (oecp) then
        write(6,'(a,a)')
     &      ' something wrong with a contraction in:ecp ',name
        write(8,'(a,a)')
     &      ' something wrong with a contraction in:ecp ',name
        call util_flush(8)
      else
        write(6,'(a,a)')
     &      ' something wrong with a contraction in:basis ',name
        write(8,'(a,a)')
     &      ' something wrong with a contraction in:basis ',name
        call util_flush(8)
      endif
      if (oecp) then
        do ip = 1,nprim
          write(6,10000)rex(ip),ee(ip),(cc(ip,ig),ig=1,ng)
        enddo
      else
        do ip = 1,nprim
          write(6,10001)ee(ip),(cc(ip,ig),ig=1,ng)
        enddo
      endif
10000 format(1x,i3,1x,30(1pd14.6))
10001 format(1x,30(1pd14.6))
00001 continue
      end
