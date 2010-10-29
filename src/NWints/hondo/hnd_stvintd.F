c
c $Id$
c
c Taken and modified from HONDO
c
      subroutine hnd_stvintd(xyzi,expi,coefi,i_nprim,i_ngen, Li, 
     1  i_cent,xyzj,expj,coefj, j_nprim, j_ngen, Lj, j_cent,
     2  xyz,zan,nat,sder,tder,vder,nint,doS,doT,doV,scr,lscr)
c
      implicit none
#include "stdio.fh"
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer nat      ! [input] number of atoms
      integer lscr     ! [input] size of scratch array
      integer nint     ! [input] size of any integral buffer
      integer i_cent   ! [input] lexical number of center i
      integer j_cent   ! [input] lexical number of center j
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,nat)  ! [input] all atom positions
      double precision zan(nat)  ! [input] charges on all atoms
      double precision scr(lscr)  ! [scratch] scratch buffers
      double precision sder(*)  ! [output] overlap integrals
      double precision tder(*)  ! [output] kinetic energy integrals
      double precision vder(*)  ! [output] potential integrals
      logical doS     ! [input] compute overlap (True/False)
      logical doT     ! [input] compute kinetic (True/False)
      logical doV     ! [input] compute potential (True/False)
      logical some, doH1
      integer nder, maxi, maxj
      integer i, index
      double precision zero
      data zero   /0.0d+00/
c
      doH1 = doT.and.doV.and..not.doS
      nder = 1
c
c     ----- calculate -s-, -t-, and -v- derivative integrals -----
c
      some = .false.
      if(some) write(luout,9999)
c
      if (doS) then
        call dcopy(nint*3*2,zero,0,sder,1)
        call hnd_sd_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1    xyzj,expj,coefj, j_nprim, j_ngen, Lj, sder, sder(nint*3+1),
     2    scr, scr, scr, nder,nint,scr,lscr)
      endif
c
      if (doV) then
        call dcopy(nint*3*nat,zero,0,vder,1)
        do i=1,nat
          call hnd_hlf_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1      xyzj,expj,coefj,j_nprim,j_ngen,Lj,i,vder((i-1)*3*nint+1),
     2      scr,zan,xyz,nder,nint,scr,lscr)
        enddo
      endif

      if (doT.or.doV) then
        if (doT) call dcopy(nint*3*2,zero,0,tder,1)
        if (doV) then
           call hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1       xyzj,expj,coefj,j_nprim,j_ngen,Lj,xyz,tder,tder(nint*3+1),
     2       scr,scr,scr,vder((i_cent-1)*3*nint+1),
     3       vder((j_cent-1)*3*nint+1),scr,scr,scr,nder,nint,zan,
     4       nat,doT,doV,i_cent,j_cent,scr,lscr)
        else
           call hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1       xyzj,expj,coefj,j_nprim,j_ngen,Lj,xyz,tder,tder(nint*3+1),
     2       scr,scr,scr,scr,scr,scr,scr,scr,nder,nint,zan,
     3       nat,doT,doV,i_cent,j_cent,scr,lscr)
        endif
      endif
c
      if(some) write(luout,9998)
      return
 9999 format(/,10x,20(1h-),/,10x,'1 electron derivative integrals',
     2       /,10x,20(1h-))
 9998 format(' ...... end of one-electron integrals ......')
      end
