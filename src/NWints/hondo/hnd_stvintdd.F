c
c $Id: hnd_stvintdd.F,v 1.3 1999-12-21 23:54:04 windus Exp $
c
c Taken and modified from HONDO
c
      subroutine hnd_stvintdd(xyzi,expi,coefi,i_nprim,i_ngen, Li, 
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
      integer nder, mini, maxi, minj, maxj, minprim, maxprim
      integer i, index
      dimension minprim(5), maxprim(5)
      double precision zero
      data zero   /0.0d+00/
      data minprim / 1, 2, 5,11,21/
      data maxprim / 1, 4,10,20,35/
c
      doH1 = doT.and.doV.and..not.doS
      nder = 2
c
c     ----- calculate -s-, -t-, and -v- derivative integrals -----
c
      some = .false.
      if(some) write(luout,9999)
c
      if (doS) then
        call dfill(nint*3*3*4,zero,sder,1)
        call hnd_sd_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1    xyzj,expj,coefj, j_nprim, j_ngen, Lj, scr, scr, sder, 
     2    sder(nint*3*3+1),sder(nint*3*3*3+1),nder,nint)
c
c to get the djdi term which is equal to the didj term which is calculated
c
        call dcopy(nint*3*3,sder(nint*3*3+1),1,sder(nint*3*3*2+1),1)
      endif
c
      if (doV) then
        call dfill(nint*3*3*nat*3,zero,vder,1)
c
c Do the Helman-Feynman part (i|d2/dxcdxcV|j)
c The first nat*nint*3*3 holds these terms
c
        do i=1,nat
          call hnd_hlf_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1      xyzj,expj,coefj,j_nprim,j_ngen,Lj,i,scr,
     2      vder((i-1)*3*3*nint+1),zan,xyz,nder,nint)
        enddo
c
c Do the mixed part (di|d/dxcV|j) and (i|d/dxcV|dj)
c These are held in the rest of vder
c
         do i=1,nat
           call hnd_hlfd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1       xyzj,expj,coefj,j_nprim,j_ngen,Lj,i,
     2       vder((i-1)*3*3*nint+nint*3*3*nat+1),
     3       vder((i-1)*3*3*nint+nint*3*3*nat*2+1),
     4       zan,xyz,nder,nint)
         enddo
c
      endif

      if (doT.or.doV) then
        if (doT) call dfill(nint*3*3*3,zero,tder,1)
c note that scr is where the didj vder ints are temporarily going
        if (doV) call dfill(nint*3*3,zero,scr,1)
        call hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1    xyzj,expj,coefj,j_nprim,j_ngen,Lj,xyz,scr,scr,tder,
     2    tder(nint*3*3+1),tder(nint*3*3*2+1),scr,scr,
     3    vder(nint*3*3*(i_cent-1)+1),scr,
     4    vder(nint*3*3*(j_cent-1)+1),nder,nint,
     5    zan,nat,doT,doV,i_cent,j_cent)
c
c add in the scr vder ints
c
         call daxpy(nint*3*3,1.0d00,scr(1),1,
     1     vder(nint*3*3*nat+nint*3*3*(j_cent-1)+1),1)
         call daxpy(nint*3*3,1.0d00,scr(1),1,
     1     vder(nint*3*3*nat*2+nint*3*3*(i_cent-1)+1),1)
      endif
c
c     ----- perform transformation to NWChem order if necessary ----
c  bring back after get everything else working
c
c     if ((Li.gt.1).or.(Lj.gt.1)) then
c       mini = minprim(Li+1)
c       maxi = maxprim(Li+1)
c       minj = minprim(Lj+1)
c       maxj = maxprim(Lj+1)
c       if (doS) then
c         index = 1
c         do i=1,6
c           call hnd_to_nw_order(scr,sder(index),minj,mini,
c    1                           maxj-minj+1,maxi-mini+1)
c           call dcopy (nint,scr,1,sder(index),1)
c           index = index + nint
c         enddo
c       endif
c       if (doT) then
c         index = 1
c         do i=1,6
c           call hnd_to_nw_order(scr,tder(index),minj,mini,
c    1                           maxj-minj+1,maxi-mini+1)
c           call dcopy (nint,scr,1,tder(index),1)
c           index = index + nint
c         enddo
c       endif
c       if (doV) then
c         call hnd_to_nw_order(scr,vder,minj,mini,maxj-minj+1,
c    1                         maxi-mini+1)
c         call dcopy (nint,scr,1,vder,1)
c       endif
c     endif
c
      if(some) write(luout,9998)
      return
 9999 format(/,10x,20(1h-),/,10x,'1 electron derivative integrals',
     2       /,10x,20(1h-))
 9998 format(' ...... end of one-electron integrals ......')
      end
