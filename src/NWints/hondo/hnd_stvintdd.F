c
c $Id: hnd_stvintdd.F,v 1.6 2000-08-16 19:48:15 windus Exp $
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
      logical some
      integer nder, mini, maxi, minj, maxj, minprim, maxprim
      integer i, index
      dimension minprim(5), maxprim(5)
      double precision zero, negone
      data zero   /0.0d+00/
      data negone /-1.0d+00/
      data minprim / 1, 2, 5,11,21/
      data maxprim / 1, 4,10,20,35/
c
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
c     write(luout,*) 'before the call to doV'
c     call flush(luout)
c
      if (doV) then
        call dfill(nint*3*3*nat*3,zero,vder,1)
c
c Do the Helman-Feynman part (i|d2/dxcdxcV|j)
c The first nat*nint*3*3 holds these terms
c
        do i=1,nat
c         write(luout,*) 'number of calls to hnd_hlf_ij2 ',i
          call hnd_hlf_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1      xyzj,expj,coefj,j_nprim,j_ngen,Lj,i,scr,
     2      vder((i-1)*3*3*nint+1),zan,xyz,nder,nint)
        enddo
c       write(6,*) 'Hel-Fey part'
c        do i=1,nat*3*3*3*nint
c          write(6,*) vder(i)
c        enddo
c
c Do the mixed part (di|d/dxcV|j) and (i|d/dxcV|dj)
c These are held in the rest of vder
c
         do i=1,nat
c         write(luout,*) 'number of calls to hnd_hlfd_ij ',i
           call hnd_hlfd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1       xyzj,expj,coefj,j_nprim,j_ngen,Lj,i,
     2       vder((i-1)*3*3*nint+nint*3*3*nat+1),
     3       vder((i-1)*3*3*nint+nint*3*3*nat*2+1),
     4       zan,xyz,nder,nint)
         enddo
c
c        write(6,*) 'before the second doV, vder is'
c        do i=1,nat*3*3*3*nint
c          write(6,*) vder(i)
c        enddo
      endif

c     write(luout,*) 'before the doT call'
c     call flush(luout)
      if (doT.or.doV) then
        call dfill(nint*3*3*3,zero,tder,1)
c note that scr is used as a dummy and tder also holds the vder terms
c       if (doV) call dfill(nint*3*3,zero,scr,1)
        call hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1    xyzj,expj,coefj,j_nprim,j_ngen,Lj,xyz,scr,scr,tder,
     2    tder(nint*3*3+1),tder(nint*3*3*2+1),scr,scr,
     3    tder,tder(nint*3*3+1),
     4    tder(nint*3*3*2+1),nder,nint,
     5    zan,nat,doT,doV,i_cent,j_cent)
c        write(6,*) 'after the second doV, vder is'
c        do i=1,nat*3*3*3*nint
c          write(6,*) vder(i)
c        enddo
      endif
c
c     ----- perform transformation to NWChem order if necessary ----
c  bring back after get everything else working
c
      if ((Li.gt.1).or.(Lj.gt.1)) then
        mini = minprim(Li+1)
        maxi = maxprim(Li+1)
        minj = minprim(Lj+1)
        maxj = maxprim(Lj+1)
        if (doS) then
          index = 1
          do i=1,36
            call hnd_to_nw_order(scr,sder(index),minj,mini,
     1                           maxj-minj+1,maxi-mini+1)
            call dcopy (nint,scr,1,sder(index),1)
            index = index + nint
          enddo
        endif
        if (doV) then
          index = 1
          do i=1,27
            call hnd_to_nw_order(scr,tder(index),minj,mini,
     1                           maxj-minj+1,maxi-mini+1)
            call dcopy (nint,scr,1,tder(index),1)
            index = index + nint
          enddo
          index = 1
          do i=1,27*nat
            call hnd_to_nw_order(scr,vder(index),minj,mini,
     1                           maxj-minj+1,maxi-mini+1)
            call dcopy (nint,scr,1,vder(index),1)
            index = index + nint
          enddo
        endif
      endif
c
      if(some) write(luout,9998)
      return
 9999 format(/,10x,20(1h-),/,10x,'1 electron derivative integrals',
     2       /,10x,20(1h-))
 9998 format(' ...... end of one-electron integrals ......')
      end
