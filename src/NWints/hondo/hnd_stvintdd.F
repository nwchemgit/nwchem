c
c $Id: hnd_stvintdd.F,v 1.1 1999-04-28 22:48:09 windus Exp $
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
      double precision zero, negone
      data zero   /0.0d+00/
      data negone /-1.0d+00/
      data minprim / 1, 2, 5,11,21/
      data maxprim / 1, 4,10,20,35/
c
      write(6,*) 'lscr is ',lscr
      write(6,*) 'i_cent and j_cent',i_cent,j_cent
      write(6,*) 'nint is ',nint
      call flush(6)
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
c       call dscal(nint*3*3,negone,sder(nint*3*3+1),1)
      endif
c
c     write(luout,*) 'before the call to doV'
c     call flush(luout)
      if (doV) then
        call dfill(nint*3*3*nat,zero,vder,1)
        do i=1,nat
          write(luout,*) 'number of calls to hnd_hlf_ij2 ',i
          call hnd_hlf_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1      xyzj,expj,coefj,j_nprim,j_ngen,Lj,i,scr,
     2      vder((i-1)*3*3*nint+1),zan,xyz,nder,nint)
        enddo
c        write(6,*) 'before the second doV, vder is'
c        do i=1,nat*3*nint
c          write(6,*) vder(i)
c        enddo
c     endif

      write(luout,*) 'before the doT call'
      call flush(luout)
      if (doT.or.doV) then
        if (doT) call dfill(nint*3*3*2,zero,tder,1)
c note that scr is where the vder ints are temporarily going
        if (doV) call dfill(nint*3*3*4,zero,scr,1)
c       call hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
c    1    xyzj,expj,coefj,j_nprim,j_ngen,Lj,xyz,scr,tder,
c    2    tder(nint*3*3+1),scr,scr,scr(nint*3*3+1),nder,nint,
c    3    zan,nat,doT,doV,i_cent,j_cent)
         write(6,*) 'after the first tvd_ij call, scr is'
         do i=1,3*nint
           write(6,*) scr(i)
         enddo
c note that the first scr is where the vder ints are going
c       call hnd_tvd_ij2(xyzj,expj,coefj,j_nprim,j_ngen,Lj,
c    1    xyzi,expi,coefi,i_nprim,i_ngen,Li,xyz,tder,scr,scr,
c    2    scr,scr(nint*3*3*2+1),scr(nint*3*3*3+1),nder,nint,
c    3    zan,nat,.false.,doV,j_cent,i_cent)
c        write(6,*) 'after the second tvd_ij call, scr is'
c        do i=3*nint+1,6*nint
c          write(6,*) scr(i)
c        enddo
c        write(6,*) 'after the second doV, vder is'
c        do i=1,nat*3*nint
c          write(6,*) vder(i)
c        enddo
        if (doT) then
          call dcopy(nint*3*3,tder,1,tder(nint*3*3+1),1)
          call dscal(nint*3*3,negone,tder(nint*3*3+1),1)
        endif
c       if (doV) then
c         if (i_cent.eq.j_cent) then
c           write(6,*) 'i equals j'
c         if (i_cent.eq.j_cent) write(6,*) 'i equals j'
c           call daxpy(nint*3,1.D0,scr,1,vder((i_cent-1)*3*nint+1),1)
c           write(6,*) 'after the first daxpy, vder is'
c           do i=1,nat*3*nint
c             write(6,*) vder(i)
c           enddo
c           call daxpy(nint*3,1.D0,scr(nint*3+1),1,
c    1                 vder((j_cent-1)*3*nint+1),1)
c           write(6,*) 'after the second daxpy, vder is'
c           do i=1,nat*3*nint
c             write(6,*) vder(i)
c           enddo
c         else
c           write(6,*) 'i not equal j'
c           call daxpy(nint*3,1.D0,scr,1,vder((i_cent-1)*3*nint+1),1)
c           call daxpy(nint*3,1.D0,scr(nint*3+1),1,
c    1                 vder((j_cent-1)*3*nint+1),1)
c         endif
        endif
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
