c
c $Id: hnd_stvintd.F,v 1.4 1999-12-29 23:30:24 windus Exp $
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
      integer nder, mini, maxi, minj, maxj, minprim, maxprim
      integer i, index
      dimension minprim(5), maxprim(5)
      double precision zero
      data zero   /0.0d+00/
      data minprim / 1, 2, 5,11,21/
      data maxprim / 1, 4,10,20,35/
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
        call dfill(nint*3*2,zero,sder,1)
        call hnd_sd_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1    xyzj,expj,coefj, j_nprim, j_ngen, Lj, sder, sder(nint*3+1),
     2    scr, scr, scr, nder,nint)
      endif
c
      if (doV) then
        call dfill(nint*3*nat,zero,vder,1)
        do i=1,nat
          call hnd_hlf_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1      xyzj,expj,coefj,j_nprim,j_ngen,Lj,i,vder((i-1)*3*nint+1),
     2      scr,zan,xyz,nder,nint)
        enddo
      endif

      if (doT.or.doV) then
        if (doT) call dfill(nint*3*2,zero,tder,1)
c note that scr is where the vder ints are temporarily going
        if (doV) call dfill(nint*3*2,zero,scr,1)
        call hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1    xyzj,expj,coefj,j_nprim,j_ngen,Lj,xyz,tder,tder(nint*3+1),
     2    scr,scr,scr,scr,scr((nint*3)+1),scr,scr,scr,nder,nint,zan,
     3    nat,doT,doV,i_cent,j_cent)
c
c Combine the orbital and operator derivative parts
c
        if (doV) then
          call daxpy(nint*3,1.D0,scr,1,vder((i_cent-1)*3*nint+1),1)
          call daxpy(nint*3,1.D0,scr(nint*3+1),1,
     1                vder((j_cent-1)*3*nint+1),1)
        endif
      endif
c
c     ----- perform transformation to NWChem order if necessary ----
c
      if ((Li.gt.1).or.(Lj.gt.1)) then
        mini = minprim(Li+1)
        maxi = maxprim(Li+1)
        minj = minprim(Lj+1)
        maxj = maxprim(Lj+1)
        if (doS) then
          index = 1
          do i=1,6
            call hnd_to_nw_order(scr,sder(index),minj,mini,
     1                           maxj-minj+1,maxi-mini+1)
            call dcopy (nint,scr,1,sder(index),1)
            index = index + nint
          enddo
        endif
        if (doT) then
          index = 1
          do i=1,6
            call hnd_to_nw_order(scr,tder(index),minj,mini,
     1                           maxj-minj+1,maxi-mini+1)
            call dcopy (nint,scr,1,tder(index),1)
            index = index + nint
          enddo
        endif
        if (doV) then
          index = 1
          do i=1,nat*3
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
