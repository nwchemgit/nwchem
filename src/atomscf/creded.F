      subroutine creded(g_dens, dt, iiloc, nbb, ospherical)
C$Id: creded.F,v 1.9 1996-08-06 05:26:05 d3g681 Exp $
      implicit none
#include "global.fh"
#include "cguess.fh"
c
c..   routine to merge atomic density matrix into molecular one
c..   is called straight after atom, so all info is still in
c..   common /cguess/ : i.e. dt,nbc
c
      integer nbb               ! [input] Leading dimension of arrays
      integer g_dens            ! [input] Handle of GA where dens is stuffed
      real *8 dt(*)             ! [input] Atomic density
      integer iiloc(nbb,4)      ! [input] First basis function in each shell
      logical ospherical        ! [input] True if spherical basis
c
c..
      real *8 dmult(145)        ! cartesian
      real *8 smult(83)         ! spherical
c..
crz   order of d's and f's in GAMESS is completely different
crz   from MOLECULE !!
c
crjh  NWCHEM is even more whacked out ... normalization differs
c
c     dmult consists of :
c     nothing for s;
c      p functions (1 line)
c      d functions (next 4 lines)
c      f functions (last lines)
c      -0.670820393 = -0.3*sqrt(5) (used for f)
c
c
      integer k, l, m
      integer lm, nbci, kdim, kmone, kpoint
      integer ll,mm
      real *8 factor

      double precision third, sixth
      parameter (third=1.0d0/3.0d0)
      parameter (sixth=-1.0d0/6.0d0)

      data dmult/ 1.0d0, 3*0.0d0, 1.0d0, 3*0.0d0, 1.0d0,

     $        1d0,   0d0,    0d0, -0.5d0,    0d0,  -0.5d0,
     $        0d0,   3d0,    0d0,    0d0,    0d0,     0d0,
     $        0d0,   0d0,    3d0,    0d0,    0d0,     0d0,
     $     -0.5d0,   0d0,    0d0,    1d0,    0d0,  -0.5d0,
     $        0d0,   0d0,    0d0,    0d0,    3d0,     0d0,
     $     -0.5d0,   0d0,    0d0, -0.5d0,    0d0,     1d0,

     x  1.0d0,   4*0.0d0, -0.670820393d0, 0.0d0,-0.670820393d0, 2*0.0d0,
     x  0.0d0,1.0d0,0.0d0,-0.670820393d0, 4*0.0d0,-0.670820393d0, 0.0d0,
     x2*0.0d0,1.0d0,0.0d0,-0.670820393d0, 0.0d0,-0.670820393d0, 3*0.0d0,
     x   0.0d0, -0.670820393d0,  0.0d0,  1.2d0, 4*0.0d0, -0.3d0,  0.0d0,
     x 2*0.0d0, -0.670820393d0,  0.0d0,  1.2d0,  0.0d0, -0.3d0, 3*0.0d0,
     x  -0.670820393d0,         4*0.0d0, 1.2d0,  0.0d0, -0.3d0, 2*0.0d0,
     x 2*0.0d0, -0.670820393d0,  0.0d0, -0.3d0,  0.0d0,  1.2d0, 3*0.0d0,
     x  -0.670820393d0,     4*0.0d0,    -0.3d0,  0.0d0,  1.2d0, 2*0.0d0,
     x   0.0d0, -0.670820393d0, 0.0d0, -0.3d0, 4*0.0d0,  1.2d0,   0.0d0,
     x 9*0.0d0,  1.0d0/

* New broken NWChem d normalization
*     $      third,   0d0,    0d0,  sixth,    0d0,   sixth,
*     $        0d0,   1d0,    0d0,    0d0,    0d0,     0d0,
*     $        0d0,   0d0,    1d0,    0d0,    0d0,     0d0,
*     $      sixth,   0d0,    0d0,  third,    0d0,   sixth,
*     $        0d0,   0d0,    0d0,    0d0,    1d0,     0d0,
*     $      sixth,   0d0,    0d0,  sixth,    0d0,   third,


* Old NWChem d normalization
*     $        1d0,   0d0,    0d0, -0.5d0,    0d0,  -0.5d0,
*     $        0d0,   3d0,    0d0,    0d0,    0d0,     0d0,
*     $        0d0,   0d0,    3d0,    0d0,    0d0,     0d0,
*     $     -0.5d0,   0d0,    0d0,    1d0,    0d0,  -0.5d0,
*     $        0d0,   0d0,    0d0,    0d0,    3d0,     0d0,
*     $     -0.5d0,   0d0,    0d0, -0.5d0,    0d0,     1d0,


c$$$      data dmult/ 1.0d0, 3*0.0d0, 1.0d0, 3*0.0d0, 1.0d0,
c$$$     x  1.0d0, -0.5d0, -0.5d0, 3*0.0d0,
c$$$     x -0.5d0,  1.0d0, -0.5d0, 3*0.0d0,
c$$$     x -0.5d0, -0.5d0,  1.0d0, 3*0.0d0,
c$$$     x  3*0.0d0, 1.0d0, 6*0.0d0, 1.0d0, 6*0.0d0, 1.0d0,
c$$$     x  1.0d0,   4*0.0d0, -0.670820393d0, 0.0d0,-0.670820393d0, 2*0.0d0,
c$$$     x  0.0d0,1.0d0,0.0d0,-0.670820393d0, 4*0.0d0,-0.670820393d0, 0.0d0,
c$$$     x2*0.0d0,1.0d0,0.0d0,-0.670820393d0, 0.0d0,-0.670820393d0, 3*0.0d0,
c$$$     x   0.0d0, -0.670820393d0,  0.0d0,  1.2d0, 4*0.0d0, -0.3d0,  0.0d0,
c$$$     x 2*0.0d0, -0.670820393d0,  0.0d0,  1.2d0,  0.0d0, -0.3d0, 3*0.0d0,
c$$$     x  -0.670820393d0,         4*0.0d0, 1.2d0,  0.0d0, -0.3d0, 2*0.0d0,
c$$$     x 2*0.0d0, -0.670820393d0,  0.0d0, -0.3d0,  0.0d0,  1.2d0, 3*0.0d0,
c$$$     x  -0.670820393d0,     4*0.0d0,    -0.3d0,  0.0d0,  1.2d0, 2*0.0d0,
c$$$     x   0.0d0, -0.670820393d0, 0.0d0, -0.3d0, 4*0.0d0,  1.2d0,   0.0d0,
c$$$     x 9*0.0d0,  1.0d0/

      data smult /
     $     1d0, 3*0d0, 1d0, 3*0d0, 1d0, ! p
     $     1d0, 5*0d0, 1d0, 5*0d0, 1d0, 5*0d0, 1d0, 5*0d0, 1d0, !d
     $     1d0, 7*0d0, 1d0, 7*0d0, 1d0, 7*0d0, 1d0, 7*0d0, 1d0, 7*0d0,
     $     1d0, 7*0d0, 1d0/      !f


c
c.......................................................................
c         among the variables used are:
c
c              nsym       - highest l-quantum no. used in atomic calc.
c              dt         - atomic density matrix
c              d          - area for final molecular density matrix.
c.......................................................................
c..
      lm = 0
      do k = 1 , nsym
         nbci = nbc(k)
         if (k.eq.1) then
c.......................................................................
c
c       s-orbitals, simple distribution of matrix elements.
c
c.......................................................................
            do l = 1 , nbci
               do m = 1 , l
                  ll = iiloc(l,1)
                  mm = iiloc(m,1)
                  lm = lm + 1
                  call ga_put(g_dens,ll,ll,mm,mm,dt(lm),1)
                  call ga_put(g_dens,mm,mm,ll,ll,dt(lm),1)
               enddo
            enddo
         else if (k.le.4) then
            factor = 1.d0/(k+k-1)
            if (.not. ospherical) then
               kdim = k*(k+1)/2
               kmone = k - 1
               kpoint = kmone*(kmone+1)*
     $              (((6*kmone+24)*kmone+26)*kmone+4)/120
               call pdfded(k,kdim,g_dens,dt,factor,dmult(kpoint),lm,
     $              nbci, iiloc(1,k))
            else
               kdim = 2*k - 1
               if (k.eq.2) then
                  kpoint = 1
               else if (k.eq.3) then
                  kpoint = 10
               else
                  kpoint = 35
               endif
               call pdfded(k,kdim,g_dens,dt,factor,smult(kpoint),lm,
     $              nbci, iiloc(1,k))
            endif
         end if
      enddo
c..
      return
      end
