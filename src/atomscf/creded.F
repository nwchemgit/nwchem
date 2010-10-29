      subroutine creded(g_dens, dt, iiloc, nbb, ospherical)
C$Id$
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

     D        1d0,   0d0,    0d0, -0.5d0,    0d0,  -0.5d0,
     D        0d0,   3d0,    0d0,    0d0,    0d0,     0d0,
     D        0d0,   0d0,    3d0,    0d0,    0d0,     0d0,
     D     -0.5d0,   0d0,    0d0,    1d0,    0d0,  -0.5d0,
     D        0d0,   0d0,    0d0,    0d0,    3d0,     0d0,
     D     -0.5d0,   0d0,    0d0, -0.5d0,    0d0,     1d0,

     F   0.12501d00,-0.00949d00,-0.00677d00, 0.74987d00, 0.02698d00,
     F  -1.12488d00, 0.00266d00, 0.00674d00, 0.00150d00, 0.00001d00,
     F  -0.00949d00,11.99983d00, 0.01345d00, 0.01873d00, 0.02677d00,
     F   0.00973d00,-3.00018d00, 0.02706d00,-2.99928d00,-0.01350d00,
     F  -0.00677d00, 0.01345d00,10.12512d00, 0.02007d00, 0.02497d00,
     F   0.00025d00, 0.02703d00,-6.74884d00,-0.09453d00,-1.12542d00,
     F   0.74987d00, 0.01873d00, 0.02007d00, 4.49894d00, 0.16215d00,
     F  -6.74852d00,-0.00542d00, 0.00003d00,-0.00247d00,-0.00670d00,
     F   0.02698d00, 0.02677d00, 0.02497d00, 0.16215d00, 0.00600d00,
     F  -0.24309d00, 0.00007d00,-0.01564d00,-0.02698d00,-0.00311d00,
     F  -1.12488d00, 0.00973d00, 0.00025d00,-6.74852d00,-0.24309d00,
     F  10.12319d00,-0.00257d00,-0.02025d00,-0.00203d00, 0.00667d00,
     F   0.00266d00,-3.00018d00, 0.02703d00,-0.00542d00, 0.00007d00,
     F  -0.00257d00, 1.99961d00, 0.05401d00,-2.99864d00,-0.02701d00,
     F   0.00674d00, 0.02706d00,-6.74884d00, 0.00003d00,-0.01564d00,
     F  -0.02025d00, 0.05401d00, 4.50381d00,-0.18909d00, 0.74835d00,
     F   0.00150d00,-2.99928d00,-0.09453d00,-0.00247d00,-0.02698d00,
     F  -0.00203d00,-2.99864d00,-0.18909d00,11.99521d00, 0.09454d00,
     F   0.00001d00,-0.01350d00,-1.12542d00,-0.00670d00,-0.00311d00,
     F   0.00667d00,-0.02701d00, 0.74835d00, 0.09454d00, 0.12570d00/

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
