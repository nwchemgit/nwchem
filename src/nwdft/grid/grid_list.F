      Subroutine grid_list(rtdb, qsh, nqlen)
c
c$Id: grid_list.F,v 1.3 1999-06-14 22:39:22 d3h449 Exp $
c
      implicit none
c      
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "geom.fh"
#include "sym.fh"
c
      integer rtdb, qsh, nqlen
      dimension qsh(3,nqlen)
      double precision ktheta
      Parameter (ktheta = 5.d0)
      double precision eps
      parameter (eps = 1.d-20)
      logical leb, oprint, oprint_grid, no_prune
      double precision asg1(4,3)
      integer iti, ictr, irsh,
     &        iang, iang_mx, iqsh, nr, nserie, izat, ir,
     &        me, nq_ctr_tot, nq_tot
      double precision R_MAX, r, r1, r2, r3, 
     &                 r4, RR, q1

      Data asg1/0.25  , 0.5, 1.0, 4.5,
     &          0.1667, 0.5, 0.9, 3.5,
     &          0.1   , 0.4, 0.8, 2.5/
c
c     determine node id
c
      me = ga_nodeid()
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
      leb = .false.
      if (iangquad.ne.0) leb = .true.
      R_MAX = dble(irqmax)
c
c     check for no pruning; no_prune
c
      if (.not. rtdb_get(rtdb, 'dft:no_prune', mt_log, 1,
     &   no_prune))then
         no_prune = .false.
      endif
c
      iqsh = 0
      nqmax = 0
      nq_tot = 0
c
c     grid list is generated here (qsh); with or without pruning
c
      do 420 ictr = 1, ncenters
         iti = iatype(ictr)
         if (bsrad_atom_type(iti).le.EPS) goto 420
c
c        check if this atom center is unique
c
         if (oskel)then
            if (.not.sym_atom(geom, ictr, q1))goto 420
         endif
         if (ga_nodeid().eq.0.and.oprint_grid)then
            write(*,*)' ncenters, ictr, iti, bsrad_atom_type(iti)',
     &                  ncenters, ictr, iti, bsrad_atom_type(iti)
         endif
c
c        parameter for pruning
c
         if (leb) then
            izat = nint(znuc_atom_type(iti))
            RR = bsrad_atom_type(iti)
            if (izat.lt.3)then
               nserie = 1
            elseif (izat.gt.2.and.izat.lt.11) then
               nserie = 2
            else
               nserie = 3
            endif
            r1 = asg1(1,nserie)*RR
            r2 = asg1(2,nserie)*RR
            r3 = asg1(3,nserie)*RR
            r4 = asg1(4,nserie)*RR
         endif
c
         nr = nrad(iti)
c
         nq_ctr_tot = 0
c
         if (ga_nodeid().eq.0.and.oprint_grid)then
            write(*,*)' nr, iti, nrad(iti), nang(iti): ', 
     &                  nr, iti, nrad(iti), nang(iti)
         endif
c
         do 410 ir = 1, nr
            irsh = nr - (ir - 1)
            r = rpts(irsh,iti)
c
               iqsh = iqsh + 1
c
               if (iqsh.gt.MXqshells)then
                  write(LuOut,*)
                  write(luout,*)'grid_list: Too many quadrature shells.'
                  call errquit('grid_list: raise MXqshells to ',iqsh)
               endif
c
c              iang =  iangquad for lebedev
c              iang =  nagrid (or nang) for user specified gauss-legendre quadrature
c              iang = grid_params (or nang) for all other gauss-legendre quadratures
c
               if (leb)then
                  iang_mx = iangquad
               else
                  iang_mx = nang(iti)
               endif
               iang = iang_mx
               if (no_prune)goto 405
c            
               if (leb) then ! lebedev angular quadrature
                  if (iangquad.eq.1) then
                     iang = iangquad
                  elseif (iangquad.eq.2) then 
                     if (r.lt.r1) then
                        iang = 1
                     elseif (r.ge.r1.and.r.lt.r4) then
                        iang = 2
                     else
                        iang = 1
                     endif
                  else
                     if (r.lt.r1) then
                        iang = 1
                     elseif (r.ge.r1.and.r.lt.r2) then
                        iang = 2
                     elseif (r.ge.r2.and.r.lt.r3) then
                        iang = max(iangquad-1,1)
                     elseif (r.ge.r3.and.r.lt.r4) then
                        iang = iangquad
                     else
                        iang = max(iangquad-1,1)
                     endif
                  endif
               else ! gauss-legendre angular quadrature
                  rr = ktheta*r/bsrad_atom_type(iti)
                  iang = min(nang(iti), int(dble(nang(iti))*rr))
                  iang = max(iang,4)
               endif
c
c              check iang value
c
               if (iang .gt. iang_mx) then
                  call errquit('Grid gen error in grid_list.',911)
               endif 
  405          continue
c
c              temp continue to fill common iqsh to keep code working
c
               iqshell(1,iqsh) = irsh
               iqshell(2,iqsh) = 1
               iqshell(3,iqsh) = ictr
               iqshell(4,iqsh) = iang
c
               qsh(1,iqsh) = irsh
               qsh(2,iqsh) = ictr
               qsh(3,iqsh) = iang
c               write(*,*)' qsh(1-3,iqsh): ',
c     &         qsh(1,iqsh), qsh(2,iqsh), 
c     &         qsh(3,iqsh)
c
c              determine nqmax
c
               if (leb)then
                  nqmax = nang_leb_pts
               else
                  nqmax = max(nqmax,ntheta(iang)*nphi(iang))
               endif
c
            nq_ctr_tot = nq_ctr_tot + ntheta(iang)*nphi(iang)
c
  410    continue
c
         if (ga_nodeid().eq.0.and.oprint_grid)then
            write(6,*)' ictr, iti, nq_ctr_tot: ',
     &                  ictr, iti, nq_ctr_tot
         endif
c
         nq_tot = nq_tot + nq_ctr_tot      
c
  420 continue
c
      nqshells = iqsh
      nqmax = nqmax*nquad_task
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
c
         write(6,*)' nq_tot = ', nq_tot
c
         write(6,*)' nqshells, nqmax, nqlen: ',
     &               nqshells, nqmax, nqlen
c
         write(6,*)' iqshell: '
         call i_output(iqshell, 1, 4, 1, nqshells, 4, mxqshells, 1)
c
         write(6,*)' qsh: '
         call i_output(qsh, 1, 3, 1, nqshells, 3, nqlen, 1)
c
      endif
      return
      end
