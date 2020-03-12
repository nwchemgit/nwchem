      Subroutine grid_list(rtdb, qsh, nqlen,nq_tot)
c
c$Id$
c
      implicit none
#include "errquit.fh"
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
      integer rtdb, qsh
      integer nqlen ! [in]
      integer nq_tot ! [out]
      dimension qsh(3,nqlen)
      double precision ktheta
      Parameter (ktheta = 5.d0)
      double precision eps
      parameter (eps = 1.d-20)
      logical  oprint, oprint_grid, no_prune, isecp
      double precision asg1(4,3)
      integer iti, ictr, irsh,
     &        iang, iang_mx, iqsh, nr, nserie, izat, ir,
     &        me, nq_ctr_tot
      double precision r, r1, r2, r3, 
     &                 r4, RR, q1
      integer xc_leb_nrd2acc
      external xc_leb_nrd2acc
      integer xc_dau_nrd2acc
      external xc_dau_nrd2acc
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
         isecp=geom_ecp_get(geom,ictr)
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
                  call errquit('grid_list: raise MXqshells to ',iqsh,
     &       INPUT_ERR)
               endif
c
c              iang =  iangquad for lebedev
c              iang =  nagrid (or nang) for user specified gauss-legendre quadrature
c              iang = grid_params (or nang) for all other gauss-legendre quadratures
c
               if (leb)then
                  iang_mx = xc_leb_nrd2acc(nang(iti))
               else
                  iang_mx = nang(iti)
               endif
               iang = iang_mx
               if (no_prune)goto 405
c            
               if (leb) then ! lebedev angular quadrature
                  if (iang_mx.eq.1) then
                     iang = iang_mx
                  elseif (iang_mx.eq.2) then 
                     if (r.lt.r1) then
                        iang = 1
                     elseif (r.ge.r1.and.r.lt.r4) then
                        iang = 2
                     else
                        iang = 1
                     endif
                  else
                     if (r.lt.r1) then
                        iang = 2
                     elseif (r.ge.r1.and.r.lt.r2) then
                        iang = 4
                        if(isecp) iang=iang_mx
c                       if(isecp) iang=max(14,iang_mx-3)
                     elseif (r.ge.r2.and.r.lt.r3) then
                        iang = max(iang_mx-1,4)
                     elseif (r.ge.r3.and.r.lt.r4) then
                        iang = iang_mx
                     else
                        iang = max(iang_mx-1,1)
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
cedo               if (iang .gt. iang_mx) then
cedo                  write(0,*) ' iang iang_mx ',iang,iang_mx
cedo                  call errquit('Grid gen error in grid_list.',911)
cedo               endif 
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
                  nqmax = max(nqmax,nang_leb_pts)
                  nq_ctr_tot = nq_ctr_tot + nang_leb_pts
               else
                  nq_ctr_tot = nq_ctr_tot + ntheta(iang)*nphi(iang)
                  nqmax = max(nqmax,ntheta(iang)*nphi(iang))
               endif
c
c
  410    continue
c
         if (ga_nodeid().eq.0.and.oprint_grid)then
            write(LuOut,*)' ictr, iti, nq_ctr_tot: ',
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
         write(LuOut,*)' nq_tot = ', nq_tot
c
         write(LuOut,*)' nqshells, nqmax, nqlen: ',
     &               nqshells, nqmax, nqlen
c
c
      endif
      return
      end
