      subroutine grid_finpts(dosympts,
     .        rad_sh,xyz, ictr, 
     .        irsh, iang, iang_done,nq_left,
     ,        aqbuf, nq, fnxyz,nxyz,
     ,        gridscr,nops,symops,gops,
     .        iscr_loc,
     .     ncenters, rij,
     .     rq,bp,expo,
     .     oskel,symfact)

      implicit none
c
      logical dosympts
      double precision rad_sh
      double precision xyz(3,*)
      integer ictr,irsh,iang,iang_done,nq_left
      double precision aqbuf(4,*)
      integer nq,fnxyz,nxyz
      double precision  gridscr(*)
      integer nops
      double precision symops(*),gops(*)
      integer iscr_loc(*)
c
      integer ncenters
      double precision rij
      double precision rq(nq,*),bp(nq,*),expo(*),symfact
      logical oskel

c     
c     Generate all angular quadrature points for this radial shell
c     and keep only those inside cube
c     
      if(dosympts) then
         call grid_sympts(dosympts,
     .        rad_sh,xyz(1,ictr), ictr, 
     .        irsh, iang, iang_done,nq_left,
     ,        aqbuf, nq, fnxyz,nxyz,
     ,        gridscr,nops,symops,gops,
     .        iscr_loc)
      else
         call grid_qpts(rad_sh,xyz(1,ictr), ictr, 
     ,        irsh, iang, 
     .        aqbuf, nq, 
     *        fnxyz,nxyz)
      endif
c     
c     modify the weight to include the spatial weight (and symmetry factor)
c     
      if(ncenters.gt.1) then
         call grid_setspac(xyz, aqbuf,
     ,        rij, 
     &        rq, bp, ictr, nq,
     ,        expo, rad_sh)
      endif
      if(nq.ne.0) then
         if(oskel) call dscal(
     /        nq, symfact, aqbuf(4,1), 4)
      endif

      return
      end
c $Id$
