      subroutine argos_main()
c
c $Id: argos_main.F,v 1.216 2008-07-03 18:05:48 d3j191 Exp $
c
      implicit none
c
#include "argos.fh"
c
      if(ntype.eq.0) then
c
c     single energy
c     -------------
c
      if(nftri.eq.0) then
      call argos_sp()
      else
      call argos_spi()
      endif
c
      elseif(ntype.eq.1) then
c
c     energy minimization
c     -------------------
c
      call argos_em()
c
      elseif(ntype.eq.2) then
c
c     molecular dynamics
c     ------------------
c
      call argos_md()
c
      elseif(ntype.eq.3) then
c
c     free energy simulation
c     ----------------------
c
      call argos_ti()
c
      else
      call md_abort('Unknown calculation type',ntype)
      endif
c
      return
      end
