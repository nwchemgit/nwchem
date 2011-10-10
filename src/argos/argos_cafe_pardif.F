      subroutine argos_cafe_pardif(wgt,vdw,chg,iwat,iwqt,mwb,
     + nwb,nbp,bndw,
     + mwh,nwh,nhp,angw,mwd,nwd,ndp,dihw,mwo,nwo,nop,oopw,
     + msb,npb,bnds,msh,nph,angs,msd,npd,dihs,mso,npo,oops)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 wgt(mat,mset),vdw(mat,mat,map,mset),chg(mqt,mqp,mset)
      integer iwat(mwa),iwqt(mwa)
      integer mwb,mwh,mwd,mwo,nbp,nhp,ndp,nop
      integer nwb,nwh,nwd,nwo
      integer msb,msh,msd,mso,npb,nph,npd,npo
      real*8 bndw(mwb,nbp,6),angw(mwh,nhp,6)
      real*8 dihw(mwd,ndp,6),oopw(mwo,nop,6)
      real*8 bnds(msb,npb,6),angs(msh,nph,6)
      real*8 dihs(msd,npd,6),oops(mso,npo,6)
c
      return
      end
c $Id$
