      subroutine argos_fopen(lclose)
c
      implicit none
c
#include "argos.fh"
c
      logical lclose
c
      if(nfnewf.gt.0) then
      ibatch=ibatch+1
      write(filtrj,1000) root(1:index(root,' ')-1),
     + ibatch,'.trj'
      write(filhop,1000) root(1:index(root,' ')-1),
     + ibatch,'.hop'
      write(filprp,1000) root(1:index(root,' ')-1),
     + ibatch,'.prp'
      write(filpmf,1000) root(1:index(root,' ')-1),
     + ibatch,'.pmf'
      write(filtim,1000) root(1:index(root,' ')-1),
     + ibatch,'.tim'
      write(rfile,1000) root(1:index(root,' ')-1),
     + ibatch,'.rst'
 1000 format(a,i3.3,a)
      endif
c
      if(ntype.ne.3) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      if(lclose) close(lfntrj)
      open(unit=lfntrj,file=filtrj(1:index(filtrj,' ')-1),
     + form='formatted',status='unknown')
      call argos_cafe_trjhdr(lfntrj)
      endif
      if(nfhop.gt.0) then
      if(lclose) close(lfnhop)
      open(unit=lfnhop,file=filhop(1:index(filhop,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nfprop.gt.0) then
      if(lclose) close(lfnprp)
      open(unit=lfnprp,file=filprp(1:index(filprp,' ')-1),
     + form='formatted',status='unknown')
      call argos_prop_header()
      endif
      endif
      if(ntype.eq.2.and.iprpmf.ne.0) then
      if(lclose) close(lfnpmf)
      open(unit=iabs(lfnpmf),file=filpmf(1:index(filpmf,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nftime.gt.0) then
      if(lclose) close(lfntim)
      open(unit=lfntim,file=filtim(1:index(filtim,' ')-1),
     + form='formatted',status='unknown')
      call argos_hdrtim()
      endif
c
      return
      end
