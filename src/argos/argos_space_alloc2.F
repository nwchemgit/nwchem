      subroutine argos_space_alloc2
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer isize
c
      if(nbget.ne.0) then
      msa=max(1,msa,(mbox-30+nrempr+madbox)*((nwm*nwa+nsa)/nbtot+1)+1)
      mwm=max(1,mwm,msa/mwa+1)
      msa=min(msa,2*nsa+1)
      mwm=min(mwm,2*nwm+1)
      endif
c
      if(util_print('distribution',print_default)) then
      if(me.eq.0) then
      if(nbget.ne.0) then
      write(lfnout,2005) nrempr
 2005 format(' Number of remote cell pairs: ',i7)
      if(nbget.gt.0) then
      write(lfnout,2006) nbget
 2006 format(' Number of prefetch cells:    ',i7)
      endif
      endif
      write(lfnout,2007) mwm,mwmg
 2007 format(' Dimension solvent local:     ',i7,', global:',i7)
      write(lfnout,2008) msa,msag
 2008 format(' Dimension solute local:      ',i7,', global:',i7)
      endif
      endif
c
      if(.not.ga_create(mt_int,np*mwmg,miw2,'iw',mwmg,miw2,ga_iw))
     + call md_abort('Failed to create global array iw',0)
      isize=6+12*mwa
      if(lpola) isize=6+18*mwa
      if(lpola.and.lfree) isize=6+30*mwa
      if(.not.ga_create(mt_dbl,np*mwmg,isize,'w',mwmg,isize,ga_w))
     + call md_abort('Failed to create global array w',0)
      if(.not.ga_create(mt_int,np*msag,mis2,'is',msag,mis2,ga_is))
     + call md_abort('Failed to create global array is',0)
      isize=39
      if(lpola) isize=45
      if(lpola.and.lfree) isize=57
      if(.not.ga_create(mt_dbl,np*msag,isize,'s',msag,isize,ga_s))
     + call md_abort('Failed to create global array s',0)
c
      if(.not.ga_create(mt_int,np*mwmg,1,'iwz',mwmg,1,ga_iwz))
     + call md_abort('Failed to create global array iwz',0)
      if(.not.ga_create(mt_int,np*msag,1,'isz',msag,1,ga_isz))
     + call md_abort('Failed to create global array isz',0)
c
      return
      end
c $Id$
