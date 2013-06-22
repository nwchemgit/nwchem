      subroutine tce_mrcc_calcexcit(erank,efrom,eto,espin,jr,ir,norb,
     1rtdb,iocc0)
c iphnum deleted

      implicit none
#include "tce_mrcc.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"

      integer norb,erank,efrom,eto,espin,jr,ir,i,j,k,is
      integer ndoexcit,iexfrom,iexto,iexspin,erankfrom,erankto,ioccwork
c      integer iphnum(maxorb,maxref,2)
c      logical mrccdebug
      integer rtdb
      logical nodezero
      integer iocc0(maxorb,2)
c      logical lusescffv

      dimension efrom(maxexcit), eto(maxexcit),espin(maxexcit)
      dimension norb(2),ioccwork(maxorb,2)

      erank=0
      if(ir.eq.jr) return 

      nodezero = (ga_nodeid().eq.0) 
c       if (.not.rtdb_get(rtdb,'mrcc:debug',mt_log,1,mrccdebug))
c     1 mrccdebug = .false.
c       if (.not.rtdb_get(rtdb,'mrcc:usescffermiv',mt_log,1,lusescffv))
c     1 lusescffv = .false.
c
c -------------------
c  lusescffv=true
c -------------------
c
      if(lusescffv.and.(ir.eq.0)) then

        do is=1,2
         erankfrom=erank
         erankto=erank

      do i=1,norb(is)
      if(iocc(i,jr,is).ne.iocc0(i,is)) then
          if(iocc0(i,is) .eq. 1 ) then
              erankfrom=erankfrom+1
c      write(6,*)'EFROM',iocc(i,jr,is),iocc0(i,is),is
              if(erankfrom.gt.maxexcit) then
                 write(6,*)'too high excitation found'
                 stop 'too high excitation found'
              endif
              efrom(erankfrom) = iphnum(i,maxref,is)
          else
              erankto=erankto+1
c      write(6,*)'ETO',iocc(i,jr,is),iocc0(i,is),is
              if(erankto.gt.maxexcit) then
                write(6,*)'too high excitation found'
                 stop 'too high excitation found'
              endif
              eto(erankto) = iphnum(i,maxref,is)
          endif
      endif
      enddo

      if(erankfrom.ne.erankto) then
         write(6,*) 'inconsistent number of electrons in ref. ',jr
         stop 'inconsistent number of electrons in references'
      endif

      if(erankfrom-erank.gt.1) then
         call tce_mrcc_sortint(erankfrom-erank,efrom(erank+1))
         call tce_mrcc_sortint(erankfrom-erank,eto(erank+1))
      endif
      do i=erank+1,erankfrom
      espin(i)=is
      enddo
      erank=erankfrom
      enddo

       return
      endif
c
c -------------------
c Occ for actual ref
c -------------------
c
      do is=1,2
      do i=1,norb(is)
        ioccwork(i,is)=iocc(i,jr,is)
      enddo
      enddo
c
c --------------------
c Search excit. level
c --------------------
c
      do is=1,2
      erankfrom=erank
      erankto=erank

      do i=1,norb(is)
      if(iocc(i,ir,is).ne.ioccwork(i,is)) then
          if(iocc(i,ir,is) .eq. 1 ) then
              erankfrom=erankfrom+1
              if(erankfrom.gt.maxexcit) then
                 write(6,*)'too high excitation found'
                 stop 'too high excitation found'
              endif
              efrom(erankfrom) = iphnum(i,ir,is)
          else
              erankto=erankto+1
              if(erankto.gt.maxexcit) then
                write(6,*)'too high excitation found'
                 stop 'too high excitation found'
              endif
              eto(erankto) = iphnum(i,ir,is)
          endif
      endif
      enddo

      if(erankfrom.ne.erankto) then
         write(6,*) 'inconsistent number of electrons in refs. ',ir,jr
         stop 'inconsistent number of electrons in references'
      endif
cjp sort indices; the order must correspond to the amplitude storage order,
cjp at least what the overall sign concerns in order to get correct sign
cjp of the rs couplings
      if(erankfrom-erank.gt.1) then
         call tce_mrcc_sortint(erankfrom-erank,efrom(erank+1))
         call tce_mrcc_sortint(erankfrom-erank,eto(erank+1))
      endif
      do i=erank+1,erankfrom
      espin(i)=is
      enddo
      erank=erankfrom
      enddo

      return 
      end

c $Id$
