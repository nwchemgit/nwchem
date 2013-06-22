      subroutine tce_mrcc_fillindexarray(rtdb)
        implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"

      integer rtdb
      integer iref
      integer ib
      integer ir
      integer orbspin
      integer orbindex
c      
c     ======================
c     Initialize parameters
c     ======================
c
c       if (.not.rtdb_get(rtdb,'mrcc:debug',mt_log,1,mrccdebug))
c     1 mrccdebug = .false.
c      if (.not.rtdb_get(rtdb,'bwcc:nref', mt_int, 1, nref))
c     1  call ERRQUIT('tce_mrcc_fillindexarray',1,RTDB_ERR)
c      
c     ======================
c           Main loop
c     ======================
c
      do iref=1,nref

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)
         k_active = k_active_tmpm(iref)

         do ib=1,noab+nvab
c            write(LuOut,*)ib,'Offset:',int_mb(k_offset+ib-1),'Range:',
c     1 int_mb(k_range+ib-1)
            do ir=1,int_mb(k_range+ib-1)
               orbspin = int_mb(k_spin+ib-1)-1
               orbindex = (1-orbspin+int_mb(k_mo_indexm(iref)+
     1 int_mb(k_offset+ib-1)+ir-1))/2
               orbinblck(orbindex,orbspin+1,iref) = ib
               offsetinblck(orbindex,orbspin+1,iref) =  ir - 1
c     1 int_mb(k_offset+ib-1) + ir - 1
c               orbfromblck(ib,ir-1,iref) = orbindex
            enddo
         enddo
       enddo

      return
      end


c $Id$
