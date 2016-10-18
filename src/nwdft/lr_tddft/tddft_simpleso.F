c
c     ----------------------------------------
c     Get the spin-orbit corrected eigenvalues
c     ----------------------------------------
c
      subroutine tddft_simpleso(rtdb,ipol,nbf_ao,evl)
c
      implicit none
c
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "errquit.fh"
#include "inp.fh"
c
      integer rtdb
      integer ipol
      integer nbf_ao 
      double precision evl(nbf_ao,2)
c
      integer i,j
      logical lsimpleso
      character*255 fn_simpleso
      logical status
      integer fn
      integer nmo,imo
      double precision evlso(nbf_ao,2)
c
      integer ddblsize
c
      logical util_io_unit
      external util_io_unit
c
c     start processing simple so file
c     -------------------------------
c
      if (.not.rtdb_cget(rtdb,'tddft:simpleso',1,fn_simpleso))
     &  call errquit('tddft_simpleso:failed to get tddft:simpleso',0, 
     &   RTDB_ERR)
c
      if(.not.util_io_unit(80,90,fn))
     &  call errquit('tddft_simpleso:cannot get free unit', 0,0)
c
      if (ga_nodeid().eq.0) then
c
c       initialize
        do i = 1,ipol
          do j = 1,nbf_ao
            evlso(j,i) = 0.d0
          end do ! j
        end do ! i
c
        open(unit=fn,status="unknown",form="formatted",file=fn_simpleso)
        write(luout,*)
        write(luout,9982) fn_simpleso(1:inp_strlen(fn_simpleso))
        write(luout,*)
        read(fn,*) nmo
        do i = 1,nmo
            read(fn,*) imo,evlso(imo,1)
        end do
        read(fn,*) nmo
        do i = 1,nmo
            read(fn,*) imo,evlso(imo,2)
        end do
        close(fn)
c
c       Assign spin-orbit corrected eigenvalues
        do i = 1,ipol
          do j = 1,nbf_ao
            if (dabs(evlso(j,i)).gt.0.d0) evl(j,i) = evlso(j,i)
          end do ! j
        end do ! i
c
      end if ! ga_nodeid()
c
c     Broadcast status to other nodes
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
      call ga_brdcst(Msg_Vec_Stat,evl,nbf_ao*2*ddblsize,0)
      call ga_sync()
c
 9982 format(' Reading spin-orbit corrected eigenvalues: ',a)
c
      return
      end
