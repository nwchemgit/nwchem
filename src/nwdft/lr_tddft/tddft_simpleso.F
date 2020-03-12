c
c     ----------------------------------------
c     Get the spin-orbit corrected eigenvalues
c     ----------------------------------------
c
      subroutine tddft_simpleso(rtdb,ipol,nbf_ao)
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
c
      integer i,j,ok
      logical lsimpleso
      character*255 fn_simpleso
      logical status
      integer fn
      integer nmo,imo
      double precision evlso(nbf_ao,2)
c
      logical util_io_unit
      external util_io_unit
c
      call ga_sync()
      ok = 0
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
      call util_file_name_resolve(fn_simpleso, .false.)
c
      if (ga_nodeid().eq.0) then
c
c      initialize
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
c       if (nmo.ne.nbf_ao)  
c     &  call errquit('tddft_simpleso:nmo is not equal to nbf_ao', 0,0)
       read(fn,*)
       do i = 1,nmo
          if (ipol.eq.1) then
           read(fn,*) imo,evlso(imo,1)
          else 
           read(fn,*) imo,evlso(imo,1),evlso(imo,2)
          end if
       end do
       close(fn)
       ok = 1
      end if  ! ga_nodeid().eq.0
c
c     Broadcast status to other nodes
      call ga_brdcst(Msg_Vec_Stat, ok, ma_sizeof(mt_int,1,mt_byte),0)
      call ga_sync()

      if(.not.rtdb_put(rtdb,'simpleso:evlso',mt_dbl,2*nbf_ao,evlso))
     $     call errquit('tddft_simpleso: rtdb put evlso  ',914,
     &         rtdb_err)
c
 9982 format(' Reading spin-orbit corrected eigenvalues: ',a)
c
      return
      end
