      subroutine dft_mulwrap(me,g_dens,g_s)
c
c     wrapper to mul_pop
c
      implicit none
#include "cdft.fh"
      integer me,g_dens(2),g_s ! [in]
c

      if (me.eq.0)
     &     call dft_header
     &     (' Total Density - Mulliken Population Analysis')
      call mull_pop(geom,ao_bas_han,g_dens(1),g_s, 'total')
      if (ipol.eq.2)then
c     
c     analysis of spin density
c     
         if (me.eq.0)call dft_header
     &        (' Spin Density - Mulliken Population Analysis')
         call ga_dadd(1d0,g_dens(1),-2.d0,g_dens(2),g_dens(2))
         call mull_pop(geom,ao_bas_han,g_dens(2),g_s,'spin')
         call ga_dadd(1d0,g_dens(1),-1.d0,g_dens(2),g_dens(2))
         call ga_dscal(g_dens(2),0.5d0)
      endif
      return
      end
c
      subroutine dft_prevals(when,me,oprint_eval,oprint_vecs,
     ,  g_movecs,evals,occ)
      implicit none
#include "util.fh"
#include "cdft.fh"
#include "stdio.fh"
      logical oprint_eval,oprint_vecs
      character*5 when
      integer me
      integer g_movecs(2)
      double precision evals(*),occ(*)
c
      if(me.eq.0.and.oprint_eval) then
         if(when.eq.'inter') then
            call util_print_centered(LuOut,'alpha eigenvalues',
     &           20,.true.)
         else
            call util_print_centered(LuOut,'Final alpha eigenvalues',
     &           20,.false.)
         endif
         call util_print_centered(LuOut,
     &        '(all occupied plus 10 virtual)',20,.true.)
         call output(evals, 1, min(noc(1)+10,nbf_ao),
     &        1, 1, nbf_ao, 1, 1)
         call util_print_centered(LuOut,'alpha foccs',
     &        20,.true.)
         call output(occ, 1, min(noc(1)+10,nbf_ao),
     &        1, 1, nbf_ao, 1, 1)
         if (ipol.eq.2)then
            write(LuOut,*)
            call util_print_centered(LuOut,'beta eigenvalues',
     &           20,.true.)
            call output(evals(1+nbf_ao), 1, min(noc(1)+10,nbf_ao),
     &           1, 1, nbf_ao, 1, 1)
            call util_print_centered(LuOut,'beta foccs',
     &           20,.true.)
            call output(occ(nbf_ao+1),1,min(noc(1)+10,nbf_ao),
     &           1, 1, nbf_ao, 1, 1)
         endif
         call util_flush(LuOut)
      endif
      if (oprint_vecs)then
         if (me .eq. 0)then
            write(LuOut,*)
            if(when.eq.'inter') then
               call util_print_centered(LuOut,
     &              'Intermediate MO vectors',40,.true.)
            else
               call util_print_centered(
     &              LuOut,'Final MO vectors',40,.true.)
            endif
            write(LuOut,*)
            call util_flush(LuOut)
         endif
         call ga_print(g_movecs)
         if (ipol.eq.2)call ga_print(g_movecs(2))
      endif
      return 
      end
c
      subroutine dft_prevals2(me,g_movecs,evals,occ)
      implicit none
#include "util.fh"
#include "cdft.fh"
#include "stdio.fh"
      logical oprint_eval,oprint_vecs
      character*5 when
      integer me
      integer g_movecs(2)
      double precision evals(*),occ(*)
c
      if(me.eq.0) then
            call util_print_centered(LuOut,'Alpha eigenvalues',
     &           20,.false.)
         call util_print_centered(LuOut,
     &        '(all occupied plus 10 virtual)',20,.true.)
         call output(evals, 1,noc(1)+nbf_ao,
     &        1, 1, nbf_ao, 1, 1)
         call util_print_centered(LuOut,'Alpha occupations',
     &        20,.true.)
         call output(occ, 1, noc(1)+nbf_ao,nbf_ao,
     &        1, 1, nbf_ao, 1, 1)
         if (ipol.eq.2)then
            write(LuOut,*)
            call util_print_centered(LuOut,'Beta eigenvalues',
     &           20,.true.)
            call output(evals(1+nbf_ao), 1, noc(1)+nbf_ao,
     &           1, 1, nbf_ao, 1, 1)
            call util_print_centered(LuOut,'Beta occupations',
     &           20,.true.)
            call output(occ(nbf_ao+1),1,noc(1)+nbf_ao,
     &           1, 1, nbf_ao, 1, 1)
         endif
         call util_flush(LuOut)
      endif
c
      return 
      end
c $Id$
c
c     check if fractional occupation is present
      logical function has_frac_occ(rtdb)
c
      implicit none
c
      integer rtdb
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "msgids.fh"
c
      logical frac_occ
      integer ifocc
      logical fon
c
c     checks instances where fractional occupation is used
c
      has_frac_occ = .false. ! initialization
c
      if (.not.rtdb_get(rtdb, 'dft:fractional_occup', mt_log, 1,
     & frac_occ )) frac_occ=.false.
c
      if (.not.rtdb_get(rtdb,'focc:occup_switch',mt_int,1,ifocc))
     &             ifocc = 0
c
      if (.not.rtdb_get(rtdb, 'dft:fon', mt_log, 1, fon))
     &             fon = .false.
c
      if (fon.or.(ifocc.eq.1).or.frac_occ) has_frac_occ = .true.
c
c     write(luout,*) "fon: ",fon
c     write(luout,*) "ifocc: ",ifocc
c     write(luout,*) "frac_occ: ",frac_occ
c     write(luout,*) "has_frac_occ: ",has_frac_occ
c
      return
      end
c
c     Write out the alpha/beta split eigenvalues
c
      subroutine dft_print_evals_ab(nbf_ao, occ, evals,ipol)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "util.fh"
c
      integer nbf_ao            ! [input] No. of basis functions
      integer ipol ! [in]
      double precision occ(2*nbf_ao)   ! [input] Occupation numbers 
      double precision evals(2*nbf_ao) ! [input] Eigenvalues
c
      integer i,ok
      character*255 evalfile
      integer uniteval
c
      logical util_io_unit
      external util_io_unit
c
      call util_file_name('evals',.false.,.false.,evalfile)
      if(.not.util_io_unit(80,90,uniteval))
     &  call errquit('print_evals_ab:cannot get free unit', 0,0)
c
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
       open(uniteval, status='unknown', form='formatted',file=evalfile)
       write(uniteval,*) nbf_ao
       write(uniteval,*)
       if(ipol.eq.1) then
          do i = 1,nbf_ao
             write(uniteval,*) i, evals(i)
          end do
       else
          do i = 1,nbf_ao
             write(uniteval,*) i, evals(i), evals(nbf_ao+i)
          end do
       endif
       close(uniteval)
       ok = 1
      endif ! ga_nodeid() = 0
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok,
     C     MA_sizeof(MT_INT,1,MT_BYTE), 0) ! Propagate status
c
      call ga_sync()
c
      return
      end
