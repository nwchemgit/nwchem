      block data block_int1e
      implicit none
#include "cint1cache.fh"
      data ocache /.false./
      end
      subroutine int_1e_ga(ibas, jbas, g, integ_type, oskel)
C$Id: int_1e_ga.F,v 1.26 2000-05-20 22:21:55 mdupuis Exp $
      implicit none
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
c     
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     and ADD them into the given global array.
c     This version computes the full square of integrals and should work
c     OK even if ibas != jbas.
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... requires that ibas = jbas.
c     
c     arguments
c     
      integer ibas, jbas        ! [input] bra and ket basis sets
      integer g(*)              ! [output] GA handle to array 
      character*(*) integ_type  ! [input] Name of integrals to compute
      logical oskel             ! [input] If true generate symmetry unique list
c
c     local variables
c
      integer type
c
      call ga_sync()
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: use of symmetry requires ibas=jbas', ibas)
      end if
c     
      if (inp_compare(.false., integ_type, 'kinetic')) then
         type = 1
      else if (inp_compare(.false., integ_type, 'potential')) then
         type = 2
      else if (inp_compare(.false., integ_type, 'overlap')) then
         type = 3
      else if (inp_compare(.false., integ_type, 'pvp')) then
         type = 4
      else if (inp_compare(.false., integ_type, 'so'))then 
         type = 5 
      else if (inp_compare(.false., integ_type, 'cos_chg_pot'))then 
         type = 6 
      else
         write(6,*) ' integ_type = ', integ_type
         call errquit('int_1e_ga: unknown integ_type', 0)
      end if
c
c     Check if Douglas-Kroll is required, if so:
c     getting the Douglas-Kroll kinetic and potential energy 
c     integrals
c
      if (doug_kroll .and. (type.le.2)) then
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: Douglas-Kroll requires ibas=jbas', ibas)
         call int_1edk_ga(ibas,g,type,oskel)
         call ga_sync()            ! So that no nasty races can result
         return
      endif
c
      call int_1e_oldga(ibas, jbas, g, integ_type, oskel)
c
      end
c
      subroutine int_1e_oldga(ibas, jbas, g, integ_type, oskel)
      implicit none
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
c     
c     This is the original routine of int_1e_ga. It is still needed by the 
c     Douglas-Kroll routines to avoid recursive routine calling.
c
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     and ADD them into the given global array.
c     This version computes the full square of integrals and should work
c     OK even if ibas != jbas.
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... requires that ibas = jbas.
c     
c     arguments
c     
      integer ibas, jbas        ! [input] bra and ket basis sets
      integer g(3)              ! [output] GA handle to array 
      character*(*) integ_type  ! [input] Name of integrals to compute
      logical oskel             ! [input] If true generate symmetry unique list
c     
c     local variables
c     
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
      logical odoit
      double precision q2
      external block_int1e      ! For T3D
      integer i, noffset
c
      call ga_sync()
      if (oscfps) call pstat_on(ps_int_1e)
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: use of symmetry requires ibas=jbas', ibas)
      end if
c     
      if (inp_compare(.false., integ_type, 'kinetic')) then
         type = 1
      else if (inp_compare(.false., integ_type, 'potential')) then
         type = 2
      else if (inp_compare(.false., integ_type, 'overlap')) then
         type = 3
      else if (inp_compare(.false., integ_type, 'pvp')) then
         type = 4
      else if (inp_compare(.false., integ_type, 'so'))then 
         type = 5 
      else if (inp_compare(.false., integ_type, 'cos_chg_pot'))then 
         type = 6 
      else
         write(6,*) ' integ_type = ', integ_type
         call errquit('int_1e_ga: unknown integ_type', 0)
      end if
c
c     In-core caching
c
c      ocache = .false.
      if (ocache .and. (ibas.eq.bas1) .and. (jbas.eq.bas1)) then
*         if (ga_nodeid().eq.0) then
*           call ga_summarize(1)
*           call util_flush(6)
*         endif
*         call ga_sync()
         call ga_dadd(1.0d0, g_cache(type), 1.0d0, g, g)
         if (oscfps) call pstat_off(ps_int_1e)
	 return
      endif
c     
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for ibas', ibas)
      if (.not. bas_numcont(jbas, nshell_j)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for jbas', jbas)
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c     
c     k_* are the offsets corrsponding to the l_* handles
c     
cso
      call int_mem_1e(max1e, mem1)
      if (type .eq. 5) then  
         max1e = 3*max1e 
         mem1 = 3*mem1 
      endif 
cso
      if (.not. MA_push_get(MT_DBL,max1e,'int_1e_ga:buf',l_buf,k_buf))
     $     call errquit('int_1e_ga: ma failed', max1e)
      if (.not. MA_push_get(MT_DBL, mem1,'int_1e_ga:scr',l_scr,k_scr))
     $     call errquit('int_1e_ga: ma failed', mem1)
c     
c     Loop thru shells with static parallel work decomposition
c     
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      q2 = 1.0d0
      do jshell = 1, nshell_j
         do ishell = 1, nshell_i
c
            if (mod(ijshell, nproc) .eq. iproc) then
               odoit = .true.
               if (oskel) 
     $              odoit = sym_shell_pair(ibas, ishell, jshell, q2)
c
               if (odoit) then
                  if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) 
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', ibas)
                  if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', jbas)
                  idim = ihi - ilo + 1
c     
c     Generate the integrals
c     
                  if (type .eq. 1) then
                     call int_1eke (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 2) then
                     call int_1epe (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 3) then
                     call int_1eov (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 4) then
                     call int_1epvpe (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 5) then 
                     call intso_1e (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 6) then 
                     call int_1epe (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else
                     call errquit('int_1e_ga: invalid type?', type)
                  end if
c     
c     Add the integrals into the global array
c    
                  if (inp_compare(.false., integ_type, 'so')) then
                     do i = 1, 3
                        noffset = (ihi-ilo+1)*(jhi-jlo+1)*(i-1)   
                        call ga_acc(g(i), ilo, ihi, jlo, jhi, 
     $                    dbl_mb(k_buf+noffset),
     $                    idim, q2)
                     enddo 
                  else  
                     call ga_acc(g, ilo, ihi, jlo, jhi, dbl_mb(k_buf),
     $                    idim, q2) 
                  end if
               end if
            endif
            ijshell = ijshell + 1
         end do
      end do
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_1e_ga: pop failed', 0)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_1e_ga: pop failed', 0)
c     
      call ga_sync()            ! So that no nasty races can result
c     
      if (oscfps) call pstat_off(ps_int_1e)
c     
      end
      subroutine int_1e_cache_ga(basis, oskel)
      implicit none
#include "cint1cache.fh"
#include "bas.fh"
      integer basis
      logical oskel
c
      integer geom
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      if (ocache) call int_1e_uncache_ga()
c
      if (.not. bas_geom(basis, geom)) call errquit
     $   ('int_1e_cache_ga: basis corrupt?', 0)
c
c	The use of ga_dadd to do the copy requires all uses of
c	these integrals to be allocated via ga_create_atom_blocked
c
*      write(6,*) ' Enabling caching of 1-e integrals '
c
      g_cache(1) = ga_create_atom_blocked(geom, basis,'int1e: t')
      g_cache(2) = ga_create_atom_blocked(geom, basis,'int1e: v')
      g_cache(3) = ga_create_atom_blocked(geom, basis,'int1e: s')
      call ga_zero(g_cache(1))
      call ga_zero(g_cache(2))
      call ga_zero(g_cache(3))
      call int_1e_ga(basis, basis, g_cache(1), 'kinetic', oskel)
      call int_1e_ga(basis, basis, g_cache(2), 'potential', oskel)
      call int_1e_ga(basis, basis, g_cache(3), 'overlap', .false.)
      bas1 = basis
      ocache = .true.
c
      end
      subroutine int_1e_uncache_ga()
      implicit none
#include "cint1cache.fh"
      logical ga_destroy
c
*      write(6,*) ' Disabling caching of 1-e integrals '
c
      if (.not. ocache) return
      if (.not. ga_destroy(g_cache(1))) call errquit('i1ega: gad?',0)
      if (.not. ga_destroy(g_cache(2))) call errquit('i1ega: gad?',0)
      if (.not. ga_destroy(g_cache(3))) call errquit('i1ega: gad?',0)
      bas1 = -1
      ocache = .false.
      end
