      block data block_int1e
      implicit none
#include "cint1cache.fh"
      data ocache /.false./
      end
c
C> \brief Compute 1-electron integrals and add them to a global array
C>
C> This routine computes 1-electron integrals of the specified kind and
C> adds them to a global array. The routine always computes a
C> rectangular matrix. It is even capable of using two different 
C> basis sets. Different kinds of integrals can be computed, they are:
C>
C> * "kinetic": the kinetic energy integrals
C>
C> * "potential": the nuclear attraction
C>
C> * "overlap": the overlap integrals
C>
C> * "pvp": ???
C>
C> * "so": spin-orbit integrals
C>
C> * "cos_chg_pot": the integrals with the COSMO charges
C>
C> * "bq_pot": the integrals with BQ-charges
C>
cc AJL/Begin
C> * "potential_beta": the nuclear attraction with ECPs for the
C>                     beta channel
cc AJL/End
C>
C> The routine can exploit symmetry. If oskel is true then the 
C> "petite-list" symmetry will be used, although this requires both
C> basis sets to be the same.
C
      subroutine int_1e_ga(ibas, jbas, g, integ_type, oskel)
C$Id$
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
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
      integer ibas             !< [Input] bra basis sets
      integer jbas             !< [Input] ket basis sets
      integer g(*)             !< [Output] GA handle to array 
      character*(*) integ_type !< [Input] Name of integrals to compute
      logical oskel            !< [Input] If true generate symmetry unique list
c
c     local variables
c
      integer type
      logical dobq
      character*255 integ_type1
c
      call ga_sync()
c
      dobq = geom_extbq_on()
      integ_type1 = integ_type
c
      if (inp_compare(.false., integ_type1, 'potential0')) then
        integ_type1='potential'
        dobq=.false.
cc AJL/Begin/SPIN-POLARISED ECPs
      elseif (inp_compare(.false., integ_type1, 'potential_beta0')) then
        integ_type1='potential_beta'
        dobq=.false.
cc AJL/End
      end if
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: use of symmetry requires ibas=jbas', ibas,
     &       BASIS_ERR)
      end if
c     
      if (inp_compare(.false., integ_type1, 'kinetic')) then
         type = 1
      else if (inp_compare(.false., integ_type1, 'potential')) then
         type = 2
      else if (inp_compare(.false., integ_type1, 'overlap')) then
         type = 3
      else if (inp_compare(.false., integ_type1, 'pvp')) then
         type = 4
      else if (inp_compare(.false., integ_type1, 'so'))then 
         type = 5 
      else if (inp_compare(.false., integ_type1, 'cos_chg_pot'))then 
         type = 6 
      else if (inp_compare(.false., integ_type1, 'bq_pot'))then
         type=7
         if(.not.dobq) return
cc AJL/Begin/SPIN-POLARISED ECPs
      else if (inp_compare(.false., integ_type1, 'potential_beta')) then
         type = 8
cc AJL/End
      else 
         write(6,*) ' integ_type1 = ', integ_type1
         call errquit('int_1e_ga: unknown integ_type', 0, INT_ERR)
      end if
       if(type.lt.0.or.type.gt.10) call
     E errquit(' bogus int1 type',0,INT_ERR)
c
c     Check if Douglas-Kroll is required, if so:
c     getting the Douglas-Kroll kinetic and potential energy 
c     integrals
c
cc AJL/Begin/SPIN-POLARISED ECPs
c      if (doug_kroll .and. (type.le.2)) then
      if (doug_kroll .and. ((type.le.2) .or. (type.eq.8))) then
cc AJL/End
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: Douglas-Kroll requires ibas=jbas', ibas,
     &       BASIS_ERR)
         call int_1edk_ga(ibas,g,type,oskel)
         call ga_sync()            ! So that no nasty races can result
         return
      endif
c
c     Doing Douglas-Kroll Spin-Orbit terms
c
      if (doug_kroll .and. (type.eq.5)) then
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: Douglas-Kroll requires ibas=jbas', ibas,
     &       BASIS_ERR)
         call int_1edk_so_ga(ibas,g,type,oskel)
         call ga_sync()            ! So that no nasty races can result
         return
      endif
c
      if (dobq) then
         call int_1e_ooldga(ibas, jbas, g, integ_type1, oskel)
cc AJL/Begin/SPIN-POLARISED ECPs
cc         if(inp_compare(.false., integ_type1, 'potential'))
cc     &      call int_1e_ooldga(ibas, jbas, g, 'bq_pot', oskel)
         if (type.eq.2.or.type.eq.8) ! potential or potential_beta
     &      call int_1e_ooldga(ibas, jbas, g, 'bq_pot', oskel)
cc AJL/End
      else
         call int_1e_oldga(ibas, jbas, g, integ_type1, oskel)
      end if
c
      end
c
      subroutine int_1e_ooldga(ibas, jbas, g, integ_type, oskel)
      implicit none
#include "errquit.fh"
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
      integer g(3)              ! [output] GA handle to array. g(1..3) are used only
                                ! for spin-orbit calculations. All other calculations use
                                ! g(1) only 
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
      integer i, noffset,g_loc
c
      logical ocache_save
c
      logical odbug
      logical osome
c
      odbug=.false.
      osome=.false.
      osome=osome.or.odbug
      odbug=odbug.and.(ga_nodeid().eq.0)
      osome=osome.and.(ga_nodeid().eq.0)
      if(osome) then
         write(6,*) 'in -int_1e_ooldga- ... integ_type = ',
     $              integ_type,ga_nodeid()
         call util_flush(6)
      endif
c
      call ga_sync()
      if (oscfps) call pstat_on(ps_int_1e)
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: use of symmetry requires ibas=jbas', ibas,
     &       BASIS_ERR)
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
      else if (inp_compare(.false., integ_type, 'bq_pot'))then 
         type = 7
cc AJL/Begin/SPIN ECPs 
      else if (inp_compare(.false., integ_type, 'potential_beta')) then
         type = 8
cc AJL/End
      else
         write(6,*) ' integ_type = ', integ_type,ga_nodeid()
         call errquit('int_1e_ga: unknown integ_type', 0, INT_ERR)
      end if
c
c     ----- save ocache logical variable -----
c
      if(type.eq.6 .or. type.eq.7 ) then
         ocache_save=ocache
         ocache     =.false.
      endif
c
c     In-core caching
c
c      ocache = .false.
      if(osome) then
         write(6,*) 'ocache = ',ocache,ga_nodeid()
         call util_flush(6)
      endif
      if (ocache .and. (ibas.eq.bas1) .and. (jbas.eq.bas1)) then
*         if (ga_nodeid().eq.0) then
*           call ga_summarize(1)
*           call util_flush(6)
*         endif
*         call ga_sync()
c
cc AJL/Begin/SPIN-POLARISED ECPs
         if (type.ne.8) then
           call ga_dadd(1.0d0, g_cache(type), 1.0d0, g, g)
         else ! accomodate for potential_beta, stored in g_cache(4)
           call ga_dadd(1.0d0, g_cache(4), 1.0d0, g, g)
         endif
cc AJL/End
c
         if (oscfps) call pstat_off(ps_int_1e)
	 return
      endif
c     
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for ibas', ibas,
     &       BASIS_ERR)
      if (.not. bas_numcont(jbas, nshell_j)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for jbas', jbas,
     &       BASIS_ERR)
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
     $     call errquit('int_1e_ga: ma failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'int_1e_ga:scr',l_scr,k_scr))
     $     call errquit('int_1e_ga: ma failed', mem1, MA_ERR)
c     
c     Loop thru shells with static parallel work decomposition
c     
      if (.not.inp_compare(.false., integ_type, 'so')) then
         if(.not.ga_duplicate(g,g_loc,'local g')) call
     .        errquit('int1e: dupl failed',0, GA_ERR)
         call ga_zero(g_loc)
      endif
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
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
                  if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', jbas,
     &       BASIS_ERR)
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
                     if(odbug) then
                        write(6,*) 'type = 6 ... potential ...',
     $                             ga_nodeid()
                        call util_flush(6)
                     endif
                     call int_1epot (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 7) then 
                     call int_1epot1 (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
cc AJL/Begin/SPIN ECPs
                 else if (type .eq. 8) then
                     call int_1epe_beta (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
cc AJL/End
                  else
                     call errquit('int_1e_ga: invalid type?', type,
     &       GA_ERR)
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
                     if(odbug) then
                        write(6,*) 'ga_acc = ...',ga_nodeid()
                        call util_flush(6)
                     endif
cedo                     call ga_acc(g, ilo, ihi, jlo, jhi, dbl_mb(k_buf),
cedo     $                    idim, q2) 
                    if(q2.ne.1) then
                      call dscal((ihi-ilo+1)*(jhi-jlo+1),
     ,            q2,dbl_mb(k_buf),1)
                    endif
                     call ga_put(g_loc, ilo, ihi, jlo, jhi,
     ,             dbl_mb(k_buf),  idim) 
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
     $     ('int_1e_ga: pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_1e_ga: pop failed', 0, GA_ERR)
c
c     ----- restore ocache -----
c
      if(type.eq.6 .or. type.eq.7) then
         ocache=ocache_save
      endif
c     
      if (.not.inp_compare(.false., integ_type, 'so')) then
         call ga_dadd(1.0d0, g_loc, 1.0d0, g, g)
         if (.not. ga_destroy(g_loc)) call errquit('i1ega: gad?',0,
     &       GA_ERR)
      endif
      call ga_sync()            ! So that no nasty races can result
c     
      if (oscfps) call pstat_off(ps_int_1e)
c     
      end
      subroutine int_1e_cache_ga(basis, oskel)
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
#include "bas.fh"
      integer basis
      logical oskel
c
      integer geom
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
cc AJL/Begin/SPIN ECPs
#include "global.fh"
      integer ecp_handle
cc AJL/End
c
      if (ocache) call int_1e_uncache_ga()
c
      if (.not. bas_geom(basis, geom)) call errquit
     $   ('int_1e_cache_ga: basis corrupt?', 0, BASIS_ERR)
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
      call int_1e_ga(basis, basis, g_cache(2), 'potential0', oskel)
      call int_1e_ga(basis, basis, g_cache(3), 'overlap', .false.)
c
cc AJL/Begin/SPIN ECPs
      if (bas_get_ecp_handle(basis,ecp_handle)) then
        if (.not.ecp_get_high_chan(ecp_handle,channels)) channels = 1
        if (channels.gt.1) then
          if (ga_nodeid().eq.0) then
            write(6,*) 'ECP Channels : ', channels
            call util_flush(6)
          endif
          g_cache(4) = 
     &            ga_create_atom_blocked(geom, basis,'int1e: v_beta')
          call ga_zero(g_cache(4))
          call int_1e_ga(basis, basis, g_cache(4),
     &                   'potential_beta0', oskel)
        end if
      end if
cc AJL/End
c
      bas1 = basis
      ocache = .true.
c
      end
      subroutine int_1e_uncache_ga()
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
      logical ga_destroy
c
*      write(6,*) ' Disabling caching of 1-e integrals '
c
      if (.not. ocache) return
      if (.not. ga_destroy(g_cache(1))) call errquit('i1ega: gad?',0,
     &       BASIS_ERR)
      if (.not. ga_destroy(g_cache(2))) call errquit('i1ega: gad?',0,
     &       BASIS_ERR)
      if (.not. ga_destroy(g_cache(3))) call errquit('i1ega: gad?',0,
     &       BASIS_ERR)
c
cc AJL/Begin/SPIN ECPs
      if (channels.gt.1) then
        if (.not. ga_destroy(g_cache(4))) call errquit('i1ega: gad?',0,
     &                       BASIS_ERR)
      end if
cc AJL/End
c
      bas1 = -1
      ocache = .false.
      end
      subroutine int_1e_oldga(ibas, jbas, g, integ_type, oskel)
      implicit none
#include "inp.fh"
#include "errquit.fh"
#include "global.fh"
      integer ibas, jbas        ! [input] bra and ket basis sets
      integer g(3)              ! [output] GA handle to array. g(1..3) are used only
                                ! for spin-orbit calculations. All other calculations use
                                ! g(1) only 
      character*(*) integ_type  ! [input] Name of integrals to compute
      logical oskel             ! [input] If true generate symmetry unique list
c
      external int_1eke,int_1epe,int_1eov,int_1epvpe,
     I     int_1epot,int_1epot1,
cc AJL/Begin/SPIN ECPs
     I     int_1epe_beta
cc AJL/End
      integer type
c

      if(ibas.ne.jbas) then
         call int_1e_ooldga(ibas, jbas, g,integ_type,oskel)
         return
      endif
      if (inp_compare(.false., integ_type, 'so'))then 
cso
          call int_1e_ooldga(ibas, jbas, g,integ_type,oskel)
          return
       endif
      if (inp_compare(.false., integ_type, 'kinetic')) then
         type = 1
         call int_1e_oldga0(ibas, g, type,oskel,int_1eke)
      else if (inp_compare(.false., integ_type, 'potential')) then
         type = 2
         call int_1e_oldga0(ibas, g, type,oskel,int_1epe)
      else if (inp_compare(.false., integ_type, 'overlap')) then
         type = 3
         call int_1e_oldga0(ibas, g, type,oskel,int_1eov)
      else if (inp_compare(.false., integ_type, 'pvp')) then
         type = 4
         call int_1e_oldga0(ibas, g, type,oskel,int_1epvpe)
      else if (inp_compare(.false., integ_type, 'cos_chg_pot'))then 
         type = 6 
         call int_1e_oldga0(ibas, g, type,oskel,int_1epot)
      else if (inp_compare(.false., integ_type, 'bq_pot'))then 
         type = 7 
         call int_1e_oldga0(ibas, g, type,oskel,int_1epot1)
cc AJL/Begin/SPIN ECPs
      else if (inp_compare(.false., integ_type, 'potential_beta')) then
         type = 8
         call int_1e_oldga0(ibas, g, type, oskel, int_1epe_beta)
cc AJL/End/SPIN ECPs
      else
         write(6,*) ' integ_type = ', integ_type,ga_nodeid()
         call errquit('int_1e_oldga: unknown integ_type', 0, INT_ERR)
      end if

      return
      end
      subroutine int_1e_oldga0(ibas, g, type, oskel,
     I     int_call)
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
c     
c     This is the original routine of int_1e_ga. It is still needed by the 
c     Douglas-Kroll routines to avoid recursive routine calling.
c
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     and ADD them into the given global array.
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... requires that ibas = jbas.
c     
c     arguments
c     
      integer ibas              ! [input] bra and ket basis sets
      integer g(3)              ! [output] GA handle to array. g(1..3) are used only
                                ! for spin-orbit calculations. All other calculations use
                                ! g(1) only 
      integer type              ! [input] Name of integrals to compute
      logical oskel             ! [input] If true generate symmetry unique list
c     
c     local variables
c     
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      logical odoit
      double precision q2
      external block_int1e      ! For T3D
      integer i, g_loc
      integer geom
c
      logical ocache_save
c
      logical odbug
      logical osome
      external int_call
c
      odbug=.false.
      osome=.false.
      osome=osome.or.odbug
      odbug=odbug.and.(ga_nodeid().eq.0)
      osome=osome.and.(ga_nodeid().eq.0)
      if(osome) then
         write(6,*) 'in -int_1e_oldga0- ... integ_type = ',
     $              type,ga_nodeid()
         call util_flush(6)
      endif
c
      call ga_sync()
c
      if (oscfps) call pstat_on(ps_int_1e)
c     
c     ----- save ocache logical variable -----
c
      if(type.eq.6 .or. type.eq.7 ) then
         ocache_save=ocache
         ocache     =.false.
      endif
c
c     In-core caching
c
c      ocache = .false.
      if(osome) then
         write(6,*) 'ocache = ',ocache,ga_nodeid()
         call util_flush(6)
      endif
      if (ocache .and. (ibas.eq.bas1)) then
cc AJL/Begin/SPIN ECPs
         if (type.ne.8) then ! everything except potential_beta
           call ga_dadd(1.0d0, g_cache(type), 1.0d0, g, g)
         else
           call ga_dadd(1.0d0, g_cache(4), 1.0d0, g, g)
         end if
cc AJL/End
         if (oscfps) call pstat_off(ps_int_1e)
         return
      endif
c     
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('int_1e_oldga0: bas_numcont failed for ibas', ibas,
     &       BASIS_ERR)
      nshell_j=nshell_i
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
c     
c     Loop thru shells with static parallel work decomposition
c     
      if(.not.ga_duplicate(g,g_loc,'local g')) call
     .     errquit('int_1e_oldga0: dupl failed',0, GA_ERR)
      call ga_zero(g_loc)
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      if (.not. bas_geom(ibas, geom)) call errquit
     $   ('int_1e_oldga0: basis corrupt?', 0, BASIS_ERR)
c
      call int_mem_1e(max1e, mem1)
      call ga_get2eri(ibas, g_loc, oskel, max1e, mem1, int_call)
      if (oskel) call sym_symmetrize(geom,ibas,.false.,g_loc)
c
c     ----- restore ocache -----
c
      if(type.eq.6 .or. type.eq.7) then
         ocache=ocache_save
      endif
c
      call ga_dadd(1.0d0, g_loc, 1.0d0, g, g)
      call ga_sync()            ! So that no nasty races can result
      if (.not. ga_destroy(g_loc)) call errquit('int_1e_oldga0: gad?',0,
     &     GA_ERR)
      call ga_sync()            ! So that no nasty races can result
c     
      if (oscfps) call pstat_off(ps_int_1e)
c     
      end
