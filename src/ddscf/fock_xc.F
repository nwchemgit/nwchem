C>
C> \brief Wrapper for AO-basis XC matrices without fitting
C>
      subroutine fock_xc(geom, nbf_ao, ao_bas_han,
     &   nfock, g_dens, g_xc, Exc, nExc, l3d)
c
c     $Id$
c
c     Wrapper routine for AO-basis XC matrices without fitting
c
c     BGJ - 8/98
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
c!!! BGJ
#include "bgj.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c!!! BGJ
c
      integer geom              !< [Input] Geometry handle
      integer nbf_ao,ao_bas_han !< [Input] No. of AOs and AO handle
      integer nfock             !< [Input] No. of XC matrices
c                               !<         = No. of density matx (RDFT)
c                               !<         = 4 * No. of density matx (UDFT)
      integer ndens             !< [Local] = ipol
      integer g_dens(*)         !< [Input] Array of handles to densities
      integer g_xc(*)           !< [Input] Array of handles to XC matrices
      integer calc_type         !< [Input] Type of XC matrix calculation
c                               !<         = 1 XC contrib to Fock
c                               !<         = 2 XC contrib to CPKS LHS
c                               !<         = 3 XC contrib to CPKS RHS
c                               !<         = 4 XC NMR contrib
c                               !<         = 5 XC third-derivative contrib
      logical l3d               !< [Local] Whether 3-D NGA density
c                               !<         and fock matrices are used.
c                               !<         Must be true when and only when
c                               !<         calctype = 2
      logical triplet           !< [Rtdb]  True if TDDFT triplet excitation
c                               !<         energy calculations from RDFT.
c                               !<         If not set, then .false. is
c                               !<         assumed.  For Hessian, this must
c                               !<         be false (or not set).
c
      integer nExc              !< [Input]  no. exchange-correlation energy
c                               !<          components.
      double precision Exc(nExc)!< [Output] Exchange-correlation energy
c                               !<          components
c
c     Local declarations
c
      integer rtdb
      integer g_wght_dum, g_xyz_dum, g_nq_dum,
     &        l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom
      integer g_dens_scf(2)
      logical wght_GA_dum,grid_reopen,xc_gotxc
      double precision rho_n_dum,tdum
      integer dims(3),chunk(3)
      integer igxcd,lgxcd
      integer xc_getipol
      logical oprint, debug
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer nfock_xc,i
      integer ityp, ndim
      logical lcgmin
      external xc_gotxc,xc_getipol
      integer ga_create_atom_blocked,dft_npatch
      external ga_create_atom_blocked
c
      integer nroots  !number of roots for the TDDFT gradient calculations
c
      if(.not.xc_gotxc()) return
      oprint= util_print('fock_xc',print_debug)

c ... jochen: 
      debug = .false.
c
c Note that ndens=ipol and is not necessarily
c equal to No. of density matrices.  "nfock" or "nmat"
c will be related to No. of density matrices. 
c
      ndens=1
      if(xc_getipol().eq.2) ndens=2
      rtdb = bgj_get_rtdb_handle()
      alo(2) = 1
      ahi(2) = nbf_ao
      alo(3) = 1
      ahi(3) = nbf_ao
      blo(2) = 1
      bhi(2) = nbf_ao
      blo(3) = 1
      bhi(3) = nbf_ao
c
c     Get fock_xc variables
c
      if (.not. rtdb_get(rtdb, 'fock_xc:calc_type', mt_int, 1,
     &   calc_type)) then
         if ((oprint.or.debug).and.ga_nodeid().eq.0)
     &      write(luout,*)' fockxc: calc_type not set: setting to 1'
         calc_type = 1
      endif
      if (calc_type .eq. 0) then
         if ((oprint.or.debug).and.ga_nodeid().eq.0)
     &      write(luout,*)' fockxc: calc_type not set: setting to 1'
         calc_type = 1
      endif
c
c     == need a better solution for this == NXG
c Daniel (2-16-13): This line of code is VERY dangerous.  It caused 
c difficult to locate problems in the TDDFT gradient routines.
      if (.not.rtdb_get(rtdb, 'dft:cgmin', mt_log, 1, lcgmin))
     &   lcgmin=.false.
c
c     Whether this is a triplet excitation energy calculation
c
      if (.not. rtdb_get(rtdb, 'fock_xc:triplet', mt_log, 1,
     &   triplet)) triplet=.false.
c
c     for uhf calculations, nmat is 4*nvec. The nmat passed to _quadv0_gen
c     should be the number of vectors 
c
      nfock_xc=nfock
c Daniel (2-13-13): The original way this part was written malfunctions 
c for third derivatives.  The routine, xc_rhogen, only cares about 
c whether the calculation is restricted or unrestricted (via ndens/ipol) 
c and how many Fock matrices we need to build (via nfock_xc).  We
c therefore need to divide nfock by ndens for 3rd derivatives, since 
c nfock is set to nroots*ipol in the TDDFT gradient routines.  This
c is only important for unrestricted calculations.
      if(ndens.eq.2)  then
        if (calc_type.eq.5) then
          nfock_xc=nfock/ndens
        else
          nfock_xc=nfock/4
        endif
      endif
      if ((debug.or.oprint).and.ga_nodeid().eq.0) then
         write(luout,*) '--------Entered fock_xc-------------'
         write(luout,*) ' calc_type =',calc_type
         write(luout,*) ' nfock ',nfock,' ndens ',ndens
         write(luout,*) ' nfock_xc ',nfock_xc
         write(luout,*) ' l3d ',l3d
         write(luout,*) ' triplet ',triplet
         call util_flush(luout)
       endif
       if(nfock_xc.eq.0) return
c
c     assumed l3d true AND calc_type=2
c
      if(l3d.and.calc_type.ne.2.or.
     &   (.not.l3d).and.calc_type.eq.2) then
        write (luout,*) 'l3d, calc_type =',l3d,calc_type
        call errquit ('fxc: calc_type-l3d logic wrong ',0, UNKNOWN_ERR)
      end if
c
c     Set up local copies of density matrix handles. The number of 
c     density matrices needed is ndens*nvec+2. The last two are scf
c     (unperturbed) density matrices.
c
c     g_dens_xc is the density matrix used in the XC routines
c
      if(l3d) then
c     -------------
c     calc_type = 2
c     -------------
c     the number of density matrices 
c     needed is nfock_xc+1 (typically 1 for perturbed density and 1
c     for unperturbed SCF density) for RDFT and 2*nfock_xc+2 (typically
c     2 for perturbed densities and 2 for unperturbed SCF densities) for
c     UDFT.
c
         if(ndens.eq.1) then
            dims(1) = nfock_xc+1
         else
            dims(1) = 2*nfock_xc+2
         endif
         if (.not.MA_alloc_Get(MT_int,dims(1),'XCDens',lgxcd,igxcd))
     &      call errquit('fock_xc: cannot allocate xcdens',0, MA_ERR)
c
c     alpha bit
c
         alo(1) = 1
         ahi(1) = nfock_xc
         blo(1) = 1
         bhi(1) = nfock_xc
         do i=1,nfock_xc
            int_mb(igxcd+i-1)=
     =           dft_npatch(g_dens,i,nbf_ao,geom,ao_bas_han)
         enddo
c
c     beta bit
c
         if(ndens.eq.2)  then
            alo(1) = nfock_xc*2+1
            ahi(1) = nfock_xc*3
            blo(1) = nfock_xc+1+1
            bhi(1) = nfock_xc*2+1
            do i=blo(1),bhi(1)
               int_mb(igxcd+i-1)=dft_npatch(g_dens,i+nfock_xc-1,
     N              nbf_ao,geom,ao_bas_han)
            enddo
         endif
c Daniel (1-10-13): We will need to do things differently for TDDFT
c gradients because we are (possibly) feeding in density matrices for 
c many excited states.
      elseif (calc_type.eq.5) then
c     -----------------------
c     calc_type = 5
c     -----------------------
c Daniel (1-10-13): For the moment, I want the calculation to work for
c a single root.  However, we at this point should have nfock+ndens
c density matrices (nfock = nroots*ndens, and ndens accounts for the
c spin-polarization, and hence the number of ground state density 
c matrices).  I build the MA array anticipating testing the multiple
c root situation in the future.
         if (.not.MA_alloc_Get(MT_int,nfock+ndens,'XCDens',lgxcd,igxcd))
     &      call errquit('fock_xc: cannot allocate xcdens',0, MA_ERR)
         do i = 1, nfock+ndens
            int_mb(igxcd+i-1) = g_dens(i)
         enddo
      else
c     -----------------------
c     calc_type = 1 or 3 or 4
c     -----------------------
         if (.not.MA_alloc_Get(MT_int,2,'XCDens',lgxcd,igxcd))
     &      call errquit('fock_xc: cannot allocate xcdens',0, MA_ERR)
         int_mb(igxcd)=g_dens(1)
         if(ndens.eq.2) int_mb(igxcd+1)=g_dens(3)
      endif
      if (calc_type .eq. 2) then
         if (.not. bgj_get_scf_dens(g_dens_scf))
     &        call errquit('fock_xc: cant get scf density handles',0,
     &       UNKNOWN_ERR)
         alo(1) = 1
         ahi(1) = nbf_ao
         blo(1) = nfock_xc+1
         bhi(1) = nfock_xc+1
         int_mb(igxcd+nfock_xc) =
     =        ga_create_atom_blocked(geom,ao_bas_han,'gdens+1')
         call ga_copy(g_dens_scf(1),int_mb(igxcd+nfock_xc))
         if(ndens.eq.2) then
            blo(1) = 2*nfock_xc+1+1
            bhi(1) = 2*nfock_xc+1+1
            int_mb(igxcd+2*nfock_xc+1) =
     =        ga_create_atom_blocked(geom,ao_bas_han,'gdens+2')
            call ga_copy(g_dens_scf(2),int_mb(igxcd+2*nfock_xc+1))
         endif
c Daniel (1-10-13): The SCF density matrix is stored differently than
c it is for calc_type = 1, 3, or 4
      elseif (calc_type.eq.5) then
c Use nroots to determine where the ground state density matrix is.
         nroots = nfock/ndens
         g_dens_scf(1) = int_mb(igxcd+nroots)
         if (ndens.eq.2) g_dens_scf(2) = int_mb(igxcd+2*nroots+1)
      else
         g_dens_scf(1) = int_mb(igxcd)
         if (ndens.eq.2) g_dens_scf(2) = int_mb(igxcd+1)
      endif
c
c     Prepare to call xc_quadv0
c
      g_wght_dum = -1
      g_xyz_dum = -1
      g_nq_dum = -1
      wght_GA_dum = .false.
      rho_n_dum = 0
      Exc(1) = 0.0d0
      Exc(2) = 0.0d0
      tdum = 0
      if(oprint) then
         if(ga_nodeid().eq.0) then
            write(luout,*) 'gdenscf1',g_dens_scf(1)
            call util_flush(luout)
         endif
         call ga_print(g_dens_scf(1))
         if(ndens.eq.2) then
            if(ga_nodeid().eq.0) then
               write(luout,*) 'gdenscf2',g_dens_scf(2)
               call util_flush(luout)
            endif
            call ga_print(g_dens_scf(2))
         endif
         if(ga_nodeid().eq.0) then
            write(luout,*) 'g_dens_xc1',int_mb(igxcd)
            call util_flush(luout)
         endif
         if(l3d) then
            do i=1,ndens*(nfock_xc+1)
               alo(1)=i
               ahi(1)=i
               call ga_print(int_mb(igxcd+i-1))
            enddo
         else
            call ga_print(int_mb(igxcd))
         endif
         if(ndens.eq.2.and.(.not.l3d)) then
            if(ga_nodeid().eq.0) then
               write(luout,*) 'g_dens_xc2',int_mb(igxcd+1)
               call util_flush(luout)
            endif
            call ga_print(int_mb(igxcd+1))
         endif
      endif  
      call xc_init_index(geom,ndens,nbf_ao,ao_bas_han,
     &   g_dens_scf, l_cntoce, k_cntoce,
     &   l_cntobfr, k_cntobfr, l_cetobfr, k_cetobfr,
     &   l_rdens_atom, k_rdens_atom)
c     
c     check if we can get the grid from a file
c         
c
c Daniel (2-18-13): For calc_type = 5, we don't use a 3D global array
c and therefore don't want l3d to be set equal to true. 
c      if (.not.lcgmin) then
      if (.not.lcgmin.and.calc_type.ne.5) then
       call grid_quadv0_gen(rtdb, 
     &     int_mb(igxcd),
     &     g_xc, nexc,
     &     rho_n_dum, Exc,
     &     nfock_xc, calc_type, tdum, .true., triplet)
      else
       call grid_quadv0_gen(rtdb, 
     &     int_mb(igxcd),
     &     g_xc, nexc,
     &     rho_n_dum, Exc,
     &     nfock_xc, calc_type, tdum, l3d, triplet)
      end if
c
      call dft_store_rhon(rho_n_dum)
c
c     Clean up
c
      if(l3d) then
         do i=1,dims(1)
         if (.not. ga_destroy(int_mb(igxcd+i-1)))
     &      call errquit('fock_xc: could not destroy g_dens_xc',0,
     &       GA_ERR)
         enddo
      endif
      if (.not.ma_free_heap(lgxcd))
     &   call errquit('fockxc: cannot pop stack',0, MA_ERR)
c
      ndim = ga_ndim(g_xc)
      if(ndens.eq.2.and.ndim.eq.3)  then
c
c     need to go from 2nfock ga to 4nfock (and sign flip)
c
c     move beta from block 2 to block 4
         alo(1) = nfock_xc+1
         ahi(1) = nfock_xc*2
         blo(1) = nfock_xc*3+1
         bhi(1) = nfock_xc*4
         call nga_copy_patch('N',g_xc,alo,ahi,
     &      g_xc,blo,bhi)
c     move alpha from block 1 to block 2
         alo(1) = 1
         ahi(1) = nfock_xc
         blo(1) = nfock_xc+1
         bhi(1) = nfock_xc*2
         call nga_copy_patch('N',g_xc,alo,ahi,
     &      g_xc,blo,bhi)
c     zero block 1 that was occupied by alpha
         call nga_zero_patch(g_xc, alo, ahi)
c     flip sign to make it consistent with HF K (Exchange)
         call ga_scale(g_xc,-1d0)
      endif
c
      call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr, l_rdens_atom)
c
      if (calc_type .eq. 2) then
         if (.not. ga_destroy(g_dens_scf(1)))
     &      call errquit('fock_xc: could not destroy DM',1, GA_ERR)
         if (ndens. eq. 2) then
            if (.not. ga_destroy(g_dens_scf(2)))
     &         call errquit('fock_xc: could not destroy DM',2, GA_ERR)
         endif
      endif
c
      if (oprint) then
         if(ga_nodeid().eq.0)write(luout,*)'XC matrix at end fock_xc:'
         tdum=ga_ddot(g_xc,g_xc)
         if(ga_nodeid().eq.0)write(luout,*)'g_xc dotproduct',tdum
         do i=1,nfock
            alo(1)=i
            ahi(1)=i
            tdum=nga_ddot_patch(g_xc, 'N', alo, ahi, 
     &         g_xc, 'N', alo, ahi) 
            if(ga_nodeid().eq.0) then
               write(luout,*) ' matrix no. ',i,tdum
               call util_flush(luout)
            endif
            if(abs(tdum).gt.1d-6) then
            if(ga_nodeid().eq.0) then
               call nga_inquire(g_xc, ityp, ndim, dims)
               write(luout,*) ' ndim ',ndim,' dims ',dims
               write(luout,*) ' alo ',alo
               write(luout,*) ' ahi ',ahi
               call util_flush(luout)
            endif
            call nga_print_patch(g_xc,alo,ahi,0)
            endif
         enddo
      endif
      return
      end
C> 
C> \brief Return a 2D slice of a 3D array of density matrices
C>
C> If the density matrices are stored in a 3D GA in the form
C> g_dens(1:nmat,1:nbf,1:nbf) then this routine returns a 2D GA in form
C> g_dens1(1:nbf,1:nbf) that equals g_dens(imat:imat,1:nbf,1:nbf).
C>
      integer function dft_npatch(g_dens,imat,nbf,geom,basis)
      implicit none
      integer g_dens ! N-D GA
      integer imat
      integer g_dens1 ! 2-D GA
      integer nbf
      integer geom,basis
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer alo(3), ahi(3), blo(2), bhi(2)
c
c     2-d GA dens1 
c
      g_dens1 = ga_create_atom_blocked(geom,basis,'gdens1')
      alo(1) = imat
      ahi(1) = imat
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens1,blo,bhi)
      dft_npatch=g_dens1
      return
      end
      subroutine dft_store_rhon(dft_rho_n)
      implicit none
#include "dft_conv.fh"
      double precision dft_rho_n
      rho_n = dft_rho_n
      end
