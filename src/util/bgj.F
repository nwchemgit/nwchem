c
c     Miscellaneous utility routines used by BGJ
c     Some will evolve into generally useful routines
c
c     BGJ - 8/98
c
c     $Id: bgj.F,v 1.3 1999-05-08 13:00:40 bjohnson Exp $
c

c
c     Returns print level for BGJ debug writes
c
      function bgj_print()
c
      implicit none
      integer bgj_print
#include "bgj_common.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer print_level
      save print_level
      data print_level / -1 /
c
      if (print_level .lt. 0) then
         if (.not. rtdb_get(bgj_rtdb, 'bgj:print', mt_int,
     &        1, print_level)) then
            print_level = 0
         endif
      endif
      bgj_print = print_level
c
      return
      end

c
c     Returns current rtdb handle
c
      function bgj_get_rtdb_handle()
c
      implicit none
      integer bgj_get_rtdb_handle
#include "bgj_common.fh"
c
c      write(*,*)'*** bgj_get_rtdb_handle called',bgj_rtdb
      bgj_get_rtdb_handle = bgj_rtdb
c
      return
      end

c
c     Returns either 'HF', 'DFT' or 'HYB' in name
c
      function bgj_get_scf_method(name)
c
      implicit none
      logical bgj_get_scf_method
      character*3 name
c
#include "rtdb.fh"
#include "mafdecls.fh"
c!!!
#include "bgj_common.fh"
c!!!
c
      integer itype
c
      if (.not. rtdb_get(bgj_rtdb, 'bgj:scf_type', mt_int,
     &     1, itype)) then
         call errquit('bgj_get_scf_method: error reading rtdb',1)
      endif
      if (itype.eq.1) then
         name = 'HF'
      else if (itype.eq.2) then
         name = 'DFT'
      else if (itype.eq.3) then
         name = 'HYB'
      else
         call errquit('bgj_get_scf_method: illegal type in rtdb',1)
      endif
c      write(*,*)'*** bgj_get_scf_method called: ',name
      bgj_get_scf_method = .true.
c
      return
      end

c
c     Returns whether there is an XC contribution is this calculation
c
      function bgj_have_xc()
c
      implicit none
      logical bgj_have_xc
c
      external bgj_get_scf_method
      logical bgj_get_scf_method
c
      character*3 name
c
      if (.not.bgj_get_scf_method(name)) then
         call errquit('bgj_have_xc: error getting method name',1)
      endif
      bgj_have_xc = name .eq. 'DFT' .or. name .eq. 'HYB'
c      write(*,*)'*** bgj_have_xc called: ',name,' ',bgj_have_xc
c
      return
      end

c
c     Returns whether J fitting is involved is this calculation
c
      function bgj_have_j_fit()
c
      implicit none
      logical bgj_have_j_fit
c
#include "../nwdft/include/cdft.fh"
c
      external bgj_have_xc
      logical bgj_have_xc
c
      bgj_have_j_fit = bgj_have_xc()
      if (bgj_have_j_fit) then
         bgj_have_j_fit = CDFIT
      endif
c      write(*,*)'*** bgj_have_j_fit called: ',name,' ',bgj_have_j_fit
c
      return
      end

c
c     Returns CD basis handle, or -1 if there's not one
c
      function bgj_cd_bas_han()
c
      implicit none
      integer bgj_cd_bas_han
c
#include "../nwdft/include/cdft.fh"
c
      external bgj_have_xc
      logical bgj_have_xc
c
      bgj_cd_bas_han = -1
      if (bgj_have_xc()) then
         if (CDFIT) then
            bgj_cd_bas_han = CD_bas_han
         endif
      endif
c
      return
      end

c
c     Returns whether the calculation is restricted
c
      function bgj_restricted()
c
      implicit none
      logical bgj_restricted
c
c !!! BGJ
#include "bgj_common.fh"
#include "rtdb.fh"
#include "../nwdft/include/cdft.fh"
c
      external bgj_have_xc
      logical bgj_have_xc
c
      character*10 scftype
c
      if (.not.bgj_have_xc()) then   ! HF case
         if (.not. rtdb_cget(bgj_rtdb, 'scf:scftype', 1, scftype))
     $        call errquit('bgj_restricted: problem getting scftype',0)
         bgj_restricted = scftype .ne. 'UHF'
      else                           ! DFT case
c         write(*,*)'*** bgj_restricted: ipol ',ipol
         bgj_restricted = ipol .eq. 1
         if (ipol .ne. 1 .and. ipol .ne. 2)
     $        call errquit('bgj_restricted: problem with ipol',0)
      endif
c      write(*,*)'*** bgj_restricted called: ',bgj_restricted
c
      return
      end

c
c     Returns handles to AO-basis SCF density matrices
c
      function bgj_get_scf_dens(g_dens)
c
      implicit none
      logical bgj_get_scf_dens
      integer g_dens(2)
c
c !!! BGJ
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "../nwdft/include/cdft.fh"
#include "bgj_common.fh"
      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c !!! BGJ
      integer rtdb, ispin, noc(2), nmo(2), g_vecs(2), nsets, nbf_vecs,
     &     l_eval, k_eval, l_occ, k_occ
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      double precision rhffact
      logical status
c
      rtdb = bgj_rtdb
c
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('bgj_get_scf_dens: rtdb_get of noc failed', 0) 
c
c     allocate and initialize global and local memory
c
c     mo-vectors
c
      if (ipol .eq. 1)then
         g_vecs(1) = ga_create_atom_blocked (geom, ao_bas_han,
     &               'bgj_get_scf_dens:mo vectors')
      else 
         g_vecs(1) = ga_create_atom_blocked (geom, ao_bas_han,
     &               'bgj_get_scf_dens:alpha mo vectors')
         g_vecs(2) = ga_create_atom_blocked (geom, ao_bas_han,
     &               'bgj_get_scf_dens:beta mo vectors')
      endif 
c
c     global density
c
      if (ipol .eq. 1)then
         g_dens(1) = ga_create_atom_blocked (geom, ao_bas_han,
     &                 'bgj_get_scf_dens:density')
      else 
         g_dens(1) = ga_create_atom_blocked (geom, ao_bas_han,
     &                 'bgj_get_scf_dens:alpha density')
         g_dens(2) = ga_create_atom_blocked (geom, ao_bas_han,
     &                 'bgj_get_scf_dens:beta density')
      endif 
c
c     occupation numbers (not used, but necessary for movecs_read)
c****
c**** should do k_occ for both spins, in case used at some point...
c****
      status = MA_push_get(MT_DBL, nbf_ao, 'occ. numbers',
     &         l_occ, k_occ)
      if (.not.status) 
     &   call errquit('bgj_get_scf_dens: could not allocate l_occ',1)
      if(.not.ma_push_get(mt_dbl,nbf_ao,'evals',
     $     l_eval,k_eval))
     $     call errquit('bgj_get_scf_dens: problem with ma for evals',1)
c
c     get MO vectors from file
c
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     &   call errquit('bgj_get_scf_dens: DFT MO vectors not defined',0)
      status = movecs_read_header(movecs_in, title_vecs, basis_vecs,
     &         scftype_vecs, nbf_vecs, nsets, nmo, 2)
c
      if (ipol .ne. nsets)then ! 1 for RHF, 2 for UHF
         write (6,*) 'bgj_get_scf_dens:  ERROR ipol, nsets:',ipol,nsets
         call errquit('bgj_get_scf_dens:  ERROR ipol, nsets disagree',2)
      endif 
c
c     Should check much more info than just nbf for consistency
c
c
c     get mo eigevectors
c
      if (nbf_ao .ne. nbf_vecs)then
         write(6,*)'bgj_get_scf_dens movecs output = ',movecs_in
         call errquit('bgj_get_scf_dens: could not read mo vectors',911)
      else 
         status = .true.
         do ispin = 1, ipol
c****
c****       get orbitals (i.e. fock eigenvectors) and fock eigenvalues
c****
            status = status .and.
     &               movecs_read(movecs_in, ispin, dbl_mb(k_occ),
     &               dbl_mb(k_eval), g_vecs(ispin))
         enddo 
      endif
c
      if (.not.status)then
         write(6,*)'bgj_get_scf_dens movecs output = ',movecs_in
         call errquit('bgj_get_scf_dens: could not read mo vectors',917)
      endif
c
      if(.not.ma_pop_stack(l_eval))
     $     call errquit('bgj_get_scf_dens: problem popping evals',1)
      if(.not.ma_pop_stack(l_occ))
     $     call errquit('bgj_get_scf_dens: problem popping occ',1)
c      
      if (ipol .eq. 1)then
         rhffact = 2.0d0
      else 
         rhffact = 1.0d0
      endif 
c
      do ispin = 1, ipol
c
c        dens = vecs*vecs
c
         call ga_dgemm('n', 't', nbf_ao, nbf_ao, noc(ispin), rhffact, 
     &        g_vecs(ispin), g_vecs(ispin), 0.0d0, g_dens(ispin))
         call ga_symmetrize(g_dens(ispin))
      enddo 

c
c     Free temporary storage
c
      if (ipol .eq. 1)then
         if (.not. ga_destroy(g_vecs(1)))
     &      call errquit('bgj_get_scf_dens: destroying mos',0)
      else 
         if (.not. ga_destroy(g_vecs(1)))
     &      call errquit('bgj_get_scf_dens: destroying alpha mos',0)
         if (.not. ga_destroy(g_vecs(2)))
     &      call errquit('bgj_get_scf_dens: destroying beta mos',0)
      endif 
c      write(*,*)'*** bgj_get_scf_dens: g_dens',g_dens
c      call ga_print(g_dens(1))
c      if (ipol.eq.2) call ga_print(g_dens(2))
c
      bgj_get_scf_dens = .true.
c
      return
      end
