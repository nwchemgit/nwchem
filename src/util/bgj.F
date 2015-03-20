c
c     Miscellaneous utility routines used by BGJ
c     Some will evolve into generally useful routines
c
c     BGJ - 8/98
c
c     $Id$
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
c     Returns either 'hf', 'dft' or 'hyb' in name
c
      function bgj_get_scf_method(name)
c
      implicit none
#include "errquit.fh"
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
         call errquit('bgj_get_scf_method: error reading rtdb',bgj_rtdb,
     &       RTDB_ERR)
      endif
      if (itype.eq.1) then
         name = 'hf'
      else if (itype.eq.2) then
         name = 'dft'
      else if (itype.eq.3) then
         name = 'hyb'
      else
         call errquit('bgj_get_scf_method: illegal type in rtdb',1,
     &       RTDB_ERR)
      endif
c      write(*,*)'*** bgj_get_scf_method called: ',name
      bgj_get_scf_method = .true.
c
      return
      end

c
c     Returns whether there is an XC contribution in this calculation
c
      function bgj_have_xc()
c
      implicit none
#include "errquit.fh"
      logical bgj_have_xc
c
      external bgj_get_scf_method
      logical bgj_get_scf_method
c
      character*3 name
c
      if (.not.bgj_get_scf_method(name)) then
         call errquit('bgj_have_xc: error getting method name',1,
     &       UNKNOWN_ERR)
      endif
      bgj_have_xc = name .eq. 'dft' .or. name .eq. 'hyb'
c      write(*,*)'*** bgj_have_xc called: ',name,' ',bgj_have_xc
c
      return
      end

c
c     Returns coefficient of HF exchange for the current job
c
      function bgj_kfac()
c
      implicit none
      double precision bgj_kfac
c
#include "../nwdft/include/cdft.fh"
c
      external bgj_have_xc
      logical bgj_have_xc
c
      if (bgj_have_xc()) then
         bgj_kfac = xfac(1)
      else
         bgj_kfac = 1.d0
      endif
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
#include "errquit.fh"
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
     $        call errquit('bgj_restricted: problem getting scftype',0,
     &       RTDB_ERR)
         bgj_restricted = scftype .ne. 'UHF'
      else                           ! DFT case
c         write(*,*)'*** bgj_restricted: ipol ',ipol
         bgj_restricted = ipol .eq. 1
         if (ipol .ne. 1 .and. ipol .ne. 2)
     $        call errquit('bgj_restricted: problem with ipol',0,
     &       INPUT_ERR)
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
#include "errquit.fh"
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
      integer rtdb, ispin, nmo(2), g_vecs(2), nsets, nbf_vecs,
     &     l_eval, k_eval, l_occ, k_occ
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      double precision rhffact
      logical status
c
      rtdb = bgj_rtdb
c
c     get MO vectors from file
c
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
cDEBUG
     &then
c       if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
c    &    call errquit(
c    &         'bgj_get_scf_dens: SCF/DFT MO vectors not defined',0,
c    &         RTDB_ERR)
         call errquit('bgj_get_scf_dens: DFT MO vectors not defined',0,
     &                RTDB_ERR)
      endif
cDEBUG
      status = movecs_read_header(movecs_in, title_vecs, basis_vecs,
     &         scftype_vecs, nbf_vecs, nsets, nmo, 2)
c
c     ipol  1  closed-shell (RKS)
c     ipol  2  open-shell (both ROKS and UKS)
c     nsets 1  RKS and ROKS
c     nsets 2  UKS
c
      if (ipol .ne. nsets)then
        if (ipol.eq.1.and.nsets.eq.2) then
          write (6,*) 'bgj_get_scf_dens: ERROR ipol, nsets:',ipol,nsets
          call errquit('bgj_get_scf_dens: ERROR ipol, nsets disagree',2,
     &         INPUT_ERR)
        endif
      endif
c
c     Should check much more info than just nbf for consistency
c
c     allocate and initialize global and local memory
c
c     mo-vectors
c
      if (nsets .eq. 1)then
         g_vecs(1) = ga_create_atom_blocked (geom, ao_bas_han,
     &               'bgj_get_scf_dens:mo vectors')
         if (ipol .eq. 2) then
           g_vecs(2) = g_vecs(1)
         endif
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
c     get mo eigenvectors
c
      if (nbf_ao .ne. nbf_vecs)then
         write(6,*)'bgj_get_scf_dens movecs output = ',movecs_in
         call errquit('bgj_get_scf_dens: could not read mo vectors',911,
     &       DISK_ERR)
      else 
         status = .true.
c
c        occupation numbers (not used, but necessary for movecs_read)
c
c...     should do k_occ for both spins, in case used at some point...
c
         if(.not.MA_push_get(MT_DBL, nbf_ao, 'occ. numbers',
     &                       l_occ, k_occ))
     &      call errquit('bgj_get_scf_dens: could not allocate l_occ',1,
     &                   MA_ERR)
         if(.not.ma_push_get(mt_dbl,nbf_ao,'evals',
     $                       l_eval,k_eval))
     $      call errquit('bgj_get_scf_dens: problem with ma for evals',
     $                   1, MA_ERR)
         do ispin = 1, nsets
c
c           get orbitals (i.e. fock eigenvectors) and fock eigenvalues
c
            status = status .and.
     &               movecs_read(movecs_in, ispin, dbl_mb(k_occ),
     &               dbl_mb(k_eval), g_vecs(ispin))
         enddo
         if(.not.ma_pop_stack(l_eval))
     $        call errquit('bgj_get_scf_dens: problem popping evals',1,
     &                      MA_ERR)
         if(.not.ma_pop_stack(l_occ))
     $        call errquit('bgj_get_scf_dens: problem popping occ',1,
     &                      MA_ERR)
      endif
c
      if (.not.status)then
         write(6,*)'bgj_get_scf_dens movecs output = ',movecs_in
         call errquit('bgj_get_scf_dens: could not read mo vectors',917,
     &       DISK_ERR)
      endif
c
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
      if (nsets .eq. 1)then
         if (.not. ga_destroy(g_vecs(1)))
     &      call errquit('bgj_get_scf_dens: destroying mos',0, GA_ERR)
      else
         if (.not. ga_destroy(g_vecs(1)))
     &      call errquit('bgj_get_scf_dens: destroying alpha mos',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_vecs(2)))
     &      call errquit('bgj_get_scf_dens: destroying beta mos',0,
     &       GA_ERR)
      endif
c
      bgj_get_scf_dens = .true.
c
      return
      end
c
c     Checks for unrestricted 2nd derivative calculation on equal
c     number of alpha and beta electrons with functional having a
c     restricted/unrestricted discontinuity in that case
c
      subroutine bgj_check_xc_u2d()
c
      implicit none
#include "errquit.fh"
#include "bgj_common.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "../nwdft/include/cdft.fh"
      double precision c
      logical dontcare
c
      if (ipol .eq. 2) then
         if(.not.rtdb_get(bgj_rtdb,'dft:dontcare',mt_log, 1,dontcare))
     .        dontcare=.false.
         if(dontcare) return
         if (noc(1) .eq. noc(2)) then
c           Add up coefficents for problematic functionals
            c = abs(cfac(1)) + abs(cfac(7)) + abs(cfac(8))
     $        + abs(cfac(9)) + abs(cfac(10)) + abs(cfac(11))
     $        + abs(cfac(6)) + abs(cfac(3))
            if (c .gt. 1d-8) then
      write(*,*)'*** For a closed-shell system the unrestricted version'
      write(*,*)'*** of this XC functional yields different second'
      write(*,*)'*** derivative results than the restricted form due to'
      write(*,*)'*** a non-zero second derivative of the interpolation'
      write(*,*)'*** function at rhoa = rhob.  You may run this job as'
      write(*,*)'*** a restricted calculation if you wish.'
      call errquit('Problem with XC functional 2nd derivative',0,
     &       CALC_ERR)
            endif
         endif
      endif
c
      return
      end
