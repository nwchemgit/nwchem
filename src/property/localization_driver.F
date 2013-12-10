      subroutine localization_driver (rtdb, basis, geom, loctype)

c     $Id$

c     ==================================================================
c
c     purpose: perform orbital localization, HF and DFT


c     input:   rtdb    - nwchem run time database 
c              basis   - basis set handle
c              geom    - geometry handle
c              loctype - localization type
c                        1: Pipek-Mezey
c                        2: Boys (not yet interfaced with this driver)
c
c     driver witten by by J. Autschbach, jochena@buffalo.edu

c     ==================================================================

      implicit none

#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "util.fh"


c     subroutine arguments:

      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer loctype

c     local GA handles:
      integer g_dipel, g_smat, g_sc, g_tran
      integer g_movecs(2), g_cmo(2), g_temp

c     MA variables:
      integer l_c, k_c, l_sc, k_sc, l_eval, k_eval, l_occ, k_occ

c     other local variables:

      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo, nloc, i

      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)

      double precision atn, tol2e, val
      character*3 scftyp

      logical     oskel, status
      data tol2e   /1.0d-10/

      integer maxnloc
      parameter (maxnloc = 10000)
      integer iloc(maxnloc)

      double precision dummy(3)
      double precision origin(3)
      data origin/0d0,0d0,0d0/

      character*(256) lmotrans
      logical file_write_ga
      external file_write_ga

      logical debug
      logical oprint

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

c     ==================================================================

      debug = .true. .and. ga_nodeid().eq.0 ! .true. during development

      if (debug) write(luout,*) 'hello from localization_driver'

      if (loctype.lt.1 .or. loctype.gt.2)
     &   call errquit('loc_driver: loctype out of range',0, RTDB_ERR)

      oprint = util_print('information', print_low)
      oprint = oprint .and. (ga_nodeid().eq.0)
      
c     ----------------------
c     start-up: print header
c     ----------------------
      
      if (oprint) then
        write (LuOut,*) 'Entering MO Localization driver routine'
      endif
      
      if (oprint) then
        write (luout,9000)
c       write (luout,9010) 
        if (loctype.eq.1) then
          call util_print_centered(luout,
     &       'Pipek-Mezey orbital localization', 32, .true.)
        elseif (loctype.eq.2) then
          call util_print_centered(luout,
     &       'Boys orbital localization', 25, .true.)
        end if
      end if
      

c     -------------------------------------------------------
c     assemble some data:
c     MO coeffs, orbital energies, density matrix, occupation
c     numbers, basis set info, matrix elements, etc.
c     -------------------------------------------------------

      if (.not. bas_numbf(basis,nbf)) call
     &   errquit('loc_driver: could not get nbf',0, BASIS_ERR)
c     allocate dbl_mb(k_occ) = fractional occupation numbers:
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('loc_driver: ma_push_get failed k_occ',0,MA_ERR)
c     allocate dbl_mb(leval) = orbital energies:
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('loc_driver: ma_push_get failed k_eval',0,MA_ERR)

      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,g_movecs,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)

      if (oprint .and. debug)
     &   write (luout,*) 'loc driver: nclosed,nopen,nvirt',nclosed(1),
     &   nopen(1), nvirt(1)
      
c     Exit gracefully if we have an open-shell system:
      if (scftyp.eq.'UHF') then
        if (oprint) write(luout,9020)
        goto 7000
c       Note: ndens = 1 means closed shell
c             ndens = 3 open shell, then g_dens has 3 components
c                       with 1=alpha, 2=beta, 3=total
      endif
      
c     perform some sanity checks for the orbital occupations:
      
c     for RHF system there should be no open shell orbitals
      if (nopen(1).ne.0) call
     &   errquit('loc_driver: nopen.ne.0',0, CALC_ERR)
c     in this case, nclosed and nvirt should add up to nmo
      if ((nclosed(1)+nvirt(1)).ne.nmo) call
     &   errquit('loc_driver: wrong no. of orbitals',0, CALC_ERR)


c     maxnloc is hard coded. We should allocate array iloc
c     dynamically but for now let's make sure we don't get out of
c     bounds:

      if (nmo.gt.maxnloc) call
     &   errquit('loc_driver: maxnloc too small',0, BASIS_ERR)

      
c     --------------------------------------------
c     setup complete. Now proceed to the localization
c     --------------------------------------------

c
c     Integral initialization
c     
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)

      if (loctype.eq.1) then

c       ========================
c       Pipek-Mezey localization
c       ========================
        

c       AO Overlap Matrix S:

        g_smat  = ga_create_atom_blocked(geom, basis, 'loc:smat')
        call ga_zero(g_smat)
        call int_1e_ga(basis, basis, g_smat, 'overlap', .false.)

c       PM localization needs S*C: 
 
        if (.not. ga_create(MT_DBL, nbf, nmo, 'loc:sc',
     $        nbf, 0, g_sc)) call errquit('loc_driver: sc',0, GA_ERR)
         call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $        1.0d0, g_smat, g_movecs(1), 0.0d0, g_sc)

c        allocate some memory used in the localization routine:

         if (.not. ma_push_get(mt_dbl, 2*nbf, 'sc', l_sc, k_sc))
     $        call errquit('ma for sc', 0, MA_ERR)
         if (.not. ma_push_get(mt_dbl, 2*nbf, 'c', l_c, k_c))
     $        call errquit('ma for c', 0, MA_ERR)

c        localize core and occupied orbitals

         do i = 1, nclosed(1)
            iloc(i) = i
         end do
 0       nloc = nclosed(1)
c        jochen: comment:
c        the localization routines were already available
c        in nwchem
         call localizePM(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $        nloc, iloc, nbf, nmo, g_movecs, g_sc)


c        clean up some temp arrays
         if (.not. ga_destroy(g_sc)) call errquit(
     &      'loc_driver: error destroying g_sc',0, GA_ERR)
         if (.not. ma_pop_stack(l_c)) call errquit(
     &      'loc_driver: error MA pop c',0, MA_ERR)
         if (.not. ma_pop_stack(l_sc)) call errquit(
     &      'loc_driver: error MA pop sc',0, MA_ERR)

c        calculate the LMO transformation matrix 
c        from C(dagger,canonical) S C(locorb)

c        first, read the CMOs again -> g_cmo

         call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &      nvirt,scftyp,g_cmo,dbl_mb(k_occ),
     &      dbl_mb(k_eval),nmo)

         if (.not. ga_create(MT_DBL, nclosed(1), nbf, 'loc:temp',
     &      nbf, 0, g_temp))
     &      call errquit('loc_driver: temp',0, GA_ERR)
         call ga_zero(g_temp)
         call ga_dgemm('t', 'n', nclosed(1), nbf, nbf, 
     $      1.0d0, g_cmo(1), g_smat, 0.0d0, g_temp)

         if (.not. ga_destroy(g_smat)) call errquit(
     &      'loc_driver: error destroying g_smat',0, GA_ERR)
         if (.not. ga_create(MT_DBL, nclosed(1), nclosed(1),
     &      'loc:smat',nbf, 0, g_smat))
     &      call errquit('loc_driver: smat',0, GA_ERR)
         call ga_zero(g_smat)
         call ga_dgemm('n', 'n', nclosed(1), nclosed(1), nbf, 
     $      1.0d0, g_temp, g_movecs(1), 0.0d0, g_smat)

         if (.not. ga_destroy(g_temp)) call errquit(
     &      'loc_driver: error destroying g_temp',0, GA_ERR)
         if (.not. ga_destroy(g_cmo(1))) call errquit(
     &      'loc_driver: error destroying g_cmo(1)',0, GA_ERR)

c        loc transform is in g_smat. Write to file

         call util_file_name('lmotrans',.true.,.true.,lmotrans)
         if(.not.file_write_ga(lmotrans,g_smat)) call errquit
     $      ('loc_driver: could not write lmotrans',0, DISK_ERR)

c        smat, temp not needed anymore
         if (.not. ga_destroy(g_smat)) call errquit(
     &      'loc_driver: error destroying g_smat',0, GA_ERR)

      else if (loctype.eq.2) then

c       =================
c       Boys localization
c       =================

        write (luout,*) 'Boys not yet implemented. Will do nothing :-)'
        continue ! nothing to do right now

      end if                    ! loctype

c     -----------------------------------
c     localization done. Store MO vectors
c     -----------------------------------

      if (debug) write (luout,*) 'calling hnd_vec_write'
      call hnd_vec_write(rtdb,geom,basis,nbf,nclosed,nopen,
     &   nvirt,scftyp,g_movecs,dbl_mb(k_occ),
     &   dbl_mb(k_eval),nmo, 'locorb.movecs')


c     --------
c     clean up
c     --------


      if (.not.ga_destroy(g_movecs(1))) call 
     &    errquit('loc_driver: ga_destroy failed g_movecs',0,GA_ERR)

      if (.not.ma_chop_stack(l_occ)) call
     &   errquit('loc_driver: ma_chop_stack failed k_occ',l_occ,MA_ERR)

      call schwarz_tidy()
      call int_terminate()

 7000 continue

c     ---------------------------------------
c     Localization done. return
c     ---------------------------------------
      
      if (oprint) then
        write (LuOut,*) 'Exiting Localization driver routine'
      endif

c     ==================================================================

      return



 9000 format(//80(1h-)//,
     &   10x,20(1h*),/,10x,'*** LOCALIZATION ***',/,10x,20(1h*)/)
 9010 format(' Response module for NWChem and dynamic CPKS solver'/
     &       ' developed by J. Autschbach and coworkers, SUNY Buffalo'/
     &       ' The methodology used in this program is described in '/
     &       ' J. Chem. Phys. 123 (2005), 114103'/
     &       ' J. Chem. Phys. 122 (2005), 224115'/
     &       ' J. Chem. Phys. 122 (2005), 074105'/
     &       ' Comp. Lett. 3 (2007), 131-150 (contact JA for a copy)'/
     &       ' Please cite this work in publications based on results'/
     &       ' obtained with this code. Thank you!'/)
 9020 format(/,10x,'Open shell calculations not yet implemented',
     &      ' for the Localization module. Exiting')
      end
