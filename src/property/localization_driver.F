      subroutine localization_driver (rtdb, basis, geom, loctype)

c     $Id$

c     ==================================================================
c
c     purpose: perform orbital localization, HF and DFT


c     input:   rtdb    - nwchem run time database 
c              basis   - basis set handle
c              geom    - geometry handle
c              loctype - localization type
c                        1: Pipek-Mezey (PM)
c                        2: Boys
c                        3: IAO-IBO
c     
c
c     driver witten by by J. Autschbach, jochena@buffalo.edu
c     localization routines are in src/ddscf/localize.F

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
      integer g_uc(4), g_smat, g_sc, g_tran
      integer g_movecs(2), g_cmo(2), g_temp, g_tmp1, g_tmp2
      integer g_s2, g_s12, g_p12, g_p21
      integer g_ctilde, g_iao, g_mo

c     MA variables:
      integer l_c, k_c, l_sc, k_sc, l_eval, k_eval, l_occ, k_occ
      integer l_dip(3), k_dip(3)

c     other local variables:

      integer nclosed(2), nopen(2), nvirt(2), ncore,
     &  ndens, nbf, nmo, nloc, i

      integer info, dbg, ispin
      
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)

      double precision atn, tol2e, val, rtmp, sumd2, dx, dy, dz, rtemp
      character*3 scftyp

      logical     oskel, status
      data tol2e   /1.0d-10/

      integer maxnloc
      parameter (maxnloc = 10000)
      integer iloc(maxnloc)

      double precision dummy(3)
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      double precision angstrom
      parameter (angstrom = 0.529177d0)

      character*(256) lmotrans
      logical file_write_ga
      external file_write_ga

      logical debug, master
      logical oprint

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      logical int_normalize
      external int_normalize

c     variables for IAO basis
      integer
     &  nbf_iao, nbf_iao_mxprim, nbf_iao_mxang, nbf_iao_mxcont,
     &  nbf_iao_mxnbf_cn, nbf_iao_mxnbf_ce, nshells_iao,
     &  minbas      

      character*(19) pname

c     ==================================================================

      pname = 'localization_driver'

      dbg = 0
      master =  ga_nodeid().eq.0
      debug = (dbg>0) .and. master ! .true. during development

      if (debug) write(luout,*) 'entering '//pname

      if (loctype.lt.1 .or. loctype.gt.3)
     &   call errquit('loc_driver: loctype out of range',0, RTDB_ERR)

      oprint = util_print('information', print_low)
      oprint = oprint .and. ga_nodeid().eq.0
      
c     ----------------------
c     start-up: print header
c     ----------------------
      
      if (oprint) then
        write (luout,9000) ! general header
        if (loctype.eq.1) then
          call util_print_centered(luout,
     &       'Pipek-Mezey localization', 32, .true.)
        elseif (loctype.eq.2) then
          call util_print_centered(luout,
     &      'Boys localization', 25, .true.)
        elseif (loctype.eq.3) then
          call util_print_centered(luout,
     &       'IAO-IBO localization', 20, .true.)
        else
          write (luout,*) 'loctype =',loctype
          write (luout,*) 'LocType value not supported. Exiting.'
          return
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

      if (debug)
     &   write (luout,*) 'loc driver: nclosed,nopen,nvirt',nclosed(1),
     &   nopen(1), nvirt(1)
      
c     Skip localization if we have an open-shell system
c     (to be extended eventually, localizing alpha and
c     beta spin orbitals separately):
      
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
     &  errquit('loc_driver: nopen.ne.0',0, CALC_ERR)
      
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

c     =================================================================

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
     $        call errquit('loc_driver: ma for sc', 0, MA_ERR)
         if (.not. ma_push_get(mt_dbl, 2*nbf, 'c', l_c, k_c))
     $        call errquit('loc_driver: ma for c', 0, MA_ERR)

c        localize core and occupied orbitals

         ispin = 1

         do i = 1, nclosed(ispin)
            iloc(i) = i
         end do
         nloc = nclosed(ispin)
         
c        jochen: comment:
c        the PM localization routine was already available
c        in nwchem
         call localizePM(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     &        nloc, iloc, nbf, nmo, g_movecs, g_sc)


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

c        smat not needed anymore
         if (.not. ga_destroy(g_smat)) call errquit(
     &     'loc_driver: error destroying g_smat',0, GA_ERR)

c     =================================================================

      else if (loctype.eq.2) then

c       =================
c       Boys localization
c       =================

        do i = 1,4
          if (.not. ga_create(MT_DBL, nbf, nbf, 'uc',
     $       nbf, 0, g_uc(i))) call errquit('loc_driver:uc'
     &       ,i,GA_ERR)
          call ga_zero(g_uc(i))
         end do

         if (debug) write(luout,*) 'g_uc created'

c        dipole moment AO matrices -> uc(1) - uc(3)
         call int_dip_ga(basis, basis, g_uc(1), g_uc(2), g_uc(3))

         if (debug) write(luout,*) 'dipole done'

c        calculate dipole-AO times C (mo-coeffs), store in uc(i)
         if (.not. ga_create(MT_DBL, nbf, nmo, 'sc',
     $        nbf, 0, g_sc)) call errquit('loc_driver: sc',0, GA_ERR)
         do i = 1, 3
            call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $           1.0d0, g_uc(i), g_movecs, 0.0d0, g_sc)
            call ga_copy_patch('n',g_sc,1,nbf,1,nmo,g_uc(i),1,nbf,1,nmo)
         end do

         if (debug) write(luout,*) 'g_uc complete'

c       AO Overlap Matrix S:

        g_smat  = ga_create_atom_blocked(geom, basis, 'loc:smat')
        call ga_zero(g_smat)
        call int_1e_ga(basis, basis, g_smat, 'overlap', .false.)

        if (debug) write(luout,*) 'smat done'

c       store S C in array uc(4)
        call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $     1.0d0, g_smat, g_movecs(1), 0.0d0, g_uc(4))

        if (debug) write(luout,*) 'g_uc(4) done'
c     
         if (.not. ma_push_get(mt_dbl, 8*nbf, 'sc', l_sc, k_sc))
     $        call errquit('loc_driver: ma for sc', 0, MA_ERR)
         if (.not. ma_push_get(mt_dbl, 8*nbf, 'c', l_c, k_c))
     $        call errquit('ma for c', 0, MA_ERR)

         if (debug) write(luout,*) 'MA c, sc complete'

c        localize core and occupied orbitals:

         do i = 1, nclosed(1)
            iloc(i) = i
         end do
         nloc = nclosed(1)

c        calculate orbital centroids and determine B2 = the sum of 
c        square distances of the centroids from the coordinate origin.
c        The Boys localization maximizes B2. 
c        we can use g_sc as a temp array as it is not further used

         if (.not. ga_destroy(g_sc)) call errquit(
     &      'loc_driver: error destroying g_sc',0, GA_ERR)
         if (.not. ga_create(MT_DBL, nmo, nmo, 'sc',
     $      nbf, 0, g_sc)) call errquit('loc_driver: sc',0, GA_ERR)

         do i = 1,3
           if (.not. ma_push_get(mt_dbl, nmo, 'sc', l_dip(i), k_dip(i)))
     $        call errquit('loc_driver: ma for dip', i, MA_ERR)
         end do

         do i = 1,3
           call ga_zero(g_sc)
           call ga_dgemm('t','n',nmo,nmo,nbf, 1.0d0, g_movecs, g_uc(i),
     &        0.0d0, g_sc)
           call ga_get_diagonal(g_sc, dbl_mb(k_dip(i)))
         end do

         if (ga_nodeid().eq.0) then
           write (luout,'(/t6,a,t15,a,t27,a,t39,a,t49,a)')
     &        '#','<x>','<y>','<z>','<r>**2'
           write (luout,'(1x,53(''-''))')
         end if
         sumd2 = 0.0d0
         do i = 1,nmo
           dx = dbl_mb(k_dip(1) + i-1) 
           dy = dbl_mb(k_dip(2) + i-1) 
           dz = dbl_mb(k_dip(3) + i-1)
           rtmp = dx**2 + dy**2 + dz**2
           sumd2 = sumd2 + rtmp
           if (ga_nodeid().eq.0) then
             write (luout,'(1x,i5,4f12.4)') i,dx,dy,dz,rtmp
           end if
         end do
c         if (ga_nodeid().eq.0) then
c           write (luout,'(1x,a,2f12.4)') 'B2 before localization, av:',
c     &        sumd2, sumd2/float(nloc)
c         end if

         do i = 3,1,-1
           if (.not. ma_pop_stack(l_dip(i))) call errquit(
     &        'loc_driver: error MA pop dip',i, MA_ERR)
         end do

c        jochen: comment:
c        the Boys localization routine was already available
c        in nwchem
         call localizeFB(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $        nloc, iloc, nbf, nmo, g_movecs, g_uc)

c        calculate orbital centroids again and print information

         do i = 1,3
           if (.not. ma_push_get(mt_dbl, nmo, 'sc', l_dip(i), k_dip(i)))
     $        call errquit('loc_driver: ma for dip', i, MA_ERR)
         end do

         do i = 1,3
           call ga_zero(g_sc)
           call ga_dgemm('t','n',nmo,nmo,nbf, 1.0d0, g_movecs,g_uc(i),
     &        0.0d0, g_sc)
           call ga_get_diagonal(g_sc, dbl_mb(k_dip(i)))
         end do

         if (ga_nodeid().eq.0) then
           write (luout,'(/t6,a,t15,a,t27,a,t39,a,t49,a)')
     &        '#','<x>','<y>','<z>','<r>**2'
           write (luout,'(1x,53(''-''))')
         end if
         sumd2 = 0.0d0
         do i = 1,nmo
           dx = dbl_mb(k_dip(1) + i-1)
           dy = dbl_mb(k_dip(2) + i-1)
           dz = dbl_mb(k_dip(3) + i-1)
           rtmp = dx**2 + dy**2 + dz**2
           sumd2 = sumd2 + rtmp
           if (ga_nodeid().eq.0) then
             write (luout,'(1x,i5,4f12.4)') i,dx,dy,dz,rtmp
           end if
         end do
c         if (ga_nodeid().eq.0) then
c           write (luout,'(1x,a,2f12.4)') 'B2 after localization, av:',
c     &        sumd2, sumd2/float(nloc)
c         end if

         do i = 3,1,-1
           if (.not. ma_pop_stack(l_dip(i))) call errquit(
     &        'loc_driver: error MA pop dip',i, MA_ERR)
         end do

c        clean up  temp arrays:

         if (.not. ga_destroy(g_sc)) call errquit(
     &      'loc_driver: error destroying g_sc',0, GA_ERR)

         if (.not. ma_pop_stack(l_c)) call errquit(
     &      'loc_driver: error MA pop c',0, MA_ERR)

         if (.not. ma_pop_stack(l_sc)) call errquit(
     &      'loc_driver: error MA pop sc',0, MA_ERR)

         do i = 1, 4
            if (.not. ga_destroy(g_uc(i)))
     &        call errquit('loc_driver: error destroying uc',i, GA_ERR)
         end do

c        smat not needed anymore
         if (.not. ga_destroy(g_smat)) call errquit(
     &     'loc_driver: error destroying g_smat',0, GA_ERR)

c      ================================================================

       else if (loctype.eq.3) then

c       ====================
c       IAO-IBO localization
c       ====================
        
c       AO Overlap Matrix S1 -> g_smat:

        g_smat  = ga_create_atom_blocked(geom, basis, 'loc:smat')
        call ga_zero(g_smat)
        call int_1e_ga(basis, basis, g_smat, 'overlap', .false.)

c       auxiliary basis overlap S2 -> g_s2 basis needs to defined in the
c       input as "iao basis". we will now create the basis here and then
c       calculate the overlap and the mixed ao-iao basis overlap S12
        
        if (.not. bas_create(minbas, 'iao basis'))
     &    call errquit(pname//': cannot create iao bas', 86, BASIS_ERR)
        
        if (.not. bas_rtdb_load(rtdb, geom, minbas, 'iao basis'))
     &    call errquit(pname//': iao basis not present', 86, BASIS_ERR)
        
        if (.not. int_normalize(rtdb,minbas))
     &    call errquit(pname//': int_normalize failed iao', 86, INT_ERR)
        
        if (.not. bas_numbf(minbas, nbf_iao))
     &    call errquit(pname//': bas error nbf_iao', 86, BASIS_ERR)
        
        if (.not. bas_nprim_cn_max(minbas,nbf_iao_mxprim))
     &    call errquit(pname//': bas error iao mxprim', 86, BASIS_ERR)
        
        if (.not. bas_high_angular(minbas,nbf_iao_mxang))
     &    call errquit(pname//': bas error iao mxang', 86, BASIS_ERR)
        
        if (.not. bas_ncontr_cn_max(minbas,nbf_iao_mxcont))
     &    call errquit(pname//': bas error iao mxcont', 86, BASIS_ERR)
        
        if (.not. bas_nbf_cn_max(minbas,nbf_iao_mxnbf_cn))
     &    call errquit(pname//': bas error iao mcnbf_cn', 86, BASIS_ERR)
        
        if (.not. bas_nbf_ce_max(minbas,nbf_iao_mxnbf_ce))
     &    call errquit(pname//': bas error iao mxnbf_ce', 86, BASIS_ERR)
        
        if (.not. bas_numcont(minbas,nshells_iao))
     &    call errquit(pname//': bas error iao nshells', 86, BASIS_ERR)

c       create overlap for minbas, and the mixed basis-minbas
c       overlap S12.
c       Then we calculate the projectors P12 and P21
c       via solving S1 P12 = S12 for P12 etc.
        
        g_s2  = ga_create_atom_blocked(geom, minbas, 'loc:s2')
        call ga_zero(g_s2)
        call int_1e_ga(minbas, minbas, g_s2, 'overlap', .false.)

        if (.not. ga_create(MT_DBL, nbf, nbf_iao, 'loc:s12',
     &    nbf, 0, g_s12)) call errquit(pname//': S12',61, GA_ERR)

        call ga_zero(g_s12)
        call int_1e_ga(basis, minbas, g_s12, 'overlap', .false.)

        if (.not.ga_duplicate(g_s12, g_p12, 'g_p12')) call
     &    errquit(pname//': ga_dup P12',62,GA_ERR)

        call ga_copy(g_s12, g_p12)
        
#ifdef SCALAPACK
        info = ga_solve(g_smat, g_p12) 
        if (debug) write(luout,*)
     &    'info from ga_solve P12:', info
#else
        if (debug) write(luout,*) 'call 1 ma_solve' 
        call ma_solve(g_smat, g_p12)
#endif

        if (.not. ga_create(MT_DBL, nbf_iao, nbf, 'loc:p21',
     &    nbf_iao, 0, g_p21)) call errquit(pname//': P21',63, GA_ERR)

        call ga_transpose(g_s12, g_p21)
        if (debug) write(luout,*) 'P12 assembled' 

c       S12 no longer needed
        if (.not. ga_destroy(g_s12))
     &    call errquit(pname//': ga_destroy failed g_s12',61, GA_ERR)  

c : scalapack throws a 'ga_to_SL: rows/cols error' so we're using
c       ma_solve for the time being
c#ifdef SCALAPACK
c        info = ga_solve(g_s2, g_p21) 
c        if (debug) write(luout,*)
c     &    'info from ga_solve P21:', info
c#else
        if (debug) write(luout,*) 'call 2 ma_solve' 
        call ma_solve(g_s2, g_p21)
c#endif

c       g_s2 no longer needed
        if (.not. ga_destroy(g_s2)) call errquit(
     &    pname//': error destroying g_s2',0, GA_ERR)

        if (dbg>1) then
          call ga_print(g_p12)
          call ga_print(g_p21)
        end if

c       store P12 * P21 in array g_temp; we'll need it soon

        if (.not.ga_duplicate(g_smat, g_temp, 'g_temp')) call
     &    errquit(pname//': ga_dup temp', 71, GA_ERR)
        
        call ga_dgemm('n', 'n', nbf, nbf, nbf_iao,
     &    1.0d0, g_p12, g_p21, 0.0d0, g_temp)

         if (debug) write (luout,*) 'p12 p21 -> temp'        

c        P21 is no longer needed
         
         if (.not. ga_destroy(g_p21))
     &     call errquit(pname//': ga_destroy failed g_p21',61, GA_ERR)
        
c        we will localize the doubly occupied orbitals, assuming a
c        closed-shell calculation and no designated core
c        corbitals. Eventually, this will need to be extended to a loop
c        over spins, for unrestricted SCF, and we may want to add
c        localization of the virtual space, too

        if (.not. geom_num_core(rtdb, geom, 'ddscf', ncore)) ncore = 0

        if (debug) write (luout,*) 'ncore = ',ncore

        ispin = 1
        
         do i = 1, nclosed(ispin)
            iloc(i) = i
         end do
         nloc = nclosed(ispin)

c        copy the relevant CMOs to g_mo
         
         if (.not. ga_create(MT_DBL, nbf, nloc, 'loc:g_mo',
     &     nbf, 0, g_mo)) call errquit(pname//': g_mo',0, GA_ERR)

         call ga_copy_patch('n',
     &     g_movecs(ispin), 1, nbf, 1, nloc,
     &     g_mo,            1, nbf, 1, nloc)

         if (debug) write (luout,*) 'movecs -> mo'         

c        create C-tilde from Appendix C of Knizia's IBO paper.  g_temp
c        holds P12 * P21; we won't need it after the next matrix
c        multiplication
         
         if (.not. ga_create(MT_DBL, nbf, nloc, 'loc:ctilde',
     &     nbf, 0, g_ctilde)) call errquit('loc_driver: sc',0, GA_ERR)

         call ga_dgemm('n', 'n', nbf, nloc, nbf, 
     &     1.0d0, g_temp, g_mo, 0.0d0, g_ctilde)

         if (.not. ga_destroy(g_temp))
     &     call errquit(pname//': ga_destroy failed g_temp',71, GA_ERR)

         if (debug) write (luout,*) 'pre sym_orth 1'
         call sym_orth(g_ctilde, g_smat, nbf, nloc)

         if (dbg>1) call ga_print(g_ctilde)

c        now assemble the matrix products to form matrix A from Appendix
c        C of the IBO paper. The result will be accumulated in array
c        g_iao which are the IAOs in basis 1. The shape of the array is
c        the same as P12

         if (.not.ga_duplicate(g_smat, g_tmp1, 'g_tmp1')) call
     &     errquit(pname//': ga_dup tmp1', 72, GA_ERR)
         
         if (.not.ga_duplicate(g_smat, g_tmp2, 'g_tmp2')) call
     &     errquit(pname//': ga_dup tmp2', 72, GA_ERR)

         if (.not.ga_duplicate(g_smat, g_temp, 'g_temp')) call
     &     errquit(pname//': ga_dup temp', 72, GA_ERR)

         if (.not.ga_duplicate(g_p12, g_iao, 'g_iao')) call
     &     errquit(pname//': ga_dup iao', 72, GA_ERR)


c        C * C(trans) 
         call ga_dgemm('n', 't', nbf, nbf, nloc, 
     &     1.0d0, g_mo, g_mo, 0.0d0, g_temp)

c        C * C(trans) * S1 -> tmp1
         call ga_dgemm('n', 'n', nbf, nbf, nbf, 
     &     1.0d0, g_temp, g_smat, 0.0d0, g_tmp1)

c        tilde-C * tilde-C(trans) 
         call ga_dgemm('n', 't', nbf, nbf, nloc, 
     &     1.0d0, g_ctilde, g_ctilde, 0.0d0, g_temp)

c        tilde-C * tilde-C(trans) * S1 -> tmp2
         call ga_dgemm('n', 'n', nbf, nbf, nbf, 
     &     1.0d0, g_temp, g_smat, 0.0d0, g_tmp2)         

c        C * C(trans) * S1  * tilde-C * tilde-C(trans) * S1

         call ga_dgemm('n', 'n', nbf, nbf, nbf, 
     &     1.0d0, g_tmp1, g_tmp2, 0.0d0, g_temp)

c        multiply prev. result by P12 and store in g_iao

         call ga_dgemm('n', 'n', nbf, nbf_iao, nbf, 
     &     1.0d0, g_temp, g_p12, 0.0d0, g_iao)
         
c        form C * C(trans) -1
         call ga_add_todiag(g_tmp1, -1.0d0)
         
c        form tilde-C * tilde-C(trans) -1
         call ga_add_todiag(g_tmp2, -1.0d0)

c        multiply the two prev. results with each other, then multiply
c        the result with P12 and add to g_iao
         
         call ga_dgemm('n', 'n', nbf, nbf, nbf, 
     &     1.0d0, g_tmp1, g_tmp2, 0.0d0, g_temp)

         call ga_dgemm('n', 'n', nbf, nbf_iao, nbf, 
     &     1.0d0, g_temp, g_p12, 1.0d0, g_iao)

c        g_ctilde no longer needed
         
         if (.not. ga_destroy(g_ctilde))
     &     call errquit(pname//': ga_destroy failed ctilde',62, GA_ERR)

c        delete temp arrays before calling sym_orth, which allocates its
c        own temp memory

         if (.not. ga_destroy(g_tmp1))
     &     call errquit(pname//': ga_destroy failed g_tmp1',71, GA_ERR)
         
         if (.not. ga_destroy(g_tmp2))
     &     call errquit(pname//': ga_destroy failed g_tmp2',71, GA_ERR)
         
         if (.not. ga_destroy(g_temp))
     &     call errquit(pname//': ga_destroy failed g_temp',71, GA_ERR)         
c        last step: orthonormalize g_iao

         if (debug) write (luout,*) 'pre sym_orth 2'
         call sym_orth(g_iao, g_smat, nbf, nbf_iao)

         call ga_sync
         
         if (dbg>1) call ga_print(g_iao)

c        check orthogonality of IAOs w.r.t. metric S1
c        g_p12 is used for scratch

        if (.not. ga_create(MT_DBL, nbf_iao, nbf_iao, 'g_temp',
     &    nbf_iao, 0, g_temp)) call errquit(pname//': temp',73, GA_ERR)
         
         call ga_dgemm('n', 'n', nbf, nbf_iao, nbf, 
     &     1.0d0, g_smat, g_iao, 0.0d0, g_p12)

         call ga_dgemm('t', 'n', nbf_iao, nbf_iao, nbf, 
     &     1.0d0, g_iao, g_p12, 0.0d0, g_temp)

         if (dbg>1) call ga_print(g_temp)

         call ga_add_todiag(g_temp, -1.0d0)

         call ga_maxelt(g_temp, rtemp)
         if (master) then
           write(luout,'(/1x,a,1x,f18.12)')
     &       'IBO loc: largest element in C(iao, trans) S C(iao) -1: ',
     &     rtemp
           write(luout,'(1x,a/1x,a/)')
     &       'Significant deviations from zero may indicate',
     &       'elevated numerical noise in the IAO generation'
         end if

         if (.not. ga_destroy(g_temp))
     &     call errquit(pname//': ga_destroy failed g_temp',73, GA_ERR)

c        g_p12 no longer needed
         
         if (.not. ga_destroy(g_p12))
     &     call errquit(pname//': ga_destroy failed g_p12',61, GA_ERR)
         
c        allocate some memory used in the localization routine:

         if (.not. ma_push_get(mt_dbl, 2*nbf_iao, 'c', l_c, k_c))
     &        call errquit(pname//': ma for c', 0, MA_ERR)


c        jochen: comment
c        IBO localization was modeled on the pre-existing PM
c        localization code
         call localizeIBO(basis, minbas,
     &     dbl_mb(k_c),
     &     nloc, iloc, nbf, nmo, g_mo, g_smat, g_iao, nbf_iao,
     &     dbl_mb(k_eval+(ispin-1)*nbf), dbl_mb(k_occ+(ispin-1)*nbf))

c        g_mo now holds the localized MOs in the AO basis

c        assuming occupied MO localization, we copy the LMOs into the
c        relevant part of movecs

         call ga_copy_patch('n',
     &     g_mo,        1, nbf, 1, nloc,
     &     g_movecs(1), 1, nbf, 1, nloc)


c        deallocate remaining arrays that are no longer needed
         
         if (.not. ma_pop_stack(l_c)) call errquit(
     &      pname//': error MA pop c',0, MA_ERR)

         if (.not. ga_destroy(g_iao)) call errquit(
     &     pname//': error destroying g_iao',0, GA_ERR)

         if (.not. ga_destroy(g_mo)) call errquit(
     &     pname//': error destroying g_mo',0, GA_ERR)         
         
c        smat not needed anymore
         if (.not. ga_destroy(g_smat)) call errquit(
     &     pname//': error destroying g_smat',0, GA_ERR)

c        destroy minimal basis (iao basis)
         if (.not.bas_destroy(minbas))
     &     call errquit(pname//'iao bas_destroy failed',0,BASIS_ERR)


       end if ! loctype

c      ================================================================

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
     &  10x,28(1h*),/,10x,'*** ORBITAL LOCALIZATION ***',/,
     &  10x,28(1h*)/)

 9020 format(/,10x,'Open shell calculations not yet implemented',
     &      ' for the Localization module. Exiting')
      end

c     =================================================================

      subroutine sym_orth(g_c, g_s, n1, n2)

      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"      
#include "numerical_constants.fh"      

c     Given input vectors collected in matrix C, and a metric S,
c     calculate
c     C' = C times [transpose(C) S C]**(-1/2) such that
c     transpose(C') S C = unit matrix

c     matrix S is square symmetric, dimension n1, n1
c     matrix C has dimension n1, n2
      
c     in the calling routine, matrix C may be allocated with more than
c     n2 columns, which is the reason we use ga_copy_patch below
      
      integer g_c !  in/out
      integer g_s !  input

      integer n1, n2  ! input

      integer g_sc, g_csc

      character*(8) pname

      logical debug

c     -----------------------------------------------------------------

      pname = 'sym_orth'

      debug = .false. .and. ga_nodeid().eq.0 !during development

      if (debug) write(luout,*) 'entering '//pname

       if(.not.ga_create(mt_dbl,n1,n2,'g_sc',0,0,g_sc))
     &  call errquit(pname//': error creating g_sc',141, GA_ERR)

       if(.not.ga_create(mt_dbl,n2,n2,'g_csc',0,0,g_csc))
     &    call errquit(pname//': error creating g_csc',142, GA_ERR) 

       call ga_dgemm('n','n',n1, n2, n1, 1.0d0,
     &   g_s, g_c, 0.0d0, g_sc)

       call ga_dgemm('t','n',n2, n2, n1, 1.0d0,
     &   g_c, g_sc, 0.0d0, g_csc)

       call ga_matpow(g_csc, -FP_Half, 1d-12)

       call ga_copy_patch('n', 
     &   g_c,  1, n1, 1, n2,
     &   g_sc, 1, n1, 1, n2)

       call ga_dgemm('n','n',n1, n2, n2, 1.0d0,
     &   g_sc, g_csc, 0.0d0, g_c)       


       if (.not. ga_destroy(g_sc))
     &   call errquit(pname//': error destroy g_sc',141, GA_ERR)

       if (.not. ga_destroy(g_csc))
     &   call errquit(pname//': error destroy g_csc',142, GA_ERR)       

       if (debug) write(luout,*) 'leaving '//pname
       
       return

       end

      
