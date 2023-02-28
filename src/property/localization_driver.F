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
#include "nwc_const.fh"

c     subroutine arguments:
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer loctype

c     local GA handles:
      integer g_uc(4), g_smat, g_sc, g_tran
      integer g_movecs(2), g_cmo(2), g_temp, g_tmp1, g_tmp2

c     MA variables:
      integer l_c, k_c, l_sc, k_sc, l_eval, k_eval, l_occ, k_occ
      integer l_dip(3), k_dip(3)
      integer l_pop, k_pop, l_list, k_list

c     other local variables:

      integer loc_opt

      integer nclosed(2), nopen(2), nvirt(2), ncore, nocc, nvir,
     &  ndens, nbf, nmo, nloc, i, natoms

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

      character*(19) pname
      character*(3) ltyp

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
     &       'Pipek-Mezey localization', 25, .true.)
        elseif (loctype.eq.2) then
          call util_print_centered(luout,
     &      'Boys localization', 25, .true.)
        elseif (loctype.eq.3) then
          call util_print_centered(luout,
     &      'IAO-IBO localization', 25, .true.)
        else
          write (luout,*) 'loctype =',loctype
          write (luout,*) 'LocType value not supported. Exiting.'
          return
        end if
      end if

c     retrieve localization option from rtdb if it exists

      if (.not. rtdb_get(rtdb, 'prop:loc_opt', mt_int, 1,
     &                               loc_opt  ))
     &  loc_opt = 0

c     option is currently only implemented for IBOs:

      if (loc_opt.gt.0 .and. (loctype.ne.3)) then
        if (oprint) write(luout,'(/1x,a/1x,a)')
     &    'localization option not implemented for chosen loc. type',
     &    'setting option to 0 ...'
        loc_opt = 0
      end if

      if (oprint) then
        if (loc_opt .eq.0) then
          call util_print_centered(luout,
     &      '(occupied orbitals)', 25, .false.)
        else if (loc_opt.eq.1) then
          call util_print_centered(luout,
     &      '(virtual orbitals)', 25, .false.)
        else if (loc_opt.eq.2) then
          call util_print_centered(luout,
     &      '(occ and virt orbitals)', 25, .false.)
        else
          call errquit(pname//': localization option out of range',
     &      loc_opt, INPUT_ERR)
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
     &  nopen(1), nvirt(1)

      natoms = 0

      if (.not. geom_ncent(geom, natoms)) call errquit
     &  (pname//': geom',0, GEOM_ERR)

      if (natoms.lt.1) call errquit
     &  (pname//': could not get natoms from rtdb', 911, UNKNOWN_ERR)

      if (natoms.gt.nw_max_atom) call errquit
     &  (pname//': nw_max_atom too small ', 911, UNKNOWN_ERR)

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
     &  errquit('loc_driver: maxnloc too small',0, BASIS_ERR)

c     for the time being, we set these variables here:

      ispin = 1
      nocc = nclosed(1)
      nvir = nvirt(ispin)


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


c      ===============================================================

       else if (loctype.eq.3) then

c        ===================================================
c        IAO-IBO localization (occupied or virtual orbitals)
c        ===================================================


         ispin = 1
         nocc = nclosed(ispin)
         nvir = nvirt(ispin)

c        allocate MA array for pairs of MOs. The actual dimensions
c        needed are less or equal to 2*nbf

         if (.not. ma_push_get(mt_dbl, 2*nbf, 'c mo', l_c, k_c))
     &     call errquit(pname//': ma for c', 0, MA_ERR)

         if (.not. ma_push_get(mt_dbl, natoms, 'pop', l_pop, k_pop))
     &     call errquit(pname//': loc:pop', 0, MA_ERR)

         if (.not. ma_push_get(mt_int, natoms, 'list', l_list, k_list))
     &     call errquit(pname//': loc:list', 0, MA_ERR)

         if (loc_opt.eq.0) then
           ltyp = 'occ'
           call ibo_localization(rtdb, geom, ltyp, basis,
     &       g_movecs(ispin),nocc,nvir, nmo, nbf, natoms,
     &       dbl_mb(k_eval+(ispin-1)*nbf),
     &       dbl_mb(k_occ+(ispin-1)*nbf),
     &       dbl_mb(k_c),
     &       dbl_mb(k_pop), int_mb(k_list))

         else if (loc_opt.eq.1) then
           ltyp = 'vir'
           call ibo_localization(rtdb, geom, ltyp, basis,
     &       g_movecs(ispin),nocc,nvir, nmo, nbf, natoms,
     &       dbl_mb(k_eval+(ispin-1)*nbf),
     &       dbl_mb(k_occ+(ispin-1)*nbf),
     &       dbl_mb(k_c),
     &       dbl_mb(k_pop), int_mb(k_list))

         else if(loc_opt.eq.2) then
           ltyp = 'occ'
           call ibo_localization(rtdb, geom, ltyp, basis,
     &       g_movecs(ispin),nocc,nvir, nmo, nbf, natoms,
     &       dbl_mb(k_eval+(ispin-1)*nbf),
     &       dbl_mb(k_occ+(ispin-1)*nbf),
     &       dbl_mb(k_c),
     &       dbl_mb(k_pop), int_mb(k_list))

           ltyp = 'vir'
           call ibo_localization(rtdb, geom, ltyp, basis,
     &       g_movecs(ispin),nocc,nvir, nmo, nbf, natoms,
     &       dbl_mb(k_eval+(ispin-1)*nbf),
     &       dbl_mb(k_occ+(ispin-1)*nbf),
     &       dbl_mb(k_c),
     &       dbl_mb(k_pop), int_mb(k_list))

         else
           call errquit(pname//': loc_opt out of range',loc_opt,
     &       UNKNOWN_ERR)

         end if ! loc_opt

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

 9020 format(/,10x,'Spin-unrestricted calculations not yet implemented',
     &      ' for the Localization module. Exiting')
      end
