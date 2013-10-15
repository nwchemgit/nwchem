      logical function scf(rtdb)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
#include "case.fh"
#include "frozemb.fh"
c
c     ROHF module.
c
c     Currently performs closed-shell and high-spin open-shell
c     restricted Hatree-Fock calculations
c
      integer rtdb              ! [input] database handle
      logical converged
      logical status, olocal
      double precision energy
      logical oprint
      character*255 name, ri_trans
      character*180 filename
c
      integer nbf_per_ir(20),i,j
      integer k,l,ij,kl,ijkl,ni,nj,nk,nl,nij,nkl,nijkl
      integer ism_lu,jsm_lu,ksm_lu,lsm_lu,smij,smkl,ijsm_lu,klsm_lu
      integer fn
      double precision int1e,int2e,enrep
      integer natoms, Nff, molo, mohi
      integer g_dens(2) ! dummy global array handles for initial guess
      integer maxireps
      parameter (maxireps=20)
      integer nop,nir,class_dim(maxireps)
      integer ipnt(maxireps),ipnt2(maxireps)
      character*8 zir(maxireps),zclass(maxireps)
      character*8 zname
      double precision chars(maxireps*maxireps)
      integer ioff(maxireps),icount(maxireps)
      integer irs,reord
c     integer k_reord,l_reord,k_reord2,l_reord2
      integer g_moreo,l_veci,k_veci
      integer g_hcore,noper,g_coul,g_exch
      logical rohf, uhf, riscf_init, movecs_converged
      external rohf, uhf, riscf_init, movecs_converged
      logical int_normalize
      external int_normalize
      logical cphf_poliz, do_poliz
      external cphf_poliz
      integer necp,itype
      integer blen
      integer minchunk
      parameter (minchunk = -1)
c
      integer ga_create
      external ga_create
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ga_create_JKblocked
      external ga_create_JKblocked
      logical util_io_unit
      external util_io_unit
c
      data blen/25/
c
c     Make sure database is parallel, 
c     and initialize performance statistics
c
      status = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'scf')
      call ecce_print_module_entry('scf')
      oprint = util_print('information', print_low)
c
c     Extract high level info from the data-base setting defaults
c     return parameters here and also shoving stuff into commons
c
c !!! BGJ
c !!! Can't place the following in scf_init since scf_init is called
c !!! by cphf for either HF or DFT
c     Store SCF hamiltonian type as HF for use in BGJ routines
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 1))
     $     call errquit('scf: put of bgj:scf_type failed',0, RTDB_ERR)
c !!! BGJ
      call scf_init(rtdb)
c
c     Initialize performance statistics
c
      call scf_pstat_init(rtdb)
c
c     Print info
c 
      if (ga_nodeid().eq.0 .and. oprint) then
         call util_print_centered(6, 'NWChem SCF Module', 40, .true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
c
         if(.not. geom_systype_get(geom,itype))
     &        itype = 0
         if(itype.ne.0) call errquit(
     A     'scf: code working only for molecules',itype, INPUT_ERR)

         if (.not. geom_ncent(geom, natoms)) call errquit
     $        ('scf: geom_ncent failed', 0, GEOM_ERR)
c
         write(6,1) trans(1:inp_strlen(trans)), nbf
 1       format(/
     $        '  ao basis        = "',a,'"'/
     $        '  functions       = ', i5)
         if (nriscf .ne. 0) then
            status = bas_numbf (riscf_basis, Nff)
            if (.not. bas_name(riscf_basis, name, ri_trans))
     $           call errquit('scf: bas_name?', 0,
     &       BASIS_ERR)
            write(6,11) ri_trans(1:inp_strlen(ri_trans)), nff
 11         format(
     $           '  ri basis        = "',a,'"'/
     $           '  ri functions    = ', i5)
         endif
         write(6,12) natoms
 12      format('  atoms           = ', i5)
         if (scftype .eq. 'UHF') then
            write(6,121) nalpha, nbeta
 121        format('  alpha electrons = ', i5/
     $           '  beta  electrons = ', i5)
         else
            write(6,122) nclosed, nopen
 122        format('  closed shells   = ', i5/
     $           '  open shells     = ', i5)
         endif
         write(6,123)
     $        charge, scftype,
     $        movecs_in(1:inp_strlen(movecs_in)),
     $        movecs_out(1:inp_strlen(movecs_out)), oskel, oadapt
 123     format(
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/
     $        '  input vectors   = ', a/
     $        '  output vectors  = ', a/
     $        '  use symmetry    = ', l1/
     $        '  symmetry adapt  = ', l1)
         if (olock) write(6,124) olock
 124     format('  lock orbitals   = ', l1)
         if (util_print('basis summary', print_default)) then
            if (.not.bas_summary_print(basis))
     &           call errquit('scf: basis summary print failed',911,
     &       BASIS_ERR)
         else
            write(6,*)
         endif
         if (oadapt .and. util_print('char table', print_high)) then
            call sym_print_char_table(geom)
         endif
         if (oadapt .and. util_print('basis sym', print_default)) then
            call sym_bas_irreps(basis, .true., nbf_per_ir)
         endif
         call util_flush(6)
      endif
c
c     Print out other info as requested by user
c
      if (ga_nodeid() .eq. 0) then
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $        call errquit('scf: geom_print', 0, GEOM_ERR)
            write(6,*)
         endif
         if (util_print('symmetry', print_debug)) then
            call sym_print_all(geom, .true., .true., .true., 
     $           .true., .true.)
            write(6,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('scf: bas_print', 0, BASIS_ERR)
            write(6,*)
         endif
         if (util_print('basis labels',print_high) .and.
     $        ga_nodeid().eq.0) then
            call bas_print_labels(basis)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(basis)) 
     $           call errquit('scf:gmp', 0, BASIS_ERR)
            write(6,*)
         endif
         call util_flush(6)
      endif
c
c     Get and check the Coulomb Attenuation Method (CAM) parameters
      if (.not.rtdb_get(rtdb, 'scf:cam_exch', mt_log, 1,
     &   cam_exch)) cam_exch=.false.
      if (.not.rtdb_get(rtdb, 'scf:cam_omega', mt_dbl, 1,
     &   cam_omega)) cam_omega=0.d0
      if (.not.rtdb_get(rtdb, 'scf:cam_alpha', mt_dbl, 1,
     &   cam_alpha)) cam_alpha=0.d0
      if (.not.rtdb_get(rtdb, 'scf:cam_beta', mt_dbl, 1,
     &   cam_beta)) cam_beta=0.d0
c
      if (.not.int_normalize(rtdb,basis))
     $           call errquit('scf:int_normalize failed', 0, INT_ERR)
c
c     If input and output movecs are same then might be converged MOs
c
      if (movecs_in .eq. movecs_out) then
         converged = movecs_converged(rtdb, basis, 'scf', movecs_in)
         if (converged) then
            if (.not. rtdb_get(rtdb, 'scf:energy', mt_dbl, 1, energy))
     $           call errquit('scf: failed getting converged energy',0,
     &       RTDB_ERR)
            if (ga_nodeid().eq.0 .and. oprint) then
               write(6,1101) energy
 1101          format(/'  The SCF is already converged '//,
     $              '         Total SCF energy =', f20.12/)
               call util_flush(6)
            endif
            goto 3131
         endif
      endif
c
c     initialize for RI and calculate 3c integrals if necessary
c
      if (nriscf.gt.0) then
        status = riscf_init(rtdb, geom, basis, riscf_basis, tol2e)
        if (.not. status) call errquit('scf: riscf_init failed', 0,
     &       UNKNOWN_ERR)
        if (ga_nodeid().eq.0 .and. oprint) then
          if (nriscf.eq.1) then
            call util_print_centered(6, 'RI hessian', 40, .true.)
          else if (nriscf.eq.2) then
            call util_print_centered(6, 'RISCF', 40, .true.)
          else if (nriscf.eq.3) then
            call util_print_centered(6, 'preconverge with RISCF', 40, 
     $           .true.)
          endif
          write(6,*)
        endif
      endif
c
c     Initialize integral library, schwarz screening, cached 1-e integs
c
      call int_init(rtdb,1, basis)
c
      call print_integrals(basis,oskel) ! print integrals, if requested
c 
      call schwarz_init(geom, basis)
      if ((nbf*nbf*3*8)/ga_nnodes() .lt. 2000000) ! Cache 1e if < 2MB/node
     $     call int_1e_cache_ga(basis, oskel)
c
c     Form intial guess vectors
c
      if (ga_nodeid().eq.0 .and. oprint) write(6,711) util_wallsec()
 711  format(/' Forming initial guess at ',f9.1,'s'/)
c      call scf_vectors_guess(rtdb)
c
c     prev call replaced with new call
c
      call scf_vectors_guess(rtdb, tol2e, geom, basis, trans, 
     &                       movecs_in, movecs_out, movecs_guess,
     &                       scftype, nclosed, nopen, nbf, nmo,
     &                       nalpha, nbeta,  k_eval, k_occ, k_irs,
     &                       g_movecs, g_dens, 'vectors', 
     &                       'scf', title, oskel, oadapt, 
     &                       oinitialized) 
c
c     Here branch to the appropriate wavefunction type ... returns
c     true if happily converged, false otherwise (e.g., if restart
c     is needed or problems encountered).
c
      if (ga_nodeid().eq.0 .and. oprint) write(6,712) util_wallsec()
 712  format(/' Starting SCF solution at ',f9.1,'s'/)
c
c     frozen embedding 
      frozemb = .false.
      if (.not.rtdb_get(rtdb, 'scf:frozemb', mt_log, 1, frozemb))
     &   frozemb=.false.
      if (frozemb) then
       g_frozemb = ga_create_atom_blocked(geom, basis, 'rohf: frozemb')
       call ga_zero(g_frozemb)
       call scf_frozemb(rtdb, geom, basis, nbf, g_frozemb)
      end if   ! frozemb
c
      if (scftype .eq. 'UHF') then
         converged = uhf(rtdb, energy)
      else
         converged = rohf(rtdb, energy)
      endif
c
c     reconverge with exact energy/gradient if requested
c
      if (nriscf.eq.3) then
        if (ga_nodeid().eq.0 .and. oprint) then
          call util_print_centered(6,
     $         'switching to exact energy/gradient', 40, .true.)
          write(6,*)
        endif
        nriscf = 1
        if (scftype .eq. 'UHF') then
          converged = uhf(rtdb, energy)
        else
          converged = rohf(rtdb, energy)
        endif
      endif
c
c     Store energy and convergence status ... not sure if better to
c     do this here or have each wavefunction type do it
c
      if (.not. rtdb_put(rtdb, 'scf:energy', MT_DBL, 1, energy))
     $     call errquit
     $     ('scf: failed to store energy in rtdb', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:converged', MT_LOG, 1, converged)) 
     $     call errquit
     $     ('scf: failed to store converged in rtdb', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:alpha irreps', MT_INT, nbf,
     &   int_mb(k_irs)))
     &   call errquit('scf: failed to MO irreps in rtdb', 0, RTDB_ERR)
      if (scftype.eq.'UHF') then
      if (.not. rtdb_put(rtdb, 'scf:beta irreps', MT_INT, nbf,
     &   int_mb(k_irs+nbf)))
     &   call errquit('scf: failed to MO irreps in rtdb', 0, RTDB_ERR)
      endif
c
c     Must store the final MOs AFTER writing the convergence flag
c     since the time it is put into the RTDB is used to verify
c     that the two correspond.
c
      if (.not. rtdb_get(rtdb, 'scf:localize', mt_log, 1, olocal))
     $     olocal = .false.
      if (olocal) CALL TESTLOCALIZE(rtdb)
c
      call scf_movecs_write(rtdb)
c
c     With some better data compression in the ECCE output
c     we can now output all occupied plus the top 20 virtuals.
c     Also, ECCE may currently be confused if molo != 1
c
      molo = 1
      mohi = nmo ! min(nmo,nalpha+20)
c
      if (scftype .eq. 'UHF') then
         call movecs_ecce(nbf, nmo, molo, mohi, dbl_mb(k_eval),
     $        dbl_mb(k_occ), int_mb(k_irs), 
     $        g_movecs(1), scftype, 'alpha')
         call movecs_ecce(nbf, nmo, molo, mohi, dbl_mb(k_eval+nbf),
     $        dbl_mb(k_occ+nbf), int_mb(k_irs+nmo),
     $        g_movecs(2), scftype, 'beta')
      else
         call movecs_ecce(nbf, nmo, molo, mohi, dbl_mb(k_eval),
     $        dbl_mb(k_occ), int_mb(k_irs),
     $        g_movecs(1), scftype, ' ')
      endif
c.. dongxia added June 2013, transform AO ints to MO ints
c.. 
c..  reorder 
       if(.not.ma_push_get(mt_int,nbf,'reord',l_reord,k_reord))
     $    call errquit('scf:memory?',nbf,MA_ERR)
c..
       if(.not.ma_push_get(mt_int,nbf,'reord2',l_reord2,k_reord2))
     $    call errquit('scf:memory?',nbf,MA_ERR)
c..
      if (ga_nodeid().eq.0) then
       call ifill(nir,0,ipnt,1)
       call ifill(nir,0,ipnt2,1)
       call ifill(maxireps,0,ioff,1)
       call ifill(maxireps,0,icount,1)
       call sym_group_name(geom,zname)
       if(.not.sym_mol(zname,nop,nir,class_dim,zir,zclass,chars,ipnt))
     $    call errquit('scf: no char table available ',geom,unknown_err)
       do i=1,nir
         ipnt2(ipnt(i))=i
         write(6,*)'MOLCAS pointer',ipnt(i),ipnt2(i)
       end do
       ioff(1)=1
       do i=2,nir
         ioff(i) = ioff(i-1)+nbf_per_ir(ipnt2(i-1))
       end do
       do i=1,nmo
         irs=ipnt(int_mb(k_irs+i-1))
         write(6,*)'irs',ga_nodeid(),irs
         icount(irs)=icount(irs)+1
         reord = ioff(irs)+icount(irs)-1
         int_mb(k_reord+i-1)=reord
         int_mb(K_reord2+reord-1)=i
       end do
         write(6,*)'reorder arrays:'
         do i=1,nmo
          write(6,*)int_mb(k_reord+i-1),int_mb(k_reord2+i-1)
         end do
       call util_flush(6)
      end if
c..
c.. reorder movecs and save to g_moreo
c..
        g_moreo = ga_create_atom_blocked(geom,basis,'reorder MO')
        call ga_zero(g_moreo)
        call ga_zero(g_moreo)
        if(.not. ma_push_get(mt_dbl,nbf,'reo',l_veci,k_veci))
     $    call errquit('scf reo: no scratch space',nbf,MA_ERR)
c..
      if(ga_nodeid().eq.0)then
        do i=1,nmo
         call ga_get(g_movecs(1),1,nbf,i,i,dbl_mb(k_veci),1)
         j=int_mb(k_reord+i-1)
         write(6,*)'new index:',j
         call ga_put(g_moreo,1,nbf,j,j,dbl_mb(k_veci),1)
        end do
      endif
        call ga_print(g_moreo)
c       call util_flush(6)
c..
        if(.not.ma_pop_stack(l_veci))
     $    call errquit('scf: ma pop failed',0,0)
c..
       if(.not.ma_pop_stack(l_reord2))
     $    call errquit('scf: ma pop failed',0,0)
c..
       if(.not.ma_pop_stack(l_reord))
     $    call errquit('scf: ma pop failed',0,0)
c..  open file to write
       fn = 90
       filename = 'moints'
       if(.not.util_io_unit(90,99,fn))
     &    call errquit("cannot get io unit",0,0)
       call util_file_name_resolve(filename,.false.)
       open(fn,file=filename,form='formatted',status='unknown')
       rewind (fn)
c!!! note!!!
c.. dimension. Need to find all the arguments used below. Not done yet.
c.. nsxfsm and ndxfsm are functions in lucia.f.
         nint1 = nsxfsm(nsmob,mxpobs,ntoobs,ntoobs,itssx,adsxa,1,iprnt)
         nint2 = ndxfsm(nsmob,nsmsx,mxpobs,ntoobs,ntoobs,ntoobs,ntoobs,
     &           itsdx,adsxa,sxdxsx,i12s,i34s,i1234s,iprnt)
c.. generate GA to store 1-e ints for LUCIA
         if(.not.ga_create(mt_dbl,1,nint1,'lucia_1e',minchunk,1,g_lu1e))
     &       call errquit('scf: cannot create lucia_1e',0,GA_ERR)
         if(.not.ga_create(mt_dbl,1,nint2,'lucia_2e',minchunk,1,g_lu2e))
     &       call errquit('scf: cannot create lucia_2e',0,GA_ERR)
c..  1-e integrals
         g_hcore = ga_create_atom_blocked(geom,basis, 'rohf: hcore')
         call ga_zero(g_hcore)
         call int_1e_ga(basis, basis, g_hcore, 'kinetic', oskel)
         call int_1e_ga(basis, basis, g_hcore, 'potential', oskel)
         call moints_1e(nbf, basis, g_moreo, g_hcore)
         call ga_print(g_hcore)
c..  2-e integrals
         noper = nbf*(nbf+1)/2
         g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb Oper')
         g_exch = ga_create_JKblocked(noper,nbf,nbf,'Exch Oper')
         call ga_zero(g_coul)
         call ga_zero(g_exch)
         call moints_build_2x(basis,.true.,oskel,1,1,nbf,1,nbf,
     $        g_moreo,g_coul,.true.,g_exch,.true.,blen,.true.)
c        call ga_print(g_coul)
         call util_flush(6)
c        call ga_print(g_exch)
c.. To write in LUCIA format
         if(ga_nodeid().eq.0)then
c.. 1e ints
         ind1e=0
         ind2e=0
         do ism_lu=1,nir
          ni = nbf_per_ir(ipnt2(ism_lu))
          nij = ni*(ni+1)/2
          do ij = 1,nij
           ind1e = ind1e + 1
           i = int(sqrt(2.0*ij)+0.5)
           j = ij - i*(i-1)/2
c.. add offset
           i = i + ioff(ism_lu) -1
           j = j + ioff(ism_lu) -1
         call ga_get(g_hcore,i,i,j,j,int1e,1)
c        write(6,'(E22.15)')int1e
c        write(fn,'(E22.15)')int1e
         call ga_put(g_lu1e,1,1,ind1e,ind1e,int1e,1)
          end do
         end do 
c.. 2e ints
         do ism_lu=1,nir
          do jsm_lu=1,ism_lu
           smij=(ism_lu-1)*nir+jsm_lu
           ijsm_lu=ieor(ism_lu-1,jsm_lu-1)+1
           do ksm_lu=1,ism_lu
            do lsm_lu=1,ksm_lu
             smkl=(ksm_lu-1)*nir+lsm_lu
             klsm_lu=ieor(ksm_lu-1,lsm_lu-1)+1
             if ((ijsm_lu.eq.klsm_lu).and.(smij.ge.smkl)) then
              ni = nbf_per_ir(ipnt2(ism_lu))
              nj = nbf_per_ir(ipnt2(jsm_lu))
              nk = nbf_per_ir(ipnt2(ksm_lu))
              nl = nbf_per_ir(ipnt2(lsm_lu))
              if(ism_lu.eq.jsm_lu)then
               nij = ni*(ni+1)/2
               nkl = nk*(nk+1)/2
              else
               nij = ni*nj
               nkl = nk*nl
              endif
              nijkl = nij*nkl
              if (smij.eq.smkl) nijkl = nij*(nij+1)/2  
              do ijkl =1, nijkl
               ind2e = ind2e + 1
               if(smij.eq.smkl)then
                 kl = nkl+1-int(sqrt(2.0*(nijkl-ijkl+1))+0.5)
                 ij = ijkl + kl*(kl-1)/2 -nij*(kl-1)       
               else
                 kl = int((ijkl-0.5)/nij)+1
                 ij = ijkl-(kl-1)*nij 
               end if
               if(ism_lu.eq.jsm_lu)then
                 k = int(sqrt(2.0*kl)+0.5) 
                 l = kl - k*(k-1)/2 
                 i = int(sqrt(2.0*ij)+0.5) 
                 j = ij - i*(i-1)/2 
               else
c.. tra_rou = 1
c                k = int((kl-0.5)/nl)+1
c                l = kl - (k-1)*nl 
c                i = int((ij-0.5)/nj)+1
c                j = ij - (i-1)*nj 
c.. tra_rou = 2
                 l = int((kl-0.5)/nk)+1
                 k = kl - (l-1)*nk 
                 j = int((ij-0.5)/ni)+1
                 i = ij - (j-1)*ni 
               end if      
c.. add on the offset for each irep
                 i = i+ioff(ism_lu)-1
                 j = j+ioff(jsm_lu)-1
                 k = k+ioff(ksm_lu)-1
                 l = l+ioff(lsm_lu)-1
c.. look for index in ga_coul
                  ij = (j-1)*nbf+i 
                  kl = k*(k-1)/2 + l
         call ga_get(g_coul,ij,ij,kl,kl,int2e,1)
c        write(6,'(E22.15)')int2e
c        write(fn,'(E22.15)')int2e
         call ga_put(g_lu2e,1,1,ind2e,ind2e,int2e,1)
               call util_flush(6)
              end do
             endif
            end do
           end do
          end do
         end do
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $  call errquit('unable to get nuclear rep energy',geom, GEOM_ERR)
         write(fn,'(E22.15)')enrep
         end if
         close(fn)
c..
c
c     terminate integrals 
c
      call int_1e_uncache_ga()
      call schwarz_tidy()
      call int_terminate
c 
c     clean up frozen embedding
      if (frozemb) then
        if (.not. ga_destroy(g_frozemb)) call errquit
     $     ('rohf_fock: ga_destroy failed',0, GA_ERR)
      end if
c
c     SCF is done ... tidy up memory etc.
c
      if (nriscf.ne.0) then
        call riscf_tidy ( geom, basis, riscf_basis )
      endif
 3131 call scf_tidy(rtdb)
c
      call scf_pstat_print
c
c     Nothing else for now ... return convergence status
c
      call util_print_pop
c
      if (converged) then
         call ecce_print_module_exit('scf', 'ok')
      else
         call ecce_print_module_exit('scf', 'failed')
      endif
c
      scf = converged
c
c !!! BGJ
      if (.not. rtdb_get(rtdb, 'bgj:poliz', mt_log,
     &     1, do_poliz)) then
         do_poliz = .false.
      endif
      if (do_poliz) then
         write(*,*)'*** scf: calling cphf_poliz'
         if (.not. cphf_poliz(rtdb)) ! Never executed.
     $        call errquit(' cphf_poliz: failed from scf !',0, RTDB_ERR)
      endif
c !!! BGJ
c
      end
