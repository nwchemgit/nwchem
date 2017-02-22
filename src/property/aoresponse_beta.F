
      subroutine aoresponse_beta (rtdb, basis, geom, nfreq,
     &   response_order, frequencies, g_smat0, g_smat1, g_sket1, g_dipel
     &   , g_quadel, g_dipmag, g_vectors, epst, froct, g_dens, nbf, nbfx
     &   , nmo, nocct, nvirt, deps, lgiao, lquad, lvelocity, lmagpert,
     &   lifetime, gamwidth, alfare, alfaim, betare, betaim)

c $Id$
      
c     =================================================================
      
c     purpose: perform response calculation of order "response_order"

c     ALTERNATE VERSION of aoresponse.F: CALCULATE BETA DIRECTLY using 
c     frequency dependent version of Amos equation as published
c     by J. Autschbach, Comp. Lett. 3 (2007), 131
c     This allows to compute the static limit directly instead of
c     using -G'/omega for small omega which is numerically unstable.
c     Moreover, we can use a static magnetic field which makes it
c     a little bit easier to implement GIAO optical rotation. 

c     called from: aoresponse_driver

c     output: alfare, alfaim - electric-electric response matrices
c             betare, betaim - electric-magnetic response matrices

c     NOTE: the GIAO functionality is not yet functioning for 
c           HF and hybrid DFT, only for pure DFT (Jan 15, 2011)

c     NOTE: With A-matrix it is meant the matrix that is used to
c           parameterize the response of the MO coefficients, i.e.
c           perturbed MO coeffs = C*A where C is the coefficient
c           matrix for the unperturbed MOs.
c           the matrix A can be one for linear, quadratic, etc.
c           response.

c     =================================================================

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

c     ---------------------
c     subroutine arguments:
c     ---------------------

      integer rtdb    ! [input] run-time database handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle

      integer g_smat0, g_dens(3), g_vectors(2), g_dipel, g_epsdif,
     &   g_dipmag, g_smat1, g_sket1, g_quadel
      integer nfreq, response_order, nbf, nbfx, nmo
      integer nocct(2), nvirt(2)
      double precision froct(nbf,2), epst(nbf,2)
      double precision frequencies(nfreq)
      double precision deps(nvirt(1),nocct(1))
      double precision alfare(3,3), alfaim(3,3)
      double precision betare(3,3), betaim(3,3)
      double precision gamwidth
      logical lgiao, lquad, lvelocity, lmagpert, lifetime


c     ----------------
c     local variables:
c     ----------------

c     global array handles:
      
      integer g_amat(2), 
     &   g_work, g_occ, g_rhs(2),
     &   g_rhs_beta, g_amat_beta, g_vecE1, g_vecB1, g_temp

c     GA arrays for treating imaginary parts of response
c     (only used if lifetime.eq.true.):

      integer g_amat_im(2), 
     &        g_rhs_im(2), g_rhs_beta_im, g_amat_beta_im,
     +        g_pmat_beta_im

c     other local variables: 

      integer nmot(2), nocvir(2)

      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      
      integer LCTensor(3,3,3)
      double precision origin(3)
      data origin/0.,0.,0./
      logical oskel
      parameter (oskel = .false.)

      character*(256) cstemp

      character*(1) direction(3)
      data direction/'x','y','z'/
      
      integer ispin, nspin, ncomp
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, iresp, idir
      logical debug, dbgmat, 
     &   lzora, lantisym, lstatic
      double precision omega, sum, scaling
      double precision tenm8, one, two, zero, half
      parameter (tenm8=1d-8, one=1d0, two=2d0, zero=0d0, half=one/two)

c     nwchem file access character strings:

      character*256 cphf_rhs(2), cphf_sol(2), cphf_rhs_im(2),
     &   cphf_sol_im(2)

c     external functions:

      double precision ga_trace_diag
      external ga_trace_diag

      logical  cphf2, cphf3, file_write_ga, file_read_ga, cphf
      external cphf2, cphf3, file_write_ga, file_read_ga, cphf

c  ====================================================================

      debug = .false. .and. ga_nodeid().eq.0 ! .true. during development
      dbgmat = .false. .and. ga_nodeid().eq.0 ! debug large matrices

      if (debug) write (luout,*) 'hello from aoresponse_beta'

c     make sure lvelocity.ne.T., we do not support that in this
c     subroutine to keep the clutter at a manageable level.
c     same for lmagpert

      if (lvelocity) 
     &   call errquit ('aor_beta: lvelocity set',1,INPUT_ERR)

      if (lmagpert) 
     &   call errquit ('aor_beta: lmagpert set',1,INPUT_ERR)

c     -------------------------
c     define Levi-Civita tensor for quadrupole additions
c     -------------------------
      LCtensor(:,:,:) = 0      
      LCtensor(1,2,3) = 1
      LCtensor(2,3,1) = 1
      LCtensor(3,1,2) = 1      
      LCtensor(2,1,3) = -1
      LCtensor(1,3,2) = -1
      LCtensor(3,2,1) = -1        

c     -------------------
c     determine frequency
c     -------------------

      omega = frequencies(response_order)
      lstatic = (abs(omega).lt.tenm8) ! static response or not
      ncomp = 2                 ! no. of Fourier components to treat
      if (lstatic) ncomp = 1    ! treat only one component for static
      
c     -------------------------
c     set-up for CPKS procedure
c     -------------------------

      if (ga_nodeid().eq.0) write (LuOut,'(a,i3,a/2x,a,E15.7,a/)')
     &   ' Performing order ',response_order, ' CPKS',
     &   ' with frequency omega = ',omega,' a.u.'
      if (lstatic .and. ga_nodeid().eq.0) write (LuOut,'(a/)')
     &   ' STATIC response' 

c     set parameters that control the various computational options
c     (later we will set most of this by input)
      nspin      =  1           ! assume closed shell
      lzora      = .false.      ! not yet available here 

c     check if we should calculate OR beta directly or G'

      if (ga_nodeid().eq.0) write (luout,*)
     & 'Optical Rotation: will calculate Beta directly, not G'''

      if (debug) write (luout,*) 'giao, velocity, magpert',
     &    lgiao, lvelocity, lmagpert   


c     -----------------------------------------
c     determine number of occ * virt orbitals
c     and nmot(1:2) and fix froct, if necessary
c     -----------------------------------------

      do ispin = 1,nspin
        nocvir(ispin) = nocct(ispin) * nvirt(ispin)
        nmot(ispin) = nmo
        if (nmo .lt.nbf) then
          do imo = nmo+1,nbf
            froct(imo,ispin) = 0d0
          enddo
        endif
      end do

c     ----------------------------------------------------
c     for each spin: load perturbing integrals into array
c     work, transform work to mo representation, store
c     vir-occ block in amat -> right hand side of CPKS
c     (it is NOT divided by (e_a - e_i -/+ omega), this
c     will be considered in the CPKS solver, in the precon-
c     ditioner and the 1e part of the "product" routine)
c     ----------------------------------------------------
      
      do ispin = 1, nspin
        
        nmo1 = nmot(ispin)      ! total no.of MOs for this spin
        nocc = nocct(ispin)     ! occupied MOs
        nvir = nvirt(ispin)     ! virtual MOs
        nocv = nocvir(ispin)    ! nocc * nvir

        chunk(1) = nbf
        chunk(2) = -1
        chunk(3) = -1
        dims(1) = nbf
        dims(2) = nbf
        dims(3) = 3 

c       ------------------------------
c       allocate some temp. work space
c       ------------------------------
        
        write(cstemp,'(a)') 'work'
        if (.not.nga_create(MT_DBL,3,dims,cstemp(1:4),chunk,
     &     g_work)) call 
     &     errquit('aoresponse: nga_create failed: '//cstemp(1:4),
     &     0,GA_ERR)     
        call ga_zero (g_work)
        
c       -----------------------------------------------------
c       GA-allocate components of A-matrices (+/-, Re and Im)
c       -----------------------------------------------------
        
        chunk(1) = nmo          ! distribution over processors
        chunk(2) = -1
        chunk(3) = -1
        dims(1) = nmo           ! dimensions of the array 
        dims(2) = nocc
        dims(3) = 3

c       allocate A-matrices for main response (dynamic, Re/Im)

        do ipm = 1,ncomp
          write(cstemp,'(a,i1)') 'amat_',ipm
          if (debug) write (luout,*) cstemp(1:6)
          if (.not.nga_create(MT_DBL,3,dims,cstemp(1:6),chunk,
     &       g_amat(ipm))) call 
     &       errquit('aoresponse: nga_create failed'//cstemp(1:6),
     &       0,GA_ERR)
          call ga_zero(g_amat(ipm))

          if (lifetime) then
            write(cstemp,'(a,i1)') 'amat_im_',ipm
            if (debug) write (luout,*) cstemp(1:9)
            if (.not.nga_create(MT_DBL,3,dims,cstemp(1:9),chunk,
     &         g_amat_im(ipm))) call 
     &         errquit('aoresponse: nga_create failed'//cstemp(1:9),
     &         0,GA_ERR)
            call ga_zero(g_amat_im(ipm))
          endif                 ! lifetime
        enddo                   ! ipm = 1,ncomp        


c       2. allocate additional A-matrices for the case that we want to 
c       calculate the OR parameter beta directly and we need 
c       the static response w.r.t. the B-field (same code as for
c       lzero.eq.true. in aoresponse.F)

       
        write(cstemp,'(a)') 'amat_beta'
        if (debug) write (luout,*) cstemp(1:9)
        if (.not.nga_create(MT_DBL,3,dims,cstemp(1:9),chunk,
     &     g_amat_beta)) call 
     &     errquit('aoresponse: nga_create failed'//cstemp(1:9),
     &     0,GA_ERR)
        call ga_zero(g_amat_beta)
        
        if (lifetime) then
          write(cstemp,'(a)') 'amat_zero_im'
          if (debug) write (luout,*) cstemp(1:12)
          if (.not.nga_create(MT_DBL,3,dims,cstemp(1:12),chunk,
     &       g_amat_beta_im)) call 
     &       errquit('aoresponse: nga_create failed'//cstemp(1:12),
     &       0,GA_ERR)
          call ga_zero(g_amat_beta_im)
        endif                   ! lifetime

        
c       -------------------------------------------------------
c       ipm counts the plus and minus combinations according to
c       Santry's notation. ipm=1: plus, ipm=2: minus
c       ipm = 1 corresponds to F(k,i)(+)
c       ipm = 2 corresponds to F(k,i)(-) = F(i,k)(+)*
c       NOTE: WE DON'T CONSIDER THE * HERE YET SO THAT MIGHT BE
C       A SOURCE OF ERROR FOR THE GIAO CASE
c       -------------------------------------------------------        
          
c       fill perturbing field integrals into work array   
        
        alo(1) = 1
        ahi(1) = nbf
        alo(2) = 1
        ahi(2) = nbf
        alo(3) = 1
        ahi(3) = 3
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nbf
        blo(3) = 1              ! 3 components of the
        bhi(3) = 3              ! perturbing field       
        
c       electric field: electric dipole integrals
        call nga_copy_patch('n',g_dipel,blo,bhi,g_work,alo,ahi)
    
        
c       ------------------------------------------------
c       transform to MO basis and extract vir-occ block,
c       store in g_amat's vir - occ block
c       ------------------------------------------------     
        
        call giao_aotomo(g_work,g_vectors(1),nocct(1),
     &     nvirt(1),1,3,nbf)   
        
        alo(1) = nocc+1
        ahi(1) = nmo
        alo(2) = 1
        ahi(2) = nocc
        alo(3) = 1
        ahi(3) = 3
        blo(1) = nocc + 1
        bhi(1) = nmo
        blo(2) = 1
        bhi(2) = nocc
        blo(3) = 1
        bhi(3) = 3
        
        do ipm = 1,ncomp
          
          call nga_copy_patch('n',g_work,alo, ahi,
     &       g_amat(ipm),blo,bhi)
          
        enddo                   ! ipm = 1,ncomp


c       need to redo the transformation and copy to A-matrix for
c       the B-field perturbed MOs. Same dimensioning as above.

        alo(1) = 1
        ahi(1) = nbf
        alo(2) = 1
        ahi(2) = nbf
        alo(3) = 1
        ahi(3) = 3
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nbf
        blo(3) = 1
        bhi(3) = 3
c       magnetic field: magnetic dipole or GIAO derivative of h
        call nga_copy_patch('n',g_dipmag,blo,bhi,g_work,alo,ahi)     
        call giao_aotomo(g_work,g_vectors(1),nocct(1),
     &     nvirt(1),1,3,nbf)             
        alo(1) = nocc+1 
        ahi(1) = nmo
        alo(2) = 1
        ahi(2) = nocc
        alo(3) = 1
        ahi(3) = 3
        blo(1) = nocc + 1
        bhi(1) = nmo
        blo(2) = 1
        bhi(2) = nocc
        blo(3) = 1
        bhi(3) = 3          
        call nga_copy_patch('n',g_work,alo, ahi,
     &     g_amat_beta,blo,bhi)  
        
c       call ga_print(g_amat_beta)
        
        
        if (debug) write (luout,*) 'amat ao2mo complete'
       
        
c       Create CPHF rhs arrays of proper dimension : (nocc*nvirt,3)
c       
c       TYPE = MT_DBL
c       nvir*nocc = dimension 1
c       3 = dimension 2, number of perturbing field components
c       -1 = equal distribution over processors
c       g_rhs = handle name (integer)
        do ipm = 1,ncomp
          write(cstemp,'(a,i1)') 'cphf_rhs_',ipm
          if(.not.ga_create(MT_DBL,nvir*nocc,3,cstemp(1:10),
     &       -1,-1,g_rhs(ipm)))
     &       call errquit
     &       ('aoresponse: ga_create failed '//cstemp(1:10),0,GA_ERR)
          call ga_zero(g_rhs(ipm))
          if (lifetime) then
            write(cstemp,'(a,i1)') 'cphf_rhs_im_',ipm
            if(.not.ga_create(MT_DBL,nvir*nocc,3,cstemp(1:13),
     &         -1,-1,g_rhs_im(ipm)))
     &         call errquit
     &         ('aoresponse: ga_create failed '//cstemp(1:13),0,GA_ERR)
            call ga_zero(g_rhs_im(ipm))
          endif                 ! lifetime
        end do

c       we also need the corresponding rhs matrices for the B-field
c       perturbed MOs

        if(.not.ga_create(MT_DBL,nvir*nocc,3,'rhs_beta',
     &     -1,-1,g_rhs_beta))
     &     call errquit
     &     ('aoresponse: ga_create failed rhs_beta',0,GA_ERR)
        call ga_zero(g_rhs_beta)
        if (lifetime) then
          if(.not.ga_create(MT_DBL,nvir*nocc,3,'rhs_beta_im',
     &       -1,-1,g_rhs_beta_im))
     &       call errquit
     &       ('aoresponse: ga_create failed rhs_beta_im',0,GA_ERR)
          call ga_zero(g_rhs_beta_im)
        endif

c       
c       In a GIAO computation with the magnetic-field perturbing,
c       first evaluate the constant
c       part of the rhs other than the derivatives of the one-
c       electron part of the Fock operator (now in amat).
c       These components are the same for ipm = 1,2 so we just
c       do this once and duplicate, if necessary.
c
c       NOTE: THIS IS NOT YET DEBUGGED AT ALL
c
c       Update 01/11: we probably need a call like the one below 
c       for the static magnetic (B) field perturbation if 
c       lgiao is set
c       
c$$$        if (lmagpert .and. lgiao) then
c$$$          call aoresponse_giao_rhs (rtdb, basis, geom, ncomp, g_rhs,
c$$$     &       g_vectors(ispin), g_smat1, g_sket1, g_dens(1), g_amat,
c$$$     &       nocct(ispin), nvirt(ispin), nmot(ispin), nbf,
c$$$     &       froct(1,ispin), epst(1,ispin), omega, lstatic)
c$$$        endif

c       01/11: the code above is disabled because we assume an electric field
c       perturbation for amat. However, there might be GIAO terms in
c       the magnetic field perturbation. We call it with lstatic = .true.
c       and a single component and zero frequency for the "beta" terms

        if (lgiao) then
          call aoresponse_giao_rhs (rtdb, basis, geom, 1, g_rhs_beta,
     &       g_vectors(ispin), g_smat1, g_sket1, g_dens(1), g_amat_beta,
     &       nocct(ispin), nvirt(ispin), nmot(ispin), nbf,
     &       froct(1,ispin), epst(1,ispin), 0d0, .true.)
        end if! lgiao
        
c       add amat vir-occ block to CPKS rhs, regular A-matrices first
        alo(1) = nocc + 1
        ahi(1) = nmo
        alo(2) = 1
        ahi(2) = nocc
        alo(3) = 1
        ahi(3) = 3
        blo(1) = 1
        bhi(1) = nocc*nvir
        blo(2) = 1
        bhi(2) = 3
        do ipm = 1,ncomp
          call nga_add_patch(1d0,g_amat(ipm),alo,ahi,
     &       1d0, g_rhs(ipm),blo,bhi, g_rhs(ipm), blo, bhi)
        end do

c       set up rhs matrices for opposite perturbation, zero freq
c       (for OR beta)
    
        call nga_add_patch(1d0,g_amat_beta,alo,ahi,
     &     1d0, g_rhs_beta,blo,bhi, g_rhs_beta, blo, bhi)
        

c       we don't need this work array for a while
        
        if (.not.ga_destroy(g_work))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_work',
     &     0,GA_ERR)
        
c       there is a division by -4 somewhere in the
c       CPKS solver so we scale the rhs by -4 in order to
c       get the correct magnitude of the result back
        
        do ipm = 1,ncomp
          call ga_scale (g_rhs(ipm), -4.0d0)
        enddo

        call ga_scale (g_rhs_beta, -4.0d0)
 
c       
c       Write ga_rhs to disk (+/- and Re/Im if applicable)
c         
        call cphf_fname('cphf_rhs',cphf_rhs(1))
        call cphf_fname('cphf_sol',cphf_sol(1))
        if(.not.file_write_ga(cphf_rhs(1),g_rhs(1))) call errquit
     $     ('aoresponse: could not write cphf_rhs 1',0, DISK_ERR)

        if (ncomp.gt.1) then
        call cphf_fname('cphf_rhs_2',cphf_rhs(2))
        call cphf_fname('cphf_sol_2',cphf_sol(2))
          if(.not.file_write_ga(cphf_rhs(2),g_rhs(2))) call errquit
     $       ('aoresponse: could not write cphf_rhs 2',0, DISK_ERR)
        endif

        if (lifetime) then
c
c         mainly, there are only zeros written to file here. i.e.
c         we start up the CPKS with zero imaginary parts (the
c         perturbation itself is always treated as real). 
c         
          call cphf_fname('cphf_rhs_im',cphf_rhs_im(1))
          call cphf_fname('cphf_sol_im',cphf_sol_im(1))
          if(.not.file_write_ga(cphf_rhs_im(1),g_rhs_im(1)))
     +       call errquit
     $       ('aoresponse: could not write cphf_rhs_im 1',0, DISK_ERR)
          
          if (ncomp.gt.1) then
            call cphf_fname('cphf_rhs_im_2',cphf_rhs_im(2))
            call cphf_fname('cphf_sol_im_2',cphf_sol_im(2))
            if(.not.file_write_ga(cphf_rhs_im(2),g_rhs_im(2)))
     +         call errquit
     $         ('aoresponse: could not write cphf_rhs_im 2',0, DISK_ERR)
          endif
        endif                   ! lifetime

        if (dbgmat) call ga_print(g_rhs(1))
        
c       next step is necessary is we use nwchem CPKS solver
        call schwarz_tidy()
        call int_terminate()
        
c       We do need to tell the CPHF that the density matrix
c       is real for the electric field perturbation.
c       the dynamic solver handles the +/- components of Pmat explicitly.
c       Done via rtdb, put cphf:skew .true. on rtdb and later
c       remove it. This comes from hnd_giaox.F where for some reason 
c       that I don't understand cphf:skew = .false. means that
c       the density matrix *IS* skew symmetric...
c       Here, we assume an electric field, though
c      
        if (.not.rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.true.)) call 
     &       errquit('aoresponse: failed to write skew ', 0, RTDB_ERR) 

c       we also need to tell the dynamic CPKS solver if the perturbation
c       is imaginary (it is not) and how many components to treat
        
        if (.not. rtdb_put(rtdb, 'cphf:imag', mt_log, 1, .false.)) call
     &     errquit('aoresponse: failed to write limag ', 0, RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'cphf:ncomp', mt_int, 1, ncomp)) call
     &     errquit('aoresponse: failed to write ncomp ', 0, RTDB_ERR)        
c       
c       ----------------------------------------------------
c       Call the CPKS solver. It will solve the sets of
c       equations for each component of the perturbing field
c       simultaneously. Dynamic CPKS for electric field
c       ----------------------------------------------------
c       
        if (debug) write (luout,*) 'aoresponse: calling CPKS' 

        if (lstatic .and. .not.lifetime) then
c         ... static CPKS solver:         
          if (.not.cphf2(rtdb)) call errquit
     &       ('aoresponse: failure in cphf2 ',0, RTDB_ERR)
        else       
c         ... dynamic CPKS solver:
          if (debug) write (luout,*)
     &       'calling cphf3',omega,lifetime,gamwidth
          if (.not.cphf3(rtdb, omega, lifetime, gamwidth))
     +       call errquit
     $       ('aoresponse: failure in cphf3 ',0, RTDB_ERR)
        endif

c       ---------------------------
c       CPKS done. Retrieve results
c       ---------------------------
        
        if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $     errquit('aoresponse: failed to delete skew', 0, RTDB_ERR)
        if (.not. rtdb_delete(rtdb, 'cphf:imag')) call
     &     errquit('aoresponse: failed to delete limag', 0, RTDB_ERR)
        if (.not. rtdb_delete(rtdb, 'cphf:ncomp')) call
     &     errquit('aoresponse: failed to delete ncomp', 0, RTDB_ERR)

        if (debug) write (6,*) 'aoresponse: back from CPKS'
        
c       Occ-virt blocks are the solution pieces of the CPHF
c       Read solution vector from disk and put solutions in amat's
c       vir - occ block. Any existing occ-occ block (GIAO) is still
c       intact from the call to aoresponse_giao_rhs
c       
        do ipm = 1,ncomp
          call ga_zero(g_rhs(ipm))
          if(.not.file_read_ga(cphf_sol(ipm),g_rhs(ipm)))
     &       call errquit
     $       ('aoresponse: could not read cphf solution',
     &       ipm, DISK_ERR)
          if (dbgmat) call ga_print(g_rhs(ipm))
          if (lifetime) then
            if(.not.file_read_ga(cphf_sol_im(ipm),g_rhs_im(ipm)))
     &         call errquit
     $         ('aoresponse: could not read cphf solution Im',
     &         ipm, DISK_ERR)
            if (dbgmat) call ga_print(g_rhs_im(ipm))
          endif                 ! lifetime
        enddo
        
c       alo, ahi and blo, bhi are the same as above when g_amat
c       was copied into g_rhs

        do ipm = 1,ncomp
          call nga_copy_patch('n',g_rhs(ipm),blo,bhi,
     &       g_amat(ipm),alo,ahi)
          if (dbgmat) call nga_print_patch(g_amat(ipm),alo,ahi,1)
          if (lifetime) then
            call nga_copy_patch('n',g_rhs_im(ipm),blo,bhi,
     &         g_amat_im(ipm),alo,ahi)
          endif                 ! lifetime
        enddo
        
c       deallocate CPKS right hand side vector and delete CPKS
c       scratch files
        
        do ipm = 1,ncomp
          write(cstemp,'(a,i1)') 'cphf_rhs_',ipm
          if (.not.ga_destroy(g_rhs(ipm))) call 
     &       errquit('hnd_giaox: ga_destroy failed '//cstemp(1:10),
     &       0,GA_ERR)
          if (lifetime) then
            write(cstemp,'(a,i1)') 'cphf_rhs_im_',ipm
            if (.not.ga_destroy(g_rhs_im(ipm))) call 
     &         errquit('hnd_giaox: ga_destroy failed '//cstemp(1:13),
     &         0,GA_ERR)
          endif
        enddo
        call util_file_unlink(cphf_rhs(1))
        call util_file_unlink(cphf_sol(1))
        if (lifetime) then
          call util_file_unlink(cphf_rhs_im(1))
          call util_file_unlink(cphf_sol_im(1))
        endif
        if (ncomp.gt.1) then
          call util_file_unlink(cphf_rhs(2))
          call util_file_unlink(cphf_sol(2))
          if (lifetime) then
            call util_file_unlink(cphf_rhs_im(2))
            call util_file_unlink(cphf_sol_im(2))
          endif
        endif

c       --------------------------------------------------------------
c       redo the CPKS for the B-field, i.e. we need the CPKS solved
c       for the second perturbation at omega = 0. The procedure
c       is essentially the same as for the lzero case in
c       aoresponse.F.
c
c       *** jochen 01/11: we might have to double check if the results
c       are affected by the GIAO occ-occ block !!! ***
c       --------------------------------------------------------------

        call cphf_fname('cphf_rhs',cphf_rhs(1))
        call cphf_fname('cphf_sol',cphf_sol(1))
        if(.not.file_write_ga(cphf_rhs(1),g_rhs_beta)) call errquit
     $     ('aoresponse: could not write cphf_rhs beta',0, DISK_ERR)
        
        if (lifetime) then
          call cphf_fname('cphf_rhs_im',cphf_rhs_im(1))
          call cphf_fname('cphf_sol_im',cphf_sol_im(1))
          if(.not.file_write_ga(cphf_rhs_im(1),g_rhs_beta_im))
     &       call errquit
     $       ('aoresponse: could not write cphf_rhs beta Im',0,
     &       DISK_ERR)
        endif                   ! lifetime
        
c       The B-field perturbation is imaginary  and antisymmetric) 
        
        if (.not.rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.))
     +     call 
     &     errquit('aoresponse: failed to write skew', 0, RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'cphf:imag', mt_log, 1, .true.))
     +     call
     &     errquit('aoresponse: failed to write limag', 0, RTDB_ERR)
        
        
c       also, we need to let the dynamic solver know that this is
c       a one-component situation (needed if lifetime.eq.true.):
        
        if (.not.rtdb_put(rtdb, 'cphf:ncomp', mt_int, 1,1)) call 
     &     errquit('aoresponse: failed to write ncomp ', 0, RTDB_ERR)
        
c       ---------------------------
c       run static CPKS for B-field
c       ---------------------------
        
c       call static CPKS solver if no lifetime key is present,
c       otherwise use dynamic solver but set omega to zero:
        
        if (.not.lifetime) then
          if (.not.cphf2(rtdb)) call errquit
     &       ('aoresponse: failure in cphf2 beta',0, RTDB_ERR)
        else
          if (.not.cphf3(rtdb, 0d0, lifetime, gamwidth)) call errquit
     &       ('aoresponse: failure in cphf3 beta',0, RTDB_ERR)
        endif
        
c       ---------------------------
c       CPKS done. Retrieve results
c       ---------------------------
        
        if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $     errquit('aoresponse: failed to delete skew', 0, RTDB_ERR)
        
        call ga_zero(g_rhs_beta)
        if(.not.file_read_ga(cphf_sol(1),g_rhs_beta))
     &     call errquit
     $     ('aoresponse: could not read cphf -beta- solution',
     &     0, DISK_ERR)          
        call nga_copy_patch('n',g_rhs_beta,blo,bhi,
     &     g_amat_beta,alo,ahi)
        
c       
c       call ga_print(g_amat_beta)
        
        if (lifetime) then
          call ga_zero(g_rhs_beta_im)
          if(.not.file_read_ga(cphf_sol_im(1),g_rhs_beta_im))
     &       call errquit
     $       ('aoresponse: could not read cphf -beta- solution Im',
     &       0, DISK_ERR)          
          call nga_copy_patch('n',g_rhs_beta_im,blo,bhi,
     &       g_amat_beta_im,alo,ahi)
        endif                   ! lifetime
        
        write(cstemp,'(a)') 'cphf_rhs_beta'
        if (.not.ga_destroy(g_rhs_beta)) call 
     &     errquit('aoresponse: ga_destroy failed '//cstemp(1:13),
     &     0,GA_ERR)
        
        if (lifetime) then
          write(cstemp,'(a)') 'cphf_rhs_beta_im'
          if (.not.ga_destroy(g_rhs_beta_im)) call 
     &       errquit('aoresponse: ga_destroy failed '//cstemp(1:15),
     &       0,GA_ERR)
        endif
        
        call util_file_unlink(cphf_rhs(1))
        call util_file_unlink(cphf_sol(1))
        if (lifetime) then
          call util_file_unlink(cphf_rhs_im(1))
          call util_file_unlink(cphf_sol_im(1))
        endif

        
        call int_init(rtdb,1,basis)
        call schwarz_init(geom,basis)
        call hnd_giao_init(basis,1)
        
        call ga_sync()          ! when do we actually need that ???
        if (debug) write (luout,*) 'CPKS complete, results in g_amat'

c       ---------------------------------------------------------
c       solution of CPKS is now in g_amat, including the division
c       by orbital energy diffs -/+ frequency.
c       Now allocate some work space again:
c       we will use g_rhs_* as scratch space to store the field
c       components of amat for calculation of the perturbed
c       density matrices:
c       ---------------------------------------------------------

        chunk(1) = nbf
        chunk(2) = -1
        dims(1) = nbf
        dims(2) = nbf        
        write(cstemp,'(a)') 'work'
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:4),chunk,
     &     g_work)) call 
     &     errquit('aoresponse: nga_create failed: '//cstemp(1:4),
     &     0,GA_ERR)     
        call ga_zero (g_work)


        ! storage for electric field A-matrices:
        ! nvir * nocc

        dims(1) = nvir
        dims(2) = nocc
        chunk(1) = dims(1)
        chunk(2) = -1
        do ipm = 1,ncomp
          write(cstemp,'(a,i1)') 'g_rhs_',ipm
          if (debug) write (luout,*) cstemp(1:7)
          if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &       g_rhs(ipm))) call 
     &       errquit('aoresponse: nga_create failed'//cstemp(1:7),
     &       0,GA_ERR)
          call ga_zero (g_rhs(ipm))          
          if (lifetime) then
            write(cstemp,'(a,i1)') 'g_rhs_im_',ipm
            if (debug) write (luout,*) cstemp(1:10)
            if (.not.nga_create(MT_DBL,2,dims,cstemp(1:10),chunk,
     &         g_rhs_im(ipm))) call 
     &         errquit('aoresponse: nga_create failed'//cstemp(1:10),
     &         0,GA_ERR)
            call ga_zero (g_rhs_im(ipm))
          endif                 ! lifetime
        enddo

        ! storage for magnetic-field A-matrices:
        ! nmo * nocc for giao, otherwise nvir * nocc

        dims(1) = nvir
        if (lgiao) dims(1) = nmo
        dims(2) = nocc
        chunk(1) = dims(1)
        chunk(2) = -1

        write(cstemp,'(a)') 'g_rhs_beta'
        if (debug) write (luout,*) cstemp(1:10)
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:10),chunk,
     &     g_rhs_beta)) call 
     &     errquit('aoresponse: nga_create failed'//cstemp(1:10),
     &     0,GA_ERR)
        call ga_zero (g_rhs_beta)          
        if (lifetime) then
          write(cstemp,'(a)') 'g_rhs_beta_im'
          if (debug) write (luout,*) cstemp(1:13)
          if (.not.nga_create(MT_DBL,2,dims,cstemp(1:13),chunk,
     &       g_rhs_beta_im)) call 
     &       errquit('aoresponse: nga_create failed'//cstemp(1:13),
     &       0,GA_ERR)
          call ga_zero (g_rhs_beta_im)
        endif                   ! lifetime
        
        
c       ------------------------------------------------------
c       memory allocation complete, now proceed to calculating
c       the perturbed density matrices
c       ------------------------------------------------------
        
        do idir = 1,3           ! direction of the perturbing field
          

          ! electric-field perturbed MOs:

          alo(1) = nocc + 1
          ahi(1) = nmo
          alo(2) = 1
          ahi(2) = nocc
          alo(3) = idir
          ahi(3) = idir
          blo(1) = 1
          bhi(1) = nmo - (alo(1) -1)
          blo(2) = 1
          bhi(2) = nocc
          
c         copy this component of A-matrix to g_rhs* :
          
          do ipm = 1,ncomp
            call nga_copy_patch ('n',g_amat(ipm), alo, ahi,
     &         g_rhs(ipm), blo, bhi)
            if (lifetime) then
              call nga_copy_patch ('n',g_amat_im(ipm), alo, ahi,
     &           g_rhs_im(ipm), blo, bhi)
            endif
          enddo

c         -------------------------------------------------
c         calculate electric-field perturbed MOs -> g_vecE1
c         -------------------------------------------------

c         create g_vecE1:

          write (cstemp,'(a)') 'aor:MO-pert efield'
          if(.not.ga_create(MT_DBL,nbf, nocc, trim(cstemp),
     &       -1,-1,g_vecE1))
     &       call errquit (trim(cstemp),0,GA_ERR)
          call ga_zero(g_vecE1)

c         we need to add the Fourier components for the E-field pert.:

          if (ncomp.gt.1) then
            call ga_add(1d0, g_rhs(1), 1d0,  g_rhs(2),
     &         g_rhs(1))
          endif

          alo(1) = 1
          ahi(1) = nbf
          alo(2) = nocc+1
          ahi(2) = nmo
          blo(1) = 1
          bhi(1) = nmo - (alo(2) -1)
          blo(2) = 1
          bhi(2) = nocc
          clo(1) = 1
          chi(1) = nbf
          clo(2) = 1
          chi(2) = nocc
            
c         perturbed MOs: MOs * A-matrix

          call ga_zero(g_work)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_vectors(ispin),alo,ahi,
     &       g_rhs(1),blo,bhi,
     &       g_vecE1,clo,chi)

          if (debug) write (luout,*) 'g_vecE1 created'

          
c         restore +/- components of A-matrix
          if (ncomp.gt.1)
     &       call ga_add(1d0, g_rhs(1), -1d0,  g_rhs(2), g_rhs(1))
          
c         ----------------------------------------------------------
c... DEBUG: calculate beta tensor from results
c         in the MO basis AT PRESENT ONLY FOR ELECTRIC FIELD PERT.
c         ----------------------------------------------------------

c         update 01/11: changed this to calculate this in the
c         AO basis set instead from C(E+/-) S C(B). 
c         We can also calculate the polarizability from
c         C(E+/-) S(0) C(0) for debug purposes
          
c         start loop over responding field direction
          
          do iresp = 1,3
            
c           step 2: copy amat for B-field perturbation into rhs_beta
c           dimensioning: nmo * nocc with GIAOs, otherwise nvir * nocc
            
            alo(1) = nocc + 1
            if (lgiao) alo(1) = 1 
            ahi(1) = nmo
            alo(2) = 1
            ahi(2) = nocc
            alo(3) = iresp
            ahi(3) = iresp
            blo(1) = 1
            bhi(1) = nmo - (alo(1) -1)
            blo(2) = 1
            bhi(2) = nocc
            
            call nga_copy_patch ('n',g_amat_beta, alo, ahi,
     &         g_rhs_beta, blo, bhi)
            
            if (lifetime)
     &         call nga_copy_patch ('n',g_amat_beta_im, alo, ahi,
     &         g_rhs_beta_im, blo, bhi)

c         -------------------------------------------------
c         calculate magnetic-field perturbed MOs -> g_vecB1
c         -------------------------------------------------

c         create g_vecB1:

          write (cstemp,'(a)') 'aor_beta :MO-pert bfield'
          if(.not.ga_create(MT_DBL,nbf, nocc, trim(cstemp),
     &       -1,-1,g_vecB1))
     &       call errquit (trim(cstemp),0,GA_ERR)
          call ga_zero(g_vecB1)

          alo(1) = 1
          ahi(1) = nbf
          alo(2) = nocc+1
          if (lgiao) alo(2) = 1 ! full set of vectors if we have GIAOs
          ahi(2) = nmo
          blo(1) = 1
          bhi(1) = nmo - (alo(2) -1)
          blo(2) = 1
          bhi(2) = nocc
          clo(1) = 1
          chi(1) = nbf
          clo(2) = 1
          chi(2) = nocc
            
c         perturbed MOs: MOs * A-matrix

          call ga_zero(g_vecB1)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_vectors(ispin),alo,ahi,
     &       g_rhs_beta,blo,bhi,
     &       g_vecB1,clo,chi)


c         it was confirmed for H2O2 RHF/STO-3G
c         that this MO vector is identical with that
c         produced by hnd_giaox for NMR shielding

          if (debug) write (luout,*) 'g_vecB1 created, direction ',iresp
c          if (debug) call ga_print(g_vecB1)


c         ------------------------------------------------------
c         allocate intermediate vector for matrix multiplication
c         that create the final results
c         ------------------------------------------------------
c         
          write (cstemp,'(a)') 'aor_beta: temp1'
          if(.not.ga_create(MT_DBL, nbf, nocc, trim(cstemp),
     &       -1,-1,g_temp))
     &       call errquit (trim(cstemp),0,GA_ERR)
          if (debug) write (luout,*) 'g_temp allocated'
          
c         ----------------------------------
c         calculate beta from C(E) S(0) C(B)
c         ----------------------------------
          
          alo(1) = 1
          ahi(1) = nbf
          alo(2) = 1
          ahi(2) = nbf
          alo(3) = 1 ! smat0 was created with 3 components. fix that
          ahi(3) = 1
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nocc
          clo(1) = 1
          chi(1) = nbf
          clo(2) = 1
          chi(2) = nocc
            

          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_smat0,alo,ahi,
     &       g_vecB1,blo,bhi,
     &       g_temp,clo,chi)

          if (debug) write (luout,*)
     &       'beta: S(0) C(B) intermediate complete'


          alo(1) = 1
          ahi(1) = nocc
          alo(2) = 1
          ahi(2) = nbf
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nocc
          clo(1) = 1
          chi(1) = nocc
          clo(2) = 1
          chi(2) = nocc

          call ga_zero(g_work)
          call nga_matmul_patch('t','n',1d0,0d0,
     &       g_vecE1,alo,ahi,
     &       g_temp,blo,bhi,
     &       g_work,clo,chi)

          if (debug) write (luout,*) 'beta: C(E) S(0) C(B) complete'

          
c         the factor of two is for the orbital occupations,
c         assuming that ispin is never equal to two

          sum = 2d0 * ga_trace_diag(g_work)

          betare(idir,iresp) = sum
          betaim(idir,iresp) = 0d0


          if (lgiao) then

c           --------------------------------------
c           if we use GIAOs there is a second term
c           in beta which is C(E) S(1ket) C(0)
c           --------------------------------------
            
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = iresp      ! pick the correct sket1 direction
            ahi(3) = iresp
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nocc
            clo(1) = 1
            chi(1) = nbf
            clo(2) = 1
            chi(2) = nocc
            
            
            call ga_zero(g_temp)
            call nga_matmul_patch('n','n',1d0,0d0,
     &         g_sket1,alo,ahi,
     &         g_vectors(ispin),blo,bhi,
     &         g_temp,clo,chi)
            
            if (debug) write (luout,*)
     &         'beta: S(ket1) C(0) intermediate complete'
            
            
            alo(1) = 1
            ahi(1) = nocc
            alo(2) = 1
            ahi(2) = nbf
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nocc
            clo(1) = 1
            chi(1) = nocc
            clo(2) = 1
            chi(2) = nocc
            
            call ga_zero(g_work)
            call nga_matmul_patch('t','n',1d0,0d0,
     &         g_vecE1,alo,ahi,
     &         g_temp,blo,bhi,
     &         g_work,clo,chi)
            
            if (debug) write (luout,*)
     &         'beta: C(E) S(ket1) C(0) complete'
            
            
c           the factor of two is for the orbital occupations,
c           assuming that ispin is never equal to two
            
c           originally the factor below was 2.0
            sum = 2d0 * ga_trace_diag(g_work)
            
            betare(idir,iresp) = betare(idir,iresp) + sum
c           betaim(idir,iresp) = 0d0

          end if                ! lgiao
            
c         ----------------------------------
c         calculate alfa from C(E) h(E) C(0)
c         ----------------------------------
          
          alo(1) = 1
          ahi(1) = nbf
          alo(2) = 1
          ahi(2) = nbf
          alo(3) = iresp ! pick direction iresp for g_dipel
          ahi(3) = iresp
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nocc
          clo(1) = 1
          chi(1) = nbf
          clo(2) = 1
          chi(2) = nocc
            

          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_dipel,alo,ahi,
     &       g_vectors(ispin),blo,bhi,
     &       g_temp,clo,chi)

          if (debug) write (luout,*)
     &       'alfa: h(E) C(0) intermediate complete'


          alo(1) = 1
          ahi(1) = nocc
          alo(2) = 1
          ahi(2) = nbf
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nocc
          clo(1) = 1
          chi(1) = nocc
          clo(2) = 1
          chi(2) = nocc

          call ga_zero(g_work)

          call nga_matmul_patch('t','n',1d0,0d0,
     &       g_vecE1,alo,ahi,
     &       g_temp,blo,bhi,
     &       g_work,clo,chi)

          if (debug) write (luout,*) 'alfa C(E) h(E) C(0) complete'

          
c         the factor of two is for the orbital occupations,
c         assuming that ispin is never equal to two

          sum = 2d0 * ga_trace_diag(g_work)

          alfare(idir,iresp) = -sum
          betaim(idir,iresp) = 0d0

c         ----------------------------------
c         destroy temporary workspace arrays
c         ----------------------------------
          
          if (.not.ga_destroy(g_vecB1))
     &       call errquit
     &       ('aor_beta: ga_destroy failed g_vecB1',
     &       0,GA_ERR)
          
          if (.not.ga_destroy(g_temp))
     &       call errquit
     &       ('aor_beta: ga_destroy failed g_temp',
     &       0,GA_ERR)
            
          enddo                 ! iresp = 1,3


c         --------------------------------------------------------
c         g_vecE1 should be deleted at this point because we might
c         enter the idir loop again
c         --------------------------------------------------------

          if (.not.ga_destroy(g_vecE1))
     &       call errquit
     &       ('aor_beta: ga_destroy failed g_vecE1',
     &       0,GA_ERR)
          
          

        end do                  ! idir = 1,3

c       -------------------------------------------
c       end loop over perturbing E-field components
c       -------------------------------------------

c       ----------------------
c       deallocate some memory
c       ----------------------
        
        do ipm = 1,ncomp
          write(cstemp,'(a,i1)') 'cphf_rhs_',ipm
          if (.not.ga_destroy(g_rhs(ipm))) call 
     &       errquit('aoresponse: ga_destroy failed '//cstemp(1:10),
     &       0,GA_ERR)
          if (lifetime) then
            write(cstemp,'(a,i1)') 'cphf_rhs_im_',ipm
            if (.not.ga_destroy(g_rhs_im(ipm))) call 
     &         errquit('aoresponse: ga_destroy failed '//cstemp(1:13),
     &         0,GA_ERR)
          endif                 ! lifetime
        enddo
        
        
        do ipm = 1,ncomp
          if (.not.ga_destroy(g_amat(ipm)))
     &       call 
     &       errquit('aoresponse: ga_destroy failed g_amat',
     &       0,GA_ERR)
          if (lifetime) then
            if (.not.ga_destroy(g_amat_im(ipm)))
     &         call 
     &         errquit('aoresponse: ga_destroy failed g_amat_im',
     &         0,GA_ERR)
          endif                 ! lifetime
        enddo                   ! ipm = 1,ncomp

            
      enddo                     ! ispin = 1,2 from way above
              
c     ---------------------------------------------------------------
c     end loop over spin components (which we don't use right now
c     since nspin is forced to be 1 at the beginning of this routine)
c     ---------------------------------------------------------------
      
c     ------------------------------------------------
c     deallocate remaining allocated memory and return
c     ------------------------------------------------
            

      if (.not.ga_destroy(g_work))
     &   call 
     &   errquit('aoresponse: ga_destroy failed g_work',
     &   0,GA_ERR)
      
      if (.not.ga_destroy(g_rhs_beta))
     &   call 
     &   errquit('aoresponse: ga_destroy failed g_rhs_beta',
     &   0,GA_ERR)
      
      if (.not.ga_destroy(g_amat_beta))
     &   call 
     &   errquit('aoresponse: ga_destroy failed g_amat_beta',
     &   0,GA_ERR)
      
      if (lifetime) then
        
        if (.not.ga_destroy(g_rhs_beta_im))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_rhs_beta_im',
     &     0,GA_ERR)
        
        if (.not.ga_destroy(g_amat_beta_im))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_amat_beta_im',
     &     0,GA_ERR)
      endif                     ! lifetime

c     it seems that if we use GIAOs everything is off by a factor of
c     two, so we need to scale betare, betaim. If we have static
c     response then there is a factor of two missing everywhere
c     because we don't add C(+) and C(-) for the electric field.

      if (lgiao) then
        scaling = half
        do idir = 1,3
          do iresp = 1,3
            betare(idir, iresp) = betare(idir, iresp) * scaling
            betaim(idir, iresp) = betaim(idir, iresp) * scaling
          end do
        end do
      end if                    ! lstatic

      if (lstatic) then
        scaling = two
        do idir = 1,3
          do iresp = 1,3
            alfare(idir, iresp) = alfare(idir, iresp) * scaling
            alfaim(idir, iresp) = alfaim(idir, iresp) * scaling
            betare(idir, iresp) = betare(idir, iresp) * scaling
            betaim(idir, iresp) = betaim(idir, iresp) * scaling
          end do
        end do
      end if                    ! lstatic


c     for comparison, let's print the
c     transpose of the negative OR tensor

      if (debug) then
      write (luout,*) 'DEBUG: OR BETA as in some other code'
      do idir = 1,3
        write (luout,'(5x,i1,1x,E13.6,1x,E13.6,1x,E13.6)')
     &    idir, -betare(1:3,idir)
      end do
      write (luout,*) 'iso: ', (betare(1,1)+betare(2,2)+betare(3,3))/3d0
      end if
     
      
      
c     ==================================================================
      
      return
      
      end

