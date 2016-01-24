      subroutine aor_r1_beta_anl(
     &      rtdb,basis,geom,    ! in : handles
     &      omega,              ! in :
     &      lstatic,            ! in :
     &      ncomp,              ! in :  
     &      g_smat0,            ! in :
     &      g_sket1,            ! in :
     &      g_vecB1,            ! in : magn -dipole mom AO
     &      g_dipel,            ! in : elect-dipole mom AO
     &      g_quadel,           ! in : quadrupole   AO
     &      g_vectors,          ! in : MOs
     &      froct,              ! in : set of occupations
     &      epst,               ! in : set of energy values
     &      nbf, nmo,           ! in : nr. basis, nr. MOs
     &      npol,               ! in : nr. polarizations
     &      nocct, nvirt,       ! in : nocc,nvirt
     &      lgiao, lquad,       ! in : logical vars
     &      lanalyze,           ! in : logical vars
     &      lvelocity,          ! in : logical vars
     &      lifetime,           ! in : logical vars
     &      g_vecE1,g_vecE1_im) ! out: 
c $Id$     

c     =================================================================
c     purpose: analyze beta tensor. See routine
c     aor_r1_beta.F for additional comments. 
c     Use a molecular orientation in which the chiral
c     response tensor is diagonal!
c     called from: aoresponse_driver_new
c     =================================================================

c     Written by J. Autschbach, SUNY Buffalo
c     Extension to spin-unrestricted case 
c     by F. Aquino, Northwestern University 03-15-12

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
      integer npol  
c     These are all input, too
      integer g_smat0,g_dipel, 
     &        g_sket1,g_quadel,
     &        g_vectors(npol),
     &        g_vecB1(npol),
     &        g_vecE1(2,2),g_vecE1_im(2,2)
      integer nfreq, response_order, nbf, nmo, ncomp
      integer nocct(npol),nvirt(npol)
      double precision froct(nbf,2),epst(nbf,2)
      double precision gamwidth, omega
      logical lgiao, lquad, lanalyze, lvelocity, lifetime,
     &   lstatic
c     ----------------
c     local variables:
c     ----------------
c     global array handles:     
      integer g_anl, g_work,
     &        g_temp, g_tmpanl,g_tran,g_vectmp(2)
      integer l_diag, k_diag
c     other local variables: 
      integer nmot(npol),nocvir(npol),nopen(npol)
c      data nopen(1),nopen(2)/0,0/
      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
c     dipole-quadrupole polarizability, cartesian rep.:
      double precision dipquadre(3,6), dipquadim(3,6)
c     traceless dipole-quadrupole tensor, full storage
      double precision dqpol(3,3,3)
      integer qindex(3,3)
      double precision tmpmat(3,3)
      character*(256) cstemp
      character*(1) direction(3)
      data direction/'x','y','z'/     
      integer ispin, nspin
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, idir, iresp,
     &        i,j,k,l
      logical debug, dbgmat, 
     &        lzora, lantisym, lmo, status, oprint
      double precision sum, scaling
      double precision tenm8, one, two, three, zero, half, third, four
      parameter (tenm8=1d-8, one=1d0, two=2d0, three=3d0,
     &           zero=0d0, half=one/two,
     &           third=one/three, four=4.0d0)
c     external functions:
      character*(256) lmotrans
      logical file_read_ga
      external file_read_ga
      double precision trace,ga_trace_diag,sum_lmo
      external ga_trace_diag,get_C1MC,
     &         get_trace_ganl,
     &         trace,
     &         get_tracelessQtensor,write_vects,
     &         update_g_anl_quad,
     &         update_g_anl_elmag,
     &         debug_update_g_anl_quad1
c  ====================================================================

      debug  = .false. ! .true. during development
      dbgmat = .false. ! debug large matrices
      oprint = ga_nodeid().eq.0

      if (debug) write (luout,*) 'hello from aor_r1_beta_anl'

      dipquadre(:,:) = 0
      dipquadim(:,:) = 0

c     make sure lvelocity.ne.T., we do not support that in this
c     subroutine to keep the clutter at a manageable level.

      if (lvelocity) 
     &   call errquit ('aor_beta: lvelocity set',1,INPUT_ERR)
      
c     define translation table for quadrupole indices in
c     packed storage
c     XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6
      
      qindex(1,1) = 1
      qindex(1,2) = 2
      qindex(2,1) = 2
      qindex(1,3) = 3
      qindex(3,1) = 3
      qindex(2,2) = 4
      qindex(2,3) = 5
      qindex(3,2) = 5
      qindex(3,3) = 6
      
c     set parameters that control the various computational options
c     (later we will set most of this by input)
      
      lzora      = .false.      ! not yet available here 

      if (debug) write (luout,*) 'giao, velocity',
     &  lgiao, lvelocity
      
c     -----------------------------------------
c     determine number of occ * virt orbitals
c     and nmot(1:2) and fix froct, if necessary
c     -----------------------------------------
      
      do ispin=1,npol
        nocvir(ispin) = nocct(ispin) * nvirt(ispin)
        nmot(ispin) = nmo
        if (nmo .lt.nbf) then
          do imo = nmo+1,nbf
            froct(imo,ispin) = 0d0
          enddo
        endif
      enddo                     ! ispin
      
c     ------------------------------
c     allocate some temp. work space
c     ------------------------------
      
      chunk(1) = nbf
      chunk(2) = -1
      dims(1)  = nbf
      dims(2)  = nbf
      
      write(cstemp,'(a)') 'work'
      if (.not.nga_create(MT_DBL,2,dims,cstemp(1:4),chunk,
     &  g_work)) call 
     &  errquit('aoresponse: nga_create failed: '//cstemp(1:4),
     &  1,GA_ERR)     
      call ga_zero (g_work)
      
      chunk(1) = nbf
      chunk(2) = -1
      dims(1)  = nbf
      dims(2)  = nbf
      
      write(cstemp,'(a)') 'vectmp(1)'
      if (.not.nga_create(MT_DBL,2,dims,cstemp(1:4),chunk,
     &  g_vectmp(1))) call 
     &  errquit('aoresponse: nga_create failed: '//cstemp(1:9),
     &  0,GA_ERR)
      
c     the factor of two is for the orbital occupations,
c     assuming that ispin is never equal to two
      
      if (npol.eq.1) then        
        scaling = two
        if (lstatic) scaling = four
      else  if (npol.eq.2) then     
        scaling = 1.0d0
        if (lstatic) scaling = 2.0d0
      endif

c     ---------------------
c     start loop over spins
c     ---------------------
      
      do ispin=1,npol
        
        nmo1 = nmot(ispin)      ! total no.of MOs for this spin
        nocc = nocct(ispin)     ! occupied MOs
        nvir = nvirt(ispin)     ! virtual MOs
        nocv = nocvir(ispin)    ! nocc * nvir
        
c       allocate intermediate vector for matrix multiplications
c       used to create the final results
        
        write (cstemp,'(a)') 'aor_beta: temp1'
        if(.not.ga_create(MT_DBL, nbf, nocc, trim(cstemp),
     &    -1,-1,g_temp))
     &    call errquit (trim(cstemp),0,GA_ERR)
        if (debug) write (luout,*) 'g_temp allocated'
        
c       ----------------------------------------------
c       check if we have localized MOs on file. If yes
c       read them, assuming nspin.eq.1
c       ----------------------------------------------
        
        i = 1
        status=rtdb_get(rtdb,'prop:pmlocalization',MT_INT,1,i)
        lmo = (i.eq.0)
        if (lmo) then
          if (oprint) write (luout,*) 'analysis: LMO switch found'
          write (cstemp,'(a)') 'aor_beta: g_tran'
          if(.not.ga_create(MT_DBL,nocc,nocc,trim(cstemp),-1,-1,
     &      g_tran))            ! create g_tran
     &      call errquit (trim(cstemp),0,GA_ERR)
          if (debug) write (luout,*) 'g_tran allocated'
          call util_file_name('lmotrans',.true.,.true.,lmotrans)
          if(.not.file_read_ga(lmotrans,g_tran)) call errquit
     $      ('aor_r1_beta_anl: could not read lmotrans',0, DISK_ERR)
        end if                  ! lmo
        
c       allocate matrix that accumulates the analysis data
        
        write (cstemp,'(a)') 'aor_beta: g_anl'
        if(.not.ga_create(MT_DBL,nocc,nocc,trim(cstemp),-1,-1,
     &    g_anl))               ! create g_anl
     &    call errquit (trim(cstemp),0,GA_ERR)
        if (debug) write (luout,*) 'g_anl allocated'
        call ga_zero(g_anl)  
        
c       diagonal elements of the last matrix
        
        if (.not. ma_push_get(mt_dbl,nocc,'diag', l_diag, k_diag))
     &    call errquit('error alloc MA diag', 0, MA_ERR)
        
c       lmos: debug transformation
        
        if (lmo .and. debug) then
          call ga_print(g_tran)
          call ga_dgemm('t', 'n', nocc, nocc, nocc, 
     $      1.0d0, g_tran, g_tran, 0.0d0, g_anl)
          call ga_print(g_anl)
        end if
        
c       ---------------------------------------------------------
c       solution of CPKS is in g_vecE1. Below we need
c       only the sum of the +/- components so we add them here
c       and store them in g_vecE1(1)
c       ---------------------------------------------------------
        
        if (ncomp.gt.1) then
          call ga_add(1d0, g_vecE1(ispin,1), 
     &      1d0, g_vecE1(ispin,2),
     &      g_vecE1(ispin,1))
          if (lifetime) then
            call ga_add(1d0, g_vecE1_im(ispin,1), 
     &        1d0, g_vecE1_im(ispin,2),
     &        g_vecE1_im(ispin,1))
          end if
        endif
        
c       ------------------------------------------------
c       for Buckingham-Dunn tensor we need the traceless
c       quadrupole tensor
c       ------------------------------------------------
        
        if (lquad) then 
          call get_tracelessQtensor(g_quadel, nbf, g_work)          
        end if                  ! lquad    
        
c       ---------------------------------------------------------
c       start loop over the components of the response tensor and 
c       calculate the final results
c       ---------------------------------------------------------
        
        do idir = 1,3           ! direction of the perturbing field
          
c         g_anl is going to accumulate the results
          call ga_zero(g_anl)
          
          if (oprint)
     &      write (luout,'(1x,40(''-'')/1x,a,2i1)')
     &      'MO analysis of OR tensor component ',idir,idir                                 
c         -------------------------------------------------------
c         (A) calculate optical rotation beta from C(E) S(0) C(B)
c         -------------------------------------------------------
          
          
          call update_g_anl_elmag(
     &      g_anl,              ! in/out:
     &      g_vecE1(ispin,1),   ! in    : C(E) pert MO vect.
     &      g_smat0,            ! in    : 0th/1st overlap deriv.
     &      g_vecB1(ispin),     ! in    : MO vect or C(B) MO vect
     &      scaling,            ! in    : scaling factor
     &      idir,               ! in    : direction (x,y or z)
     &      3,                  ! in    : 3 or 4
     &      nbf,                ! in    : nr. basis functions
     &      nocc,               ! in    : nr. occupied MOs
     &      lstatic,            ! in    : flag for static calc.
     &      debug,              ! in    : flag for debugging
     &      g_work,             ! in    : scratch GA array
     &      g_temp)             ! in    : scratch GA array
          
          
          if (lgiao) then
            
c           --------------------------------------
c           if we use GIAOs there is a second term
c           in beta which is C(E) S(1ket) C(0)
c           --------------------------------------  
            
            
            call update_g_anl_elmag(
     &        g_anl,            ! in/out:
     &        g_vecE1(ispin,1), ! in    : C(E) pert MO vect.
     &        g_sket1,          ! in    : 0th/1st overlap deriv.
     &        g_vectors(ispin), ! in    : MO vect or C(B) MO vect
     &        scaling,          ! in    : scaling factor
     &        idir,             ! in    : direction (x,y or z)
     &        4,                ! in    : 3 or 4
     &        nbf,              ! in    : nr. basis functions
     &        nocc,             ! in    : nr. occupied MOs
     &        lstatic,          ! in    : flag for static calc.
     &        debug,            ! in    : flag for debugging
     &        g_work,           ! in    : scratch GA array
     &        g_temp)           ! in    : scratch GA array

          end if ! lgiao
          
          
c         ----------------------------------------------------
c         if requested by input, add to OR beta the quadrupole
c         polarizability terms
c         ----------------------------------------------------                  
          if (lquad) then   
            
            call update_g_anl_quad(
     &        g_anl,            ! in/out:
     &        g_quadel,         ! in    : quadrupolo AO integral
     &        g_vecE1(ispin,1), ! in    : 1st-order elect. pert. MO vector
     &        g_vectors(ispin), ! in    : MO vector
     &        scaling,          ! in    : scaling factor
     &        g_temp,           ! in    : scratch GA array
     &        g_work,           ! in    : scratch GA array
     &        idir,             ! in    : =1,2,3=x,y,z
     &        nocc,             ! in    : nr. occupied MOs      
     &        nbf,              ! in    : nr. basis functions    
     &        debug,            ! in    : debugging flag
     &        lstatic)          ! in    : static    flag
            
          endif                 ! lquad
          
c         ---------------------
c         Canonical MO analysis
c         ---------------------
          
          if (oprint) write (luout,
     &      '(/t12,a,t26,a/t11,6(''-''),t22,12(''-''))')
     &      'CMO #','contrib.'
          
          call ga_get_diagonal(g_anl, dbl_mb(k_diag) )
          
          sum = zero
          do i = 1,nocc
            sum = sum + dbl_mb(k_diag+i-1)
            if (oprint) write (luout,'(t11,i6,t22,f12.4)')
     &        i,dbl_mb(k_diag+i-1)
          end do     
          
          if (oprint)
     &      write (luout,'(1x,a,2i1,a,f12.4)') 'Component ',idir,idir,
     &      ': Sum = ',sum
          
          if (debug) then
            sum = ga_trace_diag(g_anl)
            if (oprint) write (luout,*) 'sum from ga_trace: ',sum
          end if                ! debug
          
          if (oprint) write (luout,'(1x,40(''-''))')
          
c         ---------------------
c         Localized MO analysis
c         ---------------------
          
          if (lmo) then
            
            call get_trace_ganl(g_anl, ! in/out:
     &        g_tran,           ! in:
     &        idir,             ! in: =1,2,3=x,y,z
     &        oprint,           ! in: logical var
     &        nocc)             ! in: nr. occupied MOs    
            
          end if                ! lmo
          
        end do                  ! idir = 1,3
        
c       -------------------------------------------
c       end loop over perturbing E-field components
c       -------------------------------------------
        
c       -----------------
c       deallocate memory
c       -----------------
        
        if (.not. ma_pop_stack(l_diag))
     &    call errquit('error deloc MA diag',0, MA_ERR)
        
        
        call write_vects(
     &    rtdb,geom,basis,      ! in: handles
     &    g_vecE1(ispin,1),     ! in: C^(1,E)
     &    g_vecB1(ispin),       ! in: C^(1,B)
     &    g_tran,               ! in:
     &    g_vectmp,             ! in: scratch GA array
     &    npol,nocc,            ! in: nr. polariz,occ MOs
     &    nocct,nvirt,          ! in: nr. occ,virt MOs
     &    nopen,nmot,           ! in: nr. open shells,nmot=nocc*nvirt
     &    froct,epst,           ! in:
     &    nbf,                  ! in: nr. basis functions
     &    lmo,                  ! in: logical flag
     &    debug)                ! in: logical flag for debugging
        
        
c       -------------------------------------------------
c       un-add the frequency components in vec_E1 in case
c       we reuse these arrays:
c       -------------------------------------------------
        
        if (ncomp.gt.1) then
          call ga_add(1d0, g_vecE1(ispin,1), 
     &      -1d0, g_vecE1(ispin,2),
     &      g_vecE1(ispin,1))
          if (lifetime) then
            call ga_add(1d0, g_vecE1_im(ispin,1), 
     &        -1d0, g_vecE1_im(ispin,2),
     &        g_vecE1_im(ispin,1))
          end if
        endif                   ! ncomp>1
        
        if (.not.ga_destroy(g_anl))  call errquit
     &    ('aor_beta: ga_destroy failed g_anl',0,GA_ERR)  
        if (lmo) then
          if (.not.ga_destroy(g_tran)) call errquit
     &      ('aor_beta: ga_destroy failed g_tran',0,GA_ERR)
        endif     
        if (.not.ga_destroy(g_temp)) call errquit
     &    ('aor_beta: ga_destroy failed g_temp',0,GA_ERR)
        
      enddo                     ! ispin = 1,2 from way above
      
      if (.not.ga_destroy(g_vectmp(1))) call errquit 
     &  ('aor_beta: ga_destroy failed g_vectmp(1)',0,GA_ERR)   
      if (.not.ga_destroy(g_work)) call errquit
     &  ('aoresponse: ga_destroy failed g_work',0,GA_ERR)
      
      return     
      end
