      subroutine aor_r1_tensor(
     &  rtdb,basis,geom,        ! in : handles
     &  omega,                  ! in :
     &  lstatic,                ! in :
     &  ncomp,                  ! in :  
     &  g_smat0,                ! in : AO overlap
     &  g_dipmag,               ! in : magn -dipole mom AO
     &  g_dipel,                ! in : elect-dipole mom AO
     &  g_quadel,               ! in : quadrupole   AO
     &  g_vectors,              ! in : MOs
     &  froct,                  ! in : orbital occupations 
     &  nbf, nmo,               ! in : nr basis, nr MOs
     &  npol,                   ! in : nr. polarizations
     &  nocct, nvirt,           ! in : nocc,nvirt
     &  lgiao, lquad,           ! in : logical vars
     &  lvelocity,              ! in : logical vars
     &  lifetime,               ! in : logical vars
     &  g_vecE1,g_vecE1_im,     ! in :
     &  g_dipel_i,
     &  alfare,alfaim,          ! out: electric-electric response matrices
     &  betare,betaim)          ! out: electric-magnetic response matrices       
c $Id$     

c     =================================================================
c     
c     purpose: calculate linear response tensors
c     We assume that perturbed MO coefficients have already
c     been computed elsewhere. 
c     called from: aoresponse_driver_new
c     output: alfare, alfaim - field-electric response matrices
c             betare, betaim - field-magnetic response matrices
c     =================================================================
c
c  Written by J. Autschbach, SUNY Buffalo
c  Extension to spin-unrestricted case 
c          by F. Aquino,     Northwestern University 
c          03-15-12

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
      integer nbf
      integer npol,nocct(npol), nvirt(npol)
      double precision froct(nbf,npol)
      
c     GA handles:
      
      integer g_smat0, g_vectors(npol), g_dipel, g_quadel, g_dipmag,
     &  g_vecE1(2,2),g_vecE1_im(2,2)
      
      integer  nmo, ncomp
      double precision gamwidth, omega
      logical lgiao, lquad, lvelocity, lifetime,
     &  lstatic
      
      double precision sum_my   ! Added by FA
      
c     output:
      
      double precision alfare(3,3), alfaim(3,3)
      double precision betare(3,3), betaim(3,3)
      
c     ----------------
c     local variables:
c     ----------------
      
c     global array handles:
      
      integer g_work,g_temp
      
c     other local variables:
      
      integer nmot(2), nocvir(2)
      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)

c     dipole-quadrupole polarizability, cartesian rep.:
      double precision dipquadre(3,6), dipquadim(3,6)
c     traceless dipole-quadrupole tensor, full storage
      double precision dqpol(3,3,3)

      integer LCTensor(3,3,3)
      integer qindex(3,3)
      double precision tmpmat(3,3)

c F. Krausbeck, January 2016, ETHZ: needed for -omega*betaim   
      double precision tmpmat2(3,3)

      
      character*(256) cstemp
      character*(1) direction(3)
      data direction/'x','y','z'/
      
      integer ispin
      
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, 
     &  iresp, idir, ndir, i, j, k, l
      
      logical debug, dbgmat, dbgsub,
     &  lzora, lantisym, lroa, me
      
      double precision sum, scaling
      double precision tenm8, one, two, three, zero, half, third

      parameter (tenm8=1d-8, one=1d0, two=2d0, three=3.0d0,
     &  zero=0d0, half=one/two,
     &  third=one/three)
      
c     external functions:
      
      double precision ga_trace_diag,coeffre,coeffim
      external ga_trace_diag
      
      external get_alfaorbeta_reim
      integer g_dipel_i
c  ====================================================================

      ndir = 3                  ! nr of field directions (x,y,z)
      
      me = (ga_nodeid().eq.0)      
      
      debug  = .false. .and. me ! .true. during development
      dbgmat = .false. .and. me ! debug large matrices
      dbgsub = .false.          ! debug subroutines that are called

      if (debug) write (luout,*) 'hello from aor_r1_tensor'

c     check for ROA keyword
      lroa = .false.
      if (rtdb_get(rtdb, 'aor:roadata', mt_log, 1,
     &   lroa)) continue

      if (lroa .and. me) then
        write (luout,*) 'aor_r1_tensor: roadata output requested'
      end if

      if (lroa .and. lvelocity) 
     &   call errquit (
     &  'aor_r1_tensor: roa option not compatible with velocity.',
     &  1,INPUT_ERR)      

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

c     define translation table for quadrupole incices in
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
          
      dipquadre(:,:) = 0
      dipquadim(:,:) = 0
      
c     set parameters that control the various computational options
c     (later we will set most of this by input)
      
      lzora  = .false.      ! not yet available here 

      if (debug) write (luout,*) 'giao, velocity, magpert',
     &    lgiao, lvelocity
      lantisym = lvelocity ! antisymm. perturbation
      
c     -----------------------------------------
c     determine number of occ * virt orbitals
c     and nmot(1:2) and fix froct, if necessary
c     -----------------------------------------
      
      do ispin = 1,npol
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
      
      if (npol.eq.1) then
        coeffre=-2.0d0
        coeffim=+2.0d0   
      else if (npol.eq.2) then
        coeffre=-0.5d0
        coeffim=+0.5d0
      else
        call errquit ('aor_tens: npol value not supported',1,INPUT_ERR)
      endif
      
c     ----------------------
c     start loop over spins
c     ----------------------
      
      do ispin = 1, npol
        
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
        
c       -------------------------------------------------------
c       (A) calculate optical rotation beta from C(E) H(B) C(0)
c       ------------------------------------------------------
c       ---------------------------------------------------------
c       solution of CPKS is in g_vecE1. For the OR from a length-gauge
c       dipole response (E-field perturbation) we need the difference of
c       the +/- components so we subtract them here and store the result
c       in g_vecE1(1). If the velocity formalism is used then the
c       perturbing field has an antisymmetric AO matrix and 
c       we have to ADD the vectors instead. 
c       ---------------------------------------------------------
        
        if (ncomp.gt.1) then
          if (lantisym) then
            call ga_add(1d0, g_vecE1(ispin,1), 
     &        1d0, g_vecE1(ispin,2),
     &        g_vecE1(ispin,1))
            if (lifetime) then
              call ga_add(1d0, g_vecE1_im(ispin,1), 
     &          1d0, g_vecE1_im(ispin,2),
     &          g_vecE1_im(ispin,1))
            end if              ! lifetime
          else                  ! lantisym ? 
            call ga_add(1d0, g_vecE1(ispin,1), 
     &        -1d0, g_vecE1(ispin,2),
     &        g_vecE1(ispin,1))
            if (lifetime) then
              call ga_add(1d0, g_vecE1_im(ispin,1),
     &          -1d0, g_vecE1_im(ispin,2),
     &          g_vecE1_im(ispin,1))
            end if              ! lifetime
          end if                ! lantisym
        endif                   ! ncomp.gt.1        
        
        do idir = 1,ndir        ! direction of the perturbing field
          
          do iresp = 1,3
            
            call get_alfaorbeta_reim(
     &        betare(idir,iresp), ! in/out: alpha or beta real part
     &        betaim(idir,iresp), ! in/out: alpha or beta im   part
     &        g_vecE1(ispin,1), ! in : 1st-order pert vec RE
     &        g_vecE1_im(ispin,1), ! in : 1st-order pert vec IM
     &        g_dipmag,         ! in : dipole electric or magnetic
     &        g_vectors(ispin), ! in : MO vectors
     &        idir,             ! in : = 1,2,3=x,y,z directions
     &        iresp,            ! in : = 1,2,3
     &        coeffre,-coeffim,1, ! in : (coeffre,coeffim,caseAO)
     &        nbf,              ! in : nr. basis functions
     &        nocc,             ! in : nr. occupied alpha or beta
     &        lifetime,         ! in : logical var for damping
     &        dbgsub,            ! in : logical var for debugging
     &        g_dipel_i,
     &        omega,
     &        g_temp)           ! in : scratch GA array
            
            
          enddo                 ! iresp (responding field components)
        enddo                   ! idir    (perturbing E-field components)
        
c       --------------------------------------
c       (B) calculate alfa from C(E) h(E) C(0)
c       --------------------------------------
c       --------------------------------------------------------- 
c       For alfa we need the sum of the +/- components no matter
c       if we use
c       length or velocity gauge so we add twice the icomp=2 component
c       back into g_vecE1(1)
c       ---------------------------------------------------------
        
        if (ncomp.gt.1) then
          if (lantisym) then
            continue
          else
            call ga_add(1d0, g_vecE1(ispin,1),
     &        2d0, g_vecE1(ispin,2),
     &        g_vecE1(ispin,1))
            if (lifetime) then
              call ga_add(1d0, g_vecE1_im(ispin,1), 
     &          2d0, g_vecE1_im(ispin,2),
     &          g_vecE1_im(ispin,1))
              
            end if              ! lifetime
          end if                ! lantisym
        endif                   ! ncomp.gt.1
        
        do idir = 1,ndir        ! direction of the perturbing field
          
          do iresp = 1,3
            
            call get_alfaorbeta_reim(
     &        alfare(idir,iresp), ! out: alpha or beta real part
     &        alfaim(idir,iresp), ! out: alpha or beta im   part
     &        g_vecE1(ispin,1), ! in : 1st-order pert vec RE
     &        g_vecE1_im(ispin,1), ! in : 1st-order pert vec IM
     &        g_dipel,          ! in : dipole electric or magnetic
     &        g_vectors(ispin), ! in : MO vectors
     &        idir,             ! in : = 1,2,3=x,y,z directions
     &        iresp,            ! in : = 1,2,3
     &        coeffre,-coeffim,1, ! in : (coeffre,coeffim,caseAO) 
     &        nbf,              ! in : nr. basis functions
     &        nocc,             ! in : nr. occupied alpha or beta
     &        lifetime,         ! in : logical var for damping
     &        dbgsub,            ! in : logical var for debugging
     &        g_dipel_i,
     &        omega,
     &        g_temp)           ! in : scratch GA array
            
            
          enddo                 ! iresp (responding field components)
          
c         if required, calculate dipole-quadrupole polarizability:
          
          if (lquad) then
            
            do iresp = 1,6
              call get_alfaorbeta_reim(
     &          dipquadre(idir,iresp), ! out: quad real part
     &          dipquadim(idir,iresp), ! out: quad im   part
     &          g_vecE1(ispin,1), ! in : 1st-order pert vec RE
     &          g_vecE1_im(ispin,1), ! in : 1st-order pert vec IM
     &          g_quadel,       ! in : dipole-quadrupole
     &          g_vectors(ispin), ! in : MO vectors
     &          idir,           ! in : = 1,2,3=x,y,z directions
     &          iresp,          ! in : = 1,2,3
     &          coeffre,-coeffim,1, ! in : (coeffre,coeffim,caseAO) 
     &          nbf,            ! in : nr. basis functions
     &          nocc,           ! in : nr. occupied alpha or beta
     &          lifetime,       ! in : logical var for damping
     &          dbgsub,          ! in : logical var for debugging
     &          g_dipel_i,
     &          omega,
     &          g_temp)         ! in : scratch GA array       
            end do              ! iresp
          end if                ! lquad
          
        enddo                   ! idir    (perturbing E-field components)
               
        if (.not.ga_destroy(g_temp))
     &    call errquit
     &    ('aor_beta: ga_destroy failed g_temp',
     &    0,GA_ERR)
        
      enddo                     ! ispin

c     -------------------------
c     done with loop over spins
c     -------------------------      

c     deallocate work array
      
      if (.not.ga_destroy(g_work)) call errquit(
     &  'aoresponse: ga_destroy failed g_work',1,GA_ERR)

c     with static response we need to multiply by 2:      
      
      if (lstatic) then
        scaling = two
        do idir = 1,ndir
          do iresp = 1,3
            alfare(idir, iresp) = alfare(idir, iresp) * scaling
            alfaim(idir, iresp) = alfaim(idir, iresp) * scaling
            betare(idir, iresp) = betare(idir, iresp) * scaling
            betaim(idir, iresp) = betaim(idir, iresp) * scaling
          end do
        end do
c       ... getting tired of that ol' f77 syntax ... 
        if (lquad) then
          dipquadre(:,:) = dipquadre(:,:) * scaling
          dipquadim(:,:) = dipquadim(:,:) * scaling
        end if                
      end if                    ! lstatic
      
c     for electric field perturbations (default), alfa is the
c     polarizability and beta is -G'. Let's fix that sign and
c     also consider the case of velocity formalisms where 
c     additional divisions by omega occur.

      betare(:,:) = -betare(:,:)
      betaim(:,:) = -betaim(:,:)

      if (.not.lstatic .and. lvelocity) then  
        do idir = 1,ndir 
          do iresp = 1,3          
            scaling = one/omega
            betare(idir,iresp) = betare(idir,iresp) *scaling
            betaim(idir,iresp) = betaim(idir,iresp) *scaling
            scaling = one/(omega**2)
            alfare(idir,iresp) = alfare(idir,iresp) * scaling
            alfaim(idir,iresp) = alfaim(idir,iresp) * scaling
          enddo
        enddo
      end if

c     ------------------------------------------------------------
c     if requested by input, construct traceless dipole-quadrupole
c     polarizability tensor and contruct Buckingham-Dunn OR tensor
c     At present we do only the real part !
c     ------------------------------------------------------------
      
      if (lquad) then

        if (me) then
          write (luout,*) 'dip-quad polarizability Re:'
          call output (dipquadre,1,3,1,6,3,6,1)
        end if
        
        dqpol(:,:,:) = zero
        
        do i = 1,3
          do j = 1,3
            do k = 1,j
              dqpol(i,j,k) = three * dipquadre(i,qindex(j,k))
              if (k.ne.j)
     &          dqpol(i,k,j) = three * dipquadre(i,qindex(k,j))
            end do
          end do
        end do
        
        do i = 1,3
c         calculate the trace of the quadrupole part
          sum = zero
          do j = 1,3
            sum = sum + dipquadre(i,qindex(j,j))
          end do
c         subtract trace from diagonal
          do j = 1,3
            dqpol(i,j,j) =  dqpol(i,j,j) - sum
          end do
        end do                  ! i
        
        dqpol(:,:,:) = dqpol(:,:,:) * half

        if (me) then
          write (luout,*) 'traceless dip-quad polarizability:'
          do i = 1,3
            call output (dqpol(i,1:3,1:3),1,3,1,3,3,3,1)
          end do
        end if

c       if we have the keyword for the ROA output, let's do that
c       before the BD tensor gets assembled:


c F. Krausbeck, January 2016: ETH Zurich: the same loop as above
c for tmpmat2 -> -omega*betaim(i,j) required for resonance VROA
c ... jochen 08/17: had to make some changes to Florian's code because
c     of the way this routine was re-worked in the meantime. In the
c     old version, 'betare' etc. was beta, but now it is G'.
        if (lroa .and. me .and. .not.lstatic) then
          do i = 1,3
            do j = 1,3
              tmpmat(i,j) = -betare(i,j)/omega
              tmpmat2(i,j) = -betaim(i,j)/omega
            end do
          end do


c F. Krausbeck, January 2016, ETH Zurich: modified for Re+Im output
c and commented out the old function call 
           call aor_write_roadata2
     &       (luout, alfare, alfaim, tmpmat, betare,  
     &        tmpmat2, betaim, 
     &        dipquadre, dipquadim, dqpol)
c          call aor_write_roadata
c     &       (luout, alfare, betare, tmpmat, dipquadre, dqpol)
        end if                  ! lroa

c       assemble BD tensor in tmpmat.

c       step 1: symmetrize beta tensor
        
        tmpmat(:,:) = half * betare(:,:)

        do i = 1,3
          do j = 1,3
            tmpmat(i,j) = tmpmat(i,j) + half*betare(j,i)
          end do
        end do

c       step 2: add quadrupole contributions

        scaling = -half * third * omega

        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                tmpmat(i,j) = tmpmat(i,j) 
     &             + scaling * LCtensor(i,k,l) * dqpol(k,l,j) 
     &             + scaling * LCTensor(j,k,l) * dqpol(k,l,i)
              end do
            end do
          end do
        end do

        betare(:,:) = tmpmat(:,:)
        
      end if                    ! lquad      
      
c     ----------------
c     all done. return
c     ----------------
      return
      end
