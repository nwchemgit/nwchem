      subroutine aor_r1_beta(
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
     &      nbf, nmo,           ! in : nr basis, nr MOs
     &      npol,               ! in : nr. polarizations
     &      nocct, nvirt,       ! in : nocc,nvirt
     &      lgiao, lquad,       ! in : logical vars
     &      lanalyze,           ! in : logical vars
     &      lvelocity,          ! in : logical vars
     &      lifetime,           ! in : logical vars
     &      lmagpert,           ! in : logical vars
     &      g_vecE1,g_vecE1_im, ! in : 
     &      alfare,alfaim,      ! out: electric-electric response matrices
     &      betare,betaim)      ! out: electric-magnetic response matrices      
c $Id$
      
c     =================================================================
      
c     purpose: calculate linear response tensors

c     CALCULATE BETA DIRECTLY using 
c     frequency dependent version of Amos equation as published
c     by J. Autschbach, Comp. Lett. 3 (2007), 131
c     This allows to compute the static limit directly instead of
c     using -G'/omega for small omega which is numerically unstable.
c     Moreover, we can use a static magnetic field which makes it
c     a little bit easier to implement GIAO optical rotation. 

c     We assume that perturbed MO coefficients have already
c     been computed elsewhere, and that g_vecE1 is for an
c     ELECTRIC field using the dipole-length perturbation operator

c     called from: aoresponse_driver_new

c     output: alfare, alfaim - electric-electric response matrices
c             betare, betaim - electric-magnetic response matrices

c     =================================================================
c
c  Written by J. Autschbach, SUNY Buffalo
c  Extension to spin-unrestricted case 
c          by F. Aquino,     Northwestern University 
c          03-15-12
c --> Experimental (not published yet)

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
      integer g_smat0,g_dipel, 
     &        g_sket1,g_quadel,
     &        g_vectors(npol),
     &        g_vecB1(npol),
     &        g_vecE1(2,2),g_vecE1_im(2,2)
      integer nfreq, response_order, nbf, nmo, ncomp
      integer nocct(npol),nvirt(npol)
      double precision froct(nbf,2)
      double precision gamwidth, omega
      logical lgiao,lquad,lanalyze,
     &        lvelocity, lifetime, lmagpert,
     &        lstatic
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
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3),
     &        clo(3),chi(3)
c     dipole-quadrupole polarizability, cartesian rep.:
      double precision dipquadre(3,6), dipquadim(3,6)
c     traceless dipole-quadrupole tensor, full storage
      double precision dqpol(3,3,3)
      integer LCTensor(3,3,3)
      integer qindex(3,3)
      double precision tmpmat(3,3)
      character*(256) cstemp
      character*(1) direction(3)
      data direction/'x','y','z'/
      integer ispin,nspin,ndir
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, iresp, idir,
     &        i,j,k,l
      logical debug, dbgmat, 
     &        lzora, lantisym, lroa, me
      double precision sum, scaling,coeffre,coeffim
      double precision tenm8, one, two, three, zero, half, third
      parameter (tenm8=1d-8, one=1d0, two=2d0, three=3d0,
     &           zero=0d0, half=one/two,
     &           third=one/three)
c     external functions:
      double precision ga_trace_diag
      external ga_trace_diag
c  ====================================================================

      me = (ga_nodeid().eq.0)

      debug  = .false. .and. me ! .true. during development
      dbgmat = .false. .and. me ! debug large matrices

c      debug=.true.

c     check for ROA keyword
      lroa = .false.
      if (rtdb_get(rtdb, 'aor:roadata', mt_log, 1,
     &   lroa)) continue

      if (lroa .and. me) then
        write (luout,*) 'aor_r1_beta: roadata output requested'
      end if

      if (debug) write (luout,*) 'hello from aor_r1_beta'
c     the main results are collected in alfare/in, betare/im.
c     initialize these and other tensors with zeros:
      ndir=3 ! nr. directions (x,y,z)
      do idir=1,ndir
       do iresp=1,3
        alfare(idir,iresp)    = 0.0d0
        alfaim(idir,iresp)    = 0.0d0
        betare(idir,iresp)    = 0.0d0
        betaim(idir,iresp)    = 0.0d0
       end do
       do iresp = 1,6
        dipquadre(idir,iresp) = 0.0d0
        dipquadim(idir,iresp) = 0.0d0
       enddo ! end-loop-iresp
      enddo ! end-loop-idir
c     make sure lvelocity.ne.T., we do not support that in this
c     subroutine to keep the clutter at a manageable level.
c     same for lmagpert
c      if (lvelocity) 
c     &   call errquit ('aor_beta: lvelocity set',1,INPUT_ERR)
      if (lmagpert) 
     &   call errquit ('aor_beta: lmagpert set',1,INPUT_ERR)
c     -------------------------
c     define Levi-Civita tensor for quadrupole additions
c     -------------------------
      LCtensor(:,:,:) =  0      
      LCtensor(1,2,3) =  1
      LCtensor(2,3,1) =  1
      LCtensor(3,1,2) =  1      
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
c     set parameters that control the various computational options
c     (later we will set most of this by input)
      lzora      = .false.      ! not yet available here 
c     check if we should calculate OR beta directly or G'
      if (ga_nodeid().eq.0) write (luout,*)
     & 'Optical Rotation: will calculate Beta directly, not G'''
      if (debug) write (luout,*) 'giao, velocity',
     &    lgiao, lvelocity
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
      enddo ! end-loop-ispin
c       ------------------------------
c       allocate some temp. work space
c       ------------------------------
        chunk(1) = nbf
        chunk(2) = -1
        dims(1)  = nbf
        dims(2)  = nbf        
        write(cstemp,'(a)') 'work'
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:4),chunk,
     &     g_work)) call 
     &     errquit('aoresponse: nga_create failed: '//cstemp(1:4),
     &     0,GA_ERR)     
        call ga_zero (g_work)
      if      (npol.eq.1) then
       coeffre =-2.0d0
       coeffim =+2.0d0   
      else if (npol.eq.2) then
       coeffre=-0.5d0
       coeffim=+0.5d0  
      endif
      do ispin=1,npol
        nmo1 = nmot(ispin)      ! total no.of MOs for this spin
        nocc = nocct(ispin)     ! occupied MOs
        nvir = nvirt(ispin)     ! virtual MOs
        nocv = nocvir(ispin)    ! nocc * nvir
c       allocate intermediate vector for matrix multiplications
c       used to create the final results
        write (cstemp,'(a)') 'aor_beta: temp1'
        if(.not.ga_create(MT_DBL, nbf, nocc, trim(cstemp),
     &     -1,-1,g_temp))
     &     call errquit (trim(cstemp),0,GA_ERR)
        if (debug) write (luout,*) 'g_temp allocated'    
c       ---------------------------------------------------------
c       solution of CPKS is in g_vecE1. Below we need
c       only the sum of the +/- components so we add them here
c       and store them in g_amat(1)
c       ---------------------------------------------------------
        if (ncomp.gt.1) then
          call ga_add(1d0, g_vecE1(ispin,1), 
     &                1d0, g_vecE1(ispin,2),
     &                     g_vecE1(ispin,1))
          if (lifetime) then
          call ga_add(1d0, g_vecE1_im(ispin,1), 
     &                1d0, g_vecE1_im(ispin,2),
     &                     g_vecE1_im(ispin,1))
          end if
        endif       
c       -------------------------------------------------------
c       now loop over the components of the response tensor and 
c       calculate the final results
c       -------------------------------------------------------      
        do idir = 1,3           ! direction of the perturbing field     
          do iresp = 1,3        ! responding operator direction          
c           -------------------------------------------------------
c           (A) calculate optical rotation beta from C(E) S(0) C(B)
c           ------------------------------------------------------    
c           if (ga_nodeid().eq.0) then
c            write(*,1) idir,iresp
c 1          format('==== compute beta(',i3,',',i3,') ====== START')
c           endif

              call get_alfaorbeta_reim(
     &            betare(idir,iresp), ! out: alpha or beta real part
     &            betaim(idir,iresp), ! out: alpha or beta im   part
     &            g_vecE1(ispin,1),   ! in : 1st-order pert vec RE
     &            g_vecE1_im(ispin,1),! in : 1st-order pert vec IM
     &            g_smat0,            ! in : dipole electric or magnetic
     &            g_vecB1(ispin),     ! in : MO vectors
     &            idir,               ! in : = 1,2,3=x,y,z directions
     &            iresp,              ! in : = 1,2,3
     &            -coeffre,-coeffim,2,! in : (coeffre,coeffim,caseAO)
     &            nbf,                ! in : nr. basis functions
     &            nocct(ispin),       ! in : nr. occupied alpha or beta
     &            lifetime,           ! in : logical var for damping
     &            debug,              ! in : logical var for debugging
     &            g_temp)             ! in : scratch GA array

c           if (ga_nodeid().eq.0) then
c            write(*,3) idir,iresp
c 3          format('==== compute beta(',i3,',',i3,') ====== END')
c           endif
           if (debug) then
            if (ga_nodeid().eq.0) then
             write(*,2) ispin,idir,iresp,
     &                 betare(idir,iresp),betaim(idir,iresp)
 2           format('(beta(',i3,',',i3,',',i3,
     &             ')=(',f15.8,',',f15.8,')')
            endif
           endif ! end-if-debug

           if (debug) write (luout,*) 'beta: C(E) S(0) C(B) complete'
           if (lgiao) then          
c             --------------------------------------
c             if we use GIAOs there is a second term
c             in beta which is C(E) S(1ket) C(0)
c             --------------------------------------  
c            if (ga_nodeid().eq.0) then
c             write(*,4) idir,iresp
c 4           format('==== compute beta-giao(',
c     &              i3,',',i3,') ====== START')
c            endif

            call get_alfaorbeta_reim(
     &            betare(idir,iresp), ! out: alpha or beta real part
     &            betaim(idir,iresp), ! out: alpha or beta im   part
     &            g_vecE1(ispin,1),   ! in : 1st-order pert vec RE
     &            g_vecE1_im(ispin,1),! in : 1st-order pert vec IM
     &            g_sket1,            ! in : dipole electric or magnetic
     &            g_vectors(ispin),   ! in : MO vectors
     &            idir,               ! in : = 1,2,3=x,y,z directions
     &            iresp,              ! in : = 1,2,3
     &            -coeffre,-coeffim,1,! in : (coeffre,coeffim,caseAO)
     &            nbf,                ! in : nr. basis functions
     &            nocct(ispin),       ! in : nr. occupied alpha or beta
     &            lifetime,           ! in : logical var for damping
     &            debug,              ! in : logical var for debugging
     &            g_temp)             ! in : scratch GA array

c             if (ga_nodeid().eq.0) then
c              write(*,5) idir,iresp
c 5            format('==== compute beta-giao(',
c     &               i3,',',i3,') ====== END')
c             endif
             if (debug) then        
              if (ga_nodeid().eq.0) then
               write(*,6) ispin,idir,iresp,
     &                   betare(idir,iresp),betaim(idir,iresp)
 6             format('(beta-giao(',i3,',',i3,',',i3,
     &               ')=(',f15.8,',',f15.8,')')
              endif
             endif ! end-if-debug
           endif ! end-if-lgiao      
c           --------------------------------------
c           (B) calculate alfa from C(E) h(E) C(0)
c           --------------------------------------  
c            if (ga_nodeid().eq.0) then
c             write(*,7) idir,iresp
c 7           format('==== compute alfa(',i3,',',i3,') ====== START')
c            endif  

            call get_alfaorbeta_reim(
     &            alfare(idir,iresp), ! out: alpha or beta real part
     &            alfaim(idir,iresp), ! out: alpha or beta im   part
     &            g_vecE1(ispin,1),   ! in : 1st-order pert vec RE
     &            g_vecE1_im(ispin,1),! in : 1st-order pert vec IM
     &            g_dipel,            ! in : dipole electric or magnetic
     &            g_vectors(ispin),   ! in : MO vectors
     &            idir,               ! in : = 1,2,3=x,y,z directions
     &            iresp,              ! in : = 1,2,3
     &            coeffre,coeffim,1,  ! in : (coeffre,coeffim,caseAO)
     &            nbf,                ! in : nr. basis functions
     &            nocct(ispin),       ! in : nr. occupied alpha or beta
     &            lifetime,           ! in : logical var for damping
     &            debug,              ! in : logical var for debugging
     &            g_temp)             ! in : scratch GA array 

c            if (ga_nodeid().eq.0) then
c             write(*,8) idir,iresp
c 8           format('==== compute alfa(',i3,',',i3,') ====== END')
c            endif      
            if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,9) ispin,idir,iresp,
     &                   alfare(idir,iresp),alfaim(idir,iresp)
 9            format('(alfa(',i3,',',i3,',',i3,
     &             ')=(',f15.8,',',f15.8,')')
             endif
            endif ! end-if-debug
          enddo                 ! iresp = 1,3
c         ------------------------------------------------
c         if requested by input, compute dipole-quadrupole
c         polarizability, and add to OR beta to obtain the 
c         Buckingham-Dunn OR tensor
c         ------------------------------------------------       
          if (lquad) then
           do iresp = 1,6 

c             if (ga_nodeid().eq.0) then
c              write(*,11) idir,iresp
c 11           format('==== compute dipquad(',i3,',',i3,') ====== START')
c             endif

             call get_alfaorbeta_reim(
     &            dipquadre(idir,iresp),! in/out: alpha or beta real part
     &            dipquadim(idir,iresp),! in/out: alpha or beta im   part
     &            g_vecE1(ispin,1),     ! in : 1st-order pert vec RE
     &            g_vecE1_im(ispin,1),  ! in : 1st-order pert vec IM
     &            g_quadel,             ! in : dipole electric or magnetic
     &            g_vectors(ispin),     ! in : MO vectors
     &            idir,                 ! in : = 1,2,3=x,y,z directions
     &            iresp,                ! in : = 1,2,3
     &            coeffre,coeffim,1,    ! in : (coeffre,coeffim,caseAO)
     &            nbf,                  ! in : nr. basis functions
     &            nocct(ispin),         ! in : nr. occupied alpha or beta
     &            lifetime,             ! in : logical var for damping
     &            debug,                ! in : logical var for debugging
     &            g_temp)               ! in : scratch GA array 

c            if (ga_nodeid().eq.0) then
c             write(*,50) idir,iresp
c 50          format('==== compute dipquad(',
c     &              i3,',',i3,') ====== END')
c            endif      
            if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,12) ispin,idir,iresp,
     &                dipquadre(idir,iresp),dipquadim(idir,iresp)
 12           format('Dipquad(',i3,',',i3,',',i3,
     &            ')=(',f15.8,',',f15.8,')') 
             endif
            endif ! end-if-debug
           enddo ! end-loop-iresp
          endif ! end-if-lquad
        enddo ! end-loop-idir (perturbing E-field components)

c ============ visualize (alfa,beta) ========== START
         if (debug) then
          if (ga_nodeid().eq.0) then
           do idir=1,ndir
            do iresp=1,3
             write(*,10) ispin,idir,iresp,
     &        alfare(idir,iresp),alfaim(idir,iresp),
     &        betare(idir,iresp),betaim(idir,iresp)
 10          format('NEW-beta:(ispin,idir,iresp)=(',i3,',',i3',',i3,')',
     &              ' alfa(re,im)=(',f15.8,',',f15.8,')',
     &              ' beta(re,im)=(',f15.8,',',f15.8,')')
            enddo ! end-loop-iresp
           enddo ! end-loop-idir
          endif
         endif ! end-if-debug
c ============ visualize (alfa,beta) ========== END
c       -------------------------------------------
c       -------------------------------------------------
c       un-add the frequency components in vec_E1 in case
c       we reuse these arrays:
c       -------------------------------------------------
        if (ncomp.gt.1) then
          call ga_add(1d0,g_vecE1(ispin,1),
     &               -1d0,g_vecE1(ispin,2),
     &                    g_vecE1(ispin,1))
          if (lifetime) then
            call ga_add(1d0,g_vecE1_im(ispin,1),
     &                 -1d0,g_vecE1_im(ispin,2),
     &                      g_vecE1_im(ispin,1))
          end if
        endif      
          if (.not.ga_destroy(g_temp))
     &       call errquit
     &       ('aor_beta: ga_destroy failed g_temp',
     &       0,GA_ERR) 
      enddo ! end-loop-ispin
c       -----------------
c       deallocate memory
c       -----------------
        if (.not.ga_destroy(g_work))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_work',
     &     0,GA_ERR)
c     it seems that if we use GIAOs everything is off by a factor of
c     two, so we need to scale betare, betaim. If we have static
c     response then there is a factor of two missing everywhere
c     because we don't add C(+) and C(-) for the electric field.

c$$$      if (lgiao) then
c$$$        scaling = half
c$$$        do idir = 1,3
c$$$          do iresp = 1,3
c$$$            betare(idir, iresp) = betare(idir, iresp) * scaling
c$$$            betaim(idir, iresp) = betaim(idir, iresp) * scaling
c$$$          end do
c$$$        end do
c$$$      end if                    ! lgiao
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
c        ... getting tired of that ol' f77 syntax ... 
        if (lquad) then
          dipquadre(:,:) = dipquadre(:,:) * scaling
          dipquadim(:,:) = dipquadim(:,:) * scaling
        end if
      end if                    ! lstatic
c     for comparison, let's print the
c     transpose of the negative OR beta tensor
      if (debug) then
      write (luout,*) 'DEBUG: OR BETA as in some other code'
      do idir = 1,ndir
        write (luout,'(5x,i1,1x,E13.6,1x,E13.6,1x,E13.6)')
     &    idir, -betare(1:3,idir)
      end do
      write (luout,*) 'iso: ', (betare(1,1)+betare(2,2)+betare(3,3))/3d0
      end if
c     ------------------------------------------------------------
c     if requested by input, construct traceless dipole-quadrupole
c     polarizability tensor and contruct Buckingham-Dunn OR tensor
c     At present we do only the real part !
c     ------------------------------------------------------------
      if (lquad) then
        if (debug) then
          write (luout,*) 'dip-quad polarizability:'
          call output (dipquadre,1,3,1,6,3,6,1)
        end if
        dqpol(:,:,:) = zero
        do i = 1,3
          do j = 1,3
            do k = 1,j
              dqpol(i,j,k) = three * dipquadre(i,qindex(j,k))
              if (k.ne.j)
     &           dqpol(i,k,j) = three * dipquadre(i,qindex(k,j))
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
        end do
        dqpol(:,:,:) = dqpol(:,:,:) * half
        if (debug) then
          write (luout,*) 'traceless dip-quad polarizability:'
          do i = 1,3
            call output (dqpol(i,1:3,1:3),1,3,1,3,3,3,1)
          end do
        end if

c       if we have the keyword for the ROA output, let's do that
c       before the BD tensor gets assembled:

        if (lroa .and. me) then
          do i = 1,3
            do j = 1,3
              tmpmat(i,j) = -omega*betare(i,j)
            end do
          end do
          call aor_write_roadata
     &       (luout, alfare, betare, tmpmat, dipquadre, dqpol)
        end if                 

c       assemble BD tensor in tmpmat.
c       step 1: symmetrize beta tensor

        tmpmat(:,:) = half * betare(:,:)
        do i = 1,3
          do j = 1,3
            tmpmat(i,j) = tmpmat(i,j) + half*betare(j,i)
          end do
        end do
c       step 2: add quadrupole contributions
        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                tmpmat(i,j) = tmpmat(i,j) 
     &             + half * third * LCtensor(i,k,l) * dqpol(k,l,j) 
     &             + half * third * LCTensor(j,k,l) * dqpol(k,l,i)
              end do
            end do
          end do
        end do
        betare(:,:) = tmpmat(:,:)
      end if                    ! lquad
c     ----------------
c     all done. return
c     ----------------
c     ==================================================================
      return
      end
