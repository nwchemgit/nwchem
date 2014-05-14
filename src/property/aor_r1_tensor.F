      subroutine aor_r1_tensor(
     &      rtdb,basis,geom,    ! in : handles
     &      omega,              ! in :
     &      lstatic,            ! in :
     &      ncomp,              ! in :  
     &      g_smat0,            ! in :
     &      g_dipmag,           ! in : magn -dipole mom AO
     &      g_dipel,            ! in : elect-dipole mom AO
     &      g_quadel,           ! in : quadrupole   AO
     &      g_vectors,          ! in : MOs
     &      froct,              ! in : set of occupations
     &      nbf, nmo,           ! in : nr basis, nr MOs
     &      npol,               ! in : nr. polarizations
     &      nocct, nvirt,       ! in : nocc,nvirt
     &      lgiao, lquad,       ! in : logical vars
     &      lvelocity,          ! in : logical vars
     &      lifetime,           ! in : logical vars
     &      lmagpert,           ! in : logical vars
     &      g_vecE1,g_vecE1_im, ! in : 
     &      alfare,alfaim,      ! out: electric-electric response matrices
     &      betare,betaim)      ! out: electric-magnetic response matrices       
c $Id$     
c     =================================================================   
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
      integer npol,nocct(npol), nvirt(npol)
      double precision froct(nbf,npol)
c     These are all input, too
      integer g_smat0, g_vectors(npol), g_dipel, 
     &        g_quadel, g_dipmag, 
     &        g_vecE1(2,2),g_vecE1_im(2,2)
      integer nbf, nmo, ncomp
      double precision gamwidth, omega
      logical lgiao, lquad, lvelocity, lifetime, lmagpert, lstatic
      double precision sum_my ! Added by FA
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
      character*(256) cstemp
      character*(1) direction(3)
      data direction/'x','y','z'/
      integer ispin
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, 
     &        iresp, idir
      logical debug, dbgmat, 
     &        lzora, lantisym
      double precision sum, scaling
      double precision tenm8, one, two, zero, half
      parameter (tenm8=1d-8, one=1d0, two=2d0, zero=0d0, half=one/two)
c     external functions:
      double precision ga_trace_diag,coeffre,coeffim
      external ga_trace_diag 
c ------- Added for unrestricted calc ----- START
      integer ndir
      external get_alfaorbeta_reim
c ------- Added for unrestricted calc ----- END
      ndir=3 ! nr directions (x,y,z)
c  ====================================================================
      debug  = .false. .and. ga_nodeid().eq.0 ! .true. during development
      dbgmat = .false. .and. ga_nodeid().eq.0 ! debug large matrices

c      debug=.true.

      if (debug) write (luout,*) 'hello from aor_r1_beta'
c     the main results are collected in alfare/in, betare/im.
c     initialize with zeros:
      do idir=1,ndir
       do iresp=1,3
        alfare(idir,iresp) = 0.0d0
        alfaim(idir,iresp) = 0.0d0
        betare(idir,iresp) = 0.0d0
        betaim(idir,iresp) = 0.0d0
       enddo ! end-loop-iresp
      enddo ! end-loop-idir
c     set parameters that control the various computational options
c     (later we will set most of this by input)
      lzora  = .false.      ! not yet available here 
      if (debug) write (luout,*) 'giao, velocity, magpert',
     &    lgiao, lvelocity, lmagpert
      lantisym = (lvelocity .or. lmagpert) ! antisymm. perturbation
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
c     ----------------------
c     start loop over spins
c     ----------------------
      if      (npol.eq.1) then
       coeffre=-2.0d0
       coeffim=+2.0d0   
      else if (npol.eq.2) then
       coeffre=-0.5d0
       coeffim=+0.5d0  
      endif
      do ispin = 1, npol
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
c       -------------------------------------------------------
c       (A) calculate optical rotation beta from C(E) H(B) C(0)
c       ------------------------------------------------------
c       ---------------------------------------------------------
c       solution of CPKS is in g_vecE1. For the OR and a length-gauge
c       dipole response (E-field perturbation) we need the difference of
c       the +/- components so we subtract them here and store the result
c       in g_vecE1(1). If the velocity formalism is used then the
c       perturbing field has an antisymmetric AO matrix and 
c       we have to ADD the vectors instead. 
c       ---------------------------------------------------------
        if (debug) then
         if (ga_nodeid().eq.0) then
          write(*,100) ncomp,lantisym,lifetime
 100      format('(comp,lantisym,lifetime)=(',i3,',',L1,',',L1,')')
         endif
         if (ga_nodeid().eq.0)
     &    write(*,*) '--g_vecE1-BEF-IF ------ START'
          call ga_print(g_vecE1(ispin,1))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--g_vecE1-BEF-IF ------ END'
         if (lifetime) then
          if (ga_nodeid().eq.0)
     &      write(*,*) '--g_vecE1-im-BEF-IF ------ START'
            call ga_print(g_vecE1_im(ispin,1))
          if (ga_nodeid().eq.0)
     &      write(*,*) '--g_vecE1-im-BEF-IF ------ END'
         endif ! end-if-lifetime
        endif ! end-if-debug
        if (ncomp.gt.1) then
          if (lantisym) then
            call ga_add(1d0, g_vecE1(ispin,1), 
     &                  1d0, g_vecE1(ispin,2),
     &                       g_vecE1(ispin,1))
            if (lifetime) then
              call ga_add(1d0, g_vecE1_im(ispin,1), 
     &                    1d0, g_vecE1_im(ispin,2),
     &                         g_vecE1_im(ispin,1))
            end if              ! lifetime
          else                  ! lantisym ? 
            call ga_add(1d0, g_vecE1(ispin,1), 
     &                 -1d0, g_vecE1(ispin,2),
     &                       g_vecE1(ispin,1))
            if (lifetime) then
              call ga_add(1d0, g_vecE1_im(ispin,1),
     &                   -1d0, g_vecE1_im(ispin,2),
     &                         g_vecE1_im(ispin,1))
            end if              ! lifetime
          end if                ! lantisym
        endif                   ! ncomp.gt.1
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '--g_vecE1-AFT-IF ------ START'
          call ga_print(g_vecE1(ispin,1))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--g_vecE1-AFT-IF ------ END'
         if (lifetime) then
          if (ga_nodeid().eq.0)
     &     write(*,*) '--g_vecE1-im-AFT-IF ------ START'
           call ga_print(g_vecE1_im(ispin,1))
          if (ga_nodeid().eq.0)
     &     write(*,*) '--g_vecE1-im-AFT-IF ------ END'
         endif ! end-if-lifetime
        endif ! end-if-debug

       ndir=3 ! nr. directions (x,y,z)
       do idir = 1,ndir        ! direction of the perturbing field
        do iresp = 1,3

         call get_alfaorbeta_reim(
     &            betare(idir,iresp), ! in/out: alpha or beta real part
     &            betaim(idir,iresp), ! in/out: alpha or beta im   part
     &            g_vecE1(ispin,1),   ! in : 1st-order pert vec RE
     &            g_vecE1_im(ispin,1),! in : 1st-order pert vec IM
     &            g_dipmag,           ! in : dipole electric or magnetic
     &            g_vectors(ispin),   ! in : MO vectors
     &            idir,               ! in : = 1,2,3=x,y,z directions
     &            iresp,              ! in : = 1,2,3
     &            coeffre,coeffim,1,  ! in : (coeffre,coeffim,caseAO)
c     &            -coeffre,coeffim,1,  ! in : (coeffre,coeffim,caseAO)
     &            nbf,                ! in : nr. basis functions
     &            nocc,               ! in : nr. occupied alpha or beta
     &            lifetime,           ! in : logical var for damping
     &            debug,              ! in : logical var for debugging
     &            g_temp)             ! in : scratch GA array

         if (debug) then 
          if (ga_nodeid().eq.0) then
           write(*,1) ispin,idir,iresp,
     &                betare(idir,iresp),betaim(idir,iresp)
 1         format('beta(',i3,',',i3,',',i3,
     &            ')=(',f15.8,',',f15.8,')')
          endif
         endif ! end-if-debug

        enddo  ! end-loop-iresp (responding field components)
       enddo ! end-loop-idir    (perturbing E-field components)

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
     &                  2d0, g_vecE1(ispin,2),
     &                       g_vecE1(ispin,1))
            if (lifetime) then
              call ga_add(1d0, g_vecE1_im(ispin,1), 
     &                    2d0, g_vecE1_im(ispin,2),
     &                         g_vecE1_im(ispin,1))
            end if              ! lifetime
          end if                ! lantisym
        endif                   ! ncomp.gt.1

       ndir=3 ! nr. directions (x,y,z)
       do idir = 1,ndir        ! direction of the perturbing field
        do iresp = 1,3
         call get_alfaorbeta_reim(
     &            alfare(idir,iresp), ! out: alpha or beta real part
     &            alfaim(idir,iresp), ! out: alpha or beta im   part
     &            g_vecE1(ispin,1),   ! in : 1st-order pert vec RE
     &            g_vecE1_im(ispin,1),! in : 1st-order pert vec IM
     &            g_dipel,            ! in : dipole electric or magnetic
     &            g_vectors(ispin),   ! in : MO vectors
     &            idir,               ! in : = 1,2,3=x,y,z directions
     &            iresp,              ! in : = 1,2,3
     &            coeffre,-coeffim,1,  ! in : (coeffre,coeffim,caseAO) 
     &            nbf,                ! in : nr. basis functions
     &            nocc,               ! in : nr. occupied alpha or beta
     &            lifetime,           ! in : logical var for damping
     &            debug,              ! in : logical var for debugging
     &            g_temp)             ! in : scratch GA array

         if (debug) then 
          if (ga_nodeid().eq.0) then
           write(*,2) ispin,idir,iresp,
     &                alfare(idir,iresp),alfaim(idir,iresp)
 2         format('alfa(',i3,',',i3,',',i3,
     &            ')=(',f15.8,',',f15.8,')')
          endif
         endif ! end-if-debug

        enddo  ! end-loop-iresp (responding field components)
       enddo ! end-loop-idir    (perturbing E-field components)
c ============ visualize (alfa,beta) ========== START
         if (debug) then
          if (ga_nodeid().eq.0) then
           do idir=1,ndir
            do iresp=1,3
             write(*,10) ispin,idir,iresp,
     &        alfare(idir,iresp),alfaim(idir,iresp),
     &        betare(idir,iresp),betaim(idir,iresp)
 10          format('(ispin,idir,iresp)=(',i3,',',i3',',i3,')',
     &              ' alfa(re,im)=(',f15.8,',',f15.8,')',
     &              ' beta(re,im)=(',f15.8,',',f15.8,')')
            enddo ! end-loop-iresp
           enddo ! end-loop-idir
          endif
         endif ! end-if-debug
c ============ visualize (alfa,beta) ========== END
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

      if (lgiao) then
        scaling = half
        do idir = 1,ndir ! direction of the perturbing field (x,y,z)
          do iresp = 1,3
            betare(idir, iresp) = betare(idir, iresp) * scaling
            betaim(idir, iresp) = betaim(idir, iresp) * scaling
          end do
        end do
      end if                    ! lgiao

      if (lstatic) then
        scaling = two
        do idir = 1,ndir ! direction of the perturbing field (x,y,z)
          do iresp = 1,3
            alfare(idir, iresp) = alfare(idir, iresp) * scaling
            alfaim(idir, iresp) = alfaim(idir, iresp) * scaling
            betare(idir, iresp) = betare(idir, iresp) * scaling
            betaim(idir, iresp) = betaim(idir, iresp) * scaling
          end do
        end do
      end if                    ! lstatic

c     for magnetic perturbation, alfa is the G' tensor, and
c     beta is the paramagnetic mag-mag response. 
c     for electric field perturbations (default), alfa is the
c     polarizability and beta is G'. In the output, for nonzero
c     omega we want the beta tensor, not G'. 
c     also consider the case of velocity formalisms where 
c     additional divisions by omega occur. 

      do idir = 1,ndir ! direction of the perturbing field (x,y,z)
        do iresp = 1,3
          if (lmagpert .and. .not.lstatic) then
c           case I: magnetic perturbation. alpha = G', beta = Chi-p
              alfare(idir,iresp) = alfare(idir,iresp) / omega
              alfaim(idir,iresp) = alfaim(idir,iresp) / omega
            if (lvelocity) then
              alfare(idir,iresp) = alfare(idir,iresp) / omega
              alfaim(idir,iresp) = alfaim(idir,iresp) / omega
            endif
          else if (.not.lmagpert .and. .not.lstatic) then
c           case II: alpha = polarizability, beta = G'
              betare(idir,iresp) = betare(idir,iresp) / omega
              betaim(idir,iresp) = betaim(idir,iresp) / omega
            if (lvelocity) then
              betare(idir,iresp) = betare(idir,iresp) / omega
              alfare(idir,iresp) = alfare(idir,iresp) / omega**2
              betaim(idir,iresp) = betaim(idir,iresp) / omega
              alfaim(idir,iresp) = alfaim(idir,iresp) / omega**2
            endif
          endif
        enddo
      enddo
c     ----------------
c     all done. return
c     ----------------
      return
      end
