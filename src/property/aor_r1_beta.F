
      subroutine aor_r1_beta (rtdb, basis, geom, omega, lstatic, ncomp,
     &   g_smat0, g_sket1, g_vecB1,  g_dipel, g_quadel, g_vectors, froct
     &   ,nbf, nmo, nocct, nvirt, lgiao, lquad, lanalyze, lvelocity
     &   ,lifetime, lmagpert, g_vecE1, g_vecE1_im, alfare, alfaim,
     &   betare, betaim)

c $Id:  $
      
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

c     These are all input, too
      integer g_smat0, g_vectors(2), g_dipel, 
     &   g_sket1, g_quadel, g_vecB1, g_vecE1(2), g_vecE1_im(2)
      integer nfreq, response_order, nbf, nmo, ncomp
      integer nocct(2), nvirt(2)
      double precision froct(nbf,2)
      double precision gamwidth, omega
      logical lgiao, lquad, lanalyze, lvelocity, lifetime, lmagpert,
     &   lstatic

c     output:
      double precision alfare(3,3), alfaim(3,3)
      double precision betare(3,3), betaim(3,3)

c     ----------------
c     local variables:
c     ----------------

c     global array handles:
      
      integer
     &   g_work,
     &   g_temp

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

      character*(256) cstemp

      character*(1) direction(3)
      data direction/'x','y','z'/
      
      integer ispin, nspin
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, iresp, idir,
     &   i,j,k,l
      logical debug, dbgmat, 
     &   lzora, lantisym
      double precision sum, scaling
      double precision tenm8, one, two, three, zero, half, third
      parameter (tenm8=1d-8, one=1d0, two=2d0, three=3d0,
     &   zero=0d0, half=one/two,
     &   third=one/three)

c     external functions:

      double precision ga_trace_diag
      external ga_trace_diag

c  ====================================================================

      debug = .false. .and. ga_nodeid().eq.0 ! .true. during development
      dbgmat = .false. .and. ga_nodeid().eq.0 ! debug large matrices

      if (debug) write (luout,*) 'hello from aor_r1_beta'

c     the main results are collected in alfare/in, betare/im.
c     initialize these and other tensors with zeros:

      alfare(:,:) = 0
      alfaim(:,:) = 0
      betare(:,:) = 0
      betaim(:,:) = 0

      dipquadre(:,:) = 0
      dipquadim(:,:) = 0

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
      nspin      =  1           ! assume closed shell
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
c     start loop over spins (nspin=2 not yet supported)
c     ----------------------------------------------------
      
      do ispin = 1, nspin
        
        nmo1 = nmot(ispin)      ! total no.of MOs for this spin
        nocc = nocct(ispin)     ! occupied MOs
        nvir = nvirt(ispin)     ! virtual MOs
        nocv = nocvir(ispin)    ! nocc * nvir
        
        
        
c       ------------------------------
c       allocate some temp. work space
c       ------------------------------
        
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
        
c       allocate intermediate vector for matrix multiplications
c       used to create the final results
        
        write (cstemp,'(a)'), 'aor_beta: temp1'
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
          call ga_add(1d0, g_vecE1(1), 1d0,  g_vecE1(2),
     &       g_vecE1(1))
          if (lifetime) then
            call ga_add(1d0, g_vecE1_im(1), 1d0,  g_vecE1_im(2),
     &         g_vecE1_im(1))
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
            
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = 1 
            ahi(3) = 1
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nocc
            blo(3) = iresp      ! pick magnetic field direction
            bhi(3) = iresp 
            clo(1) = 1
            chi(1) = nbf
            clo(2) = 1
            chi(2) = nocc
            
            call ga_zero(g_temp)
            call nga_matmul_patch('n','n',1d0,0d0,
     &         g_smat0,alo,ahi,
     &         g_vecB1,blo,bhi,
     &         g_temp,clo,chi)
            
            if (debug) write (luout,*)
     &         'beta: S(0) C(B) intermediate complete'
            
            alo(1) = 1
            ahi(1) = nocc
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = idir
            ahi(3) = idir
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
                        
c           the factor of two is for the orbital occupations,
c           assuming that ispin is never equal to two
            
            sum = 2d0 * ga_trace_diag(g_work)
            
            betare(idir,iresp) = betare(idir,iresp) + sum

            if (lifetime) then

              call ga_zero(g_work)
              call nga_matmul_patch('t','n',1d0,0d0,
     &           g_vecE1_im,alo,ahi,
     &           g_temp,blo,bhi,
     &           g_work,clo,chi)              
            
              sum = 2d0 * ga_trace_diag(g_work)              
              
              betaim(idir,iresp) = betaim(idir,iresp) - sum

            end if              ! lifetime

            if (debug) write (luout,*) 'beta: C(E) S(0) C(B) complete'
            
            if (lgiao) then
              
c             --------------------------------------
c             if we use GIAOs there is a second term
c             in beta which is C(E) S(1ket) C(0)
c             --------------------------------------
              
              alo(1) = 1
              ahi(1) = nbf
              alo(2) = 1
              ahi(2) = nbf
              alo(3) = iresp    ! pick the correct sket1 direction
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
     &           g_sket1,alo,ahi,
     &           g_vectors(ispin),blo,bhi,
     &           g_temp,clo,chi)
              
              if (debug) write (luout,*)
     &           'beta: S(ket1) C(0) intermediate complete'
                            
              alo(1) = 1
              ahi(1) = nocc
              alo(2) = 1
              ahi(2) = nbf
              alo(3) = idir
              ahi(3) = idir
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
     &           g_vecE1,alo,ahi,
     &           g_temp,blo,bhi,
     &           g_work,clo,chi)
                            
c             the factor of two is for the orbital occupations,
c             assuming that ispin is never equal to two
              
              sum = 2d0 * ga_trace_diag(g_work)
              
              betare(idir,iresp) = betare(idir,iresp) + sum
              
              if (lifetime) then
                
                call ga_zero(g_work)
                call nga_matmul_patch('t','n',1d0,0d0,
     &             g_vecE1_im,alo,ahi,
     &             g_temp,blo,bhi,
     &             g_work,clo,chi)              
                
                sum = 2d0 * ga_trace_diag(g_work)              
                
                betaim(idir,iresp) = betaim(idir,iresp) - sum
                
              end if            ! lifetime
                           
              if (debug) write (luout,*)
     &           'beta: C(E) S(ket1) C(0) complete'

            end if              ! lgiao
            
c           --------------------------------------
c           (B) calculate alfa from C(E) h(E) C(0)
c           --------------------------------------
          
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = iresp      ! pick direction iresp for g_dipel
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
     &         g_dipel,alo,ahi,
     &         g_vectors(ispin),blo,bhi,
     &         g_temp,clo,chi)
            
            if (debug) write (luout,*)
     &         'alfa: h(E) C(0) intermediate complete'
            
            
            alo(1) = 1
            ahi(1) = nocc
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = idir
            ahi(3) = idir
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
                                    
c           the factor of two is for the orbital occupations,
c           assuming that ispin is never equal to two
            
            sum = 2d0 * ga_trace_diag(g_work)
            
            alfare(idir,iresp) = alfare(idir,iresp) - sum

c           repeat for imaginary part, if necessary:

            if (lifetime) then
              
              call ga_zero(g_work)
              call nga_matmul_patch('t','n',1d0,0d0,
     &           g_vecE1_im,alo,ahi,
     &           g_temp,blo,bhi,
     &           g_work,clo,chi)              
              
              sum = 2d0 * ga_trace_diag(g_work)              
              
              alfaim(idir,iresp) = alfaim(idir,iresp) + sum
              
            end if              ! lifetime     

            if (debug) write (luout,*) 'alfa C(E) h(E) C(0) complete'
            
          enddo                 ! iresp = 1,3

c         ------------------------------------------------
c         if requested by input, compute dipole-quadrupole
c         polarizability, and add to OR beta to obtain the 
c         Buckingham-Dunn OR tensor
c         ------------------------------------------------       

          if (lquad) then

            do iresp = 1,6 
              
              alo(1) = 1
              ahi(1) = nbf
              alo(2) = 1
              ahi(2) = nbf
              alo(3) = iresp    ! pick direction iresp for g_quadel
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
     &           g_quadel,alo,ahi,
     &           g_vectors(ispin),blo,bhi,
     &           g_temp,clo,chi)
              
              if (debug) write (luout,*)
     &           'quad: h(Q) C(0) intermediate complete'
              
              
              alo(1) = 1
              ahi(1) = nocc
              alo(2) = 1
              ahi(2) = nbf
              alo(3) = idir
              ahi(3) = idir
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
     &           g_vecE1,alo,ahi,
     &           g_temp,blo,bhi,
     &           g_work,clo,chi)
              
c             the factor of two is for the orbital occupations,
c             assuming that ispin is never equal to two
              
              sum = 2d0 * ga_trace_diag(g_work)
              
              dipquadre(idir,iresp) = dipquadre(idir,iresp) - sum
              
c             repeat for imaginary part, if necessary:
              
              if (lifetime) then
              
                call ga_zero(g_work)
                call nga_matmul_patch('t','n',1d0,0d0,
     &             g_vecE1_im,alo,ahi,
     &             g_temp,blo,bhi,
     &             g_work,clo,chi)              
                
                sum = 2d0 * ga_trace_diag(g_work)              
                
                dipquadim(idir,iresp) = dipquadim(idir,iresp) + sum
                
              end if            ! lifetime     
              
              if (debug) write (luout,*) 'quad C(Q) h(E) C(0) complete'
              
            end do              ! iresp

          end if                ! lquad

c         -----------------------------------------
c         end loop over responding field components
c         -----------------------------------------

          
        end do                  ! idir = 1,3

c       -------------------------------------------
c       end loop over perturbing E-field components
c       -------------------------------------------

c       -----------------
c       deallocate memory
c       -----------------

          if (.not.ga_destroy(g_temp))
     &       call errquit
     &       ('aor_beta: ga_destroy failed g_temp',
     &       0,GA_ERR)


        if (.not.ga_destroy(g_work))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_work',
     &     0,GA_ERR)

c       -------------------------------------------------
c       un-add the frequency components in vec_E1 in case
c       we reuse these arrays:
c       -------------------------------------------------

        if (ncomp.gt.1) then
          call ga_add(1d0, g_vecE1(1), -1d0,  g_vecE1(2),
     &       g_vecE1(1))
          if (lifetime) then
            call ga_add(1d0, g_vecE1_im(1), -1d0,  g_vecE1_im(2),
     &         g_vecE1_im(1))
          end if
        endif       

        
      enddo                     ! ispin = 1,2 from way above
              
c     ---------------------------------------------------------------
c     end loop over spin components (which we don't use right now
c     since nspin is forced to be 1 at the beginning of this routine)
c     ---------------------------------------------------------------
            

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
        do idir = 1,3
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
      do idir = 1,3
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

