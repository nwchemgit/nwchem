
      subroutine aor_r1_beta_anl (rtdb, basis, geom, omega, lstatic,
     &   ncomp,g_smat0, g_sket1, g_vecB1,  g_dipel, g_quadel, g_vectors,
     &   froct,nbf, nmo, nocct, nvirt, lgiao, lquad, lanalyze, lvelocity
     &   ,lifetime, lmagpert, g_vecE1, g_vecE1_im)

c $Id:  $
      
c     =================================================================
      
c     purpose: analyze beta tensor. See routine
c     aor_r1_beta.F for additional comments. 
c     Use a molecular orientation in which the chiral
c     response tensor is diagonal!

c     called from: aoresponse_driver_new

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

c     ----------------
c     local variables:
c     ----------------

c     global array handles:
      
      integer
     &   g_anl, g_work,
     &   g_temp, g_tmpanl, g_tran

      integer l_diag, k_diag

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
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, idir, iresp,
     &   i,j,k,l
      logical debug, dbgmat, 
     &   lzora, lantisym, lmo, status, oprint
      double precision sum, scaling
      double precision tenm8, one, two, three, zero, half, third, four
      parameter (tenm8=1d-8, one=1d0, two=2d0, three=3d0,
     &   zero=0d0, half=one/two,
     &   third=one/three, four=4.0d0)

c     external functions:

      character*(256) lmotrans
      logical file_read_ga
      external file_read_ga

      double precision ga_trace_diag
      external ga_trace_diag

c  ====================================================================

      debug = .false.  ! .true. during development
      dbgmat = .false. ! debug large matrices
      oprint = ga_nodeid().eq.0

      if (debug) write (luout,*) 'hello from aor_r1_beta_anl'

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

c     ----------------------------------------------
c     check if we have localized MOs on file. If yes
c     read them, assuming nspin.eq.1
c     ----------------------------------------------

      i = 1
      status=rtdb_get(rtdb,'prop:pmlocalization',MT_INT,1,i)
      lmo = (i.eq.0)
      if (lmo) then
        if (oprint) write (luout,*) 'analysis: LMO switch found'
        write (cstemp,'(a)'), 'aor_beta: g_tran'
        if(.not.ga_create(MT_DBL, nocct(1), nocct(1), trim(cstemp),
     &     -1,-1,g_tran))
     &     call errquit (trim(cstemp),0,GA_ERR)
        if (debug) write (luout,*) 'g_anl allocated'
        call util_file_name('lmotrans',.true.,.true.,lmotrans)
        if(.not.file_read_ga(lmotrans,g_tran)) call errquit
     $     ('aor_r1_beta_anl: could not read lmotrans',0, DISK_ERR)
      end if ! lmo

c     ----------------------------------------------------
c     start loop over spins (nspin=2 not yet supported !!!)
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

c       allocate matrix that accumulates the analysis data

        write (cstemp,'(a)'), 'aor_beta: g_anl'
        if(.not.ga_create(MT_DBL, nocc, nocc, trim(cstemp),
     &     -1,-1,g_anl))
     &     call errquit (trim(cstemp),0,GA_ERR)
        if (debug) write (luout,*) 'g_anl allocated'
        call ga_zero(g_anl)

c       diagonal elements of the last matrix
        if (.not. ma_push_get(mt_dbl, nocc, 'diag', l_diag, k_diag))
     &        call errquit('error alloc MA diag', 0, MA_ERR)

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
          call ga_add(1d0, g_vecE1(1), 1d0,  g_vecE1(2),
     &       g_vecE1(1))
          if (lifetime) then
            call ga_add(1d0, g_vecE1_im(1), 1d0,  g_vecE1_im(2),
     &         g_vecE1_im(1))
          end if
        endif    

c       ------------------------------------------------
c       for Buckingham-Dunn tensor we need the traceless
c       quadrupole tensor
c       ------------------------------------------------

        if (lquad) then
          
c         a) from trace -> g_work
          call ga_zero(g_work)  ! use for trace
          do iresp = 1,3
            
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = qindex(iresp,iresp) 
            ahi(3) = qindex(iresp,iresp)
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nbf
            
            call nga_add_patch(one,g_quadel,alo,ahi,
     &         one,g_work,blo,bhi,g_work,blo,bhi)
          end do                ! iresp
          
c         b) scale quadel by 3
          call ga_scale(g_quadel,three)
          
c         c) subtract trace from diagonal
          
          do iresp = 1,3
            
            alo(1) = 1
            ahi(1) = nbf
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = qindex(iresp,iresp) 
            ahi(3) = qindex(iresp,iresp)
            blo(1) = 1
            bhi(1) = nbf
            blo(2) = 1
            bhi(2) = nbf
            
            call nga_add_patch(one,g_quadel,alo,ahi,
     &         -one,g_work,blo,bhi,g_quadel,alo,ahi)
          end do                ! iresp
          
c         d) divide the result by two, then by three
c         because of the factor by which the result enters
c         the Buckingham-Dunn tensor
          
          call ga_scale(g_quadel, half*third)
          
!          if (debug) call ga_print(g_quadel)
        end if                  ! lquad
        
c       ---------------------------------------------------------
c       start loop over the components of the response tensor and 
c       calculate the final results
c       ---------------------------------------------------------
        
        do idir = 1,3           ! direction of the perturbing field

c         g_anl is going to accumulate the results
          call ga_zero(g_anl)

          if (oprint)
     &       write (luout,'(1x,40(''-'')/1x,a,2i1)')
     &       'MO analysis of OR tensor component ',idir,idir
                                   
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
            blo(3) = idir      ! pick magnetic field direction
            bhi(3) = idir 
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
            
            scaling = two
            if (lstatic) scaling = four
            call nga_add_patch(scaling,g_work,clo,chi,
     &         one,g_anl,clo,chi,g_anl,clo,chi)

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
              alo(3) = idir    ! pick the correct sket1 direction
              ahi(3) = idir
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

              scaling = two
              if (lstatic) scaling = four
              call nga_add_patch(scaling,g_work,clo,chi,
     &           one,g_anl,clo,chi,g_anl,clo,chi)
              
              if (debug) write (luout,*)
     &           'beta: C(E) S(ket1) C(0) complete'

            end if              ! lgiao
            


c           ----------------------------------------------------
c           if requested by input, add to OR beta the quadrupole
c           polarizability terms
c           ----------------------------------------------------       
            
            if (lquad) then
              
c             add the quadrupole terms to the B.-D. tensor
              do k = 1,3
                do l = 1,3
                  
                  if (LCtensor(idir,k,l).eq.0) goto 1000
                  
                  iresp = qindex(l,idir)
                  
                  alo(1) = 1
                  ahi(1) = nbf
                  alo(2) = 1
                  ahi(2) = nbf
                  alo(3) = iresp ! pick direction iresp for g_quadel
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
     &               g_quadel,alo,ahi,
     &               g_vectors(ispin),blo,bhi,
     &               g_temp,clo,chi)
                  
                  if (debug) write (luout,*)
     &               'quad: h(Q) C(0) intermediate complete'
                                
                  alo(1) = 1
                  ahi(1) = nocc
                  alo(2) = 1
                  ahi(2) = nbf
                  alo(3) = k ! not idir, see equation
                  ahi(3) = k
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
     &               g_vecE1,alo,ahi,
     &               g_temp,blo,bhi,
     &               g_work,clo,chi)
                  
c                 the factor of two is for the orbital occupations,
c                 assuming that ispin is never equal to two
                  
                  scaling = -two
                  if (lstatic) scaling = -four
                  
c                 Levi-Civita symbol:
                  scaling = scaling * LCtensor(idir,k,l)
                  
                  if (debug) write (luout,*) 'scaling=',scaling
                  call nga_add_patch(scaling,g_work,clo,chi,
     &               one,g_anl,clo,chi,g_anl,clo,chi)
                
                if (debug) write (luout,*)
     &             'quad C(Q) h(E) C(0) complete'

 1000           continue ! jump here if LCTensor(idir,k,l) is 0
                
              end do            ! l
            end do              ! k

          end if                ! lquad

c         -----------------------------------------
c         end loop over responding field components
c         -----------------------------------------

c         ---------------------
c         Canonical MO analysis
c         ---------------------

          if (oprint) write (luout,
     &       '(/t12,a,t26,a/t11,6(''-''),t22,12(''-''))')
     &       'CMO #','contrib.'

          call ga_get_diagonal(g_anl, dbl_mb(k_diag) )

          sum = zero
          do i = 1,nocc
            sum = sum + dbl_mb(k_diag+i-1)
            if (oprint) write (luout,'(t11,i6,t22,f12.4)')
     &         i,dbl_mb(k_diag+i-1)
          end do                ! i
          if (oprint)
     &       write (luout,'(1x,a,2i1,a,f12.4)') 'Component ',idir,idir,
     &       ': Sum = ',sum

          if (oprint) write (luout,'(1x,40(''-''))')

          if (debug) then
            sum = ga_trace_diag(g_anl)
            if (oprint) write (luout,*) 'sum from ga_trace: ',sum
          end if                ! debug

c         ---------------------
c         Localized MO analysis
c         ---------------------

          if (lmo) then

            write (cstemp,'(a)'), 'aor_beta: tmpanl'
            if(.not.ga_create(MT_DBL, nocc, nocc, trim(cstemp),
     &         -1,-1,g_tmpanl))
     &         call errquit (trim(cstemp),0,GA_ERR)
            if (debug) write (luout,*) 'g_tmpanl allocated'
            
            call ga_zero(g_tmpanl)
            call ga_dgemm('t', 'n', nocc, nocc, nocc, 
     $         1.0d0, g_tran, g_anl, 0.0d0, g_tmpanl)
            call ga_zero(g_anl)
            call ga_dgemm('n', 'n', nocc, nocc, nocc, 
     $         1.0d0, g_tmpanl, g_tran, 0.0d0, g_anl)

          if (.not.ga_destroy(g_tmpanl))
     &         call errquit
     &         ('aor_beta: ga_destroy failed g_tmpanl',
     &         0,GA_ERR)

            if (oprint) write (luout,
     &         '(/t12,a,t26,a/t11,6(''-''),t22,12(''-''))')
     &         'LMO #','contrib.'
            
            call ga_get_diagonal(g_anl, dbl_mb(k_diag) )
            
            sum = zero
            do i = 1,nocc
              sum = sum + dbl_mb(k_diag+i-1)
              if (oprint) write (luout,'(t11,i6,t22,f12.4)')
     &           i,dbl_mb(k_diag+i-1)
            end do              ! i
            if (oprint)
     &         write (luout,'(1x,a,2i1,a,f12.4)')
     &         'Component ',idir,idir,': Sum = ',sum

            if (oprint) write (luout,'(1x,40(''-''))')
            
            if (debug) then
              sum = ga_trace_diag(g_anl)
              if (oprint) write (luout,*) 'sum from ga_trace: ',sum
            end if              ! debug
            
            
          end if                ! lmo
          
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

          if (.not.ga_destroy(g_anl))
     &       call errquit
     &       ('aor_beta: ga_destroy failed g_anl',
     &       0,GA_ERR)

         if (.not. ma_pop_stack(l_diag))
     &       call errquit('error deloc MA diag',0, MA_ERR)

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
            

      if (lmo) then
        if (.not.ga_destroy(g_tran))
     &       call errquit
     &       ('aor_beta: ga_destroy failed g_tran',0,GA_ERR)
      end if

c     ----------------
c     all done. return
c     ----------------
                  
      
c     ==================================================================
      
      return
      
      end

