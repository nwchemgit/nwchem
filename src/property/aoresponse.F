
      subroutine aoresponse (rtdb, basis, geom, nfreq,
     &   response_order, frequencies,
     &   g_smat0, g_smat1, g_sket1, g_dipel, g_dipmag, g_vectors, epst,
     &   froct, g_dens, nbf, nbfx, nmo, nocct, nvirt, deps,
     &   lgprime, lgiao, lvelocity, lmagpert, lifetime, gamwidth,
     &   alfare, alfaim, betare, betaim)

c     $Id$
      
c     =================================================================
      
c     purpose: perform response calculation

c     called from: aoresponse_driver

c     output: alfare, alfaim - electric response matrices
c             betare, betaim - magnetic response matrices

c     NOTE: the GIAO functionality is not yet working because the
c           subroutine that calculates the perturbed density matrix
c           incorrectly assumes that S(1,0) = 0 and therefore
c           assumes that A(+)dagger = A(-) which yields incorrect
c           results.

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
#include "util.fh"
c#include "dimqm.fh"

c     ---------------------
c     subroutine arguments:
c     ---------------------

      integer rtdb    ! [input] run-time database handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle

      integer g_smat0, g_dens(3), g_vectors(2), g_dipel, g_epsdif,
     &   g_dipmag, g_smat1, g_sket1
      integer nfreq, response_order, nbf, nbfx, nmo
      integer nocct(2), nvirt(2)
      double precision froct(nbf,2), epst(nbf,2)
      double precision frequencies(nfreq)
      double precision deps(nvirt(1),nocct(1))
      double precision alfare(3,3), alfaim(3,3)
      double precision betare(3,3), betaim(3,3)
      double precision gamwidth
      logical lgprime,lgiao, lvelocity, lmagpert, lifetime


c     ----------------
c     local variables:
c     ----------------

c     global array handles:
      
      integer g_amat(2), g_pmats(2),
     &   g_pmata(2), g_work, g_occ, g_rhs(2), g_rhs_zero, 
     &   g_amat_zero, g_pmat_zero, g_work_im ! DIM/QM JEM

c     GA arrays for treating imaginary parts of response
c     (only used if lifetime.eq.true.):

      integer g_amat_im(2), g_pmats_im(2),  g_pmata_im(2),
     &        g_rhs_im(2), g_rhs_zero_im, g_amat_zero_im,
     +        g_pmat_zero_im

c     other local variables: 

      integer nmot(2), nocvir(2)

      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      
      integer LCTensor(3,3,3)
      double precision origin(3)
      data origin/0.,0.,0./
      logical oskel
      parameter (oskel = .false.)

      double precision alfa0re(3,3), alfa0im(3,3), beta0re(3,3),
     &                 beta0im(3,3)

      character*(255) cstemp

      character*(1) direction(3)
      data direction/'x','y','z'/
      
      integer ispin, nspin, ncomp
      integer ipm, nocc, nvir, nocv, imo, jmo, nmo1, iresp, idir
      logical debug, dbgmat, 
     &   lzora, lantisym, lstatic, haveocc, limag, lzero
      logical oprint
      double precision omega, sum, scaling
      
      double precision tenm8, one, two, three, zero, half, third      
      parameter (tenm8=1d-8, one=1d0, two=2d0, three=3.0d0,
     &  zero=0d0, half=one/two,
     &  third=one/three)      

c     nwchem file access character strings:

      character*256 cphf_rhs(2), cphf_sol(2), cphf_rhs_im(2),
     &   cphf_sol_im(2)

      logical lclfld,use_dimqm,ldimqm,lrsp
      integer g_dipel_i
c     external functions:

      double precision ga_trace_diag
      external ga_trace_diag

      logical  cphf2, cphf3, cphf4, file_write_ga, file_read_ga, cphf
      external cphf2, cphf3, cphf4, file_write_ga, file_read_ga, cphf

c  ====================================================================

      oprint = util_print('information', print_low)
      oprint = oprint .and. (ga_nodeid().eq.0)
      use_dimqm=util_module_avail("dimqm")

c     -------------------
c     determine frequency
c     -------------------

      omega = frequencies(response_order)
      lstatic = (abs(omega).lt.tenm8) ! static response or not
      ncomp = 2                 ! no. of Fourier componnts to treat
      if (lstatic) ncomp = 1    ! treat only one component for static
      
c     -------------------------
c     set-up for CPKS procedure
c     -------------------------

      if (oprint) write (LuOut,'(a,i3,a/2x,a,E15.7,a/)')
     &   ' Performing order ',response_order, ' CPKS',
     &   ' with frequency omega = ',omega,' a.u.'
      if (lstatic .and. oprint) write (LuOut,'(a/)')
     &   ' STATIC response' 

c     set parameters that control the various computational options
c     (later we will set most of this by input)
      nspin      =  1           ! assume closed shell
      debug      = (.false. .and. oprint) ! produce debug output
      dbgmat     = .false.      ! debug large matrices
      lzora      = .false.      ! not yet available in nwchem

      haveocc = (lgiao .and. lmagpert) ! consider occ-occ blocks

c     check if perturbation orperator in CPKS is imaginary (antisymm.)
      limag = lmagpert .or. lvelocity

c     check need for computing zero-frequency response in addition
c     to the input-requested response (only for omega.ne.0)

      lzero = (lvelocity .and. .not.lstatic)
      if (debug) write (6,*) 'lzero:',lzero
      if (lzero .and. (oprint)) write (luout,*)
     & 'velocity gauge: will subtract the zero-frequency response' 

      if (debug) write (6,*) 'giao, velocity, magpert',
     &    lgiao, lvelocity, lmagpert
      if (debug) write (6,*) 'haveocc, limag:',haveocc,limag
      
c     -------------------------
c     define Levi-Civita tensor (not needed right now)
c     -------------------------
c      LCtensor = 0      
c      LCtensor(1,2,3) = 1
c      LCtensor(2,3,1) = 1
c      LCtensor(3,1,2) = 1      
c      LCtensor(2,1,3) = -1
c      LCtensor(1,3,2) = -1
c      LCtensor(3,2,1) = -1        

c     ---------------------------------------------------------------
c     allocate memory for perturbed density matrices (+/- components)
c     ---------------------------------------------------------------

      chunk(1) = nbf
      chunk(2) = -1
      chunk(3) = -1
      dims(1) = nbf
      dims(2) = nbf
      dims(3) = 3

      do ipm = 1,ncomp
        write(cstemp,'(a,i1)') 'pmats_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmats(ipm))) call 
     &     errquit('aoresponse: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmats(ipm))
        write(cstemp,'(a,i1)') 'pmata_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &       g_pmata(ipm))) call 
     &       errquit('aoresponse: nga_create failed '//cstemp(1:7),
     &       0,GA_ERR)
        call ga_zero(g_pmata(ipm))

        if (lifetime) then
          write(cstemp,'(a,i1)') 'pmats_im_',ipm
          if (.not.nga_create(MT_DBL,2,dims,cstemp(1:10),chunk,
     &       g_pmats_im(ipm))) call 
     &       errquit('aoresponse: nga_create failed '//cstemp(1:10),
     &       0,GA_ERR)
          call ga_zero(g_pmats_im(ipm))
          write(cstemp,'(a,i1)') 'pmata_im_',ipm
          if (.not.nga_create(MT_DBL,2,dims,cstemp(1:10),chunk,
     &       g_pmata_im(ipm))) call 
     &       errquit('aoresponse: nga_create failed '//cstemp(1:10),
     &       0,GA_ERR)
          call ga_zero(g_pmata_im(ipm))
        endif                   ! lifetime
      enddo

      if (lzero) then ! allocate arrays to hold zeroth-order response
        write(cstemp,'(a)') 'pmat_zero'
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:9),chunk,
     &     g_pmat_zero)) call 
     &     errquit('aoresponse: nga_create failed '//cstemp(1:9),
     &     0,GA_ERR)
        call ga_zero(g_pmat_zero)
        if (lifetime) then
          write(cstemp,'(a)') 'pmat_zero_im'
          if (.not.nga_create(MT_DBL,2,dims,cstemp(1:12),chunk,
     &       g_pmat_zero_im)) call 
     &       errquit('aoresponse: nga_create failed '//cstemp(1:12),
     &       0,GA_ERR)
          call ga_zero(g_pmat_zero_im)
        endif                   ! lifetime
      endif                     ! lzero

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
c
c      DIM/QM JEM
c      Allocate imaginary array for DIM with local fields
c      Normally, the imaginary part of the A matrix is zero, but
c      this is not true when working with local fields
      if (use_dimqm) then
         call dimqm_used(ldimqm)
         call dimqm_getlrsp(lrsp)
         call dimqm_getlclfld(lclfld)
c jbecca: START: Botched this call in util, because I didn't like it 
c           being a global variable. Could fix this later, but no DIM
c           parts are tested against this driver anyway.
c         call dimqm_getgdipeli(g_dipel_i)
      else
         ldimqm=.false.
         lrsp=.false.
         lclfld=.false.
         g_dipel_i=999
      endif
        if(ldimqm .and. lrsp .and. lclfld .and. lifetime) then
          if (.not.nga_create(MT_DBL,3,dims,'g_work_im',chunk,
     &     g_work_im)) call
     &     errquit('aoresponse: nga_create failed: g_work_im',
     &     0,GA_ERR)
          call ga_zero (g_work_im)

        end if

        
c       -----------------------------------------------------
c       GA-allocate components of A-matrices (+/-, Re and Im)
c       -----------------------------------------------------
        
        chunk(1) = nmo          ! distribution over processors
        chunk(2) = -1
        chunk(3) = -1
        dims(1) = nmo           ! dimensions of the array 
        dims(2) = nocc
        dims(3) = 3
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

        if (lzero) then
          write(cstemp,'(a)') 'amat_zero'
          if (debug) write (luout,*) cstemp(1:9)
          if (.not.nga_create(MT_DBL,3,dims,cstemp(1:9),chunk,
     &       g_amat_zero)) call 
     &       errquit('aoresponse: nga_create failed'//cstemp(1:9),
     &       0,GA_ERR)
          call ga_zero(g_amat_zero)

          if (lifetime) then
            write(cstemp,'(a)') 'amat_zero_im'
            if (debug) write (luout,*) cstemp(1:12)
            if (.not.nga_create(MT_DBL,3,dims,cstemp(1:12),chunk,
     &         g_amat_zero_im)) call 
     &         errquit('aoresponse: nga_create failed'//cstemp(1:12),
     &         0,GA_ERR)
            call ga_zero(g_amat_zero_im)
          endif                 ! lifetime
        endif                   ! lzero
        
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
        
        if (lmagpert) then
          call nga_copy_patch('n',g_dipmag,blo,bhi,g_work,alo,ahi)
        else
          call nga_copy_patch('n',g_dipel,blo,bhi,g_work,alo,ahi)
c       DIM/QM JEM
          if(ldimqm .and. lrsp .and. lclfld .and. lifetime) then
           write(luout,*) "Copy g_dipel_i"
           call nga_copy_patch('n',g_dipel_i,blo,bhi,g_work_im,alo,ahi)
          end if
        endif
        
c       ------------------------------------------------
c       transform to MO basis and extract vir-occ block,
c       store in g_amat's vir - occ block
c       ------------------------------------------------     
        
        call giao_aotomo(g_work,g_vectors(1),nocct(1),
     &     nvirt(1),1,3,nbf)   
c   DIM/QM JEM
        if(ldimqm .and. lrsp .and. lclfld .and. lifetime) then
          call giao_aotomo(g_work_im,g_vectors(1),nocct(1),
     &       nvirt(1),1,3,nbf)
        end if

        
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
c         DIM/QM JEM
          if(ldimqm .and. lrsp .and. lclfld .and. lifetime) then
            call nga_copy_patch('n',g_work_im,alo, ahi,
     &         g_amat_im(ipm),blo,bhi)
          end if

          
        enddo                   ! ipm = 1,ncomp

        if (lzero) then
          call nga_copy_patch('n',g_work,alo, ahi,
     &       g_amat_zero,blo,bhi)
        endif
        
        if (debug) write (luout,*) 'amat ao2mo complete'
        
c       we don't need this work array for a while
        
        if (.not.ga_destroy(g_work))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_work',
     &     0,GA_ERR)
c       DIM/QM JEM
        if(ldimqm .and. lrsp .and. lclfld .and. lifetime) then
          if (.not.ga_destroy(g_work_im))
     &       call
     &       errquit('aoresponse: ga_destroy failed g_work_im',
     &       0,GA_ERR)
        end if

        
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

c       for velocity formalism we will have to subtract the zero-
c       frequency limit

        if (lzero) then
          if(.not.ga_create(MT_DBL,nvir*nocc,3,'rhs_zero',
     &       -1,-1,g_rhs_zero))
     &       call errquit
     &       ('aoresponse: ga_create failed rhs_zero',0,GA_ERR)
          call ga_zero(g_rhs_zero)
          if (lifetime) then
            if(.not.ga_create(MT_DBL,nvir*nocc,3,'rhs_zero_im',
     &         -1,-1,g_rhs_zero_im))
     &         call errquit
     &         ('aoresponse: ga_create failed rhs_zero_im',0,GA_ERR)
            call ga_zero(g_rhs_zero_im)
          endif
        end if
c       
c       In a GIAO computation with the magnetic-field perturbing,
c       first evaluate the constant
c       part of the rhs other than the derivatives of the one-
c       electron part of the Fock operator (now in amat).
c       These components are the same for ipm = 1,2 so we just
c       do this once and duplicate, if necessary.
c       
        if (lmagpert .and. lgiao) then
          call aoresponse_giao_rhs (rtdb, basis, geom, ncomp, g_rhs,
     &       g_vectors(ispin), g_smat1, g_sket1, g_dens(1), g_amat,
     &       nocct(ispin), nvirt(ispin), nmot(ispin), nbf,
     &       froct(1,ispin), epst(1,ispin), omega, lstatic)
        endif

c       debug:
        haveocc = .false.
        
c       add amat vir-occ block to CPKF rhs
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
c   DIM/QM JEM
          if(ldimqm .and. lrsp .and. lclfld .and. lifetime) then
          call nga_add_patch(1d0,g_amat_im(ipm),alo,ahi,
     &       1d0, g_rhs_im(ipm),blo,bhi, g_rhs_im(ipm), blo, bhi)
          end if
        end do

        if (lzero) then
          call nga_add_patch(1d0,g_amat_zero,alo,ahi,
     &       1d0, g_rhs_zero,blo,bhi, g_rhs_zero, blo, bhi)
        endif
        
c       there is a division by -4 somewhere in the
c       CPKS solver so we scale the rhs by -4 in order to
c       get the correct magnitude of the result back
        
        do ipm = 1,ncomp
          call ga_scale (g_rhs(ipm), -4.0d0)
c   DIM/QM JEM
          if(ldimqm .and. lrsp .and. lclfld .and. lifetime) then
            call ga_scale (g_rhs_im(ipm), -4.0d0)
          end if
        enddo
        if (lzero) call ga_scale (g_rhs_zero, -4.0d0)

c       ... jochen: TEST
c        call ga_scale(g_rhs(1),-1.0d0)
 
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
        
c       We do need to tell the CPHF that the density is
c       might be skew symmetric here (omega.ne.0 only;, the dynamic solver
c       handles the +/- components of Pmat explicitly)
c       Done via rtdb, put cphf:skew .true. on rtdb and later
c       remove it.
c       should also work for velocity perturbation
c      
        if (lmagpert .or. lvelocity) then 
          if (.not.rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     &       errquit('aoresponse: failed to write skew ', 0, RTDB_ERR)
        else
          if (.not.rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.true.)) call 
     &       errquit('aoresponse: failed to write skew ', 0, RTDB_ERR) 
        endif

c       we also need to tell the dynamic CPKS solver if the perturbation
c       is imaginary and how many components to treat
        
        if (.not. rtdb_put(rtdb, 'cphf:imag', mt_log, 1, limag)) call
     &     errquit('aoresponse: failed to write limag ', 0, RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'cphf:ncomp', mt_int, 1, ncomp)) call
     &     errquit('aoresponse: failed to write ncomp ', 0, RTDB_ERR)        
c       
c       ----------------------------------------------------
c       Call the CPKS solver. It will solve the sets of
c       equations for each component of the perturbing field
c       simultaneously.
c       ----------------------------------------------------
c       
        if (debug) write (6,*) 'aoresponse: calling CPKS' 

        if (lstatic .and. .not.lifetime) then
c         ... static CPKS solver:         
          if (.not.cphf2(rtdb)) call errquit
     &       ('aoresponse: failure in cphf2 ',0, RTDB_ERR)
        else       
c         ... dynamic CPKS solver:
          if (debug) write (6,*) 'calling cphf3',omega,lifetime,gamwidth
          if (.not.cphf3(rtdb, omega, lifetime, gamwidth))
     +       call errquit
     $       ('aoresponse: failure in cphf3 ',0, RTDB_ERR)
        endif
        
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

c       ----------------------------------------------------
c       redo the CPKS for zero frequency in dynamic response 
c       based on velocity formalism
c       ----------------------------------------------------

        if (lzero) then
          call cphf_fname('cphf_rhs',cphf_rhs(1))
          call cphf_fname('cphf_sol',cphf_sol(1))
          if(.not.file_write_ga(cphf_rhs(1),g_rhs_zero)) call errquit
     $       ('aoresponse: could not write cphf_rhs zero',0, DISK_ERR)

          if (lifetime) then
            call cphf_fname('cphf_rhs_im',cphf_rhs_im(1))
            call cphf_fname('cphf_sol_im',cphf_sol_im(1))
            if(.not.file_write_ga(cphf_rhs_im(1),g_rhs_zero_im))
     &         call errquit
     $         ('aoresponse: could not write cphf_rhs zero Im',0,
     &         DISK_ERR)
          endif                 ! lifetime
          
c         for velocity calculation the perturbation is always
c         imaginary no matter if it is electric or magnetic field
          if (.not.rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call 
     &       errquit('aoresponse: failed to write skew ', 0, RTDB_ERR)
          if (.not. rtdb_put(rtdb, 'cphf:imag', mt_log, 1, limag)) call
     &       errquit('aoresponse: failed to write limag ', 0, RTDB_ERR)

c         also, we need to let the dynamic solver know that this is
c         a one-component situation (needed if lifetime.eq.true.)
          if (.not.rtdb_put(rtdb, 'cphf:ncomp', mt_int, 1,1)) call 
     &       errquit('aoresponse: failed to write ncomp ', 0, RTDB_ERR)
          
c         call static CPKS solver if no lifetime key is present,
c         otherwise use dynamic solver but set omega to zero:
          if (.not.lifetime) then
            if (.not.cphf2(rtdb)) call errquit
     &         ('aoresponse: failure in cphf2 zero',0, RTDB_ERR)
          else
            if (.not.cphf3(rtdb, 0d0, lifetime, gamwidth)) call errquit
     &         ('aoresponse: failure in cphf3 zero',0, RTDB_ERR)
          endif
          
          if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $       errquit('aoresponse: failed to delete skew', 0, RTDB_ERR)

          call ga_zero(g_rhs_zero)
          if(.not.file_read_ga(cphf_sol(1),g_rhs_zero))
     &       call errquit
     $       ('aoresponse: could not read cphf zero solution',
     &       0, DISK_ERR)          
          call nga_copy_patch('n',g_rhs_zero,blo,bhi,
     &       g_amat_zero,alo,ahi)

          if (lifetime) then
            call ga_zero(g_rhs_zero_im)
            if(.not.file_read_ga(cphf_sol_im(1),g_rhs_zero_im))
     &         call errquit
     $         ('aoresponse: could not read cphf zero solution Im',
     &         0, DISK_ERR)          
            call nga_copy_patch('n',g_rhs_zero_im,blo,bhi,
     &         g_amat_zero_im,alo,ahi)
          endif                 ! lifetime
          
          write(cstemp,'(a)') 'cphf_rhs_zero'
          if (.not.ga_destroy(g_rhs_zero)) call 
     &       errquit('aoresponse: ga_destroy failed '//cstemp(1:13),
     &       0,GA_ERR)

          if (lifetime) then
            write(cstemp,'(a)') 'cphf_rhs_zero_im'
            if (.not.ga_destroy(g_rhs_zero_im)) call 
     &         errquit('aoresponse: ga_destroy failed '//cstemp(1:15),
     &         0,GA_ERR)
          endif

          call util_file_unlink(cphf_rhs(1))
          call util_file_unlink(cphf_sol(1))
          if (lifetime) then
            call util_file_unlink(cphf_rhs_im(1))
            call util_file_unlink(cphf_sol_im(1))
          endif
          
        endif                   ! lzero
        
c$$$        call int_init(rtdb,1,basis)
c$$$        call schwarz_init(geom,basis)
c$$$        call hnd_giao_init(basis,1)
        
        call ga_sync()          ! when do we actually need that ???
        if (debug) write (6,*) 'CPKS complete, results in g_amat'

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

        dims(1) = nvir
        if (haveocc) dims(1) = nmo ! need full range of amat
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
        if (lzero) then
          write(cstemp,'(a)') 'g_rhs_zero'
          if (debug) write (luout,*) cstemp(1:10)
          if (.not.nga_create(MT_DBL,2,dims,cstemp(1:10),chunk,
     &       g_rhs_zero)) call 
     &       errquit('aoresponse: nga_create failed'//cstemp(1:10),
     &       0,GA_ERR)
          call ga_zero (g_rhs_zero)          
          if (lifetime) then
            write(cstemp,'(a)') 'g_rhs_zero_im'
            if (debug) write (luout,*) cstemp(1:13)
            if (.not.nga_create(MT_DBL,2,dims,cstemp(1:13),chunk,
     &         g_rhs_zero_im)) call 
     &         errquit('aoresponse: nga_create failed'//cstemp(1:13),
     &         0,GA_ERR)
            call ga_zero (g_rhs_zero_im)
          endif                 ! lifetime
        endif                   ! lzero

c       ------------------------------------------------------
c       memory allocation complete, now proceed to calculating
c       the perturbed density matrices
c       ------------------------------------------------------
        
        do idir = 1,3           ! direction of the perturbing field
          
          alo(1) = nocc + 1
          if (haveocc) alo(1) = 1 ! need full range of amat
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
          if (lzero) then
            call nga_copy_patch ('n',g_amat_zero, alo, ahi,
     &         g_rhs_zero, blo, bhi)
            if (lifetime)
     &         call nga_copy_patch ('n',g_amat_zero_im, alo, ahi,
     &         g_rhs_zero_im, blo, bhi)
          endif                 ! lzero
                    
          if (dbgmat) then
            do ipm = 1,ncomp
              call ga_print(g_rhs(ipm))
            enddo
          endif

c         calculate density matrices for perturbing field component:

          if (debug) write (6,*) 'calculating density matrices'
          lantisym = .false.    ! do not symmetrize density matrices
          if (lstatic) then 
            call CalcPerturbedTDPmat1(ncomp, g_pmats(1), g_pmata(1),
     &         g_rhs(1),
     &         g_vectors(ispin),
     &         nbf, nocct(ispin),
     &         nvirt(ispin), nmot(ispin), lantisym, lstatic, limag,
     &         haveocc)
            if (lifetime)
     &         call CalcPerturbedTDPmat1(ncomp,
     &         g_pmats_im(1), g_pmata_im(1),
     &         g_rhs_im(1),
     &         g_vectors(ispin),
     &         nbf, nocct(ispin),
     &         nvirt(ispin), nmot(ispin), lantisym, lstatic, limag,
     &         haveocc)
          else
c           calculate both components
c           TEST LATER IF ONE IS SUFFICIENT. WE NEED ONLY ONE LATER
            call CalcPerturbedTDPmat1(ncomp, g_pmats, g_pmata, g_rhs,
     &         g_vectors (ispin),
     &         nbf, nocct(ispin),
     &         nvirt(ispin), nmot(ispin),  lantisym, lstatic,
     &         limag, haveocc)
            if (lifetime)
     &         call CalcPerturbedTDPmat1(ncomp, g_pmats_im,
     &         g_pmata_im, g_rhs_im,
     &         g_vectors (ispin),
     &         nbf, nocct(ispin),
     &         nvirt(ispin), nmot(ispin),  lantisym, lstatic,
     &         limag, haveocc)
c
            if (lzero) then     ! evaluate w=0 response (later subtracted)
              call CalcPerturbedTDPmat1(1, g_pmat_zero, g_pmata(1),
     &           g_rhs_zero,
     &           g_vectors(ispin),
     &           nbf, nocct(ispin),
     &           nvirt(ispin), nmot(ispin), lantisym, .true., limag,
     &           haveocc)
              if (lifetime)
     &           call CalcPerturbedTDPmat1(1, g_pmat_zero_im,
     &           g_pmata_im(1),
     &           g_rhs_zero_im,
     &           g_vectors(ispin),
     &           nbf, nocct(ispin),
     &           nvirt(ispin), nmot(ispin), lantisym, .true., limag,
     &           haveocc)
            endif               ! lzero
          end if                ! lstatic
          
          if (debug) write (6,*) 'density matrices calculated'
c         
c         P(+) is now in pmats(1), P(-) in pmats(2)
c         
          if (dbgmat) then
            call ga_print (g_pmats(1))
            if (.not.lstatic) call ga_print (g_pmats(2))
          end if
          
c         if (debug) then
c         c         calculate tr(P*S) with the perturbed density matrix
c         do ipm = 1,2
c         call ga_add(1d0, g_pmats(ipm), 1d0, g_pmata(ipm),
c         &         g_pmunu)
c         c           pmunu now holds the full P(+) or (P-)
c         call ga_matmul_patch('n','n', 2d0,0d0,
c         &         g_pmunu,1,nbf,1,nbf,
c         &         g_smat0,1,nbf,1,nbf,
c         &         g_work,1,nbf,1,nbf)
c         sum = ga_trace_diag(g_work)
c         if (ga_nodeid().eq.0)
c         &         write (luout,'(1x,a,i3,e15.7)') 'tr(P*S): ',ipm,sum
c         enddo                 ! ipm = 1,2
c         endif                   ! debug         
          
c         we use ipm = 1 _or_ ipm = 2 here; yields the same result
          
          ipm = 1
          
c         --------------------------------------------------
c         start loop over components of responding vector
c         and calculate tr(H*P), where H are the integrals
c         for the response property (e.g. dipole integrals)
c         and P is the perturbed density matrix in g_pmats/a
c         --------------------------------------------------
          
          alo(1) = 1
          ahi(1) = nbf
          alo(2) = 1
          ahi(2) = nbf
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nbf
          
          do iresp = 1,3        ! component of the response vector
            
            blo(3) = iresp      ! direction of response
            bhi(3) = iresp      ! for arrays g_dipel, g_magel
            
c           -----------------------------------------
c           calculate perturbation of electric dipole
c           -----------------------------------------

c           assign component of "Re(alfa)" matrix (output)

            call ga_zero(g_work)
            call nga_matmul_patch('n','n',1d0,0d0,
     &         g_pmats(ipm),alo,ahi,
     &         g_dipel,blo,bhi,
     &         g_work,alo,ahi)
            
            sum = -1d0 * ga_trace_diag(g_work)
            if (debug .and. oprint)
     &         write (luout,'(1x,a,i3,e15.7)') 'tr(D*P1): ',iresp,sum 
            
            alfare(idir, iresp) = sum

            if (lifetime) then

c             assign component of "Im(alfa)" matrix (output)
c             (note: sign reversed)

              call ga_zero(g_work)
              call nga_matmul_patch('n','n',1d0,0d0,
     &           g_pmats_im(ipm),alo,ahi,
     &           g_dipel,blo,bhi,
     &           g_work,alo,ahi)
              
c              sum = 1d0 * ga_trace_diag(g_work)   
              sum = -1d0 * ga_trace_diag(g_work)   ! FA-03-12-14               
              
              alfaim(idir, iresp) = sum
            else
              alfaim(idir, iresp) = 0d0
            endif               ! lifetime
            
c           determine zero-freq. response for velocity formalism
c           in a similar fasion for Re(alpha) and Im(alpha):

            if (lzero) then
              call ga_zero(g_work)
              call nga_matmul_patch('n','n',1d0,0d0,
     &           g_pmat_zero,alo,ahi,
     &           g_dipel,blo,bhi,
     &           g_work,alo,ahi)            
              sum = -1d0 * ga_trace_diag(g_work)
              alfa0re(idir, iresp) = sum

              if (lifetime) then
                call ga_zero(g_work)
                call nga_matmul_patch('n','n',1d0,0d0,
     &             g_pmat_zero_im,alo,ahi,
     &             g_dipel,blo,bhi,
     &             g_work,alo,ahi)            
                sum = 1d0 * ga_trace_diag(g_work)
                alfa0im(idir, iresp) = sum
              else
                alfa0im(idir, iresp) = 0d0
              endif             ! lifetime
            endif               ! lzero
            
c           -----------------------------------------
c           calculate perturbation of magnetic dipole
c           -----------------------------------------

c           assign component of "Re(beta)" matrix (output)
            
            call ga_zero(g_work)
            call nga_matmul_patch('n','n',1d0,0d0,
     &         g_pmats(ipm),alo,ahi,
     &         g_dipmag,blo,bhi,
     &         g_work,alo,ahi)
            
            if (dbgmat) call nga_print_patch(g_dipmag,blo,bhi,1)
            
            sum = -1d0 * ga_trace_diag(g_work)
            if (debug .and. oprint)
     &         write (luout,'(1x,a,i3,e15.7)') 'tr(M*P1): ',iresp,sum 
            
            betare(idir, iresp) = sum

            if (lifetime) then
            
c             assign component of "Im(beta)" matrix (output):
c             (note: sign reversed)

              call ga_zero(g_work)
              call nga_matmul_patch('n','n',1d0,0d0,
     &           g_pmats_im(ipm),alo,ahi,
     &           g_dipmag,blo,bhi,
     &           g_work,alo,ahi)
              
              if (dbgmat) call nga_print_patch(g_dipmag,blo,bhi,1)
              
              sum =  1d0 * ga_trace_diag(g_work)    
              betaim(idir, iresp) = sum
            else
              betaim(idir, iresp) = 0d0
            endif               ! lifetime

c           determine zero-freq. response for velocity formalism
c           for beta matrix:

            if (lzero) then
c             real part:
              call ga_zero(g_work)
              call nga_matmul_patch('n','n',1d0,0d0,
     &           g_pmat_zero,alo,ahi,
     &           g_dipmag,blo,bhi,
     &           g_work,alo,ahi)            
              sum = -1d0 * ga_trace_diag(g_work)
              beta0re(idir, iresp) = sum

              if (lifetime) then
c               imaginary part:
                call ga_zero(g_work)
                call nga_matmul_patch('n','n',1d0,0d0,
     &             g_pmat_zero_im,alo,ahi,
     &             g_dipmag,blo,bhi,
     &             g_work,alo,ahi)            
                sum =  1d0 * ga_trace_diag(g_work)
                beta0im(idir, iresp) = sum
              else
                beta0im(idir, iresp) = 0d0
              endif             ! lifetime
            endif               ! lzero
            
          end do                ! iresp = 1,3

        end do                  ! idir = 1,3

       if (debug) then
          write (6,*) 'alfare: ',alfare
          write (6,*) 'betare: ',betare
          write (6,*) 'alfaim: ',alfaim
          write (6,*) 'betaim: ',betaim
          write (6,*) 'alfa0re: ',alfa0re
          write (6,*) 'beta0re: ',beta0re
          write (6,*) 'alfa0im: ',alfa0im
          write (6,*) 'beta0im: ',beta0im
       endif
        
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
        
c       ---------------
c       deallocate amat
c       ---------------
        
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

c       ---------------------------------------------------------
c       subtract static response form dynamic in case of velocity
c       formalism
c       ---------------------------------------------------------
        
        if (lvelocity .and. .not.lstatic) then
          do idir = 1,3
            do iresp = 1,3
              alfare(idir,iresp) =alfare(idir,iresp)-alfa0re(idir,iresp)
              betare(idir,iresp) =betare(idir,iresp)-beta0re(idir,iresp)
              alfaim(idir,iresp) =alfaim(idir,iresp)-alfa0im(idir,iresp)
              betaim(idir,iresp) =betaim(idir,iresp)-beta0im(idir,iresp)
            enddo
          enddo
        endif                   ! lvelocity
        
c       for magnetic perturbation calculate beta tensor
c       (only for nonzero frequency, otherwise calculate G').
c       also consider the case of velocity formalism
        
        do idir = 1,3
          do iresp = 1,3
            if (lmagpert .and. .not.lstatic) then
c             case I: magnetic perturbation. alpha = G', beta = Chi-p
              if (.not.lgprime) then
                alfare(idir,iresp) = alfare(idir,iresp) / omega
                alfaim(idir,iresp) = alfaim(idir,iresp) / omega
              end if
              if (lvelocity) then
                alfare(idir,iresp) = alfare(idir,iresp) / omega
                alfaim(idir,iresp) = alfaim(idir,iresp) / omega
              endif
              
            else if (.not.lmagpert .and. .not.lstatic) then
              
c             case II: alpha = polarizability, beta = G'
              if (.not.lgprime .or. lvelocity) then
                scaling = one/omega
                if (lvelocity.and..not.lgprime) scaling = one/(omega**2)
                betare(idir,iresp) = betare(idir,iresp) * scaling
                betaim(idir,iresp) = betaim(idir,iresp) * scaling
              endif
              if (lvelocity) then
                alfare(idir,iresp) = alfare(idir,iresp) / omega**2
                alfaim(idir,iresp) = alfaim(idir,iresp) / omega**2
              endif
            endif ! lmagpert or not
          enddo                 ! iresp
        enddo                   ! idir
        
c       print warning in case omega=0 that we haven't divided by omega
c       where it would have been necessary:
        
        if (lstatic .and. lvelocity) then
          if (oprint) write (luout,*)
     +       'NOTE: Because of omega=0 the response functions involving'
          if (oprint) write (luout,*)
     +       '      the dipole-velocity operator have NOT been '//
     +       ' divided by omega'
        endif
        if (lstatic) then
          if (oprint) write (luout,*)
     +       'NOTE: because of omega=0 the electric-magnetic response'
          if (oprint) write (luout,*)
     +       '      is incorrect'
        endif
        
      enddo                     ! ispin = 1,2 from way above
              
c     ---------------------------------------------------------------
c     end loop over spin components (which we don't use right now
c     since nspin is forced to be 1 at the beginning of this routine)
c     ---------------------------------------------------------------
      
c     ------------------------------------------------
c     deallocate remaining allocated memory and return
c     ------------------------------------------------
            
      do ipm = 1,ncomp
        if (.not.ga_destroy(g_pmats(ipm)))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_pmats',
     &     0,GA_ERR)

        if (.not.ga_destroy(g_pmata(ipm)))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_pmata',
     &     0,GA_ERR)

        if (lifetime) then
          if (.not.ga_destroy(g_pmats_im(ipm)))
     &       call 
     &       errquit('aoresponse: ga_destroy failed g_pmats_im',
     &       0,GA_ERR)    
      
          if (.not.ga_destroy(g_pmata_im(ipm)))
     &       call 
     &       errquit('aoresponse: ga_destroy failed g_pmata_im',
     &       0,GA_ERR)
        endif                   ! lifetime
      enddo                     ! ipm = 1,ncomp

      if (.not.ga_destroy(g_work))
     &   call 
     &   errquit('aoresponse: ga_destroy failed g_work',
     &   0,GA_ERR)

      if (lzero) then
        if (.not.ga_destroy(g_pmat_zero))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_pmat_zero',
     &     0,GA_ERR)

        if (.not.ga_destroy(g_rhs_zero))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_rhs_zero',
     &     0,GA_ERR)

        if (.not.ga_destroy(g_amat_zero))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_amat_zero',
     &     0,GA_ERR)

        if (lifetime) then
        if (.not.ga_destroy(g_pmat_zero_im))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_pmat_zero_im',
     &     0,GA_ERR)

        if (.not.ga_destroy(g_rhs_zero_im))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_rhs_zero_im',
     &     0,GA_ERR)

        if (.not.ga_destroy(g_amat_zero_im))
     &     call 
     &     errquit('aoresponse: ga_destroy failed g_amat_zero_im',
     &     0,GA_ERR)
        endif                   ! lifetime

      endif                     ! lzero

c     it seems that if we use GIAOs everything is off by a factor of
c     two, so we need to scale alfare, alfaim. 

      if (lgiao .and. lmagpert) then
        scaling = half
        do idir = 1,3
          do iresp = 1,3
            alfare(idir, iresp) = alfare(idir, iresp) * scaling
            alfaim(idir, iresp) = alfaim(idir, iresp) * scaling
          end do
        end do
      end if                    ! lstatic

      
c     ==================================================================
      
      return
      
      end

