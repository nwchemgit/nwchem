      subroutine CalcPerturbedTDPmat1
     &   (ncomp, g_pmats, g_pmata, g_amat, g_vectors, naos, nocc,
     &   nvir, nmo, lantisym, lstatic, imag, haveocc)

* $Id$

c     ==================================================================

c     calculate frequency-dependent density matrix perturbation
c     (symmetric and antisymmetric part), linear response,
c     from a set of perturbed MO coefficients

c     we assume DOUBLE occupation of all occupied orbitals and
c     REAL unperturbed orbitals. The perturbation can be either
c     purely real or purely imaginary

c     THIS ROUTINE USES TOO MUCH MEMORY; IT COULD DO HE SAME
C     JOB WITH LESS TEMP SPACE. FIX THIS

c     input: 
c     ncomp      - number of components to be calculated
c     g_amat     - the perturbed MO coefficients are
c                  written as C(+/-) = C(0)*A(+/-),
c                  g_amat contains the elements of matrix A
c                  (only the virt - occ block, or nmo - occ block)
c     g_vectors  - unperturbed MO coefficients C(0)
c     lantisym   - logical switch to calculate symmetric
c                  and antisymmetric
c                  part separately or just the total density matrix
c     lstatic    - static response, assume that both components
c                  of amat are equal. assumes ncomp = 1 (!)
c     imag      - true if amat is imaginary instad of real
c     haveocc    - true if amat contains occ-occ block, too
      
c     output : g_pmats, g_pmata: symmetric and antisymmetric
c     part of perturbed density matrix, global arrays, if (lantisym),
c     otherwise the total density matrix is in pmats, and pmata=0
      
c     remark : all perturbed matrices are classified by
c     (+/-) frequency components 
c     
c     remark: the density matrix is given by
c     transpose(P) = C n C(dagger), i.e. in the end we transpose the
c     result to get the correct density matrix out

c     ==================================================================
      
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"

c     subroutine arguments:
      integer ncomp
      integer g_pmats(ncomp), g_pmata(ncomp), g_amat(*),
     &   g_vectors           ! GA
      integer naos, nocc, nvir, nmo
      logical lantisym
      logical lstatic, imag, haveocc

c     local variables:
      integer g_ptmp, g_eig1, g_work
      integer imo, imu, inu, ll, ip, im, ipm
      double precision rtemp
      double precision half, one, two
      parameter (half = 0.5d0, one = 1.0d0)
      logical debug

      integer type, dim1, dim2


c     ==================================================================

      debug = .false.

c     check range of ncomp

      if (ncomp.le.0 .or. ncomp.gt.2) then
        call errquit('CalcPerturbedTDPmat: ncomp out of range',
     &     0,CALC_ERR)
      endif

c     cowardy refuse so calculate two components of perturbed
c     density matrix if lstatic switch is set to true
 
      if (ncomp.gt.1 .and. lstatic) then
        call errquit
     &     ('CalcPerturbedTDPmat1: ncomp > 1 but lstatic.eq.true.',
     &     0,CALC_ERR)
      endif


c     assign + and - components for indexing amat:

      if (lstatic) then
        ip = 1
        im = 1
      else
        ip = 1
        im = 2
      endif

c     if (ga_nodeid().eq.0) then
c      write(*,1) ip,im,lstatic
c1     format('(ip,im,lstatic)=(',i3,',',i3,',',L1,')')
c     endif
        
      do ipm = 1,ncomp
        call ga_zero(g_pmats(ipm))
        call ga_zero(g_pmata(ipm))
      enddo

      if (debug) write (6,'(a,4i6)') 'nocc,nvir,nmo',nocc, nvir, nmo

      if (nocc+nvir .ne. nmo) call errquit
     &   ('CalcPerturbedTDPmat1: wrong no. of orbitals',0,CALC_ERR)

c     ------------------------
c     allocate workspace (GAs)
c     ------------------------
      
      if (.not. ga_create(MT_DBL, naos, naos,
     &   'CalcPerturbedTDPmat1:ptmp',
     &   0, 0, g_ptmp)) call errquit('CalcPerturbedTDPmat1:ptmp', 0,
     &   GA_ERR)

      if (.not. ga_create(MT_DBL, naos, naos,
     &   'CalcPerturbedTDPmat1:work',
     &   0, 0, g_work)) call errquit('CalcPerturbedTDPmat1:work', 0,
     &   GA_ERR)
      
      if (.not. ga_create(MT_DBL, naos, nocc,
     &   'CalcPerturbedTDPmat1:eig1',
     &   0, 0, g_eig1)) call errquit('CalcPerturbedTDPmat1:eig1', 0,
     &   GA_ERR)

      if (debug) then
c       debug array dimensions
        call ga_inquire (g_eig1,type, dim1, dim2)
        write (6,'(a,2i4)') 'g_eig1:',dim1,dim2
        call ga_inquire (g_ptmp,type, dim1, dim2)
        write (6,'(a,2i4)') 'g_ptmp:',dim1,dim2
        call ga_inquire (g_work,type, dim1, dim2)
        write (6,'(a,2i4)') 'g_work:',dim1,dim2
        call ga_inquire (g_amat(1),type, dim1, dim2)
        write (6,'(a,2i4)') 'g_amat(1):',dim1,dim2
        call ga_inquire (g_vectors,type, dim1, dim2)
        write (6,'(a,2i4)') 'g_vectors:',dim1,dim2
      endif


c     -------------------------------------------------------------
c     First we assemble P(+). Note that A(-) is assumed to be A(-)*
c     in fact (A = amat)
c     This allows us to use the same algorithm no matter if A is
c     real and symmetric or imaginary and antisymmetric
c     -------------------------------------------------------------

      call ga_zero(g_ptmp)
      call ga_zero(g_work)
      call ga_zero(g_eig1)
      call ga_sync()

c     ----------------------------
c     First step:  C n C(-,dagger)
c     ----------------------------

c     calculate C(-,dagger)
      two = 2d0
      if (imag) two = -2d0

      if (.not.haveocc) then
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-enter-no-haveocc...'
       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,nocc+1,nmo,
     &   g_amat(im),1,nvir,1     ,nocc,
     &   g_eig1    ,1,naos,1     ,nocc)

c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-----g_eig1-nohaveocc---- START'
c          call ga_print(g_eig1)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-----g_eig1-nohaveocc---- END'

      else
c         if (ga_nodeid().eq.0)
c     &    write(*,*) 'FA-enter-yes-haveocc...'

       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,1,nmo,
     &   g_amat(im),1,nmo ,1,nocc,
     &   g_eig1    ,1,naos,1,nocc)

c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-----g_eig1-yeshaveocc---- START'
c         call ga_print(g_eig1)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-----g_eig1-yeshaveocc---- END'
      endif ! .not.haveocc
      call ga_sync()

      if (debug) write (6,*) '1'

c     note: the dimensioning for array B is that of the transposed
c     matrix, not of the original matrix. 
      
c     calculate C(0)C(-,dagger), store in g_ptmp
      call ga_matmul_patch('n','t', 1d0,0d0,
     &   g_vectors,1,naos,1,nocc,
     &   g_eig1,1,nocc,1,naos,
     &   g_ptmp,1,naos,1,naos)  
      call ga_sync()

      if (debug) write (6,*) '2'   
      
c     -----------------------------
c     Second step: C(+) n C(dagger)
c     -----------------------------
      
c     calculate C(+)
      two = 2.0d0
      if (.not.haveocc) then

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-enter-no-haveocc-2...'

       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,nocc+1,nmo,
     &   g_amat(ip),1,nvir,1     ,nocc,
     &   g_eig1    ,1,naos,1     ,nocc)

c       if (ga_nodeid().eq.0)
c     &  write(*,*) '-----g_eig1-nohaveocc-2-- START'
c       call ga_print(g_eig1)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '-----g_eig1-nohaveocc-2-- END'

      else

c         if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-enter-yes-haveocc-2...'

       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,1,nmo,
     &   g_amat(ip),1,nmo ,1,nocc,
     &   g_eig1    ,1,naos,1,nocc)

c        if (ga_nodeid().eq.0)
c     &   write(*,*) '-----g_eig1-yeshaveocc-2---- START'
c        call ga_print(g_eig1)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '-----g_eig1-yeshaveocc-2---- END'

      endif
      call ga_sync()

      if (debug) write (6,*) '3'
      
c     calculate C(+)C(0,dagger), store in g_work
      call ga_matmul_patch('n','t', 1d0,0d0,
     &   g_eig1   ,1,naos,1,nocc,
     &   g_vectors,1,nocc,1,naos,
     &   g_work   ,1,naos,1,naos)
      call ga_sync()

      if (debug) write (6,*) '4'

c     add the two terms together and transpose the density matrix
c             if (ga_nodeid().eq.0)
c     &       write(*,*) '---- CC1^t-------- START'
c             call ga_print(g_ptmp)
c            if (ga_nodeid().eq.0)
c     &       write(*,*) '---- CC1^t--------- END'           
c             if (ga_nodeid().eq.0)
c     &       write(*,*) '---- C1C^t-------- START'
c             call ga_print(g_work)
c            if (ga_nodeid().eq.0)
c     &       write(*,*) '---- C1C^t--------- END'           

      call ga_add(1d0,g_ptmp,1d0,g_work,g_work)

c             if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmats-0-------- START'
c             call ga_print(g_work)
c            if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmats-0-------- END'           

      call ga_sync()
      call ga_transpose(g_work, g_ptmp)
      call ga_sync()

c     calculate symmetrized and antisymmetrized part (+ component)
c     if requested on input:

      if (lantisym) then
        call ga_transpose(g_ptmp,g_work)
        call ga_sync()
        call ga_add(half,g_ptmp,half,g_work,g_pmats(1))
        call ga_sync()
        call ga_add(half,g_ptmp,-half,g_work,g_pmata(1))
        call ga_sync()
      else
        call ga_copy(g_ptmp, g_pmats(1))
        call ga_sync()
      endif
c             if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmats-------- START'
c             call ga_print(g_pmats(1))
c            if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmats-------- END'           

      if (debug) write (6,*) '4a'

      if (lstatic .or. ncomp.eq.1) then
c       skip calculation of component 2 of the density matrix
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-Skipping calc of 2nd component'

        goto 7000
      endif  
      
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-Doing calc of 2nd component'     
 
c     -----------------------------------------
c     Next step: assemble P(-). Same as before,
c     but +/- interchanged in amat:
c     -----------------------------------------

      call ga_zero(g_ptmp)
      call ga_zero(g_work)
      call ga_zero(g_eig1)
      call ga_sync()

c     ----------------------------
c     First step:  C n C(+,dagger)
c     ----------------------------

c     calculate C(+,dagger)
      two = 2d0
      if (imag) two = -2d0
      if (.not.haveocc) then
      call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors,1,naos,nocc+1,nmo,
     &   g_amat(ip),1,nvir,1,nocc,
     &   g_eig1,1,naos,1,nocc)
      else
      call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors,1,naos,1,nmo,
     &   g_amat(ip),1,nmo,1,nocc,
     &   g_eig1,1,naos,1,nocc)
      endif
      call ga_sync()

      if (debug) write (6,*) '5'

c     calculate C(0)C(+,dagger), store in g_ptmp
      call ga_matmul_patch('n','t', 1d0,0d0,
     &   g_vectors,1,naos,1,nocc,
     &   g_eig1,1,nocc,1,naos,
     &   g_ptmp,1,naos,1,naos)
      call ga_sync()

      if (debug) write (6,*) '6'
      
c     -----------------------------
c     Second step: C(-) n C(dagger)
c     -----------------------------

c     calculate C(-)
      two = 2d0
      if (.not.haveocc) then
      call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,nocc+1,nmo,
     &   g_amat(im),1,nvir,1     ,nocc,
     &   g_eig1    ,1,naos,1     ,nocc)
      else
      call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,1,nmo,
     &   g_amat(im),1,nmo ,1,nocc,
     &   g_eig1    ,1,naos,1,nocc)
      endif
      call ga_sync()

      if (debug) write (6,*) '7'

c     calculate C(-)C(0,dagger), store in g_work
      call ga_matmul_patch('n','t', 1d0,0d0,
     &   g_eig1   ,1,naos,1,nocc,
     &   g_vectors,1,nocc,1,naos,
     &   g_work   ,1,naos,1,naos)
      call ga_sync()

      if (debug) write (6,*) '8'

c     add the two terms together and transpose
      call ga_add(1d0,g_ptmp,1d0,g_work,g_work)
      call ga_sync()
      call ga_transpose(g_work, g_ptmp)
      call ga_sync()

c     calculate symmetrized and antisymmetrized part (- component)

      if (lantisym) then
        call ga_transpose(g_ptmp,g_work)
        call ga_sync()
        call ga_add(half,g_ptmp,half,g_work,g_pmats(2))
        call ga_sync()
        call ga_add(half,g_ptmp,-half,g_work,g_pmata(2))
        call ga_sync()
      else
        call ga_copy(g_ptmp, g_pmats(2))
        call ga_sync()
      endif

      if (debug) write (6,*) '8a'
      
c     ---------------------------------------------
c     deallocate temporary arrays, sync, and return
c     ---------------------------------------------

c     jump here from above in case of static calculation
 7000 continue
      
      if (.not.ga_destroy(g_ptmp))
     &   call 
     &   errquit('CalcPerturbedTDPmat: ga_destroy failed g_ptmp',
     &   0,GA_ERR)
      
      if (.not.ga_destroy(g_work))
     &   call 
     &   errquit('CalcPerturbedTDPmat: ga_destroy failed g_work',
     &   0,GA_ERR)
      
      if (.not.ga_destroy(g_eig1))
     &   call 
     &   errquit('CalcPerturbedTDPmat: ga_destroy failed g_eig1',
     &   0,GA_ERR)

      call ga_sync()

c     ==================================================================
      
      return
      end
