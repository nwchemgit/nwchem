      subroutine CalcPerturbedTDPmat1_opt(
     &                 ncomp,    ! in : nr. components to calculate
     &                 g_pmats,  ! out: density matrix      symmetrized
     &                 g_pmata,  ! out: density matrix  antisymmetrized -NOT USED
     &                 g_amat,   ! in :   perturbed MO coefficients 
     &                 g_vectors,! in : unperturbed MO coefficients
     &                 naos,     ! in : nr. AOs
     &                 nocc,     ! in : nr. occupied MOs
     &                 nvir,     ! in : nr. virtual  MOs
     &                 nmo,      ! in : nr. MOs
     &                 lantisym, ! in : = .true. calc. (symm,antisymm)=(pmats,pmata)
     &                 lstatic,  ! in : = .true. static response, dynamic otherwise
     &                 imag,     ! in : = .true. if amat is imaginary instead of real
     &                 haveocc)  ! in : = .true. if amat contains occ-occ block

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
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c        Modifying original source code: CalcPerturbedTDPmat1()
    
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c     subroutine arguments:
      integer ncomp
      integer g_pmats(ncomp), 
     &        g_pmata(ncomp), 
     &        g_amat(ncomp),
     &        g_vectors          
      integer naos, nocc, nvir, nmo
      logical lantisym
      logical lstatic, imag, haveocc

c     local variables:
      integer ip, im, ipm
      double precision half, one, two,toscl
      parameter (half = 0.5d0, one = 1.0d0)
      logical debug
      external get_PertDens

c      debug = .true. !     allow debugging printout
      debug = .false. ! not allow debugging printout
c     check range of ncomp
      if (ncomp.le.0 .or. ncomp.gt.2) then
        call errquit('CalcPerturbedTDPmat: ncomp out of range',
     &     0,CALC_ERR)
      endif
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

      call ga_sync()

c      if (ga_nodeid().eq.0) then
c       write(*,190) ip,im,lstatic,imag
c 190     format('(ip,im,lstatic,imag)=(',i3,',',i3,',',L1,',',L1,')')
c      endif
        
      do ipm = 1,ncomp
        call ga_zero(g_pmats(ipm))
c        call ga_zero(g_pmata(ipm))
      enddo

      if (nocc+nvir .ne. nmo) call errquit
     &   ('CalcPerturbedTDPmat1: wrong no. of orbitals',0,CALC_ERR)

c     -------------------------------------------------------------
c     First we assemble P(+). Note that A(-) is assumed to be A(-)*
c     in fact (A = amat)
c     This allows us to use the same algorithm no matter if A is
c     real and symmetric or imaginary and antisymmetric
c     -------------------------------------------------------------
c     ----------------------------
c     First  step: C    n C(-,dagger)
c     Second step: C(+) n C(dagger)
c     -----------------------------
      if (debug) then
       do ipm=1,2
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- BEF get_PertDens g_amat(',ipm,')--- START'
         call ga_print(g_amat(ipm))
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- BEF get_PertDens g_amat(',ipm,')--- END'
       enddo ! end-loop-ipm 
      endif ! end-if-debug     

      call get_PertDens(
     &               g_pmats(1),!out:      symmetrized dens
     &               g_pmata(1),!out: anti symmetrized dens
     &               g_amat,    ! in: the u mat
     &               g_vectors, ! in: MO vect
     &               ncomp,     ! in: nr. components g_amat
     &               ip,im,     ! in: indices of u vect
     &               imag,      ! in: = T -> imag
     &               haveocc,   ! in: logical flag
     &               lantisym,  ! in: logical flag
     &               naos,nmo,  ! in: nr. AOs,MOs
     &               nocc,nvir, ! in: nr. (occ,virt) MOs
     &               debug)     ! in: =.true. -> show debug printouts

        if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-nw-------- START'
             call ga_print(g_pmats(1))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-nw-------- END'      
c             if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmata-nw-------- START'
c             call ga_print(g_pmata(1))
c            if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmata-nw-------- END'   
        endif ! end-if-debug     
      if (lstatic .or. ncomp.eq.1) then
c       skip calculation of component 2 of the density matrix
        if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-Skipping calc of 2nd component'

        goto 7000
      endif  
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-calculating 2nd component'
c     ----------------------------
c     First  step: C    n C(+,dagger)
c     Second step: C(-) n C(dagger)
c     -----------------------------
c     Note 1.- I tried calling get_PertDens 2nd time
c            and swapping (ip,im) --> (im,ip)
c            and I got the relationship shown below:
c     lantisym=.true. :
c     g_pmats(2)= -transpose(g_pmats(1))
c     g_pmata(2)= g_pmata(1)= 0
c     lantisym=.false. :
c     g_pmats(2)= -transpose(g_pmats(1))
c     g_pmata(2)= g_pmata(1)
c     Note 2.- Swapping (ip,im) --> (im,ip)
c              is equivalent to doing the second
c              part (e.g. calc. of P(-)
      call ga_transpose(g_pmats(1),g_pmats(2))
                toscl= 1.0d0
      if (imag) toscl=-1.0d0
c      call ga_scale(g_pmats(2),-1.0d0)
      call ga_scale(g_pmats(2),toscl)
c      call ga_copy(g_pmata(1),g_pmata(2))
c     jump here from above in case of static calculation
 7000 continue    
      return
      end

      subroutine get_PertDens(
     &               g_pmats,  !out:      symmetrized dens
     &               g_pmata,  !out: anti symmetrized dens -NOT USED
     &               g_amat,   ! in: the u mat
     &               g_vectors,! in: MO vect
     &               ncomp,    ! in: nr. components g_amat
     &               ip,im,    ! in: indices of u vect
     &               imag,     ! in: = T -> imag
     &               haveocc,  ! in: logical flag
     &               lantisym, ! in: logical flag
     &               naos,nmo, ! in: nr. AOs,MOs
     &               nocc,nvir,! in: nr. (occ,virt) MOs
     &               debug)    ! in: =.true. -> show debug printouts
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"

c     subroutine arguments:
c Note.- For (im,ip) --> (g_pmats,g_pmata)(1)
c        For (ip,im) --> (g_pmats,g_pmata)(2)
      integer ncomp,ip,im
      integer g_pmats,g_pmata,
     &        g_amat(ncomp),
     &        g_vectors           
      integer naos,nocc,nvir,nmo
      logical lantisym,imag,haveocc
c     local variables:
      integer g_ptmp,g_eig1,g_work
      double precision half, one, two
      parameter (half = 0.5d0, one = 1.0d0)
      logical debug
c     ------------------------
c     allocate workspace (GAs)
c     ------------------------    
      if (.not. ga_create(
     &    MT_DBL,naos,naos,'get_PertDens:ptmp',0,0,g_ptmp)) 
     &    call errquit('get_PertDens:ptmp',0,GA_ERR)
      if (.not. ga_create(
     &    MT_DBL,naos,naos,'get_PertDens:work',0,0,g_work)) 
     &    call errquit('get_PertDens:work',0,GA_ERR)
      if (.not. ga_create(
     &    MT_DBL,naos,naos,'get_PertDens:eig1',0,0,g_eig1)) 
     &    call errquit('get_PertDens:eig1',0,GA_ERR)
c     ----------------------------
c     First step:  C n C(-,dagger)
c     ----------------------------
c     calculate C(-,dagger)
      two = 2d0
      if (imag) two = -2d0
      if (.not.haveocc) then

       if (debug) then
         if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-enter-no-haveocc...'
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-nohaveocc---- START'
          call ga_print(g_vectors)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-nohaveocc---- END'
         if (ga_nodeid().eq.0) then
           write(*,1) im
 1         format('-----g_amat(',i3,')-nohaveocc---- START')
         endif
         call ga_print(g_amat(im))
         if (ga_nodeid().eq.0) then
           write(*,2) im
 2         format('-----g_amat(',i3,')-nohaveocc---- END')
         endif
          if (ga_nodeid().eq.0) then
           write(*,35) ip
 35         format('-----g_amat(',i3,')-yeshaveocc---- START')
          endif
          call ga_print(g_amat(ip))
          if (ga_nodeid().eq.0) then
           write(*,36) ip
 36         format('-----g_amat(',i3,')-yeshaveocc---- END')
          endif

       endif ! end-if-debug

       call ga_zero(g_eig1)
       call ga_matmul_patch('n','n', two,0d0,
     &                      g_vectors ,1,naos,nocc+1,nmo,
     &                      g_amat(im),1,nvir,1     ,nocc,
     &                      g_eig1    ,1,naos,1     ,nocc)
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_eig1-nohaveocc---- START'
          call ga_print(g_eig1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_eig1-nohaveocc---- END'
        endif ! end-if-debug
      else
        if (debug) then
          if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-yeshaveocc---- START'
          call ga_print(g_vectors)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-yeshaveocc---- END'
          if (ga_nodeid().eq.0) then
           write(*,3) im
 3         format('-----g_amat(',i3,')-yeshaveocc---- START')
          endif
          call ga_print(g_amat(im))
          if (ga_nodeid().eq.0) then
           write(*,4) im
 4         format('-----g_amat(',i3,')-yeshaveocc---- END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,33) ip
 33         format('-----g_amat(',i3,')-yeshaveocc---- START')
          endif
          call ga_print(g_amat(ip))
          if (ga_nodeid().eq.0) then
           write(*,34) ip
 34         format('-----g_amat(',i3,')-yeshaveocc---- END')
          endif

        endif ! end-if-debug

       call ga_zero(g_eig1)
       call ga_matmul_patch('n','n', two,0d0,
     &                      g_vectors ,1,naos,1,nmo,
     &                      g_amat(im),1,nmo ,1,nocc,
     &                      g_eig1    ,1,naos,1,nocc)

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) '-----g_eig1-yeshaveocc---- START'
          call ga_print(g_eig1)
          if (ga_nodeid().eq.0)
     &     write(*,*) '-----g_eig1-yeshaveocc---- END'
         endif ! end-if-debug
      endif ! .not.haveocc
c     note: the dimensioning for array B is that of the transposed
c     matrix, not of the original matrix. 
      
c     calculate C(0)C(-,dagger), store in g_ptmp

      call ga_zero(g_ptmp)
      call ga_matmul_patch('n','t', 1d0,0d0,
     &                     g_vectors,1,naos,1,nocc,
     &                     g_eig1   ,1,nocc,1,naos,
     &                     g_ptmp   ,1,naos,1,naos)   
  
c     -----------------------------
c     Second step: C(+) n C(dagger)
c     -----------------------------     
c     calculate C(+)
      two = 2.0d0
      if (.not.haveocc) then

       call ga_zero(g_eig1)
       call ga_matmul_patch('n','n', two,0d0,
     &                      g_vectors ,1,naos,nocc+1,nmo,
     &                      g_amat(ip),1,nvir,1     ,nocc,
     &                      g_eig1    ,1,naos,1     ,nocc)

       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '-----g_eig1-nohaveocc-2-- START'
        call ga_print(g_eig1)
        if (ga_nodeid().eq.0)
     &   write(*,*) '-----g_eig1-nohaveocc-2-- END'
       endif
      else

       call ga_zero(g_eig1)
       call ga_matmul_patch('n','n', two,0d0,
     &                      g_vectors ,1,naos,1,nmo,
     &                      g_amat(ip),1,nmo ,1,nocc,
     &                      g_eig1    ,1,naos,1,nocc)

       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '-----g_eig1-yeshaveocc-2---- START'
        call ga_print(g_eig1)
        if (ga_nodeid().eq.0)
     &   write(*,*) '-----g_eig1-yeshaveocc-2---- END'
       endif ! end-if-debug
      endif ! end-if-haveocc
      
c     calculate C(+)C(0,dagger), store in g_work

      call ga_zero(g_work)
      call ga_matmul_patch('n','t', 1d0,0d0,
     &                     g_eig1   ,1,naos,1,nocc,
     &                     g_vectors,1,nocc,1,naos,
     &                     g_work   ,1,naos,1,naos)

c     add the two terms together and transpose the density matrix
      if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- CC1^t-------- START'
             call ga_print(g_ptmp)
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- CC1^t--------- END'           
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- C1C^t-------- START'
             call ga_print(g_work)
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- C1C^t--------- END'     
      endif ! end-if-debug
c Doing: C(Du_q)^T + (Du_p)C^T
c     (q,p)=(im,ip) --> P(+)
c     (q,p)=(ip,im) --> P(-)
      call ga_add(1d0,g_ptmp,1d0,g_work,g_work) 
      call ga_transpose(g_work, g_ptmp)
c     calculate symmetrized and antisymmetrized part (+ component)
c     if requested on input:
      if (lantisym) then
        call ga_transpose(g_ptmp,g_work)
        call ga_add(half,g_ptmp, half,g_work,g_pmats)
      else
        call ga_copy(g_ptmp, g_pmats)
      endif
      if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-------- START'
             call ga_print(g_pmats)
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-------- END'    
c             if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmata-------- START'
c             call ga_print(g_pmata)
c            if (ga_nodeid().eq.0)
c     &       write(*,*) '---- g_pmata-------- END'           
      endif ! end-if-debug  
      if (.not.ga_destroy(g_ptmp))
     &   call errquit('get_PertDens: ga_destroy failed g_ptmp',
     &                0,GA_ERR)
      if (.not.ga_destroy(g_work))
     &   call errquit('get_PertDens: ga_destroy failed g_work',
     &                0,GA_ERR)
      if (.not.ga_destroy(g_eig1))
     &   call errquit('get_PertDens: ga_destroy failed g_eig1',
     &                0,GA_ERR) 

      return
      end    

c ++++++++++++++++++ Added for debugging ++++++++++++++++ START

      subroutine CalcPerturbedTDPmat1_deb
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
      external get_PertDens

c     ==================================================================
c ++++++++++++ WARNING -CHECKING lantisym=T +++++++ START
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'FA-WARNING: check lantisym=T'
c       lantisym=.true.
c ++++++++++++ WARNING -CHECKING lantisym=T +++++++ TRUE

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

      if (ga_nodeid().eq.0) then
       write(*,190) ip,im,lstatic,imag
 190     format('(ip,im,lstatic,imag)=(',i3,',',i3,',',L1,',',L1,')')
      endif
        
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

       if (ga_nodeid().eq.0)
     &   write(*,31) imag,two
 31       format('FA-debug (imag,two)=(',L1,',',f15.8)

      if (.not.haveocc) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-enter-no-haveocc...'
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-nohaveocc---- START'
          call ga_print(g_vectors)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-nohaveocc---- END'
          if (ga_nodeid().eq.0) then
           write(*,1) im
 1         format('-----g_amat(',i3,')-nohaveocc---- START')
          endif
          call ga_print(g_amat(im))
          if (ga_nodeid().eq.0) then
           write(*,2) im
 2         format('-----g_amat(',i3,')-nohaveocc---- END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,1070) ip
 1070         format('-----g_amat(',i3,')-nohaveocc---- START')
          endif
          call ga_print(g_amat(ip))
          if (ga_nodeid().eq.0) then
           write(*,1071) ip
 1071        format('-----g_amat(',i3,')-nohaveocc---- END')
          endif
       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,nocc+1,nmo,
     &   g_amat(im),1,nvir,1     ,nocc,
     &   g_eig1    ,1,naos,1     ,nocc)

         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_eig1-nohaveocc---- START'
          call ga_print(g_eig1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_eig1-nohaveocc---- END'

      else
         if (ga_nodeid().eq.0)
     &    write(*,*) 'FA-enter-yes-haveocc...'

         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-yeshaveocc---- START'
          call ga_print(g_vectors)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_vectors-yeshaveocc---- END'
          if (ga_nodeid().eq.0) then
           write(*,3) im
 3         format('-----g_amat(',i3,')-yeshaveocc---- START')
          endif
          call ga_print(g_amat(im))
          if (ga_nodeid().eq.0) then
           write(*,4) im
 4         format('-----g_amat(',i3,')-yeshaveocc---- END')
          endif

       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,1,nmo,
     &   g_amat(im),1,nmo ,1,nocc,
     &   g_eig1    ,1,naos,1,nocc)

         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_eig1-yeshaveocc---- START'
         call ga_print(g_eig1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-----g_eig1-yeshaveocc---- END'
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

        if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-enter-no-haveocc-2...'

       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,nocc+1,nmo,
     &   g_amat(ip),1,nvir,1     ,nocc,
     &   g_eig1    ,1,naos,1     ,nocc)

       if (ga_nodeid().eq.0)
     &  write(*,*) '-----g_eig1-nohaveocc-2-- START'
       call ga_print(g_eig1)
       if (ga_nodeid().eq.0)
     &  write(*,*) '-----g_eig1-nohaveocc-2-- END'

      else

         if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-enter-yes-haveocc-2...'

       call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,1,nmo,
     &   g_amat(ip),1,nmo ,1,nocc,
     &   g_eig1    ,1,naos,1,nocc)

        if (ga_nodeid().eq.0)
     &   write(*,*) '-----g_eig1-yeshaveocc-2---- START'
        call ga_print(g_eig1)
        if (ga_nodeid().eq.0)
     &   write(*,*) '-----g_eig1-yeshaveocc-2---- END'

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
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- CC1^t-------- START'
             call ga_print(g_ptmp)
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- CC1^t--------- END'           
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- C1C^t-------- START'
             call ga_print(g_work)
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- C1C^t--------- END'           

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
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-------- START'
             call ga_print(g_pmats(1))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-------- END'      
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-------- START'
             call ga_print(g_pmata(1))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-------- END'        
      if (ga_nodeid().eq.0)
     & write(*,*) 'FA-BEF-get_PertDens'
      if (ga_nodeid().eq.0) then
       write(*,1000) ip,im,lstatic,imag,
     &               haveocc,lantisym,ncomp,naos,nmo,nocc,nvir
 1000  format('(ip,im,lstatic,imag,haveocc,lanti)=(',
     &       i2,',',i2,',',L1,',',L1,',',L1,',',L1,') ',
     &       ' (ncomp,naos,nmo,nocc,nvir)=(',
     &       i4,',',i4,',',i4,',',i4,',',i4,')')
      endif

      call get_PertDens(
     &               g_pmats(1),!out:      symmetrized dens
     &               g_pmata(1),!out: anti symmetrized dens
     &               g_amat,    ! in: the u mat
     &               g_vectors, ! in: MO vect
     &               ncomp,     ! in: nr. components g_amat
     &               ip,im,     ! in: indices of u vect
     &               imag,      ! in: = T -> imag
     &               haveocc,   ! in: logical flag
     &               lantisym,  ! in: logical flag
     &               naos,nmo,  ! in: nr. AOs,MOs
     &               nocc,nvir, ! in: nr. (occ,virt) MOs
     &               debug)     ! in: =.true. -> show debug printouts

             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-nw-------- START'
             call ga_print(g_pmats(1))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-nw-------- END'      
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-nw-------- START'
             call ga_print(g_pmata(1))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-nw-------- END'        
      if (ga_nodeid().eq.0)
     & write(*,*) 'FA-AFT-get_PertDens'
c ====> WARNING-- START
c -- Commenting lines below to enter 2nd part of calc
      if (ga_nodeid().eq.0) write(*,*) 'SKIP lstatic-if'
c ====> WARNING-- END
c      if (lstatic .or. ncomp.eq.1) then
c       skip calculation of component 2 of the density matrix
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-Skipping calc of 2nd component'

c        goto 7000
c      endif  
      
        if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-Doing calc of 2nd component'     
 
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
     &   g_vectors ,1,naos,nocc+1,nmo,
     &   g_amat(ip),1,nvir,1,nocc,
     &   g_eig1    ,1,naos,1,nocc)
      else
      call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors ,1,naos,1,nmo,
     &   g_amat(ip),1,nmo,1,nocc,
     &   g_eig1    ,1,naos,1,nocc)
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
     &   g_vectors,1,naos,nocc+1,nmo,
     &   g_amat(im),1,nvir,1,nocc,
     &   g_eig1,1,naos,1,nocc)
      else
      call ga_matmul_patch('n','n', two,0d0,
     &   g_vectors,1,naos,1,nmo,
     &   g_amat(im),1,nmo,1,nocc,
     &   g_eig1,1,naos,1,nocc)
      endif
      call ga_sync()

      if (debug) write (6,*) '7'

c     calculate C(-)C(0,dagger), store in g_work
      call ga_matmul_patch('n','t', 1d0,0d0,
     &   g_eig1,1,naos,1,nocc,
     &   g_vectors,1,nocc,1,naos,
     &   g_work,1,naos,1,naos)
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
        call ga_add(half,g_ptmp, half,g_work,g_pmats(2))
        call ga_sync()
        call ga_add(half,g_ptmp,-half,g_work,g_pmata(2))
        call ga_sync()
      else
        call ga_copy(g_ptmp, g_pmats(2))
        call ga_sync()
      endif
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-2-------- START'
             call ga_print(g_pmats(2))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-2-------- END'      
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-2-------- START'
             call ga_print(g_pmata(2))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-2-------- END' 
      if (ga_nodeid().eq.0)
     & write(*,*) 'FA-BEF-get_PertDens-2'
      if (ga_nodeid().eq.0) then
       write(*,1001) ip,im,lstatic,imag,
     &               haveocc,lantisym,ncomp,naos,nmo,nocc,nvir
 1001  format('(ip,im,lstatic,imag,haveocc,lanti)=(',
     &       i2,',',i2,',',L1,',',L1,',',L1,',',L1,') ',
     &       ' (ncomp,naos,nmo,nocc,nvir)=(',
     &       i4,',',i4,',',i4,',',i4,',',i4,')')
      endif

      call get_PertDens(
     &               g_pmats(2),!out:      symmetrized dens
     &               g_pmata(2),!out: anti symmetrized dens
     &               g_amat,    ! in: the u mat
     &               g_vectors, ! in: MO vect
     &               ncomp,     ! in: nr. components g_amat
     &               im,ip,     ! in: indices of u vect
     &               imag,      ! in: = T -> imag
     &               haveocc,   ! in: logical flag
     &               lantisym,  ! in: logical flag
     &               naos,nmo,  ! in: nr. AOs,MOs
     &               nocc,nvir, ! in: nr. (occ,virt) MOs
     &               debug)     ! in: =.true. -> show debug printouts

             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-2-nw-------- START'
             call ga_print(g_pmats(2))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-2-nw-------- END'      
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-2-nw-------- START'
             call ga_print(g_pmata(2))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmata-2-nw-------- END'        
      if (ga_nodeid().eq.0)
     & write(*,*) 'FA-AFT-get_PertDens'

      stop
      
c     ---------------------------------------------
c     deallocate temporary arrays, sync, and return
c     ---------------------------------------------

c     jump here from above in case of static calculation
 7000 continue

      if (ga_nodeid().eq.0)
     & write(*,*) 'FA-AFT-get_PertDens-did-not-enter-2'

      stop      

      
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

c ++++++++++++++++++ Added for debugging ++++++++++++++++ END 
