*
* $Id$
*
c     ga_pzheev
c
c     interface into scalapack's Hermitian eigensolver
c     
c     ga interface to scalapack
c     
c     g_a .g_c(*,i) = eval(i).g_b(*,i)
c     
c     assume equal size blocks ...
c     
c******************************************************************
c     
      subroutine ga_pzheev(n, g_a, eval, g_b)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "slface.fh"
c     
      integer g_a               ! matrix A
      integer g_b               ! matrix B
      integer n
c
c all eigenvalues are real
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, range, uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimB1, dimB2, typeB
c
      integer mpA, nqA          ! rows/cols of A held by the processor
      integer mpB, nqB          ! rows/cols of B held by the processor
c
      integer info
      integer me
      integer nwork
      integer n,lda, ldb, elemA, elemB
      integer numroc
      
      integer nb                ! block size
      integer descA(9), descB(9), descC(9) !descriptor for scalapack
      
c     
      integer ngaps, hgap, adrgaps
      integer iclu, hclustr, adrclustr
      integer if, hfail, adrfail
      integer liwork, hiwork, adriwork
      integer lrwork, hrwork, adrrwork
      integer lcwork, hcwork, adrcwork
c
      integer nn, mq0, np0
      double precision vl, vu, abstol, orfac
      integer il, iu
      integer m, nz
      integer info

      integer iceil, iii
      double precision pdlamch
      external pdlamch
      external pzheevx
      external iceil
      intrinsic max
c     
c     processor dependent; machine dependent
c     
      data nb/2/
c     
c***  check environment
c     
      me     = ga_nodeid()
      write(*,*) ' me ', me
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a, 'ga_pzheev: a')
      call ga_check_handle(g_b, 'ga_pzheev: b')
c     
      call ga_inquire(g_a, typeA, dimA1, dimA2)
      call ga_inquire(g_b, typeB, dimB1, dimB2)
c     
c     
c***  initialize SL interface
c     
      call SLinit()
c     
c***  find SBS format parameters
c     
c     
      mpA = numroc(dimA1, nb, myrow, 0, nprow)
      nqA = numroc(dimA2, nb, mycol, 0, npcol)
c     
      mpB = numroc(dimB1, nb, myrow, 0, nprow)
      nqB = numroc(dimB2, nb, mycol, 0, npcol) 
c     
c     
      lda = max(1,mpA)
      ldb = max(1,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
      oactive = .true.
      if ( mpA*nqA .eq. 0 ) oactive = .false.
      if ( mpB*nqB .eq. 0 ) oactive = .false.
c     
      call ga_sync()
      write(*,*) ' me ', me, ' mpA*nqA , mpB*nqB ', mpA*nqA, mpB*nqB
c     
      if (oactive) then
c     
         elemA= mpA*nqA
         status = .true.
         if(elemA.ne.0)status =
     $        ma_push_get(MT_DCPL,elemA,'a',ha,adra)
         if(.not.status) 
     &        call ga_error('ga_pzheev: mem alloc failed A ', -1)
c     
         elemB= mpB*nqB
c     
         if(elemB.ne.0)status =
     $        ma_push_get(MT_DCPL,elemB,'b',hb,adrb)
         if(.not.status) 
     &        call ga_error('ga_pzheev: mem alloc failed B ', -1)
c     
         ngaps = nprow*npcol
         if(ngaps.ne.0)status =
     $        ma_push_get(MT_DBL, ngaps ,'gap',hgap,adrgaps)
         if(.not.status) 
     &        call ga_error('ga_pzheev: mem alloc failed gaps  ', -1)
c     
         iclu = 2*nprow*npcol
         iclu = 2*n
         if(iclu.ne.0)status =
     $        ma_push_get(MT_INT, iclu ,'iclustr',hclustr,adrclustr)
         if(.not.status) 
     &        call ga_error('ga_pzheev: mem alloc failed iclustr  ', -1)
c     
         if = n
         if(if.ne.0) then
            status =
     $           ma_push_get(MT_INT, if ,'ifail',hfail,adrfail)
            write(*,*) ' allocating ifail ', if
         endif
         if(.not.status)
     &        call ga_error('ga_pzheev: mem alloc failed ifail  ', -1)
         
c     
c***  copy g_a to A and g_b to B using the block cyclic scalapack format 
c     
c     write(*,*) ' ga_to_SL '
c     call ffflush(6)
c     
         write(*,*) ' before ga_to_zsl ', me
         
         call ga_to_ZSL(g_a, dimA1, dimA2, nb, nb,
     $        dcpl_mb(adrA), lda, mpA, nqA)
c     
         call ga_to_ZSL(g_b, dimB1, dimB2, nb,nb,
     $        dcpl_mb(adrB), ldb, mpB, nqB)
c     
         write(*,*) ' after ga_to_zsl ', me
c     
c     
c***  fill SCALAPACK matrix descriptors
c     
         write(*,*) ' before descinit ', me
         call descinit(descA, dimA1, dimA2, nb, nb, 0, 0,
     $        iSLctxt, lda, info)
         if(info.ne.0) call ga_error(' ga_pzheev: descinit A failed ',
     .        -info)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, 0, 0,
     .        iSLctxt, ldb, info)
         if(info.ne.0) call ga_error(' ga_pzheev: descinit B failed ',
     .        -info)
         write(*,*) ' after descinit ', me
c     
c***  pdgemm
c     
         write(*,*) ' before ga_sync ', me
         write(*,*) ' after ga_sync ', me
c     
c     call scalapack
c     
c     
         jobz = 'V'
         range = 'A'
         uplo = 'L'
c     abstol = 2.0d0*pdlamch('S')
         abstol = 0.0d0
         vl = 0.d0
         vu = 0.d0
         il = 0
         iu = 0
         nz = 0
c     
c     ability to deal with orthonormality ; let's just
c     have the regular scalapack stuff for the moment
c     
         liwork = 6*max(n, nprow*npcol+1, 4)
         if(liwork.ne.0)status =
     $        ma_push_get(MT_INT, liwork ,'iwork',hiwork,adriwork)
         if(.not.status) 
     &        call ga_error('ga_pzheev: mem alloc failed iwork  ', -1)
c     
         nn = max(n, nb, 2)
         np0 = numroc(nn, nb, 0, 0, nprow)
         mq0 = numroc(nn, nb, 0, 0, npcol)
c     
         write(*,*) ' before lrwork , me n', me, n
         call ffflush(6)
c     
         orfac = 1.d-3
c     
c     
         lcwork = n + (NP0 + MQ0 + NB)*NB
c     
         if(lcwork.ne.0)
     $        status = ma_push_get(MT_DCPL, lcwork ,
     $        'cwork',hcwork,adrcwork)
         if(.not.status) 
     &        call ga_error('ga_pzheev: mem alloc failed cwork  ', -1)
c     
         lrwork = iceil(n, nprow*npcol)*nn
         lrwork = 4*n + max(5*NN, NP0*MQ0) + lrwork + n*n
         write(*,*) ' after lrwork ', me
         call ffflush(6)
c     
         if(lrwork.ne.0)
     $        status = ma_push_get(MT_DCPL, lrwork ,
     $        'rwork',hrwork,adrrwork)
         if(.not.status) 
     &        call ga_error('ga_pzheev: mem alloc failed rwork  ', -1)
c     
         write(*,*) ' before abstol ', me, ' liwork ', liwork
         call ffflush(6)
c     
c     x      abstol = 2.d0 * pdlamch('s')
c     
         abstol = 0.0d0
c     
c     
         write(*,*) ' before pzheevx ', me, ' liwork ', liwork
         call ffflush(6)
c     
         call pzheevx(jobz, range, uplo,
     $        n, dcpl_mb(adrA), 1, 1, descA,
     $        vl, vu, il, iu,
     $        abstol, m, nz, eval, orfac,
     $        dcpl_mb(adrB), 1, 1, descB, dcpl_mb(adrcwork), lcwork,
     $        dbl_mb(adrrwork), lrwork, int_mb(adriwork),
     $        liwork, int_mb(adrfail),
     $        int_mb(adrclustr), dbl_mb(adrgaps), info)
c     
         if ( info .ne.0)then
            write(*,*) ' me ', me, 'info =', info
         endif
c     
         write(*,*) ' after pzheevx ', me
         call flush(6)
c
         if ( me . eq. 0 ) then
            do iii = 1, n
               write(*,*) iii, eval(iii)
            enddo
         endif
         
c     
         if (nz .ne. n ) then
            if ( info .ne. 0 ) then
               if ( info .gt. 0 ) then
         call ga_error(' ga_pzheev: argument is illegal ', info)
               else
         call ga_error(' ga_pzheev: eigenvectors failed to converge ',
     $                 info)
               endif
            endif
         endif
c     
c     
c     
c***  copy solution matrix back to g_c
c     
         write(*,*) ' before ga_from_zsl ', me
         call ffflush(6)
c     
         call ga_from_ZSL(g_b, dimA1, dimB2,
     $        nb, nb, dcpl_mb(adrB),
     &        mpB, ldb, nqB)
c     
         write(*,*) ' after ga_from_zsl ', me
         call ffflush(6)
c     
c     
         call blacs_gridexit(iSLctxt)
c     
c***  deallocate work/SL arrays
c     
         if ( lrwork .ne. 0 ) status = ma_pop_stack(hrwork)
         if ( lcwork .ne. 0 ) status = ma_pop_stack(hcwork)
         if ( liwork .ne. 0 ) status = ma_pop_stack(hiwork)
         if ( if .ne. 0 ) status = ma_pop_stack(hfail)
         if ( iclu .ne. 0 ) status = ma_pop_stack(hclustr)
         if ( ngaps.ne.0 ) status = ma_pop_stack(hgap)
         if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
         if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
      endif
c     
      call ga_sync()
      write(*,*) ' exiting ga_pzheev ', me
      call ffflush(6)
      return
      end
      




