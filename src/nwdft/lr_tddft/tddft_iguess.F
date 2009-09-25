      subroutine tddft_iguess(rtdb,geom,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_buffer1,g_buffer2,
     2  d_trials,g_movecs,evl,ntrials,nproducts,algorithm,targetsym,
     3  symmetry,l_irs,k_irs,tda)
c
c $Id: tddft_iguess.F,v 1.11 2003-10-17 22:58:29 carlfahl Exp $
c
c Generate simple initial guess vectors. First sorts 
c the eigenvalue differences (epsilon_a - epsilon_i).
c Then in an ascending order of the eigenvalue differences,
c each of the corresponding singly excited determinants is
c a trial vector.  The number of initial trial vectors is
c equal to nroots, unless there is degeneracy at the nroots-th
c initial trial vector.
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
c
      integer rtdb             ! Runtime database handle
      integer geom             ! Geometry handle
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer nroots           ! Number of roots sought
      integer maxvecs          ! Maximum number of trial vectors
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer g_buffer2(2)     ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      integer g_movecs(2)      ! GA handle for MO eigenvectors
      double precision evl(nbf_ao,2)   ! Eigenvalues
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
      integer algorithm        ! Algorithm
      character*4 targetsym    ! Target root irrep
      logical symmetry         ! Symmetry restriction
      integer l_irs(2),k_irs(2)! Irreps
      logical tda              ! True if Tamm-Dancoff approximation
c
      integer buf              ! Number of extra initial sorted eval diffs
      double precision tiny,huge,swap
      parameter (buf=0)        ! Right now it is turned off
      parameter (tiny=1.0d-6)  ! Tolerance for judging near degeneracy
      parameter (huge=1.0d99)  ! A number greater than eigenvalue differences
      integer nov(2)           ! Number of occupied-virtual pairs
      integer l_ediff,k_ediff  ! MA handles for work array
      integer l_sorto,k_sorto  ! MA handles for work array
      integer l_sortv,k_sortv  ! MA handles for work array
      integer l_sortp,k_sortp  ! MA handles for work array
      integer l_sorta,k_sorta  ! MA handles for work array
      integer l_t(2),k_t(2)    ! MA handles for work array
      integer l_irrep,k_irrep  ! MA handles for work array
      integer isym,jsym
      character*4 oname
      integer i,j,k,l,m,n
      logical nodezero
      integer request
      integer dummy
      logical found
C     character*80 vector
C     logical vector_exists
      integer ipol_verify
      integer nroots_verify
      integer nov_verify(2)
      logical tda_verify
      double precision r
c
      nodezero=(ga_nodeid().eq.0)
      do i=1,ipol
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
      enddo
      if (ipol.eq.1) nov(2)=0
c
c ---------------------------------------
c Read initial guess from a file (if any)
c ---------------------------------------
c
C     if (.not.rtdb_cget(rtdb,'tddft:vector',1,vector))
C    1  call errquit('tddft_iguess: failed to read vector',0)
C     inquire(file=vector,exist=vector_exists)
C     if (vector_exists) then
C       open(unit=69,file=vector,form='formatted',status='old',
C    1  err=90000,access='sequential')
C       rewind(69)
C       read(69,*,err=90000) tda_verify
C       read(69,*,err=90000) ipol_verify
C       read(69,*,err=90000) nroots_verify
C       read(69,*,err=90000) nov_verify(1),nov_verify(2)
C       if ((ipol.eq.ipol_verify).and.
C    1      (nroots.le.nroots_verify).and.
C    2      (nov(1).eq.nov_verify(1)).and.
C    3      (nov(2).eq.nov_verify(2))) then
c
c Read from a file
c
C         if (tda_verify) then
C           ntrials=nroots
C         else
C           ntrials=2*nroots
C         endif
C         do i=1,ipol
C           do n=1,ntrials
C             do m=1,nov(i)
C               read(69,*) r
C               if ((algorithm.eq.1).or.(algorithm.eq.2)) then
C                 call ga_put(g_trials(i),m,m,n,n,r,1)
C               else
C                 call ga_put(g_buffer1(i),m,m,1,1,r,1)
C               endif
C             enddo
C             if ((algorithm.eq.3).or.(algorithm.eq.4)) then
C               dummy=dra_write_section
C    1            (.false.,g_buffer1(i),1,nov(i),1,1,
C    2            d_trials(i),1,nov(i),n,n,request)
C               dummy=dra_wait(request)
C             endif
C           enddo
C         enddo
C         nproducts=0
C
C Orthonormalize
C
C         do i=1,ipol
C           do n=1,ntrials
C             if ((algorithm.eq.1).or.(algorithm.eq.2)) then
C               r=ga_ddot_patch(g_trials(i),'N',1,nov(i),n,n,
C    1                          g_trials(i),'N',1,nov(i),n,n)
C             else
C               dummy=dra_read_section(.false.,g_buffer1(i),
C    1          1,nov(i),1,1,d_trials(i),1,nov(i),n,n,request)
C               dummy=dra_wait(request)
C               r=ga_ddot_patch(g_buffer1(i),'N',1,nov(i),1,1,
C    1                          g_buffer1(i),'N',1,nov(i),1,1)
C             endif
C             if (r.eq.0.0d0) call errquit
C    1          ('tddft_iguess: unable to normalize',0)
C             r=1.0d0/dsqrt(r)
C             if ((algorithm.eq.1).or.(algorithm.eq.2)) then
C               call ga_scale_patch(g_trials(i),1,nov(i),n,n,r)
C             else
C               call ga_scale_patch(g_buffer1(i),1,nov(i),1,1,r)
C               dummy=dra_write_section(.false.,g_buffer1(i),
C    1          1,nov(i),1,1,d_trials(i),1,nov(i),n,n,request)
C               dummy=dra_wait(request)
C             endif
C             if (n.ne.ntrials) then
C               do m=n+1,ntrials
C                 if ((algorithm.eq.1).or.(algorithm.eq.2)) then
C                   r=ga_ddot_patch(g_trials(i),'N',1,nov(i),n,n,
C    1                              g_trials(i),'N',1,nov(i),m,m)
C                   call ga_add_patch(1.0d0,g_trials(i),1,nov(i),m,m,
C    1                                -r,   g_trials(i),1,nov(i),n,n,
C    2                                      g_trials(i),1,nov(i),m,m)
C                 else
C                   dummy=dra_read_section(.false.,g_buffer1(i),
C    1              1,nov(i),1,1,d_trials(i),1,nov(i),n,n,request)
C                   dummy=dra_wait(request)
C                   dummy=dra_read_section(.false.,g_buffer2(i),
C    1              1,nov(i),1,1,d_trials(i),1,nov(i),m,m,request)
C                   dummy=dra_wait(request)
C                   r=ga_ddot_patch(g_buffer1(i),'N',1,nov(i),1,1,
C    1                              g_buffer2(i),'N',1,nov(i),1,1)
C                   call ga_add_patch(1.0d0,g_buffer2(i),1,nov(i),1,1,
C    1                                -r,   g_buffer1(i),1,nov(i),1,1,
C    2                                      g_buffer2(i),1,nov(i),1,1)
C                   dummy=dra_write_section(.false.,g_buffer2(i),
C    1              1,nov(i),1,1,d_trials(i),1,nov(i),m,m,request)
C                   dummy=dra_wait(request)
C                 endif
C               enddo
C             endif
C           enddo
C         enddo
C         if (nodezero.and.util_print('initial guess',print_default))
C    1      write(LuOut,1010) vector
C         return
C       endif
C     endif
C1010 format(/,2x,'Initial guess read from ',a)
c
c --------------------
c Allocate work arrays
c --------------------
c
C90000 if (nodezero.and.util_print('initial guess',print_default)
C     1  .and.(vector.ne.'')) write(LuOut,1000) 
C 1000 format(/,2x,'Unable to retrieve initial guess from a file')
      if (.not.ma_push_get(mt_dbl,nov(1)+nov(2),
     1  'eigenvalue differences',l_ediff,k_ediff)) call errquit
     2  ('tddft_iguess: failed to allocate ediff',0, MA_ERR)
c Consider nroots+20 smallest eigenvalue differences
c Extra 20 are there in case the nroots-th vector is degenarate
c In case of degenaracy, we must increase nroots
      if (.not.ma_push_get(mt_int,nroots+buf,'occupied index',
     1  l_sorto,k_sorto)) call errquit
     2  ('tddft_iguess: failed to allocate sorto',0, MA_ERR)
      if (.not.ma_push_get(mt_int,nroots+buf,'virtual index',
     1  l_sortv,k_sortv)) call errquit
     2  ('tddft_iguess: failed to allocate sortv',0, MA_ERR)
      if (.not.ma_push_get(mt_int,nroots+buf,'spin index',
     1  l_sortp,k_sortp)) call errquit
     2  ('tddft_iguess: failed to allocate sortp',0, MA_ERR)
      if (.not.ma_push_get(mt_int,nroots+buf,'composite index',
     1  l_sorta,k_sorta)) call errquit
     2  ('tddft_iguess: failed to allocate sorta',0, MA_ERR)
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nov(i),'Work trial vector',
     1    l_t(i),k_t(i))) call errquit
     2    ('tddft_iguess: failed to allocate t',0, MA_ERR)
      enddo
      if (.not.ma_push_get(mt_int,nroots+buf,'irrep',
     1  l_irrep,k_irrep)) call errquit
     2  ('tddft_iguess: failed to allocate irrep',0, MA_ERR)
c
c ----------------------------------
c Compute orbital energy differences
c ----------------------------------
c
      l=0
      do i=1,ipol
        do j=nfc(i)+1,nocc(i)
          do k=nocc(i)+1,nmo(i)-nfv(i)
            l=l+1
            dbl_mb(k_ediff+l-1)=evl(k,i)-evl(j,i)
          enddo
        enddo
      enddo
c
c ------------------------------
c Compute the ground state irrep
c ------------------------------
c
      isym=1
      if (ipol.eq.2) then
        do i=1,ipol
c here we are just getting ground state symmetry ...
c ... frozen cores are irrelevant
          do j=1,nocc(i)
            isym=ieor(isym-1,int_mb(k_irs(i)+j-1)-1)+1
          enddo
        enddo
      endif
      call sym_irrepname(geom,isym,oname)
      if (nodezero.and.util_print('initial guess',print_debug))
     1  write(LuOut,*) 'Ground state symmetry ',oname
c
c ------------------------------------------
c Find nroot smallest eigenvalue differences
c ------------------------------------------
c
      n=min(nroots+buf,nov(1)+nov(2))
      do m=1,n
        swap=huge
        l=0
        found=.false.
        do i=1,ipol
          do j=nfc(i)+1,nocc(i)
            do k=nocc(i)+1,nmo(i)-nfv(i)
              jsym=ieor(isym-1,int_mb(k_irs(i)+j-1)-1)+1
              jsym=ieor(jsym-1,int_mb(k_irs(i)+k-1)-1)+1
              call sym_irrepname(geom,jsym,oname)
              l=l+1
              if (((.not.symmetry).or.(targetsym.eq.oname))
     1          .and.(dbl_mb(k_ediff+l-1).lt.swap)) then
                int_mb(k_sorto+m-1)=j
                int_mb(k_sortv+m-1)=k
                int_mb(k_sortp+m-1)=i
                int_mb(k_sorta+m-1)=l
                int_mb(k_irrep+m-1)=jsym
                swap=dbl_mb(k_ediff+l-1)
                found=.true.
              endif
            enddo
          enddo
        enddo
        if (.not.found) call errquit
     1    ('tddft_iguess: too many roots requested',nroots, INPUT_ERR)
        dbl_mb(k_ediff+int_mb(k_sorta+m-1)-1)=huge+huge
      enddo
c
c ---------------------------------------------------
c Print the smallest eigenvalue differences for debug
c ---------------------------------------------------
c
      if (nodezero.and.util_print('initial guess',print_debug)) then
        write(LuOut,*)
        write(LuOut,9000) n
        write(LuOut,9010)
        write(LuOut,9020)
        do m=1,n
          call sym_irrepname(geom,int_mb(k_irrep+m-1),oname)
          write(LuOut,9030) m,int_mb(k_sortp+m-1),
     1      int_mb(k_sorto+m-1),int_mb(k_sortv+m-1),oname,
     2      evl(int_mb(k_sortv+m-1),int_mb(k_sortp+m-1)),
     3      evl(int_mb(k_sorto+m-1),int_mb(k_sortp+m-1)),
     4      evl(int_mb(k_sortv+m-1),int_mb(k_sortp+m-1))-
     5      evl(int_mb(k_sorto+m-1),int_mb(k_sortp+m-1))
        enddo
        write(LuOut,9020)
      endif
 9000 format(i5,' smallest eigenvalue differences')
 9010 format('  No. Spin  Occ  Vir  Irrep   E(Vir)    E(Occ)   E(Diff)')
 9020 format('--------------------------------------------------------')
 9030 format(4i5,1x,a4,1x,3f10.5)
c
c --------------------------------------------------------
c Check if nroots is too large, if so, reduce it
c Check if the nroots-th excited state could be degenerate
c If so, increase nroots iteratively
c --------------------------------------------------------
c
      if (nroots.gt.nov(1)+nov(2)) then
        if (nodezero.and.util_print('initial guess',print_default))
     1    write(LuOut,9040) nroots,nov(1)+nov(2)
        nroots=nov(1)+nov(2)
        if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1    call errquit('tddft_iguess: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (nroots.lt.nov(1)+nov(2)) then
c       i=nroots
c  10   continue
c       if (evl(int_mb(k_sortv+i),int_mb(k_sortp+i))-
c    1      evl(int_mb(k_sorto+i),int_mb(k_sortp+i))-
c    2      evl(int_mb(k_sortv+i-1),int_mb(k_sortp+i-1))+
c    3      evl(int_mb(k_sorto+i-1),int_mb(k_sortp+i-1))
c    4      .lt.tiny) then
c         i=i+1
c         if (i.ge.n) call errquit
c    1      ('tddft_iguess: increase nroots',nroots)
c       else
c         goto 20
c       endif
c       goto 10
c  20   continue
c       if (i.gt.nroots) then
c         if (nodezero.and.util_print('initial guess',print_default))
c    1      write(LuOut,9050) nroots,i
c         nroots=i
c         if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
c    1      call errquit('tddft_iguess: failed writing to rtdb',0)
c       endif
      endif
 9040 format(10x,'Number of roots decreased from  ',i8,' to ',i4)
 9050 format(10x,'Number of roots increased from  ',i8,' to ',i4)
c
c ------------------------------
c Generate initial trial vectors
c ------------------------------
c
      do m=1,nroots
        do i=1,ipol
          do j=1,nov(i)
            dbl_mb(k_t(i)+j-1)=0.0d0
          enddo
          if (i.eq.int_mb(k_sortp+m-1))
     1      dbl_mb(k_t(i)+(int_mb(k_sorto+m-1)-nfc(i)-1)
     2        *((nmo(i)-nfv(i))-nocc(i))
     3        +(int_mb(k_sortv+m-1)-nocc(i)-1))=1.0d0
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            call ga_put(g_trials(i),1,nov(i),m,m,dbl_mb(k_t(i)),1)
          else
            call ga_put(g_buffer1(i),1,nov(i),1,1,dbl_mb(k_t(i)),1)
            dummy=dra_write_section(.false.,g_buffer1(i),1,nov(i),1,1,
     1        d_trials(i),1,nov(i),m,m,request)
            dummy=dra_wait(request)
          endif
        enddo
      enddo
      ntrials=nroots
      nproducts=0
      if (nodezero.and.util_print('initial guess',print_debug))
     1  write(LuOut,9060) ntrials
 9060 format(10x,'Number of initial trial vectors generated ',i8)
c           
c ------
c Return
c ------
c
      if (.not.ma_pop_stack(l_irrep)) call errquit 
     1  ('tddft_iguess: failed to deallocate irrep',0, MA_ERR)
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_t(i))) call errquit
     1    ('tddft_iguess: failed to deallocate t',0, MA_ERR)
      enddo
      if (.not.ma_pop_stack(l_sorta)) call errquit 
     1  ('tddft_iguess: failed to deallocate sorta',0, MA_ERR)
      if (.not.ma_pop_stack(l_sortp)) call errquit 
     1  ('tddft_iguess: failed to deallocate sortp',0, MA_ERR)
      if (.not.ma_pop_stack(l_sortv)) call errquit 
     1  ('tddft_iguess: failed to deallocate sortv',0, MA_ERR)
      if (.not.ma_pop_stack(l_sorto)) call errquit 
     1  ('tddft_iguess: failed to deallocate sorto',0, MA_ERR)
      if (.not.ma_pop_stack(l_ediff)) call errquit 
     1  ('tddft_iguess: failed to deallocate ediff',0, MA_ERR)
      return
      end
