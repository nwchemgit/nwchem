      subroutine tddft_iguess(rtdb,nbf_ao,ipol,
     1  nocc,nmo,nroots,maxvecs,g_trials,g_movecs,evl,
     2  ntrials,nproducts)
c
c $Id: tddft_iguess.F,v 1.1 2002-03-22 02:10:36 sohirata Exp $
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
c
      integer rtdb             ! Runtime database handle
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nroots           ! Number of roots sought
      integer maxvecs          ! Maximum number of trial vectors
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_movecs(2)      ! GA handle for MO eigenvectors
      double precision evl(nbf_ao,2)   ! Eigenvalues
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
c
      double precision tiny,huge,swap
      parameter (tiny=1.0d-6)  ! Tolerance for judging near degeneracy
      parameter (huge=1.0d99)  ! A number greater than eigenvalue differences
      integer nov(2)           ! Number of occupied-virtual pairs
      integer l_ediff,k_ediff  ! MA handles for work array
      integer l_sorto,k_sorto  ! MA handles for work array
      integer l_sortv,k_sortv  ! MA handles for work array
      integer l_sortp,k_sortp  ! MA handles for work array
      integer l_sorta,k_sorta  ! MA handles for work array
      integer l_t(2),k_t(2)    ! MA handles for work array
      integer i,j,k,l,m,n
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c --------------------
c Allocate work arrays
c --------------------
c
      do i=1,ipol
        nov(i)=(nmo(i)-nocc(i))*nocc(i)
      enddo
      if (ipol.eq.1) nov(2)=0
      if (.not.ma_push_get(mt_dbl,nov(1)+nov(2),
     1  'eigenvalue differences',l_ediff,k_ediff)) call errquit
     2  ('tddft_iguess: failed to allocate ediff',0)
c Consider nroots+20 smallest eigenvalue differences
c Extra 20 are there in case the nroots-th vector is degenarate
c In case of degenaracy, we must increase nroots
      if (.not.ma_push_get(mt_int,nroots+20,'occupied index',
     1  l_sorto,k_sorto)) call errquit
     2  ('tddft_iguess: failed to allocate sorto',0)
      if (.not.ma_push_get(mt_int,nroots+20,'virtual index',
     1  l_sortv,k_sortv)) call errquit
     2  ('tddft_iguess: failed to allocate sortv',0)
      if (.not.ma_push_get(mt_int,nroots+20,'spin index',
     1  l_sortp,k_sortp)) call errquit
     2  ('tddft_iguess: failed to allocate sortp',0)
      if (.not.ma_push_get(mt_int,nroots+20,'composite index',
     1  l_sorta,k_sorta)) call errquit
     2  ('tddft_iguess: failed to allocate sorta',0)
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nov(i),'Work trial vector',
     1    l_t(i),k_t(i))) call errquit
     2    ('tddft_iguess: failed to allocate t',0)
      enddo
c
c ----------------------------------
c Compute orbital energy differences
c ----------------------------------
c
      l=0
      do i=1,ipol
        do j=1,nocc(i)
          do k=nocc(i)+1,nmo(i)
            l=l+1
            dbl_mb(k_ediff+l-1)=evl(k,i)-evl(j,i)
          enddo
        enddo
      enddo
c
c --------------------------------------------
c Find 2*nroot smallest eigenvalue differences
c --------------------------------------------
c
      n=min(nroots+20,nov(1)+nov(2))
      do m=1,n
        swap=huge
        l=0
        do i=1,ipol
          do j=1,nocc(i)
            do k=nocc(i)+1,nmo(i)
              l=l+1
              if (dbl_mb(k_ediff+l-1).lt.swap) then
                int_mb(k_sorto+m-1)=j
                int_mb(k_sortv+m-1)=k
                int_mb(k_sortp+m-1)=i
                int_mb(k_sorta+m-1)=l
                swap=dbl_mb(k_ediff+l-1)
              endif
            enddo
          enddo
        enddo
        dbl_mb(k_ediff+int_mb(k_sorta+m-1)-1)=huge+huge
      enddo
c
c ---------------------------------------------------
c Print the smallest eigenvalue differences for debug
c ---------------------------------------------------
c
      if (nodezero.and.util_print('initial guess',print_high)) then
        write(LuOut,*)
        write(LuOut,9000) n
        write(LuOut,9010)
        write(LuOut,9020)
        do m=1,n
          write(LuOut,9030) m,int_mb(k_sortp+m-1),
     1      int_mb(k_sorto+m-1),int_mb(k_sortv+m-1),
     2      evl(int_mb(k_sortv+m-1),int_mb(k_sortp+m-1)),
     3      evl(int_mb(k_sorto+m-1),int_mb(k_sortp+m-1)),
     4      evl(int_mb(k_sortv+m-1),int_mb(k_sortp+m-1))-
     5      evl(int_mb(k_sorto+m-1),int_mb(k_sortp+m-1))
        enddo
        write(LuOut,9020)
      endif
 9000 format(i5,' smallest eigenvalue differences')
 9010 format('  No. Spin  Occ  Vir    E(Vir)    E(Occ)   E(Diff)')
 9020 format('--------------------------------------------------')
 9030 format(4i5,3f10.5)
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
     1    call errquit('tddft_iguess: failed writing to rtdb',0)
      else if (nroots.lt.nov(1)+nov(2)) then
        i=nroots
   10   continue
        if (evl(int_mb(k_sortv+i),int_mb(k_sortp+i))-
     1      evl(int_mb(k_sorto+i),int_mb(k_sortp+i))-
     2      evl(int_mb(k_sortv+i-1),int_mb(k_sortp+i-1))+
     3      evl(int_mb(k_sorto+i-1),int_mb(k_sortp+i-1))
     4      .lt.tiny) then
          i=i+1
          if (i.ge.n) call errquit
     1      ('tddft_iguess: increase nroots',nroots)
        else
          goto 20
        endif
        goto 10
   20   continue
        if (i.gt.nroots) then
          if (nodezero.and.util_print('initial guess',print_default))
     1      write(LuOut,9050) nroots,i
          nroots=i
          if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1      call errquit('tddft_iguess: failed writing to rtdb',0)
        endif
      endif
 9040 format(10x,'Number of roots decreased from  ',i8,' to ',i8)
 9050 format(10x,'Number of roots increased from  ',i8,' to ',i8)
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
     1      dbl_mb(k_t(i)+(int_mb(k_sorto+m-1)-1)*(nmo(i)-nocc(i))
     2                   +(int_mb(k_sortv+m-1)-nocc(i))-1)=1.0d0
          call ga_put(g_trials(i),1,nov(i),m,m,dbl_mb(k_t(i)),1)
        enddo
      enddo
      ntrials=nroots
      nproducts=0
      if (nodezero.and.util_print('initial guess',print_high))
     1  write(LuOut,9060) ntrials
 9060 format(10x,'Number of initial trial vectors generated ',i8)
c           
c ------
c Return
c ------
c
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_t(i))) call errquit
     1    ('tddft_iguess: failed to deallocate t',0)
      enddo
      if (.not.ma_pop_stack(l_sorta)) call errquit 
     1  ('tddft_iguess: failed to deallocate sorta',0)
      if (.not.ma_pop_stack(l_sortp)) call errquit 
     1  ('tddft_iguess: failed to deallocate sortp',0)
      if (.not.ma_pop_stack(l_sortv)) call errquit 
     1  ('tddft_iguess: failed to deallocate sortv',0)
      if (.not.ma_pop_stack(l_sorto)) call errquit 
     1  ('tddft_iguess: failed to deallocate sorto',0)
      if (.not.ma_pop_stack(l_ediff)) call errquit 
     1  ('tddft_iguess: failed to deallocate ediff',0)
      return
      end
