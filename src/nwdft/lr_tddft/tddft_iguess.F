      subroutine tddft_iguess(rtdb,geom,nbf_ao,ipol,
     1  nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_buffer1,g_buffer2,
     2  d_trials,g_movecs,evl,ntrials,nproducts,algorithm,targetsym,
     3  symmetry,l_irs,k_irs,tda,lecut,ecut,lowin,owstart,owend,lewin,
     4  ewinl,ewinh,vshift)
c
c $Id$
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
      integer isym,jsym,maxireps
      parameter (maxireps=20)
      character*4 oname,oname_ir(maxireps)
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
      logical lsimpleso,status
      double precision evlso(nbf_ao,2)   ! spin-orbit eigenvalues
c
      integer kmin,kmax
      integer necut,newin
      double precision ecut,evl_diff,ewinl,ewinh
      integer vshift
      logical lecut, is_lower
      logical lowin, in_owin
      logical lewin, in_ewin
      integer owstart(2), owend(2), nowin
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)
c
      nodezero=(ga_nodeid().eq.0)
c     init irreps names
      do isym=1,maxireps
         oname_ir(isym)='jnk'
      enddo
c
c Determine the length of the trial vector: Occ*Virt
      call tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
c
c --------------------
c Allocate work arrays
c --------------------
c
C90000 if (nodezero.and.util_print('initial guess',print_default)
C     1  .and.(vector.ne.'')) write(LuOut,1000) 
C 1000 format(/,2x,'Unable to retrieve initial guess from a file')
c Consider nroots+20 smallest eigenvalue differences
c Extra 20 are there in case the nroots-th vector is degenarate
c In case of degenaracy, we must increase nroots
      if (.not.ma_push_get(mt_dbl,nov(1)+nov(2),
     1  'eigenvalue differences',l_ediff,k_ediff)) call errquit
     2  ('tddft_iguess: failed to allocate ediff',0, MA_ERR)
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
c      if (ga_nodeid().eq.0) then
c          write(6,*) "lecut: ",lecut
c          write(6,*) "vshift: ",vshift
c          write(6,*) "lewin: ",lewin
c          write(6,*) "lowin: ",lowin
c          write(6,*) "nov(1)+nov(2):",nov(1)+nov(2)
c          write(6,*) "nmo(1): ",nmo(1)
c          if (ipol.eq.2) write(6,*) "nmo(2): ",nmo(2)
c      end if
c
c     initialise
      l=0
      do i=1,ipol
        do j=nfc(i)+1,nocc(i)
          do k=nocc(i)+1,nmo(i)-nfv(i)
             l=l+1   
             dbl_mb(k_ediff+l-1)=0.d0
          end do
        end do
      end do
c
c Apply simplified spin-orbit correction to eigenvalues 
      lsimpleso = .false.
      status = rtdb_get(rtdb,'tddft:lsimpleso',mt_log,1,lsimpleso)
      status = status.and.lsimpleso
      if (status) then
         call tddft_simpleso(rtdb,ipol,nbf_ao)
         if(.not.rtdb_get(rtdb,'simpleso:evlso',mt_dbl,2*nbf_ao,evlso))
     $     call errquit('tddft_simpleso: rtdb get evlso ',914,rtdb_err)
c
c        Assign spin-orbit corrected eigenvalues
         do i = 1,ipol
            do j = 1,nbf_ao
             if (dabs(evlso(j,i)).gt.0.d0) evl(j,i) = evlso(j,i)
            end do ! j
         end do ! i
c
c         do i = 1,ipol
c            do j = 1,nbf_ao
c               if (ga_nodeid().eq.0) write(luout,*) j,evl(j,i)
c            end do ! j
c         end do ! i
      end if  ! status
c
c Calculate windows
      l=0
      necut = 0
      nowin = 0
      newin = 0
      do i=1,ipol
        kmin = nocc(i)+1+vshift
        kmax = nmo(i)-nfv(i)
        do j=nfc(i)+1,nocc(i)
          do k=kmin,kmax
           if (lecut) then
              is_lower = (evl(j,i).le.ecut)
              if (.not. is_lower) goto 100
              necut=necut+1
              l=necut
           else if (lowin) then
              in_owin = ((j.ge.owstart(i)).and.(j.le.owend(i))) 
              if (.not. in_owin) goto 100
              nowin=nowin+1
              l=nowin
           else if (lewin) then
              in_ewin = ((evl(j,i).le.ewinh .and. evl(j,i).ge.ewinl)) 
              if (.not. in_ewin) goto 100
              newin=newin+1
              l=newin
           else
              l=l+1
           endif  ! lecut,lowin,lewin
c
            evl_diff = evl(k,i)-evl(j,i)
            dbl_mb(k_ediff+l-1)=evl_diff

 100       continue
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
      oname_ir(isym)=oname
      if (nodezero.and.util_print('initial guess',print_debug))
     1  write(LuOut,*) 'Ground state symmetry ',oname
c
c ------------------------------------------
c Find nroot smallest eigenvalue differences
c ------------------------------------------
c
      if (lecut) then
       if (necut.eq.0) call errquit
     1  ('tddft_iguess: no excitations in energy window',necut,CALC_ERR)
       if (nroots.gt.necut) nroots = necut
      else if (lowin) then
       if (nowin.eq.0) call errquit
     1  ('tddft_iguess: no excitations in window',nowin,CALC_ERR)
       if (nroots.gt.nowin) nroots = nowin
      else if (lewin) then
       if (newin.eq.0) call errquit
     1  ('tddft_iguess: no excitations in window',newin,CALC_ERR)
       if (nroots.gt.newin) nroots = newin
      endif
c
      n=min(nroots+buf,nov(1)+nov(2))
c
      do m=1,n
        swap=huge
        l=0
        necut=0
        nowin=0
        newin=0
        found=.false.
        do i=1,ipol
         kmin = nocc(i)+1+vshift
         kmax = nmo(i)-nfv(i) 
         do j=nfc(i)+1,nocc(i)
          do k=kmin,kmax
             if (lecut) then
              is_lower = (evl(j,i).le.ecut)
              if (.not. is_lower) goto 200
              necut=necut+1
              l=necut
             else if (lowin) then
              in_owin = ((j.ge.owstart(i)).and.(j.le.owend(i))) 
              if (.not. in_owin) goto 200
              nowin=nowin+1
              l=nowin
             else if (lewin) then
              in_ewin = ((evl(j,i).le.ewinh .and. evl(j,i).ge.ewinl)) 
              if (.not. in_ewin) goto 200
              newin=newin+1
              l=newin
             else
              l=l+1
             endif  ! lecout,lowin,lewin
c
             jsym=ieor(isym-1,int_mb(k_irs(i)+j-1)-1)+1
             jsym=ieor(jsym-1,int_mb(k_irs(i)+k-1)-1)+1
             if(oname_ir(jsym).eq.'jnk') then
                call sym_irrepname(geom,jsym,oname)
                oname_ir(jsym)=oname
             else
                oname=oname_ir(jsym)
             endif
c
             if (((.not.symmetry).or.(targetsym.eq.oname))
     1             .and.(dbl_mb(k_ediff+l-1).lt.swap)) then
                   int_mb(k_sorto+m-1)=j
                   int_mb(k_sortv+m-1)=k
                   int_mb(k_sortp+m-1)=i
                   int_mb(k_sorta+m-1)=l
                   int_mb(k_irrep+m-1)=jsym
                   swap=dbl_mb(k_ediff+l-1)
                   found=.true.
             endif
c
 200       continue
          enddo
         enddo
        enddo
c
        if (.not.found) call errquit
     1    ('tddft_iguess: too many roots requested',nroots, INPUT_ERR)
        dbl_mb(k_ediff+int_mb(k_sorta+m-1)-1)=huge+huge
      enddo  ! loop over n
c
c ---------------------------------------------------
c Print the smallest eigenvalue differences for debug
c ---------------------------------------------------
c
      if (nodezero) then
        write(LuOut,*)
        write(LuOut,9000) n
        write(LuOut,9020)
        write(LuOut,9010)
        write(LuOut,9020)
        do m=1,n
           jsym=int_mb(k_irrep+m-1)
           if(oname_ir(jsym).eq.'jnk') then
                call sym_irrepname(geom,jsym,oname)
                oname_ir(jsym)=oname
             else
                oname=oname_ir(jsym)
             endif
          write(LuOut,9030) m,int_mb(k_sortp+m-1),
     1      int_mb(k_sorto+m-1),int_mb(k_sortv+m-1),oname,
     2      evl(int_mb(k_sorto+m-1),int_mb(k_sortp+m-1)),
     3      evl(int_mb(k_sortv+m-1),int_mb(k_sortp+m-1)),
     4      (evl(int_mb(k_sortv+m-1),int_mb(k_sortp+m-1))-
     5      evl(int_mb(k_sorto+m-1),int_mb(k_sortp+m-1)))*au2ev
        enddo
        write(LuOut,9020)
        call util_flush(LuOut)
      endif
 9000 format(i5,' smallest eigenvalue differences (eV) ')
 9010 format('  No. Spin  Occ  Vir  Irrep   E(Occ)    E(Vir)   E(Diff)')
 9020 format('--------------------------------------------------------')
 9030 format(4i5,1x,a4,1x,3f10.3)
c
c 9030 format(4i5,1x,a4,1x,2f10.5,f18.2)
c
c --------------------------------------------------------
c Check if nroots is too large, if so, reduce it
c Check if the nroots-th excited state could be degenerate
c If so, increase nroots iteratively
c --------------------------------------------------------
c
      if (lecut) then
       if (nroots.gt.necut) nroots = necut
      else if (lowin) then
       if (nroots.gt.nowin) nroots = nowin
      else if (lewin) then
       if (nroots.gt.newin) nroots = newin
      else if (nroots.gt.nov(1)+nov(2)) then
        if (nodezero.and.util_print('initial guess',print_low))
     1    write(LuOut,9040) nroots,nov(1)+nov(2)
        nroots=nov(1)+nov(2)
        if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1    call errquit('tddft_iguess: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (nroots.lt.nov(1)+nov(2)) then
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
            if (nov(i).gt.0) then
               if(nodezero)
     Z         call ga_put(g_trials(i),1,nov(i),m,m,dbl_mb(k_t(i)),
     L              nov(i))
            endif
          else
               if(nodezero)
     Z       call ga_put(g_buffer1(i),1,nov(i),1,1,dbl_mb(k_t(i)),
     L            nov(i))
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
c
      return
      end
