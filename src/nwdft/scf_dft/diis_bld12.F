      subroutine odiis_bld12(toll_s, svals, iga_svecs, iga_sout, 
     &                      iga_stmp,  iwhich)
c***********************************************************************
c     sbuild:  Build S, S**(-1/2) or S**(+1/2) from S eigenvalues and
c     S eigenvectors.
c     iwhich = 1; build S
c     iwhich = 2; build S**(-1/2)
c     iwhich = 3; build S**(+1/2)
c***********************************************************************
C$Id: diis_bld12.F,v 1.12 2000-03-16 00:06:27 edo Exp $
c
      implicit none
c
      double precision svals(*) !  S evals [input]
      integer iga_svecs ! GA handle for S evecs [input]
      integer iga_sout  ! GA handle for S^(n) [ouput]
      integer iga_stmp  ! GA handle for scratch [input]
      integer iwhich    ! 2->S(-1/2) 3->S(1/2) 
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
c     
      logical LResult
      integer me, nproc, i
      integer ltmpm, itmpm
      double precision toll_s, aaz
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      if (.not. MA_Push_Get(MT_Dbl, nbf_ao, 'tmpm', ltmpm, itmpm))
     &   call errquit('diis_bld12: failed to alloc tmpm',0)
      call ga_print(iga_svecs)
c
cedo      call ga_copy(iga_svecs, iga_stmp)
c      
c        Build S.
c     
         call ga_zero(iga_stmp)
         do i = me+1, nbf_ao, nproc
            aaz = svals(i)
            if (aaz.ge.toll_s)  then
               if (iwhich.eq.2)then
                  aaz = 1.d0/sqrt(aaz)
               else
                  aaz = dsqrt(aaz)
               endif
               call get_col(iga_svecs, nbf_ao, i, DBL_MB(itmpm))
               call dscal(nbf_ao, aaz, DBL_MB(itmpm), 1)
               call put_col(iga_stmp, nbf_ao, i, DBL_MB(itmpm))
!wrong               call ga_scale_patch(iga_stmp,1,nbf_ao,i,i,aaz)
            else
c               write(LuOut,*)' Small eig of S in build of S',
c     &           i, aaz, toll_s
               call get_col(iga_svecs, nbf_ao, i, DBL_MB(itmpm))
               call put_col(iga_stmp, nbf_ao, i, DBL_MB(itmpm))
            endif
         enddo
c
c
      LResult = MA_Pop_Stack(ltmpm)
      call ga_print(iga_stmp)
c
cedo      call ga_sync
c
      call ga_dgemm('N', 'T', nbf_ao, nbf_ao, nbf_ao, 1.d0, 
     &              iga_stmp, iga_svecs, 0.d0, iga_sout)
c
      call ga_print(iga_sout)
      return
      end
      subroutine diis_bld12(toll_s, svals, iga_svecs, iga_sout, 
     &                      iga_stmp, iw)
c***********************************************************************
c     sbuild:  Build S, S**(-1/2) or S**(+1/2) from S eigenvalues and
c     S eigenvectors.
c     iw = 1; build S
c     iw = 2; build S**(-1/2)
c     iw = 3; build S**(+1/2)
c***********************************************************************
C$Id: diis_bld12.F,v 1.12 2000-03-16 00:06:27 edo Exp $
c
      implicit none
c
      double precision svals(*) !  S evals [input]
      integer iga_svecs ! GA handle for S evecs [input]
      integer iga_sout  ! GA handle for S^(n) [ouput]
      integer iga_stmp  ! GA handle for scratch [input]
      integer iw    ! 2->S(-1/2) 3->S(1/2) 
c
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
#include "bas.fh"
c     
      integer me, nproc
      double precision toll_s
c
      integer nxtask
      external nxtask
      integer i_a,i_b,l_a,l_b,icount,ilo,ihi,jlo,jhi,
     ,     k_down,l_down,n2,n1,nn2,k_up,l_up
      logical diis_scal1,diis_scal2
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c 
!      call dft_tstamp(' before bld')
      icount=min(ncenters,2)
      call ga_sync()
      n1 = 0
      n2 = nxtask(nproc,icount)

      nn2=nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_dbl,nn2,'ablock',l_down,k_down))
     &     call errquit('diis_bld12: cant allocate atblock',33)
      if (.not.MA_Push_Get(MT_dbl,nn2,'ablock',l_up,k_up))
     &     call errquit('diis_bld12: cant allocate atblock',44)
c      
c        Build S.
c     
      call ga_zero(iga_stmp)

      do i_a=1,ncenters
         if (.not. bas_ce2bfr(ao_bas_han, i_a, ilo, ihi))
     $        call errquit('diis_bld12: ce2bfr failed', 0)
         do i_b=1,i_a
            if (n1.eq.n2 )then
               call get_atom_block(iga_svecs, ao_bas_han,
     &              i_a, i_b, dbl_mb(k_down),l_a, l_b)
               if(i_a.ne.i_b)  then
                  call get_atom_block(iga_svecs, ao_bas_han,
     &                 i_b, i_a, dbl_mb(k_up), l_b, l_a)
                  if (.not. bas_ce2bfr(ao_bas_han, i_b, jlo, jhi)) 
     $                 call errquit('get_atom_block: ce2bfr failed',1)
               endif
               if(i_a.eq.i_b) then
                  if(diis_scal1(iw,svals(ilo),dbl_mb(k_down),
     ,                 l_a,toll_s))
     .                 call upd_atom_block(iga_stmp,ao_bas_han, 
     .                 i_a, i_a, dbl_mb(k_down))
               else
                  if(diis_scal2(iw,svals(jlo),svals(ilo),dbl_mb(k_down),
     .                 dbl_mb(k_up),l_a,l_b,toll_s)) then
                     call upd_atom_block(iga_stmp,ao_bas_han, 
     .                    i_a, i_b, dbl_mb(k_down))
                     call upd_atom_block(iga_stmp,ao_bas_han, 
     .                    i_b, i_a, dbl_mb(k_up))
                  endif
               endif


               n1=n1+1
               n2 = nxtask(nproc,icount)
            else
               n1=n1+1
            endif
         enddo
      enddo
      if (.not.ma_pop_stack(l_up))
     &   call errquit('diis_bld12: cannot pop stack',0)
      if (.not.ma_pop_stack(l_down))
     &   call errquit('diis_bld12: cannot pop stack',1)
      n1 = nxtask(-nproc,icount)

!      call dft_tstamp(' done bld')
      call ga_dgemm('N', 'T', nbf, nbf, nbf, 1.d0, 
     &              iga_stmp, iga_svecs, 0.d0, iga_sout)
!      call dft_tstamp(' done dgemm')
c
      return
      end
      logical function diis_scal1(iw,value,down,l_a,toll_s)
      implicit none
      integer l_a,iw
      double precision down(l_a,l_a),value(l_a),toll_s
c
      integer i
      double precision sum,thresh,azz
      parameter(thresh=1d-14)
      double precision dasum
      external dasum
c
      diis_scal1=.true.
      sum=0d0
      do i=1,l_a
         azz=value(i)
         if (azz.ge.toll_s)  then
            if (iw.eq.2)then
               azz = 1.d0/sqrt(azz)
            elseif(iw.eq.3) then
               azz = sqrt(azz)
            endif
            
            call dscal(l_a,azz,down(1,i),1)
         endif
         sum=dasum(l_a,down(1,i),1)
      enddo
      diis_scal1=sum.gt.thresh
      return
      end
      logical function diis_scal2(iw,valued,valueu,
     ,     down,up,l_a,l_b,toll_s)
      implicit none
      integer l_a,l_b,iw
      double precision down(l_a,l_b),up(l_b,l_a),
     ,     valued(l_b),valueu(l_a),toll_s
c
      integer j
      double precision sum,thresh,azz
      parameter(thresh=1d-12)
      double precision dasum
      external dasum
c
      diis_scal2=.true.
      sum=0d0
      do j=1,l_b
         azz=valued(j)
         if (azz.ge.toll_s)  then
            if (iw.eq.2)then
               azz = 1.d0/sqrt(azz)
            elseif(iw.eq.3) then
               azz = sqrt(azz)
            endif
            call dscal(l_a,azz,down(1,j),1)
         endif
         sum=dasum(l_a,down(1,j),1)
      enddo
      do j=1,l_a
         azz=valueu(j)
         if (azz.ge.toll_s)  then
            if (iw.eq.2)then
               azz = 1.d0/sqrt(azz)
            elseif(iw.eq.3) then
               azz = sqrt(azz)
            endif
            call dscal(l_b,azz,up(1,j),1)
         endif
         sum=dasum(l_b,up(1,j),1)
      enddo
      diis_scal2=sum.gt.thresh
      return
      end
      subroutine diis_transp(a,b,m,n)
      implicit none
      integer m,n
      double precision a(m,n),b(m,n)
c
      integer i,j
c
      do i=1,m
         do j=1,n
            b(j,i)=a(i,j)
         enddo
      enddo
      return
      end
      
      
