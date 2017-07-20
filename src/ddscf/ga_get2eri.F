      subroutine ga_get2eri(basis, g_2ceri,oskel,maxg_in,mscratch_in, 
     &            int_call)
c$Id$
      implicit none
#include "errquit.fh"
      integer basis
      integer g_2ceri
c
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "sym.fh"
c      
      integer  maxg_in, mscratch_in
c
      integer  me, nproc, icount,
     &        maxg, mscratch_2e2c, Nintegrals,
     &        ishp, ifirstp, ilastp, nshp,
     &        ishq, ifirstq, ilastq, nshq
      integer lbuf, ibuf, lscr, iscr
      integer geom
      logical oskel,odoit
      integer ishq1,ishq2,iat,jat,ishp1,ishp2
      integer iptr,in0,in1,jn0,jn1,nbfa,nbfb
      integer natoms 
      integer vecl
      double precision q2
      integer i0,i1,j0,j1,ilo,ihi
      integer iat1,iat2,jat1,jat2
      integer adrc,ldc,j,g_2
      integer mem_lef
      logical util_ga2atoms
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,util_ga2atoms
      external int_call
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('ga_get2eri: basis corrupt?', 0, BASIS_ERR)
      if (.not. geom_ncent(geom, natoms))
     &     call errquit('ga_get2eri: geom_ncent failed',73, GEOM_ERR)
      g_2 = ga_create_atom_blocked(geom, basis,'g2g2')

      if (util_ga2atoms(g_2, natoms, basis,
     A     iat1,iat2,jat1,jat2)) then

      call ga_distribution(g_2,
     .     ga_nodeid(), i0, i1, j0, j1)
      call ga_access(g_2, i0, i1, j0, j1, adrc, ldc)
c     
c     Compute the matrix of 2-ctr ERIs.
c     
c
c     Allocate scratch space.
c
      me  = ga_nodeid()
      nproc  = ga_nnodes()
      mscratch_2e2c = mscratch_in
      maxg = maxg_in
c
      vecl=0
c     get max no. shells on a given atom
      do iat=iat1,iat2
         if (.not. bas_ce2cnr(basis, iat, ishp1, ishp2))
     &        call errquit('ga_get2eri: bas_cn2bfr',0,
     &        BASIS_ERR)
         vecl=max(vecl,ishp2-ishp1+1)
      enddo
      do iat=jat1,jat2
         if (.not. bas_ce2cnr(basis, iat, ishp1, ishp2))
     &        call errquit('ga_get2eri: bas_cn2bfr',0,
     &        BASIS_ERR)
         vecl=max(vecl,ishp2-ishp1+1)
      enddo
      if (.not.MA_Push_Get(MT_Dbl,mscratch_2e2c,'scr',lscr,iscr))
     &   call errquit('ga_get2eri: cannot allocate scr',2, MA_ERR)
cold      vecl=vecl*vecl*maxg
      vecl=vecl*maxg
      if (.not.MA_Push_Get(MT_Dbl,vecl,'buf',lbuf,ibuf))
     &   call errquit('ga_get2eri: cannot allocate buf',0, MA_ERR)
c
c
      do iat=iat1,iat2
         if (.not. bas_ce2cnr(basis, iat, ishp1, ishp2))
     &        call errquit('ga_get2eri: bas_cn2bfr',0,
     &        BASIS_ERR)
         if (.not. bas_ce2bfr(basis, iat, in0, in1))
     &        call errquit('ga_get2eri:  bas_cn2bfr',0,
     &        BASIS_ERR)
         nbfa=in1-in0+1
         do jat=jat1,jat2
            if (.not. bas_ce2cnr(basis, jat, ishq1, ishq2))
     &           call errquit('ga_get2eri:bas_cn2bfr',0,
     &           BASIS_ERR)
            if (.not. bas_ce2bfr(basis, jat, jn0, jn1))
     &           call errquit('ga_get2eri: bas_cn2bfr',0,
     &           BASIS_ERR)
            nbfb=jn1-jn0+1
            
            do  ishq = ishq1, ishq2
               if (.not. bas_cn2bfr(basis, ishq, 
     Q                 ifirstq, ilastq))
     &              call errquit('ga_get2eri:bas_cn2bfr',0,
     &              BASIS_ERR)
c     zero because of symm
!            if(oskel) 
               nshq = ilastq - ifirstq + 1
               call dcopy(nbfa*nshq, 0.0d0,0, dbl_mb(ibuf), 1)
               iptr=ibuf
               do  ishp = ishp1,ishp2
c     
                  if (.not. bas_cn2bfr(basis, 
     I                 ishp, ifirstp, ilastp))
     &                 call errquit('ga_get2eri:_cn2bfr',0,
     &                 BASIS_ERR)
c     
                  nshp = ilastp - ifirstp + 1
c     
                  odoit=.true.
                  q2=1d0
                  if (oskel.and.iat.ne.jat) odoit=
     =                 sym_shell_pair(basis,ishq,ishp,q2)
                  if (odoit) then
                     Nintegrals=maxg
                     if(nintegrals+iptr-ibuf.gt.vecl) then
                     write(0,*) ga_nodeid(),
     ,                    'I0',iptr+nintegrals-ibuf,vecl
                      write(0,*) ga_nodeid(),
     ,                    'NN',nintegrals
                     call errquit('ga_get2eri:dim busted',0,0)
                     endif 
c     
                     call int_call(basis, ishq, basis, 
     I                    ishp,
     &                    mscratch_2e2c, Dbl_MB(iscr), Nintegrals,
     &                    Dbl_MB(iptr))
                     if(q2.gt.1d0) call dscal(
     .                       nshp*nshq,q2,dbl_mb(iptr),1)
    
                  endif
                  iptr=iptr+nshq*nshp
c     
               enddo            ! i shell loop
               call dscat_arr2(basis, i0,i1,j0,j1,iat,jat,ishq,
     V              dbl_mb(ibuf),dbl_mb(adrc))
c     
            enddo               ! j shell loop
         enddo                  ! i loop
      enddo                     ! j loop
c     
      if (.not.ma_chop_stack(lscr)) then
         if(ga_nodeid().eq.0) call ma_summarize_allocated_blocks()
         call ga_sync
         call errquit('ga_get2eri: cannot chop stack',0, MA_ERR)
      endif
      call ga_release_update(g_2, i0, i1, j0, j1)
      endif
      call ga_copy(g_2,g_2ceri)
      if (.not. ga_destroy(g_2)) call errquit('ga_get2eri: gad?',
     &     0, 0)
c
      return
      end

      subroutine dscat_arr(basis_hand,x0,x1,y0,y1,ata,atb,
     V     vecin,vecout)
      implicit none
#include "bas.fh"
#include "errquit.fh"
      integer basis_hand,ata,atb
      integer x0,y0,x1,y1
      double precision vecin(*),vecout(x0:x1,y0:y1)
c
      integer sha0,sha1,shb0,shb1
      integer sha,shb,a0,a1,b0,b1,a,b
      integer iptr,idima,idimb,idimx
      idimx=x1-x0+1
c
      if (.not. bas_ce2cnr( basis_hand, ata, sha0, sha1)) then
        call errquit('Exiting in scatlab_mat',3, BASIS_ERR) 
      endif 
      if (.not. bas_ce2cnr( basis_hand, atb, shb0, shb1)) then
        call errquit('Exiting in scatlab_mat',4, BASIS_ERR)
      endif 
c
c     rearrange fron integr. order  to matrix order
c
      iptr=1
      do shb=shb0,shb1
         if (.not. bas_cn2bfr(basis_hand, shb,b0,b1))
     &        call errquit('ga_get2eri:_cn2bfr',0, BASIS_ERR)
         idimb=b1-b0+1
         do sha=sha0,sha1
            if (.not. bas_cn2bfr(basis_hand, sha,a0,a1))
     &           call errquit('ga_get2eri:_cn2bfr',0, BASIS_ERR)
            idima=a1-a0+1
            do b=b0,b1
               do a=a0,a1
                  vecout(a,b)=vecin(iptr)
                  iptr=iptr+1
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine dscat_arr2(basis_hand,x0,x1,y0,y1,ata,atb,shb,
     V     vecin,vecout)
      implicit none
#include "bas.fh"
#include "errquit.fh"
      integer basis_hand,ata,atb
      integer shb ! [in]
      integer x0,y0,x1,y1
      double precision vecin(*),vecout(x0:x1,y0:y1)
c
      integer sha0,sha1
      integer sha,a0,a1,b0,b1,a,b
      integer iptr,idima,idimb,idimx
      idimx=x1-x0+1
c
      if (.not. bas_ce2cnr( basis_hand, ata, sha0, sha1)) then
        call errquit('Exiting in scatlab_mat',3, BASIS_ERR) 
      endif 
c
c     rearrange fron integr. order  to matrix order
c
      iptr=1
      if (.not. bas_cn2bfr(basis_hand, shb,b0,b1))
     &     call errquit('ga_get2eri:_cn2bfr',0, BASIS_ERR)
      idimb=b1-b0+1
      do sha=sha0,sha1
         if (.not. bas_cn2bfr(basis_hand, sha,a0,a1))
     &        call errquit('ga_get2eri:_cn2bfr',0, BASIS_ERR)
         idima=a1-a0+1
         do b=b0,b1
            do a=a0,a1
               vecout(a,b)=vecin(iptr)
               iptr=iptr+1
            enddo
         enddo
      enddo
      return
      end
c $Id$
      logical function util_ga2atoms(g_a, natoms, basis,
     A     iat1,iat2,jat1,jat2)
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "bas.fh"
      integer g_a,natoms, basis ![in]
      integer iat1,iat2,jat1,jat2 ![out]
c
      integer i0,i1, j0, j1,iat,ilo,ihi
c
      call ga_distribution(g_a,
     .     ga_nodeid(), i0, i1, j0, j1)
      if (i0.gt.0 .and. i0.le.i1) then
c
c     get atoms list
c
         iat1=-99
         iat2=-99
         jat1=-99
         jat2=-99
         do iat=1,natoms
            if (.not. bas_ce2bfr(basis, iat, ilo, ihi)) call errquit
     $           ('ga_get2eri bas_ce2bfr ?', 0, BASIS_ERR)
            if(ilo.eq.i0) iat1=iat
            if(ilo.eq.j0) jat1=iat
            if(ihi.eq.i1) iat2=iat
            if(ihi.eq.j1) jat2=iat
         enddo
         if(iat1.eq.-99.or.iat2.eq.-99.or.
     .        jat1.eq.-99.or.jat2.eq.-99) then
            write(6,*) ' distr ',i0,i1,j0,j1
            write(6,*) ' iats ',iat1,iat2,jat1,jat2
            call errquit
     A           ('gaget2eri: iat not found ',0,0)
         endif
         util_ga2atoms=.true.
      else
         util_ga2atoms=.false.
      endif
      
      return
      end
