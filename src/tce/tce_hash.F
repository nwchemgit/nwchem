      subroutine tce_hash_v2(hash,key,offset)
c
c $Id$
c
      implicit none
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer hash(*)
      integer key
      integer offset
      integer length,length1
      integer i,j,k,l
      integer i_start,j_start,k_start,l_start
      integer i_stop,j_stop,k_stop,l_stop
      integer pos1,pos2
      integer pos1_l,pos1_u
      integer pos2_l,pos2_u
      integer middle
      integer key_loop
      integer irow,icol,index_pair,index_point
      integer sa_ij,sa_kl
      integer ieo_kl,n
      integer indx_pairij,indx_point
      indx_point(i,j,n)=
     =     (i*(2*n+1-i))/2 -n + j
      indx_pairij(i,j)=j+(i*(i-1))/2
c
      length = hash(1)
      middle=-1
      do i=1,length
        if((hash(i+1).le.key).and.(key.le.hash(i+2))) then
            middle=i
            go to 20
        end if
      end do 
      if(middle.eq.-1) 
     &      call errquit('hash: key not found 1',0,UNKNOWN_ERR) 
   20 continue
c middle is fixed
c
      i_start=hash(2*(length+1)+middle+1)
      j_start=hash(3*(length+1)+middle+1)
      k_start=hash(4*(length+1)+middle+1)
      l_start=hash(5*(length+1)+middle+1)
c
      i_stop=hash(2*(length+1)+(middle+1)+1)
      j_stop=hash(3*(length+1)+(middle+1)+1)
      k_stop=hash(4*(length+1)+(middle+1)+1)
      l_stop=hash(5*(length+1)+(middle+1)+1)
c
      offset=hash((length+1)+middle+1)
c
      pos1_l=indx_point(i_start,j_start,noa+nva)
      pos1_u=indx_point(i_stop,j_stop,noa+nva)
c
      pos2_l=indx_point(k_start,l_start,noa+nva)
      pos2_u=indx_point(k_stop,l_stop,noa+nva)
c
      i=i_start
      j=j_start
      k=k_start
      l=l_start
c
      do pos1=pos1_l,pos1_u  !---------------------------------------
      ieo_kl=ieor(int_mb(k_sym_alpha+k-1),int_mb(k_sym_alpha+l-1))
      sa_kl=int_mb(k_spin_alpha+k-1)+int_mb(k_spin_alpha+l-1)
#if 0
       if(pos1.ne.pos1_l) then
        if(pos1.eq.indx_point(i+1,i+1,noa+nva)) then
          i=i+1
          j=i
        else
          j=j+1
        end if
       end if
#endif
       sa_ij=int_mb(k_spin_alpha+i-1)+int_mb(k_spin_alpha+j-1)
c
  100  continue
c
          IF (sa_ij.eq.sa_kl) then
          IF (ieor(int_mb(k_sym_alpha+i-1),ieor(int_mb(k_sym_alpha+j-1),
     &    ieo_kl)) .eq.
     &    irrep_v) THEN
          irow=indx_pairij(j,i)
          icol=indx_pairij(l,k)
          IF(irow.ge.icol) THEN
           key_loop= l - 1 + (noa+nva) * (k - 1 + (noa+nva) * (j-
     &               1 + (noa+nva) * (i - 1)))
           if(key_loop.eq.key) go to 103
           offset=offset+int_mb(k_range_alpha+i-1)*
     &            int_mb(k_range_alpha+j-1)*int_mb(k_range_alpha+k-1)*
     &            int_mb(k_range_alpha+l-1)
          end if
          end if
          end if
c if i eq j eq k eq l eq - wyskok
          if((i.eq.i_stop).and.(j.eq.j_stop).and.(k.eq.k_stop).and.
     &       (l.eq.l_stop)) go to 102
          if((k.eq.noa+nva).and.(l.eq.noa+nva)) go to 200
          if(l.eq.noa+nva) then
            k=k+1
            l=k
          else
            l=l+1
          end if
          sa_kl=int_mb(k_spin_alpha+k-1)+int_mb(k_spin_alpha+l-1)
          ieo_kl=ieor(int_mb(k_sym_alpha+k-1),int_mb(k_sym_alpha+l-1))
          go to 100
c 
  200 continue
      k=1
      l=1
          sa_kl=int_mb(k_spin_alpha+k-1)+int_mb(k_spin_alpha+l-1)
          ieo_kl=ieor(int_mb(k_sym_alpha+k-1),int_mb(k_sym_alpha+l-1))
        if(pos1+1.eq.indx_point(i+1,i+1,noa+nva)) then
          i=i+1
          j=i
        else
          j=j+1
        end if
      enddo                 !---- pos1 -------------------------------
c
c
  102 continue
      write(*,*)'key=',key
      call errquit('hashv2: key not found 2',0,UNKNOWN_ERR)
  103 continue
      return
      end
c
c
c
c
      subroutine tce_hash_v2_act(hash,key,offset)
c
c $Id$
c
      implicit none
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer hash(*)
      integer key
      integer offset
      integer length,length1
      integer i,j,k,l
      integer i_start,j_start,k_start,l_start
      integer i_stop,j_stop,k_stop,l_stop
      integer pos1,pos2
      integer pos1_l,pos1_u
      integer pos2_l,pos2_u
      integer middle
      integer key_loop
      integer irow,icol,index_pair,index_point
c - exascale ---
      LOGICAL is_active_4_o
c --------------
c
      length = hash(1)
      middle=-1
      do i=1,length
        if((hash(i+1).le.key).and.(key.le.hash(i+2))) then
            middle=i
            go to 20
        end if
      end do 
      if(middle.eq.-1) 
     &      call errquit('hash: key not found 1',0,UNKNOWN_ERR) 
   20 continue
c middle is fixed
c
      i_start=hash(2*(length+1)+middle+1)
      j_start=hash(3*(length+1)+middle+1)
      k_start=hash(4*(length+1)+middle+1)
      l_start=hash(5*(length+1)+middle+1)
c
      i_stop=hash(2*(length+1)+(middle+1)+1)
      j_stop=hash(3*(length+1)+(middle+1)+1)
      k_stop=hash(4*(length+1)+(middle+1)+1)
      l_stop=hash(5*(length+1)+(middle+1)+1)
c
      offset=hash((length+1)+middle+1)
c
      pos1_l=index_point(i_start,j_start,noa+nva)
      pos1_u=index_point(i_stop,j_stop,noa+nva)
c
      pos2_l=index_point(k_start,l_start,noa+nva)
      pos2_u=index_point(k_stop,l_stop,noa+nva)
c
      i=i_start
      j=j_start
      k=k_start
      l=l_start
c
      do pos1=pos1_l,pos1_u  !---------------------------------------
       if(pos1.ne.pos1_l) then
        if(pos1.eq.index_point(i+1,i+1,noa+nva)) then
          i=i+1
          j=i
        else
          j=j+1
        end if
       end if
c
  100  continue
c
c - exascale ---
      IF(.not.((i.gt.noa).and.(j.gt.noa).and.(k.gt.noa).and.
     &   (l.gt.noa).and.(.not.is_active_4_o(i,j,k,l)))) THEN
c ---
          IF (int_mb(k_spin_alpha+i-1)+int_mb(k_spin_alpha+j-1).eq.
     &    int_mb(k_spin_alpha+k-1)+int_mb(k_spin_alpha+l-1)) THEN
          IF (ieor(int_mb(k_sym_alpha+i-1),ieor(int_mb(k_sym_alpha+j-1),
     &    ieor(int_mb(k_sym_alpha+k-1),int_mb(k_sym_alpha+l-1)))) .eq.
     &    irrep_v) THEN
          irow=index_pair(j,i)
          icol=index_pair(l,k)
          IF(irow.ge.icol) THEN
           key_loop= l - 1 + (noa+nva) * (k - 1 + (noa+nva) * (j-
     &               1 + (noa+nva) * (i - 1)))
           if(key_loop.eq.key) go to 103
           offset=offset+int_mb(k_range_alpha+i-1)*
     &            int_mb(k_range_alpha+j-1)*int_mb(k_range_alpha+k-1)*
     &            int_mb(k_range_alpha+l-1)
          end if
          end if
          end if
c - exascale ---
      END IF
c ---
c if i eq j eq k eq l eq - wyskok
          if((i.eq.i_stop).and.(j.eq.j_stop).and.(k.eq.k_stop).and.
     &       (l.eq.l_stop)) go to 102
          if((k.eq.noa+nva).and.(l.eq.noa+nva)) go to 200
          if(l.eq.noa+nva) then
            k=k+1
            l=k
          else
            l=l+1
          end if
          go to 100
c
  200 continue
      k=1
      l=1
      enddo                 !---- pos1 -------------------------------
c
c
c
c
  102 continue
      write(*,*)'key=',key
      call errquit('hashv2act: key not found 2',0,UNKNOWN_ERR)
  103 continue
      return
      end
c
c
c
c
c
c
      subroutine tce_hash(hash,key,offset)
c
c $Id$
c
      implicit none
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer hash(*)
      integer key
      integer offset
      integer length
      integer i
      integer less,middle,more
c
      length = hash(1)
      less = 1
      more = length
   10 continue
        if (more-less.le.4) then
          middle = -1
          do i = less, more
            if (hash(i+1).eq.key) middle = i
          enddo
          if (middle.eq.-1) then
            write(*,*) 'key = ',key
            do i = 1,hash(1)
              write(*,*) i,hash(i+1)
            enddo
            call errquit('tce_hash: key not found',0,UNKNOWN_ERR)
          endif
          goto 20
        endif
        middle = (less + more) / 2
        if (hash(middle+1).eq.key) then
          goto 20
        else if (hash(middle+1).gt.key) then
          more = middle
        else
          less = middle
        endif
        goto 10
   20 continue
      offset = hash(length+middle+1)
      return
      end
c
c
c
c
c
c
c
      subroutine tce_hash_n(hash,key,offset)
c
c $Id$
c
      implicit none
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer hash(*)
      integer key
      integer offset
      integer length
      integer i
c
cedo      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
      length = hash(1)
cedo      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
      if(length.lt.0) then
         write(6,*) key,hash(1),hash(2)
         call errquit('tce_hash_n: negative length',length,UNKNOWN_ERR)
      endif
       offset = -1
       do i = 1,length
        if(key.eq.hash(i+1)) then
         offset=hash(length+i+1)
         goto 10
        end if
       end do
   10  continue
cedo        if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
       if(offset.eq.-1) then
#if 1
      write(6,123) ' hashn: addr ',i,' key ',key
      write(6,*) ' length ',length
      do i=1,length
         write(6,*) ' i ',i,' hash ',hash(i+1)
      enddo
 123  format(a,i12,a,i12)
#endif
         call errquit('tce_hash_n: key not found',key,UNKNOWN_ERR)
       end if
c
      return
      end
