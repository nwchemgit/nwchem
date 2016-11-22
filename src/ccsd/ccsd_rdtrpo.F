      subroutine ccsd_rdtrpo(t1,buf1,buf2,g_objo,g_objv,
     $     nocc,nvir,iprt)
!     $Id$
      implicit none
#include "errquit.fh"
!     
#include "ccsd_len.fh"
      integer g_objo,g_objv,nocc,nvir,iprt
      double precision t1(nvir,nocc),buf1(lnobj),buf2(lnobj)
!     
      integer i,j,k,a,av,bv,ad1,ad2,lbfo
      Integer Nodes, IAm
      integer g_tmp, start, end
!     
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "ccsdps.fh"
!     
      Nodes = GA_NNodes()
      IAm = GA_NodeID()

      if (occsdps) then
         call pstat_on(ps_rdtrpo)
      else
         call qenter('rdtrpo',0)
      endif

      call ga_sync
!     
!     - read in objects
!     - [io|ov], [oo|iv], t2(io:vv) / [ao|vo], t2(oo:av), t2(oo,va)
!     
      lbfo=lnovv+lnoov+lnoov
!     
!     do av ... uses 2*oov+ooo memory.  Do not use locality
!     so that we can distribute more uniformly across nodes.
!     Locality is not used in the triples driver.
!     
      do av=iam+1,nvir,nodes 
!     
!     get dint and eint
         call ga_get(g_objv,1,2*lnoov+lnooo,av,av,buf1,
     &        2*lnoov+lnooo)
!     
!     process eints for occupied indexed buffer
         do i=1,nocc
!     eintc
            ad1=0
            do j=1,nocc
               do k=1,nocc
                  ad1=ad1+1
                  buf2(ad1)=buf1(lnoov+(j-1)*lnoo+(i-1)*nocc+k)
               end do
            end do
            ad2=lnovv+(av-1)*lnoo
            call ga_put(g_objo,ad2+1,ad2+lnoo,i,i,buf2,lnoo)
!     eintx
            ad1=lnoov+(i-1)*lnoo+1
            ad2=lnovv+lnoov+(av-1)*lnoo
            call ga_put(g_objo,ad2+1,ad2+lnoo,i,i,buf1(ad1),lnoo)
         end do
!     
!     process dints for virtual indexed buffer
         ad1=0
         do i=1,nocc
            do j=1,nocc
               do bv=1,nvir
                  ad1=ad1+1
                  ad2=(bv-1)*lnoo+(i-1)*nocc+j
                  buf2(ad1)=buf1(ad2)
               end do
            end do
         end do
!     
!     process t2 amplitudes for virtual indexed buffer
!     t2c
         do i=1,nocc
            do j=1,nocc
               do bv=1,nvir
                  ad1=ad1+1
                  ad2=lnoov+lnooo+(bv-1)*lnoo+(i-1)*nocc+j
                  buf2(ad1)=buf1(ad2)
               end do
            end do
         end do
!     t2x
         do i=1,nocc
            do j=1,nocc
               do bv=1,nvir
                  ad1=ad1+1
                  ad2=lnoov+lnooo+(bv-1)*lnoo+(j-1)*nocc+i
                  buf2(ad1)=buf1(ad2)
               end do
            end do
         end do
         if (iprt.gt.50)then
            print *,'t2c new ',iam,av
            write(6,4859)(buf2(j),j=lnoov+1,2*lnoov)
            print *,'t2x new ',iam,av
            write(6,4859)(buf2(j),j=2*lnoov+1,3*lnoov)
            print *,'dint new ',iam,av
            write(6,4859)(buf2(j),j=1,lnoov)
 4859       format(1x,5e14.4)
         end if
         call ga_put(g_objv,1,ad1,av,av,buf2,ad1)
!     
!     process t2 amplitudes for occupied indexed buffer
!     t2j
         do i=1,nocc
            ad1=0
            do j=1,nocc
               do bv=1,nvir
                  ad1=ad1+1
                  buf2(ad1)=
     $                 buf1(lnoov+lnooo+(bv-1)*lnoo+(i-1)*nocc+j)
               end do
            end do
            ad2=(av-1)*lnov
            call ga_put(g_objo,ad2+1,ad2+lnov,i,i,buf2,lnov)
         end do
!     
      end do
!     
      call ga_sync()            ! ga_sync ensures all GA ops complete
!     
!     sort t1 then broadcast
!     
      if (iam.eq.0)then
         call ga_get(g_objv,lnoov+lnooo+1,lnoov+lnooo+lnov,
     $        nvir+1,nvir+1, buf1,lnov)
         ad1=0
         do a=1,nvir
            do i=1,nocc
               ad1=ad1+1
               t1(a,i)=buf1(ad1)
            end do
         end do
      end if
      call ga_brdcst(msg_cc_t1b,t1,mdtob(lnov),0)
!
!     Next phase transposes bits withing g_objo.  To avoid OVV local
!     memory, duplicate g_objo and do some smaller operations
!     to accomplish the sort in O(V) memory.
!     
      call ga_sync
!
      if (.not. ga_duplicate(g_objo, g_tmp, 'rtdpro')) call errquit
     $     ('ccsd_rdtpro: failed to allocate tmp GA',0, GA_ERR)
      call ga_copy(g_objo, g_tmp)
!
      do i=1+iam,nocc,nodes
!     t2j
!*         start = 1
!*         end   = lnovv
!*         call ga_get(g_tmp,start,end,i,i,buf1,lnovv)
!*         ad1 = 0
         do j=1,nocc
            do av=1,nvir
!*               do bv=1,nvir
!*                  ad1=ad1+1
!*                  buf2(ad1)=buf1((av-1)*lnov+(j-1)*nvir+bv)
!*               end do
               start = (av-1)*lnov+(j-1)*nvir+1
               end   = start + nvir - 1
               call ga_get(g_tmp,start,end,i,i,buf1,nvir)
!
               start = (av-1)*nvir+(j-1)*nvir*nvir+1
               end   = start + nvir - 1
               call ga_put(g_objo,start,end,i,i,buf1,nvir)
            end do
         end do
!*         call ga_put(g_objo,start,end,i,i,buf2,lnovv)
!
!     eintc
         start = lnovv + 1
         end   = start + 2*lnoov - 1
         call ga_get(g_tmp,start,end,i,i,buf1,2*lnoov)
         ad1=0
         do j=1,nocc
            do av=1,nvir
               do k=1,nocc
                  ad1=ad1+1
                  ad2=(av-1)*lnoo+(j-1)*nocc+k
                  buf2(ad1)=buf1(ad2)
!     buf1(ad1)=gtint(i,a,j,k)
               end do
            end do
         end do
!     - eintx
         do j=1,nocc
            do av=1,nvir
               do k=1,nocc
                  ad1=ad1+1
                  ad2=lnoov+(av-1)*lnoo+(j-1)*nocc+k
                  buf2(ad1)=buf1(ad2)
!     buf1(ad1)=gtint(i,k,j,a)
               end do
            end do
         end do
!
         call ga_put(g_objo,start,end,i,i,buf2,2*lnoov)
!     
      end do
!     
      call ga_sync()
!
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('ccsd_trdpro: failed to destroy GA',0, GA_ERR)
!
      if (occsdps) then
         call pstat_off(ps_rdtrpo)
      else
         call qexit('rdtrpo',0)
      endif
!     
      end
