      subroutine ccsd_rdtrpo(t1,buf1,buf2,g_objo,g_objv,
     $     nocc,nvir,iprt)
C     $Id: ccsd_rdtrpo.F,v 2.8 2003-10-17 22:54:25 carlfahl Exp $
      implicit none
#include "errquit.fh"
c     
#include "ccsd_len.fh"
      integer g_objo,g_objv,nocc,nvir,iprt
      double precision t1(nvir,nocc),buf1(lnobj),buf2(lnobj)
c     
      integer i,j,k,a,av,bv,ad1,ad2,lbfo
      Integer Nodes, IAm
      integer g_tmp, start, end
C     
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
C     
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
      call ga_sync
c     
c     - read in objects
c     - [io|ov], [oo|iv], t2(io:vv) / [ao|vo], t2(oo:av), t2(oo,va)
c     
      lbfo=lnovv+lnoov+lnoov
c     
c     do av ... uses 2*oov+ooo memory.  Don't use locality
c     so that we can distribute more uniformly across nodes.
c     Locality is not used in the triples driver.
c     
      do av=iam+1,nvir,nodes 
c     
c     get dint and eint
         call ga_get(g_objv,1,2*lnoov+lnooo,av,av,buf1,
     &        2*lnoov+lnooo)
c     
c     process eints for occupied indexed buffer
         do i=1,nocc
c     eintc
            ad1=0
            do j=1,nocc
               do k=1,nocc
                  ad1=ad1+1
                  buf2(ad1)=buf1(lnoov+(j-1)*lnoo+(i-1)*nocc+k)
               end do
            end do
            ad2=lnovv+(av-1)*lnoo
            call ga_put(g_objo,ad2+1,ad2+lnoo,i,i,buf2,lnoo)
c     eintx
            ad1=lnoov+(i-1)*lnoo+1
            ad2=lnovv+lnoov+(av-1)*lnoo
            call ga_put(g_objo,ad2+1,ad2+lnoo,i,i,buf1(ad1),lnoo)
         end do
c     
c     process dints for virtual indexed buffer
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
c     
c     process t2 amplitudes for virtual indexed buffer
c     t2c
         do i=1,nocc
            do j=1,nocc
               do bv=1,nvir
                  ad1=ad1+1
                  ad2=lnoov+lnooo+(bv-1)*lnoo+(i-1)*nocc+j
                  buf2(ad1)=buf1(ad2)
               end do
            end do
         end do
c     t2x
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
c     
c     process t2 amplitudes for occupied indexed buffer
c     t2j
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
c     
      end do
c     
      call ga_sync()            ! ga_sync ensures all GA ops complete
c     
c     sort t1 then broadcast
c     
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
c
c     Next phase transposes bits withing g_objo.  To avoid OVV local
c     memory, duplicate g_objo and do some smaller operations
c     to accomplish the sort in O(V) memory.
c     
      call ga_sync
c
      if (.not. ga_duplicate(g_objo, g_tmp, 'rtdpro')) call errquit
     $     ('ccsd_rdtpro: failed to allocate tmp GA',0, GA_ERR)
      call ga_copy(g_objo, g_tmp)
c
      do i=1+iam,nocc,nodes
c     t2j
**         start = 1
**         end   = lnovv
**         call ga_get(g_tmp,start,end,i,i,buf1,lnovv)
**         ad1 = 0
         do j=1,nocc
            do av=1,nvir
**               do bv=1,nvir
**                  ad1=ad1+1
**                  buf2(ad1)=buf1((av-1)*lnov+(j-1)*nvir+bv)
**               end do
               start = (av-1)*lnov+(j-1)*nvir+1
               end   = start + nvir - 1
               call ga_get(g_tmp,start,end,i,i,buf1,nvir)
c
               start = (av-1)*nvir+(j-1)*nvir*nvir+1
               end   = start + nvir - 1
               call ga_put(g_objo,start,end,i,i,buf1,nvir)
            end do
         end do
**         call ga_put(g_objo,start,end,i,i,buf2,lnovv)
c
c     eintc
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
c     buf1(ad1)=gtint(i,a,j,k)
               end do
            end do
         end do
c     - eintx
         do j=1,nocc
            do av=1,nvir
               do k=1,nocc
                  ad1=ad1+1
                  ad2=lnoov+(av-1)*lnoo+(j-1)*nocc+k
                  buf2(ad1)=buf1(ad2)
c     buf1(ad1)=gtint(i,k,j,a)
               end do
            end do
         end do
c
         call ga_put(g_objo,start,end,i,i,buf2,2*lnoov)
c     
      end do
c     
      call ga_sync()
c
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('ccsd_trdpro: failed to destroy GA',0, GA_ERR)
c     
      end
