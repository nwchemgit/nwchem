      subroutine ccsd_rdtrpo(t1,buf1,buf2,g_objo,g_objv,
     $                       nocc,nvir,iprt)
      implicit none
c
      integer lnoo,lnov,lnvv,lnoov,lnooo,lnobj,lnovv
      common/len/lnoo,lnov,lnvv,lnoov,lnooo,lnobj,lnovv
c
      integer g_objo,g_objv,nocc,nvir,iprt
      double precision t1(nvir,nocc),buf1(lnobj),buf2(lnobj)
c
      integer il,ih,jl,jh,i,j,k,a,av,bv,ad1,ad2,
     &        lbfo,lbfv
C
#include "msgids.fh"
#include "tcgmsg.fh"
#include "para1.fh"
c
c - read in objects
c - [io|ov], [oo|iv], t2(io:vv) / [ao|vo], t2(oo:av), t2(oo,va)
c
      lbfo=lnovv+lnoov+lnoov
      lbfv=lnoov+lnoov+lnoov
      call ga_distribution(g_objv,iam,jl,jh,il,ih)
      do av=1,nvir
       if (av.ge.il.and.av.le.ih)then
c
c get dint and eint
        call ga_get(g_objv,1,2*lnoov+lnooo,av,av,buf1,
     &              2*lnoov+lnooo)
c
c process eints for occupied indexed buffer
        do i=1,nocc
c eintc
         ad1=0
         do j=1,nocc
          do k=1,nocc
           ad1=ad1+1
           buf2(ad1)=buf1(lnoov+(j-1)*lnoo+(i-1)*nocc+k)
          enddo
         enddo
         ad2=lnovv+(av-1)*lnoo
         call ga_put(g_objo,ad2+1,ad2+lnoo,i,i,buf2,lnoo)
c eintx
         ad1=lnoov+(i-1)*lnoo+1
         ad2=lnovv+lnoov+(av-1)*lnoo
         call ga_put(g_objo,ad2+1,ad2+lnoo,i,i,buf1(ad1),lnoo)
        enddo
c
c process dints for virtual indexed buffer
        ad1=0
        do i=1,nocc
         do j=1,nocc
          do bv=1,nvir
           ad1=ad1+1
           ad2=(bv-1)*lnoo+(i-1)*nocc+j
           buf2(ad1)=buf1(ad2)
          enddo
         enddo
        enddo
c
c process t2 amplitudes for virtual indexed buffer
c t2c
        do i=1,nocc
         do j=1,nocc
          do bv=1,nvir
           ad1=ad1+1
           ad2=lnoov+lnooo+(bv-1)*lnoo+(i-1)*nocc+j
           buf2(ad1)=buf1(ad2)
          enddo
         enddo
        enddo
c t2x
        do i=1,nocc
         do j=1,nocc
          do bv=1,nvir
           ad1=ad1+1
           ad2=lnoov+lnooo+(bv-1)*lnoo+(j-1)*nocc+i
           buf2(ad1)=buf1(ad2)
          enddo
         enddo
        enddo
        if (iprt.gt.50)then
         print *,'t2c new ',iam,av
         write(6,4859)(buf2(j),j=lnoov+1,2*lnoov)
         print *,'t2x new ',iam,av
         write(6,4859)(buf2(j),j=2*lnoov+1,3*lnoov)
         print *,'dint new ',iam,av
         write(6,4859)(buf2(j),j=1,lnoov)
 4859    format(1x,5e14.4)
        endif
        call ga_put(g_objv,1,ad1,av,av,buf2,ad1)
c
c process t2 amplitudes for occupied indexed buffer
c t2j
        do i=1,nocc
         ad1=0
         do j=1,nocc
          do bv=1,nvir
           ad1=ad1+1
           buf2(ad1)=buf1(lnoov+lnooo+(bv-1)*lnoo+(i-1)*nocc+j)
          enddo
         enddo
         ad2=(av-1)*lnov
         call ga_put(g_objo,ad2+1,ad2+lnov,i,i,buf2,lnov)
        enddo
c
       endif
      enddo
c
c sort t1 then broadcast which also serves to synchronize processes
      if (iam.eq.0)then
       call ga_get(g_objv,lnoov+lnooo+1,lnoov+lnooo+lnov,nvir+1,nvir+1,
     &             buf1,lnov)
       ad1=0
       do a=1,nvir
        do i=1,nocc
         ad1=ad1+1
         t1(a,i)=buf1(ad1)
        enddo
       enddo
      endif
      call ga_brdcst(msg_cc_t1b,t1,mdtob(lnov),0)
c
      call ga_distribution(g_objo,iam,jl,jh,il,ih)
      do i=1,nocc
       if (i.ge.il.and.i.le.ih)then
        call ga_get(g_objo,1,lbfo,i,i,buf1,lbfo)
c t2j
        ad1=0
        do j=1,nocc
         do av=1,nvir
          do bv=1,nvir
           ad1=ad1+1
           buf2(ad1)=buf1((av-1)*lnov+(j-1)*nvir+bv)
          enddo
         enddo
        enddo
c eintc
        do j=1,nocc
         do av=1,nvir
          do k=1,nocc
           ad1=ad1+1
           ad2=lnovv+(av-1)*lnoo+(j-1)*nocc+k
           buf2(ad1)=buf1(ad2)
c          buf1(ad1)=gtint(i,a,j,k)
          enddo
         enddo
        enddo
c - eintx
        do j=1,nocc
         do av=1,nvir
          do k=1,nocc
           ad1=ad1+1
           ad2=lnovv+lnoov+(av-1)*lnoo+(j-1)*nocc+k
           buf2(ad1)=buf1(ad2)
c          buf1(ad1)=gtint(i,k,j,a)
          enddo
         enddo
        enddo

        if (iprt.gt.50)then
         print *,'t2j new ',iam,i
         write(6,4859)(buf2(j),j=1,lnovv)
         print *,'eintc new ',iam,i
         write(6,4859)(buf2(j),j=lnovv+1,lnovv+lnoov)
         print *,'eintx new ',iam,i
         write(6,4859)(buf2(j),j=lnovv+lnoov+1,lnovv+2*lnoov)
        endif
        call ga_put(g_objo,1,ad1,i,i,buf2,ad1)
c
       endif
      enddo

      call ga_sync()
c      print *,'end rdtrpo2 '
      return
      end
