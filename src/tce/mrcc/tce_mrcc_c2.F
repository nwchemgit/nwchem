c      SUBROUTINE tce_mrcc_c2(size_t2m,k_t1_offsetm,d_t1m,nref,iref,
c     1 K_T2_OFFSETM,d_t2m,d_r2m,k_r2_offsetm,size_t1m)
      SUBROUTINE tce_mrcc_c2(size_t2m,k_t1_offsetm,d_t1m,iref,
     1 K_T2_OFFSETM,d_t2m,d_r2m,k_r2_offsetm,size_t1m)
c C2 + t2(mu/nu) - t2(mu)
c          call tce_mrcc_c2(size_t2m,k_t1_offsetm,d_t1m,nref,iref,
c     1 k_t2_offsetm,d_t2m,d_r2m)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "util.fh"
#include "msgids.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_mrcc.fh"
#include "stdio.fh"

      integer inoab,invab
      integer size_t2m(maxref),size_t1m(maxref)
      integer k_t1_offsetm(maxref)
      integer d_t1m(maxref)
      integer iref
      integer size
      integer dima,dimb
      integer p1b,p2b,h3b,h4b
      integer l_t2,k_t2
      integer l_t1a1,k_t1a1
      integer l_t1b1,k_t1b1
      integer l_t1a2,k_t1a2
      integer l_t1b2,k_t1b2
      integer l_tmp,k_tmp
      integer k_t2_offsetm(maxref),d_t2m(maxref)
      integer d_r2m(maxref)
      integer k_t2_2,l_t2_2
      integer l_r2,k_r2
      integer l_c,k_c
      integer irefnew
      integer p1,p2,h3,h4
      integer i
      character*255 filename
      integer d_c2tmp,d_t1c
      integer k_r2_offsetm(maxref)
      double precision coef
      integer d_rc2
      double precision sign
      integer ijblock(nref,noab+nvab)
      integer aux1,aux2
      logical nodezero
      double precision itcpu,itwall

      nodezero=(ga_nodeid().eq.0)

      if(lusesub) then

          call tce_filename('c2_tmp',filename)
          call pgcreatefile(filename,d_c2tmp,size_t2m(iref),mypgid)

          call tce_filename('rc2',filename)
          call pgcreatefile(filename,d_rc2,size_t2m(iref),mypgid)

          call tce_filename('t1c',filename)
          call pgcreatefile(filename,d_t1c,size_t1m(iref),mypgid)

      else

          call tce_filename('c2_tmp',filename)
          call createfile(filename,d_c2tmp,size_t2m(iref))

          call tce_filename('rc2',filename)
          call createfile(filename,d_rc2,size_t2m(iref))

          call tce_filename('t1c',filename)
          call createfile(filename,d_t1c,size_t1m(iref))

      endif

      do irefnew=1,nref
       do i=1,noab+nvab

           aux1 = (1 - (int_mb(k_spinm(iref)+i-1)-1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+i-1)))/2

        aux1 = moindexes(aux1,int_mb(k_spinm(iref)+i-1),iref)
        aux1 = moindexes(aux1,int_mb(k_spinm(iref)+i-1),irefnew)

         ijblock(irefnew,i) = orbinblck(aux1,int_mb(k_spinm(iref)+i-1),
     1 irefnew)

       enddo
      enddo

      do irefnew=1,nref
      if(irefnew.ne.iref) then

      coef = dbl_mb(k_heff+iref-1+(irefnew-1)*nref)*
     1 dbl_mb(k_sqc+(irefnew-1)*nref+mkroot-1)
      if(dbl_mb(k_sqc+(iref-1)*nref+mkroot-1).lt.0.0d0)
     1 coef=-1.0d0*coef

c     2dbl_mb(k_sqc+(mkroot-1)*nref+iref-1)
c      write(6,*)'USING',coef,
c     1dbl_mb(k_heff+iref-1+(irefnew-1)*nref),
c     2abs(dbl_mb(k_sqc+(irefnew-1)*nref+mkroot-1)),
c     3abs(dbl_mb(k_sqc+(iref-1)*nref+mkroot-1))


      call ga_zero(d_c2tmp)
      call ga_zero(d_t1c)

        itcpu = - util_cpusec()
        itwall = - util_wallsec()

      call tce_mrcc_c1_offset(k_t1_offsetm,d_t1m,iref,irefnew,d_t1c)
      if(lusesub) then
        call ga_pgroup_sync(mypgid)
      else
        call reconcilefile(d_t1c,size_t1m(iref))
      endif

          itcpu = itcpu + util_cpusec()
          itwall = itwall + util_wallsec()

c          if(nodezero) then
c            write(LuOut,9222) itcpu, itwall
c            call util_flush(LuOut)
c          endif

        itcpu = - util_cpusec()
        itwall = - util_wallsec()

      call tce_mrcc_c2_offset(k_t2_offsetm,d_t2m,iref,irefnew,d_c2tmp,
     1 ijblock,nref)
      if(lusesub) then
        call ga_pgroup_sync(mypgid)
      else
        call reconcilefile(d_c2tmp,size_t2m(iref))
      endif

          itcpu = itcpu + util_cpusec()
          itwall = itwall + util_wallsec()

          if(nodezero) then
c            write(LuOut,9223) itcpu, itwall
            call util_flush(LuOut)
          endif

c      write(6,*)'T2T2T2T2T2T2T2',iref,irefnew
c      call tce_mrcc_print_t2type(d_t2m(iref),k_t2_offsetm(iref),
c     1iref,irefnew)
      
c       write(6,*)'SELECTED',iref,irefnew
c      call tce_mrcc_print_t2type(d_c2tmp,k_t2_offsetm(iref),
c     1iref,irefnew)

c      write(6,*)'0'
c      call tce_mrcc_c1_t1t1(d_t1m(iref),k_t1_offsetm(iref),d_t1c,
c     1 k_t1_offsetm(iref),d_c2tmp,k_r2_offsetm(iref))
c      call reconcilefile(d_c2tmp,size_t2m(iref))

c      write(6,*)'1'
c      call tce_mrcc_c1_t1t1_2(d_t1c,k_t1_offsetm(iref),d_t1c,
c     1 k_t1_offsetm(iref),d_c2tmp,k_r2_offsetm(iref))
c      call reconcilefile(d_c2tmp,size_t2m(iref))
      if(lusesub) then
      call T2_C2_1x(d_rc2,d_t1c,k_t2_offsetm(iref),
     1 k_t1_offsetm(iref),coef)
      else
      call T2_C2_1a(d_rc2,d_t1c,k_t2_offsetm(iref),
     1 k_t1_offsetm(iref),coef)
      endif

      if(lusesub) then
        call ga_pgroup_sync(mypgid)
      else
        call reconcilefile(d_rc2,size_t2m(iref))
      endif

      if(lusesub) then
      call T2_pg_C2(d_rc2,d_t1m(iref),d_t2m(iref),d_t1c,d_c2tmp,
     1 k_t2_offsetm(iref),k_t1_offsetm(iref),k_t2_offsetm(iref),
     2 k_t1_offsetm(iref),k_t2_offsetm(iref),coef)
      else
      call T2_C2(d_rc2,d_t1m(iref),d_t2m(iref),d_t1c,d_c2tmp,
     1 k_t2_offsetm(iref),k_t1_offsetm(iref),k_t2_offsetm(iref),
     2 k_t1_offsetm(iref),k_t2_offsetm(iref),coef)
      endif

      if(lusesub) then
        call ga_pgroup_sync(mypgid)
      else
        call reconcilefile(d_rc2,size_t2m(iref))
      endif

c      write(6,*)'2'
c      call tce_mrcc_c1_t1t1_2(d_t1m(iref),k_t1_offsetm(iref),
c     1 d_t1m(iref),k_t1_offsetm(iref),d_c2tmp,k_r2_offsetm(iref))
c      call reconcilefile(d_c2tmp,size_t2m(iref))

c      write(6,*)'3'
c      call tce_mrcc_c2_offset(k_t2_offsetm,d_t2m,iref,irefnew,d_c2tmp)

c      call tce_mrcc_print_t2type(d_t2m(iref),k_t2_offsetm(iref),
c     1iref,irefnew)
c      call tce_mrcc_print_t2type(d_t2m(iref),k_t2_offsetm(iref),
c     1iref,irefnew)

c      call ga_add(1.0d0,d_c2tmp,-1.0d0,d_t2m(iref),d_c2tmp)
c      call tce_mrcc_mul_2(1.0d0,d_c2tmp,-1.0d0,d_t2m(iref),
c     1 k_t2_offsetm(iref),d_c2tmp,iref)
c      call reconcilefile(d_c2tmp,size_t2m(iref))

c      write(6,*)'COUPLING',iref
c      call tce_mrcc_print_t2type(d_c2tmp,k_t2_offsetm(iref),
c     1iref,irefnew)

c      call ga_add(1.0d0,d_rc2,coef,d_c2tmp,d_rc2)
c      call tce_mrcc_mul_2(1.0d0,d_rc2,coef,d_c2tmp,
c     1 k_t2_offsetm(iref),d_rc2,iref)
c      call reconcilefile(d_rc2,size_t2m(iref))

      endif 
      enddo

c       write(6,*)'M2M2M2M2M2M2M2',iref,irefnew
c      call tce_mrcc_print_t2type(d_rc2,k_t2_offsetm(iref),
c     1iref,irefnew)

c      call ga_print(d_t2m(iref))
      call util_flush(LuOut)

c      write(6,*)'INITIAL R2',iref
c      call tce_mrcc_print_t2type(d_r2m(iref),k_t2_offsetm(iref),
c     1iref,irefnew)

c      write(6,*)'USING R2mu',dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)

c abs(dbl_mbik_sqc+(mkroot-1)*nref+iref-1))
c      call ga_add(abs(dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)),d_r2m(iref),
c     1 1.0d0,d_rc2,d_r2m(iref))
      call tce_mrcc_mul_2(abs(dbl_mb(k_sqc+(iref-1)*nref+mkroot-1)),
     1 d_r2m(iref),1.0d0,d_rc2,k_t2_offsetm(iref),d_r2m(iref),iref)
      if(lusesub) then
        call ga_pgroup_sync(mypgid)
      else
        call reconcilefile(d_r2m(iref),size_t2m(iref))
      endif

c       write(6,*)'FINAL R2',iref,irefnew
c      call tce_mrcc_print_t2type(d_r2m(iref),k_t2_offsetm(iref),
c     1iref,irefnew)

      call deletefile(d_rc2)
      call deletefile(d_c2tmp)
      call deletefile(d_t1c) 

 9222 format(1x,'Corresponding T1 / sec',2f15.1)
 9223 format(1x,'Corresponding T2 / sec',2f15.1)

      return
      end
c $Id$
