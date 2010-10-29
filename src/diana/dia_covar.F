      subroutine dia_covar()
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "util.fh"
c
      integer i,j,k
c
      filpca=card(8:80)
      if(index(filpca,'.pca').eq.0)
     + filpca=filpca(1:index(filpca,' ')-1)//'.pca'
c
      if(nsel.eq.0) call md_abort('No atoms selected',0)
c
      if(active) call dia_edfinal()
      call dia_edinit()
      active=.true.
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
c      call dia_xaver(dbl_mb(ndx),nsa,min((me+1)*nfrme,nfrtot)-me*nfrdim,
c     + dbl_mb(i_xadat),nsa)
c
      call dia_cov(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
     + dbl_mb(ndx),dbl_mb(i_xadat))
c
      call ga_scale(ga_cov,1.0d0/dble(nfrtot))
c
c      call dia_pcovar(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel))
c
      if(me.eq.0) then
      if(util_print('frames_read',print_debug)) then
      write(*,'(/,a)') ' Covariance matrix generation completed'
      endif
      endif
c
      call ga_sync()
c
      call ga_diag_std(ga_cov,ga_vec,dbl_mb(i_cov))
c
      call ga_sync()
      do 19 i=1,3*nsel
      int_mb(i_ord+i-1)=i
   19 continue
      do 17 i=1,3*nsel-1
      do 18 j=i+1,3*nsel
      if(abs(dbl_mb(i_cov-1+int_mb(i_ord+i-1))).lt.
     + abs(dbl_mb(i_cov-1+int_mb(i_ord+j-1)))) then
      k=int_mb(i_ord+i-1)
      int_mb(i_ord+i-1)=int_mb(i_ord+j-1)
      int_mb(i_ord+j-1)=k 
      endif
   18 continue
   17 continue
c
      if(me.eq.0) then
      open(unit=lfnpca,file=filpca(1:index(filpca,' ')-1),
     + form='formatted',status='unknown')
      write(lfnpca,1000) 3*nsel,nsa
 1000 format(2i7)
      write(lfnpca,1001) (dbl_mb(i_cov-1+int_mb(i_ord+i-1)),i=1,3*nsel)
 1001 format(6e12.5)
      endif
c
      return
      end
      subroutine dia_edinit()
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      if(me.eq.0) then
      write(*,1000) (ma_sizeof(mt_dbl,3*nsel,mt_byte)+
     + ma_sizeof(mt_int,6*nsel,mt_byte))/104857,
     + ma_inquire_avail(mt_byte)/1048576,
     + ma_sizeof(mt_dbl,18*nsel*nsel/np,mt_byte)/1048576,
     + ga_memory_avail()/1048576
 1000 format(/,' Memory   Required    Available',//,
     + '   ma   ',2i10,' MB',/,
     + '   ga   ',2i10,' MB',/)
      endif
c      print*,'nsel=',nsel
c
c     create nsel x nsel global array covariance matrix
c
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'covar',0,3*nsel,ga_cov))
     + call md_abort('Failed to create global array covar',me)
      call ga_zero(ga_cov)
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'vector',0,3*nsel,ga_vec))
     + call md_abort('Failed to create global array vector',me)
c
c     allocate memory
c
      if(.not.ma_push_get(mt_int,3*nsel,'ndx',l_ndx,i_ndx))
     + call md_abort('Failed to allocate memory for ndx',me)
      if(.not.ma_push_get(mt_dbl,3*nsel,'cov',l_cov,i_cov))
     + call md_abort('Failed to allocate memory for cov',me)
      if(.not.ma_push_get(mt_int,3*nsel,'ord',l_ord,i_ord))
     + call md_abort('Failed to allocate memory for ord',me)
c
      if(me.eq.0) then
      write(*,1002) 3*nsel
 1002 format(' Memory allocated for covariance ',i5)
      endif
c
      return
      end
      subroutine dia_edfinal()
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_ord))
     + call md_abort('Failed to deallocate memory for ord',me)
      if(.not.ma_pop_stack(l_cov))
     + call md_abort('Failed to deallocate memory for cov',me)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Failed to deallocate memory for ndx',me)
c
c     destroy global array covariance matrix
c
      if(.not.ga_destroy(ga_vec))
     + call md_abort('Failed to destroy global array vector',me)
      if(.not.ga_destroy(ga_cov))
     + call md_abort('Failed to destroy global array covar',me)
c
      return
      end
      subroutine dia_cov(ndex,cov,isel,xtrj,xa)
c
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer ndex(*)
      real*8 cov(*)
      real*8 xtrj(ntrj,3,nfrdim),xa(nsa,3)
      integer isel(msa)
c
      integer i,j,k,ia,ja,ix,jx
      integer ilc,ihc,jlc,jhc
      integer iproc
c
c      write(*,'(a)') 'DIA_COVAR'
c
      j=0
      do 1 i=1,nsa
      if(isel(i).ne.0) then
      j=j+1
      ndex(j)=i
      endif
    1 continue
c
      do 6 iproc=0,np-1
c
      call ga_distribution(ga_cov,iproc,ilc,ihc,jlc,jhc)
c
      do 3 i=ilc,ihc
      ia=ndex((i+2)/3)
      ix=mod(i+2,3)+1
      do 5 j=jlc,jhc
      cov(j-jlc+1)=0.0d0
    5 continue
      do 2 k=1,nfrme
      do 4 j=jlc,jhc
      ja=ndex((j+2)/3)
      jx=mod(j+2,3)+1
      cov(j-jlc+1)=cov(j-jlc+1)+
     + (xtrj(ia,ix,k)-xa(ia,ix))*(xtrj(ja,jx,k)-xa(ja,jx))
    4 continue
    2 continue
      call ga_acc(ga_cov,i,i,jlc,jhc,cov,1,one)
    3 continue
    6 continue
c
      call ga_sync()
c
      return
      end
      subroutine dia_project()
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "util.fh"
c
      character*3 fmt
      character*255 fil
      integer i,ivec,nbatch
      integer ga_prj,ilp,ihp,jlp,jhp,ndxp,ldp
      integer ga_val,ilv,ihv,jlv,jhv,ndxv,ldv
      integer kvec,ip
      integer i_a,l_a,i_pcav,l_pcav
      real*8 d
c
      lesd=.true.
      read(card(8:14),'(i7)') ivec
      read(card(15:21),'(i7)') nbatch
      fil=card(22:80)
      i=index(fil,'.')
      if(i.eq.0) then
      fmt='trj'
      filprj=fil(1:index(fil,' ')-1)//'.trj'
      filval=fil(1:index(fil,' ')-1)//'.val'
      filvec=fil(1:index(fil,' ')-1)//'.vec'
      filtcl=fil(1:index(fil,' ')-1)//'.tcl'
      filecc=fil(1:index(fil,' ')-1)//'.ecce.out'
      filmin=fil(1:index(fil,' ')-1)//'_min.pdb'
      filmax=fil(1:index(fil,' ')-1)//'_max.pdb'
      filamp=fil(1:index(fil,' ')-1)//'.amp'
      i=index(filprj,'.')
      else
      filprj=fil
      fmt=filprj(i+1:i+3)
      filval=fil(1:i-1)//'.val'
      filvec=fil(1:i-1)//'.vec'
      filtcl=fil(1:i-1)//'.tcl'
      filecc=fil(1:i-1)//'.ecce.out'
      filmin=fil(1:i-1)//'_min.pdb'
      filmax=fil(1:i-1)//'_max.pdb'
      filamp=fil(1:i-1)//'.amp'
      endif
c
      if(nsel.eq.0) call md_abort('No atoms selected',0)
      if(.not.active) call md_abort('No ED matrix active',0)
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
c     allocate global array to hold the projected trajectory
c
      if(.not.ga_create(mt_dbl,nsa*3,nfrtot,'prj',nsa*3,nfrdim,ga_prj))
     + call md_abort('Failed to create ga_prj',0)
      call ga_distribution(ga_prj,me,ilp,ihp,jlp,jhp)
      call ga_access(ga_prj,ilp,ihp,jlp,jhp,ndxp,ldp)
c
c     allocate global array to hold the projection value
c
      if(.not.ga_create(mt_dbl,1,nfrtot,'val',1,nfrdim,ga_val))
     + call md_abort('Failed to create ga_val',0)
      call ga_distribution(ga_val,me,ilv,ihv,jlv,jhv)
      call ga_access(ga_val,ilv,ihv,jlv,jhv,ndxv,ldv)
c
c     retrieve the specified eigen vector
c
      kvec=int_mb(i_ord+ivec-1)
      call ga_get(ga_vec,1,3*nsel,kvec,kvec,dbl_mb(i_cov),1)
c
c     allocate temporary array a(nfr)
c
      if(.not.ma_push_get(mt_dbl,nfrdim,'a',l_a,i_a))
     + call md_abort('Error allocating a',0)
c
      call dia_proj(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
     + dbl_mb(ndx),dbl_mb(ndxp),
     + dbl_mb(i_xadat),dbl_mb(ndxv))
c
      if(.not.ma_pop_stack(l_a)) call md_abort('Error deallocating a',0)
c
c      if(me.eq.0) then
c      open(unit=lfnprj,file=filprj(1:index(filprj,' ')-1),
c     + form='formatted',status='unknown',err=9999)
c      open(unit=lfnval,file=filval(1:index(filval,' ')-1),
c     + form='formatted',status='unknown',err=9999)
c      endif
c
c      call dia_wthdr(lfnprj,fmt,byte_mb(i_snam),byte_mb(i_tag),
c     + int_mb(i_isel),.false.)
c
      call ga_sync()
c
      call dia_wrttrj(ga_prj,filprj(1:i-1),fmt,nbatch,.false.)
c
      if(me.eq.0) then
      if(.not.ma_push_get(mt_dbl,nfrtot,'ipcav',l_pcav,i_pcav))
     + call md_abort('Failed to allocate ipcav',0)
      write(lfnpca,1000) ivec,nfrtot
 1000 format(i5,i10)
c
c     write the eigenvector
c
      ip=0
      do 1 i=1,nsel
      d=sqrt(dbl_mb(i_cov+ip)**2+dbl_mb(i_cov+ip+1)**2
     + +dbl_mb(i_cov+ip+2)**2)
      write(lfnpca,1001) dbl_mb(i_cov+ip),dbl_mb(i_cov+ip+1),
     + dbl_mb(i_cov+ip+2),d
 1001 format(4f12.6)
      ip=ip+3
    1 continue
      do 2 ip=0,np-1
      call ga_distribution(ga_val,ip,ilv,ihv,jlv,jhv)
      call ga_get(ga_val,1,1,jlv,jhv,dbl_mb(i_pcav+ip*nfrdim),1)
    2 continue
      write(lfnpca,1002) (dbl_mb(i_pcav+i-1),i=1,nfrtot)
 1002 format(5e12.5)
      if(.not.ma_pop_stack(l_pcav))
     + call md_abort('Failed to deallocate ipcav',0)
      endif
c
      if(.not.ga_destroy(ga_val))
     + call md_abort('Could not destroy ga_val',0)
      if(.not.ga_destroy(ga_prj))
     + call md_abort('Could not destroy ga_prj',0)
c
      return
 9999 continue
      call md_abort(' Could not open projection files',0)
      return
      end
      subroutine dia_proj(ndex,vct,isel,xtrj,xprj,xa,val)
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ndex(nsel)
      real*8 vct(*),val(nfrdim),pnorm,vval
      real*8 xtrj(ntrj,3,nfrdim),xprj(ntrj,3,nfrdim),xa(nsa,3)
      integer isel(nsa)
c
      integer i,j,k,l,m,ivmin,ivmax,ndxval,i_valt,l_valt
      real*8 value,vmin,vmax,vming,vmaxg
c
c      write(*,'(3f12.6)') ((xtrj(5,i,k),i=1,3),k=1,5)
c      write(*,'(5i5)') (ndex(i),i=1,5)
c      write(*,'(6f12.6)') (vct(k),k=1,3*nsel)
c      write(*,'(3f12.6)') (xa(5,i),i=1,3)
      do 1 k=1,nfrme
      val(k)=0.0d0
      pnorm=0.0d0
      m=0
      do 2 i=1,nsel
      vval=0.0d0
      do 3 j=1,3
      m=m+1
      val(k)=val(k)+(xtrj(ndex(i),j,k)-xa(ndex(i),j))*vct(m)
      pnorm=pnorm+vct(m)*vct(m)
    3 continue
    2 continue
      val(k)=val(k)/sqrt(pnorm)
    1 continue
c      write(*,'(5f12.6)') (val(k),k=1,nfr)
c
      do 4 k=1,nfrme
      m=0
      do 5 i=1,nsel
      do 6 j=1,3
      m=m+1
      xprj(ndex(i),j,k)=val(k)*vct(m)+xa(ndex(i),j)
c      write(*,'(4i5,3f12.6)') k,i,j,ndex(i),val(k),vct(m),xa(ndex(i),j)
    6 continue
    5 continue
    4 continue
c
c     write the vector
c
      if(me.eq.0) then
      open(unit=lfntcl,file=filtcl(1:index(filtcl,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnamp,file=filamp(1:index(filamp,' ')-1),
     + form='formatted',status='unknown',err=9999)
      k=0
      l=0
      do 15 i=1,nsel
      value=0.0d0
      do 16 j=1,3
      value=value+(dbl_mb(i_cov+k))**2
      k=k+1
   16 continue
      l=l+1
      write(lfntcl,5555)
     + (10.0*xa(ndex(i),j),j=1,3),
     + (dbl_mb(i_cov+k-4+j),j=1,3)
 5555 format('draw vector2 {',3f12.6,'} {',3f12.6,'} 1.0 10')
      write(lfnamp,'(i5,f12.6)') i,sqrt(value)
   15 continue
      close(unit=lfntcl)
      close(unit=lfnamp)
 9999 continue
      endif
c
c     write the minimum and maximum projection structures as pdb formatted files
c
      ivmin=1
      ivmax=1
      do 17 k=2,nfrme
      if(val(k).gt.val(ivmax)) ivmax=k
      if(val(k).lt.val(ivmin)) ivmin=k
   17 continue
      vmax=val(ivmax)
      vmin=val(ivmin)
      call ga_dgop(mag_d01,vmax,1,'max')
      call ga_dgop(mag_d02,vmin,1,'min')
      if(val(ivmax).eq.vmax) then
      call dia_wrtpdb(.true.,lfnmax,filmax,byte_mb(i_snam),isel,
     + xprj(1,1,ivmax),isel,xprj,.true.,.false.,isel)
      endif
      if(val(ivmin).eq.vmin) then
      call dia_wrtpdb(.true.,lfnmin,filmin,byte_mb(i_snam),isel,
     + xprj(1,1,ivmin),isel,xprj,.true.,.false.,isel)
      endif
c
      ndxval=np*nfrme
      if(.not.ma_push_get(mt_dbl,ndxval,'val',l_valt,i_valt))
     + call md_abort('Failed to allocate val',0)
      do 18 i=1,ndxval
      dbl_mb(i_valt-1+i)=0.0d0
   18 continue
      do 19 i=1,nfrme
      dbl_mb(i_valt+me*nfrme+i-1)=val(i)
   19 continue
      call ga_dgop(mag_d03,dbl_mb(i_valt),ndxval,'+')
      if(me.eq.0) then
      open(unit=lfnval,file=filval(1:index(filval,' ')-1),
     + form='formatted',status='unknown',err=9999)
      write(lfnval,5566) (i,dbl_mb(i_valt-1+i),i=1,ndxval)
 5566 format(i8,f12.6)
      close(unit=lfnval)
      endif
      if(.not.ma_pop_stack(l_valt))
     + call md_abort('Failed to deallocate val',0)
c
      return
      end
      subroutine dia_pcovar(ndex,cov,isel)
c
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer ndex(*)
      real*8 cov(*)
      integer isel(msa)
c
      integer i,j,k,ia,ja,ix,jx
      integer iproc
c
      if(me.ne.0) return
c
      write(*,'(a)') ' Covariance Matrix'
      j=0
      do 1 i=1,nsa
      if(isel(i).ne.0) then
      j=j+1
      ndex(j)=i
      endif
    1 continue
c
      return
      do 6 iproc=0,np-1
      call ga_distribution(ga_cov,iproc,il,ih,jl,jh)
      do 3 i=il,ih
      call ga_get(ga_cov,i,i,jl,jh,cov,1)
      ia=ndex((i+2)/3)
      ix=mod(i+2,3)+1
      do 4 j=jl,jh
      ja=ndex((j+2)/3)
      jx=mod(j+2,3)+1
      write(*,'(4i5,f12.6)') i,j,ia,ja,cov(j-jl+1)
    4 continue
    3 continue
    6 continue
c
      return
      end
