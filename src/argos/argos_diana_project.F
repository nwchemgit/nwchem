      subroutine argos_diana_project()
c
      implicit none
c
#include "argos_diana_common.fh"
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
      call argos_diana_proj(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
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
c      call argos_diana_wthdr(lfnprj,fmt,byte_mb(i_snam),byte_mb(i_tag),
c     + int_mb(i_isel),.false.)
c
      call ga_sync()
c
      call argos_diana_wrttrj(ga_prj,fil(1:i-1),fmt,nbatch,.false.)
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
c $Id$
