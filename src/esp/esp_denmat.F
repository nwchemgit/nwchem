      subroutine esp_denmat(occ,ene)
c
c $Id: esp_denmat.F,v 1.7 2000-05-04 07:09:56 d3j191 Exp $
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "esp_params.fh"
#include "esp_common.fh"
#include "util.fh"
c  
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
c
      real*8 ene(nbf),occ(nbf)
c
      integer ilo,ihi,jlo,jhi,lthi,lthj,i
      integer lg_s,k_scr,l_scr,icol
c
      real*8 scale
      character*20 scftyp
c
c     allocate global arrays
c
      lg_d=ga_create_atom_blocked(igeom,ibasis,'density')
      lg_v=ga_create_atom_blocked(igeom,ibasis,'vector')
      lg_s=ga_create_atom_blocked(igeom,ibasis,'scr')
c
c     open movecs file and read header info
c
      call util_file_name('movecs',.false.,.false.,movfil)
      if(.not.movecs_read_header(movfil,title,basnam,scftyp,nbf,nset,
     + nmo,2)) call errquit('esp: movecs_read_header failed',me)
c
c     read orbital occupations, energies and vectors from movecs
c
      if(.not.movecs_read(movfil,1,occ,ene,lg_v))
     + call errquit('eps: movecs_read failed',me)
c
      nocc=0
      charge=chgnuc
      do 1 i=1,nbf
      if(occ(i).gt.small) charge=charge-occ(i)
    1 continue
c
      call ga_copy(lg_v,lg_s)
      call ga_sync()
c
      call ga_distribution(lg_s,me,ilo,ihi,jlo,jhi)
      lthi=ihi-ilo+1
      lthj=jhi-jlo+1
c
      if(lthi.ne.0.and.lthj.ne.0) then
      if(lthi.lt.0.or.lthj.lt.0) 
     + call errquit('esp_denmat: distribution error',me)
c
      if(.not.ma_push_get(mt_dbl,lthi,'vector',l_scr,k_scr))
     + call errquit('esp_denmat: failed to allocate vector',me)
c
      do 2 icol=jlo,jhi
      call ga_get(lg_s,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx
cx      write(*,6666) icol,(dbl_mb(k_scr+i),i=0,nbf-1)
cx 6666 format(//,'Vector',i5,/,(10f10.6))
cx
      scale=occ(icol)
cx
cx      write(*,'(a,i5,f12.6)') 'Scale ',icol,scale
cx
      call dscal(lthi,scale,dbl_mb(k_scr),1)
      call ga_put(lg_s,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx      call ga_get(lg_s,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx
cx      write(*,6667) (dbl_mb(k_scr+i),i=0,nbf-1)
cx 6667 format(10f10.6)
cx
    2 continue
c
      if(.not.ma_pop_stack(l_scr))
     + call errquit('esp_denmat: failed to deallocate vector',me)
c
      endif
c
      call ga_sync()
      call ga_dgemm('n','t',nbf,nbf,nbf,1.d0,lg_v,lg_s,0.d0,lg_d)
c
c     deallocate global arrays lg_s, lg_v and lg_d
c
      if(.not.ga_destroy(lg_s))
     + call errquit('esp_denmat: failed to destroy lg_s',me)
      if(.not.ga_destroy(lg_v))
     + call errquit('esp_denmat: ga_destroy lg_v failed',me)
c
      return
      end
