      subroutine argos_diana_covar()
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
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
      if(active) call argos_diana_edfinal()
      call argos_diana_edinit()
      active=.true.
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
c      call argos_diana_xaver(dbl_mb(ndx),nsa,min((me+1)*nfrme,nfrtot)-me*nfrdim,
c     + dbl_mb(i_xadat),nsa)
c
      call argos_diana_cov(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
     + dbl_mb(ndx),dbl_mb(i_xadat))
c
      call ga_scale(ga_cov,1.0d0/dble(nfrtot))
c
c      call argos_diana_pcovar(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel))
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
