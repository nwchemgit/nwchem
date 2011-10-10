      subroutine argos_cafe_lambda(lamtyp,lambda,maxlam,explam,
     + lfnout,lfnpmf,rlambd,dlambd,projct)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lamtyp,lambda,maxlam,lfnout,lfnpmf
      real*8 explam,rlambd,dlambd
      character*80 projct
      character*1 ch
c
      character*255 filnam
c
      if(me.eq.0.and.lpmf) then
      if(lfnpmf.lt.0) then
      close(unit=-lfnpmf,status='keep')
      else
      lfnpmf=-lfnpmf
      endif
      write(filnam,'(a,a,i5.5,a)') projct(1:index(projct,' ')-1),'-',
     + lambda,'.pmf '
      open(unit=-lfnpmf,file=filnam(1:index(filnam,' ')-1),
     + form='formatted')
      if(mropt.eq.3) then
    1 continue
      read(-lfnpmf,1000,end=2) ch
 1000 format(a1)
      goto 1
    2 continue
      endif
      endif
c
      call argos_cafe_lam(lamtyp,lambda,maxlam,explam,
     + lfnout,rlambd,dlambd,
     + dbl_mb(i_mas),dbl_mb(i_vdw),dbl_mb(i_chg),
     + mbt(1),numb(1),mbp(1),dbl_mb(i_bnd(1)),
     + mht(1),numh(1),mhp(1),dbl_mb(i_ang(1)),
     + mdt(1),numd(1),mdp(1),dbl_mb(i_dih(1)),
     + mit(1),numi(1),mip(1),dbl_mb(i_imp(1)),
     + mbt(2),mbp(2),dbl_mb(i_bnd(2)),mht(2),mhp(2),dbl_mb(i_ang(2)),
     + mdt(2),mdp(2),dbl_mb(i_dih(2)),mit(2),mip(2),dbl_mb(i_imp(2)),
     + lfnpmf)
c
      if(ipme.ne.0) call argos_cafe_pmelam()
c
      return
      end
c $Id$
