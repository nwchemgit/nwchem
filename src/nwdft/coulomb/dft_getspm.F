      Subroutine dft_getspm(geom,lmax,iga_dens,AO_bas_han)

C$Id: dft_getspm.F,v 1.1 1995-05-03 19:27:50 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      Logical LResult
      
c**** 
c**** nwchem handles
c**** 
      integer  geom, AO_bas_han

#include <bas.fh>
#include <geom.fh>

#include <mafdecls.h>
#include <stdio.fh>


        if (.not. geom_ncent(geom, natoms))
     &       call errquit('dft_scf: geom_ncent failed',73)      

        if(.not.Ma_Push_Get(MT_Dbl,natoms*3,'coordinates',lcoord,
     &       icoord))
     &       call errquit('dft_scf: push_get failed',74)
        if(.not.Ma_Push_Get(MT_Dbl,natoms,'charges',lcharge,
     &       icharge))
     &       call errquit('dft_scf: push_get failed',74)
        if(.not.Ma_Push_Get(MT_Byte, natoms*16, 'center tags',
     &       ltags, itags))
     &       call errquit('dft_scf: push_get failed',74)

        if (.not. geom_cart_get(geom, natoms, Byte_MB(itags),
     &       Dbl_MB(icoord), Dbl_MB(icharge)))
     &       call errquit('dft_scf: geom_cart_get failed',74)

        LResult = MA_Pop_Stack(ltags)

        if(lmax.gt.0) then
          if(.not.Ma_Push_Get(MT_Dbl,3*natoms,'dipole',ldipole,
     &         idipole))
     &         call errquit('dft_scf: push_get failed',74)
        endif

        max_at_bf = 0
        do iatom = 1, natoms
          if (.not. bas_ce2bfr(AO_bas_han, iatom, ilo, ihi))
     $         call errquit('quadvxc0: bas_ce2bfr failed', iatom)
          max_at_bf = max(max_at_bf, ihi-ilo+1)
        enddo
        max_at_bf2=max_at_bf*max_at_bf

        if(.not.Ma_Push_Get(MT_Dbl,max_at_bf2,'local DM',lPmat,
     &       iPmat))
     &       call errquit('dft_scf: push_get failed',74)

        call dft_genspm(lmax,iga_dens,AO_bas_han,
     &       natoms,dbl_MB(icoord),dbl_MB(iPmat),max_at_bf,
     &       DBL_MB(icharge),DBL_MB(idipole))

        LResult = MA_Pop_Stack(lPmat)
        if(lmax.gt.0) then
          LResult = Lresult.and.MA_Pop_Stack(ldipole)
        endif
        LResult = Lresult.and.MA_Pop_Stack(lcharge)
        LResult = Lresult.and.MA_Pop_Stack(lcoord)
        if(.not.Lresult) call errquit('dft_scf: pop_stack failed',74)

      return
      end
