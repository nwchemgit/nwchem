      Subroutine grid_setrad(rtdb)
c
c$Id: grid_setrad.F,v 1.10 1999-12-30 21:19:27 edo Exp $
c
      implicit none
c      
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "stdio.fh"
#include "cdft.fh"
c
      integer rtdb
      double precision eps
      parameter (eps = 1.d-20)
      logical oprint, oprint_grid
      logical tol_reset, mul_reset
c
      character*10 radgrid
      character*10 radgrids(6)
      data radgrids /'eumac','murakn','tarad','gausslag','gcinterv',
     ,'rlindh'/
      integer i, itype, ir, nr, irsh
      double precision r, r_max
c
      integer lexpo, iexpo, insh, lnsh
c
      integer n1,ninf,ndivide,mxnrad,znint
      double precision zz
      double precision braggmult
      double precision midpoint,range_max,range_mk
      double precision lindh_acc
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
c
c     check radial cutoffs
c
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_setrad: cannot allocate expo',0)
      if (.not.MA_Push_Get(MT_Dbl,nshells_ao,'nsh',lnsh,insh))
     &   call errquit('grid_setrad: cannot allocate nsh',0)
      if (.not. rtdb_get(rtdb, 'dft:lindh_acc', mt_dbl, 
     ,     1, lindh_acc)) lindh_acc=1d-10
      if (.not. rtdb_get(rtdb, 'dft:tol_reset', mt_log, 1,
     &   tol_reset))then
         tol_reset = .false.
      endif
      if (.not. rtdb_get(rtdb, 'dft:mul_reset', mt_log, 1,
     &   mul_reset))then
         mul_reset = .false.
      endif
      if(tol_reset)call rad_trunc(rtdb, dbl_mb(iexpo), dbl_mb(insh), 
     &     tol_reset, mul_reset,range_max)
      if (.not.ma_pop_stack(lnsh))
     &   call errquit('grid_setrad: cannot pop stack',0)
      if (.not.ma_pop_stack(lexpo))
     &   call errquit('grid_setrad: cannot pop stack',0)
c
c     Define the Euler-Maclaurin radial quadrature points and weights.
c
      if (.not. rtdb_get(rtdb,'dft:braggmult', mt_dbl, 1, braggmult)) 
     *     braggmult=2d0
      if (.not. rtdb_cget(rtdb,'dft:radgrid', 1,radgrid)) 
     *     radgrid='eumac'
      if (.not. rtdb_get(rtdb,'dft:ndivide', mt_int, 1, ndivide)) 
     *     ndivide=11
      do 100 itype = 1, ntypes
c
c        debug writes
c
         if (ga_nodeid().eq.0.and.oprint_grid)then
            write(*,*)' itype, bsrad_atom_type(itype), nrad(itype): ',
     &                  itype, bsrad_atom_type(itype), nrad(itype)
         endif
         if (bsrad_atom_type(itype).gt.EPS)then
           if(radgrid.eq.'eumac') then
            call eumac(bsrad_atom_type(itype), rpts(1,itype), 
     &            rwghts(1,itype), nrad(itype), rm)
         elseif(radgrid.eq.'rlindh') then
               if(tol_reset) then
                  range_mk=range_max
               else
                  range_mk=dble(rad_cutoff(1,itype))
               endif
               if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',
     .              lexpo,iexpo))
     &              call errquit('grid_setrad: cant alloc expo',0)
            call grid_roland(
     ,              rpts(1,itype),rwghts(1,itype),nrad(itype),
     ,              ityp2ctr(itype),dbl_mb(iexpo),
     ,              mxrad,rad_cutoff(1,itype),lindh_acc)
            if (.not.ma_pop_stack(lexpo))
     &           call errquit('grid_setrad: cannot pop stack',0)
            elseif(radgrid.eq.'tarad') then
            call xc_tarad(int(znuc_atom_type(itype)), 
     ,              rpts(1,itype), rwghts(1,itype), 
     &                 nrad(itype),0d0,'chebyshe')
          elseif(radgrid.eq.'gausslag') then
           call grid_gausslag(nrad(itype), rpts(1,itype), 
     &            rwghts(1,itype),0d0)
          elseif(radgrid.eq.'gcinterv') then
c gcgc
            n1=nrad(itype)*ndivide/16
c            n1=nrad(itype)*11/8
            midpoint=bsrad_atom_type(itype)*braggmult
c           midpoint=bsrad_atom_type(itype)*2
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'chebyshe')
            call xc_tarad(iatype(itype), rpts(n1+1,itype), 
     &       rwghts(n1+1,itype), ninf,midpoint,'chebyshe')
c           call grid_gausslag(ninf, rpts(n1+1,itype), rwghts(n1+1,itype),
c     ,          midpoint)
            write(6,*) ' itype midpoint n1',itype,midpoint,n1
          elseif(radgrid.eq.'leglag') then
            n1=nrad(itype)*ndivide/16
            midpoint=bsrad_atom_type(itype)*braggmult
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'legendre')
           call grid_gausslag(ninf,rpts(n1+1,itype),
     ,          rwghts(n1+1,itype),midpoint)
            write(6,*) ' itype midpoint n1',itype,midpoint,n1
          elseif(radgrid.eq.'gclag') then
            n1=nrad(itype)*ndivide/16
            midpoint=bsrad_atom_type(itype)*braggmult
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'chebyshe')
           call grid_gausslag(ninf, rpts(n1+1,itype), 
     ,          rwghts(n1+1,itype),midpoint)
            write(6,*) ' itype midpoint n1',itype,midpoint,n1
          elseif(radgrid.eq.'leggc') then
            n1=nrad(itype)*ndivide/16
            midpoint=bsrad_atom_type(itype)*braggmult
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'legendre')
            call xc_tarad(iatype(itype), rpts(n1+1,itype), 
     &       rwghts(n1+1,itype), ninf,midpoint,'chebyshe')
            write(6,*) ' itype midpoint n1',itype,midpoint,n1
          elseif(radgrid.eq.'murakn') then
             if(tol_reset) then
                range_mk=range_max
             else
                range_mk=dble(rad_cutoff(1,itype))
             endif
             znint=int(znuc_atom_type(itype))
             call xc_murakn(
     .            znint,
     .            rpts(1,itype), 
     &            rwghts(1,itype),nrad(itype),range_mk,
     .            tol_reset)
          elseif(radgrid.eq.'legleg') then
            n1=nrad(itype)*ndivide/16
            midpoint=bsrad_atom_type(itype)*braggmult
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'legendre')
            call xc_tarad(iatype(itype), rpts(n1+1,itype), 
     &       rwghts(n1+1,itype), ninf,midpoint,'legendre')
            write(6,*) ' itype midpoint n1',itype,midpoint,n1
         else
           write(0,*) ' radgrid eq ',radgrid
              call errquit('grid_setrad: unknown radial grid type', 0)
         endif
            if (ga_nodeid().eq.0.and.oprint)then
               write(LuOut,*)' Center: ',itype,' Radial Pts. & Wghts. '
               do 90 i = 1, nrad(itype)
                  write(LuOut,'(i5,2x,f17.9,2x,1pe22.7)') 
     &            i, rpts(i,itype), rwghts(i,itype)
   90          continue
            endif
         endif
  100 continue
c
c     we could have redefined max no rad pts
c
      mxnrad=0
      do itype=1,ntypes
         mxnrad=max(mxnrad,nrad(itype))
      enddo
      if (.not. rtdb_put(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &   call errquit('grid_setrad: rtdb_put failed', 119)
c
c     Modify number of radial shells per atom type by possible 
c     truncation.
c
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         call i_output(rad_cutoff, 1, 2, 1, ntypes, 2, ntypes, 1)
      endif
      do 420 itype = 1, ntypes
c
         if (bsrad_atom_type(itype).le.eps) goto 420
c
c        NOTE - modification of r_max
c         
         r_max = dble(rad_cutoff(1,itype))
c
         nr = nrad(itype)
c
         irsh = nr
c
         do 410 ir = nr, 1, -1
            r = rpts(ir,itype)
            if (r.gt.R_MAX)then
c
               irsh = irsh - 1
c
            else
               nrad(itype) = irsh
               goto 420
            endif
  410    continue
  420 continue
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         write(6,*)' nrad: ',(nrad(itype), itype = 1, ntypes)
         write(6,*)' radgrid: ',radgrid
         write(6,*)' tol_reset, mul_reset: ', 
     &               tol_reset, mul_reset
      endif
      return
      end
