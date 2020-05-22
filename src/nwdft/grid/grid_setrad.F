      Subroutine grid_setrad(rtdb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "cdft.fh"
c
      integer rtdb
      double precision eps
      parameter (eps = 1.d-20)
      logical oprint, oprint_grid
      logical tol_reset, llind,check_s
c
      integer i, itype, ir, nr, irsh
      double precision r, r_max
c
      integer lexpo, iexpo
cc AJL/Begin
      integer nexpo
cc AJL/End
c
      integer n1,ninf,ndivide,mxnrad,znint,me,
     .     k_2l,l_2l,k_2ex,l_2ex,k_2g,l_2g
      double precision braggmult
      double precision midpoint,range_max,range_mk
      me=ga_nodeid()
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
c
      if (.not. rtdb_get(rtdb, 'dft:tol_reset', mt_log, 1,
     &   tol_reset)) tol_reset = .true.
c
      llind=.false.
      if (.not. rtdb_get(rtdb, 'dft:check_s', mt_log, 1,
     &   check_s))then
         check_s = .false.
      endif
      if(check_s) llind=.true.
c      check_s=
c
c     define radial grid type
c
      if (.not. rtdb_cget(rtdb,'dft:radgrid', 1,radgrid)) 
     .     call errquit('grid_setrad: radgrid not defined',9, RTDB_ERR)
      if(radgrid.eq.'eumac') then
         wradgrid=1
      elseif(radgrid.eq.'murakn') then
         wradgrid=2         
      elseif(radgrid.eq.'tarad') then
         wradgrid=3
      elseif(radgrid.eq.'gaussleg') then
         wradgrid=4
      elseif(radgrid.eq.'gcinterv') then
         wradgrid=5
      elseif(radgrid.eq.'chebyshev') then
         wradgrid=7
      elseif(radgrid.eq.'legendre') then
         wradgrid=8
      else
         call errquit('grid_setrad: unknown radial grid type', 0,
     &       INPUT_ERR)
      endif
c
c     check radial cutoffs
c
      if(llind) then
         if (.not.MA_Push_Get(mt_int,ntypes,'l+l',l_2l,k_2l))
     &        call errquit('grid_setrad: cannot allocate 2l',0, MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,ntypes,'exex',l_2ex,k_2ex))
     &        call errquit('grid_setrad: cannot allocate 2e',1, MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,ntypes,'exex',l_2g,k_2g))
     &        call errquit('grid_setrad: cannot allocate 2e',2, MA_ERR)
      else
         k_2l=0
         k_2ex=0
         k_2g=0
      endif

cc AJL/Begin/FDE
      nexpo = nbf_ao_mxprim
      if (frozemb_fde) nexpo = max(nexpo, nbf_fde_mxprim)
c      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
      if (.not.MA_Push_Get(MT_Dbl,nexpo,'expo',lexpo,iexpo))
cc AJL/End
     &   call errquit('grid_setrad: cannot allocate expo',0, MA_ERR)
#if 0
      if(tol_reset) call grid_rtrunc(rtdb, dbl_mb(iexpo), 
     .     llind,dbl_mb(k_2ex),dbl_mb(k_2g),int_mb(k_2l),
     &     tol_reset, range_max,check_s)
#else
      call grid_rtrunc(rtdb, dbl_mb(iexpo),
     .     llind,dbl_mb(k_2ex),dbl_mb(k_2g),int_mb(k_2l),
     &     .true., range_max,check_s)
#endif
      if (.not.ma_pop_stack(lexpo))
     &   call errquit('grid_setrad: cannot pop stack expo',1, MA_ERR)

c
c     Define the Euler-Maclaurin radial quadrature points and weights.
c
      if (.not. rtdb_get(rtdb,'dft:braggmult', mt_dbl, 1, braggmult)) 
     *     braggmult=2d0
      if (.not. rtdb_get(rtdb,'dft:ndivide', mt_int, 1, ndivide)) 
     *     ndivide=11

      do 100 itype = 1, ntypes
c
c        debug writes
c
         if (ga_nodeid().eq.0.and.oprint_grid)then
            write(LuOut,*)
     .           ' itype, bsrad_atom_type(itype), nrad(itype): ',
     &                  itype, bsrad_atom_type(itype), nrad(itype)
         endif
         if (bsrad_atom_type(itype).gt.EPS)then
           if(wradgrid.eq.1) then
              call eumac(bsrad_atom_type(itype), rpts(1,itype), 
     &             rwghts(1,itype), nrad(itype), 2d0)
           elseif(wradgrid.eq.6) then
               if(tol_reset) then
                  range_mk=range_max
               else
                  range_mk=dble(rad_cutoff(1,itype))
               endif
               if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',
     .              lexpo,iexpo))
     &              call errquit('grid_setrad: cant alloc expo',0,
     &       MA_ERR)
            if (.not.ma_pop_stack(lexpo))
     &           call errquit('grid_setrad: cannot pop stack',2, MA_ERR)

            elseif(wradgrid.eq.3) then
               if(tol_reset) then
                range_mk=dble(rad_cutoff(1,itype))
                  call xc_tarad(int(znuc_atom_type(itype)), 
     ,                 rpts(1,itype), rwghts(1,itype), 
     &                 nrad(itype),0d0,'chebyshr',range_mk)
               else
                  call xc_tarad(int(znuc_atom_type(itype)), 
     ,                 rpts(1,itype), rwghts(1,itype), 
     &                 nrad(itype),0d0,'chebyshe',0d0)
               endif
            elseif(radgrid.eq.'gausslag') then
               call grid_gausslag(nrad(itype), rpts(1,itype), 
     &              rwghts(1,itype),0d0)
          elseif(wradgrid.eq.5) then
c gcgc
            n1=nrad(itype)*ndivide/16
c            n1=nrad(itype)*11/8
            midpoint=bsrad_atom_type(itype)*braggmult
c           midpoint=bsrad_atom_type(itype)*2
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'chebyshe')
            call xc_tarad(iatype(itype), rpts(n1+1,itype), 
     &       rwghts(n1+1,itype), ninf,midpoint,'chebyshe',0d0)
c           call grid_gausslag(ninf, rpts(n1+1,itype), rwghts(n1+1,itype),
c     ,          midpoint)
            write(LuOut,*) ' itype midpoint n1',itype,midpoint,n1
          elseif(wradgrid.eq.7) then
c
c     Chebyshev with endpoint defined from grid_rtrunc
c
             range_mk=dble(rad_cutoff(1,itype))
             n1=nrad(itype)
             call xc_interv(0d0,range_mk,
     ,            rpts(1,itype), rwghts(1,itype),n1,'chebyshe')
          elseif(wradgrid.eq.8) then
c
c     Legendre with endpoint defined from grid_rtrunc
c
             range_mk=dble(rad_cutoff(1,itype))
             n1=nrad(itype)
             call xc_interv(0d0,range_mk,
     ,            rpts(1,itype), rwghts(1,itype),n1,'legendre')
          elseif(radgrid.eq.'leglag') then
            n1=nrad(itype)*ndivide/16
            midpoint=bsrad_atom_type(itype)*braggmult
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'legendre')
           call grid_gausslag(ninf,rpts(n1+1,itype),
     ,          rwghts(n1+1,itype),midpoint)
            write(LuOut,*) ' itype midpoint n1',itype,midpoint,n1
          elseif(radgrid.eq.'gclag') then
            n1=nrad(itype)*ndivide/16
            midpoint=bsrad_atom_type(itype)*braggmult
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'chebyshe')
           call grid_gausslag(ninf, rpts(n1+1,itype), 
     ,          rwghts(n1+1,itype),midpoint)
            write(LuOut,*) ' itype midpoint n1',itype,midpoint,n1
          elseif(radgrid.eq.'leggc') then
            n1=nrad(itype)*ndivide/16
            midpoint=bsrad_atom_type(itype)*braggmult
            ninf=nrad(itype)-n1
           call xc_interv(0d0,midpoint,
     ,           rpts(1,itype), rwghts(1,itype),n1,'legendre')
            call xc_tarad(iatype(itype), rpts(n1+1,itype), 
     &       rwghts(n1+1,itype), ninf,midpoint,'chebyshe',0d0)
            write(LuOut,*) ' itype midpoint n1',itype,midpoint,n1

          elseif(wradgrid.eq.2) then
             if(tol_reset) then
                range_mk=dble(rad_cutoff(1,itype))
             else
                range_mk=bsrad_atom_type(itype)
             endif
             znint=int(znuc_atom_type(itype))
cedo             if(ga_nodeid().eq.0.and.oprint) write(Luout,*)
cedo     .            ' Murakn itype cutoff ',itype,range_mk
             call xc_murakn(
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
     &       rwghts(n1+1,itype), ninf,midpoint,'legendre',0d0)
            write(LuOut,*) ' itype midpoint n1',itype,midpoint,n1
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

      if(llind) then
         if (.not.ma_pop_stack(l_2g))
     &        call errquit('grid_setrad: cannot pop stack',5, MA_ERR)
         if (.not.ma_pop_stack(l_2ex))
     &        call errquit('grid_setrad: cannot pop stack',3, MA_ERR)
         if (.not.ma_pop_stack(l_2l))
     &        call errquit('grid_setrad: cannot pop stack',4, MA_ERR)
      endif
c
c     we could have redefined max no rad pts
c
      mxnrad=0
      do itype=1,ntypes
         mxnrad=max(mxnrad,nrad(itype))
      enddo
      if (.not. rtdb_put(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &   call errquit('grid_setrad: rtdb_put failed', 119, RTDB_ERR)
c
c     Modify number of radial shells per atom type by possible 
c     truncation.
c
c
c     debug writes
c
c      if (ga_nodeid().eq.0.and.oprint_grid)then
c         call i_output(rad_cutoff, 1, 2, 1, ntypes, 2, ntypes, 1)
c      endif
      do 420 itype = 1, ntypes
c
         if (bsrad_atom_type(itype).le.eps) goto 420
c
c        NOTE - modification of r_max
c         
         r_max = dble(rad_cutoff(1,itype))
c
         nr = nrad(itype)

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
         write(LuOut,*)' nrad: ',(nrad(itype), itype = 1, ntypes)
         write(LuOut,*)' radgrid: ',radgrid
         write(LuOut,*)' tol_reset,  ', 
     &               tol_reset
      endif
      return
      end
