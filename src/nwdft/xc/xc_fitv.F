      Subroutine xc_fitv(rtdb,Tmat,nTrows,nTcols,g_vxc,g_xcinv,
     $                    IOLGC)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c
      integer nTrows, nTcols
      double precision Tmat(nTrows, nTcols)
      integer  g_vxc(4), g_xcinv
      Logical IOLGC
      integer rtdb                       ! Qin
c
#include "bas.fh"
#include "mafdecls.fh"
#include "rtdb.fh"                    
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "cdft.fh"
#include "sym.fh"
c
      Logical LResult
c
c Compute the fitted XC energy and potential. The fitting coefficients are
c obtained by solving the normal equations for a standard linear least
c squares problem.
c
c                          _
c               ~   i     \          /      
c              (Vxc)   =  /_   C     | chi  chi  chi  dr
c                   pq       s  s,i  /    p    q    s
c                        _
c                       \     -1
c               C    =  /_   S    T
c                s,i       r  s,r  r,i
c
c                       /       i
c               T    =  | chi  V    dr
c                r,i    /    r  xc
c
c                                INDICES
c                               ---------
c                      pq --> product of AO basis functions (1,...,ntri_ao)
c                  r or s --> XC LSQ fitting functions (1,...,nbf_xc)
c                       i --> energy and potential index (1,...,nTcols)
c
c
c                      Energy and Potential Index Key
c
c              Value of     |     Definition of index "i"
c            ipol     nExc  |    1        2        3       4
c           --------------------------------------------------
c              1        1   |   Vxc      Exc   
c              2        1   |   Vxc^up   Vxc^dw   Exc
c              1        2   |   Vxc      Ex       Ec
c              2        2   |   Vxc^up   Vxc^dw   Ex       Ec    
c
c           nTcols = ipol + nExc
c
c
c Qin
      integer  nxtask
      external nxtask
      integer me, nproc, MXnshbf_ao, MXnshbf_XC
      integer ndbl
      integer lxc_coef, ixc_coef, lscr, iscr, lmiss
      integer ia, ii, j, k, next, istart, n3OI, l3OI, i3OI, lvvxc, 
     &        ivvxc, lvec2, ivec2, lscr2, iscr2
      integer n1, n2, ishb, ifirstb, ilastb, nshbfb,
     &        ishc, ifirstc, ilastc, nshbfc, Nao2, iixc,
     &        isha, ifirsta, ilasta, nshbfa, nscr2
      integer Nintegrals, npint, mag, Max_scr, nxc,nint_api
      double precision fact, q1,toll
c******************************************************************************
      toll=10d0**(-itol2e)
      me=ga_nodeid()
      nproc = ga_nnodes()
c
c     Determine the characteristics of the AO and XC Gaussian basis sets.
c

      LResult = .true.
      LResult = LResult .and. bas_numcont( AO_bas_han, nshells_ao) 
      LResult = LResult .and. bas_numbf( AO_bas_han, nbf_ao)
      LResult = LResult .and. bas_nbf_cn_max( AO_bas_han, MXnshbf_ao)
      LResult = LResult .and. bas_numcont( XC_bas_han, nshells_XC) 
      LResult = LResult .and. bas_numbf( XC_bas_han, nbf_XC)
      LResult = LResult .and. bas_nbf_cn_max( XC_bas_han, MXnshbf_XC)
      if ( .not. LResult) then
        call errquit('Exit from xc_fitv',1, UNKNOWN_ERR)
      endif 
c
c     Compute the XC fitting coefficients.
c
      ndbl = nbf_xc*nTcols
c
      if(.not.MA_Push_Get(MT_Dbl,ndbl,'xc_coef',lxc_coef,ixc_coef))
     &  call errquit('xc_fitv: cannot allocate xc_coef',0, MA_ERR)
      call dfill(ndbl,0.d0,DBL_MB(ixc_coef),1)
c
      if(.not.MA_Push_Get(MT_Dbl,nbf_xc,'scratchinv',lscr,iscr))
     &  call errquit('xc_fitv: cannot allocate scratchinv',0, MA_ERR)
c
c      if (IOLGC)call gacrea(geom, XC_bas_han, nbf_xc, nbf_xc, 
c     &                     'XCinv', g_xcinv, 'row')
      if (IOLGC)then
*ga:1:0
        if (.not. ga_create(mt_dbl, nbf_xc, nbf_xc, 'XCinv', 
     &                      nbf_xc, 0, g_xcinv)) 
     &     call errquit('xc_fitv: error creating g_xcinv',0, GA_ERR)
      endif
c
      if(me.eq.0.and.IOLGC) then
        if (ADFT) then
          call dft_invio('CDI',g_xcinv,nbf_xc,'READ',lmiss)
        else
          call dft_invio('XC',g_xcinv,nbf_xc,'READ',lmiss)
        endif
      endif
      call ga_sync
c
      next=nxtask(nproc, 1)+1
c
c
      do k=1,nbf_xc
         if (k.eq.next) then
            next=nxtask(nproc, 1)+1
c            call ga_get(g_xcinv, 1, nbf_xc, k, k, Dbl_MB(iscr), 1)
            call get_col(g_xcinv,nbf_xc,k,Dbl_MB(iscr))
            istart=ixc_coef
            do j = 1, nTcols
               call daxpy(nbf_xc, Tmat(k,j), DBL_MB(iscr), 1, 
     &                    DBL_MB(istart), 1)
               istart=istart+nbf_xc
            enddo
         endif
      enddo
c
      if (.not.ma_pop_stack(lscr))
     &   call errquit('xc_fitv: cannot pop stack',0, MA_ERR)
c
      next=nxtask(-nproc, 1)
c
      call ga_dgop(Msg_XC_Coef, DBL_MB(ixc_coef), ndbl, '+')

c
c      if (IOLGC)call gadest(g_xcinv)
      if (IOLGC)then
         if (.not. ga_destroy(g_xcinv)) call errquit
     &      ('xc_fitv: could not destroy g_xcinv', 0, GA_ERR)
      endif

      ! Copy the XC fitting coefficients to the beginning of
      ! Tmat and exit
      if (ADFT) then
        call dcopy(nTrows*nTcols,dbl_mb(ixc_coef),1,tmat,1)
        if(.not.ma_pop_stack(lxc_coef))
     &    call errquit('xc_fitv: cannot pop stack',0, MA_ERR)
        return
      endif
c
c
c     Compute the fitted Vxc.
c
c     Allocate scratch.
c
      n3OI = (MXnshbf_ao**2)*MXnshbf_xc
c
      call int_mem_3ov(mag,Max_scr)
c
      if (.not.MA_Push_Get(MT_Dbl, mag, '3OI', l3OI, i3OI))
     &   call errquit('xc_fitv: cannot allocate 3OI',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl, nTcols*mag, 'vvxc', lvvxc, ivvxc))
     &   call errquit('xc_fitv: cannot allocate vvxc',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl, mag, 'vec2', lvec2, ivec2))
     &   call errquit('xc_fitv: cannot allocate vec2',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl, Max_scr, 'scr2', lscr2, iscr2))
     &   call errquit('xc_fitv: cannot allocate scr2',0, MA_ERR)
c
c     Compute 3-ctr OIs between a given pair of AO shells and a set of XC shells.
c     
      call ga_sync
c
      n1 = nxtask(nproc, 1)
      n2 = 0
c
      do 240 ishb = 1, nshells_ao
         if (.not. bas_cn2bfr(AO_bas_han, ishb, ifirstb, ilastb))
     &       call errquit('Exiting in fitvxc.',10, BASIS_ERR)
         nshbfb=ilastb-ifirstb+1
c
         do 230 ishc = 1, ishb
c
            if (n1 .eq. n2)then
               if (.not. bas_cn2bfr(AO_bas_han, ishc, ifirstc, ilastc))
     &            call errquit('Exiting in fitvxc.',10, BASIS_ERR)
               nshbfc=ilastc-ifirstc+1
               Nao2 = nshbfb*nshbfc
               iixc=0
c
               do ii=0,nTcols*Nao2-1
                  DBL_MB(ivvxc+ii)=0.d0
               enddo
c
               do 210 isha = 1, nshells_xc
                  if (.not.bas_cn2bfr(XC_bas_han, isha, ifirsta, 
     &               ilasta))call errquit('Exiting in fitvxc.',12,
     &       BASIS_ERR)
                  nshbfa=ilasta-ifirsta+1
                  Nintegrals = Nao2*nshbfa
                  nint_api=mag
                  nscr2=max_scr
c                   
c                 check for use of symmetry
c
                  if (oskel) then
                     if (sym_shell(xc_bas_han, isha, q1))then
                        call int_1e3ov(XC_bas_han, isha, AO_bas_han, 
     &                       ishb, AO_bas_han, ishc, nscr2, 
     &                       Dbl_MB(iscr2), Nint_api, Dbl_MB(i3OI))
                        call dscal(Nintegrals, q1, Dbl_MB(i3OI), 1)
c      
c                       Multiply the matrix of 3-ctr OIs by the XC 
c                       expansion coefficients.
c
                        do ii = 0, nTcols-1
                           nxc=ixc_coef+ii*nbf_xc+iixc
                           npint=i3OI
                           istart=ivvxc+ii*Nao2
                           do ia = 0, nshbfa-1
                              fact=DBL_MB(nxc+ia)
                              if(abs(fact).gt.toll)
     /                         call daxpy(Nao2, fact, DBL_MB(npint), 1, 
     &                             DBL_MB(istart), 1)
                              npint=npint+Nao2
                           enddo
                        enddo
                     else
                        call dcopy(nintegrals,0d0,0,Dbl_MB(i3OI), 1)
                     endif
                  else
                     call int_1e3ov(XC_bas_han, isha, AO_bas_han, ishb,
     &                    AO_bas_han, ishc, nscr2, Dbl_MB(iscr2),
     &                    Nint_api, Dbl_MB(i3OI))
c      
c                    Multiply the matrix of 3-ctr OIs by the XC 
c                    expansion coefficients.
c
                     do ii = 0, nTcols-1
                        nxc=ixc_coef+ii*nbf_xc+iixc
                        npint=i3OI
                        istart=ivvxc+ii*Nao2
                        do ia = 0, nshbfa-1
                           fact=DBL_MB(nxc+ia)
                           if(abs(fact).gt.toll)
     .                          call daxpy(Nao2, fact, DBL_MB(npint), 1, 
     &                          DBL_MB(istart), 1)
                           npint=npint+Nao2
                        enddo
                     enddo
                  endif
                  iixc=iixc+nshbfa
  210          continue
c      
c              Scatter the computed values into Vxc.
c      
               do 220 k = 1, nTcols
                  call sca_lab_mat(AO_bas_han, ishb, ishc, nshbfb, 
     &                 nshbfc, g_vxc(k), DBL_MB(ivvxc+(k-1)*Nao2),
     &                 dbl_mb(ivec2),'acc')
  220          continue
c              
               n1 = nxtask(nproc, 1)
c              
            endif
c            
            n2 = n2 + 1
c            
  230    continue
  240 continue
      call ga_sync
c
      if(.not.ma_pop_stack(lscr2))
     &  call errquit('xc_fitv: cannot pop stack',0, MA_ERR)
c
      if(.not.ma_pop_stack(lvec2))
     &  call errquit('xc_fitv: cannot pop stack',0, MA_ERR)
c
      if(.not.ma_pop_stack(lvvxc))
     &  call errquit('xc_fitv: cannot pop stack',0, MA_ERR)
c
      if(.not.ma_pop_stack(l3OI))
     &  call errquit('xc_fitv: cannot pop stack',0, MA_ERR)
c
      if(.not.ma_pop_stack(lxc_coef))
     &  call errquit('xc_fitv: cannot pop stack',0, MA_ERR)
c
      n1 = nxtask(-nproc, 1)
c
      return
      end
