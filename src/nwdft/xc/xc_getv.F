      Subroutine xc_getv( rtdb, geom, AO_bas_han, 
     &     XC_bas_han, Exc,ipol,nExc,
     &     iVxc_opt,g_xcinv,
     &     g_dens,g_vxc,IOLGC,iAOacc,
     &     g_wght,g_xyz,wght_GA)
C$Id: xc_getv.F,v 1.4 1995-12-19 20:21:54 og845 Exp $
      implicit none
      
      integer ipol,nExc
      integer iVxc_opt
      integer g_xcinv,g_dens(*),g_vxc(*),g_wght,g_xyz
      integer iAOacc ! [input]

      Logical LResult,IOLGC,wght_GA
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han,  XC_bas_han
c**** 


c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"

      double precision jfac(4),kfac(4),xfac(10),cfac(10)
      integer d(4),f(4)
c
c--> XC Energy
c
      double precision Exc(nExc)
c
c This driver routine solves for the XC energy and potential (Vxc) via 
c numerical quadrature methods. The results are obtained either by direct 
c numerical integration or by means of a LSQ fit of the Vxc to a set of 
c Gaussian functions. This fitted function can be used to evaluate Vxc
c via a summation of a series of 3-center overlap integrals (3OIs). The
c algorithms are formulated in terms of matrix products. See subsequent
c subroutines for further explanation.
c 
c              XC Energy and Potential Index Key, Vxc(pq,i)
c
c              Value of     |     Definition of index "i"
c            ipol     nExc  |    1        2        3       4
c           --------------------------------------------------
c              1        1   |   Vxc      
c              2        1   |   Vxc^up   Vxc^dw   
c              1        2   |   Vxc      
c              2        2   |   Vxc^up   Vxc^dw   
c
c           nTcols = ipol 
c
c
      integer g_scr,me,nproc,i,nTrows,nbf_xc,nbf_ao,nTcols
      integer lTmat,iTmat,itol2e
      double precision zero,one,onem
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e,tot
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
      Exc(1)=0.d0
      Exc(2)=0.d0
      me=ga_nodeid()
      nproc=ga_nnodes()
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, 10, cfac))
     $     call errquit('xc_getv: rtdb_put failed', 0)
      if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, 10, xfac))
     $     call errquit('xc_getv:rtdb_put failed', 0)
      if( .not. bas_numbf(AO_bas_han,nbf_ao) )then
        call errquit('Exiting in getvxc.',1)
      end if

      if( abs(xfac(1)).gt.1e-8 )then
        
c       
c       Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
        if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1, itol2e))
     $       call errquit('fitcd: rtdb_get failed', 0)
        tol2e=10.d0**(-itol2e)
        if(ipol.eq.1) then
          jfac(1)=0.
          kfac(1) = -0.5d0*xfac(1)
          call fock_2e(geom, AO_bas_han, 1,
     &         jfac,kfac,
     &         tol2e,.false.,
     &         g_dens(1), g_vxc(1))
        else
C
C         bogus Ga
C
          call gacrea(geom,AO_bas_han,
     &         nbf_ao,nbf_ao,'ga scratch',g_scr,'atom')
          call ga_zero(g_scr)
          jfac(1)=0.d0
          jfac(2)=0.d0
          jfac(3)=0.d0
          jfac(4)=0.d0
          kfac(1)=0.0d0
          kfac(2)=-1.0d0*xfac(1)
          kfac(3)=0.0d0
          kfac(4)=-1.0d0*xfac(1)
          d(1) = g_dens(1)
          d(2) = g_dens(1)
          d(3) = g_dens(2)
          d(4) = g_dens(2)
C          f(1) = g_a_coul
          f(1) = g_scr
          f(2) = g_vxc(1)
C          f(3) = g_b_coul
          f(3) = g_scr
          f(4) = g_vxc(2)
          call fock_2e(geom, AO_bas_han,4, 
     &         jfac,kfac,
     &         tol2e,.false.,
     &         d, f)
          call gadest(g_scr)
        endif
        call ga_sync
c       
c       Compute the exact exchange energy.
c       
        if( ipol.eq.1 )then
          Exc(1) = Exc(1)+0.5d0*
     *         ga_ddot(g_dens(1),g_vxc(1))
        else
          Exc(1) = Exc(1)+0.5d0*
     *         (ga_ddot(g_dens(1),g_vxc(1))+
     &         ga_ddot(g_dens(2),g_vxc(2)))
        end if
      endif
      tot=-xfac(1)
      do i=1,10
        tot=tot+xfac(i)+cfac(i)
      enddo
      if(abs(tot).gt.1e-8) then
        if( iVxc_opt.eq.0 )then
c         
c         Option A:  Compute via direct numerical quadrature.
c         
          call xc_quadv0(rtdb,geom, AO_bas_han, 
     &         ipol,nExc,
     &     g_vxc,g_dens,iAOacc,
     &     g_wght,g_xyz,wght_GA,
     &     Exc)
 

        else if( iVxc_opt.eq.1 )then
c         
c         Option B:  Compute via a LSQ fitting procedure.
c         
c         Define the dimensions of the "T" matrix.
c         
          if( .not. bas_numbf(XC_bas_han,nbf_xc) )then
            call errquit('Exiting in getvxc.',1)
          end if

          nTrows = nbf_xc
          nTcols = ipol
c         
c         Allocate scratch space for the "T" matrix.
c         
          LResult = 
     &         MA_Push_Get(MT_Dbl,nTrows*nTcols,'Tmat',lTmat,iTmat)
          call xc_quadv1(rtdb,geom, AO_bas_han, XC_bas_han,
     &         Dbl_MB(iTmat),nbf_xc,ipol,nExc,
     &         g_dens,Exc,iAOacc)
          call xc_fitv(AO_bas_han, XC_bas_han,geom,
     &         Dbl_MB(iTmat),nTrows,nTcols,
     &         1,g_vxc,g_xcinv,IOLGC)

          LResult = MA_Pop_Stack(lTmat)

        endif

      end if
      return
      end
