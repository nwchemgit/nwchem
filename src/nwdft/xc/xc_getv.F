      Subroutine xc_getv( rtdb, geom, AO_bas_han, 
     &     XC_bas_han, Exc,ipol,nExc,
     &     iExc,iVxc_opt,iga_xcinv,
     &     iga_dens,iga_vxc,IOLGC,iAOacc,
     &     iga_wght,iga_xyz,wght_GA)
C$Id: xc_getv.F,v 1.1 1995-04-26 03:40:13 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical LResult,IOLGC,wght_GA
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han,  XC_bas_han
c**** 


c
#include <mafdecls.h>
#include <rtdb.fh>
#include <bas.fh>
#include <global.fh>
#include <tcgmsg.fh>

c removed include bas_dft.fh
c
c--> Matrix Elements of the XC Potential & Energy
c
      Dimension iga_vxc(*), iga_dens(2)
c
c--> XC Energy
c
      Dimension Exc(nExc)
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
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
      Exc(1)=0.d0
      Exc(2)=0.d0
      me=ga_nodeid()
      nproc=ga_nnodes()
      zero = 0.0d0
      one = 1.0d0
      onem = -1.0d0
      iEx = iExc/10
      iEc = mod(iExc,10)
      if( .not. bas_numbf(AO_bas_han,nbf_ao) )then
        call errquit('Exiting in getvxc.',1)
      end if

      if( iEx.eq.0 )then
        
c       
c       Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
        if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1, itol2e))
     $       call errquit('fitcd: rtdb_get failed', 0)
        tol2e=10.d0**(-itol2e)
        call rhf_fock_2e(geom, AO_bas_han, iga_dens(1), iga_vxc(1), 
     &       tol2e,.false., .true., .false.)
        if(ipol.eq.2) then
          call rhf_fock_2e(geom, AO_bas_han, iga_dens(2), iga_vxc(2), 
     &         tol2e,.false., .true., .false.)
        endif
        call ga_sync
c       
c       Compute the exact exchange energy.
c       
        if( ipol.eq.1 )then
          Exc(1) = Exc(1)+0.5d0*ga_ddot(iga_dens(1),iga_vxc(1))
        else
          Exc(1) = Exc(1)+(ga_ddot(iga_dens(1),iga_vxc(1))+
     &         ga_ddot(iga_dens(2),iga_vxc(2)))
          call ga_dscal(iga_vxc(1),2.d0)
          call ga_dscal(iga_vxc(2),2.d0)
        end if
      endif
      if(iEx.ne.0.or.iEc.ne.0) then
        if( iVxc_opt.eq.0 )then
c         
c         Option A:  Compute via direct numerical quadrature.
c         
          call xc_quadv0(rtdb,geom, AO_bas_han, 
     &         ipol,nExc,iExc,
     &     iga_vxc,iga_dens,iAOacc,
     &     iga_wght,iga_xyz,wght_GA,
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
     &         Dbl_MB(iTmat),nbf_xc,ipol,nExc,iExc,
     &         iga_dens,Exc,iAOacc)
          call xc_fitv(AO_bas_han, XC_bas_han,geom,
     &         Dbl_MB(iTmat),nTrows,nTcols,
     &         1,iga_vxc,iga_xcinv,IOLGC)

          LResult = MA_Pop_Stack(lTmat)

        endif

      end if
      return
      end
