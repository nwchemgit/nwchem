      Subroutine dftg_getxc( rtdb, geom, AO_bas_han, CD_bas_han, 
     &                    XC_bas_han, ipol, natom, iEx, iEc,
     &                    iga_dens, force)
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

c****
c**** compute XC component of energy gradient.  
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han, XC_bas_han
c**** 
#include <mafdecls.h>
#include <rtdb.fh>
#include <bas.fh>
#include <global.fh>
#include <tcgmsg.fh>

c removed include bas_dft.fh
c
c--> Matrix Elements of the XC Potential & Energy
c
      Dimension iga_dens(2)
c
c--> XC Energy
c
      Dimension force(3*natom)
c
c****
c**** comments from getvxc:
c****
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
c
c****
c**** currently, only the "exact" quadrature (no fit of Coulomb or
c**** XC terms used) form of gradient available
c****
c**** Exc(1)=0.d0
c**** Exc(2)=0.d0
      me=ga_nodeid()
      nproc=ga_nnodes()
      zero = 0.0d0
      one = 1.0d0
      onem = -1.0d0
c**** iEx = iExc/10
c**** iEc = mod(iExc,10)
      if( .not. bas_numbf(AO_bas_han,nbf_ao) )then
        call errquit('Exiting in getdexc.',1)
      end if
c     
c     Option A:  Compute via direct numerical quadrature.
c     
c**** 
c**** This is the only option presently active
c**** 
      call dftg_quadv0(rtdb,geom, AO_bas_han, 
     &     ipol,iEx,iEc,
     &     iga_dens,force,natom) ! debug


      return
      end
