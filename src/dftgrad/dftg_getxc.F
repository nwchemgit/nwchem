      Subroutine dftg_getxc( rtdb, geom, AO_bas_han, CD_bas_han, 
     &                    XC_bas_han, ipol, natom, 
     &                    iga_dens, force,
     &                    oactive, nactive)
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
c****
c**** compute XC component of energy gradient.  
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han, XC_bas_han
c**** 
#include <mafdecls.fh>
#include <rtdb.fh>
#include <bas.fh>
#include <global.fh>
#include <tcgmsg.fh>
      
      parameter(zero = 0.0d0,one = 1.0d0,onem = -1.0d0)

      integer nactive
      logical status
      logical oactive(natom)

      double precision xfac(10),cfac(10)
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
      nproc=ga_nnodes()
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, 10, cfac))
     $     call errquit('xc_getv: rtdb_put failed', 0)
      if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, 10, xfac))
     $     call errquit('xc_getv:rtdb_put failed', 0)
      if( .not. bas_numbf(AO_bas_han,nbf_ao) )then
        call errquit('Exiting in getdexc.',1)
      end if
      if( abs(xfac(1)).gt.1e-8 )then
        
c       
c       Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
        if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1, itol2e))
     $       call errquit('fitcd: rtdb_get failed', 0)
        tol2e=10.d0**(-itol2e)
        status = bas_nbf_ce_max(AO_bas_han,max_at_bf)
        call int_mem ( max1e, max2e, mscratch_1e, mscratch_2e )
        lbuf = max ( max1e, max2e )
        lscratch = max ( mscratch_1e, mscratch_2e )
        status = MA_Push_get(MT_DBL, lbuf, 'deriv buffer', l_buf,
     &     k_buf)

        if (.not.status) call errquit('could not allocate buffer',1)
        status = MA_Push_get(MT_DBL, lscratch, 'deriv scratch', 
     &       l_scr, k_scr)

C     allocate local density matrix block

        lsqa = max_at_bf * max_at_bf

        status = status.and.MA_Push_get(MT_DBL,lsqa,
     &       'local_w_density', l_ik, i_ik)
        status = status.and.MA_Push_get(MT_DBL,lsqa,
     &       'local_w_density', l_il, i_il)
        status = status.and.MA_Push_get(MT_DBL,lsqa,
     &       'local_w_density', l_jk, i_jk)
        status = status.and.MA_Push_get(MT_DBL,lsqa,
     &       'local_w_density', l_jl, i_jl)
        if (.not.status) 
     &       call errquit('force_dft: could not allocate scratch',1)
c
c     store total DM in ga_dens(1)
c
      if ( ipol .eq. 2 ) then
        call ga_dadd ( one,iga_dens(1),one,iga_dens(2),iga_dens(1) )
      endif 
          call dftg_hf_x(xfac(1),
     &         dbl_mb(i_ik), dbl_mb(i_jl),
     &         dbl_mb(i_il), dbl_mb(i_jk),
     &         force, iga_dens,
     &         max_at_bf, AO_bas_han, nproc,
     &         natom, dbl_mb(k_scr), lscratch, dbl_mb(k_buf),
     &         lbuf, tol2e,
     &         oactive, nactive )

        status = MA_Pop_Stack ( l_jl )
        status = MA_Pop_Stack ( l_jk )
        status = MA_Pop_Stack ( l_il )
        status = MA_Pop_Stack ( l_ik )
        status = MA_Pop_Stack ( l_scr )
        status = MA_Pop_Stack ( l_buf )
c
c     restore alpha DM in g_dens(1)
c
      if ( ipol .eq. 2 ) then
        call ga_dadd ( one,iga_dens(1),onem,iga_dens(2),iga_dens(1) )
      endif 

        if (.not.status) 
     &       call errquit(' dftg_getcx: failed to deallocate',11)

      endif
      call ga_sync
      tot=-xfac(1)
      do i=1,10
        tot=tot+xfac(i)+cfac(i)
      enddo

      if(abs(tot).gt.1e-8) then
c     
c       Compute via direct numerical quadrature.
c     
        call dftg_quadv0(rtdb,geom, AO_bas_han, 
     &       ipol,
     &       iga_dens, force,natom, oactive, nactive) 
      endif

      return
      end
