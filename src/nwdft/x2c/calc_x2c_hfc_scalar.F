      subroutine calc_x2c_hfc_scalar (rtdb, needu2c, nexc,
     &  g_dscf, g_dens_at, g_u, g_u2c, g_c2u)

c     -----------------------------------------------------------------
c     purpose: calculate one-electron property using the 
c     two-component exact-decoupling Hamiltonian (scalar version)
c     here: hyperfine coupling

c     Integrals are supposed to be initialized by the calling
c     routine.

c     X2C Driver (c) 2012,2022 J. Autschbach, SUNY, jochena@buffalo.edu
c     HFC code (c) 2016,2022 J. Autschbach, SUNY, jochena@buffalo.edu

c     The code was adapted to the github nwchem version in 2022.
c     See further comments in the header of file
c     calc_x2c_1e_scalar_ga.F

c     Approach [1]: The transformation matrices g_u come from
c     calc_x2c_1e_scalar_ga if certain options are set. The integrals
c     for the nuclear magnetic hyperfine interaction come from the Rys
c     quadrature code in the NWints/hondo collection, adapted for this
c     routine. See also src/NWints/hondo/hnd_ops.F
c
c     [1] J. Autschbach, 'Relativistic effects on electron-nucleus
c     hyperfine coupling studied with an "exact 2-component"
c     (X2C) method', J. Chem. Theory Comput. 2017, 13, 710–718.
c     https://doi.org/10.1021/acs.jctc.6b01014
c
c     Because of yet unspecified issues in the branch used for the
c     calculations in [1], which was based on a 2012 developer's version
c     of nwchem where the x2c one-electron Hamiltonian code was added
c     first, and the hyperfine code was added in 2015/16, the point
c     nucleus calculations with high-exponent basis functions apparently
c     did not converge as well as they should have. The x2c port in the
c     2022 version of nwchem yields slightly lower energies and slightly
c     better hyperfine integrals. As a result, the QA jobs added during
c     the development of the 2022 port of the code to the github version
c     of nwchem give slightly different, and almost certainly more
c     accurate, results compared to those reported in [1].
c     
c     -----------------------------------------------------------------

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "dftpara.fh"
#include "cdft.fh"
#include "geom.fh"
#include "rel_consts.fh"
#include "cosmo.fh"
#include "cint1cache.fh"

c     subroutine arguments:

      integer rtdb              ! [in] runtime database GA handle
      integer nexc              ! [in] no. of XC terms
      integer g_dscf(2)         ! [in] density matrix from SCF
      integer g_dens_at(2)      ! [in] atom densities for model pot.
      integer g_u(2)            ! [in] decoupling matrices
      integer g_u2c             ! [in] basis transformation (optional)
      integer g_c2u             ! [in] basis transformation (optional)
      logical needu2c           ! [in] option for u2c

c     local variables:

      integer ubas, cbas  
      integer nubas, nubas2, ncbas
      integer me, type
      double precision sum
      character*2  symbol
      character*16 element, at_tag

c     GA handles used for temp storage in this routine:
      integer g_P, g_Ptmp, g_magint, g_temp, g_s

c     dimensioning arrays
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)

c$$$  c     MA handles for arrays passed to relop routine:
c$$$  integer l_dmem , k_dmem
c$$$  integer l_amP  , k_amP
c$$$  integer l_omP  , k_omP 
c$$$  integer l_amUL , k_amUL
c$$$  integer l_amUS , k_amUS

c     other MA handles
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt

c     bq variables (MV)
      logical dobq
      integer bq_ncent
      integer i_cbq
      integer i_qbq
      double precision elpotbq

      double precision ga_trace_diag
      external ga_trace_diag
      double precision tol2e, accval 
      logical master, analytic_pvp, status
      integer dbg

      logical atom_gfac
      external atom_gfac
      
      integer nat, iat, kat, ntemp, idir, jdir, kk
      
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision rr, rr5, rtemp
      double precision hfc(3,3), amat(3,3)
      double precision origin(3)

c     unit conversions:
      double precision gnuc, hbar, ge, emf, vl, auev, aiso,
     &  evmhz, gmhz, fac, betae, betan, convf, con, gnu

      integer atn, isonr

c     these are for converting the results and the same
c     as used in hnd_spind.F, so that we can compare

      data gnuc   /5.05078343d-27/ ! Nuclear magneton
      data hbar   /1.05457168d-34/ ! Planck constant over 2 pi
      data ge     /2.002319304386d+00/ ! Electron g-factor
      data emf    /1836.152701d+00/ ! Proton-electron mass ratio
      data vl     /137.0359895d+00/ ! Speed of light in au
      data auev   /27.2113961d+00/ ! Conversion from au to eV
      data evmhz  /2.41798836d+08/ ! Conversion from eV to MHz
      data gmhz   /2.8025d+00/  ! Conversion from Gauss to MHz

      double precision pi, deg

      character*14 pname, ctmp

c     =================================================================

c     ---------------
c     initializations
c     ---------------

c     options / settings

      dbg=1 ! debug level. >1 gives huge matrix output
      me=ga_nodeid()
      master = me.eq.0
      pname = "x2c_hfc_scalar"

      if (dbg>0.and.master) write(luout,*) 'entering '//pname

c     do nothing if this is not a spin-unrestricted calculation
      if (ipol.lt.2) then
        if (master) write (luout,*)
     &    pname//': not a spin-unrestricted calculation. aborting'
        return
      end if

      call ga_sync()

c     constants:

      tol2e=10.d0**(-itol2e)

      call int_acc_get(accval)

      if (dbg>0.and.master) write(luout,*) 'accval = ',accval

c     conversion factors and such. We use the very same
c     here as in routine hnd_spind, so that we can compare
c     the results directly

      pi   = acos(-1.0d0)
      fac  =(4.0d0*pi/3.0d0)
      betae=1.0d0/(2.0d0*vl)
      betan=betae/emf
      convf=auev*evmhz
      con  =ge*betae*betan*convf

c     cau = speed of light, from rel_consts.fh
      if (master.and.dbg>0) write(luout,*) 'cau = ',cau

c     the off-diagonal block of the x2c transformation matrix 
c     carries a factor of 1/(2c). The 4-component magnetic operator
c     carries a factor of c but is not included in the AO matrix
c     elements. Therefore, we scale g_u(2) by c to make up for
c     it and un-do this at the very end

      call ga_scale(g_u(2),cau)

      if (dbg>1) then
        call ga_print(g_u(1))
        call ga_print(g_u(2))
      end if      

c     process options from rtdb set via input

      if(.not.rtdb_get(rtdb,'prop:hfc:atom',MT_INT,1,kat))
     &  kat = 1


c     the ncbas, nubas gymnastics is usually needed because
c     the x2c operator is always calculated in an uncontracted
c     basis

      ncbas = 0
      nubas = 0

      if (needu2c) then
        ubas = x2c_bas
        if (.not. bas_numbf(ubas,nubas))
     &    call errquit(pname//': nbas x2c error',nubas, BASIS_ERR)
        cbas = ao_bas_han
        if (.not. bas_numbf(cbas,ncbas))
     &    call errquit(pname//': nbas ao error',nubas, BASIS_ERR)
        if (master) write(luout,*)
     &    pname//': working with uncontracted auxiliary basis'
c       sanity check: do we have the correct u2c array?
        call ga_inquire_name(g_u2c, ctmp)
        if (ctmp(1:3) .ne. 'u2c') call errquit(
     &    pname//': g_u2c inconsistency', UNKNOWN_ERR)
      else
        ubas = ao_bas_han
        if (.not. bas_numbf(ubas,nubas))
     &    call errquit(pname//': nbas ao error',nubas, BASIS_ERR)
        cbas = ubas
        ncbas = nubas
        if (master) write(luout,*)
     &    pname//': the basis seems to be fully uncontracted ...'
      end if

c     make sure the input GAs for the transformation have the correct
c     dimension

      call ga_inquire(g_u(1),type,ntemp,ntemp)
      if (nubas.ne.ntemp) then
        call errquit(pname//': inconsistent dim. for g_u',nubas,GA_ERR)
      end if

      if (master.and.dbg>0) write (luout,*) pname//': nubas =',nubas
      nubas2 = nubas * nubas

c     -------------
c     Calculate HFC
c     -------------

c     initialize integrals that were not already initialized

      call hnd_giao_init(ubas,1)

c     ----- define points for calculation -----

      status=geom_ncent(geom,nat)
      if (master.and.dbg>0) write (luout,*) 'x2c: nat =',nat
c     
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &  call errquit(pname//': ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &  call errquit(pname//': ma failed',911,MA_ERR)
c     
      do  iat=1,nat
        status=geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &    dbl_mb(k_zanpt+iat-1))
      end do

      alo(1) = nubas
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nubas
      ahi(2) = nubas   
      ahi(3) = 9

      if (.not.nga_create(MT_DBL,3,ahi,'g_mag',alo,g_magint)) call
     &  errquit(pname//': nga_create failed g_magint',11,GA_ERR)
      call ga_zero(g_magint)

      if (dbg>1) then
        if (.not.nga_create(MT_DBL,3,ahi,'temp',alo,g_temp)) call
     &     errquit(pname//': nga_create failed g_temp',14,GA_ERR)
        call ga_zero(g_temp)
      end if

c     calculate magnetic operator integrals in the mixed
c     upper- / lower-component basis

      origin(1) = dbl_mb(k_xyzpt  +3*(kat-1))
      origin(2) = dbl_mb(k_xyzpt+1+3*(kat-1))
      origin(3) = dbl_mb(k_xyzpt+2+3*(kat-1))

      if (dbg>0.and.master) write (luout,*) pname//': origin=',origin(:)

      call int_giao_1ega(ubas,ubas,g_magint,'ops',
     &  origin,1, oskel)
      
c     debug integral matrices:
      
      if (dbg>1) then
        alo(1) = 1
        ahi(1) = nubas
        alo(2) = 1
        ahi(1) = nubas
        kk = 0
        do jdir = 1,3
          do idir = 1,3
            kk = kk+1
            alo(3) = kk
            ahi(3) = kk
            if (master.and.dbg>1) write(luout,*)
     &         'magnetic integrals block ',idir,jdir
            call nga_print_patch(g_magint,alo,ahi,1)
            call ga_zero(g_temp)
          end do
        end do
      end if

      if (master.and.dbg>0) write(luout,*) pname//': AO integrals done'

c     for sanity check, contract the spin density with the
c     overlap matrix! This is done in the contracted basis
      
      alo(1) = ncbas
      alo(2) = -1
      alo(3) = -1
      ahi(1) = ncbas
      ahi(2) = ncbas
      ahi(3) = 1
      
      if (.not.nga_create(MT_DBL,2,ahi,'smat',alo,g_s)) call 
     &  errquit(pname//': nga_create failed g_s',11,GA_ERR)
      call ga_zero(g_s)
      if (ocache) call int_1e_uncache_ga()
      call int_1e_ga(cbas, cbas ,g_s, 'overlap'   ,oskel)
      call ga_symmetrize(g_s)
      if (.not.nga_create(MT_DBL,2,ahi,'temp',alo,g_temp)) call
     &  errquit(pname//': nga_create failed g_temp',14,GA_ERR)
      call ga_zero(g_temp)

c     dimensioning for matrix multiplications used below
      
      alo(1) = 1                ! range of A, first dimension 
      ahi(1) = ncbas
      alo(2) = 1                ! range of A, second dimension
      ahi(2) = ncbas
      blo(1) = 1                ! range of B, first dimension, etc.
      bhi(1) = ncbas
      blo(2) = 1
      bhi(2) = ncbas
      clo(1) = 1
      chi(1) = ncbas
      clo(2) = 1
      chi(2) = ncbas

      call ga_zero(g_temp)
      call nga_matmul_patch('n','n',1d0,0d0,
     &  g_dscf(1),alo,ahi,
     &  g_s,blo,bhi,
     &  g_temp,clo,chi)
      sum = ga_trace_diag(g_temp)
      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': tr[P(alfa) * S] : ',sum

      call ga_zero(g_temp)
      call nga_matmul_patch('n','n',1d0,0d0,
     &  g_dscf(2),alo,ahi,
     &  g_s,blo,bhi,
     &  g_temp,clo,chi)
      sum = ga_trace_diag(g_temp)
      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': tr[P(beta) * S] : ',sum

c     subtract alpha and beta density,
c     to get the spin density. Note, g_dscf(1) contains the SUM

      call ga_dadd(1.d0, g_dscf(1), -2.d0, g_dscf(2), g_dscf(1))

      call ga_zero(g_temp)
      call nga_matmul_patch('n','n',1d0,0d0,
     &  g_dscf(1),alo,ahi,
     &  g_s,blo,bhi,
     &  g_temp,clo,chi)
      sum = ga_trace_diag(g_temp)
      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': tr[P(spin) * S] : ',sum

      if (master .and. dbg>1) then
        write (luout,*)
     &    pname//': Spin density matrix:'
      end if
      if (dbg>1) call ga_print(g_dscf(1))

      if (.not. ga_destroy(g_temp)) call errquit
     &  (pname//': temp ga corrupt?',14, GA_ERR)
      if (.not. ga_destroy(g_s)) call errquit
     &  (pname//': temp ga corrupt?',11, GA_ERR)


      if (dbg>0 .and.master) write (luout,*) pname//' before hfc step'

c     --------------------------
c     calculate HFC for atom kat
c     --------------------------
      
      if (dbg>0 .and.master) write (luout,*) pname//': kat=',kat
      
      if (.not. geom_cent_tag(geom,kat,at_tag)) call
     &  errquit(pname//': geom_cent_tag failed',0,GEOM_ERR)
      if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) call
     &  errquit(pname//': geom_tag_to_element failed',0,GEOM_ERR)
      if (.not. atom_gfac(dbl_mb(k_zanpt+kat-1),gnu,isonr)) call
     &  errquit(pname//': atom_gfac failed',0, UERR)
      xp  = dbl_mb(k_xyzpt  +3*(kat-1))
      yp  = dbl_mb(k_xyzpt+1+3*(kat-1))
      zp  = dbl_mb(k_xyzpt+2+3*(kat-1))      

c     aiso = 2.0d0*fac*con*gnu 
c     fac is 4 pi / 3, we have that already in the AO integrals 
      aiso = 2.0d0*con*gnu 
      if (dbg>0.and.master) write(luout,*) pname//': con,gnu: ',con,gnu
                
c     ---------------------------------------------------
c     calculate picture-change corrected HFC operator
c     for the selected reference point (nuclear position)
c     ---------------------------------------------------
      
c     ------------------------------------------------------------
c     in a loop over the HFC components, transform the 
c     operator to two-component form, and calculate the electronic
c     HFC component again
c     ------------------------------------------------------------
      
      kk = 0 ! index for field-spin in integral matrix array
      do jdir = 1,3
        do idir = 1,3
          kk = kk+1

          if(.not.ga_create(mt_dbl, nubas, nubas,
     &      'g_Ptmp',1,1,g_Ptmp)) call errquit(
     &      pname//': error creating g_Ptmp',1, GA_ERR)
          
          if(.not.ga_create(mt_dbl, nubas, nubas,
     &      'g_P',1,1,g_P)) call errquit(
     &      pname//': error creating g_P',1, GA_ERR)
          
          call ga_zero(g_P)     ! property integral arrray
          call ga_zero(g_Ptmp)

          alo(1) = 1          
          ahi(1) = nubas
          alo(2) = 1           
          ahi(2) = nubas
          
          blo(1) = 1
          bhi(1) = nubas
          blo(2) = 1           
          bhi(2) = nubas
          blo(3) = 1
          chi(3) = 1
          
          clo(1) = 1
          chi(1) = nubas
          clo(2) = 1
          chi(2) = nubas
          clo(3) = 1
          chi(3) = 1
          
          alo(3) = kk
          ahi(3) = kk

c         x2c-transform the operator 

          if(.not.ga_create(mt_dbl, nubas, nubas,
     &      'x2ctmp',1,1,g_temp)) call errquit(
     &      pname//': error creating g_temp',1, GA_ERR)
          
          call ga_zero(g_P)

          call nga_copy_patch('n',g_magint,alo,ahi,g_Ptmp,blo,bhi)
          call ga_dgemm('n','n',nubas,nubas,nubas,1.0d0,
     &      g_Ptmp,g_u(2), 0.0d0,g_temp)
          call ga_dgemm('t','n',nubas,nubas,nubas,1.0d0,
     &      g_u(1),g_temp, 0.0d0,g_P)

          call nga_copy_patch('t',g_magint,alo,ahi,g_Ptmp,blo,bhi)
          call ga_dgemm('n','n',nubas,nubas,nubas,1.0d0,
     &      g_Ptmp,g_u(1), 0.0d0,g_temp)
          call ga_dgemm('t','n',nubas,nubas,nubas,1.0d0,
     &      g_u(2),g_temp, 1.0d0,g_P)
          
          if (.not. ga_destroy(g_temp)) call errquit(
     &      pname//': error destroying g_temp',1, GA_ERR)
          
c         might need to transform the operator matrix
c         to the contracted basis

          if (needu2c) then
c           save g_P in g_Ptmp and re-dimension g_P
            call ga_copy(g_P, g_Ptmp)
            if (.not. ga_destroy(g_P)) call errquit(
     &        pname//': error destroying g_P',1, GA_ERR)
            if(.not.ga_create(mt_dbl, ncbas, ncbas,
     &        'g_P',1,1,g_P)) call errquit(
     &        pname//': error re-creating g_P',1, GA_ERR)
            
            if(.not.ga_create(mt_dbl, nubas, ncbas,
     &        'g_temp',1,1,g_temp)) call errquit(
     &        pname//': error creating g_temp',1, GA_ERR)
            call ga_dgemm('n','n',nubas,ncbas,nubas,1.0d0,
     &        g_Ptmp,g_u2c, 0.0d0,g_temp)
            call ga_dgemm('t','n',ncbas,ncbas,nubas,1.0d0,
     &        g_u2c,g_temp, 0.0d0,g_P)
            if (.not. ga_destroy(g_temp)) call errquit(
     &        pname//': error destroying g_temp',1, GA_ERR)
            
          end if                ! needu2c

          if (master .and. dbg>1) then
            write (luout,*)
     &        'MAG after transformation, kat, idir, jdir',kat,idir,jdir
          end if
          if (dbg>1) call ga_print(g_P)          

c         -----------------------------------------
c         calculate electronic HFC again, this time 
c         using the transformed operator
c         -----------------------------------------

          if(.not.ga_create(mt_dbl, ncbas, ncbas,
     &      'g_temp',1,1,g_temp)) call errquit(
     &      pname//': error creating g_temp',1, GA_ERR)
          
c         tr[P*HFC]
          alo(1) = 1           
          ahi(1) = ncbas
          alo(2) = 1
          ahi(2) = ncbas
          alo(3) = 1
          ahi(3) = ncbas
          blo(1) = 1           
          bhi(1) = ncbas
          blo(2) = 1
          bhi(2) = ncbas
          blo(3) = 1
          bhi(3) = ncbas
          clo(1) = 1           
          chi(1) = ncbas
          clo(2) = 1
          chi(2) = ncbas
          clo(3) = 1
          chi(3) = ncbas

          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &      g_dscf(1),alo,ahi,
     &      g_P, blo,bhi,
     &      g_temp,clo,chi)
          sum = ga_trace_diag(g_temp)

          hfc(idir,jdir) = sum

          if (.not. ga_destroy(g_temp)) call errquit(
     &      pname//': error destroying g_temp',1, GA_ERR)
          if (.not. ga_destroy(g_P)) call errquit(
     &      pname//': error destroying g_P',1, GA_ERR)
          if (.not. ga_destroy(g_Ptmp)) call errquit(
     &      pname//': error destroying g_P',1, GA_ERR)
          
        end do                  ! jdir
      end do                    ! idir

c     spin-field reduction:

      amat(1,1) = hfc(2,2) + hfc(3,3)
      amat(2,2) = hfc(1,1) + hfc(3,3)
      amat(3,3) = hfc(1,1) + hfc(2,2)
      do idir=1,3
        do jdir = 1,3
          if (idir.ne.jdir) amat(idir,jdir) = -hfc(jdir,idir)
        end do
      end do
      
c     print results:
      
      if (master) then
        if (dbg>1) then
          write(luout,*) 'hfc:'
          call output(aiso*hfc(:,:),1,3,1,3,3,3,1)
        end if
        write (luout,*)
        write (luout,*) 'X2C HFC (au and MHz):'
        write(luout,9995) kat,isonr,symbol,xp,yp,zp
        write (luout,'(3(1x,3e20.10/))') -amat(:,:)
        write (luout,'(3(1x,3e20.10/))') -aiso*amat(:,:)

        rtemp = 0.0d0 
        do idir = 1,3
          rtemp = rtemp + amat(idir,idir)
        end do
        write(luout,9993)
     &    'X2C HFC isotropic =',(-1.0d0/3.0d0) * rtemp, 'au'
        write(luout,9993)
     &    '                  =',(-aiso/3.0d0) * rtemp, 'MHz'
        write (luout,*)
      end if

c       -------------------------
c       done with HFC calculation
c       -------------------------

      if (.not. ga_destroy(g_magint)) call errquit(
     &  pname//': ga corrupt?',11, GA_ERR)

c     ------- Deallocate MA memory ------
c     
      if (.not.ma_pop_stack(l_zanpt)) call errquit(
     &  pname//': ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit(
     &  pname//': ma_pop_stack of l_xyzpt failed',911,MA_ERR)

c     undo subtraction of alpha and beta density matrix:

      call ga_dadd(1.d0, g_dscf(1), 2.d0, g_dscf(2), g_dscf(1))

c     undo scaling of off-diagonal x2c transformation matrix

      call ga_scale(g_u(2),(1.0d0/cau))

c     --------
c     all done
c     --------

      call ga_sync()

      if (dbg>0.and.master) write(luout,*) 'leaving '//pname

c     =================================================================

 9994 format(/,' Total Spin Density (Fermi Contact Term)',/,
     1  ' ------------------ --------------------')
 9997 format(5x,'Atom ',11x,'x',9x,'y',9x,'z',7x,
     1  'Density (a.u.)',6x,'Aiso(MHz)',5x,'Aiso(Gauss)')
 9995 format(1x,i5,1x,i3,'-',a2,2x,3f10.5,2x,f13.6,5x,f13.6,2x,f13.6)
 9984 format(/,3x,'Spin-Dipolar term (a.u.)',/,8x,'xx',13x,'yy',
     1  13x,'zz',13x,'xy',13x,'xz',13x,'yz',/,1x,90(1h-))
 9983 format(6f15.6)
 9992 format(1x,60(1h-),/,5x,'Atom ',11x,'x',9x,'y',9x,'z',/,
     1  1x,60(1h-))
 9991 format(1x,i5,1x,i3,'-',a2,2x,3f10.5)
 9996 format(/,/,10x,' -------------------------------------------',
     1  /,10x,' Anisotropic Interaction (Spin-Dipolar Term)',
     2  /,10x,' -------------------------------------------',/)
 9993 format (1x,a,1x,f20.2,1x,a)


      return
      end

