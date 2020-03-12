C>
C> \ingroup cosmo
C> @{
C>
C> \file cosmo_def_radii.F
C> Load default COSMO radii
C>
C> \brief Load the COSMO default atom radii
C>
C> Set the default atom radii values for the COSMO model. Currently
C> the radii are returned in the unit Angstrom.
C>
C> For the regular COSMO model and where available the radii from
C> Ref.[1] are used. Where Ref.[1] does not provide radii they are
C> derived by multiplying the van der Waals radii of Ref.[2] with
C> a factor 1.17 (for heavy elements the factor is 1.17*1.9).
C>
C> For the SMD model the radii were taken from Ref.[3] except for
C> Hydrogen. For Hydrogen the radius by Bondi [2] is used.
C> 
C> ### References ###
C>
C> [1] A. Klamt, V. Jonas, T. B&uuml;rger, J.C.W. Lohrenz,
C>     "Refinement and parametrization of COSMO-RS",
C>     <i>J. Phys. Chem. A</i> (1998) <b>102</b>, pp 5074-5085, DOI:
C>     <a href="https://doi.org/10.1021/jp980017s">
C>     10.1021/jp980017s</a>.
C>
C> [2] A. Bondi,
C>     "van der Waals volumes and radii",
C>     <i>J. Phys. Chem.</i> (1964) <b>68</b>, pp 441-451, DOI:
C>     <a href="https://doi.org/10.1021/j100785a001">
C>     10.1021/j100785a001</a>.
C>
C> [3] M. Mantina, R. Valero, C.J. Cramer, D.G. Truhlar,
C>     "Atomic radii of the elements", In CRC handbook of chemistry
C>     and physics, 91st Edition, 2010-2011; W.M. Haynes, Ed.;
C>     CRC Press, Boca Raton, FL, 2010; pp. 9-49 -- 9-50,
C>     ISBN-13: 978-1439820773.
C>
      subroutine cosmo_def_radii(rtdb,geom,nat,radius)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
c
      integer rtdb
      integer geom !< [Input] The geometry handle
      integer nat  !< [Input] The number of atoms
c
      double precision radius(nat) !< [Output] The atom radii
c
      logical ostatus
      integer mxelm
      parameter(mxelm = 102)
      double precision vander(mxelm),smdrad(mxelm)
      double precision rad
      integer l_coszan, k_coszan
      integer l_costag, k_costag
      integer l_coscoor, k_coscoor
      integer gn, iat
      integer h_crad, i_crad
      integer h_ctag, i_ctag
      integer mtype,melem,nmap
      character*26 mdate
      character*16 ctag(203), aname
      double precision crad(203)
      integer i
      character*2 symb
      character*16 elem
      integer atn
      logical isbq
c
c MN solvation models -->
C
      logical do_cosmo_smd, status
      double precision sola
c
c <-- MN solvation models
C
C     ----- van der waals radii in angstroms	-----
C
      double precision def, zero
      data def  /1.80D+00/
      data zero /0.00D+00/
c
c      data vdwr /
c     1   0.80,0.49,0.00,0.00,0.00,1.65,1.55,1.50,1.50,0.00,
c     2   2.30,1.70,2.05,2.10,1.85,1.80,1.80,0.00,2.80,2.75,
c     3   0.00,0.00,1.20,0.00,0.00,0.00,2.70,0.00,0.00,0.00,
c     4   0.00,0.00,0.00,1.90,1.90,0.00,0.00,0.00,0.00,1.55,
c     5   0.00,1.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     6   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     7   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     8   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     9   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.65,
c     2   0.00,0.00,0.00/
c
c    default radii provided by Andreas Klamt (Cosmologic)
c    vdw radii: 1.17 (+/- 0.02) * Bondi radius 
c    (Bondi, J. Phys. Chem., 68, 441, 1964)
c
c    optimal vdw radii for H, C, N, O, F, S, Cl, Br, I
c    (Klamt et al, J. Phys. Chem. A, 102, 5074 (1998)
c
c    for heavy elements: 1.17*1.9
c
      data (vander(iat),iat=1,mxelm)
     1 / 1.300,1.638,1.404,1.053,2.0475,2.00,
     2   1.830,1.720,1.720,1.8018,1.755,1.638,
     3   1.404,2.457,2.106,2.160,2.05,2.223,
     4   2.223,2.223,2.223,2.223,2.223,2.223,
     5   2.223,2.223,2.223,2.223,2.223,2.223,
     6   2.223,2.223,2.223,2.223,2.160,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223,
     8   2.223,2.223,2.223,2.223,2.223,2.223,
     9   2.223,2.223,2.223,2.223,2.320,2.223,
     1   2.223,2.223,2.223,2.223,2.223,2.223,
     2   2.223,2.223,2.223,2.223,2.223,2.223,
     3   2.223,2.223,2.223,2.223,2.223,2.223,
     4   2.223,2.223,2.223,2.223,2.223,2.223,
     5   2.223,2.223,2.223,2.223,2.223,2.223,
     6   2.223,2.223,2.223,2.223,2.223,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223/
c
c     vdw radii from Merz/Kollman/Singh
c
c      data (vander(i),i=1,36)
c     1                  /1.20D+00,1.20D+00,1.37D+00,1.45D+00,
c     2 1.45D+00,1.50D+00,1.50D+00,1.40D+00,1.35D+00,1.30D+00,
c     3 1.57D+00,1.36D+00,1.24D+00,1.17D+00,1.80D+00,1.75D+00,
c     4 1.70D+00,19*0.0D+00/
c
c     VDW RADII FROM GAVEZZOTTI (J.AM.CHEM.SOC. 105, 5220 (1983))
c     SUPPLEMENTED BY VALUES FOR BORON AND THIRD AND FOURTH ROWS
c     FROM M.A.SPACKMAN (J.CHEM.PHYS. 85, 6579 (1986))
c
c     DATA (VANDER(I),I=1,36)
c    1                  /1.20D+00,0.00D+00,0.00D+00,0.00D+00,
c    2 1.85D+00,1.50D+00,1.50D+00,1.40D+00,1.35D+00,0.00D+00,
c    3 0.00D+00,0.00D+00,2.07D+00,2.05D+00,1.96D+00,1.89D+00,
c    4 1.80D+00,19*0.0D+00/
c
c     VDW RADII FROM BRENEMAN & WIBERG
c
c     DATA (VANDER(I),I=1,36)
c    1                  /1.45D+00,1.45D+00,1.50D+00,1.50D+00,
c    2 1.50D+00,1.50D+00,1.70D+00,1.70D+00,1.70D+00,1.70D+00,
c    3 2.00D+00,2.00D+00,2.00D+00,2.00D+00,2.00D+00,2.00D+00,
c    4 2.00D+00,19*0.0D+00/
c
c MN solvation models -->
c
      data (smdrad(iat),iat=1,mxelm) /
     $ 1.20,         ! H  opt
     $ 1.40,         ! He
     $ 1.82,         ! Li
     $ 1.53,         ! Be
     $ 1.92,         ! B
     $ 1.85,         ! C  opt
     $ 1.89,         ! N  opt
     $ 1.52,         ! O  opt
     $ 1.73,         ! F  opt
     $ 1.54,         ! Ne
     $ 2.27,         ! Na
     $ 1.73,         ! Mg
     $ 1.84,         ! Al
     $ 2.47,         ! Si opt
     $ 2.12,         ! P  opt
     $ 2.49,         ! S  opt
     $ 2.38,         ! Cl opt
     $ 1.88,         ! Ar
     $ 2.75,         ! K
     $ 2.31,         ! Ca
     $ 2.16,         ! Sc
     $ 1.87,         ! Ti
     $ 1.79,         ! V
     $ 1.89,         ! Cr
     $ 1.97,         ! Mn
     $ 1.94,         ! Fe
     $ 1.92,         ! Co
     $ 1.84,         ! Ni
     $ 1.86,         ! Cu
     $ 2.10,         ! Zn
     $ 1.87,         ! Ga
     $ 2.11,         ! Ge
     $ 1.85,         ! As
     $ 1.90,         ! Se
     $ 3.06,         ! Br opt
     $ 2.02,         ! Kr
     $ 3.03,         ! Rb
     $ 2.49,         ! Sr
     $ 2.19,         ! Y
     $ 1.86,         ! Zr
     $ 2.07,         ! Nb
     $ 2.09,         ! Mo
     $ 2.09,         ! Tc
     $ 2.07,         ! Ru
     $ 1.95,         ! Rh
     $ 2.02,         ! Pd
     $ 2.03,         ! Ag
     $ 2.30,         ! Cd
     $ 1.93,         ! In
     $ 2.17,         ! Sn
     $ 2.06,         ! Sb
     $ 2.06,         ! Te
     $ 1.98,         ! I
     $ 2.16,         ! Xe
     $ 3.43,         ! Cs
     $ 2.68,         ! Ba
     $ 2.40,         ! La
     $ 2.35,         ! Ce
     $ 2.39,         ! Pr
     $ 2.29,         ! Nd
     $ 2.36,         ! Pm
     $ 2.29,         ! Sm
     $ 2.33,         ! Eu
     $ 2.37,         ! Gd
     $ 2.21,         ! Tb
     $ 2.29,         ! Dy
     $ 2.16,         ! Ho
     $ 2.35,         ! Er
     $ 2.27,         ! Tm
     $ 2.42,         ! Yb
     $ 2.21,         ! Lu
     $ 2.12,         ! Hf
     $ 2.17,         ! Ta
     $ 2.10,         ! W
     $ 2.17,         ! Re
     $ 2.16,         ! Os
     $ 2.02,         ! Ir
     $ 2.09,         ! Pt
     $ 2.17,         ! Au
     $ 2.09,         ! Hg
     $ 1.96,         ! Tl
     $ 2.02,         ! Pb
     $ 2.07,         ! Bi
     $ 1.97,         ! Po
     $ 2.02,         ! At
     $ 2.20,         ! Rn
     $ 3.48,         ! Fr
     $ 2.83,         ! Ra
     $ 2.60,         ! Ac
     $ 2.37,         ! Th
     $ 2.43,         ! Pa
     $ 2.40,         ! U
     $ 2.21,         ! Np
     $ 2.43,         ! Pu
     $ 2.44,         ! Am
     $ 2.45,         ! Cm
     $ 2.44,         ! Bk
     $ 2.45,         ! Cf
     $ 2.45,         ! Es
     $ 2.45,         ! Fm
     $ 2.46,         ! Md
     $ 2.46/         ! No
c
      if (.not. 
     $ rtdb_get(rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $ call errquit(
     $ 'cosmo_def_radii: cannot get do_cosmo_smd from rtdb',
     $ 0,rtdb_err)
      status = .false.
      if (do_cosmo_smd) 
     $ status = rtdb_get(rtdb,'cosmo:sola',mt_dbl,1,sola)
      if (.not.status) sola = -1.d0
c     use SMD radii if do_cosmo_smd = .true.
      if (do_cosmo_smd.and.sola.ge.0.D0) then      
       if (sola.lt.0.43D0) then
        smdrad(8)=1.52D0+1.80D0*(0.43D0-SOLA)
       else
        smdrad(8)=1.52D0
       endif
       do iat=1,mxelm
        vander(iat)=smdrad(iat)
       enddo 
      endif
c
c <-- MN solvation models
c
      if(rtdb_get_info(rtdb,'cosmo:crad',mtype,nmap,mdate)) then
        if(.not.rtdb_get(rtdb,'cosmo:crad',mt_dbl,nmap,crad))
     $   call errquit('cosmo_def_radii: rtdb get crad',917,rtdb_err)
        if(.not.rtdb_cget(rtdb,'cosmo:ctag',nmap,ctag))
     $   call errquit('cosmo_def_radii: rtdb get ctag',917,rtdb_err)
      else
        nmap = 0
      end if

      if(.not.ma_push_get(mt_dbl,nat*3,'coord',l_coscoor,k_coscoor))
     & call errquit('cosmo_def_radii malloc k_coscoor failed',
     &              911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nat,'cosmo z',l_coszan,k_coszan)) call
     &  errquit('cosmo_def_radii malloc k_coszan failed',911,MA_ERR)
      if(.not.ma_push_get(mt_byte,nat*16,'tags',l_costag,k_costag)) call
     &  errquit('cosmo_def_radii malloc k_costag failed',911,MA_ERR)
c
      if(.not.geom_cart_get(geom,nat,byte_mb(k_costag),
     &   dbl_mb(k_coscoor),dbl_mb(k_coszan))) call errquit
     $    (' cosmo_def_radii: geom_cart_get failed.',911, GEOM_ERR)
c
      atomloop: do iat=1,nat
         call util_set_ma_char_arrays(16,
     $                                byte_mb(k_costag+16*(iat-1)),
     $                                aname)
         do i=1,nmap
            if(inp_compare(.false.,
     &                     aname,
     &                     ctag(i))) then
              radius(iat) = crad(i)
              if(ga_nodeid().eq.0) 
     &            write(*,9977) ctag(i),crad(i)
              cycle atomloop
             
            endif
         enddo
         if(aname(1:1).eq.'X'.or.aname(1:1).eq.'x') then
c     minuscule radius for X                   
            radius(iat)=1d-9
            cycle atomloop
         endif
         if(.not.geom_tag_to_element(aname, symb, elem, atn)) then
            if (symb .ne. 'bq') then
               call errquit
     &              ('cosmodefradii: center is neither atom nor bq'
     +           //aname, 0,0)
            else
#if 0               
               isbq=.true.
               if(abs(atn).gt.1d-6) then
c     point-charge with non zero charge
                  call errquit
     &                 ('cosmodefradii: not ready for charged bqs'
     +                 //aname, 0,0)
               else
c bq with zero charge. extract element
                  if (.not. geom_tag_to_element(aname(3:), symb, 
     &                 elem, atn)) 
     C                 call errquit
     &                 ('cosmodefradii: failed extract bq elem name'
     +                 //aname, 0,0)
               endif
#else
            if(ga_nodeid().eq.0) then
               write(6,*) 
               write(6,*) ' Solvation models not compatible with bqs '
               write(6,*) ' Please remove bqs from geometry '
               write(6,*) 
            endif
            call errquit
     &           ('cosmodefradii: not ready for '
     +           //aname, 0,0)
#endif
            endif
         endif
         do i=1,nmap
            if(inp_compare(.false.,
     &                     symb,
     &                     ctag(i))) then
              radius(iat) = crad(i)
              if(ga_nodeid().eq.0) 
     &            write(*,9978) ctag(i),symb, crad(i)
              cycle atomloop
             
            endif
         enddo
         rad=vander(nint(dbl_mb(k_coszan+iat-1)))
         if (geom_any_ecp(geom).and.geom_ecp_get(geom,iat)) then
            if (.not.geom_tag_to_element(byte_mb(k_costag+16*(iat-1)), 
     &         symb, elem, gn))
     &          call errquit('cosmo_def_radii:geom_tag_to_element',0,0)
            rad=vander(gn)
         endif
         if(rad.eq.zero) then
            rad=def
         endif
         radius(iat)=rad
      enddo atomloop

      if(.not.ma_pop_stack(l_costag)) call
     &  errquit('cosmo_def_radii dealloc k_costag failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_coszan)) call
     &  errquit('cosmo_def_radii dealloc k_coszan failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_coscoor)) call
     &  errquit('cosmo_def_radii dealloc k_coscoor failed',911,MA_ERR)
9977  format(' setting custom COSMO radius based on name match: ', 
     &        a6,f7.3)
9978  format(' setting custom COSMO radius based on element match: ',
     &        a6,a4,f7.3)
c
      end
C>
C> @}
c $Id$
