      program testbasis
c $Id: testbasis.F,v 1.17 1999-07-08 19:10:51 d3e129 Exp $
      implicit none
c
#include "bas.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "inp.fh"
#include "stdio.fh"
c
c::friends functions
      logical bas_ucontinfo
      logical bas_getu_exponent
      logical bas_getu_coeff
      logical bas_setu_exponent
      logical bas_setu_coeff
      external bas_ucontinfo
      external bas_getu_exponent
      external bas_getu_coeff
      external bas_setu_exponent
      external bas_setu_coeff
c
      integer rtdb, geom, basis, ngen,nprim, iang
*      integer ncenters
      integer sphcart, i,j
      integer nbf,ibf,iat,icont, type
      character*256 namename, nametrans
*      character*16 drivtags(20)
*      double precision coords(3,20), charge(20)
      integer lexcf
      parameter (lexcf=400)
      double precision exp(lexcf), coeff(lexcf)
      logical status
      double precision expnt_new(3), coeff_new(4,3)
      logical bas_add_ucnt, bas_geobas_build
      external bas_add_ucnt, bas_geobas_build
      data expnt_new/1., 2., 3./
      data coeff_new/-1., -2., -3., 0.0, -4., -5., -6., 0.0, -7., -8.,
     $     -9., 0.0/
c
      if (.not. ma_init(MT_DBL, 400 000, 400 000))
     &      call errquit('testbasis: ma_init failed',911)

      status = rtdb_open('testbasis.rtdb','empty', rtdb)
      if (.not.status) call errquit
     &      ('testbasis rtdb open failed',911)
c
c      write(luout,*)' rtdb handle ', rtdb
*old      status = bas_321g_load(rtdb)
c
      open(unit=luin,file='testbasis.nw',
     &      form='formatted',
     &      access='sequential',
     &      status='old',
     &      err=99565)
      call input_parse(rtdb)
c
      if(.not.geom_create(geom,'geometry')) then
        write(6,*)' error getting geometry handle '
        stop ' error '
      endif
c
      if(.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &      call errquit('error loading geometry',911)
      status=geom_print(geom)
c
      basis = 0
      if(.not.bas_create(basis,'xy basis')) then
        write(6,*)' error getting basis handle '
        stop ' error '
      endif
      basis = 0
      if(.not.bas_create(basis,'wz basis')) then
        write(6,*)' error getting second basis handle '
        stop ' error '
      endif
      basis = 0
      if(.not.bas_create(basis,'ao basis')) then
        write(6,*)' error getting third basis handle '
        stop ' error '
      endif
c
c
      status = bas_rtdb_load(rtdb,geom,basis,'ao basis')
      status = bas_print(basis)
      status = gbs_map_print(basis)
c
      status = bas_continfo(basis,1,iang,nprim,ngen,sphcart)
      write(6,*)' user:query: status  cont 1 ',status
      write(6,*)' user:query: type    cont 1 ',iang
      write(6,*)' user:query: nprim   cont 1 ',nprim
      write(6,*)' user:query: ngen    cont 1 ',ngen
      write(6,*)' user:query: sphcart cont 1 ',sphcart
      if (nprim*ngen.gt.lexcf) call errquit
     &      (' lexcf too small ',(nprim*ngen))
      call dfill(lexcf,0.0d00,exp  ,1)
      call dfill(lexcf,0.0d00,coeff,1)
      status = bas_get_exponent(basis,1,exp)
      status = bas_get_coeff(basis,1,coeff)
      write(6,*)' user exponenents and coefficients '
      do 00100 i=1,nprim
        write(6,*)exp(i),(coeff(i+(j-1)*nprim),j=1,ngen)
00100 continue
c
      status = bas_ucontinfo(basis,1,type,nprim,ngen,sphcart)
      write(6,*)' unique:query: status  cont 1 ',status
      write(6,*)' unique:query: type    cont 1 ',type
      write(6,*)' unique:query: nprim   cont 1 ',nprim
      if (nprim*ngen.gt.lexcf) call errquit
     &      (' lexcf too small ',(nprim*ngen))
      write(6,*)' unique:query: ngen    cont 1 ',ngen
      write(6,*)' unique:query: sphcart cont 1 ',sphcart
      call dfill(lexcf,0.0d00,exp  ,1)
      call dfill(lexcf,0.0d00,coeff,1)
      status = bas_getu_exponent(basis,1,exp)
      status = bas_getu_coeff(basis,1,coeff)
      write(6,*)' unique exponenents and coefficients '
      do 00200 i=1,nprim
        write(6,*)exp(i),(coeff(i+(j-1)*nprim),j=1,ngen)
00200 continue
c
      exp(1) = 565.6589
      coeff(1) = 6.021023
      status = bas_setu_exponent(basis,1,exp,(nprim+1))
      status = bas_setu_coeff(basis,1,coeff,(nprim*ngen)+1)
      status = bas_setu_exponent(basis,1,exp,nprim)
      status = bas_setu_coeff(basis,1,coeff,nprim*ngen)
c
c
      status = bas_ucontinfo(basis,1,type,nprim,ngen,sphcart)
      write(6,*)' modified '
      write(6,*)' unique:query: status  cont 1 ',status
      write(6,*)' unique:query: type    cont 1 ',type
      write(6,*)' unique:query: nprim   cont 1 ',nprim
      write(6,*)' unique:query: ngen    cont 1 ',ngen
      write(6,*)' unique:query: sphcart cont 1 ',sphcart
      if (nprim*ngen.gt.lexcf) call errquit
     &      (' lexcf too small ',(nprim*ngen))
      status = bas_getu_exponent(basis,1,exp)
      status = bas_getu_coeff(basis,1,coeff)
      write(6,*)' exponenents and coefficients '
      do 00300 i=1,nprim
        write(6,*)exp(i),(coeff(i+(j-1)*nprim),j=1,ngen)
00300 continue
c
c
c     Try adding new contractions on an existing center
c
      write(6,*) ' adding 3*3 d function on H'
      if (.not. bas_add_ucnt(basis, 'H', 2, 3, 3, expnt_new,
     &    expnt_new, coeff_new, 4, 'none'))
     &    write(6,*) ' basis_add_ucnt failed'
      if (.not. bas_print(basis)) write(6,*) ' print ?'
c
c     Try adding new contractions on a new center
c
      write(6,*) ' adding 2*3 g function on Cl'
      if (.not. bas_add_ucnt(basis, 'Cl', 4, 2, 3, expnt_new,
     $     coeff_new, 4, 'none'))
     &    write(6,*) ' basis_add_ucnt failed'
      if (.not. bas_print(basis)) write(6,*) ' print ?'
c
c
      write(6,'(///,1x,a)')' bas_print_all '
      status = bas_print_all()
      call bas_err_info(' testbasis bas_err_info info ')
c
      status = bas_high_angular(basis,iang)
      write(6,*)' high angular momentum ', iang
      status = bas_version()
      namename = ' '
      nametrans = ' '
      status = bas_name(basis,namename, nametrans)
      write(6,*)' handle         :  ',basis, status
      i = inp_strlen(namename)
      write(6,'(a,a,a)')' name           : <',
     &       namename(1:i),'>'
      i = inp_strlen(nametrans)
      write(6,*)' translated name: <',
     &       nametrans(1:i),'>'
c
      write(6,*)' basis function to ce/cn check '
      status = bas_numbf(basis,nbf)
      do 00400 ibf=1,nbf
        status = bas_bf2ce(basis,ibf,iat)
        status = bas_bf2cn(basis,ibf,icont)
        write(6,*)' basis function ',ibf,
     &            ' is on center ',iat,
     &         ' and in contraction',icont
00400 continue
c
      status = bas_geobas_build(basis)
      write(6,*)' basis handle ',basis
      call bas_geomap_check(basis)
c
      status = bas_nbf_cn_max(basis,nbf)
      write(6,*)' max block of nbf on a shell for basis is:',nbf
c
      status = bas_nbf_ce_max(basis,nbf)
      write(6,*)' max block of nbf on an atom for basis is:',nbf
c
      status = bas_ncoef_cn_max(basis,nbf)
      write(6,*)' max num coefs in any contraction ',nbf
c
      status = bas_nprim_cn_max(basis,nbf)
      write(6,*)' max num prims in any contraction ',nbf
c
c
      stop ' testbasis done '
99565 continue
      stop ' error opening testbasis.nw'
      end
      subroutine bas_geomap_check(basisin)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
#include "geom.fh"
c
      integer basisin
c
      logical status
      integer basis, geom
      integer ncont, icont, ucont, ucont2
      integer cent_api, cent_sf, cent_mb
      integer inat, nat, atom_sf, atom_mb
      integer bflo_api, bflo_sf, bflo_mb
      integer bfhi_api, bfhi_sf, bfhi_mb
      integer cnlo_api, cnlo_sf, cnlo_mb
      integer cnhi_api, cnhi_sf, cnhi_mb
c
#include "bas_ibs_sfn.fh"
c
      status = bas_geom(basisin,geom)
      basis = basisin + BASIS_HANDLE_OFFSET
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '
      write(6,*)' basis geometry mapping checks '
      write(6,*)' basis handle is    :',basisin
      write(6,*)' geometry handle is :',geom
      write(6,*)' '
      write(6,*)' current memory heap pointers'
      write(6,*)' ibs_cn2ucn ma handle :',ibs_cn2ucn(H_ibs,basis)
      write(6,*)' ibs_cn2ucn ma pointer:',ibs_cn2ucn(K_ibs,basis)
      write(6,*)' ibs_cn2ucn size      :',ibs_cn2ucn(SZ_ibs,basis)
      write(6,*)' '
      write(6,*)' ibs_cn2ce ma handle  :',ibs_cn2ce(H_ibs,basis)
      write(6,*)' ibs_cn2ce ma pointer :',ibs_cn2ce(K_ibs,basis)
      write(6,*)' ibs_cn2ce size       :',ibs_cn2ce(SZ_ibs,basis)
      write(6,*)' '
      write(6,*)' ibs_ce2uce ma handle :',ibs_ce2uce(H_ibs,basis)
      write(6,*)' ibs_ce2uce ma pointer:',ibs_ce2uce(K_ibs,basis)
      write(6,*)' ibs_ce2uce size      :',ibs_ce2uce(SZ_ibs,basis)
      write(6,*)' '
      write(6,*)' ibs_cn2bfr ma handle :',ibs_cn2bfr(H_ibs,basis)
      write(6,*)' ibs_cn2bfr ma pointer:',ibs_cn2bfr(K_ibs,basis)
      write(6,*)' ibs_cn2bfr size      :',ibs_cn2bfr(SZ_ibs,basis)
      write(6,*)' '
      write(6,*)' ibs_ce2cnr ma handle :',ibs_ce2cnr(H_ibs,basis)
      write(6,*)' ibs_ce2cnr ma pointer:',ibs_ce2cnr(K_ibs,basis)
      write(6,*)' ibs_ce2cnr size      :',ibs_ce2cnr(SZ_ibs,basis)
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '
c
      status = bas_numcont(basisin,ncont)
      status = geom_ncent(geom,nat)
      write(6,*)' total number of contractions is ',ncont
      write(6,*)' total number of centers      is ',nat
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '
      write(6,*)' cont 2 uniq cont check '
      do icont = 0,ncont
        ucont = 23234
        ucont2 = ucont
        ucont = sf_ibs_cn2ucn(icont,basis)
        ucont2 = int_mb(mb_ibs_cn2ucn(icont,basis))
        write(6,'(a,i5,a,i5,i5)')
     &        ' contraction ',icont,
     &        ' maps to unique contraction ',ucont,ucont2
      enddo

      write(6,*)' '
      write(6,*)' '
      write(6,*)' cont 2 center check '
      do icont = 0,ncont
        status = bas_cn2ce(basisin,icont,cent_api)
        cent_sf = sf_ibs_cn2ce(icont,basis)
        cent_mb = int_mb(mb_ibs_cn2ce(icont,basis))
        write(6,'(a,i5,a,3i5,a)')
     &        ' contraction ',icont,' is on center ',
     &        cent_api,cent_sf,cent_mb,
     &        ' ... (sould be 3 identical numbers)'
      enddo
c
      write(6,*)' '
      write(6,*)' '
      write(6,*)' center 2 unique center check '
      do inat = 1,nat
        atom_sf = sf_ibs_ce2uce(inat,basis)
        atom_mb = int_mb(mb_ibs_ce2uce(inat,basis))
        write(6,'(a,i5,a,2i5,a)')
     &        ' center ',inat,
     &        ' maps to unique center/tag ',atom_sf,atom_mb,
     &        ' ... (sould be 2 identical numbers)'
      enddo
c
      write(6,*)' '
      write(6,*)' '
      write(6,*)' contraction 2 basis function range '
      do icont = 1,ncont
        status = bas_cn2bfr(basisin,icont,bflo_api,bfhi_api)
        bflo_sf = sf_ibs_cn2bfr(1,icont,basis)
        bfhi_sf = sf_ibs_cn2bfr(2,icont,basis)
        bflo_mb = int_mb(mb_ibs_cn2bfr(1,icont,basis))
        bfhi_mb = int_mb(mb_ibs_cn2bfr(2,icont,basis))
        write(6,'(a,i5,a,6i5,a)')
     &        ' contraction ',icont,
     &        ' maps to range of basis functions ',
     &        bflo_api,bfhi_api,bflo_sf,bfhi_sf,bflo_mb,bfhi_mb,
     &        ' ... (sould be 3 sets identical numbers)'
      enddo
c
      write(6,*)' '
      write(6,*)' '
      write(6,*)' center 2 contraction range '
      do inat = 1,nat
        status = bas_ce2cnr(basisin,inat,cnlo_api,cnhi_api)
        cnlo_sf = sf_ibs_ce2cnr(1,inat,basis)
        cnhi_sf = sf_ibs_ce2cnr(2,inat,basis)
        cnlo_mb = int_mb(mb_ibs_ce2cnr(1,inat,basis))
        cnhi_mb = int_mb(mb_ibs_ce2cnr(2,inat,basis))
        write(6,'(a,i5,a,6i5,a)')
     &        ' center ',inat,
     &        ' maps to range of contractons',
     &        cnlo_api,cnhi_api,cnlo_sf,cnhi_sf,cnlo_mb,cnhi_mb,
     &        ' ... (sould be 3 sets identical numbers)'
      enddo
c
      end
