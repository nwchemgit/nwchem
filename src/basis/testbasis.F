      program testbasis
c $Id: testbasis.F,v 1.9 1994-06-01 19:03:38 d3g681 Exp $
      implicit none
c
#include "bas.fh"
#include "rtdb.fh"
#include "mafdecls.h"
#include "geom.fh"
#include "inp.fh"
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
      integer ncenters, sphcart, i,j
      integer nbf,ibf,iat,icont
      character*256 namename, nametrans
      character*16 drivtags(20)
      double precision coords(3,20), charge(20)
      double precision exp(400), coeff(400)
      logical status
      double precision expnt_new(3), coeff_new(4,3)
      logical bas_add_ucnt
      external bas_add_ucnt
      data expnt_new/1., 2., 3./
      data coeff_new/-1., -2., -3., 0.0, -4., -5., -6., 0.0, -7., -8.,
     $     -9., 0.0/
c
      if (.not. ma_init(MT_DBL, -1, -1)) stop 99
      status = rtdb_open('shit.rtdb','unknown', rtdb)
c
      write(6,*)' rtdb handle ', rtdb
      status = bas_321g_load(rtdb)
c
      if(.not.geom_create(geom,'321g:1-20')) then
        write(6,*)' error getting geometry handle '
        stop ' error '
      endif

c      supply geom information.  
      ncenters = 3
c..  oxygen
      drivtags(1) = 'O'
      coords(1,1) = 0.0d00
      coords(2,1) = 0.0d00
      coords(3,1) = 0.0d00
      charge(1)   = 8.0d00
c.. hydrogen 1
      drivtags(2) = 'H'
      coords(1,2) = 1.0d00
      coords(2,2) = 1.0d00
      coords(3,2) = 1.0d00
      charge(2)   = 1.0d00
c.. hydrogen 2
      drivtags(3) = 'H'
      coords(1,3) = -1.0d00
      coords(2,3) = -1.0d00
      coords(3,3) = -1.0d00
      charge(3)   = 1.0d00
c
      if(.not.geom_cart_set(geom,ncenters,drivtags,coords,charge))then
        write(6,*)' geom_cart_set fail '
        stop ' error'
      else
        status=geom_print(geom)
      endif

      basis = 0
      if(.not.bas_create(basis,'3211-20')) then
        write(6,*)' error getting basis handle '
        stop ' error '
      endif
      basis = 0
      if(.not.bas_create(basis,'321g1-20')) then
        write(6,*)' error getting second basis handle '
        stop ' error '
      endif
      basis = 0
      if(.not.bas_create(basis,'321g:1-20')) then
        write(6,*)' error getting third basis handle '
        stop ' error '
      endif
c
      status = bas_rtdb_load(rtdb,geom,basis,'321g:1-20')
      status = bas_print(basis)
      status = gbs_map_print(basis)
c
      status = bas_continfo(basis,1,iang,nprim,ngen,sphcart)
      write(6,*)' f:query: nprim   cont 1 ',nprim
      write(6,*)' f:query: ngen    cont 1 ',ngen
      write(6,*)' f:query: sphcart cont 1 ',sphcart
      status = bas_get_exponent(basis,1,exp)
      status = bas_get_coeff(basis,1,coeff)
      write(6,*)' exponenents and coefficients '
      do 00100 i=1,nprim
        write(6,*)exp(i),(coeff(i+(j-1)*nprim),j=1,ngen)
00100 continue
c
      status = bas_ucontinfo(basis,1,nprim,ngen,sphcart)
      write(6,*)' t:query: nprim   cont 1 ',nprim
      write(6,*)' t:query: ngen    cont 1 ',ngen
      write(6,*)' t:query: sphcart cont 1 ',sphcart
      status = bas_getu_exponent(basis,1,exp)
      status = bas_getu_coeff(basis,1,coeff)
      write(6,*)' exponenents and coefficients '
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
      status = bas_ucontinfo(basis,1,nprim,ngen,sphcart)
      write(6,*)' modified '
      write(6,*)' t:query: nprim   cont 1 ',nprim
      write(6,*)' t:query: ngen    cont 1 ',ngen
      write(6,*)' t:query: sphcart cont 1 ',sphcart
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
     $     coeff_new, 4)) write(6,*) ' basis_add_ucnt failed'
      if (.not. bas_print(basis)) write(6,*) ' print ???'
c
c     Try adding new contractions on a new center
c
      write(6,*) ' adding 2*3 g function on Cl'
      if (.not. bas_add_ucnt(basis, 'Cl', 4, 2, 3, expnt_new,
     $     coeff_new, 4)) write(6,*) ' basis_add_ucnt failed'
      if (.not. bas_print(basis)) write(6,*) ' print ???'
c
c
      write(6,'(///,1x,a)')' bas_print_all '
      status = bas_print_all()
      call bas_err_info(' who who who ')
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
      status = bas_nbf_cn_max(basis,nbf)
      write(6,*)' max block of nbf on a shell for basis is:',nbf
c
      status = bas_nbf_ce_max(basis,nbf)
      write(6,*)' max block of nbf on an atom for basis is:',nbf
c
      stop ' testbasis done '
      end

