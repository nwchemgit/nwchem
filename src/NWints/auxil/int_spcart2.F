      subroutine spcart_get_fullcs_ma(basis,ma_hp,ncart,nsph)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "spcartP.fh"
#include "global.fh"
#include "stdio.fh"
*:functions
      logical spcart_init
      external spcart_init
*
* generate a matrix in local "ma" memory (heap)
* matrix(cartesian,spherical)
*
*:passed::
      integer basis ! [input] basis set handle
      integer ma_hp ! [output] ma handle of full matrix
      integer ncart ! [output] cartesian dimension
      integer nsph  ! [output] spherical dimension
*:local::
      integer ma_hk ! index into ma for array
      logical o_loc_debug
      integer highang
      integer ishell
      integer numcont
      integer type
      integer np, ng, sphcar
c
      o_loc_debug = (.true.).and. (ga_nodeid().eq.0)
      ma_hp = 0
      ncart = 0
      nsph  = 0
c
      if (.not.bas_check_handle(basis,'spcart_get_fullcs_ma'))
     &      call errquit
     &      ('spcart_get_fullcs_ma:invalid basis handle=>',basis)
c
      if (.not.bas_high_angular(basis,highang)) call errquit
     &      ('spcart_get_fullcs_ma:spcart_init failed',911)
      call defNxyz(highang)
      if (.not.spcart_init(highang,.true.,.false.))call errquit
     &      ('spcart_get_fullcs_ma:spcart_init failed',911)
c
      if (.not.bas_numcont(basis,numcont)) call errquit
     &      ('spcart_get_fullcs_ma:bas_numcont failed',911)
      do ishell = 1,numcont
        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
     &        call errquit
     &      ('spcart_get_fullcs_ma:bas_continfo failed on shell',
     &        ishell)
        ncart = ncart + ng*((type+1)*(type+2)/2)
        nsph  = nsph  + ng*(2*type + 1)
      enddo
c
      if (o_loc_debug) then
        write(luout,*)' cartesian dimension :',ncart
        write(luout,*)' spherical dimension :',nsph
        call util_flush(luout)
      endif
c
      if (.not.ma_alloc_get(mt_dbl,(ncart*nsph),
     &      'fullma cart,sph matrix',ma_hp,ma_hk))
     &      call errquit
     &      ('spcart_get_fullcs_ma:ma_alloc_get failed for size',
     &      (ncart*nsph))
c
      call spcart_fill_fullcs_ma(numcont,ncart,nsph,dbl_mb(ma_hk),basis,
     &      o_loc_debug)
c
      if (o_loc_debug) call output
     &      (dbl_mb(ma_hk),1,ncart,1,nsph,ncart,nsph,1)
c
      end
      subroutine spcart_fill_fullcs_ma(numcont,ncart,nsph,matrx,
     &      basis,debug)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "spcartP.fh"
*::passed:
      integer numcont
      integer ncart
      integer nsph
      double precision matrx(ncart,nsph)
      integer basis
      logical debug
*::local:
      integer ishell
      integer type
      integer np, ng, sphcar
      integer off_cart, off_sph
      integer l2c, l2s
      integer ig, is, iss, ic
      integer ims, imc
*      
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
      call dfill((ncart*nsph),0.0d00,matrx,1)
c
      off_cart = 0
      off_sph  = 0
      do ishell = 1,numcont
        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
     &        call errquit
     &      ('spcart_fill_fullcs_ma:bas_continfo failed on shell',
     &        ishell)
        l2c = (type+1)*(type+2)/2
        l2s = 2*type + 1
        do ig = 1,ng
          do is = 1,l2s
            iss = is-type-1
            ims = off_sph  + is
            do ic = 1,l2c
              imc = off_cart + ic
              matrx(imc,ims) = Dtrans(ic,iss,type)
            enddo
          enddo
          off_cart = off_cart + l2c
          off_sph  = off_sph  + l2s
        enddo
      enddo
      end

