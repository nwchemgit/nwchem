      logical function spcart_get_full_cs(matrix_type,handle,ncart,nsph,
     &    use_ma,use_ga,basis)
      implicit none
*
* $Id$
*
*
*::include files::
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "spcartP.fh"
#include "bas.fh"
#include "util.fh"
#include "errquit.fh"
*
*::functions::
      logical spcart_init, spcart_terminate
      external spcart_init, spcart_terminate
*
*::passed::
      character*(*) matrix_type ! [input] cart/sper, transpose, or inverse
* matrix_type = 'N' means return the normal cartesian/spherical tranformation
*                   matrix for the given basis set: -> matrix(cartesian,spherical)
* matrix_type = 'I' means give me the true inverse matrix  -> matrix(spherical,cartesian)
*                   e.g., matrix('I')*matrix('N') = unit matrix
* matrix_type = 'T' means give me the transpose of the normal matrix('N')
*
      integer handle ! [output] the MA(heap) or GA handle that holds the matrix
      integer ncart  ! [output] the dimension of the cartesian space of the matrix
      integer nsph   ! [output] the dimension of the spherical space of the matrix
      logical use_ma ! [input] return the matrix in MA based memory (heap)
      logical use_ga ! [input] return the matrix in GA based memory
      integer basis  ! [input] the basis set used to determine the matrix
*
*::local::
      logical regular       ! compute regular matrix
      logical inverse       ! compute inverse matrix
      logical transpose     ! compute transpose of regular matrix
      logical o_debug       ! local debug variable
      integer highang       ! highest angular momentum in the given basis set
      integer numcont       ! number of shells/contractions in the basis set
      integer ishell        ! loop index over shells
      integer type          ! generic shell type
      integer np            ! generic shell number of primitives
      integer ng            ! generic shell number of general contractions in the shell
      integer sphcar        ! 1 == shell is in a spherical basis set; 0== cartesian
      integer ma_hk         ! index for handle when handle is from ma
      integer ga_reg_handle ! compute regular matrix with ga then ga_transpose
      character*255 errstring      ! temporary error string
      logical do_spcart_terminate  ! variable to keep track if spcart was initialized
*
* initialize locals
      inverse   = inp_compare(.false.,'i',matrix_type(1:1))
      transpose = inp_compare(.false.,'t',matrix_type(1:1))
      regular   = inp_compare(.false.,'n',matrix_type(1:1))
      regular   = regular.or.inp_compare(.false.,'r',matrix_type(1:1))
      if (.not.(regular.or.transpose.or.inverse)) then
        write(luout,*)
     &      ' spcart_get_full_cs: called with a matrix type of <<< ',
     &      matrix_type,' >>>'
        write(luout,*)
     &      '    which is not understood.  Producing by default the '
        write(luout,*)
     &      '    normal or regular matrix '
        regular = .true.
      endif
      o_debug   = .false.
      do_spcart_terminate = sph_cart_init.ne.SPH_CART_INIT_VALUE
* initialize outputs      
      handle = 0
      ncart  = 0
      nsph   = 0
      spcart_get_full_cs = .true.
*
      if (.not.bas_check_handle(basis,'spcart_get_fullcs_ma'))
     &      call errquit
     &      ('spcart_get_full_cs:invalid basis handle=>',basis,
     &       BASIS_ERR)
c
      if (.not.bas_high_angular(basis,highang)) call errquit
     &      ('spcart_get_full_cs:bas_high_angular failed',911,
     &       BASIS_ERR)
      if (.not.spcart_init(highang,.true.,.false.))call errquit
     &      ('spcart_get_full_cs:spcart_init failed',911,
     &       UNKNOWN_ERR)
c
      if (.not.bas_numcont(basis,numcont)) call errquit
     &      ('spcart_get_full_cs:bas_numcont failed',911, BASIS_ERR)
      do ishell = 1,numcont
        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
     &        call errquit
     &      ('spcart_get_full_cs:bas_continfo failed on shell',
     &        ishell, BASIS_ERR)
        if (type.ge.0) then
          ncart = ncart + ng*((type+1)*(type+2)/2)
          nsph  = nsph  + ng*(2*type + 1)
        else if (type.eq.-1) then
          ncart = ncart + 4
          nsph  = nsph + 4
        else
          spcart_get_full_cs = .false.
          call errquit
     &      ('spcart_get_full_cs: unrecognized type',type, INT_ERR)
        endif
      enddo
c
      if (o_debug) then
        write(luout,*)' number of shells    :1:',numcont
        write(luout,*)' cartesian dimension :1:',ncart
        write(luout,*)' spherical dimension :1:',nsph
        call util_flush(luout)
      endif
*
      if (use_ma.and.use_ga) then
        errstring = 'spcart_get_full_cs: FATAL ERROR: '
     &      //'use_ma and use_ga are both true!'
        call errquit(errstring,911, MEM_ERR)
      else if (use_ma) then
        spcart_get_full_cs = spcart_get_full_cs.and.
     &      ma_alloc_get(mt_dbl,(ncart*nsph),
     &      'full ma-based cart,sph matrix',handle,ma_hk)
        if (.not.spcart_get_full_cs) call errquit(
     &      'spcart_get_full_cs: ma_alloc_get failed',911, MA_ERR)
        if (regular) then
          call spcart_fill_fullcs_ma(numcont,ncart,nsph,
     &        dbl_mb(ma_hk),basis,o_debug)
          if (o_debug) then
            write(luout,*) ' regular cartesian/spherical matrix:ma'
            call output
     &          (dbl_mb(ma_hk),1,ncart,1,nsph,ncart,nsph,1)
          endif
        else if (inverse) then
          call spcart_fill_fullinvrscs_ma(numcont,ncart,nsph,
     &        dbl_mb(ma_hk),basis,o_debug)
          if (o_debug) then
            write(luout,*) ' inverse cartesian/spherical matrix:ma'
            call output
     &          (dbl_mb(ma_hk),1,nsph,1,ncart,nsph,ncart,1)
          endif
        else if (transpose) then
          call spcart_fill_fulltranscs_ma(numcont,ncart,nsph,
     &        dbl_mb(ma_hk),basis,o_debug)
          if (o_debug) then
            write(luout,*) ' transpose of regular cartesian/',
     &          'spherical matrix:ma'
            call output
     &          (dbl_mb(ma_hk),1,nsph,1,ncart,nsph,ncart,1)
          endif
        else
          errstring = 'spcart_get_full_cs: '
     &        //'invalid matrix type at MA alloc point'
          call errquit(errstring,911, MA_ERR)
        endif
      else if (use_ga) then
        if (regular) then
          spcart_get_full_cs = spcart_get_full_cs.and.
     &        ga_create(mt_dbl,ncart,nsph,
     &        'regular ga-based cart/sph matrix',0,0,handle)
          if (.not.spcart_get_full_cs) call errquit(
     &        'spcart_get_full_cs: ga_create failed ',911, GA_ERR)
          call spcart_fill_cs_ga(numcont,ncart,nsph,
     &        handle,basis,.false.,o_debug)
          if (o_debug) then
            write(luout,*) ' regular cartesian/spherical matrix:ga'
            call ga_print(handle)
          endif
        else if (inverse) then
          spcart_get_full_cs = spcart_get_full_cs.and.
     &        ga_create(mt_dbl,nsph,ncart,
     &        'inverse ga-based cart/sph matrix',0,0,handle)
          if (.not.spcart_get_full_cs) call errquit(
     &        'spcart_get_full_cs: ga create failed',911, GA_ERR)
          call spcart_fill_cs_ga(numcont,ncart,nsph,
     &        handle,basis,.true.,o_debug)
          if (o_debug) then
            write(luout,*) ' inverse cartesian/spherical matrix:ga'
            call ga_print(handle)
          endif
        else if (transpose) then
          spcart_get_full_cs = spcart_get_full_cs.and.
     &        ga_create(mt_dbl,ncart,nsph,
     &        'regular for transpose ga-based cart/sph matrix',
     &        0,0,ga_reg_handle)
          call ga_zero(ga_reg_handle)
          spcart_get_full_cs = spcart_get_full_cs.and.
     &        ga_create(mt_dbl,nsph,ncart,
     &        'regular for transpose ga-based cart/sph matrix',
     &        0,0,handle)
          call ga_zero(handle)
          if (.not.spcart_get_full_cs) call errquit(
     &        'spcart_get_full_cs:one of two ga_create calls failed',
     &        911, GA_ERR)
          call spcart_fill_cs_ga(numcont,ncart,nsph,
     &        ga_reg_handle,basis,.false.,o_debug)
          call ga_transpose(ga_reg_handle,handle)
          if (o_debug) then
            write(luout,*) ' for transpose the regular cartesian/sp',
     &          'herical matrix:ga'
            call ga_print(ga_reg_handle)
            write(luout,*) ' transpose of regular cartesian/',
     &          'spherical matrix:ma'
            call ga_print(handle)
          endif
          call ga_sync()
          if (.not.ga_destroy(ga_reg_handle))call errquit(
     &        'spcart_get_full_cs:ga_destroy for reg_handle failed',
     &        911, GA_ERR)
        else
          errstring = 'spcart_get_full_cs: '
     &        //'invalid matrix type at GA create point'
          call errquit(errstring,911, GA_ERR)
        endif
      else
        errstring = 'spcart_get_full_cs: FATAL ERROR: '
     &      //'use_ma and use_ga are both false'
        call errquit(errstring,911, MEM_ERR)
      endif
c
      if (do_spcart_terminate) then
        if (.not.spcart_terminate()) call errquit
     &      ('spcart_get_full_cs:spcart_terminate failed ',911,
     &       UNKNOWN_ERR)
      endif
      end
      subroutine spcart_fill_fullcs_ma(numcont,ncart,nsph,matrx,
     &      basis,debug)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "spcartP.fh"
#include "stdio.fh"
#include "errquit.fh"
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
      integer cart_offset, sph_offset
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
      call dcopy((ncart*nsph),0.0d00,0,matrx,1)
c
      if (debug) then
        write(luout,*)' number of shells    :2:',numcont
        write(luout,*)' cartesian dimension :2:',ncart
        write(luout,*)' spherical dimension :2:',nsph
        call util_flush(luout)
      endif
      cart_offset = 0
      sph_offset  = 0
      do ishell = 1,numcont
        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
     &        call errquit
     &      ('spcart_fill_fullcs_ma:bas_continfo failed on shell',
     &        ishell, BASIS_ERR)
        if (type.ge.0) then
          l2c = (type+1)*(type+2)/2
          l2s = 2*type + 1
        else if(type.eq.-1) then
          l2c = 4
          l2s = 4
          ng  = 1
        else
          call errquit
     &      ('spcart_fill_fullcs_ma: unrecognized type',type,
     &       INT_ERR)
        endif
        do ig = 1,ng
          if (debug) then
            write(luout,*)' matrx cart dim ',(cart_offset+1),
     &          ':',(cart_offset+l2c),'(l2c=',l2c,')'
            write(luout,*)' matrx sphr dim ',(sph_offset+1),
     &          ':',(sph_offset+l2s),'(l2s=',l2s,')'
          endif
          do is = 1,l2s
            iss = is-type-1
            ims = sph_offset  + is
            do ic = 1,l2c
              imc = cart_offset + ic
              if (type.ge.2) then
                matrx(imc,ims) = Dtrans(ic,iss,type)
              else if (type.eq.-1.or.type.eq.0.or.type.eq.1) then
*..................................         then unit matrix added
                if (ic.eq.is) matrx(imc,ims) = 1.0d00
                if (ic.ne.is) matrx(imc,ims) = 0.0d00
              else
                call errquit(
     &              'spcart_fill_fullcs_ma: fatal type error',911,
     &       MA_ERR)
              endif
*rak:              if (debug)
*rak:     &            write(luout,*) 'matrx <',imc,',',ims,'> = ',
*rak:     &            matrx(imc,ims)
            enddo
          enddo
          cart_offset = cart_offset + l2c
          sph_offset  = sph_offset  + l2s
        enddo
      enddo
      end
      subroutine spcart_fill_fulltranscs_ma(numcont,ncart,nsph,matrx,
     &      basis,debug)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "spcartP.fh"
#include "stdio.fh"
#include "errquit.fh"
*::passed:
      integer numcont
      integer ncart
      integer nsph
      double precision matrx(nsph,ncart)
      integer basis
      logical debug
*::local:
      integer ishell
      integer type
      integer np, ng, sphcar
      integer cart_offset, sph_offset
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
      call dcopy((ncart*nsph),0.0d00,0,matrx,1)
c
      if (debug) then
        write(luout,*)' number of shells    :2:',numcont
        write(luout,*)' cartesian dimension :2:',ncart
        write(luout,*)' spherical dimension :2:',nsph
        call util_flush(luout)
      endif
      cart_offset = 0
      sph_offset  = 0
      do ishell = 1,numcont
        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
     &        call errquit
     &      ('spcart_fill_fulltranscs_ma:bas_continfo failed:shell=',
     &        ishell, BASIS_ERR)
        if (type.ge.0) then
          l2c = (type+1)*(type+2)/2
          l2s = 2*type + 1
        else if(type.eq.-1) then
          l2c = 4
          l2s = 4
          ng  = 1
        else
          call errquit
     &      ('spcart_fill_fulltranscs_ma: unrecognized type',type,
     &       INT_ERR)
        endif
        do ig = 1,ng
          if (debug) then
            write(luout,*)' matrx cart dim ',(cart_offset+1),
     &          ':',(cart_offset+l2c),'(l2c=',l2c,')'
            write(luout,*)' matrx sphr dim ',(sph_offset+1),
     &          ':',(sph_offset+l2s),'(l2s=',l2s,')'
          endif
          do ic = 1,l2c
            imc = cart_offset + ic
            do is = 1,l2s
              iss = is-type-1
              ims = sph_offset  + is
              if (type.ge.2) then
                matrx(ims,imc) = Dtrans(ic,iss,type)
              else if (type.eq.-1.or.type.eq.0.or.type.eq.1) then
*..................................         then unit matrix added
                if (ic.eq.is) matrx(ims,imc) = 1.0d00
                if (ic.ne.is) matrx(ims,imc) = 0.0d00
              else
                call errquit(
     &              'spcart_fill_fulltranscs_ma: fatal type error',
     &              911, INT_ERR)
              endif
*rak:              if (debug)
*rak:     &            write(luout,*) 'matrx <',ims,',',imc,'> = ',
*rak:     &            matrx(ims,imc)
            enddo
          enddo
          cart_offset = cart_offset + l2c
          sph_offset  = sph_offset  + l2s
        enddo
      enddo
      end
      subroutine spcart_fill_fullinvrscs_ma(numcont,ncart,nsph,matrx,
     &      basis,debug)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "spcartP.fh"
#include "stdio.fh"
#include "errquit.fh"
*::passed:
      integer numcont
      integer ncart
      integer nsph
      double precision matrx(nsph,ncart)
      integer basis
      logical debug
*::local:
      integer ishell
      integer type
      integer np, ng, sphcar
      integer cart_offset, sph_offset
      integer l2c, l2s
      integer ig, is, iss, ic
      integer ims, imc
*      
c::statement function ----- start
      integer iic,iis,iil
      double precision Dinvtrans
      Dinvtrans(iis,iic,iil) =
     &    dbl_mb((int_mb(k_invsp2c_lindx+iil))+
     &    (iic-1)*(2*iil+1)+iis+iil)
c::statement function ----- end
c
      call dcopy((ncart*nsph),0.0d00,0,matrx,1)
c
      if (debug) then
        write(luout,*)' number of shells    :2:',numcont
        write(luout,*)' cartesian dimension :2:',ncart
        write(luout,*)' spherical dimension :2:',nsph
        call util_flush(luout)
      endif
      cart_offset = 0
      sph_offset  = 0
      do ishell = 1,numcont
        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
     &        call errquit
     &      ('spcart_fill_fullinvrscs_ma:bas_continfo failed:shell=',
     &        ishell, BASIS_ERR)
        if (type.ge.0) then
          l2c = (type+1)*(type+2)/2
          l2s = 2*type + 1
        else if(type.eq.-1) then
          l2c = 4
          l2s = 4
          ng  = 1
        else
          call errquit
     &      ('spcart_fill_fullinvrscs_ma: unrecognized type',type,
     &       INT_ERR)
        endif
        do ig = 1,ng
          if (debug) then
            write(luout,*)' matrx cart dim ',(cart_offset+1),
     &          ':',(cart_offset+l2c),'(l2c=',l2c,')'
            write(luout,*)' matrx sphr dim ',(sph_offset+1),
     &          ':',(sph_offset+l2s),'(l2s=',l2s,')'
          endif
          do ic = 1,l2c
            imc = cart_offset + ic
            do is = 1,l2s
              iss = is-type-1
              ims = sph_offset  + is
              if (type.ge.2) then
                matrx(ims,imc) = Dinvtrans(iss,ic,type)
              else if (type.eq.-1.or.type.eq.0.or.type.eq.1) then
*..................................         then unit matrix added
                if (ic.eq.is) matrx(ims,imc) = 1.0d00
                if (ic.ne.is) matrx(ims,imc) = 0.0d00
              else
                call errquit(
     &              'spcart_fill_fullinvrscs_ma: fatal type error',
     &              911, INT_ERR)
              endif
*rak:              if (debug)
*rak:     &            write(luout,*) 'matrx <',ims,',',imc,'> = ',
*rak:     &            matrx(ims,imc)
            enddo
          enddo
          cart_offset = cart_offset + l2c
          sph_offset  = sph_offset  + l2s
        enddo
      enddo
      end
      subroutine spcart_fill_cs_ga(numcont,ncart,nsph,
     &    ga_matrx,basis,inverse,debug)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "spcartP.fh"
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
*::passed:
      integer numcont
      integer ncart
      integer nsph
      integer ga_matrx
      integer basis
      logical debug
      logical inverse
*::local:
      double precision ident4(4,4), ident3(3,3)
      integer ishell
      integer type
      integer np, ng, sphcar
      integer cart_offset, sph_offset
      integer l2c, l2s
      integer ig, ilo, ihi, jlo, jhi, ld
      integer ma_index_ptr
*      
      call dcopy((4*4),0.0d00,0,ident4,1)
      call dcopy((3*3),0.0d00,0,ident3,1)
      do ig = 1,3
        ident4(ig,ig)=1.0d00
        ident3(ig,ig)=1.0d00
      enddo
      ident4(4,4) = 1.0d00
c
      call ga_zero(ga_matrx)
c
      if (debug) then
        write(luout,*)' number of shells    :2:',numcont
        write(luout,*)' cartesian dimension :2:',ncart
        write(luout,*)' spherical dimension :2:',nsph
        call util_flush(luout)
      endif
      cart_offset = 0
      sph_offset  = 0
      do ishell = 1,numcont
        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
     &        call errquit
     &      ('spcart_fill_cs_ga:bas_continfo failed on shell',
     &        ishell, BASIS_ERR)
        if (type.ge.0) then
          l2c = (type+1)*(type+2)/2
          l2s = 2*type + 1
        else if(type.eq.-1) then
          l2c = 4
          l2s = 4
          ng  = 1
        else
          call errquit
     &      ('spcart_fill_cs_ga: unrecognized type',type, GA_ERR)
        endif
        do ig = 1,ng
          if (debug) then
            write(luout,*)' matrx cart dim ',(cart_offset+1),
     &          ':',(cart_offset+l2c),'(l2c=',l2c,')'
            write(luout,*)' matrx sphr dim ',(sph_offset+1),
     &          ':',(sph_offset+l2s),'(l2s=',l2s,')'
          endif
          if (inverse) then
            ma_index_ptr = int_mb(k_invsp2c_lindx+type)
            ilo = sph_offset+1
            ihi = sph_offset+l2s
            jlo = cart_offset+1
            jhi = cart_offset+l2c
            ld  = l2s
          else
            ma_index_ptr = int_mb(k_sp2c_lindx+type)
            ilo = cart_offset+1
            ihi = cart_offset+l2c
            jlo = sph_offset+1
            jhi = sph_offset+l2s
            ld  = l2c
          endif
          if (type.ge.2) then
            if (ga_nodeid().eq.0)
     &          call ga_put(ga_matrx,
     &          ilo,ihi,jlo,jhi,
     &          dbl_mb(ma_index_ptr),ld)
          else if (type.eq.-1) then
            if (ga_nodeid().eq.0)
     &          call ga_put(ga_matrx,
     &          ilo,ihi,jlo,jhi,
     &          ident4,ld)
          else if (type.eq.0.or.type.eq.1) then
            if (ga_nodeid().eq.0)
     &          call ga_put(ga_matrx,
     &          ilo,ihi,jlo,jhi,
     &          ident3,ld)
          else
            call errquit(
     &          'spcart_fill_cs_ga: fatal type error',911, GA_ERR)
          endif
          if (debug)
     &        call ga_print_patch(ga_matrx,
     &        ilo,ihi,jlo,jhi,1)
*
          cart_offset = cart_offset + l2c
          sph_offset  = sph_offset  + l2s
        enddo
      enddo
      end
************************************************************************************ obsolete below
*obsolete:      subroutine spcart_get_fullcs_ma(basis,ma_hp,ncart,nsph)
*obsolete:      implicit none
*obsolete:#include "mafdecls.fh"
*obsolete:#include "bas.fh"
*obsolete:#include "spcartP.fh"
*obsolete:#include "global.fh"
*obsolete:#include "stdio.fh"
*obsolete:#include "errquit.fh"
*obsolete:*:functions
*obsolete:      logical spcart_init, spcart_terminate
*obsolete:      external spcart_init, spcart_terminate
*obsolete:*
*obsolete:* generate a matrix in local "ma" memory (heap)
*obsolete:* matrix(cartesian,spherical)
*obsolete:*
*obsolete:*:passed::
*obsolete:      integer basis ! [input] basis set handle
*obsolete:      integer ma_hp ! [output] ma handle of full matrix
*obsolete:      integer ncart ! [output] cartesian dimension
*obsolete:      integer nsph  ! [output] spherical dimension
*obsolete:*:local::
*obsolete:      integer ma_hk ! index into ma for array
*obsolete:      logical o_loc_debug
*obsolete:      integer highang
*obsolete:      integer ishell
*obsolete:      integer numcont
*obsolete:      integer type
*obsolete:      integer np, ng, sphcar
*obsolete:      logical do_terminate
*obsolete:c
*obsolete:      o_loc_debug = .false.
*obsolete:      ma_hp = 0
*obsolete:      ncart = 0
*obsolete:      nsph  = 0
*obsolete:c
*obsolete:      do_terminate = sph_cart_init.ne.SPH_CART_INIT_VALUE
*obsolete:c
*obsolete:      if (.not.bas_check_handle(basis,'spcart_get_fullcs_ma'))
*obsolete:     &      call errquit
*obsolete:     &      ('spcart_get_fullcs_ma:invalid basis handle=>',basis)
*obsolete:c
*obsolete:      if (.not.bas_high_angular(basis,highang)) call errquit
*obsolete:     &      ('spcart_get_fullcs_ma:spcart_init failed',911)
*obsolete:      if (.not.spcart_init(highang,.true.,.false.))call errquit
*obsolete:     &      ('spcart_get_fullcs_ma:spcart_init failed',911)
*obsolete:c
*obsolete:      if (.not.bas_numcont(basis,numcont)) call errquit
*obsolete:     &      ('spcart_get_fullcs_ma:bas_numcont failed',911)
*obsolete:      do ishell = 1,numcont
*obsolete:        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
*obsolete:     &        call errquit
*obsolete:     &      ('spcart_get_fullcs_ma:bas_continfo failed on shell',
*obsolete:     &        ishell)
*obsolete:        if (type.ge.0) then
*obsolete:          ncart = ncart + ng*((type+1)*(type+2)/2)
*obsolete:          nsph  = nsph  + ng*(2*type + 1)
*obsolete:        else if (type.eq.-1) then
*obsolete:          ncart = ncart + 4
*obsolete:          nsph  = nsph + 4
*obsolete:        else
*obsolete:          call errquit
*obsolete:     &      ('spcart_get_fullcs_ma: unrecognized type',type)
*obsolete:        endif
*obsolete:      enddo
*obsolete:c
*obsolete:      if (o_loc_debug) then
*obsolete:        write(luout,*)' number of shells    :1:',numcont
*obsolete:        write(luout,*)' cartesian dimension :1:',ncart
*obsolete:        write(luout,*)' spherical dimension :1:',nsph
*obsolete:        call util_flush(luout)
*obsolete:      endif
*obsolete:c
*obsolete:      if (.not.ma_alloc_get(mt_dbl,(ncart*nsph),
*obsolete:     &      'fullma cart,sph matrix',ma_hp,ma_hk))
*obsolete:     &      call errquit
*obsolete:     &      ('spcart_get_fullcs_ma:ma_alloc_get failed for size',
*obsolete:     &      (ncart*nsph))
*obsolete:c
*obsolete:      call spcart_fill_fullcs_ma(numcont,ncart,nsph,dbl_mb(ma_hk),
*obsolete:     &      basis,o_loc_debug)
*obsolete:c
*obsolete:      if (o_loc_debug) call output
*obsolete:     &      (dbl_mb(ma_hk),1,ncart,1,nsph,ncart,nsph,1)
*obsolete:c
*obsolete:      if (do_terminate) then
*obsolete:        if (.not.spcart_terminate()) call errquit
*obsolete:     &      ('spcart_get_fullcs_ma:spcart_terminate failed ',911)
*obsolete:      endif
*obsolete:      end
*obsolete:      subroutine spcart_get_fullinvrscs_ma(basis,ma_hp,ncart,nsph)
*obsolete:      implicit none
*obsolete:#include "mafdecls.fh"
*obsolete:#include "bas.fh"
*obsolete:#include "spcartP.fh"
*obsolete:#include "global.fh"
*obsolete:#include "stdio.fh"
*obsolete:#include "errquit.fh"
*obsolete:*:functions
*obsolete:      logical spcart_init, spcart_terminate
*obsolete:      external spcart_init, spcart_terminate
*obsolete:*
*obsolete:* generate a matrix in local "ma" memory (heap)
*obsolete:* matrix(spherical,cartesian)
*obsolete:*
*obsolete:*:passed::
*obsolete:      integer basis ! [input] basis set handle
*obsolete:      integer ma_hp ! [output] ma handle of full matrix
*obsolete:      integer ncart ! [output] cartesian dimension
*obsolete:      integer nsph  ! [output] spherical dimension
*obsolete:*:local::
*obsolete:      integer ma_hk ! index into ma for array
*obsolete:      logical o_loc_debug
*obsolete:      integer highang
*obsolete:      integer ishell
*obsolete:      integer numcont
*obsolete:      integer type
*obsolete:      integer np, ng, sphcar
*obsolete:      logical do_terminate
*obsolete:c
*obsolete:      o_loc_debug = .false.
*obsolete:      ma_hp = 0
*obsolete:      ncart = 0
*obsolete:      nsph  = 0
*obsolete:c
*obsolete:      do_terminate = sph_cart_init.ne.SPH_CART_INIT_VALUE
*obsolete:c
*obsolete:      if (.not.bas_check_handle(basis,'spcart_get_fullinvrscs_ma'))
*obsolete:     &      call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma:invalid basis handle=>',basis)
*obsolete:c
*obsolete:      if (.not.bas_high_angular(basis,highang)) call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma:spcart_init failed',911)
*obsolete:      if (.not.spcart_init(highang,.true.,.false.))call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma:spcart_init failed',911)
*obsolete:c
*obsolete:      if (.not.bas_numcont(basis,numcont)) call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma:bas_numcont failed',911)
*obsolete:      do ishell = 1,numcont
*obsolete:        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
*obsolete:     &        call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma:bas_continfo failed on shell',
*obsolete:     &        ishell)
*obsolete:        if (type.ge.0) then
*obsolete:          ncart = ncart + ng*((type+1)*(type+2)/2)
*obsolete:          nsph  = nsph  + ng*(2*type + 1)
*obsolete:        else if (type.eq.-1) then
*obsolete:          ncart = ncart + 4
*obsolete:          nsph  = nsph + 4
*obsolete:        else
*obsolete:          call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma: unrecognized type',type)
*obsolete:        endif
*obsolete:      enddo
*obsolete:c
*obsolete:      if (o_loc_debug) then
*obsolete:        write(luout,*)' number of shells    :1:',numcont
*obsolete:        write(luout,*)' cartesian dimension :1:',ncart
*obsolete:        write(luout,*)' spherical dimension :1:',nsph
*obsolete:        call util_flush(luout)
*obsolete:      endif
*obsolete:c
*obsolete:      if (.not.ma_alloc_get(mt_dbl,(ncart*nsph),
*obsolete:     &      'fullma cart,sph matrix',ma_hp,ma_hk))
*obsolete:     &      call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma:ma_alloc_get failed for size',
*obsolete:     &      (ncart*nsph))
*obsolete:c
*obsolete:      call spcart_fill_fullinvrscs_ma(numcont,ncart,nsph,
*obsolete:     &      dbl_mb(ma_hk),basis,o_loc_debug)
*obsolete:c
*obsolete:      if (o_loc_debug) call output
*obsolete:     &      (dbl_mb(ma_hk),1,ncart,1,nsph,ncart,nsph,1)
*obsolete:c
*obsolete:      if (do_terminate) then
*obsolete:        if (.not.spcart_terminate()) call errquit
*obsolete:     &      ('spcart_get_fullinvrscs_ma:spcart_terminate failed ',911)
*obsolete:      endif
*obsolete:      end
*obsolete:      subroutine spcart_get_fulltranscs_ma(basis,ma_hp,ncart,nsph)
*obsolete:      implicit none
*obsolete:#include "mafdecls.fh"
*obsolete:#include "bas.fh"
*obsolete:#include "spcartP.fh"
*obsolete:#include "global.fh"
*obsolete:#include "stdio.fh"
*obsolete:#include "errquit.fh"
*obsolete:*:functions
*obsolete:      logical spcart_init, spcart_terminate
*obsolete:      external spcart_init, spcart_terminate
*obsolete:*
*obsolete:* generate a matrix in local "ma" memory (heap)
*obsolete:* matrix(spherical,cartesian)
*obsolete:*
*obsolete:*:passed::
*obsolete:      integer basis ! [input] basis set handle
*obsolete:      integer ma_hp ! [output] ma handle of full matrix
*obsolete:      integer ncart ! [output] cartesian dimension
*obsolete:      integer nsph  ! [output] spherical dimension
*obsolete:*:local::
*obsolete:      integer ma_hk ! index into ma for array
*obsolete:      logical o_loc_debug
*obsolete:      integer highang
*obsolete:      integer ishell
*obsolete:      integer numcont
*obsolete:      integer type
*obsolete:      integer np, ng, sphcar
*obsolete:      logical do_terminate
*obsolete:c
*obsolete:      o_loc_debug = .false.
*obsolete:      ma_hp = 0
*obsolete:      ncart = 0
*obsolete:      nsph  = 0
*obsolete:c
*obsolete:      do_terminate = sph_cart_init.ne.SPH_CART_INIT_VALUE
*obsolete:c
*obsolete:      if (.not.bas_check_handle(basis,'spcart_get_fulltranscs_ma'))
*obsolete:     &      call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma:invalid basis handle=>',basis)
*obsolete:c
*obsolete:      if (.not.bas_high_angular(basis,highang)) call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma:spcart_init failed',911)
*obsolete:      if (.not.spcart_init(highang,.true.,.false.))call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma:spcart_init failed',911)
*obsolete:c
*obsolete:      if (.not.bas_numcont(basis,numcont)) call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma:bas_numcont failed',911)
*obsolete:      do ishell = 1,numcont
*obsolete:        if (.not.bas_continfo(basis,ishell,type,np,ng,sphcar))
*obsolete:     &        call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma:bas_continfo failed on shell',
*obsolete:     &        ishell)
*obsolete:        if (type.ge.0) then
*obsolete:          ncart = ncart + ng*((type+1)*(type+2)/2)
*obsolete:          nsph  = nsph  + ng*(2*type + 1)
*obsolete:        else if (type.eq.-1) then
*obsolete:          ncart = ncart + 4
*obsolete:          nsph  = nsph + 4
*obsolete:        else
*obsolete:          call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma: unrecognized type',type)
*obsolete:        endif
*obsolete:      enddo
*obsolete:c
*obsolete:      if (o_loc_debug) then
*obsolete:        write(luout,*)' number of shells    :1:',numcont
*obsolete:        write(luout,*)' cartesian dimension :1:',ncart
*obsolete:        write(luout,*)' spherical dimension :1:',nsph
*obsolete:        call util_flush(luout)
*obsolete:      endif
*obsolete:c
*obsolete:      if (.not.ma_alloc_get(mt_dbl,(ncart*nsph),
*obsolete:     &      'fullma cart,sph matrix',ma_hp,ma_hk))
*obsolete:     &      call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma:ma_alloc_get failed for size',
*obsolete:     &      (ncart*nsph))
*obsolete:c
*obsolete:      call spcart_fill_fulltranscs_ma(numcont,ncart,nsph,
*obsolete:     &      dbl_mb(ma_hk),basis,o_loc_debug)
*obsolete:c
*obsolete:      if (o_loc_debug) call output
*obsolete:     &      (dbl_mb(ma_hk),1,ncart,1,nsph,ncart,nsph,1)
*obsolete:c
*obsolete:      if (do_terminate) then
*obsolete:        if (.not.spcart_terminate()) call errquit
*obsolete:     &      ('spcart_get_fulltranscs_ma:spcart_terminate failed ',911)
*obsolete:      vendif
*obsolete:      end
