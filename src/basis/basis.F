      logical function bas_version()
c $Id: basis.F,v 1.55 1997-03-25 01:42:26 d3e129 Exp $
c
c: Routine that calclulates the size of the common block structures 
c  used in the basis set object and the mapped representation object.
c:input none
c:output always true.
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
c
      integer cdata,idata,rdata
      integer mapidata, total4, total8
c
c.. character data 
      cdata = 256*2*nbasis_bsmx + 256*nbasis_rtdb_mx
      cdata = cdata + 16*ntags_bsmx*nbasis_bsmx 
c      
c.. real data 
c      rdata = mxbs_exndcf*nbasis_bsmx + 1
      rdata =  1 
      rdata = 8*rdata
c
c.. integer data in basis set object common
      idata = 3*nbasis_bsmx
      idata = idata + ndbs_head*nbasis_bsmx
      idata = idata + ndbs_tags*ntags_bsmx*nbasis_bsmx
      idata = idata + ndbs_ucont*(nucont_bsmx+1)*nbasis_bsmx
      idata = idata + nbasis_bsmx*6 + 1   ! assume logical and integer same size
      idata = idata + nbasis_rtdb_mx
      idata = 4*idata
c
c.. integer data in the mapped object.  
      mapidata =            5*3*nbasis_bsmx
      mapidata = mapidata + 4*nbasis_bsmx
      mapidata = 4*mapidata
c
c.. total space 
      total4 = idata + mapidata
      total8 = 2*total4 + rdata + cdata 
      total4 = total4 + rdata + cdata 
c
      write(6,'(////1x,a,f5.2,a)')
     &      ' **** basis set version ',bsversion,' ****'
      write(6,'(1x,a,i20,a)')
     &      '   character data in-core ',cdata,' bytes'
      write(6,'(1x,a,i20,a)')
     &      '   real      data in-core ',rdata,' bytes'
      write(6,'(1x,a,i20,a)')
     &      '   integer*4 data in-core ',idata,' bytes'
      write(6,'(1x,a,i20,a)')
     &      'or integer*8 data in-core ',(2*idata),' bytes'
      write(6,'(1x,a,i20,a)')
     &      '   integer*4 mapping data in-core ',
     &      mapidata,' bytes'
      write(6,'(1x,a,i20,a/)')
     &      'or integer*8 mapping data in-core ',
     &      (2*mapidata),' bytes'
      write(6,*)' total(4)   = ',total4,' bytes'
      write(6,*)' total(8)   = ',total8,' bytes'
c
c.. convert to kilobytes
c
      cdata    = (cdata    + 999) / 1000
      rdata    = (rdata    + 999) / 1000
      idata    = (idata    + 999) / 1000
      mapidata = (mapidata + 999) / 1000
      total4   = (total4   + 999) / 1000
      total8   = (total8   + 999) / 1000
      write(6,'(///1x,a,f5.2,a)')
     &      ' **** basis set version ',bsversion,' ****'
      write(6,'(1x,a,i20,a)')
     &      '   character data in-core ',cdata,' Kbytes'
      write(6,'(1x,a,i20,a)')
     &      '   real      data in-core ',rdata,' Kbytes'
      write(6,'(1x,a,i20,a)')
     &      '   integer*4 data in-core ',idata,' Kbytes'
      write(6,'(1x,a,i20,a)')
     &      'or integer*8 data in-core ',(2*idata),' Kbytes'
      write(6,'(1x,a,i20,a)')
     &      '   integer*4 mapping data in-core ',
     &      mapidata,' Kbytes'
      write(6,'(1x,a,i20,a/)')
     &      'or integer*8 mapping data in-core ',
     &      (2*mapidata),' Kbytes'
      write(6,*)' total(4)   = ',total4,' Kbytes'
      write(6,*)' total(8)   = ',total8,' Kbytes'
      write(6,'(////)')
c
      bas_version = .true.
      end
*.....................................................................
      logical function bas_create(basis,name)
c
c creates a handle and marks it active in the in-core data structure
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
c::functions
c ifill from util
c dfill from util
c::passed
c
      integer basis             ! [output] returned handle
      character*(*)name         ! [input] name of basis set.  
c
*:debug:      integer iiii
c
      external basis_data  ! This for the T3D linker
c
#include "bas_exndcf_sfn.fh"
c
      do 00100 basis=1,nbasis_bsmx
        if (.not.bsactive(basis)) goto 01000
00100 continue
c
      write(6,*)' bas_create: no free basis handles for ',name
      bas_create = .false.
      return
c
01000 continue
c
c store some information in basis data structure 
c (NOTE: name discarded in LOAD operation)
c      
      bs_name(basis) = name
      len_bs_name(basis) = inp_strlen(name)
c
c Initialize basis info to be empty
c
      bs_trans(basis) = ' '
      call ifill(ndbs_head, 0, infbs_head(1,basis), 1)
      exndcf(H_exndcf ,basis) = 0 ! handle 
      exndcf(K_exndcf ,basis) = 0 ! index
      exndcf(SZ_exndcf,basis) = 0 ! allocated size
      call ifill(ndbs_tags*ntags_bsmx, 0,
     $      infbs_tags(1,1,basis), 1)
      call ifill(ndbs_ucont*nucont_bsmx, 0,
     $      infbs_cont(1,1,basis), 1)
      ecp_name_assoc(basis) = '                                     '
      ecp_handle_assoc(basis) = 0
c
c Initialize geo-basis info to empty
c
      call ifill(3, 0, ibs_cn2ucn(1,basis), 1)
      call ifill(3, 0, ibs_cn2ce (1,basis), 1)
      call ifill(3, 0, ibs_ce2uce(1,basis), 1)
      call ifill(3, 0, ibs_cn2bfr(1,basis), 1)
      call ifill(3, 0, ibs_ce2cnr(1,basis), 1)
      ncont_tot_gb(basis) = 0
      nprim_tot_gb(basis) = 0
      nbf_tot_gb(basis)   = 0
      ibs_geom(basis)     = 0
      bas_norm_id(basis)  = BasNorm_UN
c
c Mark basis as active and return info
c
      bsactive(basis) = .true.
      basis = basis - BASIS_HANDLE_OFFSET 
      bas_create = .true.
c
c debug print all basis sets that are active
c
*:debug:      do iiii = 1,nbasis_bsmx
*:debug:        if (bsactive(iiii)) then
*:debug:          write(6,*)'bas_create:',bs_name(iiii),' ',iiii
*:debug:        endif
*:debug:      enddo
c
      end
*.....................................................................
      logical function bas_destroy(basisin)
      implicit none
c::functions
      logical bas_get_ecp_handle, ecp_check_handle, bas_do_destroy
      external bas_get_ecp_handle, ecp_check_handle, bas_do_destroy
c::passed
      integer basisin
c::local
      logical ignore
      integer ecpid
c
      ignore = bas_get_ecp_handle(basisin,ecpid)
      if (ecp_check_handle(ecpid,'bas_destroy')) then
        bas_destroy = bas_do_destroy(ecpid)
      else
        bas_destroy = .true.
      endif
      bas_destroy = bas_destroy .and. bas_do_destroy(basisin)
      end
*.....................................................................
      logical function bas_do_destroy(basisin)
c
c destroys information about an active incore basis
c and the associated mapping arrays.
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "bas_exndcf_dec.fh"
c::functions used
c::bas
      logical bas_check_handle
      logical gbs_map_clear
      external gbs_map_clear
      external bas_check_handle
c::passed
      integer basisin ![input] basis set handle to be destroyed
c::local
      integer basis
      integer h_tmp
c
#include "bas_exndcf_sfn.fh"
c
      bas_do_destroy = bas_check_handle(basisin,'bas_do_destroy')
      if (.not. bas_do_destroy) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      if(.not.gbs_map_clear(basisin)) then
        write(6,*)' error clearing map '
        bas_do_destroy = .false.
        return
      endif

c
      angular_bs(basis)    = -565
      bas_norm_id(basis)   = -565
      bsactive(basis)      = .false.
      bas_spherical(basis) = .false.
c
      if (exndcf(SZ_exndcf,basis).gt.0) then
        h_tmp = exndcf(H_exndcf,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_do_destroy: error freeing heap',911)
      endif
      exndcf(H_exndcf ,basis) = 0
      exndcf(K_exndcf ,basis) = 0
      exndcf(SZ_exndcf,basis) = 0
      bas_do_destroy = .true.
c
      return
      end
*.....................................................................
      logical function bas_check_handle(basisin,msg)
c
c Checks to see if a basis set handle is valid
c
      implicit none
#include "basP.fh"
c:passed
      integer basisin           ! [input] handle
      character*(*) msg         ! [input] error message
c::local
      integer basis
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_check_handle = basis.gt.0 .and. basis.le.nbasis_bsmx
      if (bas_check_handle)
     &      bas_check_handle = bas_check_handle .and. bsactive(basis)
c
* user's responsibility to deal with status
*      if (.not. bas_check_handle) then
*        write(6,*)msg,': basis handle is invalid '
*        write(6,*)'basis_check_handle: lexical handle ',basis
*        write(6,*)'basis_check_handle:         handle ',basisin
*      endif
      return
      end
*.....................................................................
      logical function ecp_check_handle(ecpidin,msg)
c
c Checks to see if a basis set handle is valid
c
      implicit none
#include "basdeclsP.fh"
#include "basP.fh"
#include "ecp_decP.fh"
#include "inp.fh"
c:functions
      logical bas_check_handle
      external bas_check_handle
c:passed
      integer ecpidin           ! [input] handle
      character*(*) msg         ! [input] error message
c::local
      character*255 newmsg    
c
#include "ecp_sfnP.fh"
c
      newmsg(1:17) = 'ecp_check_handle:'
      newmsg(18:) = msg 
      ecp_check_handle = bas_check_handle(ecpidin,newmsg)
      ecp_check_handle = ecp_check_handle .and.
     &    Is_ECP_in(ecpidin)
      return
      end
*.....................................................................
      logical function bas_321g_load(rtdb)
      implicit none
c
c routine to load the rtdb with 321g basis sets for atoms 1-20
c H to Ca
#include "mafdecls.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
c
c function declarations
c
c::function
      logical bas_create
      logical bas_print
      logical bas_rtdb_do_store
      logical bas_destroy
      logical bas_add_ucnt_init
      logical bas_add_ucnt
      external bas_print
      external bas_create
      external bas_rtdb_do_store
      external bas_destroy
      external bas_add_ucnt_init
      external bas_add_ucnt
c:: passed
      integer rtdb              ! [input] valid rtdb handle
c:: local
      integer basis, usebas
*      integer itag, icont
      character*2 a_tag
      integer a_lval, a_ngen, a_nprim
      double precision a_exp(3), a_cof(3)
      logical a_status
c
#include "bas_exndcf_sfn.fh"
c
      if (.not.bas_create(basis,'321g:1-20')) then
        write(6,*)' basis set handle not created '
        bas_321g_load = .false.
        return
      endif
c
      usebas = basis + BASIS_HANDLE_OFFSET
c     bs_name(usebas) = '321g:1-20'
c
c.. only h and o for now
c
c
c.. hydrogen and oxygen
c
c
      call print_exndcf(' b4 ucnt init ')
      if (.not. bas_add_ucnt_init(basis)) call errquit
     &      ('bas_321g_load: bas_add_ucnt_init failed',911)
      call print_exndcf(' after ucnt init')
c 1s on H
      a_tag   = 'H '
      a_lval  = 0 
      a_ngen  = 1
      a_nprim = 2
      a_exp (1) = 5.44717800d00
      a_exp (2) = 0.82454700d00
      a_cof (1) = 0.15628500d00
      a_cof (2) = 0.90469100d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c 2s on H
      a_lval  = 0 
      a_ngen  = 1
      a_nprim = 1
      a_exp (1) = 0.18319200d00
      a_cof (1) = 1.00000000d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c 1s on He
      a_tag   = 'He'
      a_lval  = 0 
      a_ngen  = 1
      a_nprim = 2
      a_exp (1) = 13.62670000d00
      a_exp (2) =  1.99935000d00
      a_cof (1) =  0.17523000d00
      a_cof (2) =  0.89348300d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c 2s on He
      a_lval  = 0 
      a_ngen  = 1
      a_nprim = 1
      a_exp (1) = 0.38299300d00
      a_cof (1) = 1.00000000d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c 1s on O
      a_tag   = 'O '
      a_lval  = 0 
      a_ngen  = 1
      a_nprim = 3
      a_exp (1) =  322.03700000d00 
      a_exp (2) =   48.43080000d00 
      a_exp (3) =   10.42060000d00 
      a_cof (1) =    0.05923940d00  
      a_cof (2) =    0.35150000d00  
      a_cof (3) =    0.70765800d00  
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c 1sp on O
c add s first
      a_tag   = 'O '
      a_lval  = 0
      a_ngen  = 1
      a_nprim = 2
      a_exp (1) =    7.40294000d00
      a_exp (2) =    1.57620000d00
      a_cof (1) =   -0.40445300d00
      a_cof (2) =    1.22156000d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c add p now
      a_tag   = 'O '
      a_lval  = 1
      a_ngen  = 1
      a_nprim = 2
      a_exp (1) =   7.40294000d00
      a_exp (2) =   1.57620000d00
      a_cof (1) =   0.24458600d00
      a_cof (2) =   0.85395500d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c 2sp on O
c add s first
      a_tag   = 'O '
      a_lval  = 0
      a_ngen  = 1
      a_nprim = 1
      a_exp (1) =  0.37368400d00
      a_cof (1) =  1.00000000d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
c add p now
      a_tag   = 'O '
      a_lval  = 1
      a_ngen  = 1
      a_nprim = 1
      a_exp (1) =  0.37368400d00
      a_cof (1) =  1.00000000d00
      a_status = bas_add_ucnt (basis,
     &      a_tag, a_lval, a_ngen, a_nprim, a_exp, a_cof, a_nprim)
      if (.not. a_status) call errquit
     &      ('bas_321g_load: bas_add_ucnt failed',911)
*_old_checking_code:      bs_tags(1,usebas) = 'H'
*_old_checking_code:      bs_tags(2,usebas) = 'He'
*_old_checking_code:      bs_tags(3,usebas) = 'O'
*_old_checking_code:      infbs_head(HEAD_NTAGS,usebas) = 3
*_old_checking_code:      infbs_head(HEAD_NCONT,usebas) = 9
*_old_checking_code:      infbs_head(HEAD_NPRIM,usebas) = 30
*_old_checking_code:      infbs_head(HEAD_NCOEF,usebas) = 30
*_old_checking_code:c.. hydrogen
*_old_checking_code:      itag = 1
*_old_checking_code:      infbs_tags(TAG_NCONT,itag,usebas) = 2
*_old_checking_code:      infbs_tags(TAG_NPRIM,itag,usebas) = 3
*_old_checking_code:      infbs_tags(TAG_NCOEF,itag,usebas) = 3
*_old_checking_code:      infbs_tags(TAG_FCONT,itag,usebas) = 1
*_old_checking_code:      infbs_tags(TAG_LCONT,itag,usebas) = 2
*_old_checking_code:c.. hydrogen cont 1 bas 1
*_old_checking_code:      icont = 1
*_old_checking_code:      infbs_cont(CONT_TYPE ,icont,usebas)  = 0
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas)  = 2
*_old_checking_code:      infbs_cont(CONT_NGEN ,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP ,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_ICFP ,icont,usebas)  = 3
*_old_checking_code:      infbs_cont(CONT_TAG  ,icont,usebas)  = itag
*_old_checking_code:      exndcf(1,usebas) =  5.44717800d00
*_old_checking_code:      exndcf(2,usebas) =  0.82454700d00
*_old_checking_code:      exndcf(3,usebas) =  0.15628500d00
*_old_checking_code:      exndcf(4,usebas) =  0.90469100d00
*_old_checking_code:c.. hydrogen cont 2 bas 2
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_cont(CONT_TYPE ,icont,usebas)  = 0
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_NGEN ,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP ,icont,usebas)  = 5
*_old_checking_code:      infbs_cont(CONT_ICFP ,icont,usebas)  = 6
*_old_checking_code:      infbs_cont(CONT_TAG  ,icont,usebas)  = itag
*_old_checking_code:      exndcf(5,usebas) =  0.18319200d00
*_old_checking_code:      exndcf(6,usebas) =  1.00000000d00
*_old_checking_code:c.. helium 
*_old_checking_code:      itag = itag + 1
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_tags(TAG_NCONT,itag,usebas) = 2
*_old_checking_code:      infbs_tags(TAG_NPRIM,itag,usebas) = 3
*_old_checking_code:      infbs_tags(TAG_NCOEF,itag,usebas) = 3
*_old_checking_code:      infbs_tags(TAG_FCONT,itag,usebas) = icont
*_old_checking_code:      infbs_tags(TAG_LCONT,itag,usebas) = icont + 2 - 1
*_old_checking_code:c.. helium cont 1 bas 3
*_old_checking_code:      infbs_cont(CONT_TYPE,icont,usebas)  = 0
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas) = 2
*_old_checking_code:      infbs_cont(CONT_NGEN,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP,icont,usebas)  = 7
*_old_checking_code:      infbs_cont(CONT_ICFP,icont,usebas)  = 9
*_old_checking_code:      infbs_cont(CONT_TAG,icont,usebas)   = itag
*_old_checking_code:      exndcf(7,usebas) = 13.62670000
*_old_checking_code:      exndcf(8,usebas) =  1.99935000
*_old_checking_code:      exndcf(9,usebas) =  0.17523000
*_old_checking_code:      exndcf(10,usebas)=  0.89348300
*_old_checking_code:c.. helium cont 2 bas 4
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_cont(CONT_TYPE,icont,usebas)  = 0
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas) = 1
*_old_checking_code:      infbs_cont(CONT_NGEN,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP,icont,usebas)  = 11
*_old_checking_code:      infbs_cont(CONT_ICFP,icont,usebas)  = 12
*_old_checking_code:      infbs_cont(CONT_TAG,icont,usebas)   = itag
*_old_checking_code:      exndcf(11,usebas) =  0.38299300
*_old_checking_code:      exndcf(12,usebas) =  1.00000000
*_old_checking_code:c.. oxygen
*_old_checking_code:      itag = itag + 1
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_tags(TAG_NCONT,itag,usebas) = 5
*_old_checking_code:      infbs_tags(TAG_NPRIM,itag,usebas) = 9
*_old_checking_code:      infbs_tags(TAG_NCOEF,itag,usebas) = 9
*_old_checking_code:      infbs_tags(TAG_FCONT,itag,usebas) = icont
*_old_checking_code:      infbs_tags(TAG_LCONT,itag,usebas) = icont + 5 - 1
*_old_checking_code:c.. oxygen cont 1 bas 5
*_old_checking_code:      infbs_cont(CONT_TYPE,icont,usebas)  = 0
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas) = 3
*_old_checking_code:      infbs_cont(CONT_NGEN,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP,icont,usebas)  = 13
*_old_checking_code:      infbs_cont(CONT_ICFP,icont,usebas)  = 16
*_old_checking_code:      infbs_cont(CONT_TAG,icont,usebas)   = itag
*_old_checking_code:      exndcf(13,usebas)  =  322.03700000d00
*_old_checking_code:      exndcf(14,usebas)  =   48.43080000d00
*_old_checking_code:      exndcf(15,usebas)  =   10.42060000d00
*_old_checking_code:      exndcf(16,usebas) =    0.05923940d00
*_old_checking_code:      exndcf(17,usebas) =    0.35150000d00
*_old_checking_code:      exndcf(18,usebas) =    0.70765800d00  
*_old_checking_code:c.. oxygen cont 2 bas 6
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_cont(CONT_TYPE,icont,usebas)  = 0
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas) = 2
*_old_checking_code:      infbs_cont(CONT_NGEN,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP,icont,usebas)  = 19
*_old_checking_code:      infbs_cont(CONT_ICFP,icont,usebas)  = 21
*_old_checking_code:      infbs_cont(CONT_TAG,icont,usebas)   = 2
*_old_checking_code:      exndcf(19,usebas) =   7.40294000d00
*_old_checking_code:      exndcf(20,usebas) =   1.57620000d00
*_old_checking_code:      exndcf(21,usebas) =  -0.40445300d00
*_old_checking_code:      exndcf(22,usebas) =   1.22156000d00
*_old_checking_code:c.. oxygen cont 3 bas 7      
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_cont(CONT_TYPE,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas) = 2
*_old_checking_code:      infbs_cont(CONT_NGEN,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP,icont,usebas)  = 23
*_old_checking_code:      infbs_cont(CONT_ICFP,icont,usebas)  = 25
*_old_checking_code:      infbs_cont(CONT_TAG,icont,usebas)   = 2
*_old_checking_code:      exndcf(23,usebas) =   7.40294000d00
*_old_checking_code:      exndcf(24,usebas) =   1.57620000d00
*_old_checking_code:      exndcf(25,usebas) =   0.24458600d00
*_old_checking_code:      exndcf(26,usebas) =   0.85395500d00
*_old_checking_code:c.. oxygen cont 4 bas 8
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_cont(CONT_TYPE,icont,usebas)  = 0
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas) = 1
*_old_checking_code:      infbs_cont(CONT_NGEN,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP,icont,usebas)  = 27
*_old_checking_code:      infbs_cont(CONT_ICFP,icont,usebas)  = 28
*_old_checking_code:      infbs_cont(CONT_TAG,icont,usebas)   = 2
*_old_checking_code:      exndcf(27,usebas) = 0.37368400d00
*_old_checking_code:      exndcf(28,usebas) = 1.00000000d00
*_old_checking_code:c.. oxygen cont 5 bas 9
*_old_checking_code:      icont = icont + 1
*_old_checking_code:      infbs_cont(CONT_TYPE,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_NPRIM,icont,usebas) = 1
*_old_checking_code:      infbs_cont(CONT_NGEN,icont,usebas)  = 1
*_old_checking_code:      infbs_cont(CONT_IEXP,icont,usebas)  = 29
*_old_checking_code:      infbs_cont(CONT_ICFP,icont,usebas)  = 30
*_old_checking_code:      infbs_cont(CONT_TAG,icont,usebas)   = 2
*_old_checking_code:      exndcf(29,usebas) = 0.37368400d00
*_old_checking_code:      exndcf(30,usebas) = 1.00000000d00
*_old_checking_code:c
      bas_321g_load = bas_rtdb_do_store(rtdb,
     &       bs_name(usebas), bs_tags(1,usebas), 
     &       infbs_head(1,usebas), infbs_tags(1,1,usebas),
     &       infbs_cont(1,1,usebas),
     &      dbl_mb(mb_exndcf(1,usebas)),3,9,30,' ')
c
      write(6,*)' inside 321g load '
      write(6,*)' inside 321g load '
      bas_321g_load = bas_print(basis)
      write(6,*)' inside 321g load '
      write(6,*)' inside 321g load '
c
      if (.not.bas_destroy(basis)) then
        write(6,*)' error releasing temproary basis handle '
        bas_321g_load = .false.
        return
      endif
c      
      bas_321g_load = .true.
      return
      end
*.....................................................................
      logical function bas_print(basisin)
c
c routine to print unique basis information that is in core
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
#include "geom.fh"
#include "ecp_decP.fh"
c
c function declarations
c      
      logical bas_check_handle, ecp_print
      external bas_check_handle, ecp_print
c:: passed
      integer basisin  ! [input] basis set handle
c:: local
      integer mytags, myucont, myprim, mycoef, basis
      integer i,j,k,l, ifcont, mygen, mytype, iexptr, icfptr
      integer atn, len_tag, len_ele
      character*2 symbol
      character*16 element
      character*3 ctype(0:6),cltype(2)
      character*3 shell_type
*. . . . . . . . . . . ! Room for tag+space+(+element+) = 16+1+1+16+1
      character*35 buffer  
      character*12 polynomial
c
#include "bas_exndcf.fh"
#include "ecp_sfnP.fh"
c
      if (Is_ECP_in(basisin)) then
        bas_print = ecp_print(basisin)
        return
      endif
c
      ctype(0)='S'
      ctype(1)='P'
      ctype(2)='D'
      ctype(3)='F'
      ctype(4)='G'
      ctype(5)='H'
      ctype(6)='I'
      cltype(1)='SP'
      cltype(2)='SPD'
      bas_print = .true.
      basis = basisin + BASIS_HANDLE_OFFSET
c
      bas_print = bas_check_handle(basisin,'bas_print')
      if (.not. bas_print) return
c
c print basis set information
c      
      if (bas_spherical(basis)) then
         polynomial = ' (spherical)'
      else
         polynomial = ' (cartesian)'
      endif
      write(6,1)bs_name(basis)(1:inp_strlen(bs_name(basis))), 
     $     bs_trans(basis)(1:inp_strlen(bs_trans(basis))), polynomial
 1    format('                      Basis "',a,'" -> "',a,'"',a/
     $       '                      -----')
      mytags  = infbs_head(HEAD_NTAGS,basis)
      if (mytags.le.0) then
        write(6,*)' Empty basis set '
        write(6,*)
        return
      endif
c
      myucont = infbs_head(HEAD_NCONT,basis)
      myprim  = infbs_head(HEAD_NPRIM,basis)
      mycoef  = infbs_head(HEAD_NCOEF,basis)
c
      do 00100 i=1,mytags

         if (geom_tag_to_element(bs_tags(i,basis), symbol, element,
     $        atn)) then
            len_tag = inp_strlen(bs_tags(i,basis))
            len_ele = inp_strlen(element)
            write(buffer,'(a,'' ('',a,'')'')')
     $           bs_tags(i,basis)(1:len_tag), element(1:len_ele)
         else
            buffer = bs_tags(i,basis)
         endif
         len_tag = inp_strlen(buffer)
         call util_print_centered(6, buffer, len_tag/2 + 1, .true.)

         myucont = infbs_tags(TAG_NCONT,i,basis)
c        
        ifcont = infbs_tags(TAG_FCONT,i,basis)
c     
        write(6,6)
 6      format(
     $       '            Exponent  Coefficients '/
     $       '       -------------- ',57('-'))
        do 00200 j=1,myucont
          myprim = infbs_cont(CONT_NPRIM,ifcont,basis)
          mygen  = infbs_cont(CONT_NGEN,ifcont,basis)
          
          mytype = infbs_cont(CONT_TYPE, ifcont, basis)
          if (mytype.lt.0) then
            shell_type = cltype(abs(mytype))
          else
            shell_type = ctype(mytype)
          endif
          iexptr = infbs_cont(CONT_IEXP,ifcont,basis) - 1
          icfptr = infbs_cont(CONT_ICFP,ifcont,basis) - 1
          do 00300 k=1,myprim
            write(6,7) j, shell_type(1:2),
     &          sf_exndcf((iexptr+k),basis),
     &          (sf_exndcf((icfptr+k+(l-1)*myprim),basis),l=1,mygen)
 7           format(1x,i2,1x,a2,1x,f14.6,20f10.6)
00300     continue
          write(6,*)
          ifcont = ifcont + 1
00200   continue
00100 continue
*debug:      write(6,*)' associated ecp ',ecp_name_assoc(basis)
c
c  If geom is set print out the info about total basis info 
c  associated with the geometry also
c
c  ... not done yet
c
      return
      end
*.....................................................................
      logical function bas_rtdb_load(rtdb, geom, basisin, name)
      implicit none
#include "inp.fh"
c::functions
      logical  bas_rtdb_do_load, bas_get_ecp_name, bas_set_ecp_name
      logical  bas_create, bas_do_destroy, bas_set_ecp_handle
      logical  ecp_set_parent_handle
      logical  bas_name_exist_rtdb
      external bas_rtdb_do_load, bas_get_ecp_name, bas_set_ecp_name
      external bas_create, bas_do_destroy, bas_set_ecp_handle
      external ecp_set_parent_handle
      external bas_name_exist_rtdb
c::passed
      integer rtdb
      integer geom
      integer basisin
      character*(*) name
c
c::local
      character*256 ecp_name
      integer ecpid
      logical status
      logical status_bas_cr
      logical status_ecp_load
      logical status_ecp_ph
      logical status_ecp
c
      bas_rtdb_load =
     &    bas_rtdb_do_load(rtdb, geom, basisin, name)
      if (.not.inp_compare(.false.,'ao basis',name(1:8)))
     &    return
      if (.not.bas_get_ecp_name(basisin,ecp_name)) call errquit
     &    ('bas_rtdb_load: bas_get_ecp_name failed',911)
      if (ecp_name .eq. '  ') then
        ecp_name = 'ecp basis'
        if (.not.bas_set_ecp_name(basisin,ecp_name)) call errquit
     &      ('bas_rtdb_load: bas_set_ecp_name failed',911)
      endif
      status = bas_name_exist_rtdb(rtdb,ecp_name)
      if (status) then
*
* ecp_name is on rtdb so load it
*
        status_bas_cr = bas_create(ecpid,ecp_name)
        if (.not.status_bas_cr) then
          call errquit
     &        ('bas_rtdb_load: bas_create failed for ecpid',911)
        endif

        status_ecp_load = bas_rtdb_do_load(rtdb,geom,ecpid,ecp_name)
        if (.not.status_ecp_load) then
          call errquit
     &        ('bas_rtdb_load: bas_rtdb_do_load failed for ecpid',911)
        endif

        status_ecp = bas_set_ecp_handle(basisin,ecpid)
        if (.not.status_ecp) then
          call errquit
     &        ('bas_rtdb_load: bas_set_ecp_handle failed',911)
        endif
        
        status_ecp_ph = ecp_set_parent_handle(ecpid,basisin)
        if (.not.status_ecp_ph) then
          call errquit
     &        ('bas_rtdb_load: ecp_set_parent_handle failed',911)
        endif
      endif
      end
*.....................................................................
      logical function bas_rtdb_do_load(rtdb, geom, basisin, name)
      implicit none
c
c routine that loads a basis set from the rtdb and using the 
c geometry information builds the mapping arrays to contractions/
c shells, basis functions, and centers.
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "global.fh"
c
c function declarations
c
c::function
c:bas
      logical bas_check_handle
      logical bas_geobas_build
      external bas_check_handle
      external bas_geobas_build
c:: passed  
      integer rtdb        ! [input] rtdb handle      
      integer geom        ! [input] geometry handle with info loaded
      integer basisin     ! [input] basis handle
      character*(*) name  ! [input] basis set name that must be on 
*. . . . . . . . . . . . .          the rtdb
c:: local
      character*25 nameex
      integer lentmp, basis, nexcf
      character*256 tmp
      logical rtdb_status
      integer h_tmp, k_tmp
c:: statement functions
#include "bas_exndcf.fh"
c:: initalize local
c
      rtdb_status = .true.
c
c
c check geom and basis handles returns false if either is invalid
c
      bas_rtdb_do_load = geom_check_handle(geom,'bas_rtdb_do_load')
      if (.not.bas_rtdb_do_load) return
      bas_rtdb_do_load = bas_check_handle(basisin,'bas_rtdb_do_load')
      if (.not.bas_rtdb_do_load) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c store geom tag with basis map info
c
      ibs_geom(basis) = geom
c
c translate "name" to current "context"
c
      bs_name(basis) = name
      len_bs_name(basis) = inp_strlen(name)
      if (.not.context_rtdb_match(rtdb,name,bs_trans(basis)))
     &       bs_trans(basis) = name
      len_bs_trans(basis) = inp_strlen(bs_trans(basis))
c
c generate rtdb names and load information
c
      tmp = 'basis:'//bs_trans(basis)(1:len_bs_trans(basis))
      lentmp = inp_strlen(tmp) + 1
c
      tmp(lentmp:) = ':bs_tags'
      rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, ntags_bsmx, bs_tags(1,basis))
c
      tmp(lentmp:) = ':assoc ecp name'
      rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, 1, ecp_name_assoc(basis))
c
      nexcf = 0
      tmp(lentmp:) = ':number of exps and coeffs'
      rtdb_status = rtdb_status .and.
     &    rtdb_get(rtdb,tmp,mt_int,1,nexcf)
c
      write(nameex,'(a23,i2)')' basis exps and coeffs ',basis
c
      if (exndcf(H_exndcf,basis) .ne. 0) then
         if (.not. ma_free_heap(exndcf(H_exndcf,basis)))
     $        call errquit('bas_rtdb_do_load: ma is corrupted',
     $        exndcf(H_exndcf,basis))
      endif

      if (.not.ma_alloc_get(mt_dbl,(nexcf+2),nameex,
     &    h_tmp, k_tmp)) then
        write(6,*)' not enough memory'
        call errquit
     &      (' bas_rtdb_do_load: error allocating space'//
     &      ' for exndcf',911)
      else
        exndcf(H_exndcf,basis) = h_tmp
        exndcf(K_exndcf,basis) = k_tmp
        exndcf(SZ_exndcf,basis) = (nexcf+2)
      endif
c
      tmp(lentmp:) = ':exps and coeffs'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_dbl, nexcf, dbl_mb(mb_exndcf(1,basis)))
c
      tmp(lentmp:) = ':header'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int, ndbs_head, infbs_head(1,basis))
c
      tmp(lentmp:) = ':tags info'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int,
     &       ndbs_tags*ntags_bsmx, infbs_tags(1,1,basis))
c
      tmp(lentmp:) = ':contraction info'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int,
     &       ndbs_ucont*nucont_bsmx, infbs_cont(1,1,basis))
c
c read the basis now get check status of read operations
c
      if (.not.rtdb_status) then
c
c     rjh ... can be quiet now since the application should
c     whine if it really did need the basis set
c
*         if(ga_nodeid().eq.0) then
*            write(6,*)
*            write(6,*) ' bas_rtdb_do_load: basis not present "',
*     $           bs_name(basis)(1:inp_strlen(bs_name(basis))),
*     &          '" -> "',
*     $           bs_trans(basis)(1:inp_strlen(bs_trans(basis))),
*     &          '"'
*            write(6,*)
*         endif
         bas_rtdb_do_load = .false.
c.....add diagnostics later
         return
      endif
c
c compute internal information and geobas maps
c
      bas_rtdb_do_load = bas_geobas_build(basisin)
c
c Add zero exponent S function for sp code, incore only
c
      k_tmp = infbs_head(HEAD_EXCFPTR,basis)
      infbs_cont(CONT_TYPE, 0,basis) = 0
      infbs_cont(CONT_NPRIM,0,basis) = 1
      infbs_cont(CONT_NGEN, 0,basis) = 1
      infbs_cont(CONT_TAG,  0,basis) = -1
      k_tmp = k_tmp + 1
      infbs_cont(CONT_IEXP, 0,basis) = k_tmp
      dbl_mb(mb_exndcf(k_tmp,basis)) = 0.0d00
      k_tmp = k_tmp + 1
      infbs_cont(CONT_ICFP, 0,basis) = k_tmp
      dbl_mb(mb_exndcf(k_tmp,basis)) = 1.0d00
      infbs_head(HEAD_EXCFPTR,basis) = k_tmp
      end
*.....................................................................
      logical function bas_geobas_build(basisin)
      implicit none
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "global.fh"
#include "bas_ibs_dec.fh"
#include "ecp_decP.fh"
#include "stdio.fh"
c::passed
      integer basisin
c::local
      integer basis
      integer geom
      integer nat
      integer i, idum_cont, idum_at
      integer j, jstart, jend, jsize
      integer kstart, kend, ksize, lsize, icount
      integer nbf, iu_cont, myang
      integer atn
      character*16 element
      character*2  symbol
      logical status
      logical foundit
      logical is_bq
      integer int_dummy, num_elec
      integer h_tmp, k_tmp
      character*2 tag12
      character*16 name_tmp
      double precision erep_save
      integer uce
*debug:mem      integer inode
c::functions
      integer nbf_from_ucont
      external nbf_from_ucont
      logical ecp_get_num_elec
      external ecp_get_num_elec
      logical bas_match_tags
      external bas_match_tags
c
#include "bas_ibs_sfn.fh"
#include "ecp_sfnP.fh"
c
      basis = basisin + BASIS_HANDLE_OFFSET
      geom  = ibs_geom(basis)
c
*debug:mem      do inode = 0,(ga_nnodes() - 1)
*debug:mem        if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem        call ga_sync()
*debug:mem      enddo
c

      status = geom_ncent(geom, nat)
      if (nat.eq.0.or..not.status) then
        write(6,*)' bas_geobas_build: ERROR '
        write(6,*)' number of centers is zero or weird'
        write(6,*)' nat = ',nat
        bas_geobas_build = .false.
c..... add diagnostics later
        return
      endif
c.... set spherical flag
      if (infbs_head(HEAD_SPH,basis).eq.1) then
        bas_spherical(basis) = .true.
      else
        bas_spherical(basis) = .false.
      endif
c
c... clear old ibs_ce2uce if it exists
      if (ibs_ce2uce(SZ_ibs,basis).gt.0) then
*debug:mem        write(6,*)' clearing old ce2uce data ',ga_nodeid()
*debug:mem        call util_flush(6)
        h_tmp = ibs_ce2uce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_ce2uce',911)
        ibs_ce2uce(H_ibs,basis)  = 0
        ibs_ce2uce(K_ibs,basis)  = 0
        ibs_ce2uce(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_ce2uce ',basis
      if (.not.ma_alloc_get(mt_int,nat,name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_ce2uce',911)
      else
*debug:mem        write(6,*)' generating ce2uce data ',ga_nodeid()
*debug:mem        call util_flush(6)
        ibs_ce2uce(H_ibs,basis)  = h_tmp
        ibs_ce2uce(K_ibs,basis)  = k_tmp
        ibs_ce2uce(SZ_ibs,basis) = nat
*debug:mem        call ifill(nat,0,int_mb(ibs_ce2uce(K_ibs,basis)),1)
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem          call ga_sync()
*debug:mem          call util_flush(6)
*debug:mem        enddo
      endif
c
c build center to unique center map
c
      if (nat.gt.nat_mx) then
        write(6,*)' nat     = ',nat
        write(6,*)' nat max = ',nat_mx
        call errquit ('bas_geobas_build: nat.gt.nat_mx',911)
      endif
      do 00100 i=1,nat
        foundit = .false.
*before match_tags:        do 00101 j = 1,infbs_head(HEAD_NTAGS,basis)
*before match_tags:          if(bas_match_tags(tags(i,geom),bs_tags(j,basis))) then
*before match_tags:            int_mb(mb_ibs_ce2uce(i,basis)) = j
*before match_tags:            foundit = .true.
*before match_tags:            goto 00102
*before match_tags:          endif
*before match_tags:00101   continue
        if (bas_match_tags(tags(i,geom),basisin,j)) then
          int_mb(mb_ibs_ce2uce(i,basis)) = j
          foundit = .true.
          goto 00102
        endif
        if (.not. foundit .and. .not. Is_ECP(basis)) then
          if (geom_tag_to_element(tags(i,geom), symbol, element,
     $        atn)) then
            if (ga_nodeid().eq.0)
     &          write(6,10) tags(i,geom)(1:inp_strlen(tags(i,geom))),
     $          element(1:inp_strlen(element)),
     $          bs_name(basis)(1:inp_strlen(bs_name(basis)))
 10         format(/' ERROR: geometry tag ',a,' (',a,
     &          ') is an atom ',
     $          'but has no functions in basis "',a,'"'/
     $          ' ERROR: only bq* centers can have no functions')
            call util_flush(6)
            call errquit
     &          ('bas_geobas_build: basis/geometry mismatch', 0)
          else
            tag12 = tags(i,geom)(1:2)
            is_bq = inp_compare(.false.,'bq',tag12)
            if (ga_nodeid().eq.0 .and.(.not.is_bq))
     &          write(6,11) i,
     &          tags(i,geom)(1:inp_strlen(tags(i,geom))),
     $          bs_name(basis)(1:inp_strlen(bs_name(basis)))
 11         format(/'WARNING: geometry tag ',i4, ' ', a,
     $          ' not found in basis "',a,'"'/)
            int_mb(mb_ibs_ce2uce(i,basis)) = 0
            call util_flush(6)
          endif
        endif
00102   continue
00100 continue
c
c build total # of contractions 
c      
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) then
*debug:mem            call MA_summarize_allocated_blocks()
*debug:mem            if (MA_verify_allocator_stuff()) then
*debug:mem              write(6,*)' no errors'
*debug:mem            else
*debug:mem              write(6,*)' errors'
*debug:mem            endif
*debug:mem          endif
*debug:mem          call ga_sync()
*debug:mem          call util_flush(6)
*debug:mem        enddo
*debug:mem      call bas_print_allocated_info('bas_geobas_build 1')
      ncont_tot_gb(basis)  = 0
      do 10200 i=1,nat
         uce = sf_ibs_ce2uce(i,basis)
*debug:mem         write(6,*)' myuce = ',uce,i,ga_nodeid()
*debug:mem         call util_flush(6)
         if (uce.gt.0) then
            idum_cont = infbs_tags(TAG_NCONT,uce,basis)
            ncont_tot_gb(basis)   = idum_cont + ncont_tot_gb(basis)
         endif
10200 continue
c
c allocate space for center -> contraction range map
c
c... clear old ibs_ce2cnr if it exists
      if (ibs_ce2cnr(SZ_ibs,basis).gt.0) then
*debug:mem        write(6,*)' clearing old ce2cnr data ',ga_nodeid()
*debug:mem        call util_flush(6)
        h_tmp = ibs_ce2cnr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_ce2cnr',911)
        ibs_ce2cnr(H_ibs,basis)  = 0
        ibs_ce2cnr(K_ibs,basis)  = 0
        ibs_ce2cnr(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_ce2cnr ',basis
      if (.not.ma_alloc_get(mt_int,(2*nat),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_ce2cnr',911)
      else
*debug:mem        write(6,*)' generating ce2cnr data ',ga_nodeid()
*debug:mem        call util_flush(6)
        ibs_ce2cnr(H_ibs,basis)  = h_tmp
        ibs_ce2cnr(K_ibs,basis)  = k_tmp
        ibs_ce2cnr(SZ_ibs,basis) = 2*nat
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem          call ga_sync()
*debug:mem          call util_flush(6)
*debug:mem        enddo
      endif
c
c build center -> contraction range map 
c
      int_dummy = 0
      do 00200 i=1,nat
         if (sf_ibs_ce2uce(i,basis).gt.0) then
            idum_cont =
     &           infbs_tags(TAG_NCONT,sf_ibs_ce2uce(i,basis),basis)
            int_mb(mb_ibs_ce2cnr(1,i,basis)) = int_dummy + 1
            int_mb(mb_ibs_ce2cnr(2,i,basis)) = int_dummy + idum_cont
            int_dummy = idum_cont + int_dummy
         else
*. . . . . . . . . . . . . . . . . . . . ! No functions on this center
            int_mb(mb_ibs_ce2cnr(1,i,basis)) = 0 
            int_mb(mb_ibs_ce2cnr(2,i,basis)) = -1
         endif
00200 continue
c
      if (ncont_tot_gb(basis) .eq. 0) call errquit
     $    ('bas_geobas_build: no functions in basis set', 0)
      if (ncont_tot_gb(basis) .gt. ncont_mx) then
        write(6,*)' number of contractions     = ',
     &      ncont_tot_gb(basis)
        write(6,*)' number of contractions max = ',ncont_mx
        call errquit ('bas_geobas_build: ncont.gt.ncont_mx ',911)
      endif
c
c allocate space for contraction -> center map
c
c... clear old ibs_cn2ce if it exists
      if (ibs_cn2ce(SZ_ibs,basis).gt.0) then
*debug:mem        write(6,*)' clearing old cn2ce data ',ga_nodeid()
*debug:mem        call util_flush(6)
        h_tmp = ibs_cn2ce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2ce',911)
        ibs_cn2ce(H_ibs,basis)  = 0
        ibs_cn2ce(K_ibs,basis)  = 0
        ibs_cn2ce(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2ce  ',basis
      if (.not.ma_alloc_get(mt_int,(1+ncont_tot_gb(basis)),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit('bas_geobas_build: error ma_alloc ibs_cn2ce',911)
      else
*debug:mem        write(6,*)' generating cn2ce data ',ga_nodeid()
*debug:mem        call util_flush(6)
        ibs_cn2ce(H_ibs,basis)  = h_tmp
        ibs_cn2ce(K_ibs,basis)  = k_tmp
        ibs_cn2ce(SZ_ibs,basis) = 1+ncont_tot_gb(basis)
      endif
c      
c build contraction -> center map
c
      do 00300 i=1,nat
         if (sf_ibs_ce2uce(i,basis).gt.0) then
            jstart = sf_ibs_ce2cnr(1,i,basis)
            jend   = sf_ibs_ce2cnr(2,i,basis)
            do 00400 j=jstart,jend
               int_mb(mb_ibs_cn2ce(j,basis)) = i
00400       continue
         endif
00300 continue
c
c set zero element of cn2ce to something useless
c
      int_mb(mb_ibs_cn2ce(0,basis)) = -1
c
c allocate space for ibs_cn2ucn map
c
c... clear old ibs_cn2ucn if it exists
      if (ibs_cn2ucn(SZ_ibs,basis).gt.0) then
*debug:mem        write(6,*)' clearing old cn2ucn data ',ga_nodeid()
*debug:mem        call util_flush(6)
        h_tmp = ibs_cn2ucn(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2ucn',911)
        ibs_cn2ucn(H_ibs,basis)  = 0
        ibs_cn2ucn(K_ibs,basis)  = 0
        ibs_cn2ucn(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2ucn ',basis
      if (.not.ma_alloc_get(mt_int,(1+ncont_tot_gb(basis)),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_cn2ucn',911)
      else
*debug:mem        write(6,*)' generating cn2ucn data ',ga_nodeid()
*debug:mem        call util_flush(6)
        ibs_cn2ucn(H_ibs,basis)  = h_tmp
        ibs_cn2ucn(K_ibs,basis)  = k_tmp
        ibs_cn2ucn(SZ_ibs,basis) = 1+ncont_tot_gb(basis)
      endif
c
c build contraction -> unique contraction map
c
      do 00500 i=1,nat
         jstart = sf_ibs_ce2cnr(1,i,basis)
         jend   = sf_ibs_ce2cnr(2,i,basis)
         jsize  = jend - jstart + 1
         if (jsize .gt. 0) then
            idum_at = sf_ibs_ce2uce(i,basis)
            kstart = infbs_tags(TAG_FCONT,idum_at,basis)
            kend   = infbs_tags(TAG_LCONT,idum_at,basis)
            ksize  = kend - kstart + 1
            lsize  = infbs_tags(TAG_NCONT,idum_at,basis)
            if (jsize.eq.ksize.and.ksize.eq.lsize) then
               icount = 0
               do 00600 j=jstart,jend
                  int_mb(mb_ibs_cn2ucn(j,basis)) = kstart + icount
                  icount = icount + 1
00600          continue
            else
               write(6,*)' bas_geobas_build: ERROR '
               write(6,*)' contraction range size mismatch'
               write(6,*)'        cont. range (',jstart,':',jend,')'
               write(6,*)' unique cont. range (',kstart,':',kend,')'
               write(6,*)'        cont. size: ',jsize
               write(6,*)' calculated unique cont. size: ',ksize
               write(6,*)'     lookup unique cont. size: ',lsize
               bas_geobas_build = .false.
               return
            endif
         endif
00500 continue
c
c set zero element
c
      int_mb(mb_ibs_cn2ucn(0,basis)) = 0
c
c allocate space for ibs_cn2bfr map
c
c... clear old ibs_cn2bfr if it exists
      if (ibs_cn2bfr(SZ_ibs,basis).gt.0) then
*debug:mem        write(6,*)' clearing old cn2bfr data ',ga_nodeid()
*debug:mem        call util_flush(6)
        h_tmp = ibs_cn2bfr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2bfr',911)
        ibs_cn2bfr(H_ibs,basis)  = 0
        ibs_cn2bfr(K_ibs,basis)  = 0
        ibs_cn2bfr(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2bfr ',basis
      if (.not.ma_alloc_get(mt_int,(2*(1+ncont_tot_gb(basis))),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_cn2bfr',911)
      else
*debug:mem        write(6,*)' generating cn2bfr data ',ga_nodeid()
*debug:mem        call util_flush(6)
        ibs_cn2bfr(H_ibs,basis)  = h_tmp
        ibs_cn2bfr(K_ibs,basis)  = k_tmp
        ibs_cn2bfr(SZ_ibs,basis) = 2*(1+ncont_tot_gb(basis))
      endif
c
c build nprim_tot_gb, nbf_tot_gb, and 
c contraction -> basis function range map 
c find nbfmax for basis (initialized in block data statement)
c
      nbf_tot_gb(basis)   = 0
      nprim_tot_gb(basis) = 0
      do 00700 i = 1,ncont_tot_gb(basis)
        iu_cont = sf_ibs_cn2ucn(i,basis)
c
        nbf = nbf_from_ucont(iu_cont,basisin)
        nbfmax_bs(basis) = max(nbfmax_bs(basis),nbf)
c
        int_mb(mb_ibs_cn2bfr(1,i,basis)) = nbf_tot_gb(basis) + 1
        int_mb(mb_ibs_cn2bfr(2,i,basis)) = nbf_tot_gb(basis) + nbf

        nbf_tot_gb(basis) = nbf_tot_gb(basis) + nbf
        nprim_tot_gb(basis) = nprim_tot_gb(basis) +
     &         infbs_cont(CONT_NPRIM,iu_cont,basis)
00700 continue
c
c set zero elements of cn2bfr
c
      int_mb(mb_ibs_cn2bfr(1,0,basis)) = 0
      int_mb(mb_ibs_cn2bfr(2,0,basis)) = 0
c
c build high angular momentum of this loaded <basis|geom> pair
c note angular_bs(*) initialized in block data function
c
      do 00800 i = 1, infbs_head(HEAD_NCONT,basis)
        myang   = abs(infbs_cont(CONT_TYPE,i,basis))
        angular_bs(basis) = max(angular_bs(basis),myang)
00800 continue
*
      if (Is_ECP(basis)) then
* if ECP basis must modify geometry data appropriately
        do i = 1,nat
          if (ecp_get_num_elec(basisin,tags(i,geom),num_elec)) then
            oecpcent(i,geom) = .true.
            charge(i,geom) = charge(i,geom) - num_elec
          endif
        enddo
        erep_save = erep(geom)
        call geom_compute_values(geom)
        write(luout,*)
     &      ' nuclear repulsion energy without ECP: ',erep_save
        write(luout,*)
     &      ' nuclear repulsion energy with    ECP: ',erep(geom)
      endif
*
      bas_geobas_build = .true.
*
      end
*.....................................................................
      logical function bas_rtdb_store(rtdb, name, basisin)
      implicit none
c
c routine that does an incore basis set store to the rtdb
c
#include "mafdecls.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "ecp_decP.fh"
c::passed
      integer rtdb              ! [input] handle to database
      character*(*) name        ! [input] name to use when storing db
      integer basisin           ! [input] handle to basis set
c::functions
      logical bas_check_handle, bas_rtdb_do_store
      external bas_check_handle, bas_rtdb_do_store
c
c     Store basis set (not geometry) related info about specified
c     basis in into the rtdb with the given name
c
c::local
      integer basis               ! Actual index into basis set arrays
      integer size_ex
c:: statement functions
#include "bas_exndcf.fh"      
#include "ecp_sfnP.fh"
c
      bas_rtdb_store = bas_check_handle(basisin,'bas_rtdb_store')
      if (.not. bas_rtdb_store) return
      basis = basisin + BASIS_HANDLE_OFFSET
c
      if (Is_ECP(basis)) then
        size_ex = (2*infbs_head(HEAD_NPRIM,basis))+
     &      infbs_head(HEAD_NCOEF,basis)
      else
        size_ex = infbs_head(HEAD_NPRIM,basis)+
     &      infbs_head(HEAD_NCOEF,basis)
      endif
      bas_rtdb_store = bas_rtdb_do_store(rtdb, name, 
     $    bs_tags(1,basis), infbs_head(1,basis),
     $    infbs_tags(1,1,basis),
     $    infbs_cont(1,1,basis),
     &    dbl_mb(mb_exndcf(1,basis)),
     $    infbs_head(HEAD_NTAGS,basis),
     $    infbs_head(HEAD_NCONT,basis), 
     $    size_ex, ecp_name_assoc(basis))
c
      end
*.....................................................................
      logical function bas_rtdb_do_store(rtdb, name, tagsin, 
     &       head_array, tags_array, ucont_array, excfin, ntagsin, 
     &       nucontin, nexcf, ecp_name)
      implicit none
c
c stores from argument data structures the basis set information to
c the rtdb. 
c
#include "mafdecls.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "inp.fh"
#include "rtdb.fh"
c
c    This routine stores the basis set information in the appropriate 
c    data structure on the run-time-data-base (rtdb).  
c
c    This is a private routine called by the user level routine 
c    bas_rtdb_store(rtdb, name, basis)
c
c::: functions
      logical bas_rtdb_in
      logical bas_rtdb_add
      external bas_rtdb_in
      external bas_rtdb_add
c::: passed
      integer rtdb ! [input] rtdb handle
      character*(*) name ! [input] name of basis set
      integer ntagsin ! [input] number of tags
      integer nucontin ! [input] number of unique contractions
      integer nexcf ! [input] number of exponents and 
      character*16 tagsin(ntagsin) ! [input] name of tags
      integer head_array(ndbs_head) ! [input] head data
      integer tags_array(ndbs_tags,ntagsin) ! [input] tag data
      integer ucont_array(ndbs_ucont,nucontin) ! [input] unique 
*. . . . . . . . . . . . . . . . . . . . . . .   contraction data
      character*(*) ecp_name                   ! [input] ecp 
*. . . . . . . . . . . . . . . . . . . . . . .   associated with 
*. . . . . . . . . . . . . . . . . . . . . . .   normal basis
      double precision excfin(nexcf) ! [input] exponents 
c. . . . . . . . . . . . . . . . . . .         contractions coeffs.
c::: local
      character*256 tmp
      integer len_name, lentmp
      logical status
c
      bas_rtdb_do_store = .true.
c
      status = bas_rtdb_in(rtdb)
c
c generate rtdb names and store information
c      
      len_name = inp_strlen(name)
      tmp = 'basis:'//name(1:len_name)
      lentmp = inp_strlen(tmp) + 1

      status = .true.
      status = status.and.bas_rtdb_add(rtdb,name)
      tmp(lentmp:) = ':bs_tags'
      status = status .and. rtdb_cput(rtdb,tmp,ntagsin,tagsin)
      
      tmp(lentmp:) = ':assoc ecp name'
      status = status .and. rtdb_cput(rtdb,tmp,1,ecp_name)
c
      tmp(lentmp:) = ':number of exps and coeffs'
      status = status .and. rtdb_put(rtdb,tmp,mt_int,1,nexcf)

      tmp(lentmp:) = ':exps and coeffs'
      status = status .and. rtdb_put(rtdb,tmp,mt_dbl,nexcf,excfin)

      tmp(lentmp:) = ':header'
      status = status .and.
     &       rtdb_put(rtdb,tmp,mt_int,ndbs_head,head_array)

      tmp(lentmp:) = ':tags info'
      status = status .and. rtdb_put(
     &       rtdb,tmp,mt_int,(ndbs_tags*ntagsin), tags_array)

      tmp(lentmp:) = ':contraction info'
      status = status .and. rtdb_put(
     &       rtdb,tmp,mt_int,(ndbs_ucont*nucontin),ucont_array)
c
c read the basis now get check status of read operations
c
      if (.not.status) then
        write(6,*)' bas_rtdb_store: ERROR '
        write(6,*)' one or more put operations failed '
        bas_rtdb_do_store = .false.
c..... add diagnostics later
        return
      endif
      return
      end
*.....................................................................
      logical function bas_high_angular(basisin,high_angular)
      implicit none
c
c  calculate, return and store high angular momentem function
c   for given basis. 
c
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin       ! [input] basis set handle
      integer high_angular  ! [output] high angular momentum of basis
c:local
      integer basis, myucont, i, myang
c
      bas_high_angular = bas_check_handle(basisin,'bas_high_angular')
      if (.not. bas_high_angular ) then
        write(6,*) 'bas_high_angular: basis handle not valid '
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
      if (angular_bs(basis) .gt. -565) then
        high_angular = angular_bs(basis)
        bas_high_angular = .true.
        return
      endif
c
      myucont = infbs_head(HEAD_NCONT,basis)
      high_angular = -565
c
      do 00100 i = 1,myucont
        myang = abs(infbs_cont(CONT_TYPE,i,basis))
        high_angular = max(high_angular, myang)
00100 continue
c
      angular_bs(basis) = high_angular
      bas_high_angular = .true.
      return
      end
*.....................................................................
      logical function gbs_map_clear(basisin)
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c
c routine to clear online map information and basis information
c
c::functions
      logical bas_check_handle
      external bas_check_handle
c:util
c     ifill      
c::passed
      integer basisin  ! [input] basis set handle
c::local
      integer basis
      integer h_tmp
c
#include "bas_ibs_sfn.fh"
c
      gbs_map_clear = bas_check_handle(basisin,'gbs_map_clear')
      if (.not. gbs_map_clear ) then
        write(6,*) ' basis handle not valid '
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
c... clear old ibs_ce2uce if it exists
      if (ibs_ce2uce(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_ce2uce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_ce2uce',911)
        ibs_ce2uce(H_ibs,basis)  = 0
        ibs_ce2uce(K_ibs,basis)  = 0
        ibs_ce2uce(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_ce2cnr if it exists
      if (ibs_ce2cnr(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_ce2cnr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_ce2cnr',911)
        ibs_ce2cnr(H_ibs,basis)  = 0
        ibs_ce2cnr(K_ibs,basis)  = 0
        ibs_ce2cnr(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_cn2ce if it exists
      if (ibs_cn2ce(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_cn2ce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_cn2ce',911)
        ibs_cn2ce(H_ibs,basis)  = 0
        ibs_cn2ce(K_ibs,basis)  = 0
        ibs_cn2ce(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_cn2ucn if it exists
      if (ibs_cn2ucn(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_cn2ucn(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_cn2ucn',911)
        ibs_cn2ucn(H_ibs,basis)  = 0
        ibs_cn2ucn(K_ibs,basis)  = 0
        ibs_cn2ucn(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_cn2bfr if it exists
      if (ibs_cn2bfr(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_cn2bfr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_cn2bfr',911)
        ibs_cn2bfr(H_ibs,basis)  = 0
        ibs_cn2bfr(K_ibs,basis)  = 0
        ibs_cn2bfr(SZ_ibs,basis) = 0 
      endif
c
      call ifill(3,0,ibs_cn2ucn(1,basis),1)
      call ifill(3,0,ibs_cn2ce (1,basis),1)
      call ifill(3,0,ibs_cn2bfr(1,basis),1)
      call ifill(3,0,ibs_ce2uce(1,basis),1)
      call ifill(3,0,ibs_ce2cnr(1,basis),1)
      ncont_tot_gb(basis) = 0
      nprim_tot_gb(basis) = 0
      nbf_tot_gb(basis)   = 0
c
      gbs_map_clear = .true.
      return
      end
*.....................................................................
      logical function gbs_map_print(basisin)
      implicit none
c
c prints the basis set <-> geometry mapping information
c
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin  ! [input] basis set handle
c::local
      integer nat, basis, i
      integer myfirst, mylast, mysize, mycenter, myucont
      integer mygeom
      logical status
c
#include "bas_ibs_sfn.fh"
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c check geom and basis handles returns false if either is invalid
c
      mygeom = ibs_geom(basis)
      gbs_map_print=geom_check_handle(mygeom,'gbs_map_print')
      if (.not.gbs_map_print) return
      gbs_map_print=bas_check_handle(basisin,'gbs_map_print')
      if (.not.gbs_map_print) return
c
c find number of atoms
      status = geom_ncent(mygeom, nat)
c
      if (nat.eq.0.or..not.status) then
        write(6,*)' gbs_map_print: ERROR '
        write(6,*)' number of centers is zero or weird'
        write(6,*)' nat = ',nat
        gbs_map_print = .false.
c..... add diagnostics later
        return
      endif
c
c print global information
c
      write(6,*)'<<< GBS_MAP_PRINT >>>'
      write(6,*)' total number of atoms           :',nat
      write(6,*)' total number of contractions    :',
     &       ncont_tot_gb(basis)
      write(6,*)' total number of primitives      :',
     &       nprim_tot_gb(basis)
      write(6,*)' total number of basis functions :',
     &       nbf_tot_gb(basis)
c
c print center based mapping information.  
c
      write(6,*)' '
      write(6,*)'==================================================',
     &       '=============================='
      write(6,*)' center -> unique center map          <ibs_ce2uce>'
      write(6,*)'        -> contraction range map      <ibs_ce2cnr>'
      write(6,*)'==================================================',
     &       '=============================='
      do 00100 i=1,nat
        write(6,'(1x,a,i4,2x,a,i3)')
     &         'center:',i,'maps to unique center:',
     &        sf_ibs_ce2uce(i,basis)
        myfirst = sf_ibs_ce2cnr(1,i,basis)
        mylast  = sf_ibs_ce2cnr(2,i,basis)
        mysize  = mylast - myfirst + 1
        write(6,'(14x,a,i4,2x,a,i4,a,i4,a,/)')
     &         'has',mysize,'contractions   <first:',
     &         myfirst,'>  <last:',mylast,'>'
00100 continue
c
c print contraction based mapping information
c
      write(6,*)' '
      write(6,*)'==================================================',
     &       '=============================='
      write(6,*)' contraction -> center map                     ',
     &       '                <ibs_cn2ce>'
      write(6,*)'             -> unique contraction in basis set',
     &       '                <ibs_cn2ucn>'
      write(6,*)'             -> basis function range           ',
     &       '                <ibs_cn2bfr>'
      write(6,*)'==================================================',
     &       '=============================='
c
      do 00200 i=1,ncont_tot_gb(basis)
        mycenter = sf_ibs_cn2ce(i,basis)
        myucont  = sf_ibs_cn2ucn(i,basis)
        myfirst  = sf_ibs_cn2bfr(1,i,basis)
        mylast   = sf_ibs_cn2bfr(2,i,basis)
        mysize   = mylast - myfirst + 1
        write(6,'(1x,a,i4,2x,a,i3)')
     &         'contraction',i,'is on center:',mycenter
        write(6,'(18x,a,i3)')
     &         'is represented by unique contraction:',myucont
        write(6,'(18x,a,i5,a,i5,a,i5,a,/)')
     &         'has',mysize,' basis functions <first:',myfirst,
     &         '>  <last:',mylast,'>'
00200 continue
c
      gbs_map_print = .true.
      return
      end
*.....................................................................
      logical function bas_rtdb_out(rtdb)
c     
c     output to rtdb info about known basis sets
c     
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "rtdb.fh"
c     
c::passed
      integer rtdb ! [input] run time data base handle
c
      bas_rtdb_out  =
     &     rtdb_put(rtdb, 'basis:nbasis',
     &       MT_INT, 1, nbasis_rtdb)
     &     .and.
     &     rtdb_cput(rtdb, 'basis:names', nbasis_rtdb,
     &       bs_names_rtdb)
      if (.not. bas_rtdb_out) 
     &     write(6,*) ' bas_rtdb_out: rtdb is corrupt '
c     
      end
*.....................................................................
      logical function bas_rtdb_add(rtdb, name)
      implicit none
c
c add basis set name to known basis set list on rtdb
c
#include "mafdecls.fh"
#include "basP.fh"
#include "rtdb.fh"
#include "inp.fh"
c     
      integer rtdb              ! [input] run time data base handle
      character*(*) name        ! [input] name of basis set to add
      integer basis
      logical status
      integer ln
      logical bas_rtdb_in, bas_rtdb_out
      external bas_rtdb_in, bas_rtdb_out
c     
c     See if name is on the rtdb already
c     
      ln = inp_strlen(name)
      status = bas_rtdb_in(rtdb)
      bas_rtdb_add = .true.
      do 00100 basis = 1, nbasis_rtdb
        if (name(1:ln) .eq.
     &         bs_names_rtdb(basis)(1:len_bs_rtdb(basis))) return
00100 continue
c     
c     Name is not present ... add and rewrite info
c     
      if (nbasis_rtdb .eq. nbasis_rtdb_mx) then
         write(6,*) ' bas_rtdb_add: too many basis tries on rtdb ',
     &      name
         bas_rtdb_add = .false.
         return
      endif
      nbasis_rtdb = nbasis_rtdb + 1
      bs_names_rtdb(nbasis_rtdb) = name
      len_bs_rtdb(nbasis_rtdb) = ln
c     
      bas_rtdb_add = bas_rtdb_out(rtdb)
      if (.not. bas_rtdb_add) then
         write(6,*) ' bas_rtdb_add: rtdb error adding ', name(1:ln)
         return
      endif
c     
      bas_rtdb_add = .true.
c     
      end
*.....................................................................
      logical function bas_print_all()
c
c routine to print active all basis set(s) information
c
      implicit none
#include "basP.fh"
c::function
      logical bas_print
      external bas_print
c::local
      integer basis,basin
c
      bas_print_all = .true.
      do 00100 basis=1,nbasis_bsmx
        if(bsactive(basis)) then
          basin = basis - BASIS_HANDLE_OFFSET
          bas_print_all = bas_print_all .and. bas_print(basin)
        endif
00100 continue
c
      return
      end
*.....................................................................
      logical function bas_rtdb_in(rtdb)
c     
c     load in info about known basis sets ... this is more
c     for diagnostic and debugging purposes
c     
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "rtdb.fh"
#include "inp.fh"
c     
      integer rtdb     ! [input] run time data base handle
      integer bas
      bas_rtdb_in = .false.
      nbasis_rtdb = 0
      if (rtdb_get(rtdb, 'basis:nbasis', MT_INT, 1, nbasis_rtdb))
     $     then
        if (.not. rtdb_cget(rtdb,'basis:names', nbasis_rtdb_mx,
     $        bs_names_rtdb)) then
          write(6,*) 'bas_rtdb_in: rtdb corrupt'
        else
          do 00100 bas = 1, nbasis_rtdb
            len_bs_rtdb(bas) = inp_strlen(bs_names_rtdb(bas))
00100     continue
          bas_rtdb_in = .true.
        endif
      endif
c     
      return
      end
*.....................................................................
      logical function bas_cn2ce(basisin,cont,center)
      implicit none
c
c returns the center for a given mapped contraction
c
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin   ! [input] basis set handle
      integer cont      ! [input] mapped contraction index
      integer center    ! [output] center index
c::local
      integer basis
#include "bas_ibs_sfn.fh"
c
      bas_cn2ce = bas_check_handle(basisin,'bas_cn2ce')
      if(.not.bas_cn2ce) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_cn2ce = cont.ge.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2ce) then
        write(6,*)' bas_cn2ce: invalid contraction information '
        write(6,*)' contraction range is 0:',ncont_tot_gb(basis)
        write(6,*)' input contraction was: ',cont
        return
      endif
      center = sf_ibs_cn2ce(cont,basis)
c
      return
      end
*.....................................................................
      logical function bas_cn2uce(basisin,cont,ucenter)
      implicit none
c
c returns the UNIQUE center for a given mapped contraction
c
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin   ! [input] basis set handle
      integer cont      ! [input] mapped contraction index
      integer ucenter    ! [output] center index
c::local
      integer basis
#include "bas_ibs_sfn.fh"
c
      bas_cn2uce = bas_check_handle(basisin,'bas_cn2uce')
      if(.not.bas_cn2uce) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_cn2uce = cont.ge.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2uce) then
        write(6,*)' bas_cn2uce: invalid contraction information '
        write(6,*)' contraction range is 0:',ncont_tot_gb(basis)
        write(6,*)' input contraction was: ',cont
        return
      endif
c
      ucenter = sf_ibs_cn2ce(cont,basis)
      ucenter = sf_ibs_ce2uce(ucenter,basis)
c
      end
*.....................................................................
      logical function bas_cn2bfr(basisin,cont,ifirst,ilast)
c
c returns the first basis function index of a mapped contraction 
c in ifirst and the last basis function index in ilast
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin   ! [input] basis set handle
      integer cont      ! [input] mapped contraction index
      integer ifirst    ! [output] first basis function
      integer ilast     ! [output] last basis function     
c::local 
      integer basis
c
#include "bas_ibs_sfn.fh"
c
      bas_cn2bfr = .true.
#ifdef BASIS_DEBUG      
      bas_cn2bfr = bas_check_handle(basisin,'bas_cn2bfr')
      if(.not.bas_cn2bfr) return
#endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_cn2bfr = cont.ge.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2bfr) then
        write(6,*)' bas_cn2bfr: invalid contraction information '
        write(6,*)' contraction range is 0:',ncont_tot_gb(basis)
        write(6,*)' input contraction was: ',cont
        return
      endif
c
      ifirst = sf_ibs_cn2bfr(1,cont,basis)
      ilast  = sf_ibs_cn2bfr(2,cont,basis)
c
      return
      end
*.....................................................................
      logical function bas_ce2bfr(basis, icent, ibflo, ibfhi)
c
c  returns the basis function range for a given center
c
      implicit none
      integer basis ! [input] handle
      integer icent ! [input] no. of center
      integer ibflo, ibfhi ! [output] range of functions on center
c
      integer cnlo, cnhi, tmp
      logical status
      logical bas_ce2cnr, bas_cn2bfr
      external bas_ce2cnr, bas_cn2bfr
c
      status = .true.
      status = status .and. bas_ce2cnr(basis, icent, cnlo, cnhi)
      if (cnhi .gt. 0) then
         status = status .and. bas_cn2bfr(basis, cnlo, ibflo, tmp)
         status = status .and. bas_cn2bfr(basis, cnhi, tmp, ibfhi)
      else
         ibflo = 0
         ibfhi = -1
      endif
c
      bas_ce2bfr = status
c
      end
*.....................................................................
      logical function bas_ce2cnr(basisin,center,ifirst,ilast)
c
c returns the mapped contraction range on a given center
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geom.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin  ! [input] basis set handle         
      integer center   ! [input] center index 
      integer ifirst   ! [output] first mapped contraction
      integer ilast    ! [output] last mapped contraction
c::local 
      integer basis, nat
#include "bas_ibs_sfn.fh"
c
      bas_ce2cnr = .true.
#ifdef BASIS_DEBUG
      bas_ce2cnr = bas_check_handle(basisin,'bas_ce2cnr')
      if(.not.bas_ce2cnr) return
#endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_ce2cnr = geom_ncent(ibs_geom(basis),nat)
      if (nat.eq.0.or..not.bas_ce2cnr) then
        bas_ce2cnr = .false.
        write(6,*)' bas_ce2cnr: ERROR '
        write(6,*)' number of centers is zero or weird'
        write(6,*)' nat = ',nat
c..... add diagnostics later
        return
      endif

      bas_ce2cnr = center.gt.0 .and. center.le.nat
      if (.not.bas_ce2cnr) then
        write(6,*)' bas_ce2cnr: invalid center information '
        write(6,*)' center range is 1:',nat
        write(6,*)' input center was : ',center
        return
      endif
c
      ifirst = sf_ibs_ce2cnr(1, center, basis)
      ilast  = sf_ibs_ce2cnr(2, center, basis)
c
      return
      end
*.....................................................................
      logical function bas_bf2ce(basisin,testbf,center)
c
c routine to return the center of a given basis function 
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geom.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin ! [input] basis set handle         
      integer testbf  ! [input] basis function index
      integer center  ! [output] center index
c::local
      integer basis, nat, iat, ibflo, ibfhi
c
#include "bas_ibs_sfn.fh"
c
      bas_bf2ce = bas_check_handle(basisin,'bas_bf2ce')
      if (.not. bas_bf2ce) return
c
      basis = basisin + BASIS_HANDLE_OFFSET 
      bas_bf2ce = geom_ncent(ibs_geom(basis),nat)
      if (.not.bas_bf2ce .or. nat.le.0) then
        bas_bf2ce = .false.
        write(6,*)' bas_bf2ce: ERROR '
        write(6,*)' number of centers is zero or weird'
        write(6,*)' nat = ',nat
        return
      endif
c
c... linear search through atoms
c
      center = -1
      do 00100 iat = 1,nat
        ibflo = sf_ibs_cn2bfr(1,(sf_ibs_ce2cnr(1,iat,basis)),basis)
        ibfhi = sf_ibs_cn2bfr(2,(sf_ibs_ce2cnr(2,iat,basis)),basis)
        if (testbf.ge.ibflo.and.testbf.le.ibfhi) then
          center = iat
          return
        endif
00100 continue
c
      end
*.....................................................................
      logical function bas_bf2cn(basisin,testbf,cont)
c
c returns the mapped contraction index that contains the given 
c basis function index
c 
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_numcont
      logical bas_check_handle
      external bas_numcont
      external bas_check_handle
c::passed
      integer basisin ! [input] basis set handle         
      integer testbf  ! [input] basis function index
      integer cont    ! [output] mapped contraction index
c::local
      integer basis, icont, ibflo, ibfhi
      integer numcont
c
#include "bas_ibs_sfn.fh"
c
      bas_bf2cn = bas_check_handle(basisin,'bas_bf2cn')
      if (.not. bas_bf2cn) return
c
      bas_bf2cn = bas_numcont(basisin,numcont)
      if (.not.bas_bf2cn) then
        write(6,*)'bas_bf2cn: could not get number of contractions'
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET 
c
c... linear search through contractions
      cont = -1 
      do 00100 icont = 1,numcont
        ibflo = sf_ibs_cn2bfr(1,icont,basis)
        ibfhi = sf_ibs_cn2bfr(2,icont,basis)
        if(testbf.ge.ibflo.and.testbf.le.ibfhi) then
          cont = icont
          return
        endif
00100 continue
c
      end
*.....................................................................
      logical function bas_numbf(basisin,nbf)
c
c returns the total number of basis functions of the mapped basis set.
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin ! [input] basis set handle         
      integer nbf     ! [output] number of basis functions
c::local
      integer basis 
c
      nbf = -6589
      bas_numbf = bas_check_handle(basisin,'bas_numbf')
      if (.not. bas_numbf) return

      basis = basisin + BASIS_HANDLE_OFFSET 
      nbf = nbf_tot_gb(basis)
      bas_numbf = .true.
      return
      end
*.....................................................................
      logical function bas_name(basisin,basis_name,trans_name)
c
c returns the name and translated name of the basis set
c
      implicit none
#include "basP.fh"      
#include "inp.fh"
c::functions
c inp_strlen() from inp
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer       basisin ! [input] basis set handle
      character*(*) basis_name ! [output] basis set name when loaded
      character*(*) trans_name ! [output] basis set name in context 
c::local
      integer basis   ! actual offset into basis arrays
*      integer lenofit ! length of name
c
      bas_name = bas_check_handle(basisin,'bas_name')
      if (.not. bas_name) return
c
      basis = basisin + BASIS_HANDLE_OFFSET 
c
*      lenofit = inp_strlen(bs_name(basis))
*      basis_name(1:lenofit) = bs_name(basis)(1:lenofit)
*      lenofit = inp_strlen(bs_trans(basis))
*      basis_name(1:lenofit) = bs_trans(basis)(1:lenofit)
      basis_name = bs_name(basis)
      trans_name = bs_trans(basis)
c
      end
*.....................................................................
      logical function bas_get_exponent(basisin,icont,exp)
c
c  returns the exponents of a given mapped contraction 
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin         ! [input] basis set handle
      integer icont           ! [input] mapped contraction index
      double precision exp(*) ! [output] exponents
c::local
      integer basis, myucont, icontmax
      integer myprim,myexptr
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      bas_get_exponent = bas_check_handle(basisin,'bas_get_exponent')
      if (.not.bas_get_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 
      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)
      bas_get_exponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_get_exponent)) then
        write(6,*)' bas_get_exponent: ERROR '
        write(6,*)' contraction range for basis is 0:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      call dcopy(myprim,dbl_mb(mb_exndcf(myexptr,basis)),1,exp,1)
c
      bas_get_exponent = .true.
c
      return
      end
*.....................................................................
      logical function bas_cont_tag(basisin,icont,tagout)
      implicit none
#include "basP.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c::-function
      logical bas_check_handle
      external bas_check_handle
c::-passed
      integer basisin
      integer icont
      character*(*) tagout
c::-local
      integer center, ucenter, basis
      integer len_tagout
#include "bas_ibs_sfn.fh"
c
      bas_cont_tag = bas_check_handle(basisin,'bas_cont_tag')
      if (.not.bas_cont_tag) return

      basis = basisin + BASIS_HANDLE_OFFSET 
      
      center = sf_ibs_cn2ce(icont,basis)
      ucenter = sf_ibs_ce2uce(center,basis)
      len_tagout = len(tagout)
      tagout = bs_tags(ucenter,basis)(1:len_tagout)
      end
*.....................................................................
      logical function bas_continfo(basisin,icont,
     &       type,nprimo,ngeno,sphcart)
c
c  returns the generic information about the given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin           ! [input] basis handle
      integer icont             ! [input] contraction index
      integer type              ! [output] type (sp/s/p/d/..)
      integer nprimo            ! [output] no. of primitives
      integer ngeno             ! [output] no. of contractions
      integer sphcart           ! [output] 0/1 for cartesian/shperical
c::local
      integer basis,myucont,icontmax
c
#include "bas_ibs_sfn.fh"
c
      nprimo = -123
      ngeno  = -456
      sphcart = -789
c
      bas_continfo = bas_check_handle(basisin,'bas_continfo')
      if (.not.bas_continfo) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = ncont_tot_gb(basis)
c
      if (.not.(icont.ge.0.and.icont.le.icontmax)) then
        write(6,*)' bas_continfo: ERROR '
        write(6,*)' contraction range for basis is 0:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        bas_continfo = .false.
        return
      endif
c
      myucont = sf_ibs_cn2ucn(icont,basis)
c... 
      if (bas_spherical(basis)) then
        sphcart = 1
      else
        sphcart = 0
      endif
      type    = infbs_cont(CONT_TYPE, myucont,basis)
      nprimo  = infbs_cont(CONT_NPRIM,myucont,basis)
      ngeno   = infbs_cont(CONT_NGEN, myucont,basis)
      bas_continfo=.true.
      return
      end
*.....................................................................
      logical function bas_numcont(basisin,numcont)
c
c returns the total number of mapped contractions/shells for the 
c given basis set
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin   ! [input] basis set handle
      integer numcont   ! [output] number of mapped contractions
c::local
      integer basis
c
      numcont = -6589
      bas_numcont = bas_check_handle(basisin,'bas_numcont')
      if (.not.bas_numcont) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      numcont = ncont_tot_gb(basis)

      bas_numcont = .true.
      return
      end
*.....................................................................
      logical function bas_get_coeff(basisin,icont,coeff)
c
c returns the coefficients for the given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin            ! [input] basis set handle
      integer icont              ! [input] mapped contraction index
      double precision coeff(*)  ! [output] mapped contraction coeffs.
c::local
      integer basis, myucont, icontmax
      integer mycoeffptr, myprim, mygen
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      bas_get_coeff = bas_check_handle(basisin,'bas_get_coeff')
      if (.not.bas_get_coeff) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)
      bas_get_coeff = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_get_coeff)) then
        write(6,*)' bas_get_coeff: ERROR '
        write(6,*)' contraction range for basis is 0:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      mycoeffptr = infbs_cont(CONT_ICFP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      mygen   = infbs_cont(CONT_NGEN,myucont,basis)
      call dcopy((myprim*mygen),
     &    dbl_mb(mb_exndcf(mycoeffptr,basis)),1,coeff,1)
c
      bas_get_coeff = .true.
c
      return
      end
*.....................................................................
      logical function bas_set_exponent(basisin,icont,exp,nexp)
c
c sets the exponents to the input array for a given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin ! [input] basis set handle        
      integer icont ! [input] mapped contraction index
      integer nexp ! [input] number of exponents for contraction
      double precision exp(nexp) ! [input] "new" exponents for 
*. . . . . . . . . . . . . . . .            contraction    
c::local
      integer basis, myucont, icontmax
      integer myprim,myexptr
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      bas_set_exponent = bas_check_handle(basisin,'bas_set_exponent')
      if (.not.bas_set_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)

      bas_set_exponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_set_exponent)) then
        write(6,*)' bas_set_exponent: ERROR '
        write(6,*)' contraction range for basis is 0:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      bas_set_exponent = myprim.eq.nexp
      if (.not.bas_set_exponent) then
        write(6,*)' bas_set_exponent: ERROR '
        write(6,*)' input and stored number of exponents ',
     &         '(nprim) differ '
        write(6,*)' input  nprim: ',nexp
        write(6,*)' stored nprim: ',myprim
        return
      endif
      call dcopy(nexp,exp,1,dbl_mb(mb_exndcf(myexptr,basis)),1)
c
      bas_set_exponent = .true.
c
      return
      end
*.....................................................................
      logical function bas_set_coeff(basisin,icont,coeff,ncoeff)
c
c sets the coefficients of the given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin ! [input] basis set handle                   
      integer icont ! [input] mapped contraction index           
      integer ncoeff ! [input] number of coeffs. for contraction
      double precision coeff(ncoeff) ! [input] "new" coeffs. 
*. . . . . . . . . . . . . . . . . .            for contraction    
c::local
      integer basis, myucont, icontmax
      integer mycoeffptr, myprim, mygen
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      bas_set_coeff = bas_check_handle(basisin,'bas_set_coeff')
      if (.not.bas_set_coeff) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)

      bas_set_coeff = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_set_coeff)) then
        write(6,*)' bas_set_coeff: ERROR '
        write(6,*)' contraction range for basis is 0:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      mycoeffptr = infbs_cont(CONT_ICFP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      mygen   = infbs_cont(CONT_NGEN,myucont,basis)
c
      bas_set_coeff = ncoeff .eq. (myprim*mygen)
      if(.not.bas_set_coeff) then
        write(6,*)' bas_set_coeff: ERROR '
        write(6,*)' input and stored number of coefficients ',
     &         '(nprim*ngen) differ '
        write(6,*)' input  nprim*ngen: ',ncoeff
        write(6,*)' stored nprim*ngen: ',(myprim*mygen)
        return
      endif
      call dcopy(ncoeff,coeff,1,dbl_mb(mb_exndcf(mycoeffptr,basis)),1)
c
      bas_set_coeff = .true.
c
      return
      end
*.....................................................................
      logical function bas_nbf_cn_max(basisin,nbf_max)
      implicit none
c
c  calculate, return and store maximum basis function block size
c   for all contractions in a given basis. 
c
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
      integer nbf_from_ucont
      external nbf_from_ucont
c:: passed
      integer basisin ! [input] basis set handle
      integer nbf_max ! [output] largest basis function block on any 
*. . . . . . . . . . .           shell in the given basis set.
c:local
      integer basis, myucont, i, mynbf
c
      bas_nbf_cn_max = bas_check_handle(basisin,'bas_nbf_cn_max')
      if (.not. bas_nbf_cn_max ) then
        write(6,*) 'bas_nbf_cn_max: basis handle not valid '
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
      if (nbfmax_bs(basis) .gt. -565) then
        nbf_max = nbfmax_bs(basis)
        bas_nbf_cn_max = .true.
        return
      endif
c
      myucont = infbs_head(HEAD_NCONT,basis)
      nbf_max = -565
c
      do 00100 i = 1,myucont
        mynbf = nbf_from_ucont(i,basisin)
        nbf_max = max(nbf_max, mynbf)
00100 continue
c
      nbfmax_bs(basis) = nbf_max
      bas_nbf_cn_max = .true.
      return
      end
*.....................................................................
      logical function bas_nbf_ce_max(basisin,nbf_max)
      implicit none
c
c  calculate, return and store maximum basis function block size
c   for all contractions in a given basis. 
c
#include "mafdecls.fh"
#include "basP.fh"
#include "geom.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin       ! [input] basis set handle
      integer nbf_max       ! [output] largest nbf block on an 
c                           !          atom in the given basis.
c:local
      integer basis, mynat, iat, mylo, myhi, mynbf
c
#include "bas_ibs_sfn.fh"
c
      bas_nbf_ce_max = bas_check_handle(basisin,'bas_nbf_ce_max')
      if (.not. bas_nbf_ce_max ) then
        write(6,*) 'bas_nbf_ce_max: basis handle not valid '
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      bas_nbf_ce_max = geom_ncent(ibs_geom(basis),mynat)
      if (mynat.le.0.or. .not.bas_nbf_ce_max) then
        write(6,*)' bas_nbf_ce_max: ERROR '
        write(6,*)' number of centers is zero or weird'
        write(6,*)' nat = ',mynat
        return
      endif
c
      nbf_max = -1
      do 00100 iat = 1,mynat
        mylo  = sf_ibs_cn2bfr(1,(sf_ibs_ce2cnr(1,iat,basis)),basis)
        myhi  = sf_ibs_cn2bfr(2,(sf_ibs_ce2cnr(2,iat,basis)),basis)
        mynbf = myhi - mylo + 1
        nbf_max = max(nbf_max, mynbf)
00100 continue
c
      return
      end
*.....................................................................
      logical function bas_geom(basisin,geom)	
      implicit none 
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin       ! [input] basis set handle
      integer basis
      integer geom          ! [output] geometry used to load basis set
c
      if (.not.bas_check_handle(basisin,'bas_geom'))
     & call errquit('bas_geom: handle invalid',911)
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      geom = ibs_geom(basis)
      bas_geom = .true.
      end
*.....................................................................
      logical function bas_ncoef_cn_max(basisin, ncoef_max)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
      integer basisin   ! [input] basis set handle
      integer ncoef_max ! [output] maximum number of coefs
c
      integer basis
      integer nu_cont, iu_cont, my_prim, my_gen
c
      basis = basisin + BASIS_HANDLE_OFFSET
      if (.not.bsactive(basis)) call errquit
     &    ('bas_ncoef_cn_max: basis handle invalid',911)
      nu_cont = infbs_head(HEAD_NCONT,basis)
      ncoef_max = 0
      do 00100 iu_cont = 1,nu_cont
        my_prim = infbs_cont(CONT_NPRIM,iu_cont,basis)
        my_gen  = infbs_cont(CONT_NGEN, iu_cont,basis)
        ncoef_max = max(ncoef_max,(my_prim*my_gen))
00100 continue
      bas_ncoef_cn_max = .true.
      end
*.....................................................................
      logical function bas_nprim_cn_max(basisin, nprim_max)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
      integer basisin   ! [input] basis set handle
      integer nprim_max ! [output] maximum number of prims in any 
*. . . . . . . . . . . .           contraction
c
      integer basis
      integer nu_cont, iu_cont, my_prim
c
      basis = basisin + BASIS_HANDLE_OFFSET
      if (.not.bsactive(basis)) call errquit
     &    ('bas_nprim_cn_max: basis handle invalid',911)
      nu_cont = infbs_head(HEAD_NCONT,basis)
      nprim_max = 0
      do 00100 iu_cont = 1,nu_cont
        my_prim = infbs_cont(CONT_NPRIM,iu_cont,basis)
        nprim_max = max(nprim_max,my_prim)
00100 continue
      bas_nprim_cn_max = .true.
      end
*.....................................................................
      logical function bas_norm_get(basisin,norm_id)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
c
      integer basisin ! [input] basis set handle
      integer norm_id ! [output] Normalization id type
c
      integer basis
c
      basis = basisin + BASIS_HANDLE_OFFSET
      norm_id = bas_norm_id(basis)
      bas_norm_get = norm_id.ge.BasNorm_lo.and.norm_id.le.BasNorm_hi
      end
*.....................................................................
      logical function bas_norm_set(basisin,norm_id)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
c
      integer basisin ! [input] basis set handle
      integer norm_id ! [input] Normalization id type
c
      integer basis
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      if (norm_id.ge.BasNorm_lo.and.norm_id.le.BasNorm_hi) then
        bas_norm_id(basis) = norm_id
        bas_norm_set = .true.
      else
        bas_norm_set = .false.
      endif
      end
*.....................................................................
      logical function bas_norm_print(basisin)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "stdio.fh"
c
      integer basisin ! [input] basis set handle
c
      integer basis, norm_id
c
      basis = basisin + BASIS_HANDLE_OFFSET
      norm_id = bas_norm_id(basis)
      bas_norm_print = .true.
      if      (norm_id.eq.BasNorm_UN)  then
        write(luout,*)' basis is unnormalized'
      else if (norm_id.eq.BasNorm_STD) then
        write(luout,*)' basis has standard normalization via int_norm'
      else if (norm_id.eq.BasNorm_2c)  then
        write(luout,*)
     &        ' basis has dft/fitting normalization via int_norm_2c'
      else
        write(luout,*)' basis handle: ',basisin
        write(luout,*)' norm_id = ',norm_id
        write(luout,*)' does not match a known normalization mode'
        bas_norm_print = .false.
      endif
      end
*.....................................................................
      logical function ecp_get_rexponent(ecpin,icont,exp)
c
c returns the r-exponents of a given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c::function
      logical ecp_check_handle
      external ecp_check_handle
c:blas
c     dcopy
c::passed
      integer ecpin         ! [input] ecp set handle
      integer icont           ! [input] mapped contraction index
      double precision exp(*) ! [output] exponents
c::local
      integer ecp, myucont, icontmax
      integer myprim,myexptr
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      ecp_get_rexponent = ecp_check_handle(ecpin,'ecp_get_rexponent')
      if (.not.ecp_get_rexponent) return

      ecp = ecpin + BASIS_HANDLE_OFFSET 
      icontmax = ncont_tot_gb(ecp)
      myucont = sf_ibs_cn2ucn(icont,ecp)
      ecp_get_rexponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(ecp_get_rexponent)) then
        write(6,*)' ecp_get_rexponent: ERROR '
        write(6,*)' contraction range for ecp is 0:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IREXP,myucont,ecp)
      myprim  = infbs_cont(CONT_NPRIM,myucont,ecp)
      call dcopy(myprim,dbl_mb(mb_exndcf(myexptr,ecp)),1,exp,1)
c
      ecp_get_rexponent = .true.
c
      return
      end
      
*.....................................................................
      logical function ecp_set_rexponent(ecpin,icont,exp,nexp)
c
c sets the exponents to the input array for a given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
c::function
      logical ecp_check_handle
      external ecp_check_handle
c:blas
c     dcopy
c::passed
      integer ecpin ! [input] ecp set handle        
      integer icont ! [input] mapped contraction index
      integer nexp ! [input] number of exponents for contraction
      double precision exp(nexp) ! [input] "new" r-exponents 
*. . . . . . . . . . . . . . . .            for contraction    
c::local
      integer ecp, myucont, icontmax
      integer myprim,myexptr
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      ecp_set_rexponent = ecp_check_handle(ecpin,'ecp_set_rexponent')
      if (.not.ecp_set_rexponent) return

      ecp = ecpin + BASIS_HANDLE_OFFSET 

      icontmax = ncont_tot_gb(ecp)
      myucont = sf_ibs_cn2ucn(icont,ecp)

      ecp_set_rexponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(ecp_set_rexponent)) then
        write(6,*)' ecp_set_rexponent: ERROR '
        write(6,*)' contraction range for ecp is 0:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,ecp)
      myprim  = infbs_cont(CONT_NPRIM,myucont,ecp)
      ecp_set_rexponent = myprim.eq.nexp
      if (.not.ecp_set_rexponent) then
        write(6,*)' ecp_set_rexponent: ERROR '
        write(6,*)' input and stored number of exponents ',
     &         '(nprim) differ '
        write(6,*)' input  nprim: ',nexp
        write(6,*)' stored nprim: ',myprim
        return
      endif
      call dcopy(nexp,exp,1,dbl_mb(mb_exndcf(myexptr,ecp)),1)
c
      ecp_set_rexponent = .true.
c
      return
      end
      logical function ecp_print(ecpidin)
c
c routine to print unique ecpid information that is in core
c
      implicit none
#include "mafdecls.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
#include "geom.fh"
c
c function declarations
c      
      logical ecp_check_handle
      external ecp_check_handle
c:: passed
      integer ecpidin ! [input] basis set handle
c:: local
      integer mytags, myucont, myprim, mycoef, ecpid
      integer i,j,k,l, ifcont, mygen, mytype, irexptr, iexptr, icfptr
      integer atn, len_tag, len_ele
      character*2 symbol
      character*16 element
      character*3 ctype(-1:6)
      character*3 shell_type
      character*60 buffer 
c
#include "bas_exndcf.fh"
c
      ctype(-1)='U L'
      ctype(0) ='U-s'
      ctype(1) ='U-p'
      ctype(2) ='U-d'
      ctype(3) ='U-f'
      ctype(4) ='U-g'
      ctype(5) ='U-h'
      ctype(6) ='U-i'
      ecp_print = .true.
      ecpid = ecpidin + BASIS_HANDLE_OFFSET
c
      ecp_print = ecp_check_handle(ecpidin,'ecp_print')
      if (.not. ecp_print) return
c
c print basis set information
c      
      write(6,1)bs_name(ecpid)(1:inp_strlen(bs_name(ecpid))), 
     $    bs_trans(ecpid)(1:inp_strlen(bs_trans(ecpid)))
 1    format('                 ECP Basis "',a,'" -> "',a,'"'/
     $    '                      -----')
      if (bas_spherical(ecpid)) write(6,2)
 2    format('    ECP Basis is spherical, 5d, 7f, 9g ... ')
      mytags  = infbs_head(HEAD_NTAGS,ecpid)
      if (mytags.le.0) then
        write(6,*)' Empty ECP Basis set '
        write(6,*)
        return
      endif
c
      myucont = infbs_head(HEAD_NCONT,ecpid)
      myprim  = infbs_head(HEAD_NPRIM,ecpid)
      mycoef  = infbs_head(HEAD_NCOEF,ecpid)
c
      do 00100 i=1,mytags
        
        if (geom_tag_to_element(bs_tags(i,ecpid), symbol, element,
     $      atn)) then
          len_tag = inp_strlen(bs_tags(i,ecpid))
          len_ele = inp_strlen(element)
          write (buffer,
     &        '(a,'' ('',a,'') Replaces '',i5,'' electrons'' )')
     &        bs_tags(i,ecpid)(1:len_tag), element(1:len_ele),
     &        infbs_tags(Tag_Nelec,i,ecpid)
        else
          buffer = bs_tags(i,ecpid)
        endif
        len_tag = inp_strlen(buffer)
        call util_print_centered(6, buffer, len_tag/2 + 1, .true.)
        
        myucont = infbs_tags(TAG_NCONT,i,ecpid)
c        
        ifcont = infbs_tags(TAG_FCONT,i,ecpid)
c     
        write(6,6)
 6      format(
     $      '          R-exponent    Exponent     Coefficients '/
     $      '         ------------ ',60('-'))
        do 00200 j=1,myucont
          myprim = infbs_cont(CONT_NPRIM,ifcont,ecpid)
          mygen  = infbs_cont(CONT_NGEN,ifcont,ecpid)
          
          mytype = infbs_cont(CONT_TYPE, ifcont, ecpid)
          shell_type = ctype(mytype)
          iexptr  = infbs_cont(CONT_IEXP, ifcont,ecpid) - 1
          icfptr  = infbs_cont(CONT_ICFP, ifcont,ecpid) - 1
          irexptr = infbs_cont(CONT_IREXP,ifcont,ecpid) - 1
          do 00300 k=1,myprim
            write(6,7) j, shell_type,
     &          sf_exndcf((irexptr+k),ecpid),
     &          sf_exndcf((iexptr+k),ecpid),
     &          (sf_exndcf((icfptr+k+(l-1)*myprim),ecpid),l=1,mygen)
00300     continue
          write(6,*)
          ifcont = ifcont + 1
00200   continue
00100 continue
c
c  If geom is set print out the info about total basis info 
c  associated with the geometry also
c
c  ... not done yet
c
      return
 7          format(1x,i2,1x,a3,1x,f9.2,2x,f14.6,20f15.6)
      end
*.....................................................................
      logical function bas_set_ecp_name(basisin,ecp_name)
      implicit none
#include "basP.fh"
      integer basisin
      character*(*) ecp_name
c
      ecp_name_assoc((basisin+basis_handle_offset)) =
     &    ecp_name
      bas_set_ecp_name = .true.
      end
*.....................................................................
      logical function bas_get_ecp_name(basisin,ecp_name)
      implicit none
#include "basP.fh"
      integer basisin
      character*(*) ecp_name
c
      ecp_name =
     &    ecp_name_assoc((basisin+basis_handle_offset)) 
      bas_get_ecp_name = .true.
      end
*.....................................................................
      logical function bas_set_ecp_handle(basisin,ecp_handle)
      implicit none
#include "basP.fh"
      integer basisin
      integer ecp_handle
c
      ecp_handle_assoc((basisin+basis_handle_offset)) =
     &    ecp_handle
      bas_set_ecp_handle = .true.
      end
*.....................................................................
      logical function bas_get_ecp_handle(basisin,ecp_handle)
      implicit none
#include "basP.fh"
      integer basisin
      integer ecp_handle
c
      ecp_handle =
     &    ecp_handle_assoc((basisin+basis_handle_offset)) 
      if (ecp_handle.ne.0) then
        bas_get_ecp_handle = .true.
      else
        bas_get_ecp_handle = .false.
      endif
      end
*.....................................................................
      logical function ecp_set_parent_handle(ecp_handle,basis_handle)
      implicit none
#include "basP.fh"
      integer ecp_handle
      integer basis_handle
c
      logical bas_check_handle, ecp_check_handle
      external bas_check_handle, ecp_check_handle
c      
      if (.not.bas_check_handle(basis_handle,'ecp_set_parent_handle'))
     &      call errquit
     &      ('ecp_set_parent_handle: basis_handle invalid',911)
      if (.not.ecp_check_handle(ecp_handle,'ecp_set_parent_handle'))
     &      call errquit
     &      ('ecp_set_parent_handle: ecp_handle invalid',911)
      ecp_parent_assoc(ecp_handle+basis_handle_offset) =
     &      basis_handle
      ecp_set_parent_handle = .true.
      end
*.....................................................................
      logical function ecp_get_parent_handle(ecp_handle,basis_handle)
      implicit none
#include "basP.fh"
      integer ecp_handle
      integer basis_handle
c
      logical bas_check_handle, ecp_check_handle
      external bas_check_handle, ecp_check_handle
c      
      if (.not.ecp_check_handle(ecp_handle,'ecp_get_parent_handle'))
     &      call errquit
     &      ('ecp_get_parent_handle: ecp_handle invalid',911)
      basis_handle = ecp_parent_assoc(ecp_handle+basis_handle_offset)
      if (basis_handle.ne.0) then
        ecp_get_parent_handle = .true.
      else
        ecp_get_parent_handle = .false.
      endif
      if (.not.bas_check_handle(basis_handle,'ecp_get_parent_handle'))
     &      call errquit
     &      ('ecp_get_parent_handle: stored basis_handle invalid',911)
      end
      subroutine bas_print_allocated_info(msg)
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
      character*(*) msg
      integer basis
      integer me, nproc, inode
      integer inode_use
c
#include "bas_ibs_sfn.fh"
#include "bas_exndcf_sfn.fh"
c
      me = ga_nodeid()
      nproc = ga_nnodes()
*
      do inode = 0,(nproc-1)
        inode_use = inode
        call ga_sync()
        call ga_brdcst(1234,inode_use,mitob(1),0)
        call ga_sync()
        if (inode_use.eq.me) then
          write(luout,*)' msg: ',msg,me
          write(luout,*)' basis data for node ',me
          write(luout,*)' number of possible basis sets',nbasis_bsmx
          do basis = 1,nbasis_bsmx
            write(luout,*)' active :',bsactive(basis)
            write(luout,*)' exndcf Handle/index/size :',
     &          exndcf(H_exndcf,basis),exndcf(K_exndcf,basis),
     &          exndcf(SZ_exndcf,basis)
            write(luout,*)' cn2ucn (H/K/SZ) :',
     &          ibs_cn2ucn(H_ibs,basis),
     &          ibs_cn2ucn(K_ibs,basis),
     &          ibs_cn2ucn(SZ_ibs,basis)
            write(luout,*)' cn2ce (H/K/SZ)  :',
     &          ibs_cn2ce(H_ibs,basis),
     &          ibs_cn2ce(K_ibs,basis),
     &          ibs_cn2ce(SZ_ibs,basis)
            write(luout,*)' ce2uce (H/K/SZ) :',
     &          ibs_ce2uce(H_ibs,basis),
     &          ibs_ce2uce(K_ibs,basis),
     &          ibs_ce2uce(SZ_ibs,basis)
            write(luout,*)' cn2bfr (H/K/SZ) :',
     &          ibs_cn2bfr(H_ibs,basis),
     &          ibs_cn2bfr(K_ibs,basis),
     &          ibs_cn2bfr(SZ_ibs,basis)
            write(luout,*)' ce2cnr (H/K/SZ) :',
     &          ibs_ce2cnr(H_ibs,basis),
     &          ibs_ce2cnr(K_ibs,basis),
     &          ibs_ce2cnr(SZ_ibs,basis)
          enddo
        endif
        call util_flush(luout)
      enddo
      end
      logical function bas_match_tags(tag_from_geom,basisin,btag)
      implicit none
*
* This routine matches geometry tags to basis set tags in such a way
*     that "H34" matches "H" or "hydrogen" or "h" if "H34" does not 
*     exist in the basis set specification.  
*     "H" will however not match "h" i.e., case sensitivity is 
*     preserved.
*
c::includes
#include "stdio.fh"
#include "inp.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      logical geom_tag_to_element
      external geom_tag_to_element
c::passed
      character*(*) tag_from_geom ! [input] geometry tag
      integer basisin             ! [input] basis set handle
      integer btag                ! [output] lexical basis set tag
*     bas_match_tags              ! [output] true if they match
c::local
      integer basis           ! lexical basis index
      integer nbtgs           ! number of basis set tags
      character*16 gstring    ! geometry or basis tags can only be *16
      character*16 bsmatch    ! the matched basis set tag
      integer lgstring        ! length of gstring
      integer lnotalpha
      integer i               ! loop counter
      integer ind             ! dummy counter
*rak:      logical digits_found    ! did geometry tag have any digits
      logical status          ! dummy storage
      character*2 g_sym       ! geometry tag -> symbol name
      character*16 g_elem     ! geometry tag -> element name
      integer g_atn           ! geometry tag -> atomic number
      character*2 b_sym       ! basis set tag -> symbol name
      character*16 b_elem     ! basis set tag -> element name
      integer b_atn           ! basis set tag -> atomic number
      logical debug           ! true for extra output
c
      integer na2z
      parameter (na2z = 26)
      character*1 a2z(na2z)
      data a2z /'a','b','c','d','e','f','g','h','i','j',
     &          'k','l','m','n','o','p','q','r','s','t',
     &          'u','v','w','x','y','z'/
*rak:      integer ndigits         ! number of digits
*rak:      parameter (ndigits=10)  ! set number of digits
*rak:      character*1 digits(ndigits)  ! array of character digits
*rak:c      
*rak:      data digits /'0','1','2','3','4','5','6','7','8','9'/
c      
      debug = .true.
c
      bas_match_tags = .false.
c
      basis = basisin+basis_handle_offset
c
      gstring = tag_from_geom               !  copy tag to work on it
c
c... first match full geometry tag to full basis tag list 
c          did the user specifically assign a tag ??

      nbtgs = infbs_head(HEAD_NTAGS,basis)
      do i = 1, nbtgs
        if (gstring.eq.bs_tags(i,basis)) then
          bas_match_tags = .true.
          btag = i
          bsmatch = bs_tags(i,basis)
          goto 00009
        endif
      enddo
c
c... Now check to see if the tag has non alpha chracters. 
c...  the substring prior to the first non-alpha character is the users idea
c...  of the "name" of the tag.  
c
      lnotalpha = lgstring+1
      do i = 1,lgstring
        if (.not.(inp_match(na2z,.false.,
     &      gstring(lgstring:lgstring),a2z,ind))) then ! compare character to alpha (case-less test)
          lnotalpha = ind
          goto 00001
        endif
      enddo
00001 continue
      do i = lnotalpha,lgstring
        gstring(i:i) = ' '
      enddo
*rak:
*rak:      digits_found = .false.
*rak:00001 continue
*rak:      lgstring = inp_strlen(gstring)   !  get the length of the geometry tag
*rak:      if (lgstring.eq.0) return        !  if empty string or string of digits return false
*rak:      if (inp_match(ndigits,.false.,
*rak:     &    gstring(lgstring:lgstring),digits,i)) then   ! compare last character to a digit
*rak:        gstring(lgstring:lgstring) = ' '             ! if a digit remove it
*rak:        digits_found = .true.
*rak:        goto 00001
*rak:      endif
*rak:c
*rak:c if no digits then enforce case matching between say "H" and "h" by a false return      
*rak:c
*rak:      if (.not.digits_found) return
*rak:c
*rak:c... first match numberless geometry tag to full basis tag list
c
c... match user substring to basis set tags
      do i = 1, nbtgs
        if (gstring.eq.bs_tags(i,basis)) then
          bas_match_tags = .true.
          btag = i
          bsmatch = bs_tags(i,basis)
          goto 00009
        endif
      enddo
c
c... now get symbol and element names for each tag and match those
c    geometry tag is based on users substring
c    basis set tag is based on user input to the basis object
c
      status = geom_tag_to_element(gstring,g_sym,g_elem,g_atn)
      if (.not.status)then
        if (.not.(g_sym.eq.'bq')) then
          write(luout,*)'geometry tag',gstring,
     &        ' could not be matched to an element symbol'
          call errquit('bas_match_tags: fatal error ',911)
        endif
      endif
      do i = 1, nbtgs
        status =
     &      geom_tag_to_element(bs_tags(i,basis),b_sym,b_elem,b_atn)
        if (.not.status) then
          if (.not.(b_sym.eq.'bq')) then
            write(luout,*)'basis tag',bs_tags(i,basis),
     &          ' could not be matched to an element symbol'
            call errquit('bas_match_tags: fatal error ',911)
          endif
        endif
        if (g_elem.eq.b_elem) then
          bas_match_tags = .true.
          btag = i
          bsmatch = bs_tags(i,basis)
          goto 00009
        else if (g_sym.eq.b_sym) then
          bas_match_tags = .true.
          btag = i
          bsmatch = bs_tags(i,basis)
          goto 00009
        endif
      enddo
      write(luout,*)'bas_match_tags:debug: no match for tag <',
     &    tag_from_geom,'>'
      return  ! no match
c
00009 continue
      if (debug) then
        write(luout,10000)tag_from_geom,bsmatch
      endif
      return
10000 format('bas_match_tags:debug: geometry tag ',a16,
     &    ' matched basis set tag ',a16)
      end
      logical function bas_cando_mpoles(basis)
      implicit none
      integer basis             ! [input]
c
c     Return true if it's possible to compute multipoles for this basis
c
c     The multipole code cannot handle sphericals, general contractions
c     or SP shells.
c
      integer nshell, ishell, type, nprim, ngen, sphcart
      logical bas_numcont, bas_continfo
      external bas_numcont, bas_continfo
c
      bas_cando_mpoles = .false.
c
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('multipole: basis bad?',0)
      do ishell = 1, nshell
         if (.not. bas_continfo(basis, ishell, type, 
     $        nprim, ngen, sphcart)) call errquit
     $        ('multipole: basis bad?',0)
         if (type.lt.0 .or. ngen.gt.1 .or. sphcart.ne.0) return
      enddo
c
      bas_cando_mpoles = .true.
c
      end
