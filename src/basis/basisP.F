      subroutine print_exndcf(msg)
      implicit none
c $Id$
#include "nwc_const.fh"
#include "basP.fh"
#include "util.fh"
      character *(*) msg
      integer i,j
      write(6,*)'<<<< print_exndcf: ',msg,' >>>>'
      do i = 1, nbasis_bsmx
        write(6,10000)i,(exndcf(j,i), j=1,3)
      enddo
10000 format(' basis:',i2,' handle :',i10,/
     &       '         ', ' index  :',i10,/
     &       '         ', ' size   :',i10/)
      end
*.....................................................................
      logical function bas_add_ucnt_init(basisin)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "bas_exndcf_dec.fh"
      integer basisin
c::local
      integer input_size
      parameter(input_size = 10000)
      integer basis, h_tmp, k_tmp
c
#include "bas_exndcf_sfn.fh"
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      bas_add_ucnt_init = ma_alloc_get(
     &    mt_dbl,input_size,' input for basis heap ',
     &    h_tmp,k_tmp)
      exndcf(H_exndcf,basis) = h_tmp
      exndcf(K_exndcf,basis)= k_tmp
      if (.not. bas_add_ucnt_init) call errquit
     &    ('bas_add_unct_init: error allocating input heap space',911,
     &       MEM_ERR)
      call dfill(input_size,0.0d00,dbl_mb(k_tmp),1)
      exndcf(SZ_exndcf,basis) = input_size
      end
*.....................................................................
      logical function bas_add_ucnt_tidy(basisin)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "bas_exndcf_dec.fh"
      integer basisin
c::local
      integer basis 
      integer h_tmp
c
#include "bas_exndcf_sfn.fh"
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      h_tmp = exndcf(H_exndcf,basis)
      bas_add_ucnt_tidy = ma_free_heap(h_tmp)
      if (.not.bas_add_ucnt_tidy) call errquit
     &    ('bas_add_unct_tidy: error freeing heap',911, MEM_ERR)
      exndcf(H_exndcf ,basis) = -1
      exndcf(K_exndcf ,basis) = 0
      exndcf(SZ_exndcf,basis) = 0
      bas_add_ucnt_tidy = .true.
      end
*.....................................................................
      logical function bas_set_ecp_basis(basis)
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
      logical bas_check_handle
      external bas_check_handle
c      
      integer basis
c
      integer bas
c
      bas_set_ecp_basis = bas_check_handle(basis,'bas_set_ecp_basis')
      if (.not. bas_set_ecp_basis) return
      bas = basis + BASIS_HANDLE_OFFSET
c
      infbs_head(Head_ECP,bas) = 1
      end
*.....................................................................
      logical function bas_set_so_basis(basis)
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
      logical bas_check_handle
      external bas_check_handle
c      
      integer basis
c
      integer bas
c
      bas_set_so_basis = bas_check_handle(basis,'bas_set_so_basis')
      if (.not. bas_set_so_basis) return
      bas = basis + BASIS_HANDLE_OFFSET
c
      infbs_head(Head_ECP,bas) = 2
      end
*.....................................................................
      logical function ecp_set_num_elec(ecpid,tag,num_elec,stdtag)
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
*functions::
      logical ecp_check_handle
      logical bas_add_utag
      external ecp_check_handle
      external bas_add_utag
*passed::
      integer ecpid      ! [input] basis set handle
      character*16 tag   ! [input] tag to set number of electrons for
      integer num_elec   ! [input] number of electrons ecp replaces
      character*(*) stdtag ! [input] name associated with ecp on tag
c
*local::
      integer ecp_indx
      integer itag
c
c     sets the number of electrons replaced by the ecp on the 
c     specified tag.  If the tag is not present it will also add 
c     that by calling bas_add_utag
c      
      ecp_set_num_elec = ecp_check_handle(ecpid,'ecp_set_num_elec')
      if (.not. ecp_set_num_elec) return
      ecp_indx = ecpid + BASIS_HANDLE_OFFSET
c
c     Make sure that the tag is in the list
c
      ecp_set_num_elec = bas_add_utag(ecpid, tag, stdtag, itag)
      if (.not. ecp_set_num_elec) return
      infbs_tags(Tag_Nelec,itag,ecp_indx) = num_elec
c
      end
*.....................................................................
      logical function ecp_get_num_elec(ecpid,tag,num_elec)
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "stdio.fh"
*functions::
      logical ecp_check_handle
      external ecp_check_handle
*passed::
      integer ecpid      ! [input] basis set handle
      character*16 tag   ! [input] tag to get number of electrons for
      integer num_elec   ! [input] number of electrons ecp replaces
c
*local::
      integer ecp_indx
      integer itag, ntag
c
c     gets the number of electrons replaced by the ecp on the 
c     specified tag.  If the tag is not present it will error
*
c
      ecp_get_num_elec = ecp_check_handle(ecpid,'ecp_get_num_elec')
      if (.not.ecp_get_num_elec) return
      ecp_indx = ecpid + BASIS_HANDLE_OFFSET
c
      ntag = infbs_head(Head_Ntags,ecp_indx)
      do itag = 1,ntag
        if (tag.eq.bs_tags(itag,ecp_indx)) then
          num_elec = infbs_tags(Tag_Nelec,itag,ecp_indx)
          ecp_get_num_elec = .true.
          return
        endif
      enddo
c
      num_elec = 0
      ecp_get_num_elec = .false.
c
      end
*.....................................................................
      logical function bas_add_utag(basisin, tag, stdtag, itag)
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "inp.fh"
      integer basisin           ! [input] basis handle
      character*(*) tag         ! [input] name of tag
      character*(*) stdtag      ! [input] name of basis set on tag
      integer itag              ! [output] index of tag
c
      integer basis             ! [local] index into basis arrays
      logical bas_check_handle
      external bas_check_handle
      integer tmp
c
c     Add the unique tag to the list of tags in the basis,
c     incrementing the no. of tags if necessary.
c     Return in itag the index of the unique tag
c
      bas_add_utag = bas_check_handle(basisin, 'bas_add_utag')
      if (.not. bas_add_utag) return
      basis = basisin + BASIS_HANDLE_OFFSET
c
      do itag = 1, infbs_head(HEAD_NTAGS,basis)
        if (bs_tags(itag,basis) .eq. tag) then
          if (bs_stdname(itag,basis).eq.'unknown') then
            bs_stdname(itag,basis) = stdtag
          else if (bs_stdname(itag,basis) .ne. stdtag) then
            if (.not.(bs_stdname(itag,basis)(1:9).eq.'modified:')) then
              tmp = inp_strlen(bs_stdname(itag,basis))
              bs_stdname(itag,basis) =
     &            'modified:'//bs_stdname(itag,basis)(1:tmp)
            endif
          endif
          return             
        endif
      enddo
c
c     No match found ... append new tag to the list
c
      itag = infbs_head(HEAD_NTAGS,basis) + 1
      if (itag .gt. ntags_bsmx) then
         write(6,*) 'bas_add_utag: too many tags', itag
         bas_add_utag = .false.
         return
      endif
c
      infbs_head(HEAD_NTAGS,basis) = itag
      bs_tags(itag,basis) = tag
      bs_stdname(itag,basis) = stdtag
c
      end
*.....................................................................
      subroutine bas_err_info(info)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
c     
      character*(*) info        ! [input]
      integer bas,basin
      integer nbas
      logical status
c     
c     For internal use of the basis set routines only: print out
c     info of known basis sets to aid in diagnosing a problem
c     
c::function
      logical bas_print
      external bas_print
c
      nbas = 0
      do 00100 bas = 1, nbasis_bsmx
        if (bsactive(bas)) nbas = nbas + 1
00100 continue
      write(6,'(1x,a,a,i2)')
     &       info, ': open basis sets:',nbas
c
      nbas = 0
      do 00200 bas = 1, nbasis_bsmx
        if (bsactive(bas)) then
          basin = bas - BASIS_HANDLE_OFFSET
          status = bas_print(basin)
        endif
00200 continue
c
      if (nbasis_rtdb .gt. 0) then
        write(6,'(1x,a,a,i3)')
     &         info,': basis sets in current rtdb ',nbasis_rtdb
        do 00300 bas = 1, nbasis_rtdb
          write(6,'(1x,a,1x,i3,3x,a,1x,a)')
     &           'number:',bas,
     &           'basis set name:',
     &           bs_names_rtdb(bas)(1:len_bs_rtdb(bas))
00300   continue
      endif
c     
      end
*.....................................................................
      logical function bas_ucontinfo(basisin,icont,itype,
     &       nprimo,ngeno,sphcart)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin, icont, nprimo, ngeno, sphcart, itype
c::local
      integer basis,myucont,icontmax
c
      nprimo = -123
      ngeno  = -456
      sphcart = -789
c
      bas_ucontinfo = bas_check_handle(basisin,'bas_ucontinfo')
      if (.not.bas_ucontinfo) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = infbs_head(HEAD_NCONT,basis)
c
      if (.not.(icont.gt.0.and.icont.le.icontmax)) then
        write(6,*)' bas_continfo: ERROR '
        write(6,*)' unique contraction range for basis is 1:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        bas_ucontinfo = .false.
        return
      endif
c
      myucont = icont
      if (bas_spherical(basis)) then
        sphcart = 1
      else
        sphcart = 0
      endif
      itype   = infbs_cont(CONT_TYPE,myucont,basis)
      nprimo  = infbs_cont(CONT_NPRIM,myucont,basis)
      ngeno   = infbs_cont(CONT_NGEN,myucont,basis)
      bas_ucontinfo=.true.
      return
      end
*.....................................................................
      logical function bas_unumcont(basisin,numcont)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin,numcont
c::local
      integer basis
c
      numcont = -6589
      bas_unumcont = bas_check_handle(basisin,'bas_numcont')
      if (.not.bas_unumcont) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      numcont = infbs_head(HEAD_NCONT,basis)

      bas_unumcont = .true.
      return
      end
*.....................................................................
      block data basis_data
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
c     
      data nbasis_rtdb /0/
      data bsactive /nbasis_bsmx*.false./
      data bas_spherical /nbasis_bsmx*.false./
      data angular_bs /nbasis_bsmx*-565/
      data bas_norm_id /nbasis_bsmx*-565/
      data nbfmax_bs /nbasis_bsmx*-565/
      data bsversion /5.00d00/
* version 5 includes so stuff
c     
      end
*.....................................................................
      integer function nbf_from_ucont(ucont,basisin)
c
c function that returns the number of basis functions in a contraction
c 
c  types 0->S, 1->P, 2->D, 3->F etc. -1->SP -2->SPD
c
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c:: function
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer ucont   ! [input] unique contraction 
      integer basisin ! [input] basis set handle
c:: local
      integer type
      integer basis
      integer ngen
c
      if(.not.bas_check_handle(basisin,'nbf_from_ucont'))
     &  call errquit('nbf_from_ucont: bad basis handle',basisin,
     &       BASIS_ERR)
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      type = infbs_cont(CONT_TYPE,ucont,basis)
      ngen = infbs_cont(CONT_NGEN,ucont,basis)
      if (type.ge.0) then 
        if (bas_spherical(basis)) then
          nbf_from_ucont = ngen*(2*type+1)
        else
          nbf_from_ucont = ngen*(type+1)*(type+2)/2
        endif
      else if (type.eq.-1) then
        nbf_from_ucont = ngen*2
      else if (type.eq.-2) then
        if (bas_spherical(basis)) then
          nbf_from_ucont = ngen*9/3
        else
          nbf_from_ucont = ngen*10/3
        endif
      else 
        call errquit('nbf_from_ucont: bad cont type',type, BASIS_ERR)
      endif
      end
*.....................................................................
      logical function bas_set_spherical(basisin, ospherical)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c
      integer basisin    ! [input] basis set handle
      logical ospherical ! [input] logical for spherical setting
c
      integer basis
c
      bas_set_spherical = .false.
c
      basis = basisin + BASIS_HANDLE_OFFSET
      if (ospherical) then
        infbs_head(HEAD_SPH,basis) = 1
        bas_spherical(basis) = .true.
      else
        infbs_head(HEAD_SPH,basis) = 0
        bas_spherical(basis) = .false.
      endif
      bas_set_spherical = .true.
      end
*.....................................................................
      logical function bas_get_spherical(basisin, is_spherical)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c
      integer basisin      ! [input] basis get handle
      logical is_spherical ! [output] logical for spherical getting
c
      integer basis
c
      bas_get_spherical = .false.
c
      basis = basisin + BASIS_HANDLE_OFFSET
      if (infbs_head(head_sph,basis).eq.1) then 
        is_spherical = .true.
      else
        is_spherical = .false.
      endif
      bas_get_spherical = .true.
      end
*.....................................................................
      logical function bas_name_exist_rtdb(rtdb,name)
      implicit none
#include "errquit.fh"
*
* function to determin if "name" has been stored on the 
* current rtdb in actual or translated by a set directive.
*
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
#include "inp.fh"      
#include "nwc_const.fh"
#include "basP.fh"      
c::functions
      logical  bas_rtdb_in
      external bas_rtdb_in
c::passed 
      integer rtdb ! [input] run time data base handle
      character*(*) name ! [input] test name
c::local
      integer index
      character*256 trans_name
c
      bas_name_exist_rtdb = bas_rtdb_in(rtdb)
      if (.not.bas_name_exist_rtdb) call errquit
     &    ('bas_name_exist_rtdb: bas_rtdb_in failed',911, RTDB_ERR)
c
      bas_name_exist_rtdb = .false.
      if (inp_match
     &    (nbasis_rtdb,.false.,name,bs_names_rtdb,index)) then
        bas_name_exist_rtdb = .true.
        return
      endif
      if (context_rtdb_match(rtdb,name,trans_name)) then
        if (inp_match
     &      (nbasis_rtdb,.false.,trans_name,bs_names_rtdb,index)) then
          bas_name_exist_rtdb = .true.
          return
        endif
      endif
      end
