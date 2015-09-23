      logical function bas_version()
c $Id$
c
c: Routine that calclulates the size of the common block structures 
c  used in the basis set object and the mapped representation object.
c:input none
c:output always true.
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
c
      integer cdata,idata,rdata
      integer mapidata, total4, total8
c
c.. character data 
      cdata =         256*nbasis_bsmx            ! bs_name
      cdata = cdata + 16*ntags_bsmx*nbasis_bsmx  ! bs_tags
      cdata = cdata + 80*ntags_bsmx*nbasis_bsmx  ! bs_stdname
      cdata = cdata + 256*nbasis_bsmx            ! bs_trans
      cdata = cdata + 256*nbasis_rtdb_mx         ! bs_names_rtdb
      cdata = cdata + 256*nbasis_bsmx*nbasis_assoc_max ! name_assoc
c      
c.. real data 
      rdata =  1                   ! bsversion
      rdata = 8*rdata
c
c.. integer data in basis set object common
      idata =         3*nbasis_bsmx                ! exndcf
      idata = idata + ndbs_head*nbasis_bsmx        ! infbs_head
      idata = idata +
     &    ndbs_tags*ntags_bsmx*nbasis_bsmx         ! infbs_tags
      idata = idata +
     &    ndbs_ucont*(nucont_bsmx+1)*nbasis_bsmx   ! infbs_cont
      idata = idata + nbasis_bsmx                  ! len_bs_name
      idata = idata + nbasis_bsmx                  ! len_bs_trans
      idata = idata + nbasis_rtdb_mx               ! len_bs_rtdb
      idata = idata + nbasis_bsmx                  ! bsactive (int==logical)
      idata = idata + nbasis_bsmx                  ! bas_spherical (int==logical)
      idata = idata + nbasis_bsmx                  ! bas_any_gc (int==logical)
      idata = idata + nbasis_bsmx                  ! bas_any_sp_shell (int==logical)
      idata = idata + nbasis_bsmx                  ! bas_norm_id
      idata = idata + nbasis_bsmx                  ! angular_bs
      idata = idata + nbasis_bsmx                  ! nbfmax_bs
      idata = idata + nbasis_assoc_max*nbasis_bsmx ! handle_assoc
      idata = idata + nbasis_assoc_max*nbasis_bsmx ! parent_assoc
      idata = idata + 1                            ! nbasis_rtdb
      idata = 4*idata
c
c.. integer data in the mapped object.  
      mapidata =            3*nbasis_bsmx   ! ibs_cn2ucn
      mapidata = mapidata + 3*nbasis_bsmx   ! ibs_cn2ce
      mapidata = mapidata + 3*nbasis_bsmx   ! ibs_ce2uce
      mapidata = mapidata + 3*nbasis_bsmx   ! ibs_cn2bfr
      mapidata = mapidata + 3*nbasis_bsmx   ! ibs_ce2cnr
      mapidata = mapidata + nbasis_bsmx     ! ncont_tot_gb
      mapidata = mapidata + nbasis_bsmx     ! nprim_tot_gb
      mapidata = mapidata + nbasis_bsmx     ! nbf_tot_gb
      mapidata = mapidata + nbasis_bsmx     ! ibs_geom
      mapidata = 4*mapidata 
c
c.. total space 
      total4 = idata + mapidata
      total8 = 2*total4 + rdata + cdata 
      total4 = total4 + rdata + cdata 
c
      write(LuOut,'(////1x,a,f5.2,a)')
     &      ' **** basis set version ',bsversion,' ****'
      write(LuOut,'(1x,a,i20,a)')
     &      '   character data in-core ',cdata,' bytes'
      write(LuOut,'(1x,a,i20,a)')
     &      '   real      data in-core ',rdata,' bytes'
      write(LuOut,'(1x,a,i20,a)')
     &      '   integer*4 data in-core ',idata,' bytes'
      write(LuOut,'(1x,a,i20,a)')
     &      'or integer*8 data in-core ',(2*idata),' bytes'
      write(LuOut,'(1x,a,i20,a)')
     &      '   integer*4 mapping data in-core ',
     &      mapidata,' bytes'
      write(LuOut,'(1x,a,i20,a/)')
     &      'or integer*8 mapping data in-core ',
     &      (2*mapidata),' bytes'
      write(LuOut,*)' total(4)   = ',total4,' bytes'
      write(LuOut,*)' total(8)   = ',total8,' bytes'
c
c.. convert to kilobytes
c
      cdata    = (cdata    + 999) / 1000
      rdata    = (rdata    + 999) / 1000
      idata    = (idata    + 999) / 1000
      mapidata = (mapidata + 999) / 1000
      total4   = (total4   + 999) / 1000
      total8   = (total8   + 999) / 1000
      write(LuOut,'(///1x,a,f5.2,a)')
     &      ' **** basis set version ',bsversion,' ****'
      write(LuOut,'(1x,a,i20,a)')
     &      '   character data in-core ',cdata,' Kbytes'
      write(LuOut,'(1x,a,i20,a)')
     &      '   real      data in-core ',rdata,' Kbytes'
      write(LuOut,'(1x,a,i20,a)')
     &      '   integer*4 data in-core ',idata,' Kbytes'
      write(LuOut,'(1x,a,i20,a)')
     &      'or integer*8 data in-core ',(2*idata),' Kbytes'
      write(LuOut,'(1x,a,i20,a)')
     &      '   integer*4 mapping data in-core ',
     &      mapidata,' Kbytes'
      write(LuOut,'(1x,a,i20,a/)')
     &      'or integer*8 mapping data in-core ',
     &      (2*mapidata),' Kbytes'
      write(LuOut,*)' total(4)   = ',total4,' Kbytes'
      write(LuOut,*)' total(8)   = ',total8,' Kbytes'
      write(LuOut,'(////)')
c
      bas_version = .true.
      end
*.....................................................................
C> \ingroup bas
C> @{
c
C> \brief Creates a new basis instance
c
C> \return Returns .true. if successfull, and .false. otherwise
c
      logical function bas_create(basis,name)
c
c creates a handle and marks it active in the in-core data structure
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
c::functions
c ifill from util
c dfill from util
c::passed
c
      integer basis      !< [Output] the basis handle
      character*(*) name !< [Input] the name of basis set  
c
*:debug:      integer iiii, jjjj
      integer ii ! dummy loop counter
c
      external basis_data  ! This for the T3D linker
c
#include "bas_exndcf_sfn.fh"
c
      do 00100 basis=1,nbasis_bsmx
        if (.not.bsactive(basis)) goto 01000
00100 continue
c
      write(LuOut,*)' bas_create: no free basis handles for ',name
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
      exndcf(H_exndcf ,basis) = -1 ! handle 
      exndcf(K_exndcf ,basis) = 0 ! index
      exndcf(SZ_exndcf,basis) = 0 ! allocated size
      call ifill(ndbs_tags*ntags_bsmx, 0,
     $      infbs_tags(1,1,basis), 1)
      call ifill(ndbs_ucont*nucont_bsmx, 0,
     $      infbs_cont(1,1,basis), 1)
      do ii = 1,ntags_bsmx
         bs_stdname(ii,basis) = ' '
         bs_tags(ii,basis) = ' '
      enddo
      do ii = 1,nbasis_assoc_max
        name_assoc(ii,basis) = ' '
      enddo
      call ifill(nbasis_assoc_max,0,handle_assoc(1,basis),1)
      call ifill(nbasis_assoc_max,0,parent_assoc(1,basis),1)
      bas_nassoc(basis) = 0
c
c Initialize geo-basis info to empty
c
      call ifill(3, 0, ibs_cn2ucn(1,basis), 1)
      call ifill(3, 0, ibs_cn2ce (1,basis), 1)
      call ifill(3, 0, ibs_ce2uce(1,basis), 1)
      call ifill(3, 0, ibs_cn2bfr(1,basis), 1)
      call ifill(3, 0, ibs_ce2cnr(1,basis), 1)
      nbfmax_bs(basis)    = -565
      ncont_tot_gb(basis) = 0
      nprim_tot_gb(basis) = 0
      ncoef_tot_gb(basis) = 0
      nbf_tot_gb(basis)   = 0
      ibs_geom(basis)     = 0
      bas_norm_id(basis)  = BasNorm_UN
c
c Mark basis as active and return info
c
      bsactive(basis) = .true.
      basis = basis - Basis_Handle_Offset 
      bas_create = .true.
c
c debug print all basis sets that are active
c
*:debug:      write(LuOut,*)' bas_create: active basis sets'
*:debug:      do iiii = 1,nbasis_bsmx
*:debug:        if (bsactive(iiii)) then
*:debug:          jjjj = inp_strlen(bs_name(iiii))
*:debug:          write(LuOut,*)'bas_create:',bs_name(iiii)(1:jjjj),' ',iiii
*:debug:        endif
*:debug:      enddo
c
      end
*.....................................................................
c
C> \brief Destroys a basis instance
c
C> \return Returns .true. if successfull, and .false. otherwise
c
      logical function bas_destroy(basisin)
      implicit none
c::functions
      logical bas_get_ecp_handle, ecp_check_handle, bas_do_destroy
      external bas_get_ecp_handle, ecp_check_handle, bas_do_destroy
      logical bas_get_so_handle, so_check_handle
      external bas_get_so_handle, so_check_handle
c::passed
      integer basisin !< [Input] the basis handle
c::local
      logical ignore
      integer ecpid, soid
c
      ignore = bas_get_ecp_handle(basisin,ecpid)
      ignore = bas_get_so_handle(basisin,soid)
      bas_destroy = .true.
      if (ecp_check_handle(ecpid,'bas_destroy')) then
        bas_destroy = bas_destroy.and.bas_do_destroy(ecpid)
      endif
      if (so_check_handle(soid,'bas_destroy')) then
        bas_destroy = bas_destroy.and.bas_do_destroy(soid)
      endif
      bas_destroy = bas_destroy .and. bas_do_destroy(basisin)
      end
*.....................................................................
C> @}
      logical function bas_do_destroy(basisin)
c
c destroys information about an active incore basis
c and the associated mapping arrays.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "bas_exndcf_dec.fh"
#include "ecpso_decP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::functions used
c::bas
      logical bas_check_handle
      logical gbs_map_clear
      logical geom_ncent
      logical ecp_get_num_elec
      external gbs_map_clear
      external bas_check_handle
      external geom_ncent
      external ecp_get_num_elec
c::passed
      integer basisin ![input] basis set handle to be destroyed
c::local
      integer idbstag
      integer i, nat, num_elec
      integer geom
      integer basis
      integer h_tmp
c
#include "bas_ibs_sfn.fh"
#include "bas_exndcf_sfn.fh"
#include "ecpso_sfnP.fh"
c
      bas_do_destroy = bas_check_handle(basisin,'bas_do_destroy')
      if (.not. bas_do_destroy) return

      basis = basisin + Basis_Handle_Offset 

c
c must restore active geometry data appropriately
      geom = ibs_geom(basis)
      if (active(geom)) then
        if (Is_ECP(basis)) then
          if (.not.geom_ncent(geom,nat))
     &        call errquit('bas_do_destroy:geom_ncent failed',911,
     &       BASIS_ERR)
          do i = 1,nat
            idbstag = sf_ibs_ce2uce(i,basis)
            if (ecp_get_num_elec(basisin,
     &          bs_tags(idbstag,basis),num_elec)) then
              charge(i,geom) = charge(i,geom) + dble(num_elec)
            endif
          enddo
        endif
      endif
c
      if(.not.gbs_map_clear(basisin)) then
        write(LuOut,*)' error clearing map '
        bas_do_destroy = .false.
        return
      endif

c
      angular_bs(basis)       = -565
      bas_norm_id(basis)      = -565
      nbfmax_bs(basis)        = -565
      bsactive(basis)         = .false.
      bas_spherical(basis)    = .false.
      bas_any_gc(basis)       = .false.
      bas_any_sp_shell(basis) = .false.
c
      h_tmp = exndcf(H_exndcf,basis)
      if (h_tmp .ne. -1) then
         if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_do_destroy: error freeing heap',911, MEM_ERR)
      endif
      exndcf(H_exndcf ,basis) = -1
      exndcf(K_exndcf ,basis) = 0
      exndcf(SZ_exndcf,basis) = 0
c
      bas_do_destroy = .true.
      end
*.....................................................................
C> \ingroup bas
C> @{
c
C> \brief Checks whether a given handle refers to a valid basis instance
c
C> \return Return .true. if basisin is a valid basis instance, and
C> .false. otherwise
c
      logical function bas_check_handle(basisin,msg)
c
c Checks to see if a basis set handle is valid
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "stdio.fh"
c:passed
      integer basisin   !< [Input] the basis handle
      character*(*) msg !< [Input] an error message
c::local
      integer basis
c
      basis = basisin + Basis_Handle_Offset
      bas_check_handle = basis.gt.0 .and. basis.le.nbasis_bsmx
      if (bas_check_handle)
     &      bas_check_handle = bas_check_handle .and. bsactive(basis)
c
* user's responsibility to deal with status
*      if (.not. bas_check_handle) then
*        write(LuOut,*)msg,': basis handle is invalid '
*        write(LuOut,*)'basis_check_handle: lexical handle ',basis
*        write(LuOut,*)'basis_check_handle:         handle ',basisin
*      endif
      return
      end
*.....................................................................
c
C> \brief Checks whether a given handle refers to a valid ECP basis
C> instance
c
C> \return Return .true. if ecpidin is a valid ECP basis instance, and
C> .false. otherwise
c
      logical function ecp_check_handle(ecpidin,msg)
c
c Checks to see if an ECP basis set handle is valid
c
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "ecpso_decP.fh"
#include "inp.fh"
c:functions
      logical bas_check_handle
      external bas_check_handle
c:passed
      integer ecpidin   !< [Input] the ECP basis handle
      character*(*) msg !< [Input] an error message
c::local
      character*255 newmsg    
c
#include "ecpso_sfnP.fh"
c
      newmsg(1:17) = 'ecp_check_handle:'
      newmsg(18:) = msg 
      ecp_check_handle = bas_check_handle(ecpidin,newmsg)
      if(.not.ecp_check_handle) return
      ecp_check_handle = ecp_check_handle .and.
     &    Is_ECP_in(ecpidin)
      return
      end
*.....................................................................
c
C> \brief Checks whether a given handle refers to a valid spin-orbit
C> potential instance
c
C> \return Return .true. if soidin is a valid spin-orbit potential
C> instance, and .false. otherwise
c
      logical function so_check_handle(soidin,msg)
c
c Checks to see if an ECP so basis set handle is valid
c
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "ecpso_decP.fh"
#include "inp.fh"
c:functions
      logical bas_check_handle
      external bas_check_handle
c:passed
      integer soidin    !< [Input] the spin-orbit ECP basis handle
      character*(*) msg !< [Input] an error message
c::local
      character*255 newmsg    
c
#include "ecpso_sfnP.fh"
c
      newmsg(1:17) = 'so_check_handle:'
      newmsg(18:) = msg 
      so_check_handle = bas_check_handle(soidin,newmsg)
      if(.not.so_check_handle) return
      so_check_handle = so_check_handle .and.
     &    Is_SO_in(soidin)
      return
      end
*.....................................................................
c
C> \brief Print a summary of a basis instance
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function bas_summary_print(basisin)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
#include "ecpso_decP.fh"
#include "bas_starP.fh"
c::-functions
      logical bas_check_handle
      external bas_check_handle
      integer nbf_from_ucont
      external nbf_from_ucont
c::-passed
      integer basisin !< [Input] the basis set handle
c::-local
      character*16 dum_tag
      character*255 dum_stdtag
      character*12 polynomial  ! string for spherical/cartesian
      character*80 dum_type, dum_string
      integer basis   ! lexical index
      integer i_tag   ! loop counter
      integer i_cont  ! loop counter
      integer mytags  ! dummy 
      integer myfcont ! loop range value (first)
      integer mylcont ! loop range value (last)
      integer mycont  ! number of contractions
      integer mynbf   ! number of functions
      integer tmp1, tmp2, tmp3
      integer myngen
      integer mytype
      integer jtype
      character*1 ctype(0:6)
      integer cnt_type(0:6)
*
#include "ecpso_sfnP.fh"
*
      bas_summary_print =
     &    bas_check_handle(basisin,'bas_summary_print')
*
      if (Is_ECP_in(basisin).or.Is_SO_in(basisin)) then
        bas_summary_print = .true.
        return
      endif
*
      ctype(0)='s'
      ctype(1)='p'
      ctype(2)='d'
      ctype(3)='f'
      ctype(4)='g'
      ctype(5)='h'
      ctype(6)='i'
*
      basis = basisin + Basis_Handle_Offset
*      
      if (bas_spherical(basis)) then
         polynomial = ' (spherical)'
      else
         polynomial = ' (cartesian)'
      endif
      write(luout,'(/)')
      write(luout,10000)
     &    bs_name(basis)(1:inp_strlen(bs_name(basis))), 
     &    bs_trans(basis)(1:inp_strlen(bs_trans(basis))),
     &    polynomial
      mytags = infbs_head(Head_Ntags,basis)
      do i_tag = 1,mytags
        call ifill(7,0,cnt_type,1)
        myfcont = infbs_tags(Tag_Fcont,i_tag,basis)
        mylcont = infbs_tags(Tag_Lcont,i_tag,basis)
        mycont  = mylcont - myfcont + 1
        mynbf   = 0 
        do i_cont = myfcont, mylcont
          mynbf = mynbf + nbf_from_ucont(i_cont,basisin)
          mytype = infbs_cont(Cont_Type,i_cont,basis)
          myngen = infbs_cont(Cont_Ngen,i_cont,basis)
          if (myngen.lt.1) call errquit (
     &        'bas_summary_print: fatal error myngen:',myngen,
     &       BASIS_ERR)
          if (mytype.ge.0) then
            cnt_type(mytype) = cnt_type(mytype) + myngen
          else 
            do jtype = 0,abs(mytype)
              cnt_type(jtype) = cnt_type(jtype) + 1
            enddo
          endif
        enddo
        dum_tag = bs_tags(i_tag,basis)
        tmp1 = inp_strlen(bs_stdname(i_tag,basis))
        if (tmp1 .lt. (30-1)) then
          tmp2 = (30-tmp1)/2
        else
          tmp2 = 1
        endif
        dum_stdtag = ' '
        dum_stdtag(tmp2:) = bs_stdname(i_tag,basis)
        dum_type = ' '
        dum_string = ' '
        do jtype = 0,6
          dum_type = dum_string
          tmp1 = inp_strlen(dum_type)
          if (tmp1.lt.1) tmp1 = 1
          if (cnt_type(jtype).gt.0) then
            write(dum_string,'(a,i5,a)')dum_type(1:tmp1),
     &        cnt_type(jtype),ctype(jtype)
          endif
        enddo
        tmp1 = 0
        dum_type = dum_string
        do jtype = 1,inp_strlen(dum_type)
          if (dum_type(jtype:jtype).ne.char(0).and.
     &        dum_type(jtype:jtype).ne.' ') then
            tmp1 = tmp1 + 1
            dum_string(tmp1:tmp1) = dum_type(jtype:jtype)
          endif
        enddo
        tmp1 = tmp1 + 1
        do jtype = tmp1,len(dum_string)
          dum_string(jtype:jtype) = ' '
        enddo
        tmp1 = inp_strlen(dum_string)
        write(luout,10001)
     &      dum_tag,
     &      dum_stdtag(1:30),
     &      mycont,
     &      mynbf,dum_string(1:tmp1)
      enddo
      do i_tag=1,star_nr_tags
        tmp1 = inp_strlen(star_bas_typ(i_tag))
        if (tmp1 .lt. (30-1)) then
          tmp2 = (30-tmp1)/2
        else
          tmp2 = 1
        endif
        dum_stdtag = ' '
        dum_stdtag(tmp2:) = star_bas_typ(i_tag)
        tmp1 = 1
        if (i_tag .gt. 1) tmp1 = star_nr_excpt(i_tag-1) + 1
        if ((star_nr_excpt(i_tag) - tmp1) .gt. -1) then
           write(LuOut,10002) star_tag(i_tag), 
     &           dum_stdtag(1:30), 'all atoms except',
     &           (star_excpt(tmp3)(1:inp_strlen(star_excpt(tmp3))),
     &           tmp3=tmp1,star_nr_excpt(i_tag))
        else
           write(LuOut,10002) star_tag(i_tag), 
     &           dum_stdtag(1:30), 'on all atoms '
        endif
      enddo

      write(luout,'(/)')
      call util_flush(luout)
      bas_summary_print = .true.
*
10000 format(' Summary of "',a,'" -> "',a,'"',a/
     $    ' ---------------------------------------------',
     &    '---------------------------------'/
     $    '       Tag                 Description            ',
     &    'Shells   Functions and Types'/
     $    ' ---------------- ------------------------------  ',
     &    '------  ---------------------')
10001 format(1x,a16,1x,a30,2x,i4,5x,i4,3x,a)
10002 format(1x,a16,1x,a30,2x,a17,1x,10(a))
*
      end
*.....................................................................
c
C> \brief Print the contents of a basis instance
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function bas_print(basisin)
c
c routine to print unique basis information that is in core
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
#include "geom.fh"
#include "ecpso_decP.fh"
#include "stdio.fh"
#include "bas_starP.fh"
c
c function declarations
c      
      logical  bas_check_handle, ecp_print, so_print
      external bas_check_handle, ecp_print, so_print
c:: passed
      integer basisin !< [Input] the basis set handle
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
#include "ecpso_sfnP.fh"
c
      if (Is_ECP_in(basisin)) then
        bas_print = ecp_print(basisin)
        return
      endif
      if (Is_SO_in(basisin)) then
        bas_print = so_print(basisin)
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
      basis = basisin + Basis_Handle_Offset
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
      write(LuOut,1)bs_name(basis)(1:inp_strlen(bs_name(basis))), 
     $     bs_trans(basis)(1:inp_strlen(bs_trans(basis))), polynomial
 1    format('                      Basis "',a,'" -> "',a,'"',a/
     $       '                      -----')
      mytags  = infbs_head(HEAD_NTAGS,basis)
      if (mytags.le.0) then
        write(LuOut,*)'No explicit basis set is defined !'
        write(LuOut,*)
c
c there could be star tags defined, so check that before returning
c
        goto 00010
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
         call util_print_centered(LuOut, buffer, len_tag/2 + 1, .true.)

         myucont = infbs_tags(TAG_NCONT,i,basis)
c        
        ifcont = infbs_tags(TAG_FCONT,i,basis)
c     
        write(LuOut,6)
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
            write(LuOut,7) j, shell_type(1:2),
     &          sf_exndcf((iexptr+k),basis),
     &          (sf_exndcf((icfptr+k+(l-1)*myprim),basis),l=1,mygen)
 7           format(1x,i2,1x,a2,1x,1pE14.8,0p20f10.6)
00300     continue
          write(LuOut,*)
          ifcont = ifcont + 1
00200   continue
00100 continue
c
c  Check if we have star tag definitions in the basis set
c
00010 if (star_nr_tags .gt. 0) then
         write(LuOut,*)
         write(LuOut,8) 
  8      format(' In addition, one or more string tags have been',
     &         ' defined containing a * .'/' These tags, and ',
     &         'their exceptions list are printed below.'//
     &         ' Tag ',12(' '),' Description ',18(' '), 
     &         ' Excluding '/' ',16('-'),' ',30('-'),' ',30('-'))
         do i=1,star_nr_tags
           k = 1
           if (i .gt. 1) k = star_nr_excpt(i-1) + 1
           write(LuOut,9) star_tag(i), star_bas_typ(i)(1:30),
     &          (star_excpt(j)(1:inp_strlen(star_excpt(j))),
     &           j=k,star_nr_excpt(i))
  9        format(1x,a16,1x,a30,1x,10(a)) 
         enddo
         write(LuOut,*)
         write(LuOut,*)
      endif
c
c  If geom is set print out the info about total basis info 
c  associated with the geometry also
c
c  ... not done yet
c
      return
      end
*.....................................................................
c
C> \brief Load a basis set from the RTDB
c
C> Retrieve the basis set from the RTDB using the name specified. 
C> The data is stored in the basis set instance specified. Furthermore
C> some tables are constructed linking basis functions to atoms, etc.
C> For this purpose a geometry needs to be provided as well.
c
C> \return Return .true. if the basis was loaded successfully, and
C> .false. otherwise
c
      logical function bas_rtdb_load(rtdb, geom, basisin, name)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "context.fh"
#include "stdio.fh"
#include "inp.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      logical  bas_rtdb_do_load, bas_get_ecp_name, bas_set_ecp_name
      logical  bas_create, bas_do_destroy, bas_set_ecp_handle
      logical  ecp_set_parent_handle, bas_summary_print
      logical  bas_name_exist_rtdb, bas_print
      external bas_rtdb_do_load, bas_get_ecp_name, bas_set_ecp_name
      external bas_create, bas_do_destroy, bas_set_ecp_handle
      external ecp_set_parent_handle, bas_summary_print
      external bas_name_exist_rtdb, bas_print
      logical bas_get_so_name, bas_set_so_name, bas_set_so_handle
      logical so_set_parent_handle
      external bas_get_so_name, bas_set_so_name, bas_set_so_handle
      external so_set_parent_handle
c::passed
      integer rtdb       !< [Input] the RTDB handle
      integer geom       !< [Input] the geometry handle
      integer basisin    !< [Input] the basis set handle
      character*(*) name !< [Input] the name of the basis set on the
c                        !< RTDB
c
c::local
      character*256 ecp_name, so_name
      integer ecpid, soid
      logical status
      logical status_ecp_cr, status_so_cr
      logical status_ecp_load, status_so_load
      logical status_ecp_ph, status_so_ph
      logical status_ecp, status_so
c
*:debug:      write(LuOut,*)' bas_rtdb_load:rtdb,geom,basis ',
*:debug:     &    rtdb,geom,basisin
*:debug:      write(LuOut,*)' bas_rtdb_load:name ',name
      bas_rtdb_load =
     &    bas_rtdb_do_load(rtdb, geom, basisin, name)
      if (.not.bas_rtdb_load) return
      if (.not.inp_compare(.false.,'ao basis',name(1:8)))
     &    goto 112
*
      if (.not.bas_get_ecp_name(basisin,ecp_name)) call errquit
     &    ('bas_rtdb_load: bas_get_ecp_name failed',911, RTDB_ERR)
      if (ecp_name .eq. '  ') then
        ecp_name = 'ecp basis'
      endif

      status = bas_name_exist_rtdb(rtdb,ecp_name)

      if (status) then
*
* ecp_name is on rtdb so load it
*
        status_ecp_cr = bas_create(ecpid,ecp_name)
        if (.not.status_ecp_cr) then
          call errquit
     &        ('bas_rtdb_load: bas_create failed for ecpid',911,
     &       RTDB_ERR)
        endif

        status_ecp_load = bas_rtdb_do_load(rtdb,geom,ecpid,ecp_name)
        if (status_ecp_load) then
*.... everything is okay
          if (.not.bas_set_ecp_name(basisin,ecp_name)) call errquit
     &        ('bas_rtdb_load: bas_set_ecp_name failed',911,
     &       RTDB_ERR)

          status_ecp = bas_set_ecp_handle(basisin,ecpid)
          if (.not.status_ecp) then
            call errquit
     &          ('bas_rtdb_load: bas_set_ecp_handle failed',911,
     &       RTDB_ERR)
          endif
        
          status_ecp_ph = ecp_set_parent_handle(ecpid,basisin)
          if (.not.status_ecp_ph) then
            call errquit
     &          ('bas_rtdb_load: ecp_set_parent_handle failed',911,
     &       RTDB_ERR)
          endif
        else
*... ecp exists but is not used by current geometry so destroy it !
          if (.not.bas_do_destroy(ecpid)) then
            write(luout,*)' unused ecp basis failed to be destroyed'
            call errquit('bas_rtdb_load: bas_do_destroy failed',911,
     &       RTDB_ERR)
          endif
        endif
      endif

      if (.not.bas_get_so_name(basisin,so_name)) call errquit
     &    ('bas_rtdb_load: bas_get_so_name failed',911, RTDB_ERR)
      if (so_name .eq. '  ') then
        so_name = 'so potential'
      endif

      status = bas_name_exist_rtdb(rtdb,so_name)

      if (status) then
*
* so_name is on rtdb so load it
*
        status_so_cr = bas_create(soid,so_name)
        if (.not.status_so_cr) then
          call errquit
     &        ('bas_rtdb_load: bas_create failed for soid',911,
     &       RTDB_ERR)
        endif

        status_so_load = bas_rtdb_do_load(rtdb,geom,soid,so_name)
        if (status_so_load) then
*.... everything is okay
          if (.not.bas_set_so_name(basisin,so_name)) call errquit
     &        ('bas_rtdb_load: bas_set_so_name failed',911,
     &       RTDB_ERR)

          status_so = bas_set_so_handle(basisin,soid)
          if (.not.status_so) then
            call errquit
     &          ('bas_rtdb_load: bas_set_so_handle failed',911,
     &       RTDB_ERR)
          endif
        
          status_so_ph = so_set_parent_handle(soid,basisin)
          if (.not.status_so_ph) then
            call errquit
     &          ('bas_rtdb_load: so_set_parent_handle failed',911,
     &       RTDB_ERR)
          endif
        else
*... so exists but is not used by current geometry so destroy it !
          if (.not.bas_do_destroy(soid)) then
            write(luout,*)' unused so basis failed to be destroyed'
            call errquit('bas_rtdb_load: bas_do_destroy failed',911,
     &       RTDB_ERR)
          endif
        endif
      endif
c
c Check if we need to print the basis set. The basis set might not 
c have been printed, this happens when the user asks for the basis
c to be printed but the basis set contains star tags, which only can 
c can be resolved at the task level.
c
c Reusing ecp_name and so_name and status
c
 112  if (.not.context_rtdb_match(rtdb,name,ecp_name))
     $    ecp_name = name
      so_name = 'basisprint:'//ecp_name
     $           (1:inp_strlen(ecp_name))
      if (rtdb_get(rtdb,so_name,mt_log,1,status)) then
         if (status) then
            if (ga_nodeid() .eq. 0) then 
             if (.not. bas_print(basisin))
     $          call errquit('bas_rtdb_load: print failed', 0,
     &       RTDB_ERR)
             if (.not.bas_summary_print(basisin)) call errquit
     $          ('bas_rtdb_load: basis summary print failed',911,
     &       RTDB_ERR)
            endif
            status = .false.
            if (.not. rtdb_put(rtdb,so_name,mt_log,1,status))
     $         call errquit('bas_rtdb_load: rtdb_put failed',911,
     &       RTDB_ERR)
         endif
      endif
      call bas_ecce_print_basis(basisin,'bas_input')
c
      end
*.....................................................................
C> @}
      logical function bas_rtdb_do_load(rtdb, geom, basisin, name)
      implicit none
#include "errquit.fh"
c
c routine that loads a basis set from the rtdb and using the 
c geometry information builds the mapping arrays to contractions/
c shells, basis functions, and centers.
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "global.fh"
#include "stdio.fh"
#include "bas_starP.fh"
#include "bas.fh"
c
c function declarations
c
c::function
c:bas
      logical bas_geobas_build
      logical bas_add_ucnt_init
      logical bas_add_ucnt_tidy
      external bas_geobas_build
      external bas_add_ucnt_init
      external bas_add_ucnt_tidy
c:: passed  
      integer rtdb        ! [input] rtdb handle      
      integer geom        ! [input] geometry handle with info loaded
      integer basisin     ! [input] basis handle
      character*(*) name  ! [input] basis set name that must be on 
*. . . . . . . . . . . . .          the rtdb
c:: local
      integer ecXtra      ! amount of extra space required for zero shell info
      parameter (ecXtra = 2)
      character*25 nameex
      integer lentmp, basis, nexcf
      character*256 tmp
      logical rtdb_status
      integer h_tmp, k_tmp, h_new, k_new
      integer iunique, uniquecent, istar , ntag_read
      integer istart, iexcpt, ilen, junique
      integer i_array,l_array
      logical oIsexcept
      character*16 tag_to_add, tag_in_lib
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
      basis = basisin + Basis_Handle_Offset
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
      ntag_read=0
      tmp(lentmp:) = ' '
      tmp(lentmp:) = ':bs_nr_tags'
      rtdb_status = rtdb_status .and. 
     &              rtdb_get(rtdb,tmp,mt_int,1,ntag_read)
      if (ntag_read .gt. 0) then
         tmp(lentmp:) = ' ' 
         tmp(lentmp:) = ':bs_tags'
         rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, ntags_bsmx, bs_tags(1,basis))
c
         tmp(lentmp:) = ' ' 
         tmp(lentmp:) = ':bs_stdname'
         rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, ntags_bsmx, bs_stdname(1,basis))
      endif
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':assoc ecp name'
      rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, 2, name_assoc(1,basis))
c
      nexcf = 0
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':number of exps and coeffs'
      rtdb_status = rtdb_status .and.
     &    rtdb_get(rtdb,tmp,mt_int,1,nexcf)
c
      write(nameex,'(a23,i2)')' basis exps and coeffs ',basis
c
      if (exndcf(H_exndcf,basis) .ne. -1) then
         if (.not. ma_free_heap(exndcf(H_exndcf,basis)))
     $        call errquit('bas_rtdb_do_load: ma is corrupted',
     $        exndcf(H_exndcf,basis), RTDB_ERR)
      endif

      if (.not.ma_alloc_get(mt_dbl,(nexcf+ecXtra),nameex,
     &    h_tmp, k_tmp)) then
        write(LuOut,*)' not enough memory'
        call errquit
     &      (' bas_rtdb_do_load: error allocating space'//
     &      ' for exndcf',911, RTDB_ERR)
      else
        call dfill((nexcf+ecXtra),0.0d00,dbl_mb(k_tmp),1)
        exndcf(H_exndcf,basis) = h_tmp
        exndcf(K_exndcf,basis) = k_tmp
        exndcf(SZ_exndcf,basis) = (nexcf+ecXtra)
      endif
c
      if (nexcf .gt. 0) then
         tmp(lentmp:) = ' ' 
         tmp(lentmp:) = ':exps and coeffs'
         rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_dbl, nexcf, dbl_mb(mb_exndcf(1,basis)))
      endif
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':header'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int, ndbs_head, infbs_head(1,basis))
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':tags info'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int,
     &       ndbs_tags*ntags_bsmx, infbs_tags(1,1,basis))
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':contraction info'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int,
     &       ndbs_ucont*nucont_bsmx, infbs_cont(1,1,basis))
c
c now deal with the star tags from the input. We have the geometry, so
c we extract the tags from the geometry, load the star tag info, and
c add the basis sets
c
c first check if we have any star tags to deal with at all
c account for old rtdbs without star tag info
c
      tmp(lentmp:) = ' '
      tmp(lentmp:) = ':star nr tags'
      star_nr_tags = 0
      if (rtdb_status) rtdb_status = rtdb_status .and. rtdb_get(rtdb,
     &       tmp,mt_int,1,star_nr_tags)
c
c read remaining star tag info and analyze only if we have actually a star
c tag input (hence, if star_nr_tags > 0)
c
      if (star_nr_tags .gt. 0) then
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star tag names'
        rtdb_status = rtdb_status .and. rtdb_cget(
     &         rtdb,tmp,star_nr_tags,star_tag)
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star tag_in_lib'
        rtdb_status = rtdb_status .and. rtdb_cget(
     &         rtdb,tmp,star_nr_tags,star_in_lib)
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star bas type'
        rtdb_status = rtdb_status .and. rtdb_cget(
     &         rtdb,tmp,star_nr_tags,star_bas_typ)
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star filename'
        rtdb_status = rtdb_status .and. rtdb_cget(
     &         rtdb,tmp,star_nr_tags,star_file)
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star tot except'
        rtdb_status = rtdb_status .and. rtdb_get(
     &         rtdb,tmp,mt_int,1,star_tot_excpt)
        if (star_tot_excpt .gt. 0) then
           tmp(lentmp:) = ' '
           tmp(lentmp:) = ':star except'
           rtdb_status = rtdb_status .and. rtdb_cget(
     &            rtdb,tmp,star_tot_excpt,star_excpt)
        endif
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star nr except'
        rtdb_status = rtdb_status .and. rtdb_get(
     &         rtdb,tmp,mt_int,star_nr_tags,star_nr_excpt)
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star rel'
        rtdb_status = rtdb_status .and. rtdb_get(
     &         rtdb,tmp,mt_log,star_nr_tags,star_rel)
        tmp(lentmp:) = ' '
        tmp(lentmp:) = ':star segment'
        rtdb_status = rtdb_status .and. rtdb_get(
     &         rtdb,tmp,mt_log,1,star_segment)
c
c get a list of tags from the current geometry
c
        if (.not. geom_ncent_unique(geom,uniquecent)) 
     &     call errquit('bas_rtdb_do_load: geom_ncent_unique',211,
     &       RTDB_ERR)
      if (.not.MA_alloc_Get(mt_int,uniquecent, 'iarray', 
     &     l_array, i_array))
     &     call errquit('basrtdbdoload: cannot allocate ',0, MA_ERR)
        if (.not. geom_uniquecent_get(geom,uniquecent,
     ,     int_mb(i_array))) 
     &     call errquit('bas_rtdb_do_load: geom_uniquecent_get',211,
     &       RTDB_ERR)
c
c We got the tag_to_add info from the geometry object
c
        do iunique = 1, uniquecent
           if (.not. geom_cent_tag(geom,
     ,          int_mb(i_array+iunique-1),tag_to_add))
     &        call errquit('bas_rtdb_do_load: geom_cen_tag',211,
     &       RTDB_ERR)
c
c First check if we did this geometry tag already
c
           do junique = 1, iunique - 1
              if (.not. geom_cent_tag(geom,
     ,             int_mb(i_array+junique-1),tag_in_lib))
     &           call errquit('bas_rtdb_do_load: geom_cen_tag',211,
     &       RTDB_ERR)
              if (inp_compare(.true.,tag_to_add,tag_in_lib)
     &              .and. (inp_strlen(tag_to_add) .eq. 
     &                     inp_strlen(tag_in_lib))) goto 00012
           enddo
c
c
c Now for each star tag, check if tag_to_add matches and add basis set
c
           do istar = 1, star_nr_tags
c
c For * input, directly check exceptions list
c For aa*, first check if tag_to_add matches star tag itself
c For Bq, if B* is used and geometry tag is Bq, skip, only when Bq* 
c is used we need to go on
c Skip and do not load a basis for dummy center X
c
              if (inp_contains(.true.,'*',star_tag(istar),ilen)) then
                 if (ilen .gt. 1) then
                    if (tag_to_add(1:ilen-1) .ne. 
     &                 star_tag(istar)(1:ilen-1)) goto 00011
                 endif
                 if (inp_compare(.false.,'Bq',tag_to_add(1:2)).and..not.
     &               inp_compare(.false.,'Bq',star_tag(istar)(1:2))) 
     &               goto 00011
                 if (inp_compare(.false.,'X',tag_to_add(1:1)).and..not.
     &               inp_compare(.false.,'Xe',tag_to_add(1:2))) 
     &               goto 00011
              endif
c
c There is a match between the tag_to_add and the star_tag, check for
c matches in exceptions list
c
              oIsexcept = .false.
              istart = 1
              if (istar .gt. 1) istart = star_nr_excpt(istar-1) + 1
              do iexcpt = istart, star_nr_excpt(istar)
                 ilen = inp_strlen(tag_to_add)
                 if (inp_compare(.true.,tag_to_add,star_excpt(iexcpt))
     &              .and. (inp_strlen(tag_to_add) .eq. 
     &                     inp_strlen(star_excpt(iexcpt))))
     &              oIsexcept = .true.
              enddo
c
c If not in exceptions list, add basis set to the basis object
c
              if (.not. oIsexcept) then
                 tag_in_lib = star_in_lib(istar)
                 if (inp_contains(.true.,'*',star_in_lib(istar),ilen))
     &              tag_in_lib = tag_to_add
                 call bas_tag_lib(basisin,star_segment,tag_to_add, 
     &                tag_in_lib, star_bas_typ(istar), star_file(istar), 
     &                star_rel(istar))
              endif
00011         continue
           enddo
00012   continue
        enddo
        if (.not.MA_Free_Heap(l_array)) call
     E       errquit(' basrtdbdoload: failed freeheap',0, MA_ERR)
c
c We messed around with the memory. The basis set sould have room for 
c two (2) extra variables, which are added below. Add them now.
c
        h_tmp = exndcf(H_exndcf,basis)
        k_tmp = exndcf(K_exndcf,basis) 
        ilen = exndcf(SZ_exndcf,basis)
        if (.not. ma_alloc_get(mt_dbl,ilen+2,nameex,h_new,k_new)) call
     &     errquit('bas_rtdb_do_load: ma_alloc_get failed',ilen+2,
     &       MA_ERR)
        exndcf(H_exndcf,basis) = h_new
        exndcf(K_exndcf,basis) = k_new
        exndcf(SZ_exndcf,basis) = ilen + 2
        call dcopy(ilen,dbl_mb(k_tmp),1,dbl_mb(k_new),1)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &      ('bas_rtdb_do_load: error freeing old exponents',211,
     &       MEM_ERR)
c
      endif
c
      star_nr_tags = 0
c
c read the basis now get check status of read operations
c
      if (.not.rtdb_status) then
         if (exndcf(H_exndcf,basis) .ne. -1) then
            if (.not. ma_free_heap(exndcf(H_exndcf,basis)))
     $           call errquit('bas_rtdb_load: ma free failed?',0,
     &       MA_ERR)
            exndcf(H_exndcf,basis) = -1
            exndcf(K_exndcf,basis) = 0
            exndcf(SZ_exndcf,basis)= 0
         endif
c
c     rjh ... can be quiet now since the application should
c     whine if it really did need the basis set
c
*         if(ga_nodeid().eq.0) then
*            write(LuOut,*)
*            write(LuOut,*) ' bas_rtdb_do_load: basis not present "',
*     $           bs_name(basis)(1:inp_strlen(bs_name(basis))),
*     &          '" -> "',
*     $           bs_trans(basis)(1:inp_strlen(bs_trans(basis))),
*     &          '"'
*            write(LuOut,*)
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
c Add zero exponent S function for sp code/texas interface, 
*     incore information only
c
      if (bas_rtdb_do_load) then
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
*
        call bas_ecce_print_basis(basisin,'load_basis')
      endif
      end
*.....................................................................
      logical function bas_geobas_build(basisin)
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "global.fh"
#include "bas_ibs_dec.fh"
#include "ecpso_decP.fh"
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
      integer my_gen, my_type
      integer atn
      character*16 element
      character*2  symbol
      logical status
      logical foundit
      logical found_any
      logical is_bq
      logical ecpORso, relbas
      integer int_dummy, num_elec
      integer h_tmp, k_tmp
      character*2 tag12
      character*16 name_tmp
      double precision erep_save
      integer uce, idbstag
*debug:mem      integer inode
c::functions
      logical bas_high_angular
      external bas_high_angular
      integer nbf_from_ucont
      external nbf_from_ucont
      logical ecp_get_num_elec
      external ecp_get_num_elec
      logical bas_match_tags
      external bas_match_tags
      logical basis_is_rel
      external basis_is_rel
c
#include "bas_ibs_sfn.fh"
#include "ecpso_sfnP.fh"
c
      basis = basisin + Basis_Handle_Offset
      geom  = ibs_geom(basis)
      ecpORso = Is_ECP(basis).or.Is_SO(basis)
      relbas = basis_is_rel(basisin)
c
*debug:mem      do inode = 0,(ga_nnodes() - 1)
*debug:mem        if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem        call ga_sync()
*debug:mem      enddo
c

      status = geom_ncent(geom, nat)
      if (nat.eq.0.or..not.status) then
        write(LuOut,*)' bas_geobas_build: ERROR '
        write(LuOut,*)' number of centers is zero or weird'
        write(LuOut,*)' nat = ',nat
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
c.... set flag if any general contractions are present
      bas_any_gc(basis)       = .false.
      if(.not. ecpORso) then
        do i = 1,(infbs_head(Head_Ncont,basis))
          if (.not.bas_any_gc(basis)) then
            my_gen  = infbs_cont(Cont_Ngen,i,basis)
            my_type = infbs_cont(Cont_Type,i,basis)
            if (my_gen.gt.1.and.my_type.ge.0)
     &          bas_any_gc(basis) = .true.
          endif
        enddo
      endif
c.... set flag if any sp (spd,spdf) shells are present
      bas_any_sp_shell(basis) = .false.
      if(.not. ecpORso) then
        do i = 1,(infbs_head(Head_Ncont,basis))
          if (.not.bas_any_sp_shell(basis)) then
            my_gen  = infbs_cont(Cont_Ngen,i,basis)
            my_type = infbs_cont(Cont_Type,i,basis)
            if (my_type.lt.0) then
              if (my_gen.ne.2) then
                write(luout,*)' sp shell with n_gen = ',my_gen
                call errquit ('bas_geobas_build: fatal error',911,
     &       BASIS_ERR)
              endif
              bas_any_sp_shell(basis) = .true.
            endif
          endif
        enddo
      endif
c
c... clear old ibs_ce2uce if it exists
      if (ibs_ce2uce(SZ_ibs,basis).gt.0) then
*debug:mem        write(LuOut,*)' clearing old ce2uce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_ce2uce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_ce2uce',911,
     &       BASIS_ERR)
        ibs_ce2uce(H_ibs,basis)  = 0
        ibs_ce2uce(K_ibs,basis)  = 0
        ibs_ce2uce(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_ce2uce ',basis
      if (.not.ma_alloc_get(mt_int,nat,name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_ce2uce',911,
     &       MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating ce2uce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill(nat,0,int_mb(k_tmp),1)
        ibs_ce2uce(H_ibs,basis)  = h_tmp
        ibs_ce2uce(K_ibs,basis)  = k_tmp
        ibs_ce2uce(SZ_ibs,basis) = nat
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem          call ga_sync()
*debug:mem          call util_flush(LuOut)
*debug:mem        enddo
      endif
c
c build center to unique center map
c
      if (nat.gt.nat_mx) then
        write(LuOut,*)' nat     = ',nat
        write(LuOut,*)' nat max = ',nat_mx
        call errquit ('bas_geobas_build: nat.gt.nat_mx',911, BASIS_ERR)
      endif
      found_any = .false.
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
          found_any = .true.
          goto 00102
        endif
        if (.not. foundit .and. .not. (ecpORso .or. relbas)) then
          if (geom_tag_to_element(tags(i,geom), symbol, element,
     $        atn)) then
            if (ga_nodeid().eq.0)
     &          write(LuOut,10) 
     &          tags(i,geom)(1:inp_strlen(tags(i,geom))),
     $          element(1:inp_strlen(element)),
     $          bs_name(basis)(1:len_bs_name(basis))
 10         format(/' ERROR: geometry tag ',a,' (',a,
     &          ') is an atom ',
     $          'but has no functions in basis "',a,'"'/
     $          ' ERROR: only bq* centers can have no functions')
            if (ga_nodeid().eq.0) call util_flush(LuOut)
            call errquit
     &          ('bas_geobas_build: basis/geometry mismatch', 0,
     &       BASIS_ERR)
          else
            tag12 = tags(i,geom)(1:2)
            is_bq = inp_compare(.false.,'bq',tag12) .or.
     $           (inp_compare(.false.,'X',tags(i,geom)(1:1)) .and. 
     $           (.not. inp_compare(.false.,'e',tags(i,geom)(2:2))))
            if (ga_nodeid().eq.0 .and.(.not.is_bq))
     &          write(LuOut,11) i,
     &          tags(i,geom)(1:inp_strlen(tags(i,geom))),
     $          bs_name(basis)(1:len_bs_name(basis))
 11         format(/'WARNING: geometry tag ',i4, ' ', a,
     $          ' not found in basis "',a,'"'/)
            int_mb(mb_ibs_ce2uce(i,basis)) = 0
            if (ga_nodeid().eq.0) call util_flush(LuOut)
          endif
        endif
00102   continue
00100 continue
*
      if (.not.found_any) then
        if (ecpORso) then
          bas_geobas_build = .false.
          return
        else
          if (ga_nodeid().eq.0) then
            write(luout,*)' none of the geometry tags matched any ',
     &          'basis set tag in the basis "',
     &          bs_name(basis)(1:len_bs_name(basis)),
     &          '"'
          endif
          call errquit('bas_geobas_build: fatal error',911,
     &       BASIS_ERR)
        endif
      endif
c
c build total # of contractions 
c      
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) then
*debug:mem            call MA_summarize_allocated_blocks()
*debug:mem            if (MA_verify_allocator_stuff()) then
*debug:mem              write(LuOut,*)' no errors'
*debug:mem            else
*debug:mem              write(LuOut,*)' errors'
*debug:mem            endif
*debug:mem          endif
*debug:mem          call ga_sync()
*debug:mem          call util_flush(LuOut)
*debug:mem        enddo
*debug:mem      call bas_print_allocated_info('bas_geobas_build 1')
      ncont_tot_gb(basis)  = 0
      do 10200 i=1,nat
         uce = sf_ibs_ce2uce(i,basis)
*debug:mem         write(LuOut,*)' myuce = ',uce,i,ga_nodeid()
*debug:mem         call util_flush(LuOut)
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
*debug:mem        write(LuOut,*)' clearing old ce2cnr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_ce2cnr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_ce2cnr',911,
     &       MEM_ERR)
        ibs_ce2cnr(H_ibs,basis)  = 0
        ibs_ce2cnr(K_ibs,basis)  = 0
        ibs_ce2cnr(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_ce2cnr ',basis
      if (.not.ma_alloc_get(mt_int,(2*nat),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_ce2cnr',911, MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating ce2cnr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((2*nat),0,int_mb(k_tmp),1)
        ibs_ce2cnr(H_ibs,basis)  = h_tmp
        ibs_ce2cnr(K_ibs,basis)  = k_tmp
        ibs_ce2cnr(SZ_ibs,basis) = 2*nat
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem          call ga_sync()
*debug:mem          call util_flush(LuOut)
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
     $    ('bas_geobas_build: no functions in basis set', 0,
     &       BASIS_ERR)
      if (ncont_tot_gb(basis) .gt. ncont_mx) then
        write(LuOut,*)' number of contractions     = ',
     &      ncont_tot_gb(basis)
        write(LuOut,*)' number of contractions max = ',ncont_mx
        call errquit ('bas_geobas_build: ncont.gt.ncont_mx ',911,
     &       BASIS_ERR)
      endif
c
c allocate space for contraction -> center map
c
c... clear old ibs_cn2ce if it exists
      if (ibs_cn2ce(SZ_ibs,basis).gt.0) then
*debug:mem        write(LuOut,*)' clearing old cn2ce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_cn2ce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2ce',911,
     &       MEM_ERR)
        ibs_cn2ce(H_ibs,basis)  = 0
        ibs_cn2ce(K_ibs,basis)  = 0
        ibs_cn2ce(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2ce  ',basis
      if (.not.ma_alloc_get(mt_int,(1+ncont_tot_gb(basis)),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit('bas_geobas_build: error ma_alloc ibs_cn2ce',911,
     &       MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating cn2ce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((1+ncont_tot_gb(basis)),0,int_mb(k_tmp),1)
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
*debug:mem        write(LuOut,*)' clearing old cn2ucn data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_cn2ucn(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2ucn',911,
     &       BASIS_ERR)
        ibs_cn2ucn(H_ibs,basis)  = 0
        ibs_cn2ucn(K_ibs,basis)  = 0
        ibs_cn2ucn(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2ucn ',basis
      if (.not.ma_alloc_get(mt_int,(1+ncont_tot_gb(basis)),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_cn2ucn',911,
     &       MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating cn2ucn data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((1+ncont_tot_gb(basis)),0,int_mb(k_tmp),1)
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
               write(LuOut,*)' bas_geobas_build: ERROR '
               write(LuOut,*)' contraction range size mismatch'
               write(LuOut,*)'        cont. range (',jstart,':',jend,')'
               write(LuOut,*)' unique cont. range (',kstart,':',kend,')'
               write(LuOut,*)'        cont. size: ',jsize
               write(LuOut,*)' calculated unique cont. size: ',ksize
               write(LuOut,*)'     lookup unique cont. size: ',lsize
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
*debug:mem        write(LuOut,*)' clearing old cn2bfr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_cn2bfr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2bfr',911,
     &       BASIS_ERR)
        ibs_cn2bfr(H_ibs,basis)  = 0
        ibs_cn2bfr(K_ibs,basis)  = 0
        ibs_cn2bfr(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2bfr ',basis
      if (.not.ma_alloc_get(mt_int,(2*(1+ncont_tot_gb(basis))),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_cn2bfr',911,
     &       BASIS_ERR)
      else
*debug:mem        write(LuOut,*)' generating cn2bfr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((2*(1+ncont_tot_gb(basis))),0,int_mb(k_tmp),1)
        ibs_cn2bfr(H_ibs,basis)  = h_tmp
        ibs_cn2bfr(K_ibs,basis)  = k_tmp
        ibs_cn2bfr(SZ_ibs,basis) = 2*(1+ncont_tot_gb(basis))
      endif
c
c build nprim_tot_gb, ncoef_tot_gb, nbf_tot_gb, and 
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
        ncoef_tot_gb(basis) = ncoef_tot_gb(basis) + 
     &         infbs_cont(CONT_NPRIM,iu_cont,basis)*
     &         infbs_cont(CONT_NGEN,iu_cont,basis)
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
      if (.not.bas_high_angular(basisin,myang))call errquit
     &      ('bas_geobas_build: error bas_high_angular',911,
     &       BASIS_ERR)
*
      if (Is_ECP(basis)) then
* if ECP basis must modify geometry data appropriately
        do i = 1,nat
          oecpcent(i,geom)=.false.
          idbstag = sf_ibs_ce2uce(i,basis)
          if (ecp_get_num_elec(basisin,
     &        bs_tags(idbstag,basis),num_elec)) then
            oecpcent(i,geom) = .true.
            charge(i,geom) = charge(i,geom) - dble(num_elec)
          endif
        enddo
        erep_save = erep(geom)
        call geom_compute_values(geom)
*
*     cannot print this here
*
*        write(luout,*)
*     &      ' nuclear repulsion energy without ECP: ',erep_save
*        write(luout,*)
*     &      ' nuclear repulsion energy with    ECP: ',erep(geom)
      endif
*
      bas_geobas_build = .true.
*
      end
*.....................................................................
C> \ingroup bas
C> @{
c
C> \brief Store a basis set on the RTDB under the given name
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function bas_rtdb_store(rtdb, name, basisin)
      implicit none
c
c routine that does an incore basis set store to the rtdb
c
#include "mafdecls.fh"
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "ecpso_decP.fh"
#include "bas_starP.fh"
c::passed
      integer rtdb       !< [Input] the RTDB handle
      character*(*) name !< [Input] the name to use when storing basis
      integer basisin    !< [Input] the basis set handle
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
c
cc AJL/Begin/SPIN-POLARISED ECPs
      integer channels
      logical ecp_get_high_chan
      external ecp_get_high_chan
cc AJL/End
c
c:: statement functions
#include "errquit.fh"
#include "bas_exndcf.fh"      
#include "ecpso_sfnP.fh"
c
cc AJL/Begin/SPIN-POLARISED ECPs
      channels = 1
cc AJL/End
c
      bas_rtdb_store = bas_check_handle(basisin,'bas_rtdb_store')
      if (.not. bas_rtdb_store) return
      basis = basisin + Basis_Handle_Offset
c
      if (Is_ECP(basis).or.Is_SO(basis)) then
        size_ex = (2*infbs_head(HEAD_NPRIM,basis))+
     &      infbs_head(HEAD_NCOEF,basis)
c
cc AJL/Begin/SPIN-POLARISED ECPs
cc Use this to put a note of spin_polarised_ecps in rtdb
cc For use in nwdft/dft_fockbld.F
cc
cc By using this function we deal with anomalies for all
cc applications of the spin-polarised ECPs being for both channels,
cc as ecp_get_high_chan will treat that as spin-independent
cc calculation of the fock matrices
        if (Is_ECP(basis)) then
          if (.not.ecp_get_high_chan(basisin,channels))
     &      call errquit('bas_rtdb_store error',911, BASIS_ERR)
        endif
cc AJL/End
c
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
     $    size_ex,
     &    bs_stdname(1,basis),
     &    name_assoc(1,basis),
     &    channels,  ! AJL/SPIN-POLARISED ECPs
     &    star_nr_tags, star_tag, star_in_lib, star_bas_typ,
     &    star_file, star_nr_excpt, star_excpt, star_tot_excpt,
     &    star_rel, star_segment)
c
      end
C> @}
*.....................................................................
      logical function bas_rtdb_do_store(rtdb, name, tagsin, 
     &       head_array, tags_array, ucont_array, excfin, ntagsin, 
     &       nucontin, nexcf, stdnamein, ecp_name, ecp_channels,
     &       star_nr_tags, star_tag, star_in_lib, star_bas_typ,
     &       star_file, star_nr_excpt, star_excpt, star_tot_excpt,
     &       star_rel, star_segment)
      implicit none
c
c stores from argument data structures the basis set information to
c the rtdb. 
c
#include "mafdecls.fh"
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
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
      character*80 stdnamein(ntagsin) ! [input] names of basis on a tag
      integer head_array(ndbs_head) ! [input] head data
      integer tags_array(ndbs_tags,ntagsin) ! [input] tag data
      integer ucont_array(ndbs_ucont,nucontin) ! [input] unique 
*. . . . . . . . . . . . . . . . . . . . . . .   contraction data
      character*(*) ecp_name(2)                 ! [input] ecp 
*. . . . . . . . . . . . . . . . . . . . . . .   associated with 
*. . . . . . . . . . . . . . . . . . . . . . .   normal basis
      double precision excfin(nexcf) ! [input] exponents 
c. . . . . . . . . . . . . . . . . . .         contractions coeffs.
      integer star_nr_tags ! [input] number of star containing tags
      character*16 star_tag(star_nr_tags) ! [input] star tag data
      character*16 star_in_lib(star_nr_tags) ! [input] star in lib data
      character*255 star_bas_typ(star_nr_tags) ! [input] star basis data
      character*255 star_file(star_nr_tags) ! [input] star filename data
      logical star_rel(star_nr_tags) ! [input] star relativistic data
      logical star_segment ! [input] star segment or not ?
      integer star_tot_excpt ! [input] star total number of except tags
      integer star_nr_excpt(star_nr_tags) ! [input] number of except per tag
      character*16 star_excpt(star_nr_tags) ! [input] except tag list
c
cc AJL/Begin/SPIN-POLARISED ECPs
      integer ecp_channels ! [input] number of channels for ECP (1 or 2).
c............................        *Should* be 1 in all other scenarios
cc AJL/End
c
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
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':bs_nr_tags'
      status = status .and. rtdb_put(rtdb,tmp,mt_int,1,ntagsin)
      if (ntagsin .gt. 0) then
         tmp(lentmp:) = ' ' 
         tmp(lentmp:) = ':bs_tags'
         status = status .and. rtdb_cput(rtdb,
     &            tmp,ntagsin,tagsin)

         tmp(lentmp:) = ' ' 
         tmp(lentmp:) = ':bs_stdname'
         status = status .and. rtdb_cput(rtdb,
     &            tmp,ntagsin,stdnamein)
      endif
      
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':assoc ecp name'
      status = status .and. rtdb_cput(rtdb,tmp,2,ecp_name)
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':number of exps and coeffs'
      status = status .and. rtdb_put(rtdb,tmp,mt_int,1,nexcf)

      if (nexcf .gt. 0) then
         tmp(lentmp:) = ' ' 
         tmp(lentmp:) = ':exps and coeffs'
         status = status .and. rtdb_put(rtdb,tmp,mt_dbl,nexcf,excfin)
      endif

      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':header'
      status = status .and.
     &       rtdb_put(rtdb,tmp,mt_int,ndbs_head,head_array)

      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':tags info'
      status = status .and. rtdb_put(
     &       rtdb,tmp,mt_int,(ndbs_tags*ntags_bsmx), tags_array)

      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':contraction info'
      status = status .and. rtdb_put(
     &       rtdb,tmp,mt_int,(ndbs_ucont*nucont_bsmx),ucont_array)
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':star nr tags'
      status = status .and. rtdb_put(
     &       rtdb,tmp,mt_int,1,star_nr_tags)
c
c store remaining star tag info only if we have actually a star 
c tag input (hence, if star_nr_tags > 0)
c
      if (star_nr_tags .gt. 0) then
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star tag names'
        status = status .and. rtdb_cput(rtdb,tmp,star_nr_tags,
     &                                  star_tag)
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star tag_in_lib'
        status = status .and. rtdb_cput(rtdb,tmp,star_nr_tags,
     &                                  star_in_lib)
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star bas type'
        status = status .and. rtdb_cput(rtdb,tmp,star_nr_tags,
     &                                  star_bas_typ)
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star filename'
        status = status .and. rtdb_cput(rtdb,tmp,star_nr_tags,
     &                                  star_file)
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star tot except'
        status = status .and. rtdb_put(
     &         rtdb,tmp,mt_int,1,star_tot_excpt)
        if ( star_tot_excpt .gt. 0) then
           tmp(lentmp:) = ' ' 
           tmp(lentmp:) = ':star except'
           status = status .and. rtdb_cput(rtdb,tmp,star_tot_excpt,
     &                                     star_excpt)
        endif
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star nr except'
        status = status .and. rtdb_put(
     &         rtdb,tmp,mt_int,star_nr_tags,star_nr_excpt)
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star rel'
        status = status .and. rtdb_put(
     &         rtdb,tmp,mt_log,star_nr_tags,star_rel)
        tmp(lentmp:) = ' ' 
        tmp(lentmp:) = ':star segment'
        status = status .and. rtdb_put(
     &         rtdb,tmp,mt_log,1,star_segment)
      endif
c
c reset the stag tag data array for reuse
c
      star_nr_tags = 0
c
cc AJL/Begin/SPIN-POLARISED ECPs
cc Add if we are using spin_polarised_ecps
cc There may be a better place for this information
      if (ecp_channels.gt.1) then
        status = status .and.  rtdb_put(
     &    rtdb,'dft:spin_polarised_ecps',mt_int,1,ecp_channels)
      endif
cc AJL/End
c
c read the basis now get check status of read operations
c
      if (.not.status) then
        write(LuOut,*)' bas_rtdb_store: ERROR '
        write(LuOut,*)' one or more put operations failed '
        bas_rtdb_do_store = .false.
c..... add diagnostics later
        return
      endif
      return
      end
*.....................................................................
C> \ingroup bas
C> @{
c
C> \brief Retrieve the highest angular momentum in the specified basis
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function bas_high_angular(basisin,high_angular)
      implicit none
c
c  calculate, return and store high angular momentem function
c   for given basis. 
c
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "stdio.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin      !< [Input] the basis set handle
      integer high_angular !< [Output] the maximum angular momentum
c:local
      integer basis, i, myang
      integer myutag
c
      bas_high_angular = bas_check_handle(basisin,'bas_high_angular')
      if (.not. bas_high_angular ) then
        write(luout,*) 'bas_high_angular: basis handle not valid '
        return
      endif
c
      basis = basisin + Basis_Handle_Offset
      if (angular_bs(basis) .gt. -565) then
        high_angular = angular_bs(basis)
        bas_high_angular = .true.
        return
      endif
c
      myutag =  infbs_head(head_ntags,basis)
      high_angular = -565
      do i = 1,myutag
        myang = abs(infbs_tags(tag_high_ang,i,basis))
        high_angular = max(high_angular,myang)
      enddo
*is this needed here?
      angular_bs(basis) = high_angular
c
      bas_high_angular = .true.
      return
      end
C> @}
*.....................................................................
cc AJL/Begin
C> \ingroup ecp
C> @{
c
C> \brief Retrieve the number of channels in the ecp basis
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function ecp_get_high_chan(ecpidin,channels)
      implicit none
c
c  calculate and return high channel for given ecp. 
c
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "stdio.fh"
c::functions
      logical ecp_check_handle
      external ecp_check_handle
c:: passed
      integer ecpidin      !< [Input] the ECP set handle
      integer channels     !< [Output] the number of channels
c:local
      integer ecpid, i, j, mychannel
      integer myucont, mytags
c 
      ecpid = ecpidin + Basis_Handle_Offset
c
      ecp_get_high_chan = ecp_check_handle(ecpidin,'ecp_get_high_chan')
      if (.not. ecp_get_high_chan) then
        write (6,*) 'This is not an ECP Basis'
        return
      endif
c
      mytags  = infbs_head(HEAD_NTAGS,ecpid)
      if (mytags.le.0) return 
c
      channels = 0
      do i=1,mytags
        myucont = infbs_tags(TAG_NCONT,i,ecpid)
        do j = 1,myucont
          mychannel = infbs_cont(CONT_CHANNEL,j,ecpid) 
          channels = max(channels,mychannel)
        enddo
      enddo
c If all values are 0 we have a spin independent calculation 
      if (channels.eq.0) then
        channels = 1
c Somewhere we have initialised spin in the alpha of beta channels
      else if (channels.gt.0) then 
        channels = 2
      end if

c Save value process to be added
c
      return
      end
C> @}
cc AJL/End
*.....................................................................
      logical function gbs_map_clear(basisin)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
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
        write(LuOut,*) ' basis handle not valid '
        return
      endif
c
      basis = basisin + Basis_Handle_Offset
c... clear old ibs_ce2uce if it exists
      if (ibs_ce2uce(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_ce2uce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_ce2uce',911, MEM_ERR)
        ibs_ce2uce(H_ibs,basis)  = 0
        ibs_ce2uce(K_ibs,basis)  = 0
        ibs_ce2uce(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_ce2cnr if it exists
      if (ibs_ce2cnr(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_ce2cnr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_ce2cnr',911, MEM_ERR)
        ibs_ce2cnr(H_ibs,basis)  = 0
        ibs_ce2cnr(K_ibs,basis)  = 0
        ibs_ce2cnr(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_cn2ce if it exists
      if (ibs_cn2ce(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_cn2ce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_cn2ce',911, MEM_ERR)
        ibs_cn2ce(H_ibs,basis)  = 0
        ibs_cn2ce(K_ibs,basis)  = 0
        ibs_cn2ce(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_cn2ucn if it exists
      if (ibs_cn2ucn(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_cn2ucn(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_cn2ucn',911, MEM_ERR)
        ibs_cn2ucn(H_ibs,basis)  = 0
        ibs_cn2ucn(K_ibs,basis)  = 0
        ibs_cn2ucn(SZ_ibs,basis) = 0 
      endif
c... clear old ibs_cn2bfr if it exists
      if (ibs_cn2bfr(SZ_ibs,basis).gt.0) then
        h_tmp = ibs_cn2bfr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('gbs_map_clear: error freeing ibs_cn2bfr',911, MEM_ERR)
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
      ncoef_tot_gb(basis) = 0
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
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
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
      basis = basisin + Basis_Handle_Offset
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
        write(LuOut,*)' gbs_map_print: ERROR '
        write(LuOut,*)' number of centers is zero or weird'
        write(LuOut,*)' nat = ',nat
        gbs_map_print = .false.
c..... add diagnostics later
        return
      endif
c
c print global information
c
      write(LuOut,*)'<<< GBS_MAP_PRINT >>>'
      write(LuOut,*)' total number of atoms           :',nat
      write(LuOut,*)' total number of contractions    :',
     &       ncont_tot_gb(basis)
      write(LuOut,*)' total number of primitives      :',
     &       nprim_tot_gb(basis)
      write(LuOut,*)' total number of coefficients    :',
     &       ncoef_tot_gb(basis)
      write(LuOut,*)' total number of basis functions :',
     &       nbf_tot_gb(basis)
c
c print center based mapping information.  
c
      write(LuOut,*)' '
      write(LuOut,*)'=================================================',
     &       '==============================='
      write(LuOut,*)' center -> unique center map          <ibs_ce2uce>'
      write(LuOut,*)'        -> contraction range map      <ibs_ce2cnr>'
      write(LuOut,*)'=================================================',
     &       '==============================='
      do 00100 i=1,nat
        write(LuOut,'(1x,a,i4,2x,a,i3)')
     &         'center:',i,'maps to unique center:',
     &        sf_ibs_ce2uce(i,basis)
        myfirst = sf_ibs_ce2cnr(1,i,basis)
        mylast  = sf_ibs_ce2cnr(2,i,basis)
        mysize  = mylast - myfirst + 1
        write(LuOut,'(14x,a,i4,2x,a,i4,a,i4,a,/)')
     &         'has',mysize,'contractions   <first:',
     &         myfirst,'>  <last:',mylast,'>'
00100 continue
c
c print contraction based mapping information
c
      write(LuOut,*)' '
      write(LuOut,*)'=================================================',
     &       '==============================='
      write(LuOut,*)' contraction -> center map                     ',
     &       '                <ibs_cn2ce>'
      write(LuOut,*)'             -> unique contraction in basis set',
     &       '                <ibs_cn2ucn>'
      write(LuOut,*)'             -> basis function range           ',
     &       '                <ibs_cn2bfr>'
      write(LuOut,*)'=================================================',
     &       '==============================='
c
      do 00200 i=1,ncont_tot_gb(basis)
        mycenter = sf_ibs_cn2ce(i,basis)
        myucont  = sf_ibs_cn2ucn(i,basis)
        myfirst  = sf_ibs_cn2bfr(1,i,basis)
        mylast   = sf_ibs_cn2bfr(2,i,basis)
        mysize   = mylast - myfirst + 1
        write(LuOut,'(1x,a,i4,2x,a,i3)')
     &         'contraction',i,'is on center:',mycenter
        write(LuOut,'(18x,a,i3)')
     &         'is represented by unique contraction:',myucont
        write(LuOut,'(18x,a,i5,a,i5,a,i5,a,/)')
     &         'has',mysize,' basis functions <first:',myfirst,
     &         '>  <last:',mylast,'>'
00200 continue
c
      gbs_map_print = .true.
      return
      end
*.....................................................................
C> \ingroup bas
C> @{
c
C> \brief Store the names of all current basis set instances on the RTDB
c
C> Write the number of current basis set instances and the names
C> of all basis set instances to the RTDB. No actual contents of the
C> basis sets is stored by this function. The bas_rtdb_store routine
C> can be used to store the actual basis set data instead.
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function bas_rtdb_out(rtdb)
c     
c     output to rtdb info about known basis sets
c     
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "rtdb.fh"
#include "stdio.fh"
c     
c::passed
      integer rtdb !< [Input] the RTDB handle
c
      bas_rtdb_out  =
     &     rtdb_put(rtdb, 'basis:nbasis',
     &       MT_INT, 1, nbasis_rtdb)
     &     .and.
     &     rtdb_cput(rtdb, 'basis:names', nbasis_rtdb,
     &       bs_names_rtdb)
      if (.not. bas_rtdb_out) 
     &     write(LuOut,*) ' bas_rtdb_out: rtdb is corrupt '
c     
      end
*.....................................................................
c
C> \brief Add another basis set name to the list of known basis sets on 
C> the RTDB
c
C> Just add another name to the list of known basis set instances. No
C> actual basis set data is stored.
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function bas_rtdb_add(rtdb, name)
      implicit none
c
c add basis set name to known basis set list on rtdb
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "stdio.fh"
c     
      integer rtdb       !< [Input] the RTDB handle
      character*(*) name !< [Input] the name of basis set to add
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
         write(LuOut,*) ' bas_rtdb_add: too many basis tries on rtdb ',
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
         write(LuOut,*) ' bas_rtdb_add: rtdb error adding ', name(1:ln)
         return
      endif
c     
      bas_rtdb_add = .true.
c     
      end
*.....................................................................
c
C> \brief Print contents of all current basis set instances
c
C> \return Return .true. if successful, and .false. otherwise
c
      logical function bas_print_all()
c
c routine to print active all basis set(s) information
c
      implicit none
#include "nwc_const.fh"
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
          basin = basis - Basis_Handle_Offset
          bas_print_all = bas_print_all .and. bas_print(basin)
        endif
00100 continue
c
      return
      end
c
C> \brief Print the names of all basis set instances on the RTDB
c
C> \return Return .true. if successful, and .false. otherwise
c
      subroutine bas_print_known_bases(rtdb)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb !< [Input] the RTDB handle
c     
c     Print list of known basis sets
c     
      logical bas_rtdb_in
      external bas_rtdb_in
      logical ignore
      integer basis, ma_type, natom, nelem
      character*26 date
      character*32 name32
      character*128 key
      integer head_array(ndbs_head)
c
      ignore = bas_rtdb_in(rtdb)
      if (ga_nodeid() .eq. 0) then
         write(LuOut,*)
         call util_print_centered(LuOut,'Basis sets in the database',
     $        23,.true.)
         write(LuOut,*)
         if (nbasis_rtdb .le. 0) then
            write(LuOut,*) ' There are no basis sets in the database'
            write(LuOut,*)
         else
            if (nbasis_rtdb .gt. 0) write(LuOut,3)
 3          format(
     $           1x,4x,2x,'Name',28x,2x,'Natoms',2x,
     $           'Last Modified',/,
     $           1x,5x,2x,32('-'),2x,6('-'),2x,24('-'))
            do basis = 1, nbasis_rtdb
               key = ' '
               write(key,'(''basis:'',a,'':header'')')
     $              bs_names_rtdb(basis)(1:len_bs_rtdb(basis))
               if (.not. rtdb_get(rtdb, key, mt_int, 
     $              ndbs_head, head_array)) then
                  write(LuOut,*) ' Warning: basis set ', basis, 
     $                 ' may be corrupt'
                  natom = -1
               else
                  natom = head_array(HEAD_NTAGS)
               endif
               if (.not. rtdb_get_info(rtdb, key, ma_type, 
     $              nelem, date)) then
                  write(LuOut,*) ' Warning: basis set ', basis, 
     $                 ' may be corrupt'
                  date = 'unknown'
               endif
               name32 = bs_names_rtdb(basis)(1:len_bs_rtdb(basis))
               write(LuOut,4) basis, name32, natom, date
 4             format(1x,i4,2x,a32,2x,i6,2x,a26)
            end do
            if (nbasis_rtdb .gt. 0) then
               if (.not. rtdb_cget(rtdb,'ao basis',1,key)) 
     $              key = 'ao basis'
               write(LuOut,*)
               write(LuOut,5) key(1:inp_strlen(key))
 5             format(2x,'The basis set named "',a,
     $              '" is the default AO basis for restart')
            endif
            write(LuOut,*)
            write(LuOut,*)
         endif
         call util_flush(LuOut)
      endif
c     
      end
C> @}
*.....................................................................
      logical function bas_rtdb_in(rtdb)
c     
c     load in info about known basis sets ... this is more
c     for diagnostic and debugging purposes
c     
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "stdio.fh"
c     
      integer rtdb     ! [input] run time data base handle
      integer bas
      bas_rtdb_in = .false.
      nbasis_rtdb = 0
      if (rtdb_get(rtdb, 'basis:nbasis', MT_INT, 1, nbasis_rtdb))
     $     then
        if (.not. rtdb_cget(rtdb,'basis:names', nbasis_rtdb_mx,
     $        bs_names_rtdb)) then
          write(LuOut,*) 'bas_rtdb_in: rtdb corrupt'
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
C> \ingroup bas
C> @{
c
C> \brief Retrieve the center rank of a given shell in a basis set
c
C> After the basis set has been mapped to a geometry a complete list
C> of shells (or contractions) is available. This routine retrieves
C> the center rank of a shell from the complete list of shells.
c
C> \return Return .true. if successful, and .false. otherwise.
c
      logical function bas_cn2ce(basisin,cont,center)
      implicit none
c
c returns the center for a given mapped contraction
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle
      integer cont    !< [Input] the mapped contraction index
      integer center  !< [Output] the center rank
c::local
      integer basis
#include "bas_ibs_sfn.fh"
c
      bas_cn2ce = bas_check_handle(basisin,'bas_cn2ce')
      if(.not.bas_cn2ce) return
c
      basis = basisin + Basis_Handle_Offset
      bas_cn2ce = cont.ge.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2ce) then
        write(LuOut,*)' bas_cn2ce: invalid contraction information '
        write(LuOut,*)' contraction range is 0:',ncont_tot_gb(basis)
        write(LuOut,*)' input contraction was: ',cont
        return
      endif
      center = sf_ibs_cn2ce(cont,basis)
c
      return
      end
*.....................................................................
c
C> \brief Retrieve the rank of the symmetry unique center of a shell in
C> a basis set
c
C> After the basis set has been mapped to a geometry a complete list
C> of shells (or contractions) is available. This routine retrieves
C> the center rank of the symmetry unique center of a shell from the
C> complete list of shells.
c
C> \return Return .true. if successful, and .false. otherwise.
c
      logical function bas_cn2uce(basisin,cont,ucenter)
      implicit none
c
c returns the UNIQUE center for a given mapped contraction
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle
      integer cont    !< [Input] the mapped contraction index
      integer ucenter !< [Output] the symmetry unique center rank
c::local
      integer basis
#include "bas_ibs_sfn.fh"
c
      bas_cn2uce = bas_check_handle(basisin,'bas_cn2uce')
      if(.not.bas_cn2uce) return
c
      basis = basisin + Basis_Handle_Offset
      bas_cn2uce = cont.ge.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2uce) then
        write(LuOut,*)' bas_cn2uce: invalid contraction information '
        write(LuOut,*)' contraction range is 0:',ncont_tot_gb(basis)
        write(LuOut,*)' input contraction was: ',cont
        return
      endif
c
      ucenter = sf_ibs_cn2ce(cont,basis)
      ucenter = sf_ibs_ce2uce(ucenter,basis)
c
      end
*.....................................................................
c
C> \brief Retrieve the range of basis functions corresponding to a shell
c
C> A shell (or equivalently contraction) may contain one or more basis
C> functions. Hence in the total molecular basis a given shell 
C> corresponds to a range of basis functions. This routine retrieve
C> such a range of basis functions in that it
C>
C> - stores the first basis function of a shell in ifirst
C>
C> - stores the last basis function of a shell in ilast
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function bas_cn2bfr(basisin,cont,ifirst,ilast)
c
c returns the first basis function index of a mapped contraction 
c in ifirst and the last basis function index in ilast
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle
      integer cont    !< [Input] the mapped contraction index
      integer ifirst  !< [Output] the first basis function
      integer ilast   !< [Output] the last basis function     
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
      basis = basisin + Basis_Handle_Offset
      bas_cn2bfr = cont.ge.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2bfr) then
        write(LuOut,*)' bas_cn2bfr: invalid contraction information '
        write(LuOut,*)' contraction range is 0:',ncont_tot_gb(basis)
        write(LuOut,*)' input contraction was: ',cont
        return
      endif
c
      ifirst = sf_ibs_cn2bfr(1,cont,basis)
      ilast  = sf_ibs_cn2bfr(2,cont,basis)
c
      return
      end
*.....................................................................
c
C> \brief Retrieve the range of basis functions corresponding to a
C> center
c
C> A center may contain zero or more basis functions. Hence in the total
C> molecular basis a given center corresponds to a range of basis
C> functions. This routine retrieves such a range of basis functions in
C> that it
C>
C> - stores the first basis function of a center in ibflo
C>
C> - stores the last basis function of a center in ibfhi
C>
C> In case that a center does not have any basis functions at all
C>
C> - ibflo is set to 0
C>
C> - ibfhi is set to -1
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function bas_ce2bfr(basis, icent, ibflo, ibfhi)
c
c  returns the basis function range for a given center
c
      implicit none
#include "stdio.fh"
      integer basis !< [Input] the basis set handle
      integer icent !< [Input] the center rank
      integer ibflo !< [Output] the first basis function on the center
      integer ibfhi !< [Output] the last basis function on the center
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
c
C> \brief Retrieves the range of shells for a given center
C>
C> A center may have zero or more shells (or equivalently contractions)
C> associated with it. This routine retrieves the range of shells 
C> associated with a center in that it
C>
C> - stores the rank of the first shell in ifirst
C>
C> - stores the rank of the last shell in ilast
C>
C> In case center does not have any shells (e.g. a point charge) the
C> results are undefined.
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function bas_ce2cnr(basisin,center,ifirst,ilast)
c
c returns the mapped contraction range on a given center
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geom.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle         
      integer center  !< [Input] the center rank 
      integer ifirst  !< [Output] the first mapped contraction
      integer ilast   !< [Output] the last mapped contraction
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
      basis = basisin + Basis_Handle_Offset
      bas_ce2cnr = geom_ncent(ibs_geom(basis),nat)
      if (nat.eq.0.or..not.bas_ce2cnr) then
        bas_ce2cnr = .false.
        write(LuOut,*)' bas_ce2cnr: ERROR '
        write(LuOut,*)' number of centers is zero or weird'
        write(LuOut,*)' nat = ',nat
c..... add diagnostics later
        return
      endif

      bas_ce2cnr = center.gt.0 .and. center.le.nat
      if (.not.bas_ce2cnr) then
        write(LuOut,*)' bas_ce2cnr: invalid center information '
        write(LuOut,*)' center range is 1:',nat
        write(LuOut,*)' input center was : ',center
        return
      endif
c
      ifirst = sf_ibs_ce2cnr(1, center, basis)
      ilast  = sf_ibs_ce2cnr(2, center, basis)
c
      return
      end
*.....................................................................
c
C> \brief Retrieves the center a given basis function is sited on
C>
C> Every Gaussian basis function is sited at one particular expansion
C> center. This routine routine retrieves the particular center on
C> which a given basis function is sited.
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function bas_bf2ce(basisin,testbf,center)
c
c routine to return the center of a given basis function 
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geom.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle         
      integer testbf  !< [Input] the basis function index
      integer center  !< [Output] the center rank
c::local
      integer basis, nat, iat, ibflo, ibfhi, istart, iend
c
#include "bas_ibs_sfn.fh"
c
      bas_bf2ce = bas_check_handle(basisin,'bas_bf2ce')
      if (.not. bas_bf2ce) return
c
      basis = basisin + Basis_Handle_Offset 
      bas_bf2ce = geom_ncent(ibs_geom(basis),nat)
      if (.not.bas_bf2ce .or. nat.le.0) then
        bas_bf2ce = .false.
        write(LuOut,*)' bas_bf2ce: ERROR '
        write(LuOut,*)' number of centers is zero or weird'
        write(LuOut,*)' nat = ',nat
        return
      endif
c
c... linear search through atoms
c
      center = -1
      do 00100 iat = 1,nat
        istart = sf_ibs_ce2cnr(1,iat,basis)
        iend   = sf_ibs_ce2cnr(2,iat,basis)
        if ((iend - istart + 1 ) .le. 0) goto 00100
        ibflo = sf_ibs_cn2bfr(1,istart,basis)
        ibfhi = sf_ibs_cn2bfr(2,iend,basis)
        if (testbf.ge.ibflo.and.testbf.le.ibfhi) then
          center = iat
          return
        endif
00100 continue
c
      end
*.....................................................................
c
C> \brief Retrieves the shell a given basis function is part off
C>
C> Every Gaussian basis function is part of one particular shell
C> (or equivalently contraction). This routine routine retrieves the
C> particular shell a given basis function is part off.
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function bas_bf2cn(basisin,testbf,cont)
c
c returns the mapped contraction index that contains the given 
c basis function index
c 
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_numcont
      logical bas_check_handle
      external bas_numcont
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle         
      integer testbf  !< [Input] the basis function index
      integer cont    !< [Output] the mapped shell index
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
        write(LuOut,*)'bas_bf2cn: could not get number of contractions'
        return
      endif
c
      basis = basisin + Basis_Handle_Offset 
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
c
C> \brief Retrieves the number of basis functions in the molecular
C> basis set
C>
C> The basis set itself contains only the unique definitions of the
C> combinations of primitive basis functions that generate the actual
C> basis functions. I.e. it contains only the contraction coefficients,
C> the exponents and the kind of atom it should be used for. 
C>
C> The molecular basis set is constructed from this information by
C> mapping this specification onto a particular structure. This process
C> Generates the actual specification of basis functions for the atoms.
C> In molecular calculations various properties of the molecular basis
C> set are essential.
C>
C> This routine extracts to total number of basis functions in the
C> molecular basis set.
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function bas_numbf(basisin,nbf)
c
c returns the total number of basis functions of the mapped basis set.
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle         
      integer nbf     !< [Output] the number of basis functions
c::local
      integer basis 
c
      nbf = -6589
      bas_numbf = bas_check_handle(basisin,'bas_numbf')
      if (.not. bas_numbf) return

      basis = basisin + Basis_Handle_Offset 
      nbf = nbf_tot_gb(basis)
      bas_numbf = .true.
      return
      end
*.....................................................................
c
C> \brief Retrieve the total number of primitive basis functions in the
C> molecular basis set
C>
C> A shell of basis functions consists of a number of contracted 
C> Gaussians of a particular angular momentum. The total number of basis
C> functions corresponds to the number of contracted functions. Each
C> individual Gaussian function in a contraction is referred to as
C> a primitive function. Clearly then there is a number of primitive
C> basis functions associated with the molecular basis set as well.
C> This routine retrieves the latter number of primitive basis functions
C> from the basis set.
C>
C> \return Returns .true. if successfull, and .false. otherwise.
c
      logical function bas_num_prim(basisin,nprim)
c
c returns the total number of primitives of the mapped basis set.
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle         
      integer nprim   !< [Output] the number of primitives in mapped basis
c::local
      integer basis 
c
      nprim = -6589
      bas_num_prim = bas_check_handle(basisin,'bas_num_prim')
      if (.not. bas_num_prim) return

      basis = basisin + Basis_Handle_Offset 
      nprim = nprim_tot_gb(basis)
      bas_num_prim = .true.
      return
      end
*.....................................................................
c
C> \brief Retrieve the total number of contraction coefficients of the
C> molecular basis set
C>
C> Shells consist of contracted Gaussian basis functions. With each
C> primitive Gaussian there is associated an exponent and at least
C> one contraction coefficient. 
C>
C> This routine retrieves the total number of contraction coefficients
C> in the molecular basis set.
C>
C> \return Returns .true. if successfull, and .false. otherwise.
c
      logical function bas_num_coef(basisin,ncoef)
c
c returns the total number of coeffsof the mapped basis set.
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle         
      integer ncoef   !< [Output] the number of contraction coefficients
                      !< in the mapped basis
c::local
      integer basis 
c
      ncoef = -6589
      bas_num_coef = bas_check_handle(basisin,'bas_num_coef')
      if (.not. bas_num_coef) return

      basis = basisin + Basis_Handle_Offset 
      ncoef = ncoef_tot_gb(basis)
      bas_num_coef = .true.
      return
      end
*.....................................................................
c
C> \brief Retrieve the names of a basis set instance
C>
C> Each basis set instance has two names. By one name the basis set
C> in memory is known. This name is defined at the time the basis set
C> instance is created. The other name is the one used to store the
C> basis set on the RTDB. 
C>
C> This routine retrieves both basis set names.
C>
C> \return Returns .true. if successfull, and .false. otherwise.
c
      logical function bas_name(basisin,basis_name,trans_name)
c
c returns the name and translated name of the basis set
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"      
#include "inp.fh"
c::functions
c inp_strlen() from inp
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer       basisin    !< [Input] the basis set handle
      character*(*) basis_name !< [Output] the basis set name when loaded
      character*(*) trans_name !< [Output] the basis set name in context 
c::local
      integer basis   ! actual offset into basis arrays
*      integer lenofit ! length of name
c
      bas_name = bas_check_handle(basisin,'bas_name')
      if (.not. bas_name) return
c
      basis = basisin + Basis_Handle_Offset 
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
c
C> \brief Retrieves the tag of a given shell
C>
C> \return Returns .true. if successfull, and .false. otherwise.
c
      logical function bas_cont_tag(basisin,icont,tagout)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
c::-function
      logical bas_check_handle
      external bas_check_handle
c::-passed
      integer basisin      !< [Input] the basis set handle
      integer icont        !< [Input] the shell index
      character*(*) tagout !< [Output] the shell tag
c::-local
      integer center, ucenter, basis
      integer len_tagout
#include "bas_ibs_sfn.fh"
c
      bas_cont_tag = bas_check_handle(basisin,'bas_cont_tag')
      if (.not.bas_cont_tag) return

      basis = basisin + Basis_Handle_Offset 
      
      center = sf_ibs_cn2ce(icont,basis)
      ucenter = sf_ibs_ce2uce(center,basis)
      len_tagout = len(tagout)
      tagout = bs_tags(ucenter,basis)(1:len_tagout)
      end
*.....................................................................
c
C> \brief Retrieves generic information of a given shell
C>
C> \return Returns .true. if successfull, and .false. otherwise.
c
      logical function bas_continfo(basisin,icont,
     &       type,nprimo,ngeno,sphcart)
c
c  returns the generic information about the given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis handle
      integer icont   !< [Input] the shell index
      integer type    !< [Output] the shell type (sp/s/p/d/..)
      integer nprimo  !< [Output] the number of primitives
      integer ngeno   !< [Output] the number of contractions
      integer sphcart !< [Output] 0/1 for cartesian/spherical harmonic
                      !< basis functions
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

      basis = basisin + Basis_Handle_Offset 
c
      icontmax = ncont_tot_gb(basis)
c
      if (.not.(icont.ge.0.and.icont.le.icontmax)) then
        write(LuOut,*)' bas_continfo: ERROR '
        write(LuOut,*)' contraction range for basis is 0:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
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
c
C> \brief Retrieves the total number of shells in the molecular
C> basis set
C>
C> \return Returns .true. if successfull, and .false. otherwise.
c
      logical function bas_numcont(basisin,numcont)
c
c returns the total number of mapped contractions/shells for the 
c given basis set
c
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle
      integer numcont !< [Output] the number of mapped shells
c::local
      integer basis
c
      numcont = -6589
      bas_numcont = bas_check_handle(basisin,'bas_numcont')
      if (.not.bas_numcont) return

      basis = basisin + Basis_Handle_Offset 

      numcont = ncont_tot_gb(basis)

      bas_numcont = .true.
      return
      end
*.....................................................................
c
C> \brief Retrieves the maximum number of basis functions in any shell
C>
C> This function scans through all the shells in a basis set and finds
C> the maximum number of basis functions. It deals with Cartesian and
C> spherical harmonic basis functions as well as segmented and 
C> generally contracted basis sets.
C>
C> \return Returns .true. if successful, and .false. otherwise
c
      logical function bas_nbf_cn_max(basisin,nbf_max)
      implicit none
c
c  calculate, return and store maximum basis function block size
c   for all contractions in a given basis. 
c
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "stdio.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
      integer nbf_from_ucont
      external nbf_from_ucont
c:: passed
      integer basisin !< [Input]  the basis set handle
      integer nbf_max !< [Output] the maximum number of basis functions
                      !< in any shell
c:local
      integer basis, myucont, i, mynbf
c
      bas_nbf_cn_max = bas_check_handle(basisin,'bas_nbf_cn_max')
      if (.not. bas_nbf_cn_max ) then
        write(LuOut,*) 'bas_nbf_cn_max: basis handle not valid '
        return
      endif
c
      basis = basisin + Basis_Handle_Offset
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
C>
C> \brief Retrieves the maximum number of basis functions on any center
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function bas_nbf_ce_max(basisin,nbf_max)
      implicit none
c
c  calculate, return and store maximum basis function block size
c   for all contractions in a given basis. 
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geom.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin  !< [Input] the basis set handle
      integer nbf_max  !< [Output] the maximum number of basis functions
                       !< on any center
c:local
      integer basis, mynat, iat, mylo, myhi, mynbf
      integer myclo, mychi, mynumcont
c
#include "bas_ibs_sfn.fh"
c
      bas_nbf_ce_max = bas_check_handle(basisin,'bas_nbf_ce_max')
      if (.not. bas_nbf_ce_max ) then
        write(LuOut,*) 'bas_nbf_ce_max: basis handle not valid '
        return
      endif
c
      basis = basisin + Basis_Handle_Offset
c
      bas_nbf_ce_max = geom_ncent(ibs_geom(basis),mynat)
      if (mynat.le.0.or. .not.bas_nbf_ce_max) then
        write(LuOut,*)' bas_nbf_ce_max: ERROR '
        write(LuOut,*)' number of centers is zero or weird'
        write(LuOut,*)' nat = ',mynat
        return
      endif
c
      nbf_max = -1
      do 00100 iat = 1,mynat
        myclo = (sf_ibs_ce2cnr(1,iat,basis))
        mychi = (sf_ibs_ce2cnr(2,iat,basis))
        mynumcont = mychi - myclo + 1
* Bqs with no basis functions have mynumcont < 1
        if (mynumcont.gt.0) then
          mylo  = sf_ibs_cn2bfr(1,myclo,basis)
          myhi  = sf_ibs_cn2bfr(2,mychi,basis)
          mynbf = myhi - mylo + 1
          nbf_max = max(nbf_max, mynbf)
        endif
00100 continue
c
      return
      end
*.....................................................................
c
C> \brief Retrieve the geometry that was specified when loading the
C> basis
C>
C> Simply look up the geometry that was stored with the basis set when
C> the basis set was loaded from the RTDB.
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function bas_geom(basisin,geom)	
      implicit none 
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin !< [Input]  the basis set handle
      integer basis
      integer geom    !< [Output] the geometry used to load basis set
c
      if (.not.bas_check_handle(basisin,'bas_geom'))
     & call errquit('bas_geom: handle invalid',911, BASIS_ERR)
c
      basis = basisin + Basis_Handle_Offset
c
      geom = ibs_geom(basis)
      bas_geom = .true.
      end
*...............................................
c
C> \brief Retrieve the maximum number of contraction coefficients
C> associated with any shell
C>
C> This function scans through all shells in a basis set and checks
C> the number of contraction coefficients for each. It takes into
C> account segmented and generally contracted shells. The maximum
C> number of contraction coefficients for any shell is passed back to
C> the calling routine.
C>
C> \return Return .true. if successfull, .false. otherwise.
c
      logical function bas_ncoef_cn_max(basisin, ncoef_max)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
      integer basisin   !< [Input]  the basis set handle
      integer ncoef_max !< [Output] the maximum number of contraction 
                        !< coefficients
c
      integer basis
      integer nu_cont, iu_cont, my_prim, my_gen
c
      basis = basisin + Basis_Handle_Offset
      if (.not.bsactive(basis)) call errquit
     &    ('bas_ncoef_cn_max: basis handle invalid',911, BASIS_ERR)
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
c
C> \brief Retrieve the maximum number of primitive basis functions
C> associated with any shell
C>
C> This function scans through all shells in a basis set and checks
C> the number of primitive basis functions (or equivalently the
C> number of exponents) for each. The maximum number of primitive basis
C> basis functions for any shell is passed back to the calling routine.
C>
C> \return Return .true. if successfull, .false. otherwise.
c
      logical function bas_nprim_cn_max(basisin, nprim_max)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
      integer basisin   !< [Input]  the basis set handle
      integer nprim_max !< [Output] the maximum number of primitives in
                        !< any shell
c
      integer basis
      integer nu_cont, iu_cont, my_prim
c
      basis = basisin + Basis_Handle_Offset
      if (.not.bsactive(basis)) call errquit
     &    ('bas_nprim_cn_max: basis handle invalid',911, BASIS_ERR)
      nu_cont = infbs_head(HEAD_NCONT,basis)
      nprim_max = 0
      do 00100 iu_cont = 1,nu_cont
        my_prim = infbs_cont(CONT_NPRIM,iu_cont,basis)
        nprim_max = max(nprim_max,my_prim)
00100 continue
      bas_nprim_cn_max = .true.
      end
C> @}
*.....................................................................
      logical function bas_norm_get(basisin,norm_id)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c
      integer basisin ! [input] basis set handle
      integer norm_id ! [output] Normalization id type
c
      integer basis
c
      basis = basisin + Basis_Handle_Offset
      norm_id = bas_norm_id(basis)
      bas_norm_get = norm_id.ge.BasNorm_lo.and.norm_id.le.BasNorm_hi
      end
*.....................................................................
      logical function bas_norm_set(basisin,norm_id)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c
      integer basisin ! [input] basis set handle
      integer norm_id ! [input] Normalization id type
c
      integer basis
c
      basis = basisin + Basis_Handle_Offset
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
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "stdio.fh"
c
      integer basisin ! [input] basis set handle
c
      integer basis, norm_id
c
      basis = basisin + Basis_Handle_Offset
      norm_id = bas_norm_id(basis)
      bas_norm_print = .true.
      if      (norm_id.eq.BasNorm_UN)  then
        write(luout,*)' basis is unnormalized'
      else if (norm_id.eq.BasNorm_STD) then
        write(luout,*)' basis has standard normalization via int_norm'
      else if (norm_id.eq.BasNorm_2c)  then
        write(luout,*)
     &        ' basis has dft/fitting normalization via int_norm_2c'
      else if (norm_id.eq.BasNorm_rel)  then
        write(luout,*)
     &        ' basis has relativistic normalization via int_norm'
      else
        write(luout,*)' basis handle: ',basisin
        write(luout,*)' norm_id = ',norm_id
        write(luout,*)' does not match a known normalization mode'
        bas_norm_print = .false.
      endif
      end
*.....................................................................
C> \ingroup bas
C> @{
C>
C> \brief Print unique in core ECP information
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function ecp_print(ecpidin)
c
c routine to print unique ecpid information that is in core
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
#include "geom.fh"
#include "stdio.fh"
#include "bas_starP.fh"
#include "errquit.fh"
c
c function declarations
c      
      logical ecp_check_handle
      external ecp_check_handle
c
cc AJL/Begin/SPIN-POLARISED ECPs
      logical ecp_get_high_chan
      external ecp_get_high_chan
cc AJL/End
c
c:: passed
      integer ecpidin !< [Input] the basis set handle
c:: local
      integer mytags, myucont, myprim, mycoef, ecpid
      integer i,j,k,l, ifcont, mygen, mytype, irexptr, iexptr, icfptr
      integer atn, len_tag, len_ele
c
cc AJL/Begin/SPIN-POLARISED ECPs
      integer mychannel
      integer channels
      character*5 channel_type(0:2)
cc AJL/End
c
      character*2 symbol
      character*16 element
      character*9 cartorsph
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
c
cc AJL/Begin/SPIN-POLARISED ECPs
cc Aesthetic clean up for print. Is this necessary?
C Can I just set the array without the if statement? To test.
      channels = 1
      if (.not.ecp_get_high_chan(ecpidin,channels))
     &      call errquit('ecp_print error',911, BASIS_ERR)

      channel_type(0) = 'Both'
      channel_type(1) = 'Alpha'
      channel_type(2) = 'Beta'
cc AJL/End
c
      ecp_print = .true.
      ecpid = ecpidin + Basis_Handle_Offset
c
      ecp_print = ecp_check_handle(ecpidin,'ecp_print')
      if (.not. ecp_print) return
c
c print basis set information
c      
      if (bas_spherical(ecpid)) then
         cartorsph = 'spherical'
      else
         cartorsph = 'cartesian'
      endif
      write(LuOut,1)bs_name(ecpid)(1:inp_strlen(bs_name(ecpid))), 
     $    bs_trans(ecpid)(1:inp_strlen(bs_trans(ecpid))),cartorsph
 1    format('                 ECP       "',a,'" -> "',a,'"',' (',a,')'/
     $       '                -----')
      mytags  = infbs_head(HEAD_NTAGS,ecpid)
      if (mytags.le.0) then
        write(LuOut,*)'No explicit ECP functions are defined !'
        write(LuOut,*)
c
c there could be star tags defined, so check that before returning
c
        goto 00010
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
        call util_print_centered(LuOut, buffer, len_tag/2 + 1, .true.)
        
        myucont = infbs_tags(TAG_NCONT,i,ecpid)
c        
        ifcont = infbs_tags(TAG_FCONT,i,ecpid)
c     
        write(LuOut,6)
cc AJL/Begin/SPIN-POLARISED ECPs
c 6      format(
c     $      '          R-exponent    Exponent     Coefficients '/
 6      format(13(' '),
     $      'Channel    R-exponent     Exponent     Coefficients'/
cc AJL/End
     $      '         ------------ ',57('-'))
        do 00200 j=1,myucont
          myprim = infbs_cont(CONT_NPRIM,ifcont,ecpid)
          mygen  = infbs_cont(CONT_NGEN,ifcont,ecpid)
          
          mytype = infbs_cont(CONT_TYPE, ifcont, ecpid)
          shell_type = ctype(mytype)
          iexptr  = infbs_cont(CONT_IEXP, ifcont,ecpid) - 1
          icfptr  = infbs_cont(CONT_ICFP, ifcont,ecpid) - 1
          irexptr = infbs_cont(CONT_IREXP,ifcont,ecpid) - 1
          mychannel = infbs_cont(CONT_CHANNEL,ifcont,ecpid)
c
          do 00300 k=1,myprim
            write(LuOut,7) j, shell_type,
     &          channel_type(mychannel),
     &          sf_exndcf((irexptr+k),ecpid),
     &          sf_exndcf((iexptr+k),ecpid),
     &          (sf_exndcf((icfptr+k+(l-1)*myprim),ecpid),l=1,mygen)
00300     continue
          write(LuOut,*)
          ifcont = ifcont + 1
00200   continue
00100 continue
c
c  Check if we have star tag definitions in the basis set
c
00010 if (star_nr_tags .gt. 0) then
         write(LuOut,*)
         write(LuOut,8)
  8      format(' In addition, one or more string tags have been',
     &         ' defined containing a * .'/' These tags, and ',
     &         'their exceptions list are printed below.'//
     &         ' Tag ',12(' '),' Description ',18(' '),
     &         ' Excluding '/' ',16('-'),' ',30('-'),' ',30('-'))
         do i=1,star_nr_tags
           k = 1
           if (i .gt. 1) k = star_nr_excpt(i-1) + 1
           write(LuOut,9) star_tag(i), star_bas_typ(i)(1:30),
     &          (star_excpt(j)(1:inp_strlen(star_excpt(j))),
     &           j=k,star_nr_excpt(i))
  9        format(1x,a16,1x,a30,1x,10(a))
         enddo
         write(LuOut,*)
         write(LuOut,*)
      endif
c
c  If geom is set print out the info about total basis info 
c  associated with the geometry also
c
c  ... not done yet
c
      return
cc AJL/Begin/SPIN-POLARISED ECPs
c 7          format(1x,i2,1x,a3,1x,f9.2,2x,f14.6,20f15.6)
 7          format(1x,i2,1x,a3,7x,a5,3x,f9.2,2x,f14.6,20f15.6)
cc AJL/End
      end
*.....................................................................
C>
C> \brief Print unique in core spin-orbit potential information
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function so_print(soidin)
c
c routine to print unique soid information that is in core
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
#include "geom.fh"
#include "stdio.fh"
#include "bas_starP.fh"
c
c function declarations
c      
      logical so_check_handle
      external so_check_handle
c:: passed
      integer soidin !< [Input] the basis set handle
c:: local
      integer mytags, myucont, myprim, mycoef, soid
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
      so_print = .true.
      soid = soidin + Basis_Handle_Offset
c
      so_print = so_check_handle(soidin,'so_print')
      if (.not. so_print) return
c
c print basis set information
c      
      write(LuOut,1)bs_name(soid)(1:inp_strlen(bs_name(soid))), 
     $    bs_trans(soid)(1:inp_strlen(bs_trans(soid)))
 1    format('             SO Potential "',a,'" -> "',a,'"'/
     $    '                      -----')
      if (bas_spherical(soid)) write(LuOut,2)
 2    format('    SO Basis is spherical, 5d, 7f, 9g ... ')
      mytags  = infbs_head(HEAD_NTAGS,soid)
      if (mytags.le.0) then
        write(LuOut,*)'No explicit SO ECP functions are defined !'
        write(LuOut,*)
c
c there could be star tags defined, so check that before returning
c
        goto 00010
      endif
c
      myucont = infbs_head(HEAD_NCONT,soid)
      myprim  = infbs_head(HEAD_NPRIM,soid)
      mycoef  = infbs_head(HEAD_NCOEF,soid)
c
      do 00100 i=1,mytags
        
        if (geom_tag_to_element(bs_tags(i,soid), symbol, element,
     $      atn)) then
          len_tag = inp_strlen(bs_tags(i,soid))
          len_ele = inp_strlen(element)
          write (buffer,
     &        '(a,'' ('',a,'')'' )')
     &        bs_tags(i,soid)(1:len_tag), element(1:len_ele)
        else
          buffer = bs_tags(i,soid)
        endif
        len_tag = inp_strlen(buffer)
        call util_print_centered(LuOut, buffer, len_tag/2 + 1, .true.)
        
        myucont = infbs_tags(TAG_NCONT,i,soid)
c        
        ifcont = infbs_tags(TAG_FCONT,i,soid)
c     
        write(LuOut,6)
 6      format(
     $      '          R-exponent    Exponent     Coefficients '/
     $      '         ------------ ',60('-'))
        do 00200 j=1,myucont
          myprim = infbs_cont(CONT_NPRIM,ifcont,soid)
          mygen  = infbs_cont(CONT_NGEN,ifcont,soid)
          
          mytype = infbs_cont(CONT_TYPE, ifcont, soid)
          shell_type = ctype(mytype)
          iexptr  = infbs_cont(CONT_IEXP, ifcont,soid) - 1
          icfptr  = infbs_cont(CONT_ICFP, ifcont,soid) - 1
          irexptr = infbs_cont(CONT_IREXP,ifcont,soid) - 1
          do 00300 k=1,myprim
            write(LuOut,7) j, shell_type,
     &          sf_exndcf((irexptr+k),soid),
     &          sf_exndcf((iexptr+k),soid),
     &          (sf_exndcf((icfptr+k+(l-1)*myprim),soid),l=1,mygen)
00300     continue
          write(LuOut,*)
          ifcont = ifcont + 1
00200   continue
00100 continue
c
c  Check if we have star tag definitions in the basis set
c
00010 if (star_nr_tags .gt. 0) then
         write(LuOut,*)
         write(LuOut,8)
  8      format(' In addition, one or more string tags have been',
     &         ' defined containing a * .'/' These tags, and ',
     &         'their exceptions list are printed below.'//
     &         ' Tag ',12(' '),' Description ',18(' '),
     &         ' Excluding '/' ',16('-'),' ',30('-'),' ',30('-'))
         do i=1,star_nr_tags
           k = 1
           if (i .gt. 1) k = star_nr_excpt(i-1) + 1
           write(LuOut,9) star_tag(i), star_bas_typ(i)(1:30),
     &          (star_excpt(j)(1:inp_strlen(star_excpt(j))),
     &           j=k,star_nr_excpt(i))
  9        format(1x,a16,1x,a30,1x,10(a))
         enddo
         write(LuOut,*)
         write(LuOut,*)
      endif
c
c  If geom is set print out the info about total basis info 
c  associated with the geometry also
c
c  ... not done yet
c
      return
 7    format(1x,i2,1x,a3,1x,f9.2,2x,f14.6,20f15.6)
      end
*.....................................................................
C>
C> \brief Set the name of the ECP for a basis set
C>
C> \return Return .true. when successfull, and .false. otherwise
c
      logical function bas_set_ecp_name(basisin,ecp_name)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin        !< [Input] the basis set handle
      character*(*) ecp_name !< [Input] the ECP name
c
      name_assoc(1,(basisin+Basis_Handle_Offset)) =
     &    ecp_name
      bas_set_ecp_name = .true.
      end
*.....................................................................
C>
C> \brief Set the name of the spin-orbit potential for a basis set
C>
C> \return Return .true. when successfull, and .false. otherwise
c
      logical function bas_set_so_name(basisin,so_name)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin        !< [Input] the basis set handle
      character*(*) so_name  !< [Input] the spin-orbit potential name
c
      name_assoc(2,(basisin+Basis_Handle_Offset)) =
     &    so_name
      bas_set_so_name = .true.
      end
*.....................................................................
C>
C> \brief Retrieve the name of the ECP of a basis set
C>
C> \return Return .true. when successfull, and .false. otherwise
c
      logical function bas_get_ecp_name(basisin,ecp_name)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin        !< [Input] the basis set handle
      character*(*) ecp_name !< [Output] the ECP name
c
      ecp_name =
     &    name_assoc(1,(basisin+Basis_Handle_Offset)) 
      bas_get_ecp_name = .true.
      end
*.....................................................................
C>
C> \brief Retrieve the name of the spin-orbit potential of a basis set
C>
C> \return Return .true. when successfull, and .false. otherwise
c
      logical function bas_get_so_name(basisin,so_name)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin       !< [Input] the basis set handle
      character*(*) so_name !< [Output] the spin-orbit potential name
c
      so_name =
     &    name_assoc(2,(basisin+Basis_Handle_Offset)) 
      bas_get_so_name = .true.
      end
*.....................................................................
C>
C> \brief Associate an ECP with a given basis set
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function bas_set_ecp_handle(basisin,ecp_handle)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin    !< [Input] the molecular basis set handle
      integer ecp_handle !< [Input] the EPC handle
c
      handle_assoc(1,(basisin+Basis_Handle_Offset)) =
     &    ecp_handle
      bas_set_ecp_handle = .true.
      end
*.....................................................................
C>
C> \brief Retrieve an ECP of a given basis set
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function bas_get_ecp_handle(basisin,ecp_handle)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin    !< [Input] the molecular basis set handle
      integer ecp_handle !< [Output] the ECP handle
c
      ecp_handle =
     &    handle_assoc(1,(basisin+Basis_Handle_Offset)) 
      if (ecp_handle.ne.0) then
        bas_get_ecp_handle = .true.
      else
        bas_get_ecp_handle = .false.
      endif
      end
*.....................................................................
C>
C> \brief Set the molecular basis set for an ECP
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function ecp_set_parent_handle(ecp_handle,basis_handle)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
      integer ecp_handle   !< [Input] the ECP handle
      integer basis_handle !< [Input] the molecular basis set handle
c
      logical bas_check_handle, ecp_check_handle
      external bas_check_handle, ecp_check_handle
c      
      if (.not.bas_check_handle(basis_handle,'ecp_set_parent_handle'))
     &      call errquit
     &      ('ecp_set_parent_handle: basis_handle invalid',911,
     &       BASIS_ERR)
      if (.not.ecp_check_handle(ecp_handle,'ecp_set_parent_handle'))
     &      call errquit
     &      ('ecp_set_parent_handle: not an ECP ...'//
     $     ' check for name conflict between basis and ECP',911,
     &       BASIS_ERR)
      parent_assoc(1,(ecp_handle+Basis_Handle_Offset)) =
     &      basis_handle
      bas_nassoc(basis_handle+Basis_Handle_Offset) =
     &    bas_nassoc(basis_handle+Basis_Handle_Offset) + 1 
      ecp_set_parent_handle = .true.
      end
*.....................................................................
C>
C> \brief Retrieve the molecular basis set of an ECP
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function ecp_get_parent_handle(ecp_handle,basis_handle)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
      integer ecp_handle   !< [Input] the ECP handle
      integer basis_handle !< [Output] the molecular basis set handle
c
      logical bas_check_handle, ecp_check_handle
      external bas_check_handle, ecp_check_handle
c      
      if (.not.ecp_check_handle(ecp_handle,'ecp_get_parent_handle'))
     &      call errquit
     &      ('ecp_get_parent_handle: ecp_handle invalid',911,
     &       BASIS_ERR)
      basis_handle = parent_assoc(1,ecp_handle+Basis_Handle_Offset)
      if (basis_handle.ne.0) then
        ecp_get_parent_handle = .true.
      else
        ecp_get_parent_handle = .false.
      endif
      if (.not.bas_check_handle(basis_handle,'ecp_get_parent_handle'))
     &      call errquit
     &      ('ecp_get_parent_handle: stored basis_handle invalid',911,
     &       BASIS_ERR)
      end
*.....................................................................
C>
C> \brief Set the spin-orbit potential for a molecular basis set
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function bas_set_so_handle(basisin,so_handle)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin   !< [Input] the molecular basis set handle
      integer so_handle !< [Input] the spin-orbit potential handle
c
      handle_assoc(2,(basisin+Basis_Handle_Offset)) =
     &    so_handle
      bas_set_so_handle = .true.
      end
*.....................................................................
C>
C> \brief Retrieve the spin-orbit potential of a molecular basis set
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function bas_get_so_handle(basisin,so_handle)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin   !< [Input] the molecular basis set handle
      integer so_handle !< [Output] the spin-orbit potential handle
c
      so_handle =
     &    handle_assoc(2,(basisin+Basis_Handle_Offset)) 
      if (so_handle.ne.0) then
        bas_get_so_handle = .true.
      else
        bas_get_so_handle = .false.
      endif
      end
*.....................................................................
C>
C> \brief Set the molecular basis set for a spin-orbit potential
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function so_set_parent_handle(so_handle,basis_handle)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
      integer so_handle    !< [Input] the spin-orbit potential handle
      integer basis_handle !< [Input] the molecular basis set handle
c
      logical bas_check_handle, so_check_handle
      external bas_check_handle, so_check_handle
c      
      if (.not.bas_check_handle(basis_handle,'so_set_parent_handle'))
     &      call errquit
     &      ('so_set_parent_handle: basis_handle invalid',911,
     &       BASIS_ERR)
      if (.not.so_check_handle(so_handle,'so_set_parent_handle'))
     &      call errquit
     &      ('so_set_parent_handle: so_handle invalid',911,
     &       BASIS_ERR)
      parent_assoc(2,(so_handle+Basis_Handle_Offset)) =
     &      basis_handle
      bas_nassoc(basis_handle+Basis_Handle_Offset) =
     &    bas_nassoc(basis_handle+Basis_Handle_Offset) + 1 
      so_set_parent_handle = .true.
      end
*.....................................................................
C>
C> \brief Retrieve the molecular basis set of a spin-orbit potential
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function so_get_parent_handle(so_handle,basis_handle)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
      integer so_handle    !< [Input] the spin-orbit potential handle
      integer basis_handle !< [Output] the molecular basis set handle
c
      logical bas_check_handle, so_check_handle
      external bas_check_handle, so_check_handle
c      
      if (.not.so_check_handle(so_handle,'so_get_parent_handle'))
     &      call errquit
     &      ('so_get_parent_handle: so_handle invalid',911, BASIS_ERR)
      basis_handle =
     &    parent_assoc(2,(so_handle+Basis_Handle_Offset))
      if (basis_handle.ne.0) then
        so_get_parent_handle = .true.
      else
        so_get_parent_handle = .false.
      endif
      if (.not.bas_check_handle(basis_handle,'so_get_parent_handle'))
     &      call errquit
     &      ('so_get_parent_handle: stored basis_handle invalid',911,
     &       BASIS_ERR)
      end
C
C> \brief Summarize the internal data structures
C
      subroutine bas_print_allocated_info(msg)
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "global.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
      character*(*) msg !< [Input] A message to identify the print
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
        call ga_brdcst(1234,inode_use,
     &    ma_sizeof(mt_int, 1, mt_byte),0)
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
          call util_flush(luout)
        endif
      enddo
      end
c
C> \brief Matches a tag from the geometry to a tag in the basis set
C>
C> Associated with a shell of basis functions is a tag identifying the
C> kind of center the shell should reside on. These centers should of
C> course appear in the geometry. This routine takes a tag from the
C> geometry and matches it to a tag from the basis set. The index of
C> the first tag in the basis set that matches the geometry tag is
C> retrieved.
C>
C> \return Return .true. when a match was found, and .false. otherwise
c
      logical function bas_match_tags(tag_from_geom,basisin,btag)
      implicit none
#include "errquit.fh"
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
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      logical geom_tag_to_element
      external geom_tag_to_element
c::passed
      character*(*) tag_from_geom !< [Input] the geometry tag
      integer basisin             !< [Input] the basis set handle
      integer btag                !< [Output] lexical basis set tag index
*     bas_match_tags              ! [output] true if they match
c::local
      integer basis           ! lexical basis index
      integer nbtgs           ! number of basis set tags
      character*16 gstring    ! geometry or basis tags can only be *16
      character*16 bsmatch    ! the matched basis set tag
      integer lgstring        ! length of gstring
      integer lgstring_old    ! length of gstring at first test
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
      debug = .false.
c
      bas_match_tags = .false.
      bsmatch = ' '
c
      basis = basisin+Basis_Handle_Offset
c
      gstring = tag_from_geom               !  copy tag to work on it
c
      if (gstring(1:4).eq.'bqgh') then
        gstring = gstring(5:)
      endif
c
c... first match full geometry tag to full basis tag list 
c          did the user specifically assign a tag ?

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
      lgstring = inp_strlen(gstring)
      lnotalpha = lgstring+1
      do i = 1,lgstring
        if (.not.(inp_match(na2z,.false.,
     &      gstring(i:i),a2z,ind))) then ! compare character to alpha (case-less test)
          lnotalpha = i
          goto 00001
        endif
      enddo
00001 continue
      do i = lnotalpha,lgstring
        gstring(i:i) = ' '
      enddo
      if (debug) write(luout,*)' the string:', tag_from_geom,
     &    'has the user substring of ',gstring
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
c... match user substring to basis set tags (only if it gstring has a different length)
      lgstring_old = lgstring
      lgstring = inp_strlen(gstring)
      if (lgstring_old.gt.lgstring) then
        do i = 1, nbtgs
          if (gstring.eq.bs_tags(i,basis)) then
            bas_match_tags = .true.
            btag = i
            bsmatch = bs_tags(i,basis)
            goto 00009
          endif
        enddo
      endif
c
c... now get symbol and element names for each tag and match those
c    geometry tag is based on users substring
c    basis set tag is based on user input to the basis object
c... bq's with basis functions must match from the above rules!
c
      status = geom_tag_to_element(gstring,g_sym,g_elem,g_atn)
      if (.not.status)then
        if (.not.(g_sym.eq.'bq')) then
          write(luout,*)'geometry tag <',gstring,
     &        '> could not be matched to an element symbol'
          call errquit('bas_match_tags: fatal error ',911, BASIS_ERR)
        endif
      endif
      if (g_sym.eq.'bq') goto 00009   ! bq labes with basis functions must match from above
      do i = 1, nbtgs
        status =
     &      geom_tag_to_element(bs_tags(i,basis),b_sym,b_elem,b_atn)
        if (.not.status) then
          if (.not.(b_sym.eq.'bq')) then
            write(luout,*)'basis tag',bs_tags(i,basis),
     &          ' could not be matched to an element symbol'
            call errquit('bas_match_tags: fatal error ',911, BASIS_ERR)
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
      if (debug)
     &    write(luout,*)'bas_match_tags:debug: no match for tag <',
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
c
C> \brief Work out whether multipoles can be calculated for a basis set
C>
C> The multipole code can expand the molecular potential in multipoles.
C> However, it cannot handle all basis sets. In particular it cannot
C> handle generally contracted or SP shells.
C>
C> \return Return .true. if this basis set is compatible with the 
C> multipole code, and .false. otherwise
c
      logical function bas_cando_mpoles(basis)
      implicit none
#include "errquit.fh"
      integer basis !< [Input] the basis set handle
c
c     Return true if it's possible to compute multipoles for this basis
c
c     The multipole code cannot handle general contractions or SP shells.
c
      integer nshell, ishell, type, nprim, ngen, sphcart
      logical bas_numcont, bas_continfo
      external bas_numcont, bas_continfo
c
      bas_cando_mpoles = .false.
c
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('multipole: basis bad?',0, BASIS_ERR)
      do ishell = 1, nshell
         if (.not. bas_continfo(basis, ishell, type, 
     $        nprim, ngen, sphcart)) call errquit
     $        ('multipole: basis bad?',0, BASIS_ERR)
         if (type.lt.0 .or. ngen.gt.1) return
      enddo
c
      bas_cando_mpoles = .true.
c
      end
c
C> \brief Lookup whether a basis set consists of spherical harmonic
C> functions
C>
C> At present cartesian and spherical harmonic basis sets are supported.
C> This routine looks up whether the specified basis set is a spherical
C> harmonic basis set.
C>
C> \return Return .true. if the basis set consists of spherical 
C> harmonic functions, and .false. otherwise
C
      logical function bas_is_spherical(basisin)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"      
#include "stdio.fh"
*::passed
      integer basisin !< [Input] the basis set handle
*::local
      integer basis   ! lexical index of basis set handle
*::functions
      logical bas_check_handle
      external bas_check_handle
*::code
      bas_is_spherical = bas_check_handle(basisin,'bas_is_spherical')
      if (.not.bas_is_spherical) then
        write(luout,*)'invalid/inactive basis handle:',basisin
        call errquit('bas_is_spherical: fatal error',911, BASIS_ERR)
      endif
      basis = basisin + Basis_Handle_Offset
      bas_is_spherical = bas_spherical(basis)
      return 
*      
      end
C> @}
      subroutine bas_dump_info(msg)
      implicit none
#include "stdio.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
      character*(*) msg
c
      integer llen
      integer basis
      integer ntags
      llen = inp_strlen(msg)
      write(luout,00001)
      write(luout,*)'bas_dump_info: <',msg(1:llen),'>'
      write(luout,*)'bas_dump_info: bsversion:', bsversion
      do basis = 1, nbasis_bsmx
        if (bsactive(basis)) then
          write(luout,*)'---------------------------------------------'
          write(luout,*)'---------------------------------------------'
          write(luout,*)'handle(computed) :',
     &        (basis-basis_handle_offset)
          write(luout,*)'geometry loader  :',
     &        ibs_geom(basis)
          write(luout,*)'bs_name          :',
     &        bs_name(basis)(1:len_bs_name(basis))
          write(luout,*)'bs_trans         :',
     &        bs_trans(basis)(1:len_bs_trans(basis))
          llen = inp_strlen(name_assoc(1,basis))
          write(luout,*)'ecp_name_assoc   :',
     &        name_assoc(1,basis)(1:llen)
          llen = inp_strlen(name_assoc(2,basis))
          write(luout,*)'so_name_assoc   :',
     &        name_assoc(2,basis)(1:llen)
          ntags = infbs_head(Head_Ntags,basis)
          write(luout,*)'head ntags       :',ntags
          write(luout,*)'head ncont       :',
     &        infbs_head(Head_Ncont,basis)
          write(luout,*)'head nprim       :',
     &        infbs_head(Head_Nprim,basis)
          write(luout,*)'head ncoef       :',
     &        infbs_head(Head_Ncoef,basis)
          write(luout,*)'head excfptr     :',
     &        infbs_head(Head_Excfptr,basis)
          write(luout,*)'head sph         :',
     &        infbs_head(Head_Sph,basis)
          write(luout,*)'head ecp         :',
     &        infbs_head(Head_Ecp,basis)
          write(luout,*)'bas_spherical    :',
     &        bas_spherical(basis)
          write(luout,*)'bas_any_gc       :',
     &        bas_any_gc(basis)
          write(luout,*)'bas_any_sp_shell :',
     &        bas_any_sp_shell(basis)
          write(luout,*)'bas_norm_id      :',
     &        bas_norm_id(basis)
          write(luout,*)'angular_bs       :',
     &        angular_bs(basis)
          write(luout,*)'nbfmax_bs        :',
     &        nbfmax_bs(basis)
          write(luout,*)'ecp_handle_assoc :',
     &        handle_assoc(1,basis)
          write(luout,*)'ecp_parent_assoc :',
     &        parent_assoc(1,basis)
          write(luout,*)'so_handle_assoc :',
     &        handle_assoc(2,basis)
          write(luout,*)'so_parent_assoc :',
     &        parent_assoc(2,basis)
          write(luout,*)'num of assoc    :',
     &        bas_nassoc(basis)
          write(luout,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
          write(luout,*)'---------------------------------------------'
        endif
      enddo
      write(luout,00001)
00001 format(1x,80('='))
      end
*.....................................................................
C> \ingroup bas
C> @{
C>
C> \brief Find the length of the longest contraction in a basis set
C>
C> Retrieves the length of the longest contraction (i.e. the maximum
C> number of Gaussian primitive terms) in the specified basis set,
C> for any contraction/shell, atom, etc.
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      Logical function bas_ncontr_cn_max( basis_hand, ncontr_max )
c    
c     returns the largest number of contractions, i.e., columns of 
c     contraction coefficients for a given contraction set, 
c     in the basis input, for any contraction set in the basis
c    
      implicit none
      logical bas_numcont, bas_continfo
      external bas_numcont, bas_continfo
c::args
      integer basis_hand !< [Input] the basis set handle
      integer ncontr_max !< [Output] the maximum contraction length
c::locals
      logical LResult
      integer itype, nprimo, isphcart
      integer ncontr, icontset, numcontset
c::exec
      LResult = .true.
      LResult = LResult .and. bas_numcont(basis_hand, numcontset)
c    
c     loop over contraction sets; find largest number of contractions
c     (i.e., columns of contraction coefficients)
c    
      ncontr_max = 0
      do icontset = 1,numcontset
        LResult = LResult .and.
     &       bas_continfo( basis_hand, icontset, itype, nprimo,
     &                     ncontr, isphcart)

      ncontr_max = max( ncontr_max, ncontr )
      enddo 

      bas_ncontr_cn_max = LResult

      return 
      end 
*----------------------------------------------------------------------
C> \brief Retrieve the number of centers with either ECPs or spin-orbit
C> potentials
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function ecpso_ncent(geom,soidin,ecpidin,num_cent)
      implicit none
#include "errquit.fh"
*
* return the combined number of ecp/so centers
*      
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom     !< [Input] the geometry handle
      integer soidin   !< [Input] the spin-orbit potential (so) handle
      integer ecpidin  !< [Input] the ECP handle
      integer num_cent !< [Output] the number of ECP/so centers
*::local
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_ecp         ! is tag in the ecp?
      logical inthe_so          ! is the tag in the so pot.?
*
      ecpso_ncent = geom_ncent(geom,ntotal)
      if (.not.ecpso_ncent) call errquit
     &    ('ecpso_ncent: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('ecpso_ncent: geom_cent_get failed',911, GEOM_ERR)
        inthe_ecp = bas_match_tags(gtag,ecpidin,tag_indx)
        inthe_so  = bas_match_tags(gtag,soidin,tag_indx)
        if (inthe_ecp.or.inthe_so) num_cent = num_cent + 1
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the number of centers with spin-orbit potentials
C> \return Return .true. if successfull, and .false. otherwise
      logical function so_ncent(geom,soidin,num_cent)
      implicit none
#include "errquit.fh"
*
* return the number of so centers
*      
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom     !< [Input] the geometry handle
      integer soidin   !< [Input] the spin-orbit potential (so) handle
      integer num_cent !< [Output] the number of so centers
*::local
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_so          ! is the tag in the so pot.?
*
      so_ncent = geom_ncent(geom,ntotal)
      if (.not.so_ncent) call errquit
     &    ('so_ncent: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('so_ncent: geom_cent_get failed',911, GEOM_ERR)
        inthe_so  = bas_match_tags(gtag,soidin,tag_indx)
        if (inthe_so) num_cent = num_cent + 1
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the number of centers with ECPs
C> \return Return .true. if successfull, and .false. otherwise
      logical function ecp_ncent(geom,ecpidin,num_cent)
      implicit none
#include "errquit.fh"
*
* return the number of ecp centers
*      
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom     !< [Input] the geometry handle
      integer ecpidin  !< [Input] the ECP handle
      integer num_cent !< [Output] the number of ECP centers
*::local
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_ecp         ! is the tag in the ecp pot.?
*
      ecp_ncent = geom_ncent(geom,ntotal)
      if (.not.ecp_ncent) call errquit
     &    ('ecp_ncent: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('ecp_ncent: geom_cent_get failed',911, GEOM_ERR)
        inthe_ecp  = bas_match_tags(gtag,ecpidin,tag_indx)
        if (inthe_ecp) num_cent = num_cent + 1
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the list of centers with either ECP or spin-orbit
C> potentials
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function ecpso_list_ncent(geom,soidin,ecpidin,
     &    num_cent,ecpso_list)
      implicit none
#include "errquit.fh"
*
* return the combined number of ecp/so centers and the list of centers
*      
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom          !< [Input] the geometry handle
      integer soidin        !< [Input] the spin-orbit potential handle
      integer ecpidin       !< [Input] the ECP handle
      integer num_cent      !< [Output] the number of ECP/so centers
      integer ecpso_list(*) !< [Output] the list of ECP/so centers
*::local
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_ecp         ! is tag in the ecp?
      logical inthe_so          ! is the tag in the so pot.?
*
      ecpso_list_ncent = geom_ncent(geom,ntotal)
      if (.not.ecpso_list_ncent) call errquit
     &    ('ecpso_list_ncent: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('ecpso_list_ncent: geom_cent_get failed',911, GEOM_ERR)
        inthe_ecp = bas_match_tags(gtag,ecpidin,tag_indx)
        inthe_so  = bas_match_tags(gtag,soidin,tag_indx)
        if (inthe_ecp.or.inthe_so) then
          num_cent = num_cent + 1
          ecpso_list(num_cent) = ic
        endif
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the list of centers with spin-orbit
C> potentials
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function so_list_ncent(geom,soidin,
     &    num_cent,so_list)
      implicit none
#include "errquit.fh"
*
* return the number of so centers and the list of centers
*      
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom       !< [Input] the geometry handle
      integer soidin     !< [Input] the spin-orbit potential (so) handle
      integer num_cent   !< [Output] the number of so centers
      integer so_list(*) !< [Output] the list of so centers
*::local
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_so          ! is the tag in the so pot.?
*
      so_list_ncent = geom_ncent(geom,ntotal)
      if (.not.so_list_ncent) call errquit
     &    ('so_list_ncent: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('so_list_ncent: geom_cent_get failed',911, GEOM_ERR)
        inthe_so  = bas_match_tags(gtag,soidin,tag_indx)
        if (inthe_so) then
          num_cent = num_cent + 1
          so_list(num_cent) = ic
        endif
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the list of centers with ECP 
C> potentials
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function ecp_list_ncent(geom,ecpidin,
     &    num_cent,ecp_list)
      implicit none
#include "errquit.fh"
*
* return the number of ecp centers and the list of centers
*      
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom        !< [Input] the geometry handle
      integer ecpidin     !< [Input] the ECP handle
      integer num_cent    !< [Output] the number of ECP centers
      integer ecp_list(*) !< [Output] the list of ECP centers
*::local
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_ecp         ! is tag in the ecp?
*
      ecp_list_ncent = geom_ncent(geom,ntotal)
      if (.not.ecp_list_ncent) call errquit
     &    ('ecp_list_ncent: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('ecp_list_ncent: geom_cent_get failed',911, GEOM_ERR)
        inthe_ecp = bas_match_tags(gtag,ecpidin,tag_indx)
        if (inthe_ecp) then
          num_cent = num_cent + 1
          ecp_list(num_cent) = ic
        endif
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the coordinates of all centers that have either
C> an ECP or spin-orbit potential
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      logical function ecpso_get_coords(geom,soidin,ecpidin,
     &    nxyz,xyz)
      implicit none
#include "errquit.fh"
*
* get the coordinates for the ecp/so centers
*
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom    !< [Input] the geometry handle
      integer soidin  !< [Input] the spin-orbit potential handle
      integer ecpidin !< [Input] the ecp handle
      integer nxyz    !< [Input] the size of xyz array
      double precision xyz(3,nxyz) !< [Output] the coordinates
*::local      
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer num_cent          ! number of centers returned
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_ecp         ! is tag in the ecp?
      logical inthe_so          ! is the tag in the so pot.?
*
      ecpso_get_coords = geom_ncent(geom,ntotal)
      if (.not.ecpso_get_coords) call errquit
     &    ('ecpso_get_coords: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('ecpso_get_coords: geom_cent_get failed',911, GEOM_ERR)
        inthe_ecp = bas_match_tags(gtag,ecpidin,tag_indx)
        inthe_so  = bas_match_tags(gtag,soidin,tag_indx)
        if (inthe_ecp.or.inthe_so) then
          num_cent = num_cent + 1
          if (num_cent.gt.nxyz) call errquit
     &        ('ecpso_get_coords: array passed in was too small',911,
     &       GEOM_ERR)
          xyz(1,num_cent) = cxyz(1)
          xyz(2,num_cent) = cxyz(2)
          xyz(3,num_cent) = cxyz(3)
        endif
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the coordinates of all centers that have 
C> a spin-orbit potential
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      logical function so_get_coords(geom,soidin,
     &    nxyz,xyz)
      implicit none
#include "errquit.fh"
*
* get the coordinates for the so centers
*
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom   !< [Input] the geometry handle
      integer soidin !< [Input] the spin-orbit potential handle
      integer nxyz   !< [Input] the size of xyz array 
      double precision xyz(3,nxyz) !< [Output] the coordinates
*::local      
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer num_cent          ! number of centers returned
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_so          ! is the tag in the so pot.?
*
      so_get_coords = geom_ncent(geom,ntotal)
      if (.not.so_get_coords) call errquit
     &    ('so_get_coords: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('so_get_coords: geom_cent_get failed',911, GEOM_ERR)
        inthe_so  = bas_match_tags(gtag,soidin,tag_indx)
        if (inthe_so) then
          num_cent = num_cent + 1
          if (num_cent.gt.nxyz) call errquit
     &        ('so_get_coords: array passed in was too small',911,
     &       GEOM_ERR)
          xyz(1,num_cent) = cxyz(1)
          xyz(2,num_cent) = cxyz(2)
          xyz(3,num_cent) = cxyz(3)
        endif
      enddo
      end
*----------------------------------------------------------------------
C> \brief Retrieve the coordinates of all centers that have
C> an ECP
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      logical function ecp_get_coords(geom,ecpidin,
     &    nxyz,xyz)
      implicit none
#include "errquit.fh"
*
* get the coordinates for the ecp centers
*
*::includes
#include "geom.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom    !< [Input] the geometry handle
      integer ecpidin !< [Input] the ecp handle
      integer nxyz    !< [Input] the size of xyz array
      double precision xyz(3,nxyz) !< [Output] the coordinates
*::local      
      integer ntotal            ! total number of centers in goemetry
      integer ic                ! loop counter for centers
      integer num_cent          ! number of centers returned
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
      logical inthe_ecp         ! is tag in the ecp?
*
      ecp_get_coords = geom_ncent(geom,ntotal)
      if (.not.ecp_get_coords) call errquit
     &    ('ecp_get_coords: geom_ncent failed?',911, GEOM_ERR)
      num_cent = 0
      do ic = 1,ntotal
        if (.not.geom_cent_get(geom,ic,gtag,cxyz,chg)) call errquit
     &      ('ecp_get_coords: geom_cent_get failed',911, GEOM_ERR)
        inthe_ecp = bas_match_tags(gtag,ecpidin,tag_indx)
        if (inthe_ecp) then
          num_cent = num_cent + 1
          if (num_cent.gt.nxyz) call errquit
     &        ('ecp_get_coords: array passed in was too small',911,
     &       GEOM_ERR)
          xyz(1,num_cent) = cxyz(1)
          xyz(2,num_cent) = cxyz(2)
          xyz(3,num_cent) = cxyz(3)
        endif
      enddo
      end
*----------------------------------------------------------------------
C> \brief Check whether a given center from the geometry is known to
C> the basis set by that name
C>
C> This routine looks the tag of the specified center up in the 
C> geometry and queries the basis set with that tag.
C>
C> \return Return .true. if the basis set knows of such centers, and
C> .false. otherwise
C
      logical function bas_isce(geom,basisin,center)
      implicit none
#include "errquit.fh"
*
* return true if the given center is in the designated basis set
*      
*::includes
#include "geom.fh"
#include "stdio.fh"
*::functions
      logical bas_match_tags
      external bas_match_tags
*::passed
      integer geom    !< [Input] the geometry handle
      integer basisin !< [Input] the basis set handle
      integer center  !< [Input] the lexical geometry index of center
*::local
      integer tag_indx          ! basis set unique center index for tag
      character*16 gtag         ! geometry tag 
      double precision cxyz(3)  ! coordinates
      double precision chg      ! charge
*
      bas_isce = geom_cent_get(geom,center,gtag,cxyz,chg)
      if (.not.bas_isce) call errquit
     &      ('bas_iscet: geom_cent_get failed',911, GEOM_ERR)
*      write(LuOut,*)' center = ',center,' gtag:<',gtag,'>'
      bas_isce  = bas_match_tags(gtag,basisin,tag_indx)
      end
C
C> \brief Retrieve the largest exponent any primitive function in the
C> basis set
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      logical function bas_large_exponent(basisin,exponent)
      implicit none
#include "errquit.fh"
c
c returns the largest primitive exponent in the given basis set
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c::functions
      logical  bas_check_handle
      external bas_check_handle
c::passed
      integer basisin           !< [Input] the basis set handle
      double precision exponent !< [Output] the maximum exponent
c::local
      integer basis
      integer num_cont
      integer icont, ucont
      integer myexptr
      integer mynprim
      integer iprim
      double precision test_exponent
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      bas_large_exponent =
     &    bas_check_handle(basisin,'bas_large_exponent')
      if (.not.bas_large_exponent) then
        write(luout,*)'invalid/inactive basis handle',basisin
        call errquit('bas_large_exponent: fatal error',911, BASIS_ERR)
      endif
      basis = basisin + Basis_Handle_Offset 
      num_cont = ncont_tot_gb(basis)
      exponent = -565.6589d00
      do icont = 1,num_cont
        ucont = sf_ibs_cn2ucn(icont,basis)
        myexptr = infbs_cont(CONT_IEXP,ucont,basis)
        mynprim = infbs_cont(CONT_NPRIM,ucont,basis)
        do iprim = 1, mynprim
          test_exponent = sf_exndcf((myexptr-1+iprim),basis)
          exponent = max(exponent,test_exponent)
*          write(luout,*)' exponent ',exponent,'  test ', test_exponent
        enddo
      enddo
*      write(luout,*)' largest exponent ',exponent
      end
C
C> \brief Checks whether the basis set contains any generally contracted 
C> or SP shells
C>
C> \return Return .true. if any generally contracted or SP shells are
C> present, and .false. otherwise
C
      logical function bas_any_gcorsp(basisin)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "basP.fh"
c
c returns true if the basis set has any sp functions or general contraction
c functions.
c
c::functions
      logical  bas_check_handle
      external bas_check_handle
c::passed
      integer basisin !< [Input] the basis set handle
c::local
      integer basis
      bas_any_gcorsp =
     &    bas_check_handle(basisin,'bas_any_gcorsp')
      if (.not.bas_any_gcorsp) then
        write(luout,*)'invalid/inactive basis handle',basisin
        call errquit('bas_any_gcorsp: fatal error',911, BASIS_ERR)
      endif
      basis = basisin + Basis_Handle_Offset 
      bas_any_gcorsp = bas_any_gc(basis).or.bas_any_sp_shell(basis)
      end
C
C> \brief Retrieves the largest exponent in the specified shell in the
C> basis set
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      logical function bas_cont_large_exponent(basisin,incont,exponent)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c::functions
      logical  bas_check_handle
      external bas_check_handle
c::passed
      integer basisin           !< [Input] the basis set handle
      integer incont            !< [Input] the shell rank
      double precision exponent !< [Output] the maximum exponent
c::local
      integer basis
      integer num_cont
      integer ucont
      integer myexptr
      integer mynprim
      integer iprim
      double precision test_exponent
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      bas_cont_large_exponent =
     &    bas_check_handle(basisin,'bas_cont_large_exponent')
      if (.not.bas_cont_large_exponent) then
        write(luout,*)'invalid/inactive basis handle',basisin
        call errquit('bas_cont_large_exponent: fatal error',911,
     &       BASIS_ERR)
      endif
      basis = basisin + Basis_Handle_Offset 
      num_cont = ncont_tot_gb(basis)
      exponent = -565.6589d00
      ucont = sf_ibs_cn2ucn(incont,basis)
      myexptr = infbs_cont(CONT_IEXP,ucont,basis)
      mynprim = infbs_cont(CONT_NPRIM,ucont,basis)
      do iprim = 1, mynprim
        test_exponent = sf_exndcf((myexptr-1+iprim),basis)
        exponent = max(exponent,test_exponent)
*        write(luout,*)' exponent ',exponent,'  test ', test_exponent
      enddo
*      write(luout,*)' largest exponent ',exponent
      end
*.....................................................................
C
C> \brief Checks whether a relativistic basis set exists on the RTDB
C>
C> This routine checks for a relativistic basis set by looking for
C> the name of the small component basis functions on the RTDB.
C> The name of the small component basis is pulled from a common block
C> and hence not visible to the caller.
C>
C> \return Return .true. if the basis set exists on the RTDB, and
C> .false. otherwise
c
      logical function bas_rel_exists (rtdb)
      implicit none
      integer rtdb !< [Input] the RTDB handle
#include "nwc_const.fh"
#include "rel_nwc.fh"
      logical  bas_name_exist_rtdb
      external bas_name_exist_rtdb
      bas_rel_exists = bas_name_exist_rtdb(rtdb,small_cmpt_name)
     &    .or. bas_name_exist_rtdb(rtdb,auto_small_cmpt_name)
      end
*.....................................................................
C> \brief Check whether a basis set is a (part of) a relativistic
C> basis set
C>
C> \return Return .true. if the basis set is a relativistic basis set,
C> and .false. otherwise
C
      logical function basis_is_rel (basisin)
      implicit none
      integer basisin !< [Input] the basis set handle
#include "nwc_const.fh"
#include "rel_nwc.fh"
      character*32 basis_name,trans_name
      logical  bas_name
      external bas_name
      basis_is_rel = bas_name(basisin,basis_name,trans_name)
      if (basis_is_rel) then
        basis_is_rel = (basis_name .eq. small_cmpt_name) 
     &      .or. (basis_name .eq. large_cmpt_name)
     &      .or. (basis_name .eq. auto_small_cmpt_name)
     &      .or. (basis_name .eq. auto_large_cmpt_name)
      end if
      end
*.....................................................................
C> \brief Retrieves the handle of the basis set with the name 
C> "ao basis"
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      logical function bas_get_ao_handle (basis,nbas,ao_handle)
      implicit none
      integer nbas        !< [Input] the number of basis sets
      integer basis(nbas) !< [Input] the array of basis set handles
      integer ao_handle   !< [Output] the ao basis set handle
c
c   Returns handle of "ao basis"
c
      integer ibas
      logical odum
      character*255 basis_name,trans_name
*
      logical bas_name
      external bas_name
c
      do ibas = 1,nbas
        odum = bas_name(basis(ibas),basis_name,trans_name)
        if (basis_name(1:8) .eq. 'ao basis') then
          bas_get_ao_handle = .true.
          ao_handle = basis(ibas)
          return
        end if
      end do
      bas_get_ao_handle = .false.
c
      return
      end
*.....................................................................
C> \brief Checks whether the basis set contains any relativistic shells
C>
C> \return Return .true. if there are any relativistic shell present,
C> and .false. otherwise
C
      logical function bas_any_rel_shells (basisin)
      implicit none
#include "errquit.fh"
c
      integer basisin !< [Input] the basis set handle
#include "nwc_const.fh"
#include "basdeclsP.fh"
#include "basP.fh"
c
c   Checks basis set to see if there are any relativistic shells
c
      integer ibas   ! basis set
      integer nucont ! number of unique contractions
      integer i, isum
c
      logical bas_check_handle
      external bas_check_handle
c
      if (.not. bas_check_handle(basisin,'bas_any_rel_shells'))
     &    call errquit('bas_any_rel_shells: invalid handle',99,
     &       BASIS_ERR)
      ibas = basisin+BASIS_HANDLE_OFFSET
      nucont = infbs_head(HEAD_NCONT,ibas)
      isum = 0
      do i = 1,nucont
        isum = isum+infbs_cont(CONT_RELLS,i,ibas)
      end do
      bas_any_rel_shells = isum .ne. 0
      return
      end
*.....................................................................
C> \brief Retrieve the highest angular momentum of all relativistic
C> shells in a basis set
C>
C> \return Return .true. if successfull, and .false. otherwise
C
      logical function bas_rel_high_ang(basisin,high_angular)
      implicit none
c
c  calculate and return highest angular momentem
c  for relativistic shells in given basis. 
c
#include "stdio.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin      !< [Input] the basis set handle
      integer high_angular !< [Output] the highest angular momentum
c:local
      integer basis, i, myang
      integer ucont
c
      bas_rel_high_ang = bas_check_handle(basisin,'bas_high_angular')
      if (.not. bas_rel_high_ang ) then
        write(LuOut,*) 'bas_rel_high_ang: basis handle not valid '
        return
      endif
c
      basis = basisin + Basis_Handle_Offset
      ucont =  infbs_head(head_ncont,basis)
      high_angular = -565
      do i = 1,ucont
        if (infbs_cont(cont_rells,i,basis) .ne. 0) then
          myang = abs(infbs_cont(cont_type,i,basis))
          high_angular = max(high_angular,myang)
        end if
      enddo
c
      bas_rel_high_ang = .true.
      return
      end
C
C> \brief Checks whether two tags match
C>
C> This routine matches two input tags using the geometry/basis set 
C> tags matching mechanism, in such a way that "H34" matches "H" or 
C> "hydrogen" or "h" if "H34" is one of the tags. "H" will however 
C> not match "h" i.e., case sensitivity is preserved.
C>
C> \return Return .true. if the tags match, and .false. otherwise
C
      logical function bas_do_tags_match(tag_one,tag_two)
      implicit none
#include "errquit.fh"
*
* This routine matches two input tags using the geometry/basis set 
* tags matching mechanism, in such a way that "H34" matches "H" or 
* "hydrogen" or "h" if "H34" is one of the tags.  "H" will however 
* not match "h" i.e., case sensitivity is preserved.
*
c::includes
#include "stdio.fh"
#include "inp.fh"
#include "nwc_const.fh"
c::functions
      logical geom_tag_to_element
      external geom_tag_to_element
c::passed
      character*(*) tag_one !< [Input] a geometry/basis tag
      character*(*) tag_two !< [Input] another geometry/basis tag
*     bas_do_tags_match     ! [output] true if they match
c::local
      character*16 one16      ! geometry or basis tags can only be *16
      character*16 two16      ! geometry or basis tags can only be *16
      integer lnotalpha       ! string index to non alpha character
      integer i               ! loop counter
      integer ind             ! dummy counter
      logical status          ! dummy storage
      integer len_one         ! length of tag one or substring of tag
      integer len_one_old     ! length of tag one
      integer len_two         ! length of tag two or substring of tag
      integer len_two_old     ! length of tag two
      character*2 one_sym     ! tag one -> symbol name
      character*16 one_elem   ! tag one -> element name
      integer one_atn         ! tag one -> atomic number
      character*2 two_sym     ! tag two -> symbol name
      character*16 two_elem   ! tag two -> element name
      integer two_atn         ! tag two -> atomic number
      logical debug           ! true for extra output
c
      integer na2z
      parameter (na2z = 26)
      character*1 a2z(na2z)
      data a2z /'a','b','c','d','e','f','g','h','i','j',
     &          'k','l','m','n','o','p','q','r','s','t',
     &          'u','v','w','x','y','z'/
c      
      debug = .true.
c
      bas_do_tags_match = .false.
c
      one16 = tag_one              !  copy tag to work on it
      two16 = tag_two              !  copy tag to work on it
c
* remove bq ghost info if it exists
c
      if (one16(1:4).eq.'bqgh') then
        one16 = one16(5:)
      endif
      if (two16(1:4).eq.'bqgh') then
        two16 = two16(5:)
      endif
c
c... first match full geometry tag to full basis tag list 
c          did the user specifically assign a tag ?

      if (one16.eq.two16) then
        bas_do_tags_match = .true.
        goto 00009
      endif
c
c... Now check to see if either tag has non alpha chracters. 
c...  the substring prior to the first non-alpha character is the users idea
c...  of the "name" of the tag.  
c
      len_one = inp_strlen(one16)
      len_one_old = len_one
      lnotalpha = len_one+1
      do i = 1,len_one
        if (.not.(inp_match(na2z,.false.,
     &      one16(i:i),a2z,ind))) then ! compare character to alpha (case-less test)
          lnotalpha = i
          goto 00001
        endif
      enddo
00001 continue
      do i = lnotalpha,len_one
        one16(i:i) = ' '
      enddo
      if (debug) write(luout,*)' the string:', tag_one,
     &    'has the user substring of ',one16

      len_two = inp_strlen(two16)
      len_two_old = len_two
      lnotalpha = len_two+1
      do i = 1,len_two
        if (.not.(inp_match(na2z,.false.,
     &      two16(i:i),a2z,ind))) then ! compare character to alpha (case-less test)
          lnotalpha = i
          goto 00002
        endif
      enddo
00002 continue
      do i = lnotalpha,len_two
        two16(i:i) = ' '
      enddo
      if (debug) write(luout,*)' the string:', tag_two,
     &    'has the user substring of ',two16
c
c... match user substrings 
      len_one = inp_strlen(one16)
      len_two = inp_strlen(two16)
      if ((len_one_old.gt.len_one).or.(len_two_old.gt.len_two)) then
        if (one16.eq.two16) then
          bas_do_tags_match = .true.
          goto 00009
        endif
      endif
c
c... now get symbol and element names for each tag and match those
c    use the substring for each tag to match these
c... bq's with basis functions must match from the above rules!
c
      status = geom_tag_to_element(one16,one_sym,one_elem,one_atn)
      if (.not.status)then
        if (.not.(one_sym.eq.'bq')) then
          write(luout,*)'tag <',one16,
     &        '> could not be matched to an element symbol'
          call errquit('bas_do_tags_match: fatal error ',911, BASIS_ERR)
        endif
      endif
      if (one_sym.eq.'bq') goto 00009   ! bq labels with basis functions must match from above
      status = geom_tag_to_element(two16,two_sym,two_elem,two_atn)
      if (.not.status)then
        if (.not.(two_sym.eq.'bq')) then
          write(luout,*)'tag <',two16,
     &        '> could not be matched to an element symbol'
          call errquit('bas_do_tags_match: fatal error ',911, BASIS_ERR)
        endif
      endif
      if (two_sym.eq.'bq') goto 00009   ! bq labels with basis functions must match from above
*
      if (one_elem.eq.two_elem) then
        bas_do_tags_match = .true.
        goto 00009
      else if (one_sym.eq.two_sym) then
        bas_do_tags_match = .true.
        goto 00009
      endif
      if (debug)
     &    write(luout,*)'bas_do_tags_match:debug: no match for tags <',
     &    tag_one,'> and <',tag_two,'>'
      return  ! no match
c
00009 continue
      if (debug) then
        write(luout,10000)tag_one,tag_two
      endif
      return
10000 format('bas_do_tags_match:debug: tag ',a16,
     &    ' matched this tag ',a16)
      end
C> @}
