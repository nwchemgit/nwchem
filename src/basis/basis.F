      logical function bas_version()
c $Id: basis.F,v 1.31 1994-11-10 17:58:58 d3g681 Exp $
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
      rdata = mxbs_exndcf*nbasis_bsmx + 1
      rdata = 8*rdata
c
c.. integer data in basis set object common
      idata =         ndbs_head*nbasis_bsmx
      idata = idata + ndbs_tags*ntags_bsmx*nbasis_bsmx
      idata = idata + ndbs_ucont*nucont_bsmx*nbasis_bsmx
      idata = idata + nbasis_bsmx*5 + 1
      idata = idata + nbasis_rtdb_mx
      idata = 4*idata
c
c.. integer data in the mapped object.  
      mapidata =            4*ncont_mx*nbasis_bsmx
      mapidata = mapidata + 3*nat_mx*nbasis_bsmx
      mapidata = mapidata + 4*nbasis_bsmx
      mapidata = 4*mapidata
c
c.. total space 
      total4 = idata + mapidata
      total8 = 2*total4 + rdata + cdata 
      total4 = total4 + rdata + cdata 
c
      write(6,'(////1x,a,f5.2,a)')
     &       ' **** basis set version ',bsversion,' ****'
      write(6,'(1x,a,i20,a)')
     &       '   character data in-core ',cdata,' bytes'
      write(6,'(1x,a,i20,a)')
     &       '   real      data in-core ',rdata,' bytes'
      write(6,'(1x,a,i20,a)')
     &       '   integer*4 data in-core ',idata,' bytes'
      write(6,'(1x,a,i20,a)')
     &       'or integer*8 data in-core ',(2*idata),' bytes'
      write(6,'(1x,a,i20,a)')
     &       '   integer*4 mapping data in-core ',
     &       mapidata,' bytes'
      write(6,'(1x,a,i20,a/)')
     &       'or integer*8 mapping data in-core ',
     &       (2*mapidata),' bytes'
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
      write(6,'(1x,a,f5.2,a)')
     &       ' **** basis set version ',bsversion,' ****'
      write(6,'(1x,a,i20,a)')
     &       '   character data in-core ',cdata,' Kbytes'
      write(6,'(1x,a,i20,a)')
     &       '   real      data in-core ',rdata,' Kbytes'
      write(6,'(1x,a,i20,a)')
     &       '   integer*4 data in-core ',idata,' Kbytes'
      write(6,'(1x,a,i20,a)')
     &       'or integer*8 data in-core ',(2*idata),' Kbytes'
      write(6,'(1x,a,i20,a)')
     &       '   integer*4 mapping data in-core ',
     &       mapidata,' Kbytes'
      write(6,'(1x,a,i20,a/)')
     &       'or integer*8 mapping data in-core ',
     &       (2*mapidata),' Kbytes'
      write(6,*)' total(4)   = ',total4,' Kbytes'
      write(6,*)' total(8)   = ',total8,' Kbytes'
      write(6,'(////)')
c
      bas_version = .true.
      end
c......................................................................
      logical function bas_create(basis,name)
c
c creates a handle and marks it active in the in-core data structure
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "basdeclsP.fh"
#include "util.fh"
c::functions
c ifill from util
c dfill from util
c::passed
c
      integer basis       ! [output] returned handle
      character*(*)name   ! [input] name of basis set.  
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
      call dfill(mxbs_exndcf, 0.0d0, exndcf(1,basis), 1)
      call ifill(ndbs_tags*ntags_bsmx, 0,
     $     infbs_tags(1,1,basis), 1)
      call ifill(ndbs_ucont*nucont_bsmx, 0,
     $     infbs_cont(1,1,basis), 1)
c
c Initialize geo-basis info to empty
c
      call ifill(ncont_mx, 0, ibs_cn2ucn(1,basis), 1)
      call ifill(ncont_mx, 0, ibs_cn2ce (1,basis), 1)
      call ifill(nat_mx,   0, ibs_ce2uce(1,basis), 1)
      call ifill(2*ncont_mx,0,ibs_cn2bfr(1,1,basis), 1)
      call ifill(2*nat_mx, 0, ibs_ce2cnr(1,1,basis), 1)
      ncont_tot_gb(basis) = 0
      nprim_tot_gb(basis) = 0
      nbf_tot_gb(basis)   = 0
      ibs_geom(basis)     = 0
c
c Mark basis as active and return info
c
      bsactive(basis) = .true.
      basis = basis - BASIS_HANDLE_OFFSET 
      bas_create = .true.
c
      end
c......................................................................
      logical function bas_destroy(basisin)
c
c destroys information about an active incore basis
c and the associated mapping arrays.
c
      implicit none
#include "basP.fh"
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

      bas_destroy = bas_check_handle(basisin,'bas_destroy')
      if (.not. bas_destroy) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      if(.not.gbs_map_clear(basisin)) then
        write(6,*)' error clearing map '
        bas_destroy = .false.
        return
      endif

      bsactive(basis) = .false.
      bas_destroy = .true.
c
      return
      end
c......................................................................
      logical function bas_check_handle(basisin,msg)
c
c Checks to see if a basis set handle is valid
c
      implicit none
#include "basP.fh"
c:passed
      integer basisin    ! [input] handle
      character*(*) msg  ! [input] error message
c::local
      integer basis
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_check_handle = basis.gt.0 .and. basis.le.nbasis_bsmx
      if (bas_check_handle)
     &       bas_check_handle = bas_check_handle .and. bsactive(basis)
c
      if (.not. bas_check_handle) then
        write(6,*)msg,': basis handle is invalid '
        write(6,*)'basis_check_handle: lexical handle ',basis
        write(6,*)'basis_check_handle:         handle ',basisin
      endif
      return
      end
c......................................................................
      logical function bas_321g_load(rtdb)
      implicit none
c
c routine to load the rtdb with 321g basis sets for atoms 1-20
c H to Ca
#include "basP.fh"
#include "basdeclsP.fh"
#
c
c function declarations
c
c::function
      logical bas_create
      logical bas_print
      logical bas_rtdb_do_store
      logical bas_destroy
      external bas_print
      external bas_create
      external bas_rtdb_do_store
      external bas_destroy
c:: passed
      integer rtdb  ! [input] valid rtdb handle
c:: local
      integer basis,usebas, itag, icont
c
      if (.not.bas_create(basis,'321g:1-20')) then
        write(6,*)' basis set handle not created '
        bas_321g_load = .false.
        return
      endif
c
      usebas = basis + BASIS_HANDLE_OFFSET
      bs_name(usebas) = '321g:1-20'
c
c.. only h and o for now
c
c
c.. hydrogen and oxygen
c
c
      bs_tags(1,usebas) = 'H'
      bs_tags(2,usebas) = 'He'
      bs_tags(3,usebas) = 'O'
      infbs_head(HEAD_NTAGS,usebas) = 3
      infbs_head(HEAD_NCONT,usebas) = 9
      infbs_head(HEAD_NPRIM,usebas) = 30
      infbs_head(HEAD_NCOEF,usebas) = 30
c.. hydrogen
      itag = 1
      infbs_tags(TAG_NCONT,itag,usebas) = 2
      infbs_tags(TAG_NPRIM,itag,usebas) = 3
      infbs_tags(TAG_NCOEF,itag,usebas) = 3
      infbs_tags(TAG_FCONT,itag,usebas) = 1
      infbs_tags(TAG_LCONT,itag,usebas) = 2
c.. hydrogen cont 1 bas 1
      icont = 1
      infbs_cont(CONT_TYPE ,icont,usebas)  = 0
      infbs_cont(CONT_NPRIM,icont,usebas)  = 2
      infbs_cont(CONT_NGEN ,icont,usebas)  = 1
      infbs_cont(CONT_IEXP ,icont,usebas)  = 1
      infbs_cont(CONT_ICFP ,icont,usebas)  = 3
      infbs_cont(CONT_TAG  ,icont,usebas)  = itag
      exndcf(1,usebas) =  5.44717800d00
      exndcf(2,usebas) =  0.82454700d00
      exndcf(3,usebas) =  0.15628500d00
      exndcf(4,usebas) =  0.90469100d00
c.. hydrogen cont 2 bas 2
      icont = icont + 1
      infbs_cont(CONT_TYPE ,icont,usebas)  = 0
      infbs_cont(CONT_NPRIM,icont,usebas)  = 1
      infbs_cont(CONT_NGEN ,icont,usebas)  = 1
      infbs_cont(CONT_IEXP ,icont,usebas)  = 5
      infbs_cont(CONT_ICFP ,icont,usebas)  = 6
      infbs_cont(CONT_TAG  ,icont,usebas)  = itag
      exndcf(5,usebas) =  0.18319200d00
      exndcf(6,usebas) =  1.00000000d00
c.. helium 
      itag = itag + 1
      icont = icont + 1
      infbs_tags(TAG_NCONT,itag,usebas) = 2
      infbs_tags(TAG_NPRIM,itag,usebas) = 3
      infbs_tags(TAG_NCOEF,itag,usebas) = 3
      infbs_tags(TAG_FCONT,itag,usebas) = icont
      infbs_tags(TAG_LCONT,itag,usebas) = icont + 2 - 1
c.. helium cont 1 bas 3
      infbs_cont(CONT_TYPE,icont,usebas)  = 0
      infbs_cont(CONT_NPRIM,icont,usebas) = 2
      infbs_cont(CONT_NGEN,icont,usebas)  = 1
      infbs_cont(CONT_IEXP,icont,usebas)  = 7
      infbs_cont(CONT_ICFP,icont,usebas)  = 9
      infbs_cont(CONT_TAG,icont,usebas)   = itag
      exndcf(7,usebas) = 13.62670000
      exndcf(8,usebas) =  1.99935000
      exndcf(9,usebas) =  0.17523000
      exndcf(10,usebas)=  0.89348300
c.. helium cont 2 bas 4
      icont = icont + 1
      infbs_cont(CONT_TYPE,icont,usebas)  = 0
      infbs_cont(CONT_NPRIM,icont,usebas) = 1
      infbs_cont(CONT_NGEN,icont,usebas)  = 1
      infbs_cont(CONT_IEXP,icont,usebas)  = 11
      infbs_cont(CONT_ICFP,icont,usebas)  = 12
      infbs_cont(CONT_TAG,icont,usebas)   = itag
      exndcf(11,usebas) =  0.38299300
      exndcf(12,usebas) =  1.00000000
c.. oxygen
      itag = itag + 1
      icont = icont + 1
      infbs_tags(TAG_NCONT,itag,usebas) = 5
      infbs_tags(TAG_NPRIM,itag,usebas) = 9
      infbs_tags(TAG_NCOEF,itag,usebas) = 9
      infbs_tags(TAG_FCONT,itag,usebas) = icont
      infbs_tags(TAG_LCONT,itag,usebas) = icont + 5 - 1
c.. oxygen cont 1 bas 5
      infbs_cont(CONT_TYPE,icont,usebas)  = 0
      infbs_cont(CONT_NPRIM,icont,usebas) = 3
      infbs_cont(CONT_NGEN,icont,usebas)  = 1
      infbs_cont(CONT_IEXP,icont,usebas)  = 13
      infbs_cont(CONT_ICFP,icont,usebas)  = 16
      infbs_cont(CONT_TAG,icont,usebas)   = itag
      exndcf(13,usebas)  =  322.03700000d00
      exndcf(14,usebas)  =   48.43080000d00
      exndcf(15,usebas)  =   10.42060000d00
      exndcf(16,usebas) =    0.05923940d00
      exndcf(17,usebas) =    0.35150000d00
      exndcf(18,usebas) =    0.70765800d00
c.. oxygen cont 2 bas 6
      icont = icont + 1
      infbs_cont(CONT_TYPE,icont,usebas)  = 0
      infbs_cont(CONT_NPRIM,icont,usebas) = 2
      infbs_cont(CONT_NGEN,icont,usebas)  = 1
      infbs_cont(CONT_IEXP,icont,usebas)  = 19
      infbs_cont(CONT_ICFP,icont,usebas)  = 21
      infbs_cont(CONT_TAG,icont,usebas)   = 2
      exndcf(19,usebas) =   7.40294000d00
      exndcf(20,usebas) =   1.57620000d00
      exndcf(21,usebas) =  -0.40445300d00
      exndcf(22,usebas) =   1.22156000d00
c.. oxygen cont 3 bas 7      
      icont = icont + 1
      infbs_cont(CONT_TYPE,icont,usebas)  = 1
      infbs_cont(CONT_NPRIM,icont,usebas) = 2
      infbs_cont(CONT_NGEN,icont,usebas)  = 1
      infbs_cont(CONT_IEXP,icont,usebas)  = 23
      infbs_cont(CONT_ICFP,icont,usebas)  = 25
      infbs_cont(CONT_TAG,icont,usebas)   = 2
      exndcf(23,usebas) =   7.40294000d00
      exndcf(24,usebas) =   1.57620000d00
      exndcf(25,usebas) =   0.24458600d00
      exndcf(26,usebas) =   0.85395500d00
c.. oxygen cont 4 bas 8
      icont = icont + 1
      infbs_cont(CONT_TYPE,icont,usebas)  = 0
      infbs_cont(CONT_NPRIM,icont,usebas) = 1
      infbs_cont(CONT_NGEN,icont,usebas)  = 1
      infbs_cont(CONT_IEXP,icont,usebas)  = 27
      infbs_cont(CONT_ICFP,icont,usebas)  = 28
      infbs_cont(CONT_TAG,icont,usebas)   = 2
      exndcf(27,usebas) = 0.37368400d00
      exndcf(28,usebas) = 1.00000000d00
c.. oxygen cont 5 bas 9
      icont = icont + 1
      infbs_cont(CONT_TYPE,icont,usebas)  = 1
      infbs_cont(CONT_NPRIM,icont,usebas) = 1
      infbs_cont(CONT_NGEN,icont,usebas)  = 1
      infbs_cont(CONT_IEXP,icont,usebas)  = 29
      infbs_cont(CONT_ICFP,icont,usebas)  = 30
      infbs_cont(CONT_TAG,icont,usebas)   = 2
      exndcf(29,usebas) = 0.37368400d00
      exndcf(30,usebas) = 1.00000000d00
c
      bas_321g_load = bas_rtdb_do_store(rtdb,
     &       bs_name(usebas), bs_tags(1,usebas), 
     &       infbs_head(1,usebas), infbs_tags(1,1,usebas),
     &       infbs_cont(1,1,usebas),exndcf(1,usebas),3,9,30)
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
c......................................................................
      logical function bas_print(basisin)
c
c routine to print unique basis information that is in core
c
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
#include "geom.fh"
c
c function declarations
c      
      logical bas_check_handle
      external bas_check_handle
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
      character*35 buffer  ! Room for tag+space+(+element+) = 16+1+1+16+1
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
      write(6,1)bs_name(basis)(1:inp_strlen(bs_name(basis))), 
     $     bs_trans(basis)(1:inp_strlen(bs_trans(basis)))
 1    format('                      Basis "',a,'" -> "',a,'"'/
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
     $       '          Exponent  Coefficients '/
     $       '       ------------ ',60('-'))
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
             write(6,7) j, shell_type(1:2), exndcf((iexptr+k),basis),
     &             (exndcf((icfptr+k+(l-1)*myprim),basis),l=1,mygen)
 7           format(1x,i2,1x,a2,1x,f12.6,20f10.6)
00300     continue
          write(6,*)
          ifcont = ifcont + 1
00200   continue
00100 continue
c
c  If geom is set print out the info about total basis info associated with
c  the geometry also
c
c  ... not done yet
c
      return
      end
c......................................................................
      logical function bas_rtdb_load(rtdb, geom, basisin, name)
      implicit none
c
c routine that loads a basis set from the rtdb and using the 
c geometry information builds the mapping arrays to contractions/
c shells, basis functions, and centers.
c
#include "rtdb.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "mafdecls.fh"
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
      integer nbf_from_ucont
      logical bas_check_handle
      external nbf_from_ucont
      external bas_check_handle
c:: passed  
      integer rtdb        ! [input] rtdb handle      
      integer geom        ! [input] geometry handle with info loaded
      integer basisin     ! [input] basis handle
      character*(*) name  ! [input] basis set name that must be on the rtdb
c:: local
      character*16 element
      character*2  symbol
      integer atn
      integer lentmp, basis
      character*256 tmp
      logical rtdb_status, status
      integer i,nat, idum_cont, idum_at
      integer j,jstart,jend, jsize
      integer kstart,kend,ksize,lsize, icount
      integer nbf, iu_cont, myang
      logical foundit
c:: initalize local
      rtdb_status = .true.
c
c
c check geom and basis handles returns false if either is invalid
c
      bas_rtdb_load = geom_check_handle(geom,'bas_rtdb_load')
      if (.not.bas_rtdb_load) return
      bas_rtdb_load = bas_check_handle(basisin,'bas_rtdb_load')
      if (.not.bas_rtdb_load) return
c
c store geom tag with basis map info
c
      basis = basisin + BASIS_HANDLE_OFFSET
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
      tmp(lentmp:) = ':exps and coeffs'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_dbl, mxbs_exndcf, exndcf(1,basis))
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
        write(6,*)' bas_rtdb_load: ERROR '
        write(6,*)' basis set is not there '
        write(6,*)' name requested            <|',bs_name(basis),'|>'
        write(6,*)' translated name requested <|',bs_trans(basis),'|>'
        bas_rtdb_load = .false.
c..... add diagnostics later
        return
      endif
c
      status = geom_ncent(geom, nat)
      if (nat.eq.0.or..not.status) then
        write(6,*)' bas_rtdb_load: ERROR '
        write(6,*)' number of centers is zero or weird'
        write(6,*)' nat = ',nat
        bas_rtdb_load = .false.
c..... add diagnostics later
        return
      endif
c
c build center to unique center map
c
      if (nat.gt.nat_mx) then
        write(6,*)' nat     = ',nat
        write(6,*)' nat max = ',nat_mx
        call errquit ('bas_rtdb_load: nat.gt.nat_mx',911)
      endif
      do 00100 i=1,nat
        foundit = .false.
        do 00101 j = 1,infbs_head(HEAD_NTAGS,basis)
          if(tags(i,geom) .eq. bs_tags(j,basis)) then
            ibs_ce2uce(i,basis) = j
            foundit = .true.
            goto 00102
          endif
00101   continue
        if (.not. foundit .and. ga_nodeid().eq.0) then
           if (geom_tag_to_element(tags(i,geom), symbol, element,
     $          atn)) then
              write(6,10) tags(i,geom)(1:inp_strlen(tags(i,geom))),
     $             element(1:inp_strlen(element)),
     $             bs_name(basis)(1:inp_strlen(bs_name(basis)))
 10           format(/' ERROR: geometry tag ',a,' (',a,') is an atom ',
     $             'but has no functions in basis "',a,'"'/
     $             ' ERROR: only bq* centers can have no functions')
              call util_flush(6)
              call errquit('bas_rtdb_load: basis/geometry mismatch', 0)
           else
              write(6,11) i, tags(i,geom)(1:inp_strlen(tags(i,geom))),
     $             bs_name(basis)(1:inp_strlen(bs_name(basis)))
 11           format(/'WARNING: geometry tag ',i4, ' ', a,
     $             ' not found in basis "',a,'"'/)
              call util_flush(6)
           endif
        endif
00102   continue
00100 continue
c
c build total # of contractions and center -> contraction range map 
c      
      ncont_tot_gb(basis)  = 0
      do 00200 i=1,nat
         if (ibs_ce2uce(i,basis).gt.0) then
            idum_cont = infbs_tags(TAG_NCONT,ibs_ce2uce(i,basis),basis)
            ibs_ce2cnr(1,i,basis) = ncont_tot_gb(basis) + 1
            ibs_ce2cnr(2,i,basis) = ncont_tot_gb(basis) + idum_cont
            ncont_tot_gb(basis)   = idum_cont + ncont_tot_gb(basis)
         else
            ibs_ce2cnr(1,i,basis) = 0 ! No functions on this center
            ibs_ce2cnr(2,i,basis) = -1
         endif
00200 continue
c
      if (ncont_tot_gb(basis) .eq. 0) call errquit
     $    ('bas_rtdb_load: no functions in basis set', 0)
      if (ncont_tot_gb(basis) .gt. ncont_mx) then
        write(6,*)' number of contractions     = ',
     &      ncont_tot_gb(basis)
        write(6,*)' number of contractions max = ',ncont_mx
        call errquit ('bas_rtdb_load: ncont.gt.ncont_mx ',911)
      endif
c      
c build contraction -> center map
c
      do 00300 i=1,nat
         if (ibs_ce2uce(i,basis).gt.0) then
            jstart = ibs_ce2cnr(1,i,basis)
            jend   = ibs_ce2cnr(2,i,basis)
            do 00400 j=jstart,jend
               ibs_cn2ce(j,basis) = i
00400       continue
         endif
00300 continue
c
c build contraction -> unique contraction map
c
      do 00500 i=1,nat
         jstart = ibs_ce2cnr(1,i,basis)
         jend   = ibs_ce2cnr(2,i,basis)
         jsize  = jend - jstart + 1
         if (jsize .gt. 0) then
            idum_at = ibs_ce2uce(i,basis)
            kstart = infbs_tags(TAG_FCONT,idum_at,basis)
            kend   = infbs_tags(TAG_LCONT,idum_at,basis)
            ksize  = kend - kstart + 1
            lsize  = infbs_tags(TAG_NCONT,idum_at,basis)
            if (jsize.eq.ksize.and.ksize.eq.lsize) then
               icount = 0
               do 00600 j=jstart,jend
                  ibs_cn2ucn(j,basis) = kstart + icount
                  icount = icount + 1
00600          continue
            else
               write(6,*)' bas_rtdb_load: ERROR '
               write(6,*)' contraction range size mismatch'
               write(6,*)'        cont. range (',jstart,':',jend,')'
               write(6,*)' unique cont. range (',kstart,':',kend,')'
               write(6,*)'        cont. size: ',jsize
               write(6,*)' calculated unique cont. size: ',ksize
               write(6,*)'     lookup unique cont. size: ',lsize
               bas_rtdb_load = .false.
               return
            endif
         endif
00500 continue
c
c build nprim_tot_gb, nbf_tot_gb, and 
c contraction -> basis function range map 
c find nbfmax for basis (initialized in block data statement)
c
      nbf_tot_gb(basis)   = 0
      nprim_tot_gb(basis) = 0
      do 00700 i = 1,ncont_tot_gb(basis)
        iu_cont = ibs_cn2ucn(i,basis)
c
        nbf = nbf_from_ucont(iu_cont,basisin)
        nbfmax_bs(basis) = max(nbfmax_bs(basis),nbf)
c
        ibs_cn2bfr(1,i,basis) = nbf_tot_gb(basis) + 1
        ibs_cn2bfr(2,i,basis) = nbf_tot_gb(basis) + nbf

        nbf_tot_gb(basis) = nbf_tot_gb(basis) + nbf
        nprim_tot_gb(basis) = nprim_tot_gb(basis) +
     &         infbs_cont(CONT_NPRIM,iu_cont,basis)
00700 continue
c
c build high angular momentum of this loaded <basis|geom> pair
c note angular_bs(*) initialized in block data function
c
      do 00800 i = 1, infbs_head(HEAD_NCONT,basis)
        myang   = abs(infbs_cont(CONT_TYPE,i,basis))
        angular_bs(basis) = max(angular_bs(basis),myang)
00800 continue
      return
      end
c......................................................................
      logical function bas_rtdb_store(rtdb, name, basisin)
      implicit none
c
c routine that does an incore basis set store to the rtdb
c
#include "basdeclsP.fh"
#include "basP.fh"
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
c
      bas_rtdb_store = bas_check_handle(basisin,'bas_rtdb_store')
      if (.not. bas_rtdb_store) return
      basis = basisin + BASIS_HANDLE_OFFSET
c
      bas_rtdb_store = bas_rtdb_do_store(rtdb, name, 
     $     bs_tags(1,basis), infbs_head(1,basis),
     $     infbs_tags(1,1,basis),
     $     infbs_cont(1,1,basis), exndcf(1,basis),
     $     infbs_head(HEAD_NTAGS,basis),
     $     infbs_head(HEAD_NCONT,basis), 
     $     infbs_head(HEAD_NPRIM,basis)+infbs_head(HEAD_NCOEF,basis))
c
      end
c......................................................................
      logical function bas_rtdb_do_store(rtdb, name, tagsin, head_array,
     &       tags_array, ucont_array, excfin, ntagsin, nucontin,
     &       nexcf)
      implicit none
c
c stores from argument data structures the basis set information to
c the rtdb. 
c
#include "basdeclsP.fh"
#include "basP.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
c     This routine stores the basis set information in the appropriate 
c     data structure on the run-time-data-base (rtdb).  
c
c     This is a private routine called by the user level routine 
c     bas_rtdb_store(rtdb, name, basis)
c
c::: functions
      logical bas_rtdb_in
      logical bas_rtdb_add
      external bas_rtdb_in
      external bas_rtdb_add
c::: passed
      integer rtdb                              ! [input] rtdb handle
      character*(*) name                        ! [input] name of basis set
      integer ntagsin                           ! [input] number of tags
      integer nucontin                          ! [input] number of unique contractions
      integer nexcf                             ! [input] number of exponents and 
      character*16 tagsin(ntagsin)              ! [input] name of tags
      integer head_array(ndbs_head)             ! [input] head data
      integer tags_array(ndbs_tags,ntagsin)     ! [input] tag data
      integer ucont_array(ndbs_ucont,nucontin)  ! [input] unique contraction data
      double precision excfin(nexcf)            ! [input] exponents and contraction coeffs
c                                               !         contractions coeffs.
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
c......................................................................
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
c......................................................................
      logical function gbs_map_clear(basisin)
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
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

      gbs_map_clear = bas_check_handle(basisin,'gbs_map_clear')
      if (.not. gbs_map_clear ) then
        write(6,*) ' basis handle not valid '
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
      call ifill(ncont_mx,0,ibs_cn2ucn(1,basis),1)
      call ifill(ncont_mx,0,ibs_cn2ce(1,basis),1)
      call ifill(2*ncont_mx,0,ibs_cn2bfr(1,1,basis),1)
      call ifill(nat_mx,0,ibs_ce2uce(1,basis),1)
      call ifill(2*nat_mx,0,ibs_ce2cnr(1,1,basis),1)
      ncont_tot_gb(basis) = 0
      nprim_tot_gb(basis) = 0
      nbf_tot_gb(basis)   = 0
c
      angular_bs(basis) = -565
c
      gbs_map_clear = .true.
      return
      end
c......................................................................
      logical function gbs_map_print(basisin)
      implicit none
c
c prints the basis set <-> geometry mapping information
c
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin  ! [input] basis set handle
c::local
      integer nat, basis, i, myfirst, mylast, mysize, mycenter, myucont
      integer mygeom
      logical status
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
     &         'center:',i,'maps to unique center:',ibs_ce2uce(i,basis)
        myfirst = ibs_ce2cnr(1,i,basis)
        mylast  = ibs_ce2cnr(2,i,basis)
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
        mycenter = ibs_cn2ce(i,basis)
        myucont  = ibs_cn2ucn(i,basis)
        myfirst  = ibs_cn2bfr(1,i,basis)
        mylast   = ibs_cn2bfr(2,i,basis)
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
c......................................................................
      logical function bas_rtdb_out(rtdb)
c     
c     output to rtdb info about known basis sets
c     
      implicit none
#include "basP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
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
c......................................................................
      logical function bas_rtdb_add(rtdb, name)
      implicit none
c
c add basis set name to known basis set list on rtdb
c
#include "basP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
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
         write(6,*) ' bas_rtdb_add: too many basetries on rtdb ', name
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
c......................................................................
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
c......................................................................
      logical function bas_rtdb_in(rtdb)
c     
c     load in info about known basis sets ... this is more
c     for diagnostic and debugging purposes
c     
      implicit none
#include "basP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
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
c......................................................................
      logical function bas_cn2ce(basisin,cont,center)
      implicit none
c
c returns the center for a given mapped contraction
c
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin   ! [input] basis set handle
      integer cont      ! [input] mapped contraction index
      integer center    ! [output] center index
c::local
      integer basis
c
      bas_cn2ce = bas_check_handle(basisin,'bas_cn2ce')
      if(.not.bas_cn2ce) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_cn2ce = cont.gt.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2ce) then
        write(6,*)' bas_cn2ce: invalid contraction information '
        write(6,*)' contraction range is 1:',ncont_tot_gb(basis)
        write(6,*)' input contraction was: ',cont
        return
      endif
      center = ibs_cn2ce(cont,basis)
c
      return
      end
c......................................................................
      logical function bas_cn2bfr(basisin,cont,ifirst,ilast)
c
c returns the first basis function index of a mapped contraction 
c in ifirst and the last basis function index in ilast
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
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
      bas_cn2bfr = .true.
#ifdef BASIS_DEBUG      
      bas_cn2bfr = bas_check_handle(basisin,'bas_cn2bfr')
      if(.not.bas_cn2bfr) return
#endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_cn2bfr = cont.gt.0 .and. cont.le.ncont_tot_gb(basis)
      if (.not.bas_cn2bfr) then
        write(6,*)' bas_cn2bfr: invalid contraction information '
        write(6,*)' contraction range is 1:',ncont_tot_gb(basis)
        write(6,*)' input contraction was: ',cont
        return
      endif
c
      ifirst = ibs_cn2bfr(1,cont,basis)
      ilast  = ibs_cn2bfr(2,cont,basis)
c
      return
      end
c......................................................................
      logical function bas_ce2bfr(basis, icent, ibflo, ibfhi)
c
c  returns the basis function range for a given center
c
      implicit none
      integer basis             ! [input] handle
      integer icent             ! [input] no. of center
      integer ibflo, ibfhi      ! [output] range of functions on center
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
c......................................................................
      logical function bas_ce2cnr(basisin,center,ifirst,ilast)
c
c returns the mapped contraction range on a given center
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "geom.fh"
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
      ifirst = ibs_ce2cnr(1, center, basis)
      ilast  = ibs_ce2cnr(2, center, basis)
c
      return
      end
c......................................................................
      logical function bas_bf2ce(basisin,testbf,center)
c
c routine to return the center of a given basis function 
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "geom.fh"
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
        ibflo = ibs_cn2bfr(1,ibs_ce2cnr(1,iat,basis),basis)
        ibfhi = ibs_cn2bfr(2,ibs_ce2cnr(2,iat,basis),basis)
        if (testbf.ge.ibflo.and.testbf.le.ibfhi) then
          center = iat
          return
        endif
00100 continue
c
      end
c......................................................................
      logical function bas_bf2cn(basisin,testbf,cont)
c
c returns the mapped contraction index that contains the given 
c basis function index
c 
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
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
        ibflo = ibs_cn2bfr(1,icont,basis)
        ibfhi = ibs_cn2bfr(2,icont,basis)
        if(testbf.ge.ibflo.and.testbf.le.ibfhi) then
          cont = icont
          return
        endif
00100 continue
c
      end
c......................................................................
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
c......................................................................
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
      integer       basisin    ! [input] basis set handle
      character*(*) basis_name ! [output] basis set name when loaded
      character*(*) trans_name ! [output] basis set name in context when loaded
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
c......................................................................
      logical function bas_get_exponent(basisin,icont,exp)
c
c  returns the exponents of a given mapped contraction 
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
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
c
      bas_get_exponent = bas_check_handle(basisin,'bas_get_exponent')
      if (.not.bas_get_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 
      icontmax = ncont_tot_gb(basis)
      myucont = ibs_cn2ucn(icont,basis)
      bas_get_exponent = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_get_exponent)) then
        write(6,*)' bas_get_exponent: ERROR '
        write(6,*)' contraction range for basis is 1:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      call dcopy(myprim,exndcf(myexptr,basis),1,exp,1)
c
      bas_get_exponent = .true.
c
      return
      end
c......................................................................
      logical function bas_continfo(basisin,icont,
     &       type,nprimo,ngeno,sphcart)
c
c  returns the generic information about the given mapped contraction
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
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
      if (.not.(icont.gt.0.and.icont.le.icontmax)) then
        write(6,*)' bas_continfo: ERROR '
        write(6,*)' contraction range for basis is 1:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        bas_continfo = .false.
        return
      endif
c
      myucont = ibs_cn2ucn(icont,basis)
c... no spherical yet 3/94 only cart. 
      sphcart = 0
      type    = infbs_cont(CONT_TYPE, myucont,basis)
      nprimo  = infbs_cont(CONT_NPRIM,myucont,basis)
      ngeno   = infbs_cont(CONT_NGEN, myucont,basis)
      bas_continfo=.true.
      return
      end
c......................................................................
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
c......................................................................
      logical function bas_get_coeff(basisin,icont,coeff)
c
c returns the coefficients for the given mapped contraction
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
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
c
      bas_get_coeff = bas_check_handle(basisin,'bas_get_coeff')
      if (.not.bas_get_coeff) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = ncont_tot_gb(basis)
      myucont = ibs_cn2ucn(icont,basis)
      bas_get_coeff = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_get_coeff)) then
        write(6,*)' bas_get_coeff: ERROR '
        write(6,*)' contraction range for basis is 1:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      mycoeffptr = infbs_cont(CONT_ICFP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      mygen   = infbs_cont(CONT_NGEN,myucont,basis)
      call dcopy((myprim*mygen),exndcf(mycoeffptr,basis),1,coeff,1)
c
      bas_get_coeff = .true.
c
      return
      end
c......................................................................
      logical function bas_set_exponent(basisin,icont,exp,nexp)
c
c sets the exponents to the input array for a given mapped contraction
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin             ! [input] basis set handle        
      integer icont               ! [input] mapped contraction index
      integer nexp                ! [input] number of exponents for contraction
      double precision exp(nexp)  ! [input] "new" exponents for contraction    
c::local
      integer basis, myucont, icontmax
      integer myprim,myexptr
c
      bas_set_exponent = bas_check_handle(basisin,'bas_set_exponent')
      if (.not.bas_set_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      icontmax = ncont_tot_gb(basis)
      myucont = ibs_cn2ucn(icont,basis)

      bas_set_exponent = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_set_exponent)) then
        write(6,*)' bas_set_exponent: ERROR '
        write(6,*)' contraction range for basis is 1:',
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
      call dcopy(nexp,exp,1,exndcf(myexptr,basis),1)
c
      bas_set_exponent = .true.
c
      return
      end
c......................................................................
      logical function bas_set_coeff(basisin,icont,coeff,ncoeff)
c
c sets the coefficients of the given mapped contraction
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin                 ! [input] basis set handle                   
      integer icont                   ! [input] mapped contraction index           
      integer ncoeff                  ! [input] number of coeffs. for contraction
      double precision coeff(ncoeff)  ! [input] "new" coeffs. for contraction    
c::local
      integer basis, myucont, icontmax
      integer mycoeffptr, myprim, mygen
c
      bas_set_coeff = bas_check_handle(basisin,'bas_set_coeff')
      if (.not.bas_set_coeff) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = ncont_tot_gb(basis)
      myucont = ibs_cn2ucn(icont,basis)

      bas_set_coeff = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_set_coeff)) then
        write(6,*)' bas_set_coeff: ERROR '
        write(6,*)' contraction range for basis is 1:',
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
      call dcopy(ncoeff,coeff,1,exndcf(mycoeffptr,basis),1)
c
      bas_set_coeff = .true.
c
      return
      end
c......................................................................
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
      integer basisin       ! [input] basis set handle
      integer nbf_max       ! [output] largest basis function block on any 
c                           !          shell in the given basis set.
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
c......................................................................
      logical function bas_nbf_ce_max(basisin,nbf_max)
      implicit none
c
c  calculate, return and store maximum basis function block size
c   for all contractions in a given basis. 
c
#include "basP.fh"
#include "geom.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
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
        mylo  = ibs_cn2bfr(1,ibs_ce2cnr(1,iat,basis),basis)
        myhi  = ibs_cn2bfr(2,ibs_ce2cnr(2,iat,basis),basis)
        mynbf = myhi - mylo + 1
        nbf_max = max(nbf_max, mynbf)
00100 continue
c
      return
      end
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
      logical function bas_nprim_cn_max(basisin, nprim_max)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
      integer basisin   ! [input] basis set handle
      integer nprim_max ! [output] maximum number of prims in any contraction
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
