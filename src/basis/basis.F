      block data basis_data
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
      implicit none
#include "basP.fh"
c     
      data nbasis_rtdb /0/
      data bsactive /nbasis_bsmx*.false./
      data angular_bs /nbasis_bsmx*-565/
      data bsversion /1.00d00/
c     
      end
      logical function bas_version()
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
      idata = idata + nbasis_bsmx*4 + 1
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
      logical function bas_create(basis,name)
c
c creates a handle and marks it active in the in-core data structure
c
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "basdeclsP.fh"
c::passed
c:integer  basis [output] returned handle
c:char*(*) name  [input] name of basis set.  
c
      integer basis
      character*(*)name
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
c:integer basisin [input]
      integer basisin
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
c
      logical function bas_check_handle(basisin,msg)
c
c Checks to see if a basis set handle is valid
c
      implicit none
#include "basP.fh"
c:passed
c:Integer       basisin  [input] handle
c:Character*(*) msg      [input] error message
c
      integer basisin
      character*(*) msg
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
c
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
c:I rtdb [input] valid rtdb handle
c.
      integer rtdb
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
c
c
c
      logical function bas_print(basisin)
c
c routine to print unique basis information that is in core
c
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "inp.fh"
c
c function declarations
c      
      logical bas_check_handle
      external bas_check_handle
c:: passed
      integer basisin
c:: local
      integer mytags, myucont, myprim, mycoef, basis, len_tag, empty
      integer i,j,k,l, ifcont, mygen, mytype, iexptr, icfptr
      character*3 ctype(0:6),cltype(2)
      character*3 shell_type
      character*16 blank
c
      blank   = ' '
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

*      write(6,2) mytags, myucont, myprim, mycoef  
* 2    format(
*     $     ' number of unique tags         :',i5/
*     $     ' number of unique contractions :',i5/
*     $     ' number of unique primitives   :',i5/
*     $     ' number of unique coefficients :',i5)

      do 00100 i=1,mytags
         len_tag = inp_strlen(bs_tags(i,basis))
         empty = (16 - len_tag)/2
         write(6,3) blank(1:empty), bs_tags(i,basis)
 3       format(' ',a,a/' ----------------'/)
         myucont = infbs_tags(TAG_NCONT,i,basis)
*         write(6,4) myucont,infbs_tags(TAG_NPRIM,i,basis),
*     $        infbs_tags(TAG_NCOEF,i,basis)
* 4      format(
*     &       ' number of contractions:',i5/
*     &       ' number of primitives  :',i5/
*     &       ' number of coefficients:',i5/)
c        
        ifcont = infbs_tags(TAG_FCONT,i,basis)
c     
        write(6,6)
 6      format(
     $       '          Exponent Coefficients '/
     $       '       ----------- ',60('-'))
        do 00200 j=1,myucont
          myprim = infbs_cont(CONT_NPRIM,ifcont,basis)
          mygen  = infbs_cont(CONT_NGEN,ifcont,basis)
          
          mytype = infbs_cont(CONT_TYPE, ifcont, basis)
          if (mytype.lt.0) then
            shell_type = cltype(abs(mytype))
          else
            shell_type = ctype(mytype)
          endif
*          write(6,5) j, shell_type(1:inp_strlen(shell_type)),
*     $         myprim, mygen
* 5        format(/
*     $         ' contraction           :', i5/
*     $         ' type                  :    ',a/
*     $         ' number of primitives  :', i5/
*     $         ' number of contractions:', i5/)
          iexptr = infbs_cont(CONT_IEXP,ifcont,basis) - 1
          icfptr = infbs_cont(CONT_ICFP,ifcont,basis) - 1
          do 00300 k=1,myprim
             write(6,7) j, shell_type(1:2), exndcf((iexptr+k),basis),
     &             (exndcf((icfptr+k+(l-1)*myprim),basis),l=1,mygen)
 7           format(1x,i2,1x,a2,1x,f11.4,20f10.6)
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
c
      logical function bas_rtdb_load(rtdb, geom, basisin, name)
      implicit none
c
#include "rtdb.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "mafdecls.h"
#include "inp.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
c
c function declarations
c
c::function
c:bas
      logical bas_check_handle
      external bas_check_handle
c
c
c:: arguments
c:I   rtdb    ->  valid rtdb handle      
c:I   geom    ->  valid geometry handle with info loaded
c:I   basisin ->  valid basis handle
c:C** name    ->  basis set name that must exist on the rtdb
c:: passed  
      integer rtdb, geom, basisin
      character*(*) name
c:: local
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
     &       rtdb_par_cget(rtdb, tmp, ntags_bsmx, bs_tags(1,basis))
c
      tmp(lentmp:) = ':exps and coeffs'
      rtdb_status = rtdb_status .and.
     &       rtdb_par_get(
     &       rtdb, tmp, mt_dbl, mxbs_exndcf, exndcf(1,basis))
c
      tmp(lentmp:) = ':header'
      rtdb_status = rtdb_status .and.
     &       rtdb_par_get(
     &       rtdb, tmp, mt_int, ndbs_head, infbs_head(1,basis))
c
      tmp(lentmp:) = ':tags info'
      rtdb_status = rtdb_status .and.
     &       rtdb_par_get(
     &       rtdb, tmp, mt_int,
     &       ndbs_tags*ntags_bsmx, infbs_tags(1,1,basis))
c
      tmp(lentmp:) = ':contraction info'
      rtdb_status = rtdb_status .and.
     &       rtdb_par_get(
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
c:fucked      do 00100 i = 1,nat
c:fucked        if (.not.inp_match(infbs_head(HEAD_NTAGS,basis),.true.,
c:fucked     &         tags(i,geom),
c:fucked     &         bs_tags(1,basis),ibs_ce2uce(i,basis))) then
c:fucked          write(6,*)' geom tag was not found in basis tag list rtdb'
c:fucked          write(6,*)' geom tag searched for ',tags(i,geom)
c:fucked          write(6,*)' basis tag list was from basis ',bs_name(basis)
c:fuckedc....  add further diagnostics later          
c:fucked        endif
c:fucked#if defined(BS_DEBUG)
c:fucked        write(6,*)' ibs_ce2uce(',i,') = ',ibs_ce2uce(i,basis)
c:fucked#endif
c:fucked00100 continue
c
c build center to unique center map
c
      do 00100 i=1,nat
        foundit = .false.
        do 00101 j = 1,infbs_head(HEAD_NTAGS,basis)
          if(inp_compare(.true.,tags(i,geom),bs_tags(j,basis))) then
            ibs_ce2uce(i,basis) = j
            foundit = .true.
            goto 00102
          endif
00101   continue
        if (.not.foundit) then
          write(6,*)' geom tag was not found in basis tag list rtdb'
          write(6,*)' geom tag searched for ',tags(i,geom)
          write(6,*)' basis tag list was from basis ',bs_name(basis)
        endif
00102   continue
00100 continue
c
c build total # of contractions and center -> contraction range map 
c      
      idum_cont     = infbs_tags(TAG_NCONT,ibs_ce2uce(1,basis),basis)
      ncont_tot_gb(basis)  = idum_cont
      ibs_ce2cnr(1,1,basis) = 1
      ibs_ce2cnr(2,1,basis) = idum_cont
      do 00200 i=2,nat
        idum_cont    = infbs_tags(TAG_NCONT,ibs_ce2uce(i,basis),basis)
        ncont_tot_gb(basis) = idum_cont + ncont_tot_gb(basis)
        ibs_ce2cnr(1,i,basis) = ibs_ce2cnr(2,(i-1),basis) + 1
        ibs_ce2cnr(2,i,basis) = ibs_ce2cnr(1,i,basis) + idum_cont - 1
00200 continue
c      
c build contraction -> center map
c
      do 00300 i=1,nat
        jstart = ibs_ce2cnr(1,i,basis)
        jend   = ibs_ce2cnr(2,i,basis)
        do 00400 j=jstart,jend
          ibs_cn2ce(j,basis) = i
00400   continue
00300 continue
c
c build contraction -> unique contraction map
c
      do 00500 i=1,nat
        idum_at = ibs_ce2uce(i,basis)
        jstart = ibs_ce2cnr(1,i,basis)
        jend   = ibs_ce2cnr(2,i,basis)
        jsize  = jend - jstart + 1
        kstart = infbs_tags(TAG_FCONT,idum_at,basis)
        kend   = infbs_tags(TAG_LCONT,idum_at,basis)
        ksize  = kend - kstart + 1
        lsize  = infbs_tags(TAG_NCONT,idum_at,basis)
        if (jsize.eq.ksize.and.ksize.eq.lsize) then
          icount = 0
          do 00600 j=jstart,jend
            ibs_cn2ucn(j,basis) = kstart + icount
            icount = icount + 1
00600     continue
        else
          write(6,*)' bas_rtdb_load: ERROR '
          write(6,*)' contraction range size <-> unique contraction ',
     &           'range size mismatch'
          write(6,*)'        contraction range (',jstart,':',jend,')'
          write(6,*)' unique contraction range (',kstart,':',kend,')'
          write(6,*)'                   contraction size: ',jsize
          write(6,*)' calculated unique contraction size: ',ksize
          write(6,*)'     lookup unique contraction size: ',lsize
          bas_rtdb_load = .false.
          return
        endif
00500 continue
c
c build nprim_tot_gb, nbf_tot_gb, and 
c contraction -> basis function range map 
c
      iu_cont = ibs_cn2ucn(1,basis)
      nbf     = infbs_cont(CONT_TYPE,iu_cont,basis)
      nbf     = (nbf+1)*(nbf+2)/2
      nbf     = nbf*infbs_cont(CONT_NGEN,iu_cont,basis)

      ibs_cn2bfr(1,1,basis) = 1
      ibs_cn2bfr(2,1,basis) = nbf

      nbf_tot_gb(basis)   = nbf
      nprim_tot_gb(basis) = infbs_cont(CONT_NPRIM,iu_cont,basis)
      do 00700 i = 2,ncont_tot_gb(basis)
        iu_cont = ibs_cn2ucn(i,basis)

        nbf     = infbs_cont(CONT_TYPE,iu_cont,basis)
        nbf     = (nbf+1)*(nbf+2)/2
        nbf     = nbf*infbs_cont(CONT_NGEN,iu_cont,basis)

        ibs_cn2bfr(1,i,basis) = ibs_cn2bfr(2,(i-1),basis) + 1
        ibs_cn2bfr(2,i,basis) = ibs_cn2bfr(1,i,basis) + nbf - 1

        nbf_tot_gb(basis) = nbf_tot_gb(basis) + nbf
        nprim_tot_gb(basis) = nprim_tot_gb(basis) +
     &         infbs_cont(CONT_NPRIM,iu_cont,basis)
00700 continue
c
c build high angular momentum of this loaded <basis|geom> pair
c
      do 00800 i = 1, ncont_tot_gb(basis)
        iu_cont = ibs_cn2ucn(i,basis)
        myang   = infbs_cont(CONT_TYPE,iu_cont,basis)
        angular_bs(basis) = max(angular_bs(basis),myang)
00800 continue
      return
      end
c
      logical function bas_rtdb_store(rtdb, name, basisin)
      implicit none
#include "basdeclsP.fh"
#include "basP.fh"
      integer rtdb              ! [input] handle to database
      character*(*) name        ! [input] name to use when storing db
      integer basisin           ! [input] handle to basis set
      logical bas_check_handle, bas_rtdb_do_store
      external bas_check_handle, bas_rtdb_do_store
c
c     Store basis set (not geometry) related info about specified
c     basis in into the rtdb with the given name
c
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
      logical function bas_rtdb_do_store(rtdb, name, tagsin, head_array,
     &       tags_array, ucont_array, excfin, ntagsin, nucontin,
     &       nexcf)
      implicit none
#include "basdeclsP.fh"
#include "basP.fh"
#include "inp.fh"
#include "mafdecls.h"
#include "rtdb.fh"
c
c     This routine stores the basis set information in the appropriate 
c     data structure on the run-time-data-base (rtdb).  
c
c     This is a private routine called by the user level routine 
c     bas_rtdb_store(rtdb, name, basis)
c
c::: argument description
c:C*(*)             name         -> name of the basis set
c:C*16(1:ntagsin)   tagsin       -> name of each tag
c::: functions
      logical bas_rtdb_in
      logical bas_rtdb_add
c::: passed
      integer rtdb
      integer ntagsin, nucontin, nexcf
      character*(*) name
      character*16 tagsin(ntagsin)
      integer head_array(ndbs_head)
      integer tags_array(ndbs_tags,ntagsin)
      integer ucont_array(ndbs_ucont,nucontin)
      double precision excfin(nexcf)
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
      status = status .and. rtdb_par_cput(rtdb,tmp,ntagsin,tagsin)
      
      tmp(lentmp:) = ':exps and coeffs'
      status = status .and. rtdb_par_put(rtdb,tmp,mt_dbl,nexcf,excfin)

      tmp(lentmp:) = ':header'
      status = status .and.
     &       rtdb_par_put(rtdb,tmp,mt_int,ndbs_head,head_array)

      tmp(lentmp:) = ':tags info'
      status = status .and. rtdb_par_put(
     &       rtdb,tmp,mt_int,(ndbs_tags*ntagsin), tags_array)

      tmp(lentmp:) = ':contraction info'
      status = status .and. rtdb_par_put(
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
      logical function bas_high_angular(basisin,high_angular)
      implicit none
c
c  calculate and store high angular momentem function
c   for given basis. 
c
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      logical bas_check_handle
      external bas_check_handle
c:: passed
c  basis        [input]   ! basis set handle
c  high_angular [output]  ! high angular momentum of basis
c
      integer basisin, high_angular
c:local
      integer basis, myucont, i
c
      bas_high_angular = bas_check_handle(basisin,"bas_high_angular")
      if (.not. bas_high_angular ) then
        write(6,*) ' basis handle not valid '
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
        high_angular = max(high_angular, infbs_cont(CONT_TYPE,i,basis))
00100 continue
c
      angular_bs(basis) = high_angular
      bas_high_angular = .true.
      return
      end
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
      integer basisin
c::local
      integer basis

      gbs_map_clear = bas_check_handle(basisin,"gbs_map_clear")
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
      logical function gbs_map_print(basisin)
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin
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
c
      logical function bas_rtdb_out(rtdb)
      implicit none
#include "basP.fh"
#include "rtdb.fh"
#include "mafdecls.h"
c     
c::passed
      integer rtdb
c     
c     output to rtdb info about known basis sets
c     
      bas_rtdb_out  =
     &     rtdb_par_put(rtdb, 'basis:nbasis',
     &       MT_INT, 1, nbasis_rtdb)
     &     .and.
     &     rtdb_par_cput(rtdb, 'basis:names', nbasis_rtdb,
     &       bs_names_rtdb)
      if (.not. bas_rtdb_out) 
     &     write(6,*) ' bas_rtdb_out: rtdb is corrupt '
c     
      end
      logical function bas_rtdb_add(rtdb, name)
      implicit none
#include "basP.fh"
#include "rtdb.fh"
#include "mafdecls.h"
#include "inp.fh"
c     
      integer rtdb              ! [input]
      character*(*) name        ! [input]
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
      logical function bas_print_all()
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
      subroutine bas_err_info(info)
      implicit none
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
      logical function bas_rtdb_in(rtdb)
      implicit none
#include "basP.fh"
#include "rtdb.fh"
#include "mafdecls.h"
#include "inp.fh"
c     
      integer rtdb              ! [input]
      integer bas
c     
c     load in info about known basis sets ... this is more
c     for diagnostic and debugging purposes
c     
      bas_rtdb_in = .false.
      nbasis_rtdb = 0
      if (rtdb_par_get(rtdb, 'basis:nbasis', MT_INT, 1, nbasis_rtdb))
     $     then
        if (.not. rtdb_par_cget(rtdb,'basis:names', nbasis_rtdb_mx,
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
      logical function bas_cn2ce(basisin,cont,center)
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin, cont, center
c::local
      integer basis
c
      bas_cn2ce = bas_check_handle(basisin)
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
      logical function bas_cn2bfr(basisin,cont,ifirst,ilast)
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin, cont, ifirst, ilast
c::local 
      integer basis
c
      bas_cn2bfr = bas_check_handle(basisin)
      if(.not.bas_cn2bfr) return
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
      logical function bas_ce2cnr(basisin,center,ifirst,ilast)
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "geom.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin, center, ifirst, ilast
c::local 
      integer basis, nat
c
      bas_ce2cnr = bas_check_handle(basisin)
      if(.not.bas_ce2cnr) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
      bas_ce2cnr = geom_ncent(ibs_geom(basis),nat)
      if (nat.eq.0.or..not.bas_ce2cnr) then
        write(6,*)' bas_ce2cnr: ERROR '
        write(6,*)' number of centers is zero or weird'
        write(6,*)' nat = ',nat
c..... add diagnostics later
        return
      endif

      bas_ce2cnr = center.gt.0 .and. center.le.nat
      if (.not.bas_ce2cnr) then
        write(6,*)' bas_ce2cnr: invalid center information '
        write(6,*)' contraction range is 1:',nat
        write(6,*)' input contraction was: ',center
        return
      endif
c
      return
      end
      logical function bas_get_exponent(basisin,icont,exp)
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
      integer basisin, icont 
      double precision exp(*)
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
c
      logical function bas_continfo(basisin,icont,
     &       nprimo,ngeno,sphcart)
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin, icont, nprimo, ngeno, sphcart
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
      nprimo  = infbs_cont(CONT_NPRIM,myucont,basis)
      ngeno   = infbs_cont(CONT_NGEN,myucont,basis)
      bas_continfo=.true.
      return
      end
c
      logical function bas_numcont(basisin,numcont)
      implicit none
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
      bas_numcont = bas_check_handle(basisin,'bas_numcont')
      if (.not.bas_numcont) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      numcont = ncont_tot_gb(basis)

      bas_numcont = .true.
      return
      end
c
      logical function bas_numbf(basisin,nbf)
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin,nbf
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
c      
      logical function bas_get_coeff(basisin,icont,coeff)
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
      integer basisin, icont 
      double precision coeff(*)
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
      logical function bas_set_exponent(basisin,icont,exp,nexp)
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
      integer basisin, icont, nexp
      double precision exp(nexp)
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
c
      logical function bas_set_coeff(basisin,icont,coeff,ncoeff)
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
      integer basisin, icont, ncoeff
      double precision coeff(ncoeff)
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

