      logical function bas_add_ucnt(basis, tag, l_value, ngen, nprim, 
     $     expnt, coeffs, ldc)
c $Id: basisP.F,v 1.6 1994-05-05 17:44:42 d3e129 Exp $
      implicit none
#include "basdeclsP.fh"
#include "basP.fh"
      integer basis             ! [input] basis handle
      character*(*) tag         ! [input] tag on which to contraction
      integer l_value           ! [input] type of contraction (s, p, ..., sp)
      integer ngen              ! [input] no. of contractions
      integer nprim             ! [input] no. of primitives
      integer ldc               ! [input] leading dimension of coeffs
      double precision expnt(nprim) ! [input] exponents
      double precision coeffs(ldc, 1:*) ! [input] coefficients
c
      integer ind               ! Index into basis function structures
      integer free              ! Free space pointer
      integer i, itag, jtag, iu_cont, ntags ! Locals
c
      logical bas_add_utag
      external bas_add_utag
c
      logical bas_check_handle
      external bas_check_handle
c
c     adds a new general contraction on the specified tag.  If the
c     tag is not present it will also add that by calling bas_add_utag
c      
      bas_add_ucnt = bas_check_handle(basis,'bas_add_ucnt')
      if (.not. bas_add_ucnt) return
      ind = basis + BASIS_HANDLE_OFFSET
c
c     Make sure that the tag is in the list
c
      bas_add_ucnt = bas_add_utag(basis, tag, itag)
      if (.not. bas_add_ucnt) return
c
c     Update header information about all unique contractions on all tags
c     Free points to next free word in the exndcf
c
      free = infbs_head(HEAD_NPRIM,ind) + infbs_head(HEAD_NCOEF,ind) + 1
      if ((free+nprim*ngen+nprim-1) .gt. mxbs_exndcf) then
         write(6,*) 'bas_add_ucnt: too many prims/coeffs'
         bas_add_ucnt = .false.
         return
      endif
      if (infbs_head(HEAD_NCONT,ind)+1 .gt. nucont_bsmx) then
         write(6,*) 'bas_add_ucnt: too many contractions '
         bas_add_ucnt = .false.
         return
      endif
c         
      infbs_head(HEAD_NCONT,ind) = infbs_head(HEAD_NCONT,ind) + 1
      infbs_head(HEAD_NPRIM,ind) = infbs_head(HEAD_NPRIM,ind) + nprim
      infbs_head(HEAD_NCOEF,ind) = infbs_head(HEAD_NCOEF,ind) +
     $     ngen*nprim
c
      ntags = infbs_head(HEAD_NTAGS,ind)
      if (itag .ne. ntags) then
         do jtag = ntags, itag+1, -1
c
c     Shuffle data+pointers for following tags up one contraction
c
            do iu_cont = infbs_tags(TAG_LCONT,jtag,ind),
     $           infbs_tags(TAG_FCONT,jtag,ind), -1
               do i = 1, ndbs_ucont
                  infbs_cont(i,iu_cont+1,ind) =
     $                 infbs_cont(i,iu_cont,ind)
               enddo
            enddo
c
c     Increment first and last contractions on following tags
c
            infbs_tags(TAG_FCONT,jtag,ind) =
     $           infbs_tags(TAG_FCONT,jtag,ind) + 1
            infbs_tags(TAG_LCONT,jtag,ind) =
     $           infbs_tags(TAG_LCONT,jtag,ind) + 1
         enddo
      endif
c
c     Increment basis info on this tag
c
      infbs_tags(TAG_NCONT,itag,ind) = infbs_tags(TAG_NCONT,itag,ind)
     $     + 1
      infbs_tags(TAG_NPRIM,itag,ind) = infbs_tags(TAG_NPRIM,itag,ind)
     $     + nprim
      infbs_tags(TAG_NCOEF,itag,ind) = infbs_tags(TAG_NCOEF,itag,ind)
     $     + nprim*ngen
      if (infbs_tags(TAG_FCONT,itag,ind).eq.0) then
         if (itag .ne. ntags) call errquit
     $        ('bas_add_ucnt: tag error', itag)
         infbs_tags(TAG_FCONT,itag,ind) = infbs_head(HEAD_NCONT,ind)
         infbs_tags(TAG_LCONT,itag,ind) = infbs_head(HEAD_NCONT,ind)
      else
         infbs_tags(TAG_LCONT,itag,ind) = infbs_tags(TAG_LCONT,itag,ind)
     $        + 1
      endif
c
      iu_cont = infbs_tags(TAG_LCONT,itag,ind) ! Index of new contraction
c
      infbs_cont(CONT_TYPE, iu_cont,ind) = l_value
      infbs_cont(CONT_NPRIM,iu_cont,ind) = nprim
      infbs_cont(CONT_NGEN, iu_cont,ind) = ngen
      infbs_cont(CONT_TYPE, iu_cont,ind) = l_value
      infbs_cont(CONT_IEXP, iu_cont,ind) = free
      infbs_cont(CONT_ICFP, iu_cont,ind) = free + nprim
c
c     Copy real data over
c
      call dcopy(nprim, expnt, 1, exndcf(free,ind), 1)
      free = free + nprim
      do i = 1, ngen
         call dcopy(nprim, coeffs(1,i), 1, exndcf(free,ind), 1)
         free = free + nprim
      enddo
*      write(6,*) ' expnt input ', nprim
*      call output(expnt, 1, nprim, 1, 1, nprim, 1, 1)
*      write(6,*) ' coeffs input ', nprim, ngen, ldc
*      call output(coeffs, 1, nprim, 1, ngen, ldc, ngen, 1)
c
c     Done
c
      end
      logical function bas_add_utag(basisin, tag, itag)
      implicit none
#include "basdeclsP.fh"
#include "basP.fh"
#include "inp.fh"
      integer basisin           ! [input] basis handle
      character*(*) tag         ! [input] name of tag
      integer itag              ! [output] index of tag
c
      integer basis             ! [local] index into basis arrays
      logical bas_check_handle
      external bas_check_handle
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
         if (inp_compare(.true.,bs_tags(itag,basis),tag)) return
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
c
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
c
      logical function bas_getu_coeff(basisin,icont,coeff)
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
      bas_getu_coeff = bas_check_handle(basisin,'bas_getu_coeff')
      if (.not.bas_getu_coeff) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = infbs_head(HEAD_NCONT,basis)
      myucont  = icont
c
      bas_getu_coeff = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_getu_coeff)) then
        write(6,*)' bas_getu_coeff: ERROR '
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
      bas_getu_coeff = .true.
c
      return
      end
c
c
      logical function bas_getu_exponent(basisin,icont,exp)
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
      bas_getu_exponent = bas_check_handle(basisin,'bas_getu_exponent')
      if (.not.bas_getu_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      icontmax = infbs_head(HEAD_NCONT,basis)
      myucont = icont

      bas_getu_exponent = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_getu_exponent)) then
        write(6,*)' bas_getu_exponent: ERROR '
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
      bas_getu_exponent = .true.
c
      return
      end
c
      logical function bas_setu_coeff(basisin,icont,coeff,ncoeff)
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
      bas_setu_coeff = bas_check_handle(basisin,'bas_setu_coeff')
      if (.not.bas_setu_coeff) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      icontmax = infbs_head(HEAD_NCONT,basis)
      myucont  = icont
c
      bas_setu_coeff = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_setu_coeff)) then
        write(6,*)' bas_setu_coeff: ERROR '
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
      bas_setu_coeff = ncoeff .eq. (myprim*mygen)
      if(.not.bas_setu_coeff) then
        write(6,*)' bas_setu_coeff: ERROR '
        write(6,*)' input and stored number of coefficients ',
     &         '(nprim*ngen) differ '
        write(6,*)' input  nprim*ngen: ',ncoeff
        write(6,*)' stored nprim*ngen: ',(myprim*mygen)
        return
      endif
      call dcopy(ncoeff,coeff,1,exndcf(mycoeffptr,basis),1)
c
      bas_setu_coeff = .true.
c
      return
      end
      logical function bas_setu_exponent(basisin,icont,exp,nexp)
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
      bas_setu_exponent = bas_check_handle(basisin,'bas_setu_exponent')
      if (.not.bas_setu_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      icontmax = infbs_head(HEAD_NCONT,basis)
      myucont = icont

      bas_setu_exponent = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_setu_exponent)) then
        write(6,*)' bas_setu_exponent: ERROR '
        write(6,*)' contraction range for basis is 1:',
     &         icontmax
        write(6,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      bas_setu_exponent = myprim.eq.nexp
      if (.not.bas_setu_exponent) then
        write(6,*)' bas_setu_exponent: ERROR '
        write(6,*)' input and stored number of exponents ',
     &         '(nprim) differ '
        write(6,*)' input  nprim: ',nexp
        write(6,*)' stored nprim: ',myprim
        return
      endif
      call dcopy(nexp,exp,1,exndcf(myexptr,basis),1)
c
      bas_setu_exponent = .true.
c
      return
      end
c
      logical function bas_ucontinfo(basisin,icont,
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
c... no spherical yet 3/94 only cart. 
      sphcart = 0
      nprimo  = infbs_cont(CONT_NPRIM,myucont,basis)
      ngeno   = infbs_cont(CONT_NGEN,myucont,basis)
      bas_ucontinfo=.true.
      return
      end
c
      logical function bas_unumcont(basisin,numcont)
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
      bas_unumcont = bas_check_handle(basisin,'bas_numcont')
      if (.not.bas_unumcont) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      numcont = infbs_head(HEAD_NCONT,basis)

      bas_unumcont = .true.
      return
      end
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
      integer function nbf_from_ucont(ucont,basisin)
c
c function that returns the number of basis functions in a contraction
c 
c  types 0->S, 1->P, 2->D, 3->F etc. -1->SP -2->SPD
c
      implicit none
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
     &  call errquit('nbf_from_ucont: bad basis handle',basisin)     
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      type = infbs_cont(CONT_TYPE,ucont,basis)
      ngen = infbs_cont(CONT_NGEN,ucont,basis)
      if (type.ge.0) then 
        nbf_from_ucont = ngen*(type+1)*(type+2)/2
      else if (type.eq.-1) then
        nbf_from_ucont = ngen*4/2
      else if (type.eq.-2) then
        nbf_from_ucont = ngen*10/3
      else 
        call errquit('nbf_from_ucont: bad cont type',type)
      endif
      end
