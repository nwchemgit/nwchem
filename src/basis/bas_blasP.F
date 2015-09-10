*
* routines from basisP.F that use blas calls
*
* $Id$
*
*.....................................................................
      logical function bas_add_ucnt(basis, tag, l_value, ngen, nprim, 
     $     rex, expnt, coeffs, ldc, stdtag, oshell_is_rel)
cc AJL/Begin: This is a wrapper call to the updated version of this function
cc SPIN-POLARISED ECPS
      implicit none
c#include "errquit.fh"
c#include "mafdecls.fh"
c#include "basdeclsP.fh"
c#include "nwc_const.fh"
c#include "basP.fh"
c#include "bas_exndcf_dec.fh"
c#include "ecpso_decP.fh"
c#include "stdio.fh"
      integer basis ! [input] basis handle
      character*(*) tag ! [input] tag on which to add contraction
      character*(*) stdtag ! [input] standard basis set name associated
*                                    with the tag     
      integer l_value ! [input] type of contraction (s, p, ..., sp)
      integer ngen ! [input] no. of contractions
      integer nprim ! [input] no. of primitives
      integer ldc ! [input] leading dimension of coeffs
      double precision expnt(nprim) ! [input] exponents
      double precision coeffs(ldc, 1:*) ! [input] coefficients
      double precision rex(nprim) ! [input] gaussian R exponents 
*. . . . . . . . . . . . . . . . .          (offset by 2)
      logical oshell_is_rel ! [input] flag for relativistic shell

      logical  bas_add_ucnt0
      external bas_add_ucnt0

      bas_add_ucnt =  bas_add_ucnt0(basis, tag, l_value,ngen, nprim, 
     &                rex, expnt, coeffs, ldc, stdtag, oshell_is_rel, 0)

      end
c
      logical function bas_add_ucnt0(basis, tag, l_value, ngen, nprim, 
     $     rex, expnt, coeffs, ldc, stdtag, oshell_is_rel, channel)
cc AJL/End
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "bas_exndcf_dec.fh"
#include "ecpso_decP.fh"
#include "stdio.fh"
      integer basis ! [input] basis handle
      character*(*) tag ! [input] tag on which to add contraction
      character*(*) stdtag ! [input] standard basis set name associated
*                                    with the tag     
      integer l_value ! [input] type of contraction (s, p, ..., sp)
      integer ngen ! [input] no. of contractions
      integer nprim ! [input] no. of primitives
      integer ldc ! [input] leading dimension of coeffs
      double precision expnt(nprim) ! [input] exponents
      double precision coeffs(ldc, 1:*) ! [input] coefficients
      double precision rex(nprim) ! [input] gaussian R exponents 
*. . . . . . . . . . . . . . . . .          (offset by 2)
      logical oshell_is_rel ! [input] flag for relativistic shell
c
      integer size_add ! amount to be added to exndcf array
      integer ind ! Index into basis function structures
      integer free ! Free space pointer
c
      integer i, itag, jtag, iu_cont, ntags ! Locals
      integer s_old, s_new, k_old, k_new, h_old, h_new
c
      logical oIs_ecp ! flag saying if this is ecp
      logical oIs_so ! flag saying if this is so potential
      logical bas_add_utag
      external bas_add_utag
c
      logical bas_check_handle
      external bas_check_handle
cc AJL/Begin/SPIN-POLARISED ECPS
      integer channel           ! [input] ecp channel [both/alpha/beta]
cc AJL/End
c
#include "bas_exndcf_sfn.fh"
#include "ecpso_sfnP.fh"
c
      oIs_ecp = Is_ECP_in(basis)
      oIs_so  = Is_SO_in(basis)
c
*Ul with rexp=0 :      if (oIs_ecp.and.l_value.eq.(-1)) then
*Ul with rexp=0 :        do i = 1,nprim
*Ul with rexp=0 :          if (int(rex(i)).eq.0) then
*Ul with rexp=0 :            write(luout,*)'This version of nwchem does not support ',
*Ul with rexp=0 :     &          'local potentials with an r-exponent of 0'
*Ul with rexp=0 :            call util_flush(luout)
*Ul with rexp=0 :            write(luout,*)'This should be fixed in a future release',
*Ul with rexp=0 :     &          ' of NWChem'
*Ul with rexp=0 :            call util_flush(luout)
*Ul with rexp=0 :            call errquit(' {ecp}bas_input fatal error ',911)
*Ul with rexp=0 :          endif
*Ul with rexp=0 :        enddo
*Ul with rexp=0 :      endif
c
c     adds a new general contraction on the specified tag.  If the
c     tag is not present it will also add that by calling bas_add_utag
c      
CC AJL/Begin/SPIN POLARISED ECPS
c      bas_add_ucnt = bas_check_handle(basis,'bas_add_ucnt')
c      if (.not. bas_add_ucnt) return
      bas_add_ucnt0 = bas_check_handle(basis,'bas_add_ucnt0')
      if (.not. bas_add_ucnt0) return
      ind = basis + BASIS_HANDLE_OFFSET
c
c     Make sure that the tag is in the list
c
c      bas_add_ucnt = bas_add_utag(basis, tag, stdtag, itag)
c      if (.not. bas_add_ucnt) return
      bas_add_ucnt0 = bas_add_utag(basis, tag, stdtag, itag)
      if (.not. bas_add_ucnt0) return
CC AJL/End
c
c     Update header information about all unique contractions on all 
c     tags.  Free points to next free word in the exndcf
c
*old: free = infbs_head(HEAD_NPRIM,ind)+infbs_head(HEAD_NCOEF,ind)+1
      free = infbs_head(HEAD_EXCFPTR,ind) + 1
      s_old = exndcf(SZ_exndcf,ind)
      size_add = nprim*ngen + nprim
      if (oIs_ecp.or.oIs_so) size_add = size_add + nprim
      if ((free+size_add-1) .gt. s_old) then
        h_old = exndcf(H_exndcf,ind)
        k_old = exndcf(K_exndcf,ind)
        s_new = free+size_add-1
        if (.not.ma_alloc_get(
     &      mt_dbl,s_new,' input for basis heap ',
     &      h_new, k_new)) then
          write(LuOut,*)'bas_add_ucnt0: too many prims/coeffs'
          write(LuOut,*)' allocated size for input is :',
     &        exndcf(SZ_exndcf,ind)
          write(LuOut,*)' size requested here         :',
     &        (free+size_add-1)
          bas_add_ucnt0 = .false.
          return
        endif
        call dfill(s_new,0.0d00,dbl_mb(k_new),1)
        exndcf(H_exndcf,ind) = h_new
        exndcf(K_exndcf,ind) = k_new
        exndcf(SZ_exndcf,ind) = s_new
        call dcopy(s_old,dbl_mb(k_old),1,dbl_mb(k_new),1)
        if (.not.ma_free_heap(h_old)) call errquit
     &      ('bas_add_ucnt: error freeing old exponents',911,
     &       BASIS_ERR)
      endif
*      if (infbs_head(HEAD_NCONT,ind)+1 .gt. nucont_bsmx) then
      if (infbs_head(HEAD_NCONT,ind) .gt. nucont_bsmx) then
         write(LuOut,*) 'bas_add_ucnt0: too many contractions '
         bas_add_ucnt0 = .false.
         return
      endif
c         
      infbs_head(HEAD_NCONT,ind) = infbs_head(HEAD_NCONT,ind) + 1
      infbs_head(HEAD_NPRIM,ind) = infbs_head(HEAD_NPRIM,ind) + nprim
      infbs_head(HEAD_NCOEF,ind) = infbs_head(HEAD_NCOEF,ind) +
     $     ngen*nprim
      infbs_head(HEAD_EXCFPTR,ind) =  infbs_head(HEAD_EXCFPTR,ind) +
     &      size_add
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
      infbs_tags(Tag_High_Ang,itag,ind) =
     &      max(infbs_tags(Tag_High_Ang,itag,ind),abs(l_value))
      infbs_tags(TAG_NCONT,itag,ind) = infbs_tags(TAG_NCONT,itag,ind)
     $     + 1
      infbs_tags(TAG_NPRIM,itag,ind) = infbs_tags(TAG_NPRIM,itag,ind)
     $     + nprim
      infbs_tags(TAG_NCOEF,itag,ind) = infbs_tags(TAG_NCOEF,itag,ind)
     $     + nprim*ngen
      if (infbs_tags(TAG_FCONT,itag,ind).eq.0) then
         if (itag .ne. ntags) call errquit
     $        ('bas_add_ucnt0: tag error', itag, BASIS_ERR)
         infbs_tags(TAG_FCONT,itag,ind) = infbs_head(HEAD_NCONT,ind)
         infbs_tags(TAG_LCONT,itag,ind) = infbs_head(HEAD_NCONT,ind)
      else
         infbs_tags(TAG_LCONT,itag,ind) =
     &      infbs_tags(TAG_LCONT,itag,ind) + 1
      endif
c
*. . . . . . . . . . . . . . . . . . . . . ! Index of new contraction
      iu_cont = infbs_tags(TAG_LCONT,itag,ind) 
c
      infbs_cont(CONT_TYPE, iu_cont,ind) = l_value
      infbs_cont(CONT_NPRIM,iu_cont,ind) = nprim
      infbs_cont(CONT_NGEN, iu_cont,ind) = ngen
      infbs_cont(CONT_IEXP, iu_cont,ind) = free
      infbs_cont(CONT_ICFP, iu_cont,ind) = free + nprim
cc AJL/Begin/SPIN-POLARISED ECPs
      infbs_cont(CONT_CHANNEL, iu_cont, ind) = channel
cc AJL/End
      if (oIs_ecp.or.oIs_so) then
        infbs_cont(Cont_Irexp, iu_cont, ind) =
     &      free + nprim + nprim*ngen
      else
*. . . . . . . . . . . . . . . . . . ! point to exponents for safety?
        infbs_cont(Cont_Irexp, iu_cont,ind) = free 
      endif
      if (oshell_is_rel) then
        infbs_cont(CONT_RELLS, iu_cont,ind) = 1
      else
        infbs_cont(CONT_RELLS, iu_cont,ind) = 0
      end if
c
c     Copy real data over
c
      call dcopy(nprim, expnt, 1, dbl_mb(mb_exndcf(free,ind)), 1)
      free = free + nprim
      do i = 1, ngen
         call dcopy
     &      (nprim, coeffs(1,i), 1, dbl_mb(mb_exndcf(free,ind)), 1)
         free = free + nprim
      enddo
      if (oIs_ecp.or.oIs_so)
     &    call dcopy(nprim, rex, 1, dbl_mb(mb_exndcf(free,ind)), 1)
c
c     Done
c
      end
*.....................................................................
      logical function bas_num_uce(basisin,nucent)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin, nucent 
c::local
      integer basis
c
#include "bas_exndcf_sfn.fh"
c
      bas_num_uce = bas_check_handle(basisin,'bas_getu_coeff')
      if (.not.bas_num_uce) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      nucent = infbs_head(HEAD_NTAGS,basis)
c
      bas_num_uce = .true.
c
      return
      end
*.....................................................................
      logical function bas_uce2cnr(basisin,ucenter,ifirst,ilast)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c::passed
      integer basisin, ucenter, ifirst, ilast
c::local
      integer basis
c
#include "bas_exndcf_sfn.fh"
c
      bas_uce2cnr = bas_check_handle(basisin,'bas_getu_coeff')
      if (.not.bas_uce2cnr) return

      basis = basisin + BASIS_HANDLE_OFFSET 
c
      ifirst = infbs_tags(TAG_FCONT,ucenter,basis)
      ilast = infbs_tags(TAG_LCONT,ucenter,basis)
c
      bas_uce2cnr = .true.
c
      return
      end
*.....................................................................
      logical function bas_uce_tag(basisin,ucent,tagout)
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
      integer basisin
      integer ucent
      character*(*) tagout
c::-local
      integer basis
      integer len_tagout
#include "bas_ibs_sfn.fh"
c
      bas_uce_tag = bas_check_handle(basisin,'bas_cont_tag')
      if (.not.bas_uce_tag) return

      basis = basisin + Basis_Handle_Offset

      len_tagout = len(tagout)
      tagout = bs_tags(ucent,basis)(1:len_tagout)
      bas_uce_tag = .true.
      end
*.....................................................................
      logical function bas_getu_coeff(basisin,icont,coeff)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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
#include "bas_exndcf_sfn.fh"
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
        write(LuOut,*)' bas_getu_coeff: ERROR '
        write(LuOut,*)' contraction range for basis is 1:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      mycoeffptr = infbs_cont(CONT_ICFP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      mygen   = infbs_cont(CONT_NGEN,myucont,basis)
      call dcopy ((myprim*mygen),
     &    dbl_mb(mb_exndcf(mycoeffptr,basis)),1,coeff,1)
c
      bas_getu_coeff = .true.
c
      return
      end
*.....................................................................
      logical function bas_getu_exponent(basisin,icont,exp)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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
#include "bas_exndcf_sfn.fh"
c
      bas_getu_exponent =
     &    bas_check_handle(basisin,'bas_getu_exponent')

      if (.not.bas_getu_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      icontmax = infbs_head(HEAD_NCONT,basis)
      myucont = icont

      bas_getu_exponent = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_getu_exponent)) then
        write(LuOut,*)' bas_getu_exponent: ERROR '
        write(LuOut,*)' contraction range for basis is 1:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      call dcopy(myprim,dbl_mb(mb_exndcf(myexptr,basis)),1,exp,1)
c
      bas_getu_exponent = .true.
c
      return
      end
*.....................................................................
      logical function bas_setu_coeff(basisin,icont,coeff,ncoeff)
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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
#include "bas_exndcf_sfn.fh"
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
        write(LuOut,*)' bas_setu_coeff: ERROR '
        write(LuOut,*)' contraction range for basis is 1:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      mycoeffptr = infbs_cont(CONT_ICFP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      mygen   = infbs_cont(CONT_NGEN,myucont,basis)
c
      bas_setu_coeff = ncoeff .eq. (myprim*mygen)
      if(.not.bas_setu_coeff) then
        write(LuOut,*)' bas_setu_coeff: ERROR '
        write(LuOut,*)' input and stored number of coefficients ',
     &         '(nprim*ngen) differ '
        write(LuOut,*)' input  nprim*ngen: ',ncoeff
        write(LuOut,*)' stored nprim*ngen: ',(myprim*mygen)
        return
      endif
      call dcopy(ncoeff,coeff,1,
     &    dbl_mb(mb_exndcf(mycoeffptr,basis)),1)
c
      bas_setu_coeff = .true.
c
      return
      end
*.....................................................................
      logical function bas_setu_exponent(basisin,icont,exp,nexp)
      implicit none
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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
#include "bas_exndcf_sfn.fh"
c
      bas_setu_exponent =
     &    bas_check_handle(basisin,'bas_setu_exponent')

      if (.not.bas_setu_exponent) return

      basis = basisin + BASIS_HANDLE_OFFSET 

      icontmax = infbs_head(HEAD_NCONT,basis)
      myucont = icont

      bas_setu_exponent = icont.gt.0.and.icont.le.icontmax
      if (.not.(bas_setu_exponent)) then
        write(LuOut,*)' bas_setu_exponent: ERROR '
        write(LuOut,*)' contraction range for basis is 1:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      bas_setu_exponent = myprim.eq.nexp
      if (.not.bas_setu_exponent) then
        write(LuOut,*)' bas_setu_exponent: ERROR '
        write(LuOut,*)' input and stored number of exponents ',
     &         '(nprim) differ '
        write(LuOut,*)' input  nprim: ',nexp
        write(LuOut,*)' stored nprim: ',myprim
        return
      endif
c
      call dcopy(nexp,exp,1,dbl_mb(mb_exndcf(myexptr,basis)),1)
c
      bas_setu_exponent = .true.
*.....................................................................
      end
