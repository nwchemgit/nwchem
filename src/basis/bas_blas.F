*
* routines from basis.F that use blas calls
*
* $Id: bas_blas.F,v 1.4 2004-03-15 20:28:16 edo Exp $
*
      logical function bas_get_exponent(basisin,icont,expon)
c
c  returns the exponents of a given mapped contraction 
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin         ! [input] basis set handle
      integer icont           ! [input] mapped contraction index
      double precision expon(*) ! [output] exponents
c::local
      integer basis, myucont, icontmax
      integer myprim,myexptr
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      bas_get_exponent = bas_check_handle(basisin,'bas_get_exponent')
      if (.not.bas_get_exponent) return

      basis = basisin + Basis_Handle_Offset 
      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)
      bas_get_exponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_get_exponent)) then
        write(LuOut,*)' bas_get_exponent: ERROR '
        write(LuOut,*)' contraction range for basis is 0:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      call dcopy(myprim,dbl_mb(mb_exndcf(myexptr,basis)),1,expon,1)
c
      bas_get_exponent = .true.
c
      return
      end
*.....................................................................
      logical function bas_get_coeff(basisin,icont,coeff)
c
c returns the coefficients for the given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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

      basis = basisin + Basis_Handle_Offset 
c
      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)
      bas_get_coeff = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_get_coeff)) then
        write(LuOut,*)' bas_get_coeff: ERROR '
        write(LuOut,*)' contraction range for basis is 0:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
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
      logical function bas_set_exponent(basisin,icont,expon,nexp)
c
c sets the exponents to the input array for a given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
c::function
      logical bas_check_handle
      external bas_check_handle
c:blas
c     dcopy
c::passed
      integer basisin ! [input] basis set handle        
      integer icont ! [input] mapped contraction index
      integer nexp ! [input] number of exponents for contraction
      double precision expon(nexp) ! [input] "new" exponents for 
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

      basis = basisin + Basis_Handle_Offset 

      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)

      bas_set_exponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_set_exponent)) then
        write(LuOut,*)' bas_set_exponent: ERROR '
        write(LuOut,*)' contraction range for basis is 0:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      bas_set_exponent = myprim.eq.nexp
      if (.not.bas_set_exponent) then
        write(LuOut,*)' bas_set_exponent: ERROR '
        write(LuOut,*)' input and stored number of exponents ',
     &         '(nprim) differ '
        write(LuOut,*)' input  nprim: ',nexp
        write(LuOut,*)' stored nprim: ',myprim
        return
      endif
      call dcopy(nexp,expon,1,dbl_mb(mb_exndcf(myexptr,basis)),1)
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
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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

      basis = basisin + Basis_Handle_Offset 
c
      icontmax = ncont_tot_gb(basis)
      myucont = sf_ibs_cn2ucn(icont,basis)

      bas_set_coeff = icont.ge.0.and.icont.le.icontmax
      if (.not.(bas_set_coeff)) then
        write(LuOut,*)' bas_set_coeff: ERROR '
        write(LuOut,*)' contraction range for basis is 0:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      mycoeffptr = infbs_cont(CONT_ICFP,myucont,basis)
      myprim  = infbs_cont(CONT_NPRIM,myucont,basis)
      mygen   = infbs_cont(CONT_NGEN,myucont,basis)
c
      bas_set_coeff = ncoeff .eq. (myprim*mygen)
      if(.not.bas_set_coeff) then
        write(LuOut,*)' bas_set_coeff: ERROR '
        write(LuOut,*)' input and stored number of coefficients ',
     &         '(nprim*ngen) differ '
        write(LuOut,*)' input  nprim*ngen: ',ncoeff
        write(LuOut,*)' stored nprim*ngen: ',(myprim*mygen)
        return
      endif
      call dcopy(ncoeff,coeff,1,dbl_mb(mb_exndcf(mycoeffptr,basis)),1)
c
      bas_set_coeff = .true.
c
      return
      end
*.....................................................................
      logical function ecp_get_rexponent(ecpin,icont,expon)
c
c returns the r-exponents of a given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::function
      logical ecp_check_handle
      external ecp_check_handle
c:blas
c     dcopy
c::passed
      integer ecpin         ! [input] ecp set handle
      integer icont           ! [input] mapped contraction index
      double precision expon(*) ! [output] exponents
c::local
      integer ecp, myucont, icontmax
      integer myprim,myexptr
c:: statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      ecp_get_rexponent = ecp_check_handle(ecpin,'ecp_get_rexponent')
      if (.not.ecp_get_rexponent) return

      ecp = ecpin + Basis_Handle_Offset 
      icontmax = ncont_tot_gb(ecp)
      myucont = sf_ibs_cn2ucn(icont,ecp)
      ecp_get_rexponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(ecp_get_rexponent)) then
        write(LuOut,*)' ecp_get_rexponent: ERROR '
        write(LuOut,*)' contraction range for ecp is 0:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IREXP,myucont,ecp)
      myprim  = infbs_cont(CONT_NPRIM,myucont,ecp)
      call dcopy(myprim,dbl_mb(mb_exndcf(myexptr,ecp)),1,expon,1)
c
      ecp_get_rexponent = .true.
c
      return
      end
*.....................................................................
      logical function ecp_set_rexponent(ecpin,icont,expon,nexp)
c
c sets the exponents to the input array for a given mapped contraction
c
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "basdeclsP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
c::function
      logical ecp_check_handle
      external ecp_check_handle
c:blas
c     dcopy
c::passed
      integer ecpin ! [input] ecp set handle        
      integer icont ! [input] mapped contraction index
      integer nexp ! [input] number of exponents for contraction
      double precision expon(nexp) ! [input] "new" r-exponents 
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

      ecp = ecpin + Basis_Handle_Offset 

      icontmax = ncont_tot_gb(ecp)
      myucont = sf_ibs_cn2ucn(icont,ecp)

      ecp_set_rexponent = icont.ge.0.and.icont.le.icontmax
      if (.not.(ecp_set_rexponent)) then
        write(LuOut,*)' ecp_set_rexponent: ERROR '
        write(LuOut,*)' contraction range for ecp is 0:',
     &         icontmax
        write(LuOut,*)' information requested for contraction:',icont
        return
      endif
c
      myexptr = infbs_cont(CONT_IEXP,myucont,ecp)
      myprim  = infbs_cont(CONT_NPRIM,myucont,ecp)
      ecp_set_rexponent = myprim.eq.nexp
      if (.not.ecp_set_rexponent) then
        write(LuOut,*)' ecp_set_rexponent: ERROR '
        write(LuOut,*)' input and stored number of exponents ',
     &         '(nprim) differ '
        write(LuOut,*)' input  nprim: ',nexp
        write(LuOut,*)' stored nprim: ',myprim
        return
      endif
      call dcopy(nexp,expon,1,dbl_mb(mb_exndcf(myexptr,ecp)),1)
c
      ecp_set_rexponent = .true.
c
      return
*.....................................................................
      end
