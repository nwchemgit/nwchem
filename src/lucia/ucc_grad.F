*----------------------------------------------------------------------*
      subroutine ucc_grad(ccvec1,ccvec2,civec1,civec2,
     &                    n_cc_amp,
     &                    lu_amp,lu_exrf,luhexrf)
*----------------------------------------------------------------------*
*     calculate the gradient of the Unitary CC Energy functional using
*     the Wilson formula and Gauss-Legendre quadrature
*
*     lu_amp   contains the amplitudes
*     lu_exrf  contains exp(sig)|R>
*     luhexrf  contains H exp(sig)|R>
*
*     result is returned on ccvec1
*
*     andreas, adapted from gtbce_gradE, August 2004 
*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00, npnts = 5, mx_term = -150, iopsym = -1
      real(8), parameter ::
     &     thresh = 1d-20

      integer, intent(in) ::
     &     n_cc_amp,
     &     lu_amp, lu_exrf, luhexrf
      real(8), intent(inout) ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp),
     &     civec1(*), civec2(*)

      character ::
     &     cctype*8
      integer ::
     &     lblk, mxpnts, ipnt,
     &     lusc1,lusc2,lusc3,lusc4,lusc5,lusc6,lusc7
      real(8) ::
     &     cpu0, cpu, wall0, wall,
     &     xnrm, etest, etest2, xnrm2, dltalp
      
      real(8) ::
     &     alp(npnts), wght(npnts)

      integer, external ::
     &     iopen_nus

      real(8), external ::
     &     inprod, inprdd

      call atim(cpu0,wall0)

      lblk = -1
      if (ntest.ge.5) then
        write (6,*) '====================='
        write (6,*) ' This is ucc_grad'
        write (6,*) '====================='
        write (6,*) 
        write (6,*) 'on entry: '
        write (6,*) ' lu_amp:  ', lu_amp 
        write (6,*) ' lu_exrf: ', lu_exrf
        write (6,*) ' luhexrf: ', luhexrf
      end if

      if (ntest.ge.1000) then
        write(6,*) 'on entry:'
        write(6,*) 'e^G|0> on LU_EXRF'
        call wrtvcd(civec1,lu_exrf,1,lblk)
        write(6,*) 'H e^G|0> on LUHEXRF'
        call wrtvcd(civec1,luhexrf,1,lblk)
      end if
      
      ! for I/O
      lblk = -1
      ! for expt_ref
      cctype='GEN_CC'

      lusc1 = iopen_nus('UCCGRDSCR1')
      lusc2 = iopen_nus('UCCGRDSCR2')
      lusc3 = iopen_nus('UCCGRDSCR3')
      lusc4 = iopen_nus('UCCGRDSCR4')
      lusc5 = iopen_nus('UCCGRDSCR5')
      lusc6 = iopen_nus('UCCGRDSCR6')
      lusc7 = iopen_nus('UCCGRDSCR7')

*--------------------------------------------------------------------*
* set up points and weights
*--------------------------------------------------------------------*
      call gl_weights(0d0,1d0,npnts,alp,wght)

      mxpnts=npnts
      ! if G == 0 ...
      xnrm2 = inprod(ccvec1,ccvec1,n_cc_amp)
      ! ... things are trivial and we evaluate the formula only once
      if (xnrm2.lt.thresh) then
        mxpnts=1
        wght(1)=1d0
        alp(1)=0d0
        if (ntest.ge.5) then
          write(6,*) 'Detected zero amplitudes: ',
     &               'only case alpha = 0 will be processed'
        end if
      end if

*--------------------------------------------------------------------*
* H|0tilde>
*  init lusc1
*--------------------------------------------------------------------*
      call copvcd(luhexrf,lusc1,civec1,1,lblk)

**-------------------------------------------------------------------*
* loop over quadrature points
**-------------------------------------------------------------------*
      do ipnt = 1, mxpnts
        if (ntest.ge.5) then
          write(6,*) 'info for quadrature point: ', ipnt,'/',npnts
          write(6,*) 'point, weight: ', alp(ipnt), wght(ipnt)
        end if

        if (ipnt.gt.1.and.(alp(ipnt).le.alp(ipnt-1))) then
          write(6,*) 'quadrature point should be in ascending order!'
          stop 'ucc_grad > quadrature '
        end if

        ! read amplitudes
        call vec_from_disc(ccvec1,n_cc_amp,1,lblk,lu_amp)

        if (ipnt.eq.1) then
          dltalp = alp(1)
        else
          dltalp = alp(ipnt)-alp(ipnt-1)
          call copvcd(lusc2,lusc1,civec1,1,lblk)
        end if
*--------------------------------------------------------------------*
* |a_i> = exp(a_i G^+) [H exp(G)|0>]
*       = exp((a_i-a_{i-1})G^+)[exp(a_{i-1}G^+) H exp(G)|0>]
*       = exp(-(a_i-a_{i-1})G) [exp(-a_{i-1}G) H exp(G)|0>]
*  result on lusc2
*--------------------------------------------------------------------*
        if (ntest.ge.5) then
          write(6,*)
     &         'constructing |a_i> = exp(a_i G^+) [H exp(G)|0>]'
        end if
        
        if (abs(dltalp).lt.1d-20) then
          call copvcd(lusc1,lusc2,civec1,1,lblk)
        else
          ! scale it (minus sign: G^+ = -G)
          call scalve(ccvec1,-dltalp,n_cc_amp)
          call expt_ref2(lusc1,lusc2,lusc4,lusc5,lusc6,
     &         thresh,mx_term, ccvec1, ccvec2, civec1, civec2,          
     &         n_cc_amp,cctype, iopsym)
          if (ntest.ge.5) then
            xnrm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))
            write(6,*) '|dlta G^+|, dlta = ',xnrm, dltalp,
     &                 'for alp(i) = ', alp(ipnt) 
          end if

        end if

*--------------------------------------------------------------------*
* |b_i> = exp(-a_i G)exp(G)|0> =
*       = exp(-(a_i-a_{i-1})G) [exp(-a_{i-1}G)exp(G)|0>]
*  result on lusc3
*--------------------------------------------------------------------*
        if (ipnt.eq.1) then
          call copvcd(lu_exrf,lusc1,civec1,1,lblk)
        else
          call copvcd(lusc3,lusc1,civec1,1,lblk)
        end if

        if (ntest.ge.5) then
          write(6,*) 'constructing |b_i> = exp(-a_i G) exp(G)|0>]'
        end if

        if (abs(dltalp).lt.1d-20) then
          call copvcd(lusc1,lusc3,civec1,1,lblk)          
        else
          ! scaled vector is already on ccvec1
          call expt_ref2(lusc1,lusc3,lusc4,lusc5,lusc6,
     &         thresh,mx_term, ccvec1, ccvec2, civec1, civec2,
     &         n_cc_amp,cctype, iopsym)

          if (ntest.ge.5) then
            xnrm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
            etest = inprdd(civec1,civec2,lusc3,lusc3,1,lblk)
            etest2= inprdd(civec1,civec2,lusc2,lusc3,1,lblk)
            write(6,*) '|dltaG|, dlta = ',xnrm, dltalp
            write(6,*) '<b_i|b_i> = ', etest,
     &           'for alp(i) = ', alp(ipnt) 
            write(6,*) '<a_i|b_i>  = ', etest2,
     &           'for alp(i) = ', alp(ipnt) 
          end if
        end if
        
*--------------------------------------------------------------------*
* dE_u +=   w_i (<a_i|gamma_u-gamma^+_u|b_i>)
*         = w_i (<a_i|gamma_u|b_i> - <b_i|gamma_u|a_i>)
*--------------------------------------------------------------------*
        
        if (ntest.ge.1000) then
          write(6,*) 'Before calling sigden_cc:'
          write(6,*) '|a_i> on lusc2:'
          call wrtvcd(civec1,lusc2,1,lblk)
          write(6,*) '|b_i> on lusc3:'
          call wrtvcd(civec1,lusc3,1,lblk)
        end if

        call den_gcc_s(civec1,civec2,lusc2,lusc3,ccvec1,ccvec2,-1)

        if (ipnt.gt.1) then
          call vec_from_disc(ccvec2,n_cc_amp,1,lblk,lusc7)
          call vecsum(ccvec1,ccvec1,ccvec2,2d0*wght(ipnt),1d0,n_cc_amp)
        else
          call scalve(ccvec1,2d0*wght(ipnt),n_cc_amp)
        end if
        if (ipnt.lt.mxpnts)
     &       call vec_to_disc(ccvec1,n_cc_amp,1,lblk,lusc7)

      end do ! loop over quadrature points


      if (ntest.ge.5) then
        xnrm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))
        write(6,*) ' ucc_grad > '
        write(6,*) '     n_cc_amp,norm of grad: ',n_cc_amp,xnrm
      end if
      if (ntest.ge.100) then
        call wrt_cc_vec2(ccvec1,6,'GEN_CC')
      end if

c prelim: divide by 2d0 again
      call scalve(ccvec1,0.5d0,n_cc_amp)
      
      call relunit(lusc1,'delete')
      call relunit(lusc2,'delete')
      call relunit(lusc3,'delete')
      call relunit(lusc4,'delete')
      call relunit(lusc5,'delete')
      call relunit(lusc6,'delete')
      call relunit(lusc7,'delete')

      call atim(cpu,wall)
      call prtim(6,'time in ucc_grad',cpu-cpu0,wall-wall0)

      return
      end
c $Id$
