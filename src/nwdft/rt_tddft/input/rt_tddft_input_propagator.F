
C
C     rt_tddft_input_propagator.F
C
C
C     XXX NOT USED ANYMORE

      subroutine rt_tddft_input_propagator (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: rtdb


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_input_propagator: "

      
C     == Variables ==
      character*255 method
      integer prop_method
      
      set_method = .false.
      set_predictor = .false.
      prop_method = -99
      
      if (.not. inp_a (method)) call errquit (pname//
     $     "incorrect syntax--expecting: propagator <method>")

      if (method.eq."euler") then
         prop_method = 1        ! euler
         if (set_method) call errquit (pname//
     $        "cannot set multiple propagators",0,0)

      elseif (method.eq."rk4") then
         prop_method = 2        ! 4th order runge-kutta
         if (set_method) call errquit (pname//
     $        "cannot set multiple propagators",0,0)
         
      elseif (method.eq."magnus") then
         prop_method = 3        ! 2nd order magnus
         if (set_method) call errquit (pname//
     $        "cannot set multiple propagators",0,0)

      else
         call errquit (pname//"invalid propagator: "//method)
      endif

      if (.not.rtdb_put(rtdb,'rt_tddft:prop_method',
     $     mt_int,1,prop_method))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      call errquit ("XXX", 0, 0)
      
      end subroutine 




      subroutine rt_tddft_input_propagator_OLD (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: rtdb


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_input_propagator: "

      
C     == Variables ==
      character*255 method, opt, subopt
      integer prop_method, exp_method
      logical got_opt

      logical set_method, set_predictor, set_exp

      
      set_method = .false.
      set_predictor = .false.
      set_exp = .false.

      prop_method = -99
      exp_method = -99

      if (.not. inp_a (method)) call errquit (pname//
     $     "incorrect syntax--expecting: propagator <method> [opts]")
      
      if (method.eq."rk4") then
         prop_method = 1        ! rk4
         if (set_method) call errquit (pname//
     $        "cannot set multiple propagators",0,0)
         set_method = .true.
         
      elseif (method.eq."magnus") then
         got_opt = .true.
         if (set_method) call errquit (pname//
     $        "cannot set multiple propagators",0,0)
         set_method = .true.

C     (set to default magnus options, will overwrite if keywords found)
         prop_method = 4
         exp_method = 1
         
         do while (got_opt)
            got_opt = inp_a (opt)
            if ( trim(opt) .eq. "") exit
            
            if (opt.eq."predictor") then
               if (.not. inp_a (subopt)) call errquit (pname//
     $              "expecting option after 'predictor'")
               if (subopt.eq."pc2") then
                  prop_method = 3 !magnus2-pc
                  if (set_predictor) call errquit (pname//
     $                 "cannot set multiple magnus predictors",0,0)
                  set_predictor = .true.
               elseif (subopt.eq."interpol2") then
                  prop_method = 4 !magnus2-interpolation
                  if (set_predictor) call errquit (pname//
     $                 "cannot set multiple magnus predictors",0,0)
                  set_predictor = .true.
               elseif (subopt.eq."lin2") then
                  prop_method = 2 !magnus2-linear
                  if (set_predictor) call errquit (pname//
     $                 "cannot set multiple magnus predictors",0,0)
                  set_predictor = .true.
               else
                  call errquit (pname//"unknown magnus predictor: "//
     $                 trim(subopt), 0, 0)
               endif            !end magnus predictor
               
               
            elseif (opt.eq."exp") then
               if (.not. inp_a (subopt)) call errquit (pname//
     $              "expecting option after 'exp'")
               if (subopt.eq."pseries") then
                  exp_method = 1
                  if (set_exp) call errquit (pname//
     $                 "cannot set multiple magnus exp",0,0)
                  set_exp = .true.
               elseif (subopt.eq."diag") then
                  exp_method = 2
                  call errquit (pname//"diag buggy? disabled for now",
     $                 0,0)
                  if (set_exp) call errquit (pname//
     $                 "cannot set multiple magnus exp",0,0)
                  set_exp = .true.
               elseif (subopt.eq."bch") then
                  exp_method = 3
                  if (set_exp) call errquit (pname//
     $                 "cannot set multiple magnus exp",0,0)
                  set_exp = .true.
               else
                  call errquit (pname//"unknown magnus exp: "//
     $                 trim(subopt), 0, 0)
               endif            !end magnus exp 
               
            else
               call errquit (pname//"unknown magnus option: "//
     $              trim(opt),0,0)
            endif               !end magnus options
         enddo                  !end loop over all propagator options
      endif                     !end prop method
      
      if (.not.rtdb_put(rtdb,'rt_tddft:prop_method',
     $     mt_int,1,prop_method))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      
      if (.not.rtdb_put(rtdb,'rt_tddft:exp_method',
     $     mt_int,1,exp_method))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            
      end subroutine
c $Id$
