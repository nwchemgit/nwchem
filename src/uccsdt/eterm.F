      subroutine uccsdt_eterm(qO_handles,qV_handles)
c
c$Id: eterm.F,v 1.1 2002-08-28 22:18:19 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer qO_handles(8,0:7), qV_handles(8,0:7)
C
C     GENERATION OF E-TERM
C

c (1) R(E1) aa = +(1-(ef)) t(me,na) [ f(fa) - 0.5*q_1_O - q_2_O] m=e=n=f=a=alpha
c (1) R(E1) ab =           t(ma,nf) [ f(ea) - 0.5*q_1_O - q_2_O] m=e=a=alpha n=f=beta
c (2) R(E1) bb = +(1-(ef)) t(me,na) [ f(fa) - 0.5*q_4_O - q_5_O] m=e=n=f=a=beta 
c (3) R(E2) aa = -(1-(mn)) t(me,if) [ f(in) + 0.5*q_1_V + q_2_V] m=e=n=f=i=alpha
c (3) R(E4) ab =           t(ie,nf) [ f(im) + 0.5*q_1_V + q_2_V] i=e=m=alpha n=f=beta
c (4) R(E2) bb = -(1-(mn)) t(me,if) [ f(in) + 0.5*q_4_V + q_5_V] m=e=n=f=i=beta      
c (4) R(E3) ab =           t(me,if) [ f(in) + 0.5*q_4_V + q_5_V] m=e=alpha i=n=f=beta
c (5) R(E2) ab =           t(me,na) [ f(af) - 0.5*q_4_O - q_8_O] m=e=alpha a=n=f=beta
c
c     first combine terms in [] and them "matmul" with t's
c     (1) is in qO(2,sym)
c     (2) is in qO(5,sym)
c     (3) is in qV(2,sym)
c     (4) is in qV(5,sym)
c     (5) is in qO(8,sym)
c
c     Who has the f blocks, they should be in format
c
      do symf = 0, nir-1
         syma = symf
         call ga_add(-1.0d0,qO(2,symf),-0.5d0,qO(1,symf),qO(2,symf))
         call ga_add( 1.0d0,qO(2,symf), 1.0d0,fO(a,a,symf),qO(2,symf))   # (1)
         call ga_add(-1.0d0,qO(5,symf),-0.5d0,qO(4,symf),qO(5,symf))
         call ga_add( 1.0d0,qO(5,symf), 1.0d0,fO(b,b,symf),qO(5,symf))   # (2)
         call ga_add( 1.0d0,qV(2,symf), 0.5d0,qV(1,symf),qV(2,symf))
         call ga_add( 1.0d0,qV(2,symf), 1.0d0,fV(a,a,symf),qV(2,symf))   # (3)
         call ga_add( 1.0d0,qV(5,symf), 0.5d0,qV(4,symf),qV(5,symf))
         call ga_add( 1.0d0,qV(5,symf), 1.0d0,fV(b,b,symf),qV(5,symf))   # (4)
         call ga_add(-1.0d0,qO(8,symf),-0.5d0,qO(4,symf),qO(8,symf))
         call ga_add( 1.0d0,qO(8,symf), 1.0d0,fO(a,a,symf),qO(8,symf))   # (5)
      end do
      do symnf = 0, nir-1
         symme = symnf
         symna = symnf
         call get_T(g_t2_mena_pure,symnf,1,1,1,1,x)   #note, no backtransform needed
         do symf = 0, nir-1
            syma = symf
            symn = ieor(symf,symnf)
            do loop over a
               matmul(g_t2_mena_pure,qO(2,symf)) into R(E1) aa   #(1)
            end do
            do loop over a
               matmul(g_t2_mena_pure,qV(2,symn)) into R(E2) aa   #(3)
            end do
         end do 
         if (.not. ga_destroy(g_t2_mena_pure)) call errquit()
         do symf = 0, nir-1
            syma = symf
            symn = ieor(symf,symnf)
            do loop over a
               matmul(g_t2_mena_pure,qO(5,symf)) into R(E1) bb   #(2)
            end do
            do loop over a
               matmul(g_t2_mena_pure,qV(5,symn)) into R(E2) bb   #(4)
            end do
         end do
         if (.not. ga_destroy(g_t2_mena_mix)) call errquit()
         call get_T(g_t2_mena_pure,symnf,0,0,1,1,x)   #note, no backtransform needed
         do symf = 0, nir-1
            syma = symf
            symn = ieor(symf,symnf)
            do loop over a
               matmul(g_t2_mena_mix,qO(2,syma)) into R(E1) ab   #(1)
            end do  
            do loop over a
               matmul(g_t2_mena_mix,qV(2,symn)) into R(E4) ab   #(3)
            end do
            do loop over a
               matmul(g_t2_mena_mix,qV(5,symn)) into R(E3) ab   #(4)
            end do
            do loop over a
               matmul(g_t2_mena_mix,qV(8,symn)) into R(E2) ab   #(5)
            end do
         end do
      end do

      do symf = 0, nir-1
         do inum = 1, 8
            if (qO_handle(inum,symf) .gt. -1) then
                if (.not. ga_destroy(qO_handle(inum,symf))
     &             call errquit()
            if (qV_handle(inum,symf) .gt. -1) then
                if (.not. ga_destroy(qV_handle(inum,symf))
     &             call errquit()
            end if
         end do
      end do
c
      end
