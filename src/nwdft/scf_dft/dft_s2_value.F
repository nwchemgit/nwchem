      subroutine dft_s2_value(geom,basis,symmetry,na,nb,
     &     nbf,ga_densa,ga_densb)
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
      integer geom,basis
      logical symmetry
      integer nbf,ga_densa,ga_densb,na,nb
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer iga_s,iga_tmp,iga_tmp2
      double precision ana,anb
      double precision s2,s2_exact
      double precision one,zero,two
      parameter(one=1.d0,zero=0.d0,two=2.d0)
c     
c     prod S*Pb*S*Pa
c     
c      call gacrea(geom,basis,nbf,nbf,'ovl ',iga_s,'atom')
      iga_s = ga_create_atom_blocked(geom, basis, 'S2 temp 1')

c      call gacrea(geom,basis,nbf,nbf,'ovl ',iga_tmp,'atom')
      iga_tmp = ga_create_atom_blocked(geom, basis, 'S2 temp 2')

c      call gacrea(geom,basis,nbf,nbf,'ovl ',iga_tmp2,'atom')
      iga_tmp2 = ga_create_atom_blocked(geom, basis, 'S2 temp 3')

      call ga_zero(iga_s)
      call int_1e_ga(basis,basis,iga_s,'overlap',symmetry)

      call ga_dgemm('N','N',nbf,nbf,nbf,one,iga_s,ga_densb,zero,iga_tmp)

      call ga_dgemm('N','N',nbf,nbf,nbf,one,iga_tmp,iga_s,zero,iga_tmp2)

      ana=dble(na)
      anb=dble(nb)
      s2_exact= 0.5d0*(ana-anb)*(.5d0*(ana-anb) + one)
      s2= 0.25d0*(ana-anb)**two+0.5d0*(ana+anb) - 
     &     ga_ddot(iga_tmp2,ga_densa)

c      call gadest(iga_s)
      if (.not. ga_destroy(iga_s)) call errquit
     &   ('dft_s2_value: could not destroy iga_s', 0, GA_ERR)
c      call gadest(iga_tmp)
      if (.not. ga_destroy(iga_tmp)) call errquit
     &   ('dft_s2_value: could not destroy iga_tmp', 0, GA_ERR)
c      call gadest(iga_tmp2)
      if (.not. ga_destroy(iga_tmp2)) call errquit
     &   ('dft_s2_value: could not destroy iga_tmp2', 0, GA_ERR)
c
      if(ga_nodeid().eq.0) then
        call dft_header('Expectation value of S2:  ')
        write (LuOut,12) s2,s2_exact
   12   format('      <S2> = ',f11.4,' (Exact =',f11.4,')')
        write(LuOut,*)
      endif

      return
      end
