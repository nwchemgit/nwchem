      subroutine bas_create_uncont(usertdb,basis_in,name,basis_out)
c $Id$
      implicit none
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "sym.fh"
#include "apiP.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "bgj.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
c
c     Given a contracted basis, create a new uncontracted basis on rtdb
c
c     function declarations
c
      logical bas_add_ucnt_init
      logical bas_add_ucnt_tidy
      logical bas_add_ucnt
      logical bas_set_spherical
      logical bas_get_spherical
      logical bas_ucontinfo
      logical bas_num_uce
      logical bas_uce2cnr
      logical bas_getu_exponent
      logical bas_uce_tag
      external bas_add_ucnt_init
      external bas_add_ucnt_tidy
      external bas_add_ucnt
      external bas_set_spherical
      external bas_get_spherical
      external bas_ucontinfo
      external bas_num_uce
      external bas_uce2cnr
      external bas_getu_exponent
      external bas_uce_tag
c
c     arguments
c
      integer usertdb                ! Handle of rtdb
      integer basis_in               ! Handle of input contracted basis
      character*(*) name             ! Name of the basis set for rtdb
      integer basis_out              ! Handle new uncontracted basis
c
c     local variables
c
      integer Expm_a, ExpIm_a        ! MA array for storing exponents
      integer Lm_a, LIm_a            ! MA array for storing l-values of exponents
      integer Sm_a, SIm_a            ! MA array for storing shells of exponents
      integer Am_a, AIm_a            ! MA array for storing centers of exponents
      integer ipr_mx                 ! Maximum # of primitives in basis
      integer usegeom                ! Handle to geometry basis_in
      integer nuniqatom              ! # of unique atoms
      integer iatom, ishell, iprim
      integer ifirst, ilast, icont
      integer cartsph, l_type
      integer counter, nuniqprim
      character*16 usetag
      logical spheri
c
      integer Pmnm_a,PmnIm_a         ! Following indices are used for generating 
      integer Pmxm_a,PmxIm_a         ! a fitting basis set based on exponents and
      integer l_max,nprim            ! angular momenta
      double precision maxprim,minprim,alpha
c
c     Get geometry and rtdb handles for basis_in
c
      if (.not. bas_geom(basis_in,usegeom))
     $   call errquit('bas_create_uncont: bas_geom failed',0, BASIS_ERR)
c
c     Create a new fitting basis set
c
      if (.not. bas_create(basis_out,name))
     $   call errquit('bas_create_uncont: bas_create failed',0,
     &       BASIS_ERR)
c
c     Check if the user generated a fitting basis on input called
c     "D-K basis", if so load, normalize and return
c
      if (bas_rtdb_load(usertdb,usegeom,basis_out,'D-K basis')) then
        call int_normalize(usertdb,basis_out)
        return
      end if
c
c     Create a new uncontracted basis and initialize adding function
c
      if (.not. bas_add_ucnt_init(basis_out)) call errquit
     &      ('bas_create_uncont: bas_add_ucnt_init failed',0,
     &       BASIS_ERR)
c
c     Get # unique centers 
c
      if (.not. bas_num_uce(basis_in,nuniqatom)) 
     $   call errquit ('bas_create_uncont: bas_numcont failed',0,
     &       BASIS_ERR)
c
c     Get maximum number of exponents for unique atoms
c
      ipr_mx = 0
      do 7 iatom = 1, nuniqatom
         if (.not.bas_uce2cnr(basis_in,iatom,ifirst,ilast))
     $      call errquit ('bas_create_uncont: bas_uce2cnr failed',0,
     &       BASIS_ERR)
         do 8 ishell = ifirst, ilast
            if (.not. bas_ucontinfo(basis_in,ishell,l_type,
     $                             nprim,icont,cartsph))
     $         call errquit ('bas_create_uncont: bas_ucontinfo ',0,
     &       BASIS_ERR)
            ipr_mx = ipr_mx + nprim
    8    continue
    7 continue
c
c     Make arrays for exponents plus l-types and shells of exponents
c
      if (.not. ma_push_get(MT_DBL,ipr_mx,'primitives',Expm_a,ExpIm_a))
     $     call errquit('bas_create_uncont: ma_push_get Exp array',0,
     &       MA_ERR)
      if (.not. ma_push_get(MT_INT,ipr_mx,'ltypes',Lm_a,LIm_a))
     $     call errquit('bas_create_uncont: ma_push_get L array',0,
     &       MA_ERR)
      if (.not. ma_push_get(MT_INT,ipr_mx,'shells',Sm_a,SIm_a))
     $     call errquit('bas_create_uncont: ma_push_get Exp array',0,
     &       MA_ERR)
      if (.not. ma_push_get(MT_INT,ipr_mx,'atoms',Am_a,AIm_a))
     $     call errquit('bas_create_uncont: ma_push_get Exp array',0,
     &       MA_ERR)
c
c     Extract exponents and data from basis_in for all unique atoms 
c
      counter = 0
      do 1 iatom = 1, nuniqatom
         if (.not.bas_uce2cnr(basis_in,iatom,ifirst,ilast))
     $      call errquit ('bas_create_uncont: bas_uce2cnr failed',0,
     &       BASIS_ERR)
         do 2 ishell = ifirst, ilast
            if (.not. bas_ucontinfo(basis_in,ishell,l_type,
     $                             nprim,icont,cartsph))
     $         call errquit ('bas_create_uncont: bas_ucontinfo ',0,
     &       BASIS_ERR)
            if (.not. bas_getu_exponent(basis_in,ishell,
     $                                 dbl_mb(ExpIm_a+counter)))
     $         call errquit ('bas_create_uncont: bas_getu_exponent',0,
     &       BASIS_ERR)
            do 3 iprim = 1, nprim
                 int_mb(LIm_a+counter+iprim-1) = l_type
                 int_mb(SIm_a+counter+iprim-1) = ishell
                 int_mb(AIm_a+counter+iprim-1) = iatom
    3       continue
            counter = counter + nprim
    2    continue
    1 continue
c
c   Find unique exponents and store exponents and data into basis_out
c
      call bas_unique_sort(dbl_mb(ExpIm_a),int_mb(LIm_a),int_mb(SIm_a),
     $                 int_mb(AIm_a),ipr_mx,nuniqprim)
c
c   If the option to use a more extensive, flexible fitting basis set
c   setup data to handle even tempered basis set construction
c
      if (dk_even_bas) then
         if ((ga_nodeid() .eq. 0) .and. 
     $       util_print('douglas-kroll', print_low)) write(luout,30)
   30    format(/12x, 
     $   'Extended even tempered basis constructed for Douglas-Kroll'/)

c
c        Check if sp or spd combined basis is used, they cannot be handled
c        in the flexible basis set construction code below
c
         do 9 iprim = 1, nuniqprim
            if (int_mb(LIm_a+iprim-1) .lt.0) call errquit(
     $       'Cannot use sp or spd to create flexi fitting basis DK',0,
     &       BASIS_ERR)
    9    continue
c
c        Get highest angular momentum for array definitions and create arrays
c
         if (.not. bas_high_angular(basis_in,l_max))
     $      call errquit('bas_create_uncont: bas_high_angular failed',0,
     &       BASIS_ERR)
         ipr_mx = l_max + 1
         if (.not. ma_push_get(MT_DBL,ipr_mx,'maxprim',Pmxm_a,PmxIm_a))
     $      call errquit('bas_create_uncont: ma_push_get Pmx array',0,
     &       MA_ERR)
         if (.not. ma_push_get(MT_DBL,ipr_mx,'minprim',Pmnm_a,PmnIm_a))
     $      call errquit('bas_create_uncont: ma_push_get Pmn array',0,
     &       MA_ERR)
c
c        Loop over each unique center
c
         do 17 iatom = 1, nuniqatom
c
c        Get l-values and l_max and tag for this center
c
            do 16 iprim = 0, ipr_mx - 1
               dbl_mb(PmxIm_a+iprim) = 0.0d0
               dbl_mb(PmnIm_a+iprim) = 1.0d20
   16       continue
            l_max = 0
            do 18 iprim = 1, nuniqprim
               if (int_mb(AIm_a+iprim-1) .eq. iatom) then
                  if (l_max .lt. int_mb(LIm_a+iprim-1))
     $               l_max = int_mb(LIm_a+iprim-1)
               counter = int_mb(LIm_a+iprim-1)
               if (dbl_mb(PmxIm_a+counter) .lt. dbl_mb(ExpIm_a+iprim-1)) 
     $             dbl_mb(PmxIm_a+counter) = dbl_mb(ExpIm_a+iprim-1)
               if (dbl_mb(PmnIm_a+counter) .gt. dbl_mb(ExpIm_a+iprim-1)) 
     $             dbl_mb(PmnIm_a+counter) = dbl_mb(ExpIm_a+iprim-1)
               endif
   18       continue
c
c        Get max and min exponent for each l-value, based on ratio get 
c        number of exponents (compute adjusted ratio based on integer
c        number of exponents) using the following extra rules:
c        s space has to incorporate p space
c        p space has to incorporate s and d space
c        d space has to incorporate p and f space, etc
c
c        When max and min primitive are set, compute the number of primitive
c        functions, expansion ratio for number of primitive and store in basis
c
            do 20 l_type = 0, l_max + 1
               if (l_type .eq. 0) then
                  if (l_max .gt. 0) then
                     maxprim = max(dbl_mb(PmxIm_a),
     $                             dbl_mb(PmxIm_a+1))
                     minprim = min(dbl_mb(PmnIm_a),
     $                             dbl_mb(PmnIm_a+1))
                  else
                     maxprim = dbl_mb(PmxIm_a)
                     minprim = dbl_mb(PmnIm_a)
                  endif
               elseif (l_type .eq. (l_max + 1)) then
                  maxprim = dbl_mb(PmxIm_a+l_type-1)
                  minprim = dbl_mb(PmnIm_a+l_type-1)
               elseif (l_type .eq. l_max) then
                  maxprim = max(dbl_mb(PmxIm_a+l_type),
     $                          dbl_mb(PmxIm_a+l_type-1))
                  minprim = min(dbl_mb(PmnIm_a+l_type),
     $                          dbl_mb(PmnIm_a+l_type-1))
               else
                  maxprim = max(dbl_mb(PmxIm_a+l_type),
     $                          dbl_mb(PmxIm_a+l_type-1),
     $                          dbl_mb(PmxIm_a+l_type+1))
                  minprim = min(dbl_mb(PmnIm_a+l_type),
     $                          dbl_mb(PmxIm_a+l_type-1),
     $                          dbl_mb(PmnIm_a+l_type+1))
               endif
               nprim = log(maxprim/minprim)/log(dk_ratio) + 1
               alpha = (maxprim/minprim)**(1.0d0/dble(nprim-1))
               if (.not. bas_uce_tag(basis_in,iatom, usetag))
     $          call errquit('bas_create_uncont: bas_cont_tag failed',0,
     &       BASIS_ERR)
               do 21 iprim = 1, nprim
                  if (.not. bas_add_ucnt(basis_out,usetag,l_type,1,1,
     $                0.0d0,maxprim,1.0d0,1,name, .false.)) call errquit 
     $                ('bas_create_uncont: bas_add_ucnt',0, BASIS_ERR)
                  maxprim = maxprim / alpha
   21          continue
   20       continue
   17    continue
         if (.not. ma_pop_stack(Pmnm_a))
     $     call errquit('bas_create_uncont: ma__pop_stack Pmn array',0,
     &       MA_ERR)
         if (.not. ma_pop_stack(Pmxm_a))
     $     call errquit('bas_create_uncont: ma__pop_stack Pmx array',0,
     &       MA_ERR)
      else
         do 15 iprim = 1, nuniqprim
            if (.not. bas_uce_tag(basis_in,int_mb(AIm_a+iprim-1),
     $                             usetag)) 
     $         call errquit ('bas_create_uncont: bas_cont_tag failed',0,
     &       BASIS_ERR)
            if (.not. bas_add_ucnt(basis_out,usetag,
     $                            int_mb(LIm_a+iprim-1),1,1,0.0d0,
     $                            dbl_mb(ExpIm_a+iprim-1),1.0d0,1,name,
     $                            .false.))
     $         call errquit ('bas_create_uncont: bas_add_ucnt',0,
     &       BASIS_ERR)
   15    continue
      endif
c
      if (.not. bas_get_spherical(basis_in,spheri))
     $     call errquit('bas_create_uncont: bas_get_spherical',0,
     &       BASIS_ERR)
      if (.not. bas_set_spherical(basis_out,spheri))
     $     call errquit('bas_create_uncont: bas_set_spherical',0,
     &       BASIS_ERR)
c
      if (.not. ma_pop_stack(Am_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Sm_a',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(Sm_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Sm_a',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(Lm_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Lm_a',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(Expm_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Expm_a',0,
     &       MA_ERR)
c
c     Store the uncontracted basis on the RTDB
c
      if (.not. bas_rtdb_store(usertdb,name,basis_out))
     $   call errquit('bas_create_uncont: bas_rtdb_store failed',0,
     &       RTDB_ERR)
c
c     Do a delete for bas_add_ucnt_init with the strange name tidy
c
      if (.not. bas_add_ucnt_tidy(basis_out)) call errquit
     &      ('bas_create_uncont: bas_add_ucnt_tidy failed',0,
     &       BASIS_ERR)
c
      if (.not. bas_destroy(basis_out))
     $   call errquit('bas_create_uncont: bas_destroy failed',0,
     &       BASIS_ERR)
c
c     Load the uncontracted basis from the RTDB
c
      if (.not. bas_create(basis_out,name))
     $   call errquit('bas_create_uncont: bas_create failed',0,
     &       BASIS_ERR)
      if (.not. bas_rtdb_load(usertdb,usegeom,basis_out,name))
     $   call errquit('bas_create_uncont: bas_rtdb_load failed',0,
     &       RTDB_ERR)
c
c     Normalize the uncontracted fitting basis 
c
      call int_normalize(usertdb,basis_out)
c
      return
c
      end
c
      subroutine bas_unique_sort(rd,id1,id2,id3,length,iuniq)
      implicit none
c
      integer length, iuniq
      integer id1(length), id2(length), id3(length)
      double precision rd(length)
c
      integer iloop, jloop, kloop
c
      iuniq = length
    1 continue
      do 2 iloop = 1, iuniq
         do 3 jloop = iloop+1, iuniq 
            if ((dabs(rd(iloop)-rd(jloop)).lt.1.0d-5)
     $          .and.(id1(iloop).eq.id1(jloop)) .and.
     $          (id3(iloop).eq.id3(jloop))) then
               do 4 kloop = jloop+1, iuniq
                  rd(kloop-1) = rd(kloop)
                  id1(kloop-1) = id1(kloop)
                  id2(kloop-1) = id2(kloop)
                  id3(kloop-1) = id3(kloop)
    4          continue
               iuniq = iuniq - 1
               goto 1
            endif
    3    continue
    2 continue
c
      return
c
      end
