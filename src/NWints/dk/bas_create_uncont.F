      subroutine bas_create_uncont(usertdb,basis_in,name,basis_out)
c $Id: bas_create_uncont.F,v 1.9 2000-01-12 19:02:43 bert Exp $
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
c
c     Given a contracted basis, create a new uncontracted basis on rtdb
c
c     function declarations
c
      logical bas_add_ucnt_init
      logical bas_add_ucnt_tidy
      logical bas_add_ucnt
      logical bas_set_spherical
      external bas_add_ucnt_init
      external bas_add_ucnt_tidy
      external bas_add_ucnt
      external bas_set_spherical
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
      integer nuniqatom, iuniqcent   ! # of unique atoms
      integer nshell
      integer iatom, ishell, iprim
      integer ifirst, ilast, icont
      integer cartsph, l_type
      integer counter, nuniqprim
      integer icenter(nw_max_atom)
      character*16 usetag
      logical new, spheri
c
      integer Pmnm_a,PmnIm_a         ! Following indices are used for generating 
      integer Pmxm_a,PmxIm_a         ! a fitting basis set based on exponents and
      integer l_max,nprim,loctag     ! angular momenta
      double precision maxprim,minprim,alpha
c
c     Get geometry and rtdb handles for basis_in
c
      if (.not. bas_geom(basis_in,usegeom))
     $   call errquit('bas_create_uncont: bas_geom failed',0)
c
c     Create a new fitting basis set
c
      if (.not. bas_create(basis_out,name))
     $   call errquit('bas_create_uncont: bas_create failed',0)
c
c     Check if the user generated a fitting basis on input called
c     "D-K basis", if so load and return
c
      if (bas_rtdb_load(usertdb,usegeom,basis_out,'D-K basis')) return
c
c     Create a new uncontracted basis and initialize adding function
c
      if (.not. bas_add_ucnt_init(basis_out)) call errquit
     &      ('bas_create_uncont: bas_add_ucnt_init failed',0)
c
c     Get info from basis_in in order to extract the exponents
c
      if (.not. bas_numcont(basis_in,nshell))
     $   call errquit ('bas_create_uncont: bas_numcont failed',0)
c
c     Get unique centers and their shell range
c
      nuniqatom = 0
      do 5 ishell = 1, nshell
         if (.not. bas_cn2uce(basis_in,ishell,iuniqcent))
     $      call errquit ('bas_create_uncont: bas_continfo failed',0)
         if (nuniqatom .gt. 0) then
            new = .true.
            do 6 iatom = 1, nuniqatom
               if (iuniqcent .eq. icenter(iatom)) new = .false.
    6       continue
            if (new) then
               nuniqatom = nuniqatom +1
               icenter(nuniqatom) = iuniqcent
            endif
         else
            nuniqatom = 1
            icenter(1) = iuniqcent
         endif
    5 continue
c
c     Get maximum number of exponents for unique atoms
c
      ipr_mx = 0
      do 7 iatom = 1, nuniqatom
         if (.not.bas_ce2cnr(basis_in,icenter(iatom),ifirst,ilast))
     $      call errquit ('bas_create_uncont: bas_ce2cnr failed',0)
         do 8 ishell = ifirst, ilast
            if (.not. bas_continfo(basis_in,ishell,l_type,
     $                             nprim,icont,cartsph))
     $         call errquit ('bas_create_uncont: bas_continfo ',0)
            ipr_mx = ipr_mx + nprim
    8    continue
    7 continue
c
c     Make arrays for exponents plus l-types and shells of exponents
c
      if (.not. ma_push_get(MT_DBL,ipr_mx,'primitives',Expm_a,ExpIm_a))
     $     call errquit('bas_create_uncont: ma_push_get Exp array',0)
      if (.not. ma_push_get(MT_INT,ipr_mx,'ltypes',Lm_a,LIm_a))
     $     call errquit('bas_create_uncont: ma_push_get L array',0)
      if (.not. ma_push_get(MT_INT,ipr_mx,'shells',Sm_a,SIm_a))
     $     call errquit('bas_create_uncont: ma_push_get Exp array',0)
      if (.not. ma_push_get(MT_INT,ipr_mx,'atoms',Am_a,AIm_a))
     $     call errquit('bas_create_uncont: ma_push_get Exp array',0)
c
c     Extract exponents and data from basis_in for all unique atoms 
c
      counter = 0
      do 1 iatom = 1, nuniqatom
         if (.not.bas_ce2cnr(basis_in,icenter(iatom),ifirst,ilast))
     $      call errquit ('bas_create_uncont: bas_ce2cnr failed',0)
         do 2 ishell = ifirst, ilast
            if (.not. bas_continfo(basis_in,ishell,l_type,
     $                             nprim,icont,cartsph))
     $         call errquit ('bas_create_uncont: bas_continfo ',0)
            if (.not. bas_get_exponent(basis_in,ishell,
     $                                 dbl_mb(ExpIm_a+counter)))
     $         call errquit ('bas_create_uncont: bas_get_exponent',0)
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
     $       'Cannot use sp or spd to create flexi fitting basis DK',0)
    9    continue
c
c        Get highest angular momentum for array definitions and create arrays
c
         if (.not. bas_high_angular(basis_in,l_max))
     $      call errquit('bas_create_uncont: bas_high_angular failed',0)
         ipr_mx = l_max + 1
         if (.not. ma_push_get(MT_DBL,ipr_mx,'maxprim',Pmxm_a,PmxIm_a))
     $      call errquit('bas_create_uncont: ma_push_get Pmx array',0)
         if (.not. ma_push_get(MT_DBL,ipr_mx,'minprim',Pmnm_a,PmnIm_a))
     $      call errquit('bas_create_uncont: ma_push_get Pmn array',0)
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
            loctag = 0
            do 18 iprim = 1, nuniqprim
               if (int_mb(AIm_a+iprim-1) .eq. icenter(iatom)) then
                  if (l_max .lt. int_mb(LIm_a+iprim-1))
     $               l_max = int_mb(LIm_a+iprim-1)
               counter = int_mb(LIm_a+iprim-1)
               if (dbl_mb(PmxIm_a+counter) .lt. dbl_mb(ExpIm_a+iprim-1)) 
     $             dbl_mb(PmxIm_a+counter) = dbl_mb(ExpIm_a+iprim-1)
               if (dbl_mb(PmnIm_a+counter) .gt. dbl_mb(ExpIm_a+iprim-1)) 
     $             dbl_mb(PmnIm_a+counter) = dbl_mb(ExpIm_a+iprim-1)
               loctag = int_mb(SIm_a+iprim-1)
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
               if (.not. bas_cont_tag(basis_in,loctag, usetag))
     $          call errquit('bas_create_uncont: bas_cont_tag failed',0)
               do 21 iprim = 1, nprim
                  if (.not. bas_add_ucnt(basis_out,usetag,l_type,1,1,
     $                0.0d0,maxprim,1.0d0,1,name, .false.)) call errquit 
     $                ('bas_create_uncont: bas_add_ucnt',0)
                  maxprim = maxprim / alpha
   21          continue
   20       continue
   17    continue
         if (.not. ma_pop_stack(Pmnm_a))
     $     call errquit('bas_create_uncont: ma__pop_stack Pmn array',0)
         if (.not. ma_pop_stack(Pmxm_a))
     $     call errquit('bas_create_uncont: ma__pop_stack Pmx array',0)
      else
         do 15 iprim = 1, nuniqprim
            if (.not. bas_cont_tag(basis_in,int_mb(SIm_a+iprim-1),
     $                             usetag)) 
     $         call errquit ('bas_create_uncont: bas_cont_tag failed',0)
            if (.not. bas_add_ucnt(basis_out,usetag,
     $                            int_mb(LIm_a+iprim-1),1,1,0.0d0,
     $                            dbl_mb(ExpIm_a+iprim-1),1.0d0,1,name,
     $                            .false.))
     $         call errquit ('bas_create_uncont: bas_add_ucnt',0)
   15    continue
      endif
c
      spheri = .false.
      if (cartsph.eq.1) spheri = .true.
      if (.not. bas_set_spherical(basis_out,spheri))
     $     call errquit('bas_create_uncont: bas_set_spherical',0)
c
      if (.not. ma_pop_stack(Am_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Sm_a',0)
      if (.not. ma_pop_stack(Sm_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Sm_a',0)
      if (.not. ma_pop_stack(Lm_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Lm_a',0)
      if (.not. ma_pop_stack(Expm_a))
     $   call errquit('bas_create_uncont: ma_pop_stack Expm_a',0)
c
c     Store the uncontracted basis on the RTDB
c
      if (.not. bas_rtdb_store(usertdb,name,basis_out))
     $   call errquit('bas_create_uncont: bas_rtdb_store failed',0)
c
c     Do a delete for bas_add_ucnt_init with the strange name tidy
c
      if (.not. bas_add_ucnt_tidy(basis_out)) call errquit
     &      ('bas_create_uncont: bas_add_ucnt_tidy failed',0)
c
      if (.not. bas_destroy(basis_out))
     $   call errquit('bas_create_uncont: bas_destroy failed',0)
c
c     Load the uncontracted basis from the RTDB
c
      if (.not. bas_create(basis_out,name))
     $   call errquit('bas_create_uncont: bas_create failed',0)
      if (.not. bas_rtdb_load(usertdb,usegeom,basis_out,name))
     $   call errquit('bas_create_uncont: bas_rtdb_load failed',0)
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
            if ((dabs(rd(iloop)-rd(jloop)).lt.1.0d-8)
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
