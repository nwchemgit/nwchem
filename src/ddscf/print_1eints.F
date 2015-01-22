      subroutine print_1eints(ibas, jbas, integ_type, oskel)
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
c#include "util.fh"
c#include "stdio.fh"
c#include "tcgmsg.fh"
c#include "crohf.fh"
c#include "cscf.fh"
c#include "bas.fh"
c     
c     This is the original routine of print_1eints. It is still needed by the 
c     Douglas-Kroll routines to avoid recursive routine calling.
c
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     and ADD them into the given global array.
c     This version computes the full square of integrals and should work
c     OK even if ibas != jbas.
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... requires that ibas = jbas.
c     
c     arguments
c     
      integer ibas, jbas        ! [input] bra and ket basis sets
      character*(*) integ_type  ! [input] Name of integrals to compute
      character*10 inttitle
      logical oskel             ! [input] If true generate symmetry unique list
c     
c     local variables
c     
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
      logical odoit
      double precision q2
      external block_int1e      ! For T3D
      integer i, noffset,g_loc
c
c      logical ocache_save
c
      logical odbug
      logical osome
c
c     JEFF
c
      logical more_printing, status
      integer nbf
      character*8 cartbf(-1:7,1:36)
c
c     i... <=> ibas, j... <=> jbas
c
      integer ishdim, jshdim, ibf, jbf
      integer igeom, jgeom
      integer itype, jtype
      integer inprim, jnprim
      integer ingen, jngen, igen, jgen
      integer isphcart, jsphcart
      integer iatom, jatom
      character*16 itag, jtag
      character*2 isymbol, jsymbol
      character*16 ielement, jelement
      integer iatn, jatn
      integer inat, jnat
      character*8 icbf, jcbf
c      logical geom_cent_tag
c      external geom_cent_tag
c      logical geom_tag_to_element
c      external geom_tag_to_element
c
c     arguments of bas_continfo(ibas,ishell,itype,inprim,ingen,isphcart)
c
c      integer ibas              ! [input] basis handle
c      integer ishell            ! [input] contraction index
c      integer itype             ! [output] type (sp/s/p/d/..)
c      integer inprim            ! [output] no. of primitives
c      integer ingeni            ! [output] no. of contractions
c      integer isphcart          ! [output] 0/1 for cartesian/shperical
c
c     arguments of bas_cn2ce(ibas,ishell,iatom)
c
c      integer ibas              ! [input] basis set handle
c      integer ishell            ! [input] mapped contraction index
c      integer iatom             ! [output] center index
c
c     arguments of geom_cent_tag(geom, iatom, itag)
c
c      integer geom              ! [input]
c      integer iatom             ! [input]
c      character*16 itag         ! [output]
c
c     arguments of geom_tag_to_element(itag, isymbol, ielement, iatn)
c
c      character*16 itag         ! [input]
c      character*(*) isymbol     ! [output] (The two character abbreviation)
c      character*(*) ielement    ! [output] (The complete element name)
c      integer iatn              ! [output] (Useless?)
c
c     END JEFF
c
      odbug=.false.
      osome=.false.
      osome=osome.or.odbug
      odbug=odbug.and.(ga_nodeid().eq.0)
      osome=osome.and.(ga_nodeid().eq.0)
      if(osome) then
         write(6,*) 'in -int_1e_oldga- ... integ_type = ',
     $              integ_type,ga_nodeid()
         call util_flush(6)
      endif
c
      call ga_sync()
      if (oscfps) call pstat_on(ps_int_1e)
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $        ('print_1eints: use of symmetry requires ibas=jbas', ibas,
     &       BASIS_ERR)
      end if
c     
      if (inp_compare(.false., integ_type, 'kinetic')) then
         type = 1
      else if (inp_compare(.false., integ_type, 'potential')) then
         type = 2
      else if (inp_compare(.false., integ_type, 'overlap')) then
         type = 3
      else if (inp_compare(.false., integ_type, 'pvp')) then
         type = 4
      else if (inp_compare(.false., integ_type, 'so'))then 
         type = 5 
      else if (inp_compare(.false., integ_type, 'cos_chg_pot'))then 
         type = 6 
      else if (inp_compare(.false., integ_type, 'bq_pot'))then 
         type = 7 
      else
         write(6,*) ' integ_type = ', integ_type,ga_nodeid()
         call errquit('print_1eints: unknown integ_type', 0, INT_ERR)
      end if
c
c     ----- save ocache logical variable -----
c
ccccccccccccccccccccc
c      if(type.eq.6 .or. type.eq.7 ) then
c         ocache_save=ocache
c         ocache     =.false.
c      endif
c
c     In-core caching
c
c      ocache = .false.
c      if(osome) then
c         write(6,*) 'ocache = ',ocache,ga_nodeid()
c         call util_flush(6)
c      endif
c      if (ocache .and. (ibas.eq.bas1) .and. (jbas.eq.bas1)) then
c*         if (ga_nodeid().eq.0) then
c*           call ga_summarize(1)
c*           call util_flush(6)
c*         endif
c*         call ga_sync()
c         call ga_dadd(1.0d0, g_cache(type), 1.0d0, g, g)
c         if (oscfps) call pstat_off(ps_int_1e)
c	 return
c      endif
ccccccccccccccccccccccc
c     
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('print_1eints: bas_numcont failed for ibas', ibas,
     &       BASIS_ERR)
      if (.not. bas_numcont(jbas, nshell_j)) call errquit
     $     ('print_1eints: bas_numcont failed for jbas', jbas,
     &       BASIS_ERR)
c
c     Get the geometry handle (JEFF)
c
      if (.not. bas_geom(ibas, igeom)) call errquit
     $     ('print_1eints: geom?',igeom,BASIS_ERR)
      if (.not. geom_ncent(igeom, inat)) call errquit
     $     ('print_1eints: ncent?',inat, BASIS_ERR)
c
      if (.not. bas_geom(jbas, jgeom)) call errquit
     $     ('print_1eints: geom?',jgeom,BASIS_ERR)
      if (.not. geom_ncent(jgeom, jnat)) call errquit
     $     ('print_1eints: ncent?',jnat, BASIS_ERR)
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c     
c     k_* are the offsets corrsponding to the l_* handles
c     
cso
      call int_mem_1e(max1e, mem1)
      if (type .eq. 5) then  
         max1e = 3*max1e 
         mem1 = 3*mem1 
      endif 
cso
      if (.not. MA_push_get(MT_DBL,max1e,'print_1eints:buf',
     &     l_buf,k_buf)) call errquit('print_1eints: ma failed',
     &     max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'print_1eints:scr',
     &     l_scr,k_scr)) call errquit('print_1eints: ma failed',
     &     mem1, MA_ERR)
c     
c     Loop thru shells with static parallel work decomposition
c     
ccc      if (.not.inp_compare(.false., integ_type, 'so')) then
ccc         if(.not.ga_duplicate(g,g_loc,'local g')) call
ccc     .        errquit('int1e: dupl failed',0, GA_ERR)
ccc         call ga_zero(g_loc)
ccc      endif
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      q2 = 1.0d0
c
c
      call print_1e_init_cartbf(cartbf)
c
c
      if(iproc.eq.0) then
      do jshell = 1, nshell_j
         do ishell = 1, nshell_i
c
cprog            if (mod(ijshell, nproc) .eq. iproc) then
              odoit = .true.
              if (oskel) odoit = sym_shell_pair(ibas,ishell,jshell,q2)
cnosym              if (odoit) then
                if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) 
     $           call errquit('print_1eints: bas_cn2bfr',ibas,BASIS_ERR)
                if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $           call errquit('print_1eints: bas_cn2bfr',jbas,BASIS_ERR)
c
                if (.not. bas_continfo(ibas,ishell,itype,
     &                                 inprim,ingen,isphcart))
     &            call errquit('print_1eints: continf',ishell,BASIS_ERR)
                if (.not. bas_continfo(jbas,jshell,jtype,
     &                                 jnprim,jngen,jsphcart))
     &            call errquit('print_1eints: continf',jshell,BASIS_ERR)
c
                if (.not. bas_cn2ce(ibas, ishell, iatom))
     &            call errquit('print_1eints: cn2ce',ishell,BASIS_ERR)
                if (.not. bas_cn2ce(jbas, jshell, jatom))
     &            call errquit('print_1eints: cn2ce',jshell,BASIS_ERR)
c
                status = geom_cent_tag(igeom, iatom, itag)
                status = geom_cent_tag(jgeom, jatom, jtag)
c
                status = geom_tag_to_element(itag,isymbol,ielement,iatn)
                status = geom_tag_to_element(jtag,jsymbol,jelement,jatn)
c
                ishdim = ihi - ilo + 1
                jshdim = jhi - jlo + 1
c
c
c
         do igen = 1, ingen
            do ibf = 1, ishdim
               if (isphcart .eq. 0) then
                  icbf = cartbf(itype,ibf)
               else
                  icbf = ' '
                  write(icbf,'(i3)') ibf-ishdim/2-1
               endif
c               write(6,*) 'igen = ',igen
c               write(6,*) 'ibf = ',ibf
c               write(6,*) 'icbf = ',icbf
            enddo
         enddo
c
         do jgen = 1, jngen
            do jbf = 1, jshdim
               if (jsphcart .eq. 0) then
                  jcbf = cartbf(jtype,jbf)
               else
                  jcbf = ' '
                  write(jcbf,'(i3)') jbf-jshdim/2-1
               endif
            enddo
         enddo
c
c
c
c     
c     Generate the integrals
c     
                if (type .eq. 1) then
                   inttitle = '1eke'
                   call int_1eke (jbas, jshell, ibas, ishell, 
     $                  mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                else if (type .eq. 2) then
                   inttitle = '1epe'
                   call int_1epe (jbas, jshell, ibas, ishell, 
     $                  mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                else if (type .eq. 3) then
                   inttitle = '1eov'
                   call int_1eov (jbas, jshell, ibas, ishell, 
     $                  mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                else if (type .eq. 4) then
                   inttitle = '1epvpe'
                   call int_1epvpe (jbas, jshell, ibas, ishell, 
     $                  mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                else if (type .eq. 5) then 
                   inttitle = 'so_1e'
                   call intso_1e (jbas, jshell, ibas, ishell,
     $                  mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                else if (type .eq. 6) then 
                   if(odbug) then
                      write(6,*) 'type = 6 ... potential ...',
     $                           ga_nodeid()
                      call util_flush(6)
                   endif
                   inttitle = '1epot'
                   call int_1epot (jbas, jshell, ibas, ishell,
     $                  mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                else if (type .eq. 7) then 
                   inttitle = '1epot1'
                   call int_1epot1 (jbas, jshell, ibas, ishell,
     $                  mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                else
                   call errquit('print_1eints: invalid type?', type,
     &       GA_ERR)
                end if
c     
c     Print the integrals
c    
                more_printing = .true.
                if (.not.more_printing) then
                  call print_1e_buf(ilo, ihi,jlo,jhi,
     &                 dbl_mb(k_buf),inttitle)
                else
c
c                 Printing atomic basis set information
c
                  call print_1e_buf_more(ilo,ihi,jlo,jhi,
     &                 dbl_mb(k_buf),inttitle,ishell,jshell,
     &                 iatom,jatom,isymbol,jsymbol,itype,jtype,
     &                 icbf,jcbf)
c
c                 The end
c
                endif
cnosym              end if
c            endif !proc
            ijshell = ijshell + 1
         end do
      end do
      endif ! iproc.eq.0
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('print_1eints: pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('print_1eints: pop failed', 0, GA_ERR)
c
c     ----- restore ocache -----
c
cccccccccccccccc
c      if(type.eq.6 .or. type.eq.7) then
c         ocache=ocache_save
c      endif
c     
c      if (.not.inp_compare(.false., integ_type, 'so')) then
c         call ga_dadd(1.0d0, g_loc, 1.0d0, g, g)
c         if (.not. ga_destroy(g_loc)) call errquit('i1ega: gad?',0,
c     &       GA_ERR)
c      endif
cccccccccccccccc
      call ga_sync()            ! So that no nasty races can result
c     
      if (oscfps) call pstat_off(ps_int_1e)
c     
      end
c
      subroutine print_1e_buf(ilo, ihi, jlo, jhi, buf,inttitle)
      implicit none
      integer i, j
      integer ilo, ihi, jlo, jhi
      double precision buf(ilo:ihi,jlo:jhi)
      character*10 inttitle
c
      do j = jlo, jhi
         do i = ilo, ihi
            if (abs(buf(i,j)).gt.1d-8) then
               write(6,1) inttitle,i,j,buf(i,j)
 1             format(A,2i5,f20.8)
            endif
         enddo
      enddo
c
      end
c
      subroutine print_1e_buf_more(ilo,ihi,jlo,jhi,
     &                   buf,title,ishell,jshell,iatom,jatom,
     &                   isymbol,jsymbol,itype,jtype,icbf,jcbf)

      implicit none
      integer i, j
      integer ilo, ihi, jlo, jhi
c      double precision buf(jlo:jhi,ilo:ihi)
      double precision buf(ilo:ihi,jlo:jhi)
      character*10 title
      character*8 cartbf(-1:7,1:36)
      integer ishell, jshell
      integer itype, jtype
      integer iatom, jatom
      character*16 itag, jtag
      character*2 isymbol, jsymbol
      character*8 icart, jcart
      character*8 icbf, jcbf
c
      character*1 shtypes(-1:7)
      data shtypes/'l','s', 'p', 'd', 'f', 'g', 'h', 'i', 'k'/
c

c
      call print_1e_init_cartbf(cartbf)
c
      do j = jlo, jhi
         do i = ilo, ihi
            if (abs(buf(i,j)).gt.1d-8) then
c               write(6,*) 'ishell = ',ishell
c               write(6,*) 'jshell = ',jshell
c               write(6,*) 'iatom = ',iatom
c               write(6,*) 'jatom = ',jatom
c               write(6,*) 'isymbol = ',isymbol
c               write(6,*) 'jsymbol = ',jsymbol
c               write(6,*) 'itype = ',itype
c               write(6,*) 'jtype = ',jtype
c               write(6,*) 'icbf = ',icbf
c               write(6,*) 'jcbf = ',jcbf
c               icart = '        ' ! cartbf(itype,ishell)
c               jcart = '        ' ! cartbf(jtype,jshell)
               write(6,1) title,i,isymbol,shtypes(itype),
     &                          j,jsymbol,shtypes(jtype),buf(i,j)
 1             format(A,i5,3x,a2,1x,a1,1x,
     &                  i5,3x,a2,1x,a1,1x,f20.8)
            endif
         enddo
      enddo
c
      end
c
      subroutine print_1e_init_cartbf(cartbf)
      implicit none
c
      character*8 cartbf(-1:7,1:36)
c
      cartbf(-1,1) = 's'
      cartbf(-1,2) = 'x'
      cartbf(-1,3) = 'y'
      cartbf(-1,4) = 'z'

      cartbf(0, 1) = ' '

      cartbf(1, 1) = 'x'
      cartbf(1, 2) = 'y'
      cartbf(1, 3) = 'z'

      cartbf(2, 1) = 'xx'
      cartbf(2, 2) = 'xy'
      cartbf(2, 3) = 'xz'
      cartbf(2, 4) = 'yy'
      cartbf(2, 5) = 'yz'
      cartbf(2, 6) = 'zz'

      cartbf(3, 1) = 'xxx'
      cartbf(3, 2) = 'xxy'
      cartbf(3, 3) = 'xxz'
      cartbf(3, 4) = 'xyy'
      cartbf(3, 5) = 'xyz'
      cartbf(3, 6) = 'xzz'
      cartbf(3, 7) = 'yyy'
      cartbf(3, 8) = 'yyz'
      cartbf(3, 9) = 'yzz'
      cartbf(3,10) = 'zzz'

      cartbf(4, 1) = 'xxxx'
      cartbf(4, 2) = 'xxxy'
      cartbf(4, 3) = 'xxxz'
      cartbf(4, 4) = 'xxyy'
      cartbf(4, 5) = 'xxyz'
      cartbf(4, 6) = 'xxzz'
      cartbf(4, 7) = 'xyyy'
      cartbf(4, 8) = 'xyyz'
      cartbf(4, 9) = 'xyzz'
      cartbf(4,10) = 'xzzz'
      cartbf(4,11) = 'yyyy'
      cartbf(4,12) = 'yyyz'
      cartbf(4,13) = 'yyzz'
      cartbf(4,14) = 'yzzz'
      cartbf(4,15) = 'zzzz'

      cartbf(5, 1) = 'xxxxx'
      cartbf(5, 2) = 'xxxxy'
      cartbf(5, 3) = 'xxxxz'
      cartbf(5, 4) = 'xxxyy'
      cartbf(5, 5) = 'xxxyz'
      cartbf(5, 6) = 'xxxzz'
      cartbf(5, 7) = 'xxyyy'
      cartbf(5, 8) = 'xxyyz'
      cartbf(5, 9) = 'xxyzz'
      cartbf(5,10) = 'xxzzz'
      cartbf(5,11) = 'xyyyy'
      cartbf(5,12) = 'xyyyz'
      cartbf(5,13) = 'xyyzz'
      cartbf(5,14) = 'xyzzz'
      cartbf(5,15) = 'xzzzz'
      cartbf(5,16) = 'yyyyy'
      cartbf(5,17) = 'yyyyz'
      cartbf(5,18) = 'yyyzz'
      cartbf(5,19) = 'yyzzz'
      cartbf(5,20) = 'yzzzz'
      cartbf(5,21) = 'zzzzz'

      cartbf(6, 1) = 'xxxxxx'
      cartbf(6, 2) = 'xxxxxy'
      cartbf(6, 3) = 'xxxxxz'
      cartbf(6, 4) = 'xxxxyy'
      cartbf(6, 5) = 'xxxxyz'
      cartbf(6, 6) = 'xxxxzz'
      cartbf(6, 7) = 'xxxyyy'
      cartbf(6, 8) = 'xxxyyz'
      cartbf(6, 9) = 'xxxyzz'
      cartbf(6,10) = 'xxxzzz'
      cartbf(6,11) = 'xxyyyy'
      cartbf(6,12) = 'xxyyyz'
      cartbf(6,13) = 'xxyyzz'
      cartbf(6,14) = 'xxyzzz'
      cartbf(6,15) = 'xxzzzz'
      cartbf(6,16) = 'xyyyyy'
      cartbf(6,17) = 'xyyyyz'
      cartbf(6,18) = 'xyyyzz'
      cartbf(6,19) = 'xyyzzz'
      cartbf(6,20) = 'xyzzzz'
      cartbf(6,21) = 'xzzzzz'
      cartbf(6,22) = 'yyyyyy'
      cartbf(6,23) = 'yyyyyz'
      cartbf(6,24) = 'yyyyzz'
      cartbf(6,25) = 'yyyzzz'
      cartbf(6,26) = 'yyzzzz'
      cartbf(6,27) = 'yzzzzz'
      cartbf(6,28) = 'zzzzzz'

      cartbf(7, 1) = 'xxxxxxx'
      cartbf(7, 2) = 'xxxxxxy'
      cartbf(7, 3) = 'xxxxxxz'
      cartbf(7, 4) = 'xxxxxyy'
      cartbf(7, 5) = 'xxxxxyz'
      cartbf(7, 6) = 'xxxxxzz'
      cartbf(7, 7) = 'xxxxyyy'
      cartbf(7, 8) = 'xxxxyyz'
      cartbf(7, 9) = 'xxxxyzz'
      cartbf(7,10) = 'xxxxzzz'
      cartbf(7,11) = 'xxxyyyy'
      cartbf(7,12) = 'xxxyyyz'
      cartbf(7,13) = 'xxxyyzz'
      cartbf(7,14) = 'xxxyzzz'
      cartbf(7,15) = 'xxxzzzz'
      cartbf(7,16) = 'xxyyyyy'
      cartbf(7,17) = 'xxyyyyz'
      cartbf(7,18) = 'xxyyyzz'
      cartbf(7,19) = 'xxyyzzz'
      cartbf(7,20) = 'xxyzzzz'
      cartbf(7,21) = 'xxzzzzz'
      cartbf(7,22) = 'xyyyyyy'
      cartbf(7,23) = 'xyyyyyz'
      cartbf(7,24) = 'xyyyyzz'
      cartbf(7,25) = 'xyyyzzz'
      cartbf(7,26) = 'xyyzzzz'
      cartbf(7,27) = 'xyzzzzz'
      cartbf(7,28) = 'xzzzzzz'
      cartbf(7,29) = 'yyyyyyy'
      cartbf(7,30) = 'yyyyyyz'
      cartbf(7,31) = 'yyyyyzz'
      cartbf(7,32) = 'yyyyzzz'
      cartbf(7,33) = 'yyyzzzz'
      cartbf(7,34) = 'yyzzzzz'
      cartbf(7,35) = 'yzzzzzz'
      cartbf(7,36) = 'zzzzzzz'

      end


c $Id$
