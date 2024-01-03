
      subroutine pm_localization(rtdb, geom, ltyp, basis, g_c, g_sc,
     &  g_sc0, nocc, nvir, nmo, nbf, natoms, eval, occ, c, sc, iloc,
     &  pop, list)

c     =================================================================
c     set up Pipek-Mezey localization

c     g_c, g_sc, g_sc0 and arrays occ and eval are for the MOs of a
c     given spin (alpha or beta)

c     The localization routine called from here is in
c     src/ddscf/localize.F
c     The PM localization code was already available there, but is
c     now interfaced more cleanly with the property code, and we added
c     the LMO sorting, support for spin-unrestricted SCF

c     =================================================================

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"


c     subroutine arguments:

      integer rtdb, geom, basis
      character*(3) ltyp
      integer g_c, g_sc, g_sc0
      integer nocc, nvir, nmo, nbf, natoms
      double precision eval(nbf), occ(nbf), c(nbf,2), sc(nbf,2)
      double precision pop(natoms)
      integer list(natoms)
      integer iloc(nmo)

c     local GA handles:

      integer g_t

c     local variables:

      logical master, debug
      integer dbg, info

      integer ncore, nloc

      integer s, ss, nlist, a, i, t, tt, u, bflo, bfhi

      double precision tmp, qas

      double precision small
      parameter (small=1d-8)

      double precision minval, swap

      character*(15) pname

c     =================================================================

      pname = 'pm_localization'

      dbg = 0
      master =  ga_nodeid().eq.0
      debug = (dbg>0) .and. master ! .true. during development

      if (ltyp.ne.'occ' .and. ltyp.ne.'vir')  call errquit
     &     (pname//': loc. type unknown', 0, BASIS_ERR)

      if(debug) then
        if (ltyp.eq.'occ') write (luout,*)
     &       'entering occupied PM localization'
        if (ltyp.eq.'vir') write (luout,*)
     &    'entering virtual PM localization'
      end if

      if (.not. geom_num_core(rtdb, geom, 'ddscf', ncore)) ncore = 0

      if (debug) write (luout,*) 'ncore = ',ncore

      if (ltyp.eq.'occ') then

c       localize the occupied MOs
        do i = 1, nocc
          iloc(i) = i
        end do
        nloc = nocc

        call localizePM(basis, c, sc,
     &    nloc, iloc, nbf, nmo, g_c, g_sc)


      else if (ltyp.eq.'vir') then

        do i = nocc+1, nmo
          iloc(i-nocc) = i
        end do
        nloc = nmo - nocc

        if(nloc.ne.nvir) call errquit(
     &    pname//': error nloc.ne.nvir',
     &    nloc-nvir, UNKNOWN_ERR)

        call localizePM(basis, c, sc,
     &    nloc, iloc, nbf, nmo, g_c, g_sc)


      end if ! ltyp?

c     calculate the localization transform

      if (.not. ga_create(MT_DBL, nmo, nmo, 'loc:trans',
     &  nloc, 0, g_t)) call errquit(pname//': g_t',1, GA_ERR)

      call ga_dgemm('t', 'n', nmo, nmo, nbf,
     &  1.0d0, g_sc0, g_c, 0.0d0, g_t)

c     transform MO energies, then sort the localized MOs. We assume that
c     the Fock matrix was diagonal in the basis of input MOs.

      if (ga_nodeid().eq.0) then

        do tt = 1,nloc
          t = iloc(tt)
          call ga_get(g_t, 1, nloc, t, t, c(1,1), 1)
          tmp = 0.0d0
          do ss = 1,nloc
            s = iloc(ss)
            tmp = tmp + eval(s) * c(s,1)**2
          end do ! s
          c(t,2) = tmp ! transformed MO energies
        end do ! t
        do ss = 1,nloc
          s = iloc(ss)
          eval(s) = c(s,2)
        end do

        do ss = 1, nloc -1
          s = iloc(ss)
          minval = eval(s)
c         find lowest eval(u) below eval(s)
          u = 0
          do tt = ss+1, nloc
            t = iloc(tt)
            if (eval(t).lt.minval) then
              u = t
              minval = eval(t)
            end if
          end do
c         if u > 0 we swap s and t
          if (u.ne.0) then
c           swap orbital energies
            swap = eval(s)
            eval(s) = eval(u)
            eval(u) = swap
c           swap LMOs
            call ga_get(g_c, 1, nbf, s,   s, c(1,1), 1)
            call ga_get(g_c, 1, nbf, u,   u, c(1,2), 1)
            call ga_put(g_c, 1, nbf, s,   s, c(1,2), 1)
            call ga_put(g_c, 1, nbf, u,   u, c(1,1), 1)
c           swap corresponding SC columns
            call ga_get(g_sc, 1, nbf, s,  s, c(1,1), 1)
            call ga_get(g_sc, 1, nbf, u,  u, c(1,2), 1)
            call ga_put(g_sc, 1, nbf, s,  s, c(1,2), 1)
            call ga_put(g_sc, 1, nbf, u,  u, c(1,1), 1)
c           swap columns of LMO transformation
            call ga_get(g_t, 1, nloc, s,  s, c(1,1), 1)
            call ga_get(g_t, 1, nloc, u,  u, c(1,2), 1)
            call ga_put(g_t, 1, nloc, s,  s, c(1,2), 1)
            call ga_put(g_t, 1, nloc, u,  u, c(1,1), 1)
          end if
        end do ! ss

      end if ! ga_nodeid.eq.0

      call ga_sync

      if (.not. ga_destroy(g_t)) call errquit(
     &  pname//': error destroying g_t',1, GA_ERR)

c     Analyze localization of each mo:
c     per lmo, a list of atomic populations is printed
c     in decreasing magnitude, with the populations in parentheses.

      if (ga_nodeid() .eq. 0) then
         write(luout,*)
         do ss = 1, nloc
            s = iloc(ss)
            call ga_get(g_c,  1, nbf, s, s, c(1,1), 1)
            call ga_get(g_sc, 1, nbf, s, s,sc(1,1), 1)
            nlist = 0
            do a = 1, natoms
               if (.not. bas_ce2bfr(basis, a, bflo, bfhi))
     &              call errquit(pname//': basis ', 0,
     &       BASIS_ERR)
               qas  = 0.0d0
               do u = bflo, bfhi
                  qas  = qas  + c(u,1)*sc(u,1)
               end do
               if (abs(qas) .gt. 0.01d0) then
                  nlist = nlist + 1
                  list(nlist) = a
                  pop(nlist) = qas
               end if
            end do
            do u = 1, nlist
               do t = 1, u-1
                  if (abs(pop(t)).lt.abs(pop(u))) then
                     tmp = pop(u)
                     pop(u) = pop(t)
                     pop(t) = tmp
                     tt = list(u)
                     list(u) = list(t)
                     list(t) = tt
                  end if
               end do
            end do

            write(luout,9002) s, eval(s),
     &        occ(s),(list(a), pop(a), a=1,nlist)
 9002       format(i5, 1x, f14.6,1x, f5.3, 1x,100(2x,i4,'(',f5.2,')'))
          end do
          write(luout,*)
         call util_flush(6)
      end if
c
      call ga_sync

      return

      end



