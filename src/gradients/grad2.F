      subroutine grad2(d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $  d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $  d_ij3, d_kl3, d_ik3, d_jl3, d_il3, d_jk3, 
     $  d_ij4, d_kl4, d_ik4, d_jl4, d_il4, d_jk4, 
     $  frc, g_dens, g_force, ldim,
     $  rtdb, geom, basis, 
     $  nproc, nat, scr, lscr, buf, lbuf, tol2e,
     $  s34, oactive, oskel, scftype, omp2,
     $  pdm2,pdm2a,pdm2b,pdm2c,pdm2d,coeff,maxsh,nopen,nbf)
c     $Id: grad2.F,v 1.20 1997-02-07 22:06:12 d3g681 Exp $

      implicit none
c     
c     NOTE FOR EMACS USERS ... due to the deep nesting of do/if
c     in this routine indent it using identation variables
c     IF=1, DO=1, CONTINUATION=2
c     

#include "sym.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "rtdb.fh"

C-----------------------parameters---------------------------------
      integer ldim,             ! max dimension of density block
     $  g_dens(3),              ! array of global density handles
     $  g_force,                ! total gradients(global)
     $  rtdb, basis, geom,      ! handles
     $  nproc, nat, lscr, lbuf

c     FOR MCSCF ONLY
      integer nopen, maxsh, nbf
      double precision pdm2(nopen,nopen,nopen,nopen)
      double precision pdm2a(nopen,nopen,nopen,maxsh)
      double precision pdm2b(nopen,nopen,maxsh,maxsh)
      double precision pdm2c(nopen,maxsh,maxsh,maxsh)
      double precision pdm2d(maxsh,maxsh,maxsh,maxsh)
      double precision coeff(nbf,nopen)
c     

      double precision tol2e

      double precision scr(lscr), buf(lbuf), frc(3, nat), s34(nat),
     $  d_ij(ldim,ldim),d_kl(ldim,ldim),d_ik(ldim,ldim),
     $  d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim),
     $  d_ij2(ldim,ldim), d_kl2(ldim,ldim),
     $  d_ik2(ldim,ldim), d_jl2(ldim,ldim),
     $  d_il2(ldim,ldim), d_jk2(ldim,ldim),
     $  d_ij3(ldim,ldim), d_kl3(ldim,ldim),
     $  d_ik3(ldim,ldim), d_jl3(ldim,ldim),
     $  d_il3(ldim,ldim), d_jk3(ldim,ldim),
     $  d_ij4(ldim,ldim), d_kl4(ldim,ldim),
     $  d_ik4(ldim,ldim), d_jl4(ldim,ldim),
     $  d_il4(ldim,ldim), d_jk4(ldim,ldim)


      logical oactive(nat), oskel, omp2

      double precision scale

      character*(*) scftype

C-------------------------local variables-----------------------------
      integer i, ijatom, next, icart, idatom(4), icount,
     $  iat1, iat2, iat3, iat4, ihi4, ish1, ish2, ish3, ish4,
     $  iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l, 
     $  iac1f, iac1l, iac2f, iac2l, iac3f, iac3l, iac4f, iac4l,
     $  itop3, itop4, if1, il1, if2, il2, if3, il3, if4, il4,
     $  ip1, ip2, ip3, ip4, nshbfp, nshbfq, Nint,
     $  chktime, ndone
      integer id1, id2, id3, id4

      integer nxtask, task_size
      external nxtask

      double precision dE, den2, q4, q2, s12, half, quart, eigh, sixh

      logical status, oij, okl, oikjl, oschwarz12, osym, ouhf, orohf,
     $  ochkpt, oskip, ofinish, grad_restart, omcscf
      external grad_restart

      parameter(half = 0.5D0,
     $  quart= 0.25D0,
     $  eigh = 0.125D0,
     $  sixh = 0.0625D0)

      
#if defined(REPORT_LOOP)
      logical reportit
      character*26 date_string
#endif

      if (.not. bas_geom(basis, geom)) 
     $  call errquit('grad2: bad basis',0)

c     MCSCF is the same as RHF except for the 2pdm piece

      omcscf= (scftype .eq. 'MCSCF')
      ouhf  = (scftype .eq. 'UHF')
      orohf = (scftype .eq. 'ROHF')

      q4 = 1.0d0                ! If not using symmetry

      call fock_2e_est_task(geom, oskel, tol2e, task_size)

C     is this a restart run?
      oskip = grad_restart(rtdb, geom, basis, nat, ndone, frc, 
     $  g_force)

C     checkpointing requested?
      ochkpt = rtdb_get(rtdb, 'gradients:checkpoint', mt_int, 1, 
     $  chktime)

      if (ochkpt) then
C     convert to seconds
       chktime = chktime * 60
       ofinish = .false.
       call grad_store(rtdb, geom, g_force, chktime, nat, 0, frc, 
     $   .true., ofinish, oactive, oskel, omp2)
      end if

      if (ochkpt .or. oskip) then
C     at the moment checkpointing works only with task_size == 1
C     override fock_2e_est_task
       task_size = 1
      end if

      ijatom = -1
      next = nxtask(nproc,task_size)
      do iat3 = nat, 1, -1
#if defined(REPORT_LOOP)
       if (ga_nodeid().eq.0) then
        reportit = mod(iat3,10).eq.0
        reportit = reportit.or.(iat3.eq.1)
        reportit = reportit.or.(iat3.eq.nat)
        if (reportit) then
         call util_date(date_string)
         write(6,*) 'grad2: atom =',iat3,' ',date_string
         call util_flush(6)
        end if
       end if
#endif
       do iat4 = 1, iat3
        s34(iat4) = schwarz_atom(iat3,iat4) ! Precompute
       end do
       
       do iat1 = nat, iat3, -1  ! Want canonical atom list for symmetry
        do iat2 = iat1, 1, -1

         s12 = schwarz_atom(iat1,iat2)
         oschwarz12 = (s12*schwarz_max().ge.tol2e)
C     
         if (oschwarz12 .and. oskel) then
          oschwarz12 = sym_atom_pair(geom, iat1, iat2, q2)
         end if
c     
         if (oschwarz12) then
          ihi4 = iat3
          if (iat1.eq.iat3) ihi4 = iat2
          do iat4 = 1, ihi4

C     check if all atoms are inactive
           if (.not.oactive(iat1) .and. .not.oactive(iat2) .and.
     $       .not.oactive(iat3) .and. .not.oactive(iat4))
     $       goto 1020

           if (s12 * s34(iat4) .lt. tol2e) goto 1020

           if (oskel) then
            osym = sym_atom_quartet(geom,
     $        iat1, iat2, iat3, iat4, q4)
            if (.not. osym) goto 1020
           end if

           status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2bfr',0)
           status = bas_ce2bfr(basis,iat2,iab2f,iab2l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2bfr',0)
           status = bas_ce2bfr(basis,iat3,iab3f,iab3l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2bfr',0)
           status = bas_ce2bfr(basis,iat4,iab4f,iab4l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2bfr',0)

           if (iab1f.le.0 .or. iab2f.le.0 .or. iab3f.le.0 .or.
     $       iab4f.le.0) then
c     
c     At least one center has no functions on it ... next atom
c     
            goto 1020
           end if

           status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2cnr',0)
           status = bas_ce2cnr(basis,iat2,iac2f,iac2l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2cnr',0)
           status = bas_ce2cnr(basis,iat3,iac3f,iac3l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2cnr',0)
           status = bas_ce2cnr(basis,iat4,iac4f,iac4l)
           if (.not.status) 
     $       call errquit('grad2: error with bas_ce2cnr',0)

           ijatom = ijatom + 1
C     in case of restart skip loop until last calculated atom quadruple
           if (oskip) then
            if (ijatom .eq. ndone) then
             oskip = .false.
             ijatom = -1        ! nxtval(or nxtask) starts with 0
C     compensate for ijatom starting with 0
             ndone = ndone + 1
            end if
            goto 1020
           end if

           if (next .eq. ijatom) then

C     get density matrix blocks
            if (.not. omp2) then
             if (ouhf) then
              call grad_get_dens_uhf(d_ij, d_kl, d_ik, 
     $          d_jl, d_il, d_jk,
     $          d_ik2, d_jl2, d_il2, d_jk2, ldim, 
     $          iab1f, iab1l, iab2f, iab2l, 
     $          iab3f, iab3l, iab4f, iab4l,
     $          g_dens)
             else if (orohf) then
              call grad_get_dens_rohf(d_ij, d_kl, d_ik, 
     $          d_jl, d_il, d_jk,
     $          d_ij2, d_kl2, d_ik2, 
     $          d_jl2, d_il2, d_jk2, ldim, 
     $          iab1f, iab1l, iab2f, iab2l, 
     $          iab3f, iab3l, iab4f, iab4l,
     $          g_dens)
             else
              call grad_get_dens(d_ij, d_kl, d_ik, 
     $          d_jl, d_il, d_jk, ldim,
     $          iab1f, iab1l, iab2f, iab2l, 
     $          iab3f, iab3l, iab4f, iab4l,
     $          g_dens(1))
             end if

            else                ! MP2

             if (ouhf) then
              call grad_get_dens_ump2(d_ij, d_kl, d_ik, 
     $          d_jl, d_il, d_jk,
     $          d_ij2, d_kl2, d_ik2, 
     $          d_jl2, d_il2, d_jk2, 
     $          d_ij3, d_kl3, d_ik3, 
     $          d_jl3, d_il3, d_jk3, 
     $          d_ij4, d_kl4, d_ik4, 
     $          d_jl4, d_il4, d_jk4, 
     $          ldim,
     $          iab1f, iab1l, iab2f, iab2l, 
     $          iab3f, iab3l, iab4f, iab4l,
     $          g_dens)
             else if (orohf) then
              call errquit('no RMP2 gradients yet', 0)
             else               ! rmp2
              call grad_get_dens_rmp2(d_ij, d_kl, d_ik, 
     $          d_jl, d_il, d_jk,
     $          d_ij2, d_kl2, d_ik2, 
     $          d_jl2, d_il2, d_jk2, 
     $          ldim,
     $          iab1f, iab1l, iab2f, iab2l, 
     $          iab3f, iab3l, iab4f, iab4l,
     $          g_dens)
             end if

            end if

            oij = iat1.eq.iat2
            okl = iat3.eq.iat4
            oikjl = iat1.eq.iat3 .and. iat2.eq.iat4
            do ish1 = iac1f, iac1l
             status = bas_cn2bfr(basis,ish1,if1,il1)
             if (omcscf) call grad_mcscf_a
     $         (nbf,nopen,maxsh,if1,il1,coeff,pdm2,pdm2a)
             if1 = if1 - iab1f + 1
             il1 = il1 - iab1f + 1
             id1 = il1 - if1   + 1
c     
             if (oij) iac2l = ish1
             do ish2 = iac2f, iac2l
              if (schwarz_shell(ish1,ish2)*schwarz_max()*q4
     $          .ge. tol2e) then
               status = bas_cn2bfr(basis,ish2,if2,il2)
               if (omcscf) call grad_mcscf_b
     $           (nbf,nopen,maxsh,id1,if2,il2,coeff,
     $           pdm2a,pdm2b)
               if2 = if2 - iab2f + 1
               il2 = il2 - iab2f + 1
               id2 = il2 - if2   + 1
               itop3 = iac3l
               if (iat1.eq.iat3 .and. iat2.eq.iat4) 
     $           itop3 = ish1
               do ish3 = iac3f, itop3
                status = bas_cn2bfr(basis,ish3,if3,il3)
                if (omcscf) call grad_mcscf_c
     $            (nbf,nopen,maxsh,id1,id2,if3,il3,coeff,
     $            pdm2b,pdm2c)
                if3 = if3 - iab3f + 1
                il3 = il3 - iab3f + 1
                id3 = il3 - if3   + 1
                itop4 = iac4l
                if (okl) itop4 = ish3
                if (oikjl .and. ish1.eq.ish3) itop4 = ish2
                do ish4 = iac4f, itop4
                 if (schwarz_shell(ish1,ish2)
     $             * schwarz_shell(ish3,ish4) 
     $             * q4 .ge. tol2e) 
     $             then
                  status = bas_cn2bfr(basis,ish4,if4,il4)
                  if (omcscf) call grad_mcscf_d
     $              (nbf,nopen,maxsh,id1,id2,id3,
     $              if4,il4,coeff,pdm2c,pdm2d)
                  if4 = if4 - iab4f + 1
                  il4 = il4 - iab4f + 1
                  id4 = il4 - if4   + 1

C---------------------------------------------------------------------


C     find contraction in atomic block


c     Determine the number of ERIs in this block of integrals.

                  nshbfp = id1*id2
                  nshbfq = id3*id4
                  Nint = nshbfp*nshbfq

C     2el. integral derivatives
                  call intd_2e4c(basis, ish1, ish2, basis, ish3, ish4,
     &              lscr, scr, lbuf, buf, idatom)
                  
                  scale = 1.0d0
                  if (.not.oij .or. ish1.ne.ish2) scale = scale*2.0d0
                  if (.not.okl .or. ish3.ne.ish4) scale = scale*2.0d0
                  if (.not.oikjl .or. ish1.ne.ish3 .or. 
     $              ish2.ne.ish4) scale = scale*2.0d0

                  icount = 1
                  do i = 1, 4
                   if (idatom(i).ge.1) then
                    if (.not. omp2) then
                     if (ouhf) then
                      do icart = 1, 3
                       dE = 0.D0
                       do ip1 = if1, il1
                        do ip2 = if2, il2
                         do ip3 = if3, il3
                          do ip4 = if4, il4
                           den2 = half *(d_ij(ip1,ip2) * d_kl(ip3,ip4))
     $                       - eigh  *(d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                       +  d_il(ip1,ip4) * d_jk(ip2,ip3)
     $                       +  d_ik2(ip1,ip3) * d_jl2(ip2,ip4)
     $                       +  d_il2(ip1,ip4) * d_jk2(ip2,ip3))
                           dE = dE + den2 * buf(icount)
                           icount = icount + 1
                          end do
                         end do
                        end do
                       end do
                       frc(icart,idatom(i)) = frc(icart,idatom(i)) + 
     $                   dE * q4 * scale
                      end do
                     else if (orohf) then
                      do icart = 1, 3
                       dE = 0.D0
                       do ip1 = if1, il1
                        do ip2 = if2, il2
                         do ip3 = if3, il3
                          do ip4 = if4, il4
                           den2 = half *(d_ij(ip1,ip2) * d_kl(ip3,ip4))
     $                       + quart *(d_ij(ip1,ip2) * d_kl2(ip3,ip4)
     $                       +  d_ij2(ip1,ip2) * d_kl(ip3,ip4))
     $                       + eigh *(d_ij2(ip1,ip2) * d_kl2(ip3,ip4))
     $                       - eigh *(d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                       +  d_il(ip1,ip4) * d_jk(ip2,ip3))
     $                       - sixh *(d_ik(ip1,ip3) * d_jl2(ip2,ip4)
     $                       +  d_il2(ip1,ip4) * d_jk(ip2,ip3)
     $                       +  d_ik2(ip1,ip3) * d_jl(ip2,ip4)
     $                       +  d_il(ip1,ip4) * d_jk2(ip2,ip3)
     $                       +  d_ik2(ip1,ip3) * d_jl2(ip2,ip4)
     $                       +  d_il2(ip1,ip4) * d_jk2(ip2,ip3))
                           dE = dE + den2 * buf(icount)
                           icount = icount + 1
                          end do
                         end do
                        end do
                       end do
                       frc(icart,idatom(i)) = frc(icart,idatom(i)) + 
     $                   dE * q4 * scale
                      end do
                     else if (omcscf) then
                      do icart = 1, 3
                       dE = 0.D0
                       do ip1 = if1, il1
                        do ip2 = if2, il2
                         do ip3 = if3, il3
                          do ip4 = if4, il4
                           den2 = half *(d_ij(ip1,ip2)*d_kl(ip3,ip4))
     $                       - eigh * (d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                       +  d_il(ip1,ip4) * d_jk(ip2,ip3))
     $                       + 0.5d0*pdm2d(ip4-if4+1,ip3-if3+1,
     $                       ip2-if2+1,ip1-if1+1)
                           dE = dE + den2 * buf(icount)
                           icount = icount + 1
                          end do
                         end do
                        end do
                       end do
                       frc(icart,idatom(i)) = frc(icart,idatom(i)) + 
     $                   dE * q4 * scale
                      end do
                     else       ! RHF
                      do icart = 1, 3
                       dE = 0.D0
                       do ip1 = if1, il1
                        do ip2 = if2, il2
                         do ip3 = if3, il3
                          do ip4 = if4, il4
                           den2 = half *(d_ij(ip1,ip2)*d_kl(ip3,ip4))
     $                       - eigh *(d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                       +  d_il(ip1,ip4) * d_jk(ip2,ip3))
                           dE = dE + den2 * buf(icount)
                           icount = icount + 1
                          end do
                         end do
                        end do
                       end do
                       frc(icart,idatom(i)) = frc(icart,idatom(i)) + 
     $                   dE * q4 * scale
                      end do
                     end if      ! rhf

                    else        ! MP2

                     if (ouhf) then
                      do icart = 1, 3
                       dE = 0.D0
                       do ip1 = if1, il1
                        do ip2 = if2, il2
                         do ip3 = if3, il3
                          do ip4 = if4, il4
                           den2 = quart *(
     $                       (d_ij(ip1,ip2)  + d_ij2(ip1,ip2))
     $                       *(d_kl3(ip3,ip4) + d_kl4(ip3,ip4))
     $                       +(d_ij3(ip1,ip2) + d_ij4(ip1,ip2))
     $                       *(d_kl(ip3,ip4) + d_kl2(ip3,ip4)))
     $                       - eigh *(d_il(ip1,ip4) * d_jk3(ip2,ip3) 
     $                       + d_jk(ip2,ip3) * d_il3(ip1,ip4)
     $                       + d_jl(ip2,ip4) * d_ik3(ip1,ip3)
     $                       + d_ik(ip1,ip3) * d_jl3(ip2,ip4)
     $                       + d_il2(ip1,ip4) * d_jk4(ip2,ip3) 
     $                       + d_jk2(ip2,ip3) * d_il4(ip1,ip4)
     $                       + d_jl2(ip2,ip4) * d_ik4(ip1,ip3)
     $                       + d_ik2(ip1,ip3) * d_jl4(ip2,ip4))
                           dE = dE + den2 * buf(icount)
                           icount = icount + 1
                          end do
                         end do
                        end do
                       end do
                       frc(icart,idatom(i)) = frc(icart,idatom(i)) +
     $                   dE * q4 * scale
                      end do

                     else if (orohf) then
                      call errquit('only OUMP2 grads available',110)
                     else       ! RMP2
                      do icart = 1, 3
                       dE = 0.D0
                       do ip1 = if1, il1
                        do ip2 = if2, il2
                         do ip3 = if3, il3
                          do ip4 = if4, il4
                           den2 = quart*(d_ij(ip1,ip2)*d_kl2(ip3,ip4)
     $                       +  d_ij2(ip1,ip2) * d_kl(ip3,ip4))
     $                       - sixh *(d_il(ip1,ip4) * d_jk2(ip2,ip3) 
     $                       + d_jk(ip2,ip3) * d_il2(ip1,ip4)
     $                       + d_jl(ip2,ip4) * d_ik2(ip1,ip3)
     $                       + d_ik(ip1,ip3) * d_jl2(ip2,ip4))
                           dE = dE + den2 * buf(icount)
                           icount = icount + 1
                          end do
                         end do
                        end do
                       end do
                       frc(icart,idatom(i)) = frc(icart,idatom(i)) + 
     $                   dE * q4 * scale
                      end do

                     end if      ! closed shell
                     
                    end if       ! MP2

                   else         ! idatom
                    icount = icount + 3 * Nint
                   end if
                  end do         ! i (1:4 derivative center on integral)

C---------------------------------------------------------------------
                 end if          !schwarz(1234)
                end do           ! ish4
               end do            ! ish3
              end if             !schwarz(12)
             end do              ! ish2
            end do               ! ish1

            if (ochkpt) then

             call grad_store(rtdb, geom, g_force, chktime, nat, 
     $         ijatom+ndone, frc, .false., ofinish, oactive,
     $         oskel, omp2)
             if (ofinish) goto 100
            end if
            
            next = nxtask(nproc,task_size)
           end if                ! next

 1020      continue             ! Jump here for inactive/zero
          end do                 ! iat4
         end if
        end do                   ! iat2
       end do                    ! iat1
      end do                     ! iat3

C     final assembly of the forces
      ofinish = .true.
      call grad_store(rtdb, geom, g_force, chktime, nat, 
     $  ijatom+ndone+1, frc, .false., ofinish, oactive, oskel, omp2)
 100  continue

      next = nxtask(-nproc,task_size)
      return
      end


      subroutine  grad_mcscf_a
     $     (nbf,nopen,maxsh,plo,phi,coeff,pdm2,pdm2a)
      implicit none
      integer nbf, nopen, maxsh
      integer plo, phi
      double precision coeff(nbf,nopen)
      double precision pdm2(nopen,nopen,nopen,nopen)
      double precision pdm2a(nopen,nopen,nopen,maxsh)
c     
c     pdm2a(l,k,j,p) = sum(i) C(p,i) pdm2(l,k,j,i)
c     
c     for p in given shell(plo<=p<=phi)
c     
      integer i, j, k, l        ! Active MOs
      integer p                 ! AO
      double precision cpi
c     
      do p = 1,(phi-plo+1)
         do j = 1, nopen
            do k = 1, nopen
               do l = 1, nopen
                  pdm2a(l,k,j,p) = 0.0d0
               end do
            end do
         end do
      end do
c     
      do i = 1, nopen
         do p =  1,(phi-plo+1)
            cpi = coeff(p+plo-1,i)
            do j = 1, nopen
               do k = 1, nopen
                  do l = 1, nopen
                     pdm2a(l,k,j,p) = pdm2a(l,k,j,p) + 
     $                    pdm2(l,k,j,i)*cpi
                  end do
               end do
            end do
         end do
      end do
c     
      end

      subroutine  grad_mcscf_b
     $     (nbf,nopen,maxsh,pdim,qlo,qhi,coeff,pdm2a,pdm2b)
      implicit none
      integer nbf, nopen, maxsh
      integer pdim,qlo, qhi
      double precision coeff(nbf,nopen)
      double precision pdm2a(nopen,nopen,nopen,maxsh)
      double precision pdm2b(nopen,nopen,maxsh,maxsh)
c     
c     pdm2b(l,k,q,p) = sum(j) C(q,j) pdm2a(l,k,j,p)
c     
c     for q in given shell(qlo<=q<=qhi)
c     
      integer j, k, l           ! Active MOs
      integer p, q              ! AO
      double precision cqj
c     
      do p = 1, pdim
         do q = 1,(qhi-qlo+1)
            do k = 1, nopen
               do l = 1, nopen
                  pdm2b(l,k,q,p) = 0.0d0
               end do
            end do
         end do
      end do
c     
      do p =  1, pdim
         do q = 1,(qhi-qlo+1)
            do j = 1, nopen
               cqj = coeff(q+qlo-1,j)
               do k = 1, nopen
                  do l = 1, nopen
                     pdm2b(l,k,q,p) = pdm2b(l,k,q,p) + 
     $                    pdm2a(l,k,j,p)*cqj
                  end do
               end do
            end do
         end do
      end do
c     
      end

      subroutine  grad_mcscf_c
     $     (nbf,nopen,maxsh,pdim,qdim,rlo,rhi,coeff,pdm2b,pdm2c)
      implicit none
      integer nbf, nopen, maxsh
      integer pdim,qdim,rlo,rhi
      double precision coeff(nbf,nopen)
      double precision pdm2b(nopen,nopen,maxsh,maxsh)
      double precision pdm2c(nopen,maxsh,maxsh,maxsh)
c     
c     pdm2c(l,r,q,p) = sum(k) C(r,k) pdm2b(l,k,q,p)
c     
c     for r in given shell(rlo<=r<=rhi)
c     
      integer k, l              ! Active MOs
      integer p, q, r           ! AO
      double precision crk
c     
      do p = 1, pdim
         do q = 1, qdim
            do r = 1,(rhi-rlo+1)
               do l = 1, nopen
                  pdm2c(l,r,q,p) = 0.0d0
               end do
            end do
         end do
      end do
c     
      do p =  1, pdim
         do q = 1, qdim
            do r = 1,(rhi-rlo+1)
               do k = 1, nopen
                  crk = coeff(r+rlo-1,k)
                  do l = 1, nopen
                     pdm2c(l,r,q,p) = pdm2c(l,r,q,p) + 
     $                    pdm2b(l,k,q,p)*crk
                  end do
               end do
            end do
         end do
      end do
c     
      end


      subroutine  grad_mcscf_d
     $     (nbf,nopen,maxsh,pdim,qdim,rdim,slo,shi,coeff,pdm2c,pdm2d)
      implicit none
      integer nbf, nopen, maxsh
      integer pdim,qdim,rdim,slo,shi
      double precision coeff(nbf,nopen)
      double precision pdm2c(nopen,maxsh,maxsh,maxsh)
      double precision pdm2d(maxsh,maxsh,maxsh,maxsh)
c     
c     pdm2d(s,r,q,p) = sum(l) C(s,l) pdm2c(l,r,q,p)
c     
c     for s in given shell(slo<=s<=shi)
c     
      integer l                 ! Active MOs
      integer p, q, r, s        ! AO
      integer sdim
      double precision pdm2clrqp
c     
      sdim = (shi-slo+1)
c     
      do p = 1, pdim
         do q = 1, qdim
            do r = 1, rdim
               do s = 1, sdim
                  pdm2d(s,r,q,p) = 0.0d0
               end do
            end do
         end do
      end do
c     
      do p =  1, pdim
         do q = 1, qdim
            do r = 1, rdim
               do l = 1, nopen
                  pdm2clrqp = pdm2c(l,r,q,p)
                  do s = 1, sdim
                     pdm2d(s,r,q,p) = pdm2d(s,r,q,p) + 
     $                    pdm2clrqp*coeff(s+slo-1,l)
                  end do
               end do
            end do
         end do
      end do
c     
      end

