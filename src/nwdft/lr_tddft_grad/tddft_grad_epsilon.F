c
c     Add (epsilon_a-epsilon_i) contribution to the product vector.
c     HvD 6/2006, NG 01/2013
c
      subroutine tddft_grad_epsilon(g_x,g_Ax,ipol,nroot,nbf_ao,
     +                              nocc,nfc,eps)
c
c     This subroutine is explicitly parallelized based on the 
c     data distribution of g_x.
c
c     This subroutine is based on tddft_epsilon from lr_tddft, but 
c     has been completely reimplemented to work 3D arrays and to
c     minimize communications.
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer g_x(2)           ! GA handle for trial vectors
      integer g_Ax(2)          ! GA handle for trial vectors
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nroot            ! Number of roots
      integer nbf_ao           ! Number of AO basis functions
      integer nocc(2)          ! Number of occupied orbitals
      integer nfc(2)           ! Number of frozen cores
      integer alo(3),ahi(3)    ! Distribution of g_x array
      integer ld(3)            ! Leading dimensions of buffer
c
      double precision eps(nbf_ao,2) ! Orbital energies
c
c     Local:
c
      integer iproc            ! Processor number
      integer ip               ! counter over spin components
      integer ir               ! counter over spin vectors/roots
      integer mini             ! the lowest local occupied orbital
      integer maxi             ! the highest local occupied orbital
      integer numi             ! the number of local occupied orbitals
      integer mina             ! the lowest local virtual orbital
      integer maxa             ! the highest local virtual orbital
      integer numa             ! the number of local virtual orbitals
      integer i                ! counter over occupied orbital
      integer a                ! counter over virtual orbital
      integer imo              ! current occupied orbital
      integer amo              ! current virtual orbital
      integer ihdl_x,iptr_x    ! handle and index of memory block
      logical doit
c
      character*32 pname
c
      pname="tddft_grad_epsilon: "
c
c --------------------------------------------------------
c Add (epsilon_a-epsilon_i)*trial vector diagonal elements
c --------------------------------------------------------
c
      iproc=ga_nodeid()
      call ga_sync()
      do ip=1,ipol
c
        call nga_distribution(g_x(ip),iproc,alo,ahi)
        doit = .not.(alo(1).ne.1.or.ahi(1).ne.nroot)
c
        if (doit) then
          mini = alo(2)
          maxi = ahi(2)
          numi = maxi-mini+1
          mina = alo(3)
          maxa = ahi(3)
          numa = maxa-mina+1
          if (.not.ma_push_get(mt_dbl,nroot*numi*numa,'block x',
     +                       ihdl_x,iptr_x)) 
     +    call errquit(pname//'failed to allocate block x',0,UERR)
          ld(1) = nroot
          ld(2) = numi
          ld(3) = numa
          call nga_get(g_x(ip),alo,ahi,dbl_mb(iptr_x),ld)
          do a = 0, numa-1
            amo = nocc(ip)+mina+a
            do i = 0, numi-1
              imo = nfc(ip)+mini+i
              do ir = 0, nroot-1
                 dbl_mb(iptr_x+ir+nroot*i+nroot*numi*a)
     +        = dbl_mb(iptr_x+ir+nroot*i+nroot*numi*a)
     +        * (eps(amo,ip)-eps(imo,ip))
              enddo
            enddo
          enddo
          call nga_acc(g_Ax(ip),alo,ahi,dbl_mb(iptr_x),ld,1.0d0)
          if (.not.ma_pop_stack(ihdl_x)) 
     +       call errquit(pname//'failed to deallocate block x',0,UERR)
        endif ! doit
c
      enddo ! ipol
c
      call ga_sync()
c
      return
      end
c $Id$
