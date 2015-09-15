cc AJL/Begin/FDE
C> \brief Compute 1-electron integrals and add them to a global array
C>
C> This routine computes 1-electron integrals of the specified kind and
C> adds them to a global array. The routine always computes a
C> rectangular matrix. It is even capable of using two different 
C> basis sets. Different kinds of integrals can be computed, they are:
C>
C> * "fde_pot": the potential energy integrals for the FDE nuclei
C>
C> The routine can exploit symmetry. If oskel is true then the 
C> "petite-list" symmetry will be used, although this requires both
C> basis sets to be the same.
C
      subroutine int_1e_fde(ibas, jbas, g, integ_type, oskel, fde_basis)
      implicit none
#include "errquit.fh"
c#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
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
      integer ibas             !< [Input] bra basis sets
      integer jbas             !< [Input] ket basis sets
      integer g(*)             !< [Output] GA handle to array 
      character*(*) integ_type !< [Input] Name of integrals to compute
      logical oskel            !< [Input] If true generate symmetry unique list
      integer fde_basis        !< [Input] FDE basis set handle
c
c     local variables
c
      integer type
c     
c     local variables (ooldga)
c     
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      logical odoit
      double precision q2
      integer i, noffset,g_loc
      logical odbug, osome
cc AJL/Debug
c      write(6,*) 'int_1e_fde : ',integ_type

      call ga_sync()
c
      if (ibas.ne.jbas) call errquit
     $        ('int_1e_fde: use of symmetry requires ibas=jbas', ibas,
     &       BASIS_ERR)
c     
      if (inp_compare(.false., integ_type, 'fde_pot')) then
         type = 1
      else
         write(6,*) ' integ_type = ', integ_type
         call errquit('int_1e_fde: unknown integ_type', 0, INT_ERR)
      end if

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
c     
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for ibas', ibas,
     &       BASIS_ERR)
      if (.not. bas_numcont(jbas, nshell_j)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for jbas', jbas,
     &       BASIS_ERR)

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
cso

      if (.not. MA_push_get(MT_DBL,max1e,'int_1e_fde:buf',l_buf,k_buf))
     $     call errquit('int_1e_fde: ma failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'int_1e_fde:scr',l_scr,k_scr))
     $     call errquit('int_1e_fde: ma failed', mem1, MA_ERR)
c     
c     Loop thru shells with static parallel work decomposition
c
c      call ga_print(g)

      if(.not.ga_duplicate(g,g_loc,'local g')) call
     .     errquit('int1e: dupl failed',0, GA_ERR)

      call ga_zero(g_loc)

      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      q2 = 1.0d0
      do jshell = 1, nshell_j
         do ishell = 1, nshell_i
c
            if (mod(ijshell, nproc) .eq. iproc) then
               odoit = .true.

               if (oskel) 
     $              odoit = sym_shell_pair(ibas, ishell, jshell, q2)

               if (odoit) then
                  if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) 
     $                 call errquit('int_1e_fde: bas_cn2bfr ?', ibas,
     &                 BASIS_ERR)
                  if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $                 call errquit('int_1e_fde: bas_cn2bfr ?', jbas,
     &                 BASIS_ERR)
                  idim = ihi - ilo + 1
c     
c     Generate the integrals
c     

                  if (type .eq. 1) then
                     call int_1efde1(jbas, jshell, ibas, ishell, 
     $                    fde_basis, mem1, dbl_mb(k_scr), max1e, 
     $                    dbl_mb(k_buf))
                  else
                     call errquit('int_1e_fde: invalid type?', type,
     &               GA_ERR)
                  end if
c     
c     Add the integrals into the global array
c    
                  if(odbug) then
                     write(6,*) 'ga_acc = ...',ga_nodeid()
                     call util_flush(6)
                  endif
cedo                     call ga_acc(g, ilo, ihi, jlo, jhi, dbl_mb(k_buf),
cedo     $                    idim, q2) 
                  if(q2.ne.1) then
                     call dscal((ihi-ilo+1)*(jhi-jlo+1),
     ,               q2,dbl_mb(k_buf),1)
                  endif
                  call ga_put(g_loc, ilo, ihi, jlo, jhi,
     ,            dbl_mb(k_buf),  idim) 
               end if
            endif
            ijshell = ijshell + 1
         end do
      end do
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_1e_fde: pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_1e_fde: pop failed', 0, GA_ERR)
c
c     
      call ga_dadd(1.0d0, g_loc, 1.0d0, g, g)
      if (.not. ga_destroy(g_loc)) call errquit('i1ega: gad?',0,
     &       GA_ERR)

      call ga_sync()            ! So that no nasty races can result

      end
cc AJL/End
