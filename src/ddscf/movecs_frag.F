      subroutine movecs_fragment(rtdb, geom, basis, 
     $     module, iset, nclosed, nopen, g_vecs, evals)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb              ! [input]
      integer geom, basis       ! [input]
      character *(*) module     ! [input] Name of invoking module
      integer iset              ! [input] 1 ROHF, 1/2=alpha/beta UHF
      integer nclosed           ! [input] doubly occ in each set (zero UHF)
      integer nopen             ! [input] single (active) occ in each set
      integer g_vecs            ! [output] GA containing MOs
      double precision evals(*) ! [output] Eigenvalues from the fragments
c
c     Assemble initial guess MOs from the MO vectors of
c     fragments.
c
      integer nmo, nbf, nvirt, closed_off, open_off, virt_off, bf_off
      integer g_frag, l_eval, k_eval, l_occ, k_occ, excess
      integer nclosed_frag, nopen_frag, nvirt_frag, ioff
      integer i
      integer maxfrag
      parameter (maxfrag = 100)
      character*256 fragments(maxfrag), key
      integer ma_type, nfrag
c
      logical oprint
c
      character*255 title_frag, basis_name_frag
      character*20 scftype_frag
      integer nbf_frag, nsets_frag, nmo_frag(2)
      integer frag
      logical movecs_read_header, movecs_read
      external  movecs_read_header, movecs_read
c
      call ga_sync()
      oprint = util_print('mo guess', print_default) 
     $     .and. (ga_nodeid().eq.0)
      if (oprint) then
         write(6,*)
         call util_print_centered(6,
     $        'Fragment Molecular Orbital Initial Guess', 22, .true.)
         write(6,*)
         call util_flush(6)
      endif
      key = ' '
      write(key,'(a,a)') module(1:inp_strlen(module)),':nfrag'
      if (.not. rtdb_get(rtdb, key, mt_int, 1, nfrag))
     $     call errquit('movecs_fragment: no fragment info?',0,
     &       RTDB_ERR)
      if (nfrag .gt. maxfrag) 
     $     call errquit('movecs_fragment: too many fragments', nfrag,
     &       UNKNOWN_ERR)
      write(key,'(a,a)') module(1:inp_strlen(module)),':fragments'
      if (.not. rtdb_cget(rtdb, key, nfrag, fragments))
     $     call errquit('movecs_fragment: can''t read fragments', nfrag,
     &       RTDB_ERR)
c
      call ga_inquire(g_vecs, ma_type, nbf, nmo)
      nvirt       = nmo - nclosed - nopen
      closed_off  = 1
      open_off    = nclosed + 1
      virt_off    = nclosed + nopen + 1
      bf_off      = 1
c
      call ga_zero(g_vecs)
c
      do frag = 1, nfrag
         call util_file_name_resolve(fragments(frag),.false.)
         if (.not. movecs_read_header(fragments(frag), 
     $        title_frag, basis_name_frag, scftype_frag, nbf_frag, 
     $        nsets_frag, nmo_frag,2)) call errquit
     $        ('movecs_fragment: failed reading frag header', frag,
     &       INPUT_ERR)
         if (.not. ga_create(MT_DBL, nbf_frag, nmo_frag(1), 
     $        ' frag mos',  32, 32, g_frag)) 
     $        call errquit('movecs_fragment: GA for frag MOs', nmo_frag,
     &       GA_ERR)
         if (.not. ma_push_get(mt_dbl, nbf_frag, 
     $        ' fragment eva', l_eval, k_eval)) call errquit
     $        ('movecs fragment: insufficient memory?', nbf_frag,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, nbf_frag, 
     $        ' fragment eva', l_occ, k_occ)) call errquit
     $        ('movecs fragment: insufficient memory?', nbf_frag,
     &       MA_ERR)
c     
         if (.not. movecs_read(fragments(frag), iset, dbl_mb(k_occ), 
     $        dbl_mb(k_eval), g_frag)) call errquit
     $        ('movecs fragment: failed read fragment MOs',0,
     &       INPUT_ERR)
c
**         call ga_print(g_frag)
c
         nclosed_frag = 0
         nopen_frag   = 0
         nvirt_frag   = 0
         do i = 0, nmo_frag(1)-1
            if (dbl_mb(k_occ+i).eq.2.0d0) then
               nclosed_frag = nclosed_frag + 1
            else if (dbl_mb(k_occ+i).ne.0.0d0) then
               nopen_frag = nopen_frag + 1
            else
               nvirt_frag = nvirt_frag + 1
            endif
         enddo
         if (oprint) then
            write(6,3) title_frag(1:inp_strlen(title_frag)),
     $           basis_name_frag(1:inp_strlen(basis_name_frag)),
     $           scftype_frag, nbf_frag, nsets_frag, nmo_frag(1),
     $           nclosed_frag, nopen_frag, nvirt_frag
 3          format(
     $           ' title    ', a/
     $           ' basis    ', a/
     $           ' scftype  ', a/
     $           ' nbf      ', i6/
     $           ' nsets    ', i6/
     $           ' nmo      ', i6/
     $           ' nclosed  ', i6/
     $           ' nopen    ', i6/
     $           ' nvirtual ', i6/)
            call util_flush(6)
         endif
c
         if (bf_off+nbf_frag-1 .gt. nbf) call errquit
     $        ('movecs_fragment: no. basis function mismatch', nbf,
     &       INPUT_ERR)
         ioff = 1
         if (nclosed_frag .gt. 0) then
            if (closed_off+nclosed_frag-1 .gt. nclosed) call errquit
     $           ('movecs_fragment: closed shell mismatch', nclosed,
     &       INPUT_ERR)
            call ga_copy_patch('n',
     $           g_frag, 1, nbf_frag, ioff, ioff+nclosed_frag-1,
     $           g_vecs, bf_off, bf_off+nbf_frag-1,
     $                   closed_off,closed_off+nclosed_frag-1)
            call dcopy(nclosed_frag,dbl_mb(k_eval+ioff-1),1,
     $           evals(closed_off), 1)
            closed_off = closed_off + nclosed_frag
            ioff = ioff + nclosed_frag
         endif
         if (nopen_frag .gt. 0) then
            if (open_off+nopen_frag-1 .gt. nclosed+nopen) call errquit
     $           ('movecs_fragment: open shell mismatch', nopen,
     &       INPUT_ERR)
            call ga_copy_patch('n',
     $           g_frag, 1, nbf_frag, ioff, ioff+nopen_frag-1,
     $           g_vecs, bf_off, bf_off+nbf_frag-1,
     $                   open_off,open_off+nopen_frag-1)
            call dcopy(nopen_frag,dbl_mb(k_eval+ioff-1),1,
     $           evals(open_off), 1)
            open_off = open_off + nopen_frag
            ioff = ioff + nopen_frag
         endif
c
         if (nvirt_frag .gt. 0) then
c            if (virt_off+nvirt_frag-1 .gt. nmo) call errquit
c     $           ('movecs_fragment: virt shell mismatch', nvirt)
            if (virt_off+nvirt_frag-1 .gt. nmo) then
               excess = virt_off+nvirt_frag-1 - nmo
            if (excess.LT.0) call errquit
     $         ('movecs_fragment: virt shell mismatch',excess,
     &       INPUT_ERR)
             if(oprint) then
               write(6,*) '**************WARNING*****************'
               write(6,*) ' movecs_fragment: virt shell mismatch,'
               write(6,*) ' probably due to linear dependence.'
               write(6,*) ' virtual orbitals excluded: ', excess
               write(6,*) '**************************************'
             endif
            else 
               excess=0
            endif
            call ga_copy_patch('n',
     $           g_frag, 1, nbf_frag, ioff, ioff+nvirt_frag-1-excess,
     $           g_vecs, bf_off, bf_off+nbf_frag-1,
     $                   virt_off,virt_off+nvirt_frag-1-excess)
            call dcopy(nvirt_frag-excess,dbl_mb(k_eval+ioff-1),1,
     $           evals(virt_off), 1)
            virt_off = virt_off + nvirt_frag - excess
            ioff = ioff + nvirt_frag
         endif
c
         bf_off = bf_off + nbf_frag
c
**         call ga_print(g_vecs)
c
         if (.not. ma_pop_stack(l_occ)) 
     $        call errquit('movecs_fragment: stack is bad',0, MA_ERR)
         if (.not. ma_pop_stack(l_eval)) 
     $        call errquit('movecs_fragment: stack is bad2',0, MA_ERR)
         if (.not. ga_destroy(g_frag))
     $        call errquit('movecs_fragment: ga is bad?',0, GA_ERR)
c
      enddo
c
*      call output(evals,1,nmo,1,1,nmo,1,1)
c
      if (closed_off .ne. nclosed+1) call errquit
     $     ('movecs_fragment: bad # closed shells',  closed_off-1,
     &       INPUT_ERR)
      if (open_off .ne. nclosed+nopen+1) call errquit
     $     ('movecs_fragment: bad # open shells',  open_off-closed_off,
     &       INPUT_ERR)
c      if (virt_off .ne. nmo+1) call errquit
c     $     ('movecs_fragment: bad # virtuals', virt_off-open_off)
c
      call ga_sync()
c
      end


