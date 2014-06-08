      subroutine tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +                                  nroot,facmo,facao,trans,
     +                                  g_mo,g_mmo,block,g_mao)
c
c     $Id$
c
c     Transform a quantity from MO basis to AO basis. The quantities
c     (both MO and AO basis) are assumed to be stored in 3 dimensional
c     global arrays where the first index runs as follows:
c
c     - MO quantities:
c
c       ip=1,pol, ir=1,nroot:
c
c          g_mmo(ip) index1=ir
c
c     - AO quantities
c
c          g_mao     index1=(ip-1)*nroot+ir
c
c     The second and third indeces run over basis functions.
c
c     Dirty extension: if g_mao is a 2 dimensional global array then
c     the meaning of nroot changes to mean that only the matrix for that
c     root is to be transformed to AO basis.
c
c
c     The factors facmo and facao used in the following way
c
c        g_mao = facao*g_mao + facmo*transform(g_mmo)
c
c     to allow control over whether or not g_mao should be initialised
c     or added to, and the scale factor of the contribution of the 
c     transformed quantity.
c
c     The parameter trans requests which part of g_mmo should be 
c     transformed. Valid values are:
c
c        "ij" transform the occupied-occupied block only
c        "ib" transform the occupied-virtual block only
c        "iq" transform the occupied-all block only
c        "aj" transform the virtual-occupied block only
c        "ab" transform the virtual-virtual block only
c        "aq" transform the virtual-all block only
c        "pj" transform the all-occupied block only
c        "pb" transform the all-virtual block only
c        "pq" transform the all-all block only
c
c     The parameter block refers to what is stored in g_mmo.
c     Otherwise the valid values are the same and their meaning is
c     analogous.
c
c     Written by Huub van Dam, May 2006.
c
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
c     Input:
c
      integer ipol     ! =1 (restricted), =2 (unrestricted)
      integer nao      ! the number of AOs
      integer nfc(2)   ! the number of frozen cores
      integer naoc(2)  ! the number of active occupied orbitals
      integer nocc(2)  ! the number of occupied orbitals
      integer nav(2)   ! the number of active virtual orbitals
      integer nfv(2)   ! the number of frozen virtuals
      integer nroot    ! the number of roots to be treated
      integer g_mo(2)  ! global arrays to the MOs
      integer g_mmo(2) ! global arrays to hold the MO basis quantities
c
      double precision facmo ! scale factor for transformed MO terms
      double precision facao ! scale factor for AO terms
c
      character*(*) trans    ! which block should be transformed
      character*(*) block    ! which block in MO basis do we have
c
c     Output:
c
      integer g_mao    ! global arrays to hold the AO basis quantities
c
c     Local
c
      integer g_diu    ! work space array
      integer ndim     ! dimension of the work space array
      integer idim(3)  ! dimensions
      integer ichnk(3) ! chunking
      integer alo(3)   ! lower chunk limits on A
      integer ahi(3)   ! upper chunk limits on A
      integer blo(3)   ! lower chunk limits on B
      integer bhi(3)   ! upper chunk limits on B
      integer clo(3)   ! lower chunk limits on C
      integer chi(3)   ! upper chunk limits on C
      integer ip       ! counter for ipol
      integer ir       ! counter for nroot
      integer ndimao   ! the dimension of g_mao
      integer iroot    ! the first root to be transformed
c
c     What is the dimensionality of g_mao?
c
c     call ga_sync()
      ndimao = GA_Ndim(g_mao)
      if (ndimao.eq.2) then
        iroot = nroot
      else if (ndimao.eq.3) then
        iroot = 1
      else 
        call errquit
     +  ("tddft_grad_trans_mo2ao: invalid dimension for g_mao",0,UERR)
      endif
c
c     Check if combination of trans and block are valid
c
      if (trans.ne."ij".and.trans.ne."ib".and.trans.ne."iq".and.
     +    trans.ne."aj".and.trans.ne."ab".and.trans.ne."aq".and.
     +    trans.ne."pj".and.trans.ne."pb".and.trans.ne."pq") then
        call errquit("tddft_grad_trans_mo2ao: invalid value for trans",
     +               0, UERR)
      endif
      if (block.ne."ij".and.block.ne."ib".and.block.ne."iq".and.
     +    block.ne."aj".and.block.ne."ab".and.block.ne."aq".and.
     +    block.ne."pj".and.block.ne."pb".and.block.ne."pq") then
        call errquit("tddft_grad_trans_mo2ao: invalid value for block",
     +               0, UERR)
      endif
      if (trans(1:1).eq."i") then
        if (block(1:1).eq."a") call errquit
     +    ("tddft_grad_trans_mo2ao: invalid combination trans-block",
     +     0, UERR)
      else if (trans(1:1).eq."a") then
        if (block(1:1).eq."i") call errquit
     +    ("tddft_grad_trans_mo2ao: invalid combination trans-block",
     +     0, UERR)
      else
        if (block(1:1).ne."p") call errquit
     +    ("tddft_grad_trans_mo2ao: invalid combination trans-block",
     +     0, UERR)
      endif
      if (trans(2:2).eq."j") then
        if (block(2:2).eq."b") call errquit
     +    ("tddft_grad_trans_mo2ao: invalid combination trans-block",
     +     0, UERR)
      else if (trans(2:2).eq."b") then
        if (block(2:2).eq."j") call errquit
     +    ("tddft_grad_trans_mo2ao: invalid combination trans-block",
     +     0, UERR)
      else
        if (block(2:2).ne."q") call errquit
     +    ("tddft_grad_trans_mo2ao: invalid combination trans-block",
     +     0, UERR)
      endif
c
c     Find dimension of work space global array
c
      ndim = 0
      if (trans.eq."ij") then
        do ip = 1, ipol
          ndim = max(ndim,naoc(ip))
        enddo
      else if (trans.eq."ib".or.trans.eq."aj") then
        do ip = 1, ipol
          ndim = max(ndim,naoc(ip),nav(ip))
        enddo
      else if (trans.eq."ab") then
        do ip = 1, ipol
          ndim = max(ndim,nav(ip))
        enddo
      else
        do ip = 1, ipol
          ndim = max(ndim,naoc(ip)+nav(ip))
        enddo
      endif
c
c     Create global array to hold partially transformed quantities
c
      idim(1) = ndim
      idim(2) = nao
      ichnk(1) = -1
      ichnk(2) = -1
      if (.not.nga_create(mt_dbl,2,idim,'vectors Diu',ichnk,g_diu))
     +    call errquit('tddft_grad_trans_mo2ao: failed to create g_diu',
     +                 0, GA_ERR)
c
c     The actual transformation
c
      do ip = 1, ipol
c
c       Tuv = sum_pq Cup*Tpq*Cvq
c
        do ir = iroot, nroot
c
c         Dpv = sum_q Tpq*Cvq
c
          alo(1) = ir
          ahi(1) = ir
          blo(2) = 1
          bhi(2) = nao
          clo(2) = 1
          chi(2) = nao
          if (trans(1:1).eq.'i') then
            alo(2) = 1
            ahi(2) = naoc(ip)
            clo(1) = 1
            chi(1) = naoc(ip)
          else if (trans(1:1).eq.'a') then
            alo(2) = 1
            ahi(2) = nav(ip)
            if (block(1:1).eq.'p') then
              alo(2) = alo(2)+naoc(ip)
              ahi(2) = ahi(2)+naoc(ip)
            endif
            clo(1) = 1
            chi(1) = nav(ip)
          else
            alo(2) = 1
            ahi(2) = naoc(ip)+nav(ip)
            clo(1) = 1
            chi(1) = naoc(ip)+nav(ip)
          endif
          if (trans(2:2).eq.'j') then
            alo(3) = 1
            ahi(3) = naoc(ip)
            blo(1) = nfc(ip)+1
            bhi(1) = nocc(ip)
          else if (trans(2:2).eq.'b') then
            alo(3) = 1
            ahi(3) = nav(ip)
            if (block(2:2).eq.'q') then
              alo(3) = alo(3)+naoc(ip)
              ahi(3) = ahi(3)+naoc(ip)
            endif
            blo(1) = nocc(ip)+1
            bhi(1) = nocc(ip)+nav(ip)
          else
            alo(3) = 1
            ahi(3) = naoc(ip)+nav(ip)
            blo(1) = nfc(ip)+1
            bhi(1) = nocc(ip)+nav(ip)
          endif
          call nga_matmul_patch('n','t',1.0d0,0.0d0,g_mmo(ip),alo,ahi,
     +                          g_mo(ip),blo,bhi,g_diu,clo,chi)
c
c         Tuv = sum_p Cup*Dpv
c
          alo(1) = 1
          ahi(1) = nao
          blo(2) = 1
          bhi(2) = nao
          if (ndimao.eq.3) then
            clo(1) = ir+(ip-1)*nroot
            chi(1) = ir+(ip-1)*nroot
            clo(2) = 1
            chi(2) = nao
            clo(3) = 1
            chi(3) = nao
          else
            clo(1) = 1
            chi(1) = nao
            clo(2) = 1
            chi(2) = nao
          endif
          if (trans(1:1).eq.'i') then
            alo(2) = nfc(ip)+1
            ahi(2) = nocc(ip)
            blo(1) = 1
            bhi(1) = naoc(ip)
          else if (trans(1:1).eq.'a') then
            alo(2) = nocc(ip)+1
            ahi(2) = nocc(ip)+nav(ip)
            blo(1) = 1
            bhi(1) = nav(ip)
          else
            alo(2) = nfc(ip)+1
            ahi(2) = nocc(ip)+nav(ip)
            blo(1) = 1
            bhi(1) = naoc(ip)+nav(ip)
          endif
          call nga_matmul_patch('n','n',facmo,facao,g_mo(ip),alo,ahi,
     +                          g_diu,blo,bhi,g_mao,clo,chi)
        enddo
      enddo
c
c     call ga_sync()
      if (.not.ga_destroy(g_diu)) call errquit
     +   ('tddft_grad_trans_mo2ao: failed to destroy g_diu',0, GA_ERR)
c
      end
