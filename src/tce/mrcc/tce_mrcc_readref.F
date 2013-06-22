       subroutine tce_mrcc_readref(rtdb)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "global.fh"
#include "tce_main.fh"

      integer rtdb
      integer l_occ0,k_occ0
      integer i,s,j
      character*8 t1   ! Name of variable
      character*3 t2   ! No. of ref.
      integer izero
      character*1025 ref
      integer iocc0(maxorb,2)
      character*1 o
      character*10 scftype
      integer n,n0,iu
      integer nactive(2)
      integer numactive(maxorb,2)
      integer iocc1
c      integer iphnum(maxorb,maxref,2)
c      integer invhnum(maxorb,maxref,2)
c      integer invpnum(maxorb,maxref,2)
      integer jr,ir
      integer ih,ip
c      logical mrccdebug
      logical nodezero
      integer isign,erank2
c      logical lusescffv
      character*1024 aocc0
ckbn  moved from include/tce_mrcc.fh
      integer iheffefrom(maxexcit,maxref,maxref)
      integer iheffeto(maxexcit,maxref,maxref)
      integer iheffespin(maxexcit,maxref,maxref)

      nodezero = (ga_nodeid().eq.0)

      if(nodezero)
     1write(LuOut,"(/,'MRCC module - reference analysis',/)")

c       if (.not.rtdb_get(rtdb,'mrcc:debug',mt_log,1,mrccdebug))
c     1 mrccdebug = .false.
c       if (.not.rtdb_get(rtdb,'mrcc:usescffermiv',mt_log,1,lusescffv))
c     1 lusescffv = .false.

      nodezero = (ga_nodeid().eq.0)
c
c ------------
c Input debug
c ------------
c
      if(mrccdebug.and.nodezero) then
      write(LuOut,"('Restricted: ',I1)")ipol-1
      if(ipol.eq.1) then
        nmo(2) = nmo(1)
        nocc(2) = nocc(1)
        write(LuOut,"('Orbitals, occ, virt: ',I5,I5,I5)")2*nmo(1),
     &2*nocc(1),2*(nmo(1)-nocc(1))
      else
        write(LuOut,"('nmo(1), occ(1), virt(1): ',I5,I5,I5)")nmo(1),
     &nocc(1),nmo(1)-nocc(1)
        write(LuOut,"('nmo(2), occ(2), virt(2): ',I5,I5,I5)")nmo(2),
     &nocc(2),nmo(2)-nocc(2)
      end if
      end if
c
c ----------------------
c Allocate local arrays
c ----------------------
c
c         if (.not.ma_push_get(mt_int,nmo(1)+nmo(2),'Occupancy of FV',
c     1    l_occ0,k_occ0)) call errquit('tce_mrcc_input: MA problem',
c     2    0,MA_ERR)
      do i=1, maxorb
           do j=1,maxref
              do s=1, 2
                 iocc(i,j,s) = 0
              end do
           end do
           ref(i:i) = '0'
           iocc0(i,1) = 0
           iocc0(i,2) = 0
        end do
c
c -------------
c Fermi vacuum 
c -------------
c
      do s=1,2
        do i=1,nmo(s)
          if(i.le.nocc(s)) then
            iocc0(i,s) = 1
          else
            iocc0(i,s) = 0
          end if
        end do
      end do
c
c -----------------------
c Fermi vacuum settings
c -----------------------
c
      if(lusescffv) then
        if(nodezero)
     1write(LuOut,"('Fermi vacuum according to SCF wave function 
     1 has been requested.')")
        do i=1,max(nmo(1),nmo(2))
          if((iocc0(i,1).eq.1).and.(iocc0(i,2).eq.1)) then
            aocc0(i:i)='2'
          elseif((iocc0(i,1).eq.1)) then
            aocc0(i:i)='a'
          elseif((iocc0(i,2).eq.1)) then
            aocc0(i:i)='b'
          else
            aocc0(i:i)='0'
          end if
        end do
      endif
c        if(nodezero)write(LuOut,*)aocc0(1:max(nmo(1),nmo(2)))
c
c -------------
c   Read refs
c -------------
c

      if (.not.rtdb_get(rtdb,'bwcc:nref', mt_int, 1, nref))
     1  call ERRQUIT('tce_mrcc_readref',1,RTDB_ERR)

c      if (.not.rtdb_cget(rtdb, 'scf:scftype',1,scftype))
c     1  call ERRQUIT('tce_mrcc_readref',2,RTDB_ERR)

      if(nref.lt.1) call errquit('tce_mrcc_readref',3,RTDB_ERR)

      t1 = 'bwcc:ref'
      izero = ICHAR('0')

      do j=1,nref
         write(t2,"(I3.3)")j
         if (.not.rtdb_cget(rtdb,t1//t2,1,ref))
     1     call errquit('tce_mrcc_readref',0,RTDB_ERR)
        do i=1,maxorb-1
         o = ref(i:i)
         if(o.eq.'2') then
           iocc(i,j,1)=1
           iocc(i,j,2)=1
         else
            if (o.eq.'a' .or. o.eq.'b') then
c               if (scftype.eq.'RHF') then
c                 write(LuOut,*)'Open shell references not allowed for',
c     &                       ' RHF based calculation'
c                 call errquit('tce_mrcc_readref',1,UNKNOWN_ERR)
c               else
                  if (o.eq.'a') then
                     iocc(i,j,1) = 1
                  else
                     iocc(i,j,2) = 1
                  end if
c               end if
            else
               if (o.ne.'0'.and.o.ne.' ') then
                  write(LuOut,*) 'Illegal character in configuration ',
     &                       'specification: (',ICHAR(o),')'
                  call  errquit('tce_mrcc_readref',2,UNKNOWN_ERR)
               end if
            end if
         end if
        end do
      end do

      if(mrccdebug) then
        write(LuOut,"(/,'Debug output for iocc(orb,ref,spin) array',/)")
        do j=1,nref 
         write(LuOut,"('Reference ',I4)")j
         do iu=1,2
          write(LuOut,"(200I1)")
     &    (iocc(i,j,iu),i=1,nmo(iu))
         end do
        end do
      end if

c
c ----------------------
c Check number of elecs
c ----------------------
c
      do s = 1,2
         do j = 1, nref
            n  = 0
            n0 = 0
            do i = 1, nmo(s)
               n  = n  + iocc(i,j,s)
               n0 = n0 + iocc0(i,s)
            end do
            if (n.ne.nocc(s) .or. n0.ne.nocc(s)) then
               write(LuOut,*)
     &            'Inconsistent number of electrons in bwcc input'
     1            ,' n:',n,' ne.nocc(s):',nocc(s),' n0:',n0
               call errquit('tce_mrcc_readref',4,UNKNOWN_ERR)
            end if
         end do
      end do
c
c ---------------------
c Active orbitals 
c ---------------------
c
      do iu = 1, ipol
         do i = 1, nmo(iu)
            isactive(i,iu) = .false.
c CH  
        if(lusescffv) then
            iocc1 = iocc0(i,iu)
        else
            iocc1 = iocc(i,1,iu)
        endif
c CH E
        if(.not.lusescffv) then
            do j = 2, nref
               if (iocc(i,j,iu).ne.iocc1) then
                  isactive(i,iu) = .true.
                  goto 100
               end if
            end do
        else
            do j = 1, nref
               if (iocc(i,j,iu).ne.iocc1) then
                  isactive(i,iu) = .true.
                  goto 100
               end if
            end do
        endif
 100        continue
            if (ipol.eq.1) isactive(i,2) = isactive(i,1)
         end do
c prepare nactive, numactive - count and indices of active orbitals
         nactive(iu) = 0
         do i = 1, nmo(iu)
            if (isactive(i,iu)) then
               nactive(iu) = nactive(iu)+1
               numactive(nactive(iu),iu) = i
            end if
         end do
      end do
 
      if(mrccdebug) then
        write(LuOut,"(/,'Debug output for isactive(orb,spin) array',/)")
        do iu=1,2
         write(LuOut,"(200L1)")
     & (isactive(i,iu),i=1,nmo(iu))
        end do
      end if
c
c -----------------------------------------
c Prepare renumbering vectors for later use
c -----------------------------------------
c
      do iu = 1, ipol
         do j = 1, nref
            ip = 0
            ih = 0
            do i = 1, nmo(iu)
               if (iocc(i,j,iu).gt.0) then
                  ih = ih+1
                  iphnum(i,j,iu) = ih
                  invhnum(ih,j,iu) = i
               else
                  ip = ip+1
                  iphnum(i,j,iu) = ip
                  invpnum(ip,j,iu) = i
               end if
               if (ipol.eq.1) then
                  iphnum(i,j,2)   = iphnum(i,j,1)
                  invhnum(ih,j,2) = invhnum(ih,j,1)
                  invpnum(ip,j,2) = invpnum(ip,j,1)
               end if
            end do
         end do
      end do
  
      if(lusescffv) then

c      write(6,*)'IPHNUM...'
      do iu = 1, ipol
            ip = 0
            ih = 0
            do i = 1, nmo(iu)
               if (iocc0(i,iu).gt.0) then
                  ih = ih+1
                  iphnum(i,maxref,iu) = ih
                  invhnum(ih,maxref,iu) = i
c                  write(6,"('IH: ',I4,I4)"),i,ih
               else
                  ip = ip+1
                  iphnum(i,maxref,iu) = ip
                  invpnum(ip,maxref,iu) = i
c                  write(6,"('IP: ',I4,I4)"),i,ih
               end if
               if (ipol.eq.1) then
                  iphnum(i,maxref,2)   = iphnum(i,maxref,1)
                  invhnum(ih,maxref,2) = invhnum(ih,maxref,1)
                  invpnum(ip,maxref,2) = invpnum(ip,maxref,1)
               end if
            end do
      end do

      endif
c
c --------------------------------------------------------
c Excitations between Fermi vacuum and reference configs
c --------------------------------------------------------
c
       if(lusescffv) then

         do jr = 1, nref
           call tce_mrcc_calcexcit(ihefferank0(1,jr),
     &                  iheffefrom(1,jr,maxref),
     &                  iheffeto(1,jr,maxref),iheffespin(1,jr,maxref),
     &                  jr,0,nmo,rtdb,iocc0)

         enddo

      do ir=1,nref
         erank2=ihefferank0(1,ir)+ihefferank0(1,ir)
         do i=1, ihefferank0(1,ir)
            isign=3-2*iheffespin(i,ir,maxref)
            mutualoper0(i,ir,1)=isign*
     &         invpnum(iheffeto(i,ir,maxref),maxref,
     2 iheffespin(i,ir,maxref))
            mutualoper0(erank2-i+1,ir,1)=isign*
     &         invhnum(iheffefrom(i,ir,maxref),maxref,
     1 iheffespin(i,ir,maxref))
         enddo
      enddo

      endif
c 
c ---------------------------------------------
c Excitations between reference configurations
c ---------------------------------------------
c
      do ir = 1, nref
         do jr = 1, nref
            call tce_mrcc_calcexcit(ihefferank(jr,ir),
     &                     iheffefrom(1,jr,ir),
     &                     iheffeto(1,jr,ir),iheffespin(1,jr,ir),
     &                     jr,ir,nmo,rtdb,iocc0)
c iphnum deleted
            if (ir.ne.jr .and. ihefferank(jr,ir).eq.0 ) then
                  write(LuOut,"('References no.',I3,I3,
     &                       ' are identical!')")ir,jr
            call errquit('tce_mrcc_input',12,UNKNOWN_ERR)
            end if
         end do
      end do

c       if(nodezero) then
c         write(LuOut,"(/,'Excitations between references',/)")
c         do ir = 1, nref
c            write(LuOut,'(32I2)') (ihefferank(ir,jr),jr=1,nref)
c         end do
c        end if
c
c ---------------------------------------------
c Operators 
c ---------------------------------------------
c
      do ir=1,nref
      do jr=1,nref
      if(ir.ne.jr) then
         erank2=ihefferank(ir,jr)+ihefferank(ir,jr)
         do i=1, ihefferank(ir,jr)
            isign=3-2*iheffespin(i,ir,jr)
            mutualoper(i,ir,jr)=isign*
     &         invpnum(iheffeto(i,ir,jr),jr,iheffespin(i,ir,jr))
            mutualoper(erank2-i+1,ir,jr)=isign*
     &         invhnum(iheffefrom(i,ir,jr),jr,iheffespin(i,ir,jr))
         enddo
      endif
      enddo
      enddo
c
c ---------------------------------------------
c Operators debug output
c ---------------------------------------------
c
c        if(nodezero) then
c            write(LuOut,"(/,'Debug output for mutualoper array',/)")
c            do ir = 1, nref
c            do jr = 1, nref
c               write(LuOut,"('Heff (',I3,',',I3,'): erank: ',I3)")
c     &                    ir,jr,ihefferank(ir,jr)
c               write(LuOut,"(200I4)")
c     & (mutualoper(i,ir,jr),i=1,2*ihefferank(ir,
c     1 jr))
c            end do
c            end do
c         end if
c
c
c ------------------------
c Deallocate local arrays 
c ------------------------
c
c       if (.not.ma_pop_stack(l_occ0))
c     1    call errquit("tce_mrcc_readref: MA problem",1,MA_ERR)

      if(nodezero)
     1write(LuOut,"(/,'MRCC module finished',/)")

      return
      end

c $Id$
