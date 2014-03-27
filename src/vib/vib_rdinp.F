*
* $Id$
*

      SUBROUTINE vib_RDINP(
     &    HESS,HESST,COORD,ATMASS,ZAN,OLDLBLAT,IZMAT,NELS,
     &    PROJEC,ZEROPE,HESOUT,INTERN,
     &    rtdb,first_pass)
C
C This routine reads the user input and the hessian matrix from
C      the appropriate place.  The default is to use the "default"
C      masses and to read the analytic hessian from tape 10.
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
      INTEGER NELS
      LOGICAL SETMAS,PROJEC,ZEROPE,HESOUT,INTERN
      logical dum_log
      integer natom, nat3, nhess, nhesst
      logical first_pass
      logical o_new_masses
      DOUBLE PRECISION HESS(NAT3,NAT3)  ! Hessian matrix
      DOUBLE PRECISION HESST(NHESST)    ! Lower triangle of hessian matrix
      DOUBLE PRECISION COORD(3,NATOM)   ! Atom x,y,z coordinates
      DOUBLE PRECISION ATMASS(NATOM)    ! Atoms mass
      DOUBLE PRECISION ZAN(NATOM)       ! Nuclear charge (i.e. atomic number)
      INTEGER OLDLBLAT(2,NATOM)     ! Atom character labels
      INTEGER IZMAT(NELS)
      integer num_new_masses
      integer geom_index
      double precision new_mass
      character*32 mass_rtdb_id
      character*16 mass_tag
      character*16 tag
      character*22 element
      character*2 symbol
      double precision q
      double precision xyz(3)
      COMMON /cvib_HESS/ NATOM,NAT3,NHESS,NHESST   ! Hessian information
      integer numans
      double precision ams, wave
      COMMON /cvib_SETCON/ AMS(36),WAVE,NUMANS     ! setup parameters
c
      integer h_lblb, i_lblb, geom
      integer il_cnt, iii, ijunk, jjj
c
      integer rtdb
      HESOUT = .FALSE.            ! do not write Hessian to tape 10
      o_new_masses = .false.
*
      if (first_pass) WRITE(luout,*)' Vib: Default input used '
      PROJEC = .FALSE.
      ZEROPE = .FALSE.
      HESOUT = .FALSE.
      INTERN = .FALSE.
      if (rtdb_get(rtdb,'vib:project',MT_LOG,1,dum_log)) then
        projec = dum_log
        write(luout,*)' vib:project option set to',projec
      endif
      if (rtdb_get(rtdb,'vib:zero point energy',MT_LOG,1,dum_log)) then
        zerope = dum_log
        write(luout,*)' vib:zero point energy option set to',zerope
      endif
      if (zerope) projec = .true.
*
* over ride all input to do normal and projected analysis
*
      if (first_pass) then
        projec = .false.
        zerope = .false.
      else
        projec = .true.
        zerope = .true.
      endif
*.. hessian read in vib_vib hess and hesst set there check it here
      call vib_chkhess(hess,nat3,first_pass)
* set labels, zan and coords information and mass information
c
      if (.not.ma_push_get(MT_BYTE,(2*natom),' labels for vib',
     &    h_lblb,i_lblb))
     &    call errquit('vib_rdinp ma_get for labels failed',911, MA_ERR)
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('vib_rdinp: error creating geometry',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('vib_rdinp: error loading geometry',911, RTDB_ERR)
c
c see if ANY mass is overwritten by freq/vib input block
c
      if (rtdb_get
     &    (rtdb,'vib:remass:count',mt_int,1,num_new_masses)) then
        o_new_masses = .true.
      else
        num_new_masses = 0
        o_new_masses = .false.
      endif
      il_cnt = 0
      do iii = 1,natom
        if(.not.geom_cent_get(geom,iii,tag,xyz,q))
     &      call errquit('vib_rdinp: geom_cent_get failed',911,
     &       GEOM_ERR)
        zan(iii) = q
        coord(1,iii) = xyz(1)
        coord(2,iii) = xyz(2)
        coord(3,iii) = xyz(3)
        if(.not.geom_tag_to_element(tag,symbol,element,ijunk)) then
          if (symbol.ne.'bq')
     &        call errquit('vib_rdinp: tag2elem fail',911, GEOM_ERR)
        endif
        byte_mb(il_cnt+i_lblb) = symbol(1:1)
        il_cnt = il_cnt+1
        byte_mb(il_cnt+i_lblb) = symbol(2:2)
        il_cnt = il_cnt+1
      enddo


*
* this loop structure separated out to do one pass over data in rtdb.
* set the incore geometry with any mass information from the freq/vib
* input.
*
      if (o_new_masses) then
        do jjj = 1,num_new_masses
          write(mass_rtdb_id,11111)jjj          ! lexical index
          if (rtdb_get(rtdb,mass_rtdb_id,mt_int,1,geom_index)) then
            write(mass_rtdb_id,11113)jjj
            if (.not.rtdb_get(rtdb,mass_rtdb_id,mt_dbl,1,new_mass))
     &          call errquit
     &          ('vib_rdinp: rtdb get failed for mass',911, RTDB_ERR)
            if (.not.geom_mass_set(geom,geom_index,new_mass))
     &          call errquit
     &          ('vib_rdinp: geom_mass_set failed for mass',911,
     &       GEOM_ERR)
          else  ! check for tag input 
            write(mass_rtdb_id,11112)jjj
            if (rtdb_cget(rtdb,mass_rtdb_id,1,mass_tag)) then 
              write(mass_rtdb_id,11113)jjj
              if (.not.rtdb_get(rtdb,mass_rtdb_id,mt_dbl,1,new_mass))
     &            call errquit
     &            ('vib_rdinp: rtdb get failed for mass',911, RTDB_ERR)
* tag input found loop over all geometry atoms for a possible match
              do iii = 1,natom
                if(.not.geom_cent_get(geom,iii,tag,xyz,q))
     &              call errquit('vib_rdinp: geom_cent_get failed',911,
     &       GEOM_ERR)
                if (bas_do_tags_match(mass_tag,tag)) then
                  if (.not.geom_mass_set(geom,iii,new_mass))
     &                call errquit
     &                ('vib_rdinp: geom_mass_set failed for mass',911,
     &       GEOM_ERR)
                endif
              enddo
            else
              if(.not.rtdb_print(rtdb,.true.)) call errquit('911',911,
     &       RTDB_ERR)
              write(luout,*)' expected mass input for index ',jjj
              write(luout,*)' It was neither index or tag based '
              call errquit
     &            ('vib_rdinp: fatal error for mass informaition',911,
     &       RTDB_ERR)
            endif
          endif
        enddo
      endif
*
* best place to set mass is in geometry object
* modified in-core geometry object with new masses from
* vibrational/freq input now reset's masses.
*
      setmas = geom_masses_get(geom,natom,atmass)
      if (.not.setmas)
     &    call errquit
     &    ('vib_rdinp: geom_masses_get failed',911, GEOM_ERR)
*      
      if (.not.geom_destroy(geom)) call errquit
     &    ('vib_rdinp: geom_destroy failed',911, GEOM_ERR)
C
C Write out atom information read here
C
      if (first_pass) then
        WRITE(6,10000)
C
        il_cnt = 0
        DO 00300 III = 1,NATOM
          WRITE(6,10001)
     &        byte_mb(il_cnt+i_lblb),
     &        byte_mb(il_cnt+i_lblb+1),
     &        III,
     &        COORD(1,III),COORD(2,III),COORD(3,III),
     &        ATMASS(III)
          il_cnt = il_cnt + 2
00300   CONTINUE
        WRITE(6,10002)
      endif
      if (.not.ma_pop_stack(h_lblb)) call errquit
     &    ('vib_rdinp ma_pop failed',911, MA_ERR)
      RETURN ! leave routine
10000 FORMAT(//,1X,28('-'),' Atom information ',28('-'),/,
     +5X,'atom',4X,'#',8X,'X',14X,'Y',14X,'Z',12X,'mass',/,
     +1X,74('-'))
10001 FORMAT(4x,a1,a1,3x,I5,4(1PD15.7))
10002 FORMAT(1X,74('-')///)
11111 format('vib:remass:',i4,':lexi:')
11112 format('vib:remass:',i4,':tags:')
11113 format('vib:remass:',i4,':mass:')
C
      END


