      subroutine sd_6ds_1(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h2d,h1d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
       triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_2(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
      do h1=1,h1d
      do h3=1,h3d
       triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_3(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     2 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
      do h3=1,h3d
      do h1=1,h1d
       triplesx(h1,h3,h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h1,h3,h2,p6,p5-jstart+1,p4-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_4(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h2d,h1d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
       triplesx(h3,h2,h1,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h3,h2,h1,p6,p4-jstart+1,p5-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_5(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
      do h1=1,h1d
      do h3=1,h3d
       triplesx(h3,h1,h2,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h3,h1,h2,p6,p4-jstart+1,p5-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_6(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
      do h3=1,h3d
      do h1=1,h1d
       triplesx(h1,h3,h2,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h1,h3,h2,p6,p4-jstart+1,p5-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_7(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h2d,h1d,p4d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
       triplesx(h3,h2,h1,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h3,h2,h1,p4,p6-jstart+1,p5-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_8(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p4d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h2=1,h2d
      do h1=1,h1d
      do h3=1,h3d
       triplesx(h3,h1,h2,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h3,h1,h2,p4,p6-jstart+1,p5-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_6ds_9(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p4d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h2=1,h2d
      do h3=1,h3d
      do h1=1,h1d
       triplesx(h1,h3,h2,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h1,h3,h2,p4,p6-jstart+1,p5-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end



      subroutine sd_tx1_1(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h2d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
      do h7=1,h7d
       triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)
     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_tx1_2(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p6d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
      do h1=1,h1d
      do h3=1,h3d
      do h7=1,h7d
       triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)
     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo  
      enddo
      return
      end
c
      subroutine sd_tx1_3(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d,h2d,p6d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
      do h3=1,h3d
      do h1=1,h1d
      do h7=1,h7d
       triplesx(h1,h3,h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h1,h3,h2,p6,p5-jstart+1,p4-istart+1)
     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_tx1_4(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h2d,h1d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p6=istart,istop
      do p4=jstart,jstop
      do p5=1,p5d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
      do h7=1,h7d
       triplesx(h3,h2,h1,p5,p4-jstart+1,p6-istart+1)=
     1     triplesx(h3,h2,h1,p5,p4-jstart+1,p6-istart+1)
     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_tx1_5(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p6=istart,istop
      do p4=jstart,jstop
      do p5=1,p5d
      do h2=1,h2d
      do h1=1,h1d
      do h3=1,h3d
      do h7=1,h7d
       triplesx(h3,h1,h2,p5,p4-jstart+1,p6-istart+1)=
     1     triplesx(h3,h1,h2,p5,p4-jstart+1,p6-istart+1)
     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_tx1_6(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d,h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p6=istart,istop
      do p4=jstart,jstop
      do p5=1,p5d
      do h2=1,h2d
      do h3=1,h3d
      do h1=1,h1d
      do h7=1,h7d
       triplesx(h1,h3,h2,p5,p4-jstart+1,p6-istart+1)=
     1     triplesx(h1,h3,h2,p5,p4-jstart+1,p6-istart+1)
     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_tx1_7(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h2d,h1d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p4=istart,istop
      do p6=jstart,jstop
      do p5=1,p5d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
      do h7=1,h7d
       triplesx(h3,h2,h1,p5,p6-jstart+1,p4-istart+1)=
     1     triplesx(h3,h2,h1,p5,p6-jstart+1,p4-istart+1)
     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_tx1_8(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p4=istart,istop
      do p6=jstart,jstop
      do p5=1,p5d
      do h2=1,h2d
      do h1=1,h1d
      do h3=1,h3d
      do h7=1,h7d
       triplesx(h3,h1,h2,p5,p6-jstart+1,p4-istart+1)=
     1     triplesx(h3,h1,h2,p5,p6-jstart+1,p4-istart+1)
     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_tx1_9(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d,h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      do p4=istart,istop
      do p6=jstart,jstop
      do p5=1,p5d
      do h2=1,h2d
      do h3=1,h3d
      do h1=1,h1d
      do h7=1,h7d
       triplesx(h1,h3,h2,p5,p6-jstart+1,p4-istart+1)=
     1     triplesx(h1,h3,h2,p5,p6-jstart+1,p4-istart+1)
     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_1(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
      do p7=1,p7d
       triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)
     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_2(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h1d,h3d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h3=1,h3d
      do h1=1,h1d
      do h2=1,h2d
      do p7=1,p7d
       triplesx(h2,h1,h3,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h2,h1,h3,p6,p5-jstart+1,p4-istart+1)
     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_3(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h3=1,h3d
      do h2=1,h2d
      do p7=1,p7d
       triplesx(h2,h3,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h2,h3,h1,p6,p5-jstart+1,p4-istart+1)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_4(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
      do p7=1,p7d
       triplesx(h3,h2,h1,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h3,h2,h1,p6,p4-jstart+1,p5-istart+1)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_5(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h1d,h3d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h3=1,h3d
      do h1=1,h1d
      do h2=1,h2d
      do p7=1,p7d
       triplesx(h2,h1,h3,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h2,h1,h3,p6,p4-jstart+1,p5-istart+1)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_6(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h3=1,h3d
      do h2=1,h2d
      do p7=1,p7d
       triplesx(h2,h3,h1,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h2,h3,h1,p6,p4-jstart+1,p5-istart+1)
     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_7(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p4d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h1=1,h1d
      do h2=1,h2d
      do h3=1,h3d
      do p7=1,p7d
       triplesx(h3,h2,h1,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h3,h2,h1,p4,p6-jstart+1,p5-istart+1)
     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_8(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h1d,h3d,p4d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h3=1,h3d
      do h1=1,h1d
      do h2=1,h2d
      do p7=1,p7d
       triplesx(h2,h1,h3,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h2,h1,h3,p4,p6-jstart+1,p5-istart+1)
     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
c
      subroutine sd_ty1_9(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p4d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h1=1,h1d
      do h3=1,h3d
      do h2=1,h2d
      do p7=1,p7d
       triplesx(h2,h3,h1,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h2,h3,h1,p4,p6-jstart+1,p5-istart+1)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
