C     This is a simple program to determine the units of the Fortran
c     direct access record length.  It writes 2 records of 512 units
c     each.  Simply run this program and check the size resulting file
c     'da_rec_size.test' (in bytes) using 'ls' and divide by 1024 units.
C
      Program da_rec_size
C$Id: da_rec_size.f,v 1.2 1995-02-02 23:10:44 d3g681 Exp $
      Integer RecLen
      Parameter (RecLen = 512)
      Integer A(1)
      Open(1, FILE='da_rec_size.test', ACCESS='DIRECT', RECL=RecLen)
      Do I = 1, 2
         Write (1, REC=I) A
      EndDo
      Stop
      End
