C     This is a simple program to determine the units of the Fortran
c     direct access record length.  It writes 2 records of 512 units
c     each.  Simply run this program and check the size resulting file
c     'da_rec_size.test' (in bytes) using 'ls' and divide by 1024 units.
C
      Program da_rec_size
C$Id$
      Implicit NONE
      Integer RecLen
      Parameter (RecLen = 512)
      Integer A(1), I
      Open(1, FILE='da_rec_size.test', ACCESS='DIRECT', RECL=RecLen)
      Do I = 1, 2
         Write (1, REC=I) A
      EndDo
      Write (6, *) 'This program has written two 512-unit records in'
      Write (6, *) 'the file da_rec_size.test.  To determine the size'
      Write (6, *) 'of the record length unit, divide the length of'
      Write (6, *) 'the file in bytes by 1024.  The most common sizes'
      Write (6, *) 'for record size units are bytes and integer words.'
      Stop
      End
