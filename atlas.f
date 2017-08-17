c	file ATLAS.F
c
c	Program to read Liegi Solar Spectrum data file /d1/iraf/atlas/liegi.da
ct
c	The Range of the Atlas runs from 3601 Angstroms to 9300 Anstroms, with
c	resolution of 2 MiliAngstrom.
      logical*1 buf1(2049), buf2(2049)
      real buf1a(512), buf2a(512)
      real lambda(1000), ri(1000)
      equivalence (buf1a, buf1), (buf2a, buf2)
# 12 "atlas.for"
      open(unit=9, file='/usr/local/lib/atlas/liegi/data/liegi.dat', 
     &status='old', access='direct', err=9001, recl=8196) 
   10 write(unit=*, fmt=*) '(0)Exit,(1)Read Spectrum,(2)Plot Spectrum,'
      write(unit=*, fmt=*) '(3)Output Spectrum(4)Dump Spectrum.'
      read(unit=*, fmt=*) iop
      iop = iop + 1
      goto (999, 100, 200, 300, 400), iop
c*******************************************************************
c	Read Spectrum
# 19 "atlas.for"
      goto 10
# 24 "atlas.for"
  100 write(unit=*, fmt=101) 
  101 format(35hEnter Spectrum Range in Angstrom...,$)
      read(unit=*, fmt=*) lambda1, lambda2
      ist = lambda1 - 3600
      ino = (lambda2 - lambda1) + 1
      if (ist .lt. 1) then
      write(unit=*, fmt=*) 
     &'INPUT ERROR...Spectrum Range : 3601 A <-> 9300 A'
# 31 "atlas.for"
      goto 10
      end if
      j1 = 0
      do 103 i = ist, ino, 2
      read(unit=9, rec=i, err=9002) buf1
      read(unit=9, rec=i + 1, err=9002) buf2
      do 102 j = 48, 2047
      j1 = j1 + 1
      lambda(j1) = buf1a(j)
      ri(j1) = buf2a(j)
  102 continue
  103 continue
c*******************************************************
c	Plot Spectrum
# 43 "atlas.for"
      goto 10
c*******************************************************
c	Output Spectrum
# 48 "atlas.for"
  200 continue
c*******************************************************
c	Dump Spectrum
# 53 "atlas.for"
  300 continue
# 58 "atlas.for"
  400 write(unit=*, fmt=*) 'Enter Range...'
      read(unit=*, fmt=*) ist, iend
      do 401 i = ist, iend
      write(unit=*, fmt=*) lambda(i), ri(i)
  401 continue
      goto 10
# 65 "atlas.for"
 9001 write(unit=*, fmt=*) 'Error Open File Liegi.dat...'
      goto 10
 9002 write(unit=*, fmt=*) 'Error Reading File...'
      goto 10
# 70 "atlas.for"
  999 close(unit=9) 
      end
