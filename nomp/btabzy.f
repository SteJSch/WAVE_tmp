*CMZ :  4.00/11 17/05/2021  11.08.35  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.54/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/08 31/05/95  13.30.40  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  16.00.44  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.02.31  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.57  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BTABZY(XIN,DUMY,DUMZ,BX,BY,BZ,AX,AY,AZ)
*KEEP,gplhint.
*KEND.

C     BTABZY reads By and Bz from data files and interpolates fields

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      INTEGER ICAL

      DOUBLE PRECISION XIN,DUMY,DUMZ,DUMZY,XS,XE,XS1,XS2,XE1,XE2
      DOUBLE PRECISION AX,AY,AZ,BX,BY,BZ
      DOUBLE PRECISION AXX,AYY,AZZ,BXX,BYY,BZZ
      DOUBLE PRECISION AXXX,AYYY,AZZZ,BXXX,BYYY,BZZZ

      DATA ICAL/0/
      DATA XS1,XS2,XE1,XE2/4*0.0/

      DUMZY=DUMY
      DUMZY=DUMZ

      BXX=0.
      BYY=0.
      BZZ=0.

      AXX=0.
      AYY=0.
      AZZ=0.

      BXXX=0.
      BYYY=0.
      BZZZ=0.

      AXXX=0.
      AYYY=0.
      AZZZ=0.

      XS=XSTART
      XE=XSTOP


      CALL BTAB (XIN,0.D0,0.D0, BXX, BYY, BZZ, AXX, AYY, AZZ)

      IF (ICAL.EQ.0) THEN

         XS1=XSTART
         XE1=XSTOP

      ENDIF !ICAL

      CALL BTABZ(XIN,0.D0,0.D0,BXXX,BYYY,BZZZ,AXXX,AYYY,AZZZ,XS,XE)

      IF (ICAL.EQ.0) THEN

          XS2=XSTART
          XE2=XSTOP

          IF (XS1.NE.XS2) XSTART=DMIN1(XS1,XS2)
          IF (XE1.NE.XE2)  XSTOP=DMAX1(XE1,XE2)

          ICAL=1

      ENDIF !ICAL

      BX=BXX+BXXX
      BY=BYY+BYYY
      BZ=BZZ+BZZZ

      AX=AXX+AXXX
      AY=AYY+AYYY
      AZ=AZZ+AZZZ

      RETURN
      END
