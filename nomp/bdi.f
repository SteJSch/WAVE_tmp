*CMZ :  3.04/00 05/01/2018  16.13.44  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.65/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.63/02 25/03/2008  09.35.21  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.02/00 19/12/97  17.58.18  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.07.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.47  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.43  by  Michael Scheer
*-- Author : Michael Scheer
C***********************************************************************
      SUBROUTINE BDI(XI,YI,ZI,BX,BY,BZ,IMAG)

*KEEP,gplhint.
*KEND.

      IMPLICIT NONE

      INTEGER IMAG

*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,BY0,BY1,BY2,XLEN2,XI,YI,ZI,AY1,AY2

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

c pmag(1,imag): deflection angle (rad)
c pmag(2,imag): bending radius (T)
c pmag(3,imag): Center of magnet (m)
c pmag(4,imag): Width of edge

      IF (
     &    (IWFILF.EQ.99. .OR. IMGSQF.EQ.0)
     &    .AND.
     &    PMAG(1,IMAG)*PMAG(2,IMAG).NE.0.
     &    ) THEN

        X=XI-PMAG(3,IMAG)

        Y=YI
        Z=ZI

        BX=0.
        BZ=0.

        BY0=EMOM/CLIGHT1/PMAG(2,IMAG)
        XLEN2=DABS(PMAG(2,IMAG)*sin(pmag(1,imag)/2.0d0))


        AY1=(+X-XLEN2)*PMAG(4,IMAG)
        AY2=(-X-XLEN2)*PMAG(4,IMAG)

        IF (AY1.GT.70.0D0) THEN
          BY1=0.0d0
        ELSE IF (AY1.LT.-70.) THEN
          BY1=1.0D0
        ELSE
          BY1=1.0D0/(1.0D0+DEXP(AY1))
        ENDIF

        IF (AY2.GT.70.0D0) THEN
          BY2=0.0D0
        ELSE IF (AY2.LT.-70.0D0) THEN
          BY2=1.0D0
        ELSE
          BY2=1.0D0/(1.0D0+DEXP(AY2))
        ENDIF

        BY=BY0*BY1*BY2*CORR(IMAG)

      ELSE

        BX=0.0d0
        BY=0.0d0
        BZ=0.0d0

      ENDIF

      RETURN
      END
