*CMZ :  3.04/00 04/01/2018  17.28.38  by  Michael Scheer
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
      SUBROUTINE BDh(XI,YI,ZI,BX,BY,BZ,IMAG)

*KEEP,gplhint.
*KEND.

      IMPLICIT NONE

      INTEGER IMAG

*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,BZ0,BZ1,BZ2,XLEN2,XI,YI,ZI,AZ1,AZ2

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
     &    IWFILF.ne.99.
     &    .AND.
     &    PMAG(1,IMAG)*PMAG(2,IMAG).NE.0.
     &    ) THEN

        X=XI-PMAG(3,IMAG)

        Y=YI
        Z=ZI

        BX=0.0d0
        BY=0.0d0

        BZ0=EMOM/CLIGHT1/PMAG(2,IMAG)
        XLEN2=DABS(PMAG(2,IMAG)*sin(pmag(1,imag)/2.0d0))


        AZ1=(+X-XLEN2)*PMAG(4,IMAG)
        AZ2=(-X-XLEN2)*PMAG(4,IMAG)

        IF (AZ1.GT.70.0D0) THEN
          BZ1=0.
        ELSE IF (AZ1.LT.-70.) THEN
          BZ1=1.00D0
        ELSE
          BZ1=1.0D0/(1.0D0+DEXP(AZ1))
        ENDIF

        IF (AZ2.GT.70.0D0) THEN
          BZ2=0.0D0
        ELSE IF (AZ2.LT.-70.0D0) THEN
          BZ2=1.0D0
        ELSE
          BZ2=1.0D0/(1.0D0+DEXP(AZ2))
        ENDIF

        BZ=BZ0*BZ1*BZ2*CORR(IMAG)

      ELSE

        Bx=0.0d0
        By=0.0d0
        BZ=0.0d0

      ENDIF

      RETURN
      END
