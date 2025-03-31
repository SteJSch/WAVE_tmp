*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.52/11 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.02/00 18/12/97  13.35.57  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.24.44  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  19.35.51  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.43  by  Michael Scheer
*-- Author : Michael Scheer
C***********************************************************************
      SUBROUTINE BQF(XI,YI,ZI,BX,BY,BZ,IMAG)
*KEEP,gplhint.
*KEND.
      IMPLICIT NONE
      INTEGER IMAG

*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION Y,Z,BX,BY,BZ,G,XI,YI,ZI,PHI,SPHI,CPHI,XCEN,ZCEN,DX,DZ
      DOUBLE PRECISION AY1,AY2,XLEN2,BY1,BY2,FRINGE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DATA FRINGE/1000.0D0/

      PHI=PMAG(5,IMAG)
      CPHI=DCOS(PHI)
      SPHI=DSIN(PHI)
      XCEN=PMAG(3,IMAG)
      ZCEN=PMAG(4,IMAG)

      XLEN2=PMAG(1,IMAG)/2.0D0

      IF (IWFILF.NE.99) THEN

        DX=XI-XCEN
        DZ=ZI-ZCEN
        Z=-SPHI*DX+CPHI*DZ
        Y=YI

        AY1=(DX-XLEN2)*FRINGE
        AY2=(-DX-XLEN2)*FRINGE

        IF (AY1.GT.70.) THEN
          BY1=0.0D0
        ELSE IF (AY1.LT.-70.) THEN
          BY1=1.0D0
        ELSE
          BY1=1.0D0/(1.0D0+DEXP(AY1))
        ENDIF

        IF (AY2.GT.70.) THEN
          BY2=0.0D0
        ELSE IF (AY2.LT.-70.) THEN
          BY2=1.0D0
        ELSE
          BY2=1.0D0/(1.0D0+DEXP(AY2))
        ENDIF

        G=PMAG(2,IMAG)*EMOM/CLIGHT1*BY1*BY2

        BX=0.0D0
        BY=G*Z
        BZ=G*Y

      ELSE

        BX=0.0D0
        BY=0.0D0
        BZ=0.0D0

      ENDIF !IWFILF

      RETURN
      END
