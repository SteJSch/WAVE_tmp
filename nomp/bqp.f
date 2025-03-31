*CMZ :  2.67/01 15/03/2012  10.51.33  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  1.02/00 18/12/97  13.35.57  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.24.44  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  19.35.51  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.43  by  Michael Scheer
*-- Author : Michael Scheer
C***********************************************************************
      SUBROUTINE BQP(XI,YI,ZI,BX,BY,BZ,IMAG)
*KEEP,gplhint.
*KEND.
      IMPLICIT NONE
      INTEGER IMAG

*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION Y,Z,BX,BY,BZ,G,XI,YI,ZI,PHI,SPHI,CPHI,XCEN,ZCEN,DX,DZ
      DOUBLE PRECISION DXA,DZA,DXE,DZE,XA,XE,ZA,ZE,ADUMA,ADUME

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      PHI=PMAG(5,IMAG)
      CPHI=DCOS(PHI)
      SPHI=DSIN(PHI)
      XCEN=PMAG(3,IMAG)
      ZCEN=PMAG(4,IMAG)

      XA=XCEN-CPHI*PMAG(1,IMAG)/2.
      ZA=ZCEN-SPHI*PMAG(1,IMAG)/2.
      XE=XCEN+CPHI*PMAG(1,IMAG)/2.
      ZE=ZCEN+SPHI*PMAG(1,IMAG)/2.
      DXA=XI-XA
      DZA=ZI-ZA
      DXE=XI-XE
      DZE=ZI-ZE

      ADUMA=DXA*CPHI+DZA*SPHI
      ADUME=DXE*CPHI+DZE*SPHI

      IF (IWFILF.NE.99.AND.ADUMA*ADUME.LE.0.) THEN

        DX=XI-XCEN
        DZ=ZI-ZCEN
        Z=-SPHI*DX+CPHI*DZ
        Y=YI

        G=PMAG(2,IMAG)*EMOM/CLIGHT1

        BX=0.
        BY=G*Z
        BZ=G*Y

      ELSE

        BX=0.
        BY=0.
        BZ=0.

      ENDIF !IWFILF

      RETURN
      END
