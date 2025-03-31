*CMZ :  3.07/00 11/03/2019  13.47.13  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.13/03 15/12/99  16.42.06  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.15  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.52  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE PSPLINE(ISOUR,IFREQ)
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- CALCULATES COEFFICIENTS OF CUBIC SPLINES THAT INTERPOLATE THE
C    INTENSITY INSIDE THE PINHOLE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER ISOUR,IFREQ,IOBSV,IZ,IY

      DOUBLE PRECISION S(NSPLINEP),S2(NSPLINEP)
      DOUBLE PRECISION YPP0,YPPN

C- SPLINES IN Z

      IOBSV=0
      DO IY=1,NOBSVY

        DO IZ=1,NOBSVZ
          IOBSV=IOBSV+1
          S(IZ)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
        ENDDO      !IZ

c        YPP0=1.D30
c        YPPN=1.D30

C       CALL FSPLINEZ(OBSVZ,S,NOBSVZ,YPP0,YPPN,S2)
C060793    CALL FSPLINDX(OBSVZ, S,NOBSVZ,YPP0,YPPN,S2)
        CALL FSPLINDX(OBSVDZ,S,NOBSVZ,0.D0,0.D0,S2)

        DO IZ=1,NOBSVZ
          SPCOEF(IZ+(IY-1)*NOBSVZ)=S2(IZ)
        ENDDO      !IZ
      ENDDO      !IY

      RETURN
      END
