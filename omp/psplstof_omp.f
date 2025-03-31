*CMZ :  3.07/00 15/03/2019  12.15.26  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.34.43  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.07  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE psplstof_omp(ISTOK,kfreq)
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
C    FOLDED INTENSITY INSIDE THE PINHOLE

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

      INTEGER ISTOK,kfreq,IOBSV,IZ,IY,IIZ

      DOUBLE PRECISION S(NSPLINEP),S2(NSPLINEP)

C- SPLINES IN Z

      DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY

        IIZ=0
        DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ
          IIZ=IIZ+1
          IOBSV=(IY-1)*NOBSVZ+IZ
          S(IIZ)=STOKESF(ISTOK,IOBSV+NOBSV*(kfreq-1))
        ENDDO      !IZ

        CALL FSPLINDX_omp(OBSVDZ,S,MOBSVZ,0.D0,0.D0,S2)

        IIZ=0
        DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ
          IOBSV=(IY-1)*NOBSVZ+IZ
          IIZ=IIZ+1
          SPCOEFM(IOBSV)=S2(IIZ)
        ENDDO      !IZ
      ENDDO      !IY

      RETURN
      END
