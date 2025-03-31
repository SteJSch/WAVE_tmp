*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.56.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.20  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.34  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE CONVUN
*KEEP,gplhint.
*KEND.

      IMPLICIT NONE

C--- CONVERTS FREQUENCES ON INPUT FILE FROM EV TO NM OR VICE VERSA

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IFREQ
      DOUBLE PRECISION FBUFF(NDFREQP)

      DO IFREQ=1,NFREQ
          FBUFF(IFREQ)=WTOE1/FREQ(IFREQ)
      ENDDO

      DO IFREQ=1,NFREQ
             FREQ(IFREQ)=FBUFF(NFREQ+1-IFREQ)
      ENDDO


      IF (FREQC.NE.0.0) FREQC=WTOE1/FREQC
      IF (FREQCF.NE.0.0) FREQCF=WTOE1/FREQCF

      RETURN
      END
