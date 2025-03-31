*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.16/08 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.46.45  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.11  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.01  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE COHEREN
*KEEP,gplhint.
*KEND.

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IFREQ,ISOUR
      DOUBLE PRECISION OSZILL,OSZMAX,OSZMIN
      DOUBLE PRECISION DTIME,DTIMEP,DTIM

      OSZMAX=-1.D30
      OSZMIN= 1.D30

      DO IFREQ=1,NFREQ
          DO ISOUR=1,NSOURCE

          DTIME=SOURCET(2,ISOUR)-SOURCET(1,ISOUR)
          DTIMEP=(SOURCEE(1,1,ISOUR)-SOURCEA(1,1,ISOUR))/CLIGHT1

          DTIM=DTIME-DTIMEP

         OSZILL=FREQ(IFREQ)/HBAREV1*DTIM/2.D0/PI1
         IF (OSZILL.LT.OSZMIN) OSZMIN=OSZILL
         IF (OSZILL.GT.OSZMAX) OSZMAX=OSZILL

          ENDDO   !ISOUR
      ENDDO !IFREQ

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR COHEREN:'
      WRITE(LUNGFO,*)
     &'     Minimum and maximum number of phases  between'
      WRITE(LUNGFO,*)
     &'     photons and electrons:'
      WRITE(LUNGFO,*)'     ',SNGL(OSZMIN),SNGL(OSZMAX)
      WRITE(LUNGFO,*)

      RETURN
      END
