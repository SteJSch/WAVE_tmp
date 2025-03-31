*CMZ :  3.00/00 11/03/2013  15.13.36  by  Michael Scheer
*CMZ :  2.16/08 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 10/05/2000  16.54.10  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  14.59.11  by  Michael Scheer
*CMZ :  2.13/03 18/01/2000  17.55.50  by  Michael Scheer
*CMZ :  1.03/06 25/09/98  11.35.13  by  Michael Scheer
*-- Author :    Michael Scheer   18/09/98

      SUBROUTINE PHASEANA
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,reargf90.
      include 'reargf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEND.

      INTEGER I

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** PHASEANA CALLED!! ***'
      WRITE(LUNGFO,*)'ARRAY REAIMA OVERWRITTEN!!'
      WRITE(LUNGFO,*)

      WRITE(6,*)
      WRITE(6,*)'*** PHASEANA CALLED!! ***'
      WRITE(6,*)'ARRAY REAIMA OVERWRITTEN!!'
      WRITE(6,*)

      DO I=1,NOBSV*NFREQ
          REAIMA(1,1,I)=0.D0
          REAIMA(2,1,I)=0.D0
          REAIMA(3,1,I)=0.D0
          REAIMA(1,2,I)=0.D0
          REAIMA(2,2,I)=0.D0
          REAIMA(3,2,I)=0.D0
      ENDDO
      REAIMA(2,1,NOBSV/2+1)=1.D0

      RETURN
      END
