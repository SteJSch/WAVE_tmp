*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/07 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.16/08 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.01.12  by  Michael Scheer
*CMZ :  1.00/00 08/07/97  10.51.15  by  Michael Scheer
*CMZ : 00.01/10 02/06/96  12.03.03  by  Michael Scheer
*CMZ : 00.01/06 20/02/95  16.06.18  by  Michael Scheer
*CMZ : 00.01/05 01/02/95  14.05.08  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  17.53.53  by  Michael Scheer
*CMZ : 00.00/07 24/05/94  09.48.27  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UNAME_FIT
*KEEP,gplhint.
*KEND.

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,modulator.
      include 'modulator.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      INTEGER I
      DOUBLE PRECISION VARFIT(100)

      OPEN(UNIT=99,FILE='UNAME.FIT',STATUS='OLD')
      DO I=1,100
          READ(99,*,END=99)VARFIT(I)
            VARFIT(I)=ASIN(VARFIT(I))*radgra1
      ENDDO
99    CLOSE(99)

      DO I=1,NMAGMOD

          IF (THEROT(I).EQ.11111.) THEN
         THEROT(I)=VARFIT(1)
          else if (THEROT(I).EQ.-11111.) THEN
         THEROT(I)=-VARFIT(1)

          else if (THEROT(I).EQ.22222.) THEN
         THEROT(I)=VARFIT(2)
          else if (THEROT(I).EQ.-22222.) THEN
         THEROT(I)=-VARFIT(2)

          else if (THEROT(I).EQ.33333.) THEN
         THEROT(I)=VARFIT(3)
          else if (THEROT(I).EQ.-33333.) THEN
         THEROT(I)=-VARFIT(3)

          ENDIF

      ENDDO

      RETURN
      END
