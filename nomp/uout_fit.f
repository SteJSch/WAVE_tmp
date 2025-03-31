*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.13/03 17/12/99  11.45.46  by  Michael Scheer
*CMZ :  1.00/00 10/07/97  13.57.51  by  Michael Scheer
*CMZ : 00.01/10 02/06/96  12.04.24  by  Michael Scheer
*CMZ : 00.01/06 20/02/95  16.06.00  by  Michael Scheer
*CMZ : 00.01/05 01/02/95  15.24.06  by  Michael Scheer
*CMZ : 00.01/04 26/01/95  15.55.24  by  Michael Scheer
*CMZ : 00.00/07 24/05/94  09.48.34  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UOUT_FIT
*KEEP,gplhint.
*KEND.

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.
      OPEN(UNIT=99,FILE='UOUT.FIT',STATUS='OLD')
         WRITE(99,*)SPEC(1)
         WRITE(99,*)WTRA2I
         WRITE(99,*)WTRA(2,1,NCO)
         WRITE(99,*)WTRA(3,1,NCO)
         WRITE(99,*)WTRA(2,2,NCO)/WTRA(1,2,NCO)
         WRITE(99,*)WTRA(3,2,NCO)/WTRA(1,2,NCO)
      CLOSE(99)

      RETURN
      END
