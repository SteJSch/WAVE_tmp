*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 12/11/2009  16.27.11  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.08.53  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.02.58  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.27  by  Michael Scheer
*-- Author : Michael Scheer
C****************************************************************************
      SUBROUTINE POWINT(X,Y,IWALL,IMODE,IPOL,ISTART,N)
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER IWALL,IPOL,IMODE,IWALLO,IMODEO,IPOLO,N,ISTART

      DATA IMODEO/0/
      DATA IPOLO/0/
      DATA IWALLO/0/

      REAL*4 X,Y

C--- SPLINE COEFFICIENTS

      IF (
     &      IWALLO.NE.IWALL
     &  .OR.
     &      IPOLO.NE.IPOL
     &  .OR.
     &      IMODEO.NE.IMODE
     &  ) THEN

          CALL POWSPL(IWALL,IMODE,IPOL,ISTART,N,1.D30,1.D30)
          IWALLO=IWALL
          IPOLO=IPOL
          IMODEO=IMODE

      ENDIF !ICAL

          CALL POWSPI(X,Y,IWALL,IMODE,IPOL,ISTART,N)

      RETURN
      END
