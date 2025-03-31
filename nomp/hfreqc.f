*CMZ :  4.00/13 07/12/2021  18.47.10  by  Michael Scheer
*CMZ :  3.03/02 07/12/2015  17.18.10  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/03 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.65/01 08/10/2009  09.58.11  by  Michael Scheer
*CMZ :  2.63/05 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.61/03 27/03/2007  13.23.57  by  Michael Scheer
*CMZ :  2.53/03 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/02 25/01/2005  18.15.44  by  Michael Scheer
*CMZ :  2.52/13 09/12/2004  13.09.09  by  Michael Scheer
*CMZ :  2.52/00 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.48/04 17/03/2004  13.47.13  by  Michael Scheer
*CMZ :  2.40/02 14/03/2002  15.46.20  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.16/07 01/09/2000  14.36.07  by  Michael Scheer
*CMZ :  2.16/04 28/06/2000  17.43.04  by  Michael Scheer
*CMZ :  2.16/01 15/06/2000  15.47.11  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.17.53  by  Michael Scheer
*CMZ :  2.13/04 24/01/2000  16.35.38  by  Michael Scheer
*CMZ :  2.13/03 10/01/2000  17.32.03  by  Michael Scheer
*CMZ :  2.13/00 02/12/99  13.23.58  by  Michael Scheer
*CMZ :  1.03/00 16/01/98  11.16.49  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.27  by  Michael Scheer
*CMZ : 00.01/06 14/02/95  10.59.52  by  Michael Scheer
*CMZ : 00.01/05 31/01/95  17.02.40  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  13.06.10  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.43.29  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.08  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.39  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HFREQC
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- HISTOGRAMS FOR SPECTRA OF SINGLE OBSERVATION POINTS OR PINHOLE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      IF (IUNIT.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING SR HFREQC: NO HISTOGRAMS FOR IUNIT.NE.0'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
        RETURN
      ENDIF

c      IF (IFREQ2P.LE.1) THEN
c
c        WRITE(LUNGFO,*)
c        WRITE(LUNGFO,*)'*** WARNING SR HFREQC ***'
c        WRITE(LUNGFO,*)'PARAMETER IFREQ2P .LE. 1, NO HISTOGRAM BOOKED'
c        WRITE(LUNGFO,*)
c
c      ELSE IF (IFREQ2P.EQ.2) THEN
      IF (IFREQ2P.EQ.2) THEN

        call hfreqc2

      ELSE  !IFREQ2P

        call hfreqc1

      ENDIF !IFREQ2P

      RETURN
      END
