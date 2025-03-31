*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.14/02 27/04/2000  15.53.13  by  Michael Scheer
*CMZ : 00.00/07 18/05/94  14.54.44  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.17  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.25  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPECANAF
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- FILL ARRAY SPECF WITH USER SUPPLIED FUNCTION

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
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IFREQ,IOBSV,ISOUR

        DO IFREQ=1,NFREQ
        DO IOBSV=1,NOBSV
        DO ISOUR=1,NSOURCE
            SPECF(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))=1.D0
        ENDDO
        ENDDO
        ENDDO


      RETURN
      END
