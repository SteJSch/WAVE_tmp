*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.37/02 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  15.42.20  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.19.09  by  Michael Scheer
*CMZ :  2.13/07 17/02/2000  15.11.12  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.01.43  by  Michael Scheer
*CMZ :  2.12/03 07/07/99  12.30.49  by  Michael Scheer
*CMZ :  2.10/01 25/02/99  16.08.38  by  Michael Scheer
*-- Author :    Michael Scheer   25/02/99

      SUBROUTINE POWERCIRC(ISOUR)
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

      INTEGER ISOUR,ICAL,IOBSV,IOBSVZ

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

         DO IOBSV=1,NOBSV
             SPECPOWT(IOBSV)=0.D0
         ENDDO

         DO IOBSVZ=1,NOBSVZ
             SPECPOWVT(IOBSVZ)=0.D0
         ENDDO

         SPECPOWVHT=0.D0

         ICAL=1

      ENDIF !ICAL

      IF (IRPHI.EQ.0) THEN
          CALL CIRCPOW(ISOUR,NSOURCE)
      ELSE  !(IRPHI.EQ.0)
                   CALL CIRCPOWRP(ISOUR)
      ENDIF !(IRPHI.EQ.0)

      DO IOBSV=1,NOBSV
          SPECPOWT(IOBSV)=SPECPOWT(IOBSV)+SPECPOW(ISOUR+NSOURCE*(IOBSV-1))
      ENDDO

      DO IOBSVZ=1,NOBSVZ
            ILIOBZ=ISOUR+NSOURCE*(IOBSVZ-1)
          SPECPOWVT(IOBSVZ)=SPECPOWVT(IOBSVZ)+SPECPOWV(ILIOBZ)
          SPECPOWT(IOBSVZ)=SPECPOWT(IOBSVZ)+SPECPOW(ILIOBZ)
      ENDDO

      SPECPOWVHT=SPECPOWVHT+SPECPOWVH(ISOUR)

      RETURN
      END
