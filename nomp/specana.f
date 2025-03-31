*CMZ :  3.03/01 12/11/2015  12.12.46  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.57/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.36/00 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.34/09 26/09/2001  16.25.11  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  16.21.10  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  14.24.02  by  Michael Scheer
*CMZ :  2.13/03 15/12/99  16.21.34  by  Michael Scheer
*CMZ :  2.13/00 03/12/99  14.58.39  by  Michael Scheer
*CMZ : 00.02/04 25/02/97  16.04.42  by  Michael Scheer
*CMZ : 00.02/03 30/01/97  15.29.53  by  Michael Scheer
*CMZ : 00.00/07 18/05/94  14.54.15  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.12  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.24  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPECANA
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- FILL ARRAY SPEC WITH USER SUPPLIED FUNCTION

      IMPLICIT NONE

      INTEGER ICONTROL
      INTEGER IFREQ,IOBSV,ISOUR

      DOUBLE PRECISION SIGY

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
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      DATA ICONTROL/1/
      DATA SIGY/0.0001D0/

      ICONTROL=USER(1)
      sigy=wsigy(1)

      IF (ICONTROL.EQ.1) THEN

        write(lungfo,*) ' '
        write(lungfo,*) '*** SPECANA: ICONTROL, SIG:',ICONTROL,SIGY
        write(lungfo,*) ' '

        DO IFREQ=1,NFREQ
          DO IOBSV=1,NOBSV
            DO ISOUR=1,NSOURCE
              SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))=
     &          EXP(-(OBSV(2,IOBSV)/SIGY)**2/2.D0)/SQRT(2.D0*PI1)/SIGY
              SPECPOW(ISOUR+(IOBSV-1)*NSOURCE)=
     &          SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
            ENDDO
            STOKES(1,IOBSV+NOBSV*(IFREQ-1))=
     &        SPEC(1+IOBSV-1+NOBSV*(IFREQ-1))
          ENDDO
        ENDDO

      else IF (ICONTROL.EQ.2) THEN

        write(lungfo,*) ' '
        write(lungfo,*) '*** SPECANA: ICONTROL, SIG:',ICONTROL,SIGY
        write(lungfo,*) ' '

        DO IFREQ=1,NFREQ
          DO IOBSV=1,NOBSV
            DO ISOUR=1,NSOURCE
              SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))=
     &          EXP(-(OBSV(3,IOBSV)/SIGY)**2/2.D0)/SQRT(2.D0*PI1)/SIGY
              SPECPOW(ISOUR+(IOBSV-1)*NSOURCE)=
     &          SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
            ENDDO
            STOKES(1,IOBSV+NOBSV*(IFREQ-1))=
     &        SPEC(1+IOBSV-1+NOBSV*(IFREQ-1))
          ENDDO
        ENDDO

      ELSE  !ICONTROL

        DO IFREQ=1,NFREQ
          DO IOBSV=1,NOBSV

            IF (ISTOKES.NE.0) THEN
              STOKES(1,IOBSV+NOBSV*(IFREQ-1))=1.
            ENDIF

            DO ISOUR=1,NSOURCE
              SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))=1.D0
              SPECPOW(ISOUR+(IOBSV-1)*NSOURCE)=1.D0
            ENDDO
          ENDDO
        ENDDO

      ENDIF !ICONTROL

      RETURN
      END
