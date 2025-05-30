*CMZ :  4.00/13 26/08/2021  09.23.14  by  Michael Scheer
*CMZ :  3.08/01 04/04/2019  07.58.25  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  16.48.27  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.20.42  by  Michael Scheer
*CMZ :  2.52/01 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.51/00 24/05/2004  20.37.45  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  15.24.39  by  Michael Scheer
*CMZ :  2.33/00 02/05/2001  15.58.03  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.13/03 10/01/2000  17.17.14  by  Michael Scheer
*CMZ :  2.13/00 29/11/99  12.18.03  by  Michael Scheer
*CMZ :  2.12/03 07/07/99  12.27.14  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.47  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SOUADD
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- SUMS UP THE CONTRIBUTIONS OF ALL SOURCES

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER kfreq,IOBSV,ISOUR

C--- LOOP OVER ALL FREQUENCES

      DO kfreq=1,NFREQ

C--- LOOP OVER ALL OBSERVATION POINTS

        DO IOBSV=1,NOBSV

C--- LOOP OVER ALL SOURCES

          SPECTOT(IOBSV+NOBSV*(kfreq-1))=0.0d0
          DO ISOUR=1,NSOURCE

            IF (kfreq.EQ.1.AND.IOBSV.EQ.1.AND.IPIN.NE.2) THEN
              IF (IPINCIRC.EQ.0.or.ipin.eq.3) THEN
                CALL POWER(ISOUR)
              ELSE
                CALL POWERCIRC(ISOUR)
              ENDIF
            ENDIF !(kfreq.EQ.1.AND.IOBSV.EQ.1)

            ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1))

            IF( ISPECDIP.EQ.0.AND.SPECCUT.GT.0.0) then
              if(FREQ(kfreq).GT.SPECCUT*ecdipev1*DMYENERGY**2*ECMAX(ISOUR)) THEN
                SPEC(ILIOBFR)=0.0D0
              ENDIF
            ENDIF

            IOBFR=IOBSV+NOBSV*(kfreq-1)
            SPECTOT(IOBFR)=
     &        SPECTOT(IOBFR)+SPEC(ILIOBFR)

          ENDDO !LOOP OVER ALL SOURCES

        ENDDO !LOOP OVER ALL OBSERVATION POINTS
      ENDDO !LOOP OVER ALL FREQUENCES

      RETURN
      END
