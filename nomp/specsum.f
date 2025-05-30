*CMZ :  4.00/13 12/11/2021  11.33.50  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.28  by  Michael Scheer
*CMZ : 00.02/04 26/02/97  10.39.20  by  Michael Scheer
*CMZ : 00.01/09 20/10/95  15.51.37  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.26.21  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.32  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.22  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPECSUM
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


      INTEGER IFREQ,ISOUR,IOBSV,ICAL,ISS

      DOUBLE PRECISION DZ,DY
      DOUBLE PRECISION SUMT

      DATA ICAL/0/

      ALLOCATE(SPCSMRAT(NSOURCE))
      ALLOCATE(SPCSMSUM(NSOURCE*NOBSV))


      IF (IPIN.EQ.0) GOTO 9999

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR SPECSUM:'
      WRITE(LUNGFO,*)

         IF (IPINCIRC.NE.0.AND.ICAL.NE.1) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING SR SPECSUM ***'
             WRITE(LUNGFO,*)
     &'USE OF FLAG IPINCIRC NOT RECOMMENDED FOR THIS ROUTINE SINCE'
             WRITE(LUNGFO,*)
     &'SHAPE OF CIRCULARE PINHOLE IS TAKEN INTO ACCOUNT ONLY VERY ROUGHLY'
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
C            WRITE(6,*)
C            WRITE(6,*)
C            WRITE(6,*)'*** WARNING SR SPECSUM ***'
C            WRITE(6,*)
C     &'USE OF FLAG IPINCIRC NOT RECOMMENDED FOR THIS ROUTINE SINCE'
C            WRITE(6,*)
C     &'SHAPE OF CIRCULARE PINHOLE IS TAKEN INTO ACCOUNT ONLY VERY ROUGHLY'
C            WRITE(6,*)
C            WRITE(6,*)
             ICAL=1
         ENDIF !IPINCIRC

      WRITE(LUNGFO,*)
     & '     Photon energy and flux through pinhole simply summed up:'
      IF (IW_BLEN.NE.0) THEN
         WRITE(LUNGFO,*)'     (and ratio of spline and summation results)'
         WRITE(6,*)'     Ratio of spline and summation results:'
         write(6,*)
      ENDIF
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)


         DO IFREQ=1,NFREQ

         DO ISOUR=1,NSOURCE

             SPCSMSUM(ISOUR)=0.0

         DO IOBSV=1,NOBSV

            ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
             IF (IPINCIRC.EQ.0) THEN

            IF (
     &           DABS(OBSV(2,IOBSV)-PINCEN(2))-PINH/2.D0.LT.1.D-10
     &                      .AND.
     &           DABS(OBSV(3,IOBSV)-PINCEN(3))-PINW/2.D0.LT.1.D-10
     &                      ) THEN

                  IF(DABS(
     &                           DABS(OBSV(3,IOBSV)-PINCEN(3))
     &                           -PINW/2.D0).LT.1.D-10) THEN

                DZ=OBSVDZ/2.D0

                  ELSE

                     DZ=OBSVDZ

                  ENDIF

                  IF(DABS(
     &                           DABS(OBSV(2,IOBSV)-PINCEN(2))
     &                           -PINH/2.D0).LT.1.D-10) THEN

                DY=OBSVDY/2.D0

                  ELSE

                    DY=OBSVDY

                  ENDIF

                  SPCSMSUM(ISOUR)=SPCSMSUM(ISOUR)
     &                      +SPEC( ILIOBFR)*
     &                                   DZ*DY

            ENDIF   !OBSV

             ELSE    !IPINCIR

            IF (
     &                      (OBSV(2,IOBSV)-PINCEN(2))**2
     &                     +(OBSV(3,IOBSV)-PINCEN(3))**2
     &                       -PINR**2.LT.1.D-10) THEN

                  DZ=OBSVDZ
                  DY=OBSVDY

                  SPCSMSUM(ISOUR)=SPCSMSUM(ISOUR)
     &                      +SPEC( ILIOBFR)*
     &                                   DZ*DY


                 ENDIF  !OBSV

             ENDIF   !IPINCIR

         ENDDO !NOBSV
         ENDDO !ISOUR

         SUMT=0.D0

         DO IOBSV=1,NOBSV

             IF (IPINCIRC.EQ.0) THEN
            IF (
     &           DABS(OBSV(2,IOBSV)-PINCEN(2))-PINH/2.D0.LT.1.D-10
     &                      .AND.
     &           DABS(OBSV(3,IOBSV)-PINCEN(3))-PINW/2.D0.LT.1.D-10
     &                      ) THEN

                  IF(DABS(
     &                           DABS(OBSV(3,IOBSV)-PINCEN(3))
     &                           -PINW/2.D0).LT.1.D-10) THEN
                DZ=OBSVDZ/2.D0

                  ELSE

                     DZ=OBSVDZ

                  ENDIF

                  IF(DABS(
     &                           DABS(OBSV(2,IOBSV)-PINCEN(2))
     &                           -PINH/2.D0).LT.1.D-10) THEN
                DY=OBSVDY/2.D0

                  ELSE

                    DY=OBSVDY

                  ENDIF

                  SUMT=SUMT+SPECTOT(IOBSV+NOBSV*(IFREQ-1))*DY*DZ

            ENDIF   !OBSV

             ELSE    !IPINCIR

            IF (
     &                      (OBSV(2,IOBSV)-PINCEN(2))**2
     &                     +(OBSV(3,IOBSV)-PINCEN(3))**2
     &                       -PINR**2.LT.1.D-10) THEN

                  DZ=OBSVDZ
                  DY=OBSVDY
                  SUMT=SUMT+SPECTOT(IOBSV+NOBSV*(IFREQ-1))*DY*DZ
                 ENDIF  !OBSV

             ENDIF   !IPINCIR

         ENDDO

      IF (IW_BLEN.EQ.0) THEN

         IF (IUNIT.EQ.0)       !260194
     &      WRITE(LUNGFO,*)
     &      '  ',SNGL(FREQ(IFREQ)),(SNGL(SPCSMSUM(ISOUR)),ISOUR=1,NSOURCE)
     &      ,SNGL(SUMT)

         IF (IUNIT.NE.0) !260194
     &      WRITE(LUNGFO,*)
     &      '  ',SNGL(WELLEN(IFREQ)),(SNGL(SPCSMSUM(ISOUR)),ISOUR=1,NSOURCE)
     &      ,SNGL(SUMT)

      ENDIF

        IF (IW_BLEN.NE.0) THEN
           DO ISS=1,NSOURCE
             SPCSMRAT(ISS)=9999.
             IF (SPCSMSUM(ISS).NE.0.D0) SPCSMRAT(ISS)=WFLUX(ISS+NSOURCE*(IFREQ-1))
     &                                     /SPCSMSUM(ISS)
           ENDDO !ISS
           WRITE(LUNGFO,2584)SNGL(FREQ(IFREQ)),(SPCSMRAT(ISS),ISS=1,NSOURCE)
           WRITE(6,2584)SNGL(FREQ(IFREQ)),(SPCSMRAT(ISS),ISS=1,NSOURCE)
        ENDIF   !IW_BLEN
2584     FORMAT('     ',6(1PE12.4))

         ENDDO !IFREQ

9999  CONTINUE
      DEALLOCATE(SPCSMRAT)
      DEALLOCATE(SPCSMSUM)
      RETURN
      END
