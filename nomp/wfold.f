*CMZ :  3.07/00 15/03/2019  12.50.12  by  Michael Scheer
*CMZ :  3.03/02 21/03/2016  16.06.27  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/03 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.51/02 08/10/2009  09.58.11  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.46  by  Michael Scheer
*CMZ :  2.13/10 14/04/2000  14.26.49  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  14.54.46  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.28  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.16  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WFOLD
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- CALCULATES FOLDING OF PINHOLE INTENSITIY WITH ELECTRON PHASE SPACE
C    DISTRIBUTIONS (GAUSSIAN DISTRIBUTION)

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER ISOUR,IFREQ,IZ,IY,IOBSV

C--- CALCULATE FOURIER-COEFFICIENTS OF GAUSSIAN

      CALL WGFOUR

      DO IFREQ=1,NFREQ

        DO IY=1,NOBSVY
          DO IZ=1,NOBSVZ
            IOBSV=(IY-1)*NOBSVZ+IZ
            SPECTOTF(IOBSV+NOBSV*(IFREQ-1))=0.0d0
          ENDDO   !IZ
        ENDDO   !IY

        WFLUXTF(IFREQ)=0.0

        DO ISOUR=1,NSOURCE

C--- CALCULATE 2D POLYNOMIALS OF INTENSITY DISTRIBUTION

          IF (IFOLD.EQ.-2) CALL WPOLY2(ISOUR,IFREQ)

C--- PERFORM FOLDING

          CALL WFOLINT(ISOUR,IFREQ)

C--- CALCULATE INTEGRATED INTENSITY

          IF (ISPECANAF.NE.0) THEN

            CALL SPECANAF

          ENDIF !ISPECANAF

          IF (IPINCIRC.EQ.0.OR.IPINCIRC*IRPHI.NE.0)
     &      CALL PSPLINEF(ISOUR,IFREQ)
          CALL BLENDEF(ISOUR,IFREQ)

          DO IY=1,NOBSVY
            DO IZ=1,NOBSVZ
              IOBSV=(IY-1)*NOBSVZ+IZ
              IOBFR=IOBSV+NOBSV*(IFREQ-1)
              SPECTOTF(IOBFR)=SPECTOTF(IOBFR)+
     &          SPECF(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
            ENDDO   !IZ
          ENDDO   !IY

          WFLUXTF(IFREQ)=WFLUXTF(IFREQ)+WFLUXF(ISOUR+NSOURCE*(IFREQ-1))

C--- DELETE INTENSITY IN EDGES

          DO IY=1,NOBSVY
            DO IZ=1,NOBSVZ

              IF (IPINCIRC.EQ.0) THEN

                IF (
     &              IY.LT.(NOBSVY-MOBSVY)/2+1
     &              .OR.IY.GT.(NOBSVY-MOBSVY)/2+MOBSVY
     &              .OR.IZ.LT.(NOBSVZ-MOBSVZ)/2+1
     &              .OR.IZ.GT.(NOBSVZ-MOBSVZ)/2+MOBSVZ
     &              ) THEN
                  SPECF(ISOUR+NSOURCE*(((IY-1)*NOBSVZ+IZ)-1+NOBSV*(IFREQ-1)))=0.0d0
                ENDIF

              ELSE  !IPINCIRC

                IF (
     &              (OBSVZ(IZ)-PINCEN(3))**2
     &              +(OBSVY(IY)-PINCEN(2))**2
     &              -PINR**2
     &              .GT.1.D-10
     &              ) THEN
                  SPECF(ISOUR+NSOURCE*(((IY-1)*NOBSVZ+IZ)-1+NOBSV*(IFREQ-1)))=0.0d0
                ENDIF

              ENDIF !IPINCIRC

            ENDDO !IZ
          ENDDO !IY

        ENDDO !ISOUR
      ENDDO !IFREQ

      DO IFREQ=1,NFREQ

        DO IY=1,NOBSVY
          DO IZ=1,NOBSVZ

            IF (IPINCIRC.EQ.0) THEN

              IF (
     &            IY.LT.(NOBSVY-MOBSVY)/2+1
     &            .OR.IY.GT.(NOBSVY-MOBSVY)/2+MOBSVY
     &            .OR.IZ.LT.(NOBSVZ-MOBSVZ)/2+1
     &            .OR.IZ.GT.(NOBSVZ-MOBSVZ)/2+MOBSVZ
     &            ) THEN

                SPECTOTF((IY-1)*NOBSVZ+IZ+NOBSV*(IFREQ-1))=0.0d0

              ENDIF

            ELSE  !IPINCIRC

              IF (
     &            (OBSVZ(IZ)-PINCEN(3))**2
     &            +(OBSVY(IY)-PINCEN(2))**2
     &            -PINR**2
     &            .GT.1.D-10
     &            ) THEN

                SPECTOTF((IY-1)*NOBSVZ+IZ+NOBSV*(IFREQ-1))=0.0d0

              ENDIF

            ENDIF !IPINCIRC

          ENDDO !IZ
        ENDDO !IY
      ENDDO !IFREQ

      RETURN
      END
