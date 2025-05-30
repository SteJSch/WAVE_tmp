*CMZ :  3.07/00 15/03/2019  13.09.52  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.51/02 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.17/00 03/11/2000  14.00.48  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.31.33  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.20.25  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.43  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.09  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WFOLDSTO
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

C--- CALCULATES FOLDING OF STOKES INTENSITIY WITH ELECTRON PHASE SPACE
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

      INTEGER ISTOK,KFREQ,IZ,IY

C--- CALCULATE FOURIER-COEFFICIENTS OF GAUSSIAN

C     CALL WGFOUR !ALREADY DONE, USE VALUES OF ISIGSTO

      DO kfreq=1,NFREQ

        DO ISTOK=1,4

C--- CALCULATE 2D POLYNOMIALS OF INTENSITY DISTRIBUTION


          IF (ISPECANAF.EQ.0.AND.IFOLD.EQ.-2) CALL WPOLY2ST(ISTOK,kfreq)

C--- PERFORM FOLDING

          CALL WFOLISTO(ISTOK,kfreq) !MUST ALSO BE CALLED FOR ISPECANAF

C--- CALCULATE INTEGRATED INTENSITY

          IF(ISPECANAF.NE.0) CALL SPECANAF    !MUST BE CALLED EACH TIME

          IF (IPINCIRC.EQ.0.OR.IPINCIRC*IRPHI.NE.0)
     &      CALL PSPLSTOF(ISTOK,kfreq)

          CALL BLENSTOF(ISTOK,kfreq)

        ENDDO !ISTOK

        DO ISTOK=1,4

C--- DELETE INTENSITY IN EDGES

          IF (IPINCIRC.EQ.0) THEN

            DO IY=1,NOBSVY
              DO IZ=1,NOBSVZ

                IF (
     &              IY.LT.(NOBSVY-MOBSVY)/2+1
     &              .OR.IY.GT.(NOBSVY-MOBSVY)/2+MOBSVY
     &              .OR.IZ.LT.(NOBSVZ-MOBSVZ)/2+1
     &              .OR.IZ.GT.(NOBSVZ-MOBSVZ)/2+MOBSVZ
     &              ) THEN

                  STOKESF(ISTOK,(IY-1)*NOBSVZ+IZ+NOBSV*(kfreq-1))=0.0

                ENDIF

              ENDDO !IZ
            ENDDO !IY

          ELSE  !IPINCIRC          !CIRCULARE PINHOLE


            DO IY=1,NOBSVY
              DO IZ=1,NOBSVZ

                IF (
     &              (OBSVZ(IZ)-PINCEN(3))**2
     &              +(OBSVY(IY)-PINCEN(2))**2
     &              -PINR**2
     &              .GT.1.D-10
     &              ) THEN

                  STOKESF(ISTOK,(IY-1)*NOBSVZ+IZ+NOBSV*(kfreq-1))=0.0

                ENDIF

              ENDDO !IZ
            ENDDO !IY

          ENDIF !IPINCIRC
        ENDDO !ISTOK

      ENDDO !kfreq

      RETURN
      END
