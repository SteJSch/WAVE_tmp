*CMZ :  4.01/05 19/04/2024  10.29.48  by  Michael Scheer
*CMZ :  4.01/04 14/11/2023  13.37.28  by  Michael Scheer
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
      SUBROUTINE ampfold
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEND.

C--- CALCULATES FOLDING OF FIELD AMPLITUDES WITH ELECTRON PHASE SPACE
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

      INTEGER IFREQ,IZ,IY,IOBSV,ICOMP,ireim,ieb,i1,i2

C--- CALCULATE FOURIER-COEFFICIENTS OF GAUSSIAN

c      CALL WGFOUR ! already done in WFOLD

      IF (nsource.gt.1) then
        write(6,*)' '
        write(6,*)' *** Warning in AMPFOLD: More than one source.'
        write(6,*)' Sigmas for folding procedure are taken for source',nsource/2+1
        write(6,*)' '
        write(lungfo,*)' '
        write(lungfo,*)' *** Warning in AMPFOLD: More than one source.'
        write(lungfo,*)' Sigmas for folding procedure are taken for source',nsource/2+1
        write(lungfo,*)' '
      endif

      IF (IFOLD.EQ.-2) then
        write(6,*)' '
        write(6,*)' *** Warning in AMPFOLD: Modus IFOLD.EQ.-2 not available.'
        write(6,*)' '
        write(lungfo,*)' '
        write(lungfo,*)' *** Warning in AMPFOLD: Modus IFOLD.EQ.-2 not available.'
        write(lungfo,*)' '
      endif

c14.11,2023      DO icomp=2,3
      do ieb=1,2

        if (ieb.eq.1) then
          i1=1
          i2=3
        else
          i1=6
          i2=8
        endif

        DO icomp=i1,i2
          DO ireim=1,2
            DO IFREQ=1,NFREQ
C--- PERFORM FOLDING
              CALL AFOLINT(icomp,ireim,IFREQ)
            ENDDO !IFREQ
          ENDDO !ireim
        ENDDO !icomp

      enddo !ieb=1,2

      do ieb=1,2

        if (ieb.eq.1) then
          i1=1
          i2=3
        else
          i1=6
          i2=8
        endif

        DO icomp=i1,i2
          DO ireim=1,2

            DO IFREQ=1,NFREQ

c DELETE INTENSITY IN EDGES

              DO IY=1,NOBSVY
                DO IZ=1,NOBSVZ

                  IOBSV=NOBSVZ*(IY-1)+IZ
                  IOBFR=IOBSV+NOBSV*(IFREQ-1)

                  IF (IPINCIRC.EQ.0) THEN

                    IF (
     &                  IY.LT.(NOBSVY-MOBSVY)/2+1
     &                  .OR.IY.GT.(NOBSVY-MOBSVY)/2+MOBSVY
     &                  .OR.IZ.LT.(NOBSVZ-MOBSVZ)/2+1
     &                  .OR.IZ.GT.(NOBSVZ-MOBSVZ)/2+MOBSVZ
     &                  ) THEN
                      reaima(icomp+2,ireim,iobfr)=0.0d0
                    ENDIF

                  ELSE  !IPINCIRC

                    IF (
     &                  (OBSVZ(IZ)-PINCEN(3))**2
     &                  +(OBSVY(IY)-PINCEN(2))**2
     &                  -PINR**2
     &                  .GT.1.D-10
     &                  ) THEN
                      reaima(icomp+2,ireim,iobfr)=0.0d0
                    ENDIF

                  ENDIF !IPINCIRC

                ENDDO !IZ
              ENDDO !IY

            ENDDO !IFREQ

          ENDDO !ireim
        ENDDO !icomp

      enddo !ieb=1,2

      RETURN
      END
