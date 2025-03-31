*CMZ :  4.00/07 07/04/2020  08.56.40  by  Michael Scheer
*CMZ :  4.00/04 14/05/2019  10.06.52  by  Michael Scheer
*CMZ :  3.03/02 03/12/2015  13.58.49  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.52/16 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.52/00 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.51/00 13/05/2004  12.01.48  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.31.33  by  Michael Scheer
*CMZ : 00.01/06 01/02/95  16.22.16  by  Michael Scheer
*CMZ : 00.01/04 25/01/95  16.36.51  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.19.11  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.43  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE EFOLD
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- THE ROUTINE FOLDS THE SPECTRAL FLUXES WITH A GAUSSIAN TO TAKE THE BEAM
C    ENERGY SPREAD INTO ACCOUNT

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER kfreq
      INTEGER ICEN,IPINOLD
      REAL*4 F,S0,S1,S2,S3

      INTEGER NFOLDP,kmode,iefoldold
      integer :: iwarn=1
      PARAMETER (NFOLDP=1000)
      DOUBLE PRECISION DF,F3SIG,dnsigma
      DOUBLE PRECISION S0E(NDFREQP),S2E(NDFREQP),S3E(NDFREQP),S1E(NDFREQP)
      DOUBLE PRECISION S0EF(NDFREQP),S2EF(NDFREQP),S3EF(NDFREQP),S1EF(NDFREQP)
      double precision ws1(ndfreqp),ws2(ndfreqp)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR EFOLD called to perform energy folding'
      WRITE(LUNGFO,*)
     &  '     Beam energy spread (ESPREAD) and number of sigmas:'
      WRITE(LUNGFO,*)'     ',SNGL(ESPREAD),NSIGE
      WRITE(LUNGFO,*)

      IPINOLD=IPIN

      IF (IPIN.NE.0) THEN
          IPIN=1
          ICEN=ICBRILL
      ELSE    !IPIN
          ICEN=ICBRILL
      ENDIF   !IPIN

10    CONTINUE

      IF (IPIN.EQ.1) THEN

        DO kfreq=1,NFREQ
          F=FREQ(kfreq)
          S0=WSTOKES(1,kfreq)
          S1=WSTOKES(2,kfreq)
          S2=WSTOKES(3,kfreq)
          S3=WSTOKES(4,kfreq)
          S0E(kfreq)=S0
          S1E(kfreq)=S1
          S2E(kfreq)=S2
          S3E(kfreq)=S3
        ENDDO !kfreq

      ELSE !IPIN

        DO kfreq=1,NFREQ
          F=FREQ(kfreq)
          IOBFR=ICEN+NOBSV*(kfreq-1)
          S0=STOKES(1,IOBFR)
          S1=STOKES(2,IOBFR)
          S2=STOKES(3,IOBFR)
          S3=STOKES(4,IOBFR)
          S0E(kfreq)=S0
          S1E(kfreq)=S1
          S2E(kfreq)=S2
          S3E(kfreq)=S3
        ENDDO !kfreq

      ENDIF !IPIN

1     continue
      kmode=-1
      DO kfreq=NFREQEM+1,NFREQ-NFREQEP
        DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
        F3SIG=FREQ(kfreq)*ESPREAD*2.D0
        IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &      .AND.
     &      FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
          if (iefold.ne.3) then
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S0E,F3SIG,NSIGE,FREQ(kfreq),kmode,S0EF(kfreq))
          else
            call util_fold_gauss_lin(nfreq,freq,s0e,f3sig,
     &        dble(nsige),kfreq,s0ef(kfreq),ws1,ws2)
          endif
          kmode=1
          if (s0ef(kfreq).lt.0.0d0) then
            write(lungfo,*)'*** Warning in EFOLD_GAUSS: Negative value of S0_E found.'
            write(lungfo,*)'*** Probably problems with oszillating splines.'
            write(lungfo,*)'*** Check/change or photon energy range and number of intervals.'
            write(lungfo,*)'*** Switching to IEFOLD=3 ***'
            write(6,*)'*** Warning in EFOLD_GAUSS: Negative value of S0_E found.'
            write(6,*)'*** Probably problems with oszillating splines.'
            write(6,*)'*** Check/change photon energy range and number of intervals.'
            write(6,*)'*** Switching to IEFOLD=3 ***'
            iefold=3
            goto 1
          endif
        ELSE
          S0EF(kfreq)=9.999D-11
        ENDIF
      ENDDO !kfreq

      kmode=-1
      DO kfreq=NFREQEM+1,NFREQ-NFREQEP
        DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
        F3SIG=FREQ(kfreq)*ESPREAD*2.D0
        IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &      .AND.
     &      FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
          if (iefold.ne.3) then
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S1E,F3SIG,NSIGE,FREQ(kfreq),kmode,S1EF(kfreq))
          else
            call util_fold_gauss_lin(nfreq,freq,s1e,f3sig,
     &        dble(nsige),kfreq,s1ef(kfreq),ws1,ws2)
          endif
          kmode=1
        ELSE
          S1EF(kfreq)=0.0
        ENDIF
      ENDDO !kfreq

      kmode=-1
      DO kfreq=NFREQEM+1,NFREQ-NFREQEP
        DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
        F3SIG=FREQ(kfreq)*ESPREAD*2.D0
        IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &      .AND.
     &      FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
          if (iefold.ne.3) then
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S2E,F3SIG,NSIGE,FREQ(kfreq),kmode,S2EF(kfreq))
          else
            call util_fold_gauss_lin(nfreq,freq,s2e,f3sig,
     &        dble(nsige),kfreq,s2ef(kfreq),ws1,ws2)
          endif
          kmode=1
        ELSE
          S2EF(kfreq)=0.0
        ENDIF
      ENDDO !kfreq

      kmode=-1
      DO kfreq=NFREQEM+1,NFREQ-NFREQEP
        DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
        F3SIG=FREQ(kfreq)*ESPREAD*2.D0
        IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &      .AND.
     &      FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
          if (iefold.ne.3) then
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S3E,F3SIG,NSIGE,FREQ(kfreq),kmode,S3EF(kfreq))
          else
            call util_fold_gauss_lin(nfreq,freq,s3e,f3sig,
     &        dble(nsige),kfreq,s3ef(kfreq),ws1,ws2)
          endif
          kmode=1
        ELSE
          S3EF(kfreq)=0.0
        ENDIF
      ENDDO !kfreq

      IF (IPIN.EQ.0) THEN

        DO kfreq=1,NFREQ

          IF (kfreq.GT.NFREQEM.AND.kfreq.LE.NFREQ-NFREQEP) THEN
            WSTOKESE(1,kfreq)=S0EF(kfreq)
            WSTOKESE(2,kfreq)=S1EF(kfreq)
            WSTOKESE(3,kfreq)=S2EF(kfreq)
            WSTOKESE(4,kfreq)=S3EF(kfreq)
          ELSE
            WSTOKESE(1,kfreq)=0.0
            WSTOKESE(2,kfreq)=0.0
            WSTOKESE(3,kfreq)=0.0
            WSTOKESE(4,kfreq)=0.0
          ENDIF

          STOKECE(1,kfreq)=WSTOKESE(1,kfreq)
          STOKECE(2,kfreq)=WSTOKESE(2,kfreq)
          STOKECE(3,kfreq)=WSTOKESE(3,kfreq)
          STOKECE(4,kfreq)=WSTOKESE(4,kFREQ)

          stokesE(1,kfreq)=WSTOKESE(1,kfreq)
          stokesE(2,kfreq)=WSTOKESE(2,kfreq)
          stokesE(3,kfreq)=WSTOKESE(3,kfreq)
          stokesE(4,kfreq)=WSTOKESE(4,kFREQ)

        ENDDO !kfreq

      ELSE IF (IPIN.EQ.1) THEN

        DO kfreq=1,NFREQ

          IF (kfreq.GT.NFREQEM.AND.kfreq.LE.NFREQ-NFREQEP) THEN
            WSTOKESE(1,kfreq)=S0EF(kfreq)
            WSTOKESE(2,kfreq)=S1EF(kfreq)
            WSTOKESE(3,kfreq)=S2EF(kfreq)
            WSTOKESE(4,kfreq)=S3EF(kfreq)
          ELSE
            WSTOKESE(1,kfreq)=0.0
            WSTOKESE(2,kfreq)=0.0
            WSTOKESE(3,kfreq)=0.0
            WSTOKESE(4,kfreq)=0.0
          ENDIF

        ENDDO !kfreq

      ELSE IF (IPIN.EQ.2) THEN

        DO kfreq=1,NFREQ

          IF (kfreq.GT.NFREQEM.AND.kfreq.LE.NFREQ-NFREQEP) THEN
            STOKECE(1,kfreq)=S0EF(kfreq)
            STOKECE(2,kfreq)=S1EF(kfreq)
            STOKECE(3,kfreq)=S2EF(kfreq)
            STOKECE(4,kfreq)=S3EF(kfreq)
          ELSE
            STOKECE(1,kfreq)=0.0
            STOKECE(2,kfreq)=0.0
            STOKECE(3,kfreq)=0.0
            STOKECE(4,kfreq)=0.0
          ENDIF

        ENDDO !kfreq

      ENDIF   !IPIN

      IF (IPIN.EQ.1) THEN
        IPIN=2
        GOTO 10 !JUMP BACK TO RERUN FOR CENTER OF PINHOLE
      ELSE IF (IPIN.EQ.2) THEN
        IPIN=1
      ENDIF

      IF (IFOLD.NE.0) THEN

20      CONTINUE

        IF (IPIN.EQ.1) THEN
          DO kfreq=1,NFREQ
            F=FREQ(kfreq)
            S0=WSTOKESF(1,kfreq)
            S1=WSTOKESF(2,kfreq)
            S2=WSTOKESF(3,kfreq)
            S3=WSTOKESF(4,kfreq)
            S0E(kfreq)=S0
            S1E(kfreq)=S1
            S2E(kfreq)=S2
            S3E(kfreq)=S3
          ENDDO !kfreq
        ELSE    !IPIN
          DO kfreq=1,NFREQ
            F=FREQ(kfreq)
            IOBFR=ICEN+NOBSV*(kfreq-1)
            S0=STOKESF(1,IOBFR)
            S1=STOKESF(2,IOBFR)
            S2=STOKESF(3,IOBFR)
            S3=STOKESF(4,IOBFR)
            S0E(kfreq)=S0
            S1E(kfreq)=S1
            S2E(kfreq)=S2
            S3E(kfreq)=S3
          ENDDO !kfreq
        ENDIF   !IPIN

        kmode=-1
        DO kfreq=NFREQEM+1,NFREQ-NFREQEP
          DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(kfreq)*ESPREAD*2.D0
          IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
            if (iefold.ne.3) then
              CALL EFOLD_GAUSS
     &          (NFREQ,FREQ,S0E,F3SIG,NSIGE,FREQ(kfreq),kmode,S0EF(kfreq))
            else
              call util_fold_gauss_lin(nfreq,freq,s0e,f3sig,
     &          dble(nsige),kfreq,s0ef(kfreq),ws1,ws2)
            endif
            kmode=1
          ELSE
            S0EF(kfreq)=9.999D-11
          ENDIF
        ENDDO !kfreq

        kmode=-1
        DO kfreq=NFREQEM+1,NFREQ-NFREQEP
          DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(kfreq)*ESPREAD*2.D0
          IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
            if (iefold.ne.3) then
              CALL EFOLD_GAUSS
     &          (NFREQ,FREQ,S1E,F3SIG,NSIGE,FREQ(kfreq),kmode,S1EF(kfreq))
            else
              call util_fold_gauss_lin(nfreq,freq,s1e,f3sig,
     &          dble(nsige),kfreq,s1ef(kfreq),ws1,ws2)
            endif
            kmode=1
          ELSE
            S1EF(kfreq)=0.0
          ENDIF
        ENDDO !kfreq

        kmode=-1
        DO kfreq=NFREQEM+1,NFREQ-NFREQEP
          DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(kfreq)*ESPREAD*2.D0
          IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
            if (iefold.ne.3) then
              CALL EFOLD_GAUSS
     &          (NFREQ,FREQ,S2E,F3SIG,NSIGE,FREQ(kfreq),kmode,S2EF(kfreq))
            else
              call util_fold_gauss_lin(nfreq,freq,s2e,f3sig,
     &          dble(nsige),kfreq,s2ef(kfreq),ws1,ws2)
            endif
            kmode=1
          ELSE
            S2EF(kfreq)=0.0
          ENDIF
        ENDDO !kfreq

        kmode=-1
        DO kfreq=NFREQEM+1,NFREQ-NFREQEP
          DF=FREQ(kfreq)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(kfreq)*ESPREAD*2.D0
          IF(FREQ(kfreq)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(kfreq)+DF.LE.FREQ(NFREQ)) THEN
            if (iefold.ne.3) then
              CALL EFOLD_GAUSS
     &          (NFREQ,FREQ,S3E,F3SIG,NSIGE,FREQ(kfreq),kmode,S3EF(kfreq))
            else
              call util_fold_gauss_lin(nfreq,freq,s3e,f3sig,
     &          dble(nsige),kfreq,s3ef(kfreq),ws1,ws2)
            endif
            kmode=1
          ELSE
            S3EF(kfreq)=0.0
          ENDIF
        ENDDO !kfreq


        IF (IPIN.EQ.1) THEN

          DO kfreq=1,NFREQ

            IF (kfreq.GT.NFREQEM.AND.kfreq.LE.NFREQ-NFREQEP) THEN
              WSTOKESEF(1,kfreq)=S0EF(kfreq)
              WSTOKESEF(2,kfreq)=S1EF(kfreq)
              WSTOKESEF(3,kfreq)=S2EF(kfreq)
              WSTOKESEF(4,kfreq)=S3EF(kfreq)
            ELSE
              WSTOKESEF(1,kfreq)=0.0
              WSTOKESEF(2,kfreq)=0.0
              WSTOKESEF(3,kfreq)=0.0
              WSTOKESEF(4,kfreq)=0.0
            ENDIF

          ENDDO !kfreq

        ELSE IF (IPIN.EQ.2) THEN

          DO kfreq=1,NFREQ

            IF (kfreq.GT.NFREQEM.AND.kfreq.LE.NFREQ-NFREQEP) THEN
              STOKECEF(1,kfreq)=S0EF(kfreq)
              STOKECEF(2,kfreq)=S1EF(kfreq)
              STOKECEF(3,kfreq)=S2EF(kfreq)
              STOKECEF(4,kfreq)=S3EF(kfreq)
            ELSE
              STOKECEF(1,kfreq)=0.0
              STOKECEF(2,kfreq)=0.0
              STOKECEF(3,kfreq)=0.0
              STOKECEF(4,kfreq)=0.0
            ENDIF

          ENDDO !kfreq

        ENDIF   !IPIN

        IF (IPIN.EQ.1) THEN
          IPIN=2
          GOTO 20 !JUMP BACK TO RERUN FOR CENTER OF PINHOLE
        ELSE IF (IPIN.EQ.2) THEN
          IPIN=1
        ENDIF

      ENDIF !IFOLD

      if (ipin.eq.1) CALL EFOLDPIN !20190510

      IPIN=IPINOLD

      RETURN
      END
