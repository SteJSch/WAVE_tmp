*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.33/06 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.46  by  Michael Scheer
*CMZ :  2.16/05 04/08/2000  11.03.38  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/07 17/02/2000  15.11.12  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.01.43  by  Michael Scheer
*CMZ :  2.10/01 25/02/99  15.57.16  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.41  by  Michael Scheer
*-- Author : Michael Scheer

      SUBROUTINE BLENDSPECIV(ISOUR,IZ)
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- INTEGRATES THE SPLINES THAT INTERPOLATE THE POWER INSIDE THE PINHOLE

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
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.


      INTEGER ISOUR,IY,IZ,IOBSV
      INTEGER ICAL,IWBLEN,IDUM
      INTEGER IWSOUR
      DOUBLE PRECISION DSUM

      DOUBLE PRECISION SUMZ(NDOBSVYP),S2(NDOBSVYP),SUM
      DOUBLE PRECISION SUMZP(NDOBSVYP),S2P(NDOBSVYP),SUMP,DSUMP

      DATA ICAL/0/

      IWBLEN=0

      IF (IPINCIRC.EQ.0) THEN

C--- INTEGRATION ALONG VERTICAL AXIS Y

      DO IY=1,NOBSVY
          IOBSV=IZ+(IY-1)*NOBSVZ
          SUMZ(IY)=SPECI(ISOUR+NSOURCE*(IOBSV-1))
      ENDDO

      CALL FSPLINDX(OBSVDY,SUMZ,NOBSVY,0.D0,0.D0,S2)

        IF(MOBSVY.GT.1) THEN

          SUM=0.0
          SUMP=0.0

          DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY-1

            DSUM=
     &          OBSVDY*0.5D0
     &          *(SUMZ(IY)+SUMZ(IY+1))
     &          -OBSVDY**3/24.D0
     &          *(S2(IY)+S2(IY+1))

          IF (IW_BLEN.NE.0) THEN
                 DSUMP=
     &          OBSVDY*0.5D0
     &          *(SUMZP(IY)+SUMZP(IY+1))
     &          -OBSVDY**3/24.D0
     &          *(S2P(IY)+S2P(IY+1))
          ENDIF

             IF(
     &              (IWSOUR.NE.ISOUR)
     &              .AND.
     &              DSUM.LT.0.0) THEN
         IF (IW_BLEN.EQ.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING SR BLENDSPECIV ***'
              WRITE(LUNGFO,*)
     &              'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              ENDIF

             IWSOUR=ISOUR
           IW_BLEN=1
           IWBLEN=1
           DO IDUM=1,NOBSVY
         S2P(IDUM)=S2(IDUM)
         SUMZP(IDUM)=SUMZ(IDUM)
         SUMP=SUM
         DSUMP=DSUM
           ENDDO

            ENDIF !IWSOUR

          IF (IWBLEN.NE.0) SUMP=SUMP+DABS(DSUMP)

             SUM=SUM+DSUM

          ENDDO !IY

        ELSE !MOBSVY

         SUM=OBSVDY*SUMZ(NOBSVY/2+1)

         IF (IWBLEN.NE.0) SUMP=OBSVDY*SUMZP(NOBSVY/2+1)

        ENDIF !MOBSVY

      ELSE  !IPINCIRC

      WRITE(LUNGFO,*)'*** WARNING IN BLENDSPECIV: ***'
      WRITE(LUNGFO,*)'INTEGRATION OF POWERDENSITY NOT POSSIBLE'
      WRITE(LUNGFO,*)'FOR CIRCULAR PINHOLE, SORRY'
      WRITE(LUNGFO,*)'USE PAW AND NTUPLE FOR RAW OFFLINE INTEGRATION'
      WRITE(6,*)'*** WARNING IN BLENDSPECIV: ***'
      WRITE(6,*)'INTEGRATION OF POWERDENSITY NOT POSSIBLE'
      WRITE(6,*)'FOR CIRCULAR PINHOLE, SORRY'
      WRITE(6,*)'USE PAW AND NTUPLE FOR RAW OFFLINE INTEGRATION'
      RETURN

      ENDIF !IPINCIRC

      SPECIV(ISOUR+NSOURCE*(IZ-1))=SUM

      IF (IWBLEN.NE.0) THEN
          IF (ICAL.EQ.0) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** SR BLENDSPECIV:'
         WRITE(LUNGFO,*)
     &'LINES INDICATED BY * SHOW A RAW ESTIMATE OF ERRORS DUE TO'
         WRITE(LUNGFO,*)
     &'SPLINE FAILURE IF REL. ERROR .GT. 1E-5 (FIRST NUMBER IS SOURCE)'
         WRITE(LUNGFO,*)
         ICAL=1
          ENDIF   !ICAL
          IF (SUMP.NE.0.D0) THEN
         DSUM=SUM/SUMP
          ELSE
         DSUM=-9999.
          ENDIF
          IF (DABS(DSUM-1.D0).GT.1.D-5) THEN
            WRITE(LUNGFO,*)'*',ISOUR,IZ,
     &        SNGL(SUM),SNGL(SUMP),SNGL(DSUM)
            WRITE(6,*)'*',ISOUR,IZ,
     &        SNGL(SUM),SNGL(SUMP),SNGL(DSUM)
          ENDIF
      ENDIF !IWBLEN

      RETURN
      END
