*CMZ :  4.00/15 13/02/2022  18.53.40  by  Michael Scheer
*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 07/12/2021  18.47.10  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  10.38.43  by  Michael Scheer
*CMZ :  2.64/05 02/09/2009  09.39.09  by  Michael Scheer
*CMZ :  2.41/10 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.20/01 20/11/2000  16.56.14  by  Michael Scheer
*CMZ :  2.17/00 03/11/2000  14.52.16  by  Michael Scheer
*CMZ :  2.16/08 25/10/2000  12.21.53  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.16/03 16/06/2000  12.15.50  by  Michael Scheer
*CMZ :  2.13/00 02/12/99  11.42.16  by  Michael Scheer
*CMZ :  1.03/03 27/03/98  13.56.06  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.01/06 13/02/95  10.25.06  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.22  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.31  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HSOURCE
*KEEP,gplhint.
*KEND.

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

C--- HISTOGRAM OF SOURCES

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

      INTEGER NTUPP,IROI,IROIO,nallo
      PARAMETER (NTUPP=10)
      REAL*8 TUP(NTUPP)

      CHARACTER(3) CHTAGS(NTUPP)
      data chtags/'x','y','z','bx','by','bz','is','ir','roi','bou'/

      INTEGER I,NBIN,ISOUR
      INTEGER ICYCLE
      REAL*4 XI,XE,X
      REAL*8 Z,Y,BX,BY,BZ

      nallo=nco/ihtrsmp+1024
      CALL hbookm(NIDMINI,'SOURCES',NTUPP,'//WAVE',nallo,CHTAGS)

      IROIO=0
      DO ISOUR=1,NSOURCE
        DO I=1,NCO
          X=WTRA(1,1,I)
          IF (X.GE.SOURCEA(1,1,ISOUR).AND.X.LE.SOURCEE(1,1,ISOUR)
     &        .AND.
     &        X.GE.XIANF.AND.X.LE.XIEND) THEN
            TUP(1)=WTRA(1,1,I)
            TUP(2)=WTRA(2,1,I)
            TUP(3)=WTRA(3,1,I)
            TUP(4)=WTRA(1,3,I)
            TUP(5)=WTRA(2,3,I)
            TUP(6)=WTRA(3,3,I)
            TUP(7)=ISOUR
            TUP(8)=0.0
            DO IROI=1,NROIA-1
              IF (TUP(1).GE.ROIX(IROI).AND.TUP(1).LT.ROIX(IROI+1)) THEN
                TUP(8)=IROI
                TUP(9)=ROIP(IROI)
                GOTO 10
              ENDIF
            ENDDO
10          IF (IROI.NE.IROIO) THEN
              TUP(10)=1.0d0
            ELSE
              TUP(10)=0.0d0
            ENDIF
            CALL hfm(NIDMINI,TUP)
            IROIO=IROI
          ENDIF
        ENDDO !I
      ENDDO   !ISOUR

      IF (IHTRACKM.LT.0) RETURN

      IF(NSOURCE.GT.99) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN HSOURCE ***'
        WRITE(LUNGFO,*)'TO MUCH SOURCES, HISTOGRAM IDs COLLIDE'
        WRITE(LUNGFO,*)'CHANGE HISTOGRAM IDENTIFIER AND THIS CHECK'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN HSOURCE ***'
        WRITE(6,*)'TO MUCH SOURCES, HISTOGRAM IDs COLLIDE'
        WRITE(6,*)'CHANGE HISTOGRAM IDENTIFIER AND THIS CHECK'
        WRITE(6,*)
        STOP
      ENDIF !NSOURCE

      IF (XSTARTH.EQ.9999.) THEN
        XI=WTRA(1,1,1)-DS0/2.
      ELSE
        XI=XSTARTH-DS0/2.
      ENDIF !XSTARTH

      IF (XSTOPH.EQ.9999.) THEN
        XE=WTRA(1,1,NCO)+DS0/2.
      ELSE
        XE=XSTOPH+DS0/2.
      ENDIF !XSTARTH
      NBIN=(XE-XI)/DS0
      XE=DS0*NBIN+XI

C--- LOOP OVER ALL SOURCES

      DO ISOUR=1,NSOURCE

      call hbook1m(IDTRCKSZ+ISOUR,'Z OF SOURCE',NBIN,XI,XE,VMX)
      call hbook1m(IDTRCKSY+ISOUR,'Y OF SOURCE',NBIN,XI,XE,VMX)

      call hbook1m(IDSBX+ISOUR,'Bx OF SOURCE',NBIN,XI,XE,VMX)
      call hbook1m(IDSBY+ISOUR,'By OF SOURCE',NBIN,XI,XE,VMX)
      call hbook1m(IDSBZ+ISOUR,'Bz OF SOURCE',NBIN,XI,XE,VMX)

      XI=WTRA(1,1,1)-DS0
      DO I=1,NCO
          X=DS0*I+XI
          IF (X.GE.SOURCEA(1,1,ISOUR).AND.X.LE.SOURCEE(1,1,ISOUR)) THEN

             Y=WTRA(2,1,I)
             Z=WTRA(3,1,I)
             BX=WTRA(1,3,I)
             BY=WTRA(2,3,I)
             BZ=WTRA(3,3,I)

          ELSE    !X

             Y=-9.
             Z=-9.
             BX=-99.
             BY=-99.
             BZ=-99.
          ENDIF   !X

          CALL hfillm(IDTRCKSZ+ISOUR,X,0.,Z)
          CALL hfillm(IDTRCKSY+ISOUR,X,0.,Y)
          CALL hfillm(IDSBX+ISOUR,X,0.,BX)
          CALL hfillm(IDSBY+ISOUR,X,0.,BY)
          CALL hfillm(IDSBZ+ISOUR,X,0.,BZ)

      ENDDO !I

        CALL MHROUT(IDTRCKSZ+ISOUR,ICYCLE,' ')
        CALL MHROUT(IDTRCKSY+ISOUR,ICYCLE,' ')
        CALL MHROUT(IDSBX+ISOUR,ICYCLE,' ')
        CALL MHROUT(IDSBY+ISOUR,ICYCLE,' ')
        CALL MHROUT(IDSBZ+ISOUR,ICYCLE,' ')

      CALL hdeletm(IDTRCKSZ+ISOUR)
      CALL hdeletm(IDTRCKSY+ISOUR)
      CALL hdeletm(IDSBX+ISOUR)
      CALL hdeletm(IDSBY+ISOUR)
      CALL hdeletm(IDSBZ+ISOUR)

      ENDDO   !ISOUR

      RETURN
      END
